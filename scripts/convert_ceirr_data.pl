#!/usr/bin/perl
#======================================================================
# convert_data.pl -t [human/animal/seroloy/viral]
#======================================================================
use strict;
use Data::Dumper;
use Benchmark;
use Getopt::Std;
use POSIX qw/strftime/;
use File::Basename;
use JSON;
use Text::CSV;

use Bio::BVBRC::CEIRR::Config; 

sub usage{
  my $str = <<STR;
Usage: convert_data.pl -f file -t type[human/animal/seroloy/viral]
STR

  print STDERR $str;
}

my $para;
my %options;
my $processed_file = '';
my $type = '';
my $delimiter = ',';

getopts("f:t:",\%options);
if ( !defined $options{t} ){
  usage();
  exit(1);
} else{
  $processed_file = $options{f};
  $type = $options{t};
}

my $json = JSON->new->allow_nonref;
my $fdate = strftime('%Y%m%d',localtime);

$para = new Bio::BVBRC::CEIRR::Config; 
$para->setType($type);

my $jobDir = dirname($processed_file);

my $error_log_file   = "$jobDir/error\_$type.log";
my $warning_log_file = "$jobDir/warning\_$type.log";
my $run_log_file     = "$jobDir/run\_$type.log";
my $out_file = "${type}_${fdate}.json";

my %look    = $para->getLookUp();
my %country = $para->getCountryCode();
my %country_mapping = $para->getCountryMapping();

my %header = $para->getFileHeader($type);
my %map = $para->getMapToSolr();
  
print STDERR "Convert file $processed_file to json format\n";
 
if ( -e "${processed_file}" ){
  my @records = ();
  my @attribs = ();
  my @solr    = ();
  my @values  = ();
    
  my $check = "P";

  my $f_stats = (stat("$processed_file"))[9];
  my $f_stamp = scalar localtime  $f_stats;

  my $count = 0;

  # define csv parser objext
  my $csv = Text::CSV->new({ sep_char => $delimiter});

  open(FH, '<', "${processed_file}") or die "Failed to open $processed_file: $!";
  while( my $entry = <FH> ){
    chomp $entry;
    $count++;

    if ($csv->parse($entry)) {
      if( $count==1 ){
        $entry=~s/ *\/ */_/g;
        $entry=~s/ +/_/g;

        @attribs = $csv->fields(); 

        for (my $i=0; $i<scalar @attribs; $i++){
          if ( lc $attribs[$i] ne "row" ) {
            if( $map{lc $attribs[$i]}->{'solr'} ){
              # valid attribute name
              $solr[$i] =  $map{lc $attribs[$i]}->{'solr'};
            }else{
              print STDERR " Invalid attribute name for solr: $attribs[$i] for $map{lc $attribs[$i]}->{'solr'}\n";
              $check = "F";
            } 
          }
        }
      }else{
        my $record;
        @values = $csv->fields(); 
        #next unless scalar @attribs == scalar @values;
    
        for (my $i=0; $i<scalar @attribs; $i++){
          if ( lc $attribs[$i] ne "row" ) {
            # skip if value is null or na
            $values[$i]=~ s/^ *| *$//g;
            $values[$i]=~ s/  */ /g;
            $values[$i]=~ s/ *:  */:/g;
            $values[$i]=~ s/^,|,$//g;
            $values[$i]=~ s/^\,|\,$//g;
            $values[$i]=~ s/^\,\,|\,\,$//g;
            next if $values[$i]=~/^(-* *-*|null)$/i;
            next if ( !$values[$i] );
       
            my $new_values; 
            if ( $map{lc $attribs[$i]}->{'role'} ){
              $new_values = lookup($type, $solr[$i], "$values[$i]");
              if($map{lc $attribs[$i]}->{'role'} eq "Influenza"){
                push @{$record->{$solr[$i]}}, "$new_values";
              }else{
                push @{$record->{$solr[$i]}}, $map{lc $attribs[$i]}->{'role'}.":"."$new_values";
              }
            }else{
              if (lc $attribs[$i] eq "collection_date" || lc $attribs[$i] eq "submission_date" || lc $attribs[$i] eq "embargo_end_date" ){
                if ($values[$i] and $values[$i] =~ /^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z$/){
                }else{
                   $values[$i] = '';
                }
              }
              if ($type =~ /human/ && lc $attribs[$i] eq "host_group" && $values[$i] eq "human"){
                $record->{host_group}= "Human";
                $record->{host_species} = "Homo sapiens";
                $record->{host_common_name} = "Human";
              }else{
                if ( $map{lc $attribs[$i]}->{'flag'} eq "Y" ){
                  my @arr;
                  if (lc $attribs[$i] eq "taxon_lineage_ids"){
                    @arr = split /\s/, $values[$i];
                  }else{
                    @arr = split /,/, $values[$i];
                  }
                  foreach my $av (@arr){
                    $av =~ s/^\s+//g;
                    $av =~ s/\s+$//g;
                    push @{$record->{$solr[$i]}},$av;
                  }
                }else{
                  #Ignore collection_latitude and collection_longitude if NA, since SOLR accepts only number
                  if (not ((lc $attribs[$i] eq "collection_latitude" or lc $attribs[$i] eq "collection_longitude") and $values[$i] eq "NA")) {
                    $new_values = lookup($type, $solr[$i], "$values[$i]");
                    $record->{$solr[$i]} = $new_values;

                    #Add geographic_group
                    if (lc $attribs[$i] eq "collection_country"){
                      $record->{"geographic_group"} = $country_mapping{group}->{$values[$i]};
                    }
                  }
                }
              }
            }
          }
        }
        push @records, $record;
      }
    }else{
      print STDERR "WARNING: Line could not be parsed: $entry\n";
    }
  }
  close FH;
  if( $check eq "P" ){ 
    my $out = $json->pretty->encode(\@records);
    #my $out_file = "${type}_${fdate}.json";
    $para->write_to_file($out_file,$out,"new"); 
  }
}
      
print STDERR "Finished $processed_file converting.\n";

# Return json file
print $out_file;

1;

sub lookup {
   my ($type,$attr,$data) = @_;
   
   my %common = ( 'Y'   => 'Yes',
     'N'   => 'No',
     'NON' => 'None',
     'U'   => 'Unknown',
     'NA'  => 'Not Applicable',
   );

   my $new_v; 
   if($attr eq "collection_country" ){
      if( $country{$data} ){
         $new_v = $country{$data}.",";
      }else{
         $new_v = $data.",";
      }
   }elsif($attr eq "host_age" ){
      my %age = ( 'A'   => 'Adult',
                  'ADL' => 'Adult',
                  'H'   => 'Hatch year',
                  'J'   => 'Juvenile',
                  'JUV' => 'Juvenile',
                  'PUP' => 'Pup',
                  'SAD' => 'Subadult',
                  'WNL' => 'Weanling',
                  'U'   => 'Unknown/undetermined',
                );
      if( $age{$data} ){
         $new_v = $age{$data}.",";
      }else{
         $new_v = $data.",";
      }
   }else{
      my @vals =  split /,/, $data;
      foreach (@vals){
         if( $look{$attr}->{$_} ){
            $new_v .= $look{$attr}->{$_}.",";
         }elsif( $common{$_} ){
            $new_v .= $common{$_}.",";
         }else{
            $new_v .= $_.",";
         }
      }
   } 
   chop $new_v;
   return $new_v;
}
