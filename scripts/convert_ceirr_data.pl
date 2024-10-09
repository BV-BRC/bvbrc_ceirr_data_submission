#!/usr/bin/perl
#======================================================================
# convert_ceirr_data.pl -f file -t [human/human-sh/animal/seroloy/viral] -n number_of_data
#======================================================================
use strict;
use Data::Dumper;
use Benchmark;
use Getopt::Std;
use POSIX qw/strftime/;
use File::Basename;
use File::Copy ('copy', 'move');
use JSON;
use Text::CSV;
use Fcntl qw(:flock SEEK_END);

use Bio::BVBRC::CEIRR::Config; 

sub usage{
  my $str = <<STR;
Usage: convert_ceirr_data.pl -f file -t type[human/human-sh/animal/seroloy/viral] -n number_of_data
STR

  print STDERR $str;
}

my $para;
my %options;
my $processed_file = '';
my $type = '';
my $number_of_data;
my $delimiter = ',';

getopts("f:t:n:",\%options);
if ( !defined $options{t} ){
  usage();
  exit(1);
} else{
  $processed_file = $options{f};
  $type = $options{t};
  $number_of_data = $options{n};
}

my $json = JSON->new->allow_nonref;
my $fdate = strftime('%Y%m%d',localtime);

$para = new Bio::BVBRC::CEIRR::Config; 
$para->setType($type);

my $jobDir = dirname($processed_file);

my $error_log_file   = "$jobDir/error\_$type.log";
my $warning_log_file = "$jobDir/warning\_$type.log";
my $run_log_file     = "$jobDir/run\_$type.log";

# create bvbrc accession file
my $bvbrc_accession_file = "$jobDir/BVBRC_Accession_ID.csv";
my $accession_info = "";
if( $type =~ /human/ or $type eq "animal" ){
  $accession_info .= "Row,Sample_Identifier,BVBRC_Accession_ID\n";
}else{
  $accession_info .= "Row,Sample_Identifier,BVBRC_Accession_ID,Virus_Identifier\n";
}
$para->write_to_file($bvbrc_accession_file ,"$accession_info", "new");

my $sequence_id_file = "/vol/bvbrc/production/application-backend/bvbrc_ceirr_data_submission/sequence_id";

open(FH, '+<', $sequence_id_file) or die "Failed to open $sequence_id_file: $!";
flock(FH, LOCK_EX) or die "Cannot lock $sequence_id_file: $!";
my $seq = <FH>;
seek FH, 0, 0;
truncate FH, 0;
print FH $seq + $number_of_data;
close FH;

my $out_file = "${type}_${fdate}.json";

my %look    = $para->getLookUp();
my %country = $para->getCountryCode();
my %country_mapping = $para->getCountryMapping();
my %datasource = $para->getDataSourceCode();
#To be able to lookup by inst name
my %ds_by_value = reverse %datasource;

my %header = $para->getFileHeader($type);
my %map = $para->getMapToSolr();

print STDERR "Convert file $processed_file of valid $number_of_data data to json format\n";
 
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
  my $csv = Text::CSV->new({ sep_char => $delimiter, eol => $/ });

  # Tmp file to append BVBRC Accession ID
  my $processed_file_tmp = "${processed_file}.tmp";
  open(my $temp_file, ">", $processed_file_tmp) or die "Failed to open $processed_file_tmp: $!";

  open(FH, '<', "${processed_file}") or die "Failed to open $processed_file: $!";

  my %accession_id_map = {};
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
        push(@attribs, ("BVBRC_Accession_ID"));
        $csv->print($temp_file, \@attribs);
      }else{
        my $record;
        @values = $csv->fields(); 

        my $row_number;
        my $contributing_institution;
        my $sampleid;
        my $virusid;
        #next unless scalar @attribs == scalar @values;
        for (my $i=0; $i < scalar(@attribs) - 1; $i++){ #One less to avoid extra empty column for BVBRC_Accession_ID
          if ( lc $attribs[$i] ne "row" ) {
            # skip if value is null or na
	    $values[$i]  =~ s/\"//g;
            my $value_check = $values[$i];
            $value_check =~ s/^ *| *$//g;
            $value_check =~ s/  */ /g;
            $value_check =~ s/ *:  */:/g;
            $value_check =~ s/^,|,$//g;
            $value_check =~ s/^\,|\,$//g;
            $value_check =~ s/^\,\,|\,\,$//g;
            next if $value_check =~/^(-* *-*|null)$/i;
            next if ( !$value_check );
       
            my $new_values; 
            if ( $map{lc $attribs[$i]}->{'role'} ){
              $new_values = lookup($solr[$i], "$values[$i]");
              push @{$record->{$solr[$i]}}, $map{lc $attribs[$i]}->{'role'}.":"."$new_values";
            }else{
              if (lc $attribs[$i] eq "collection_date" || lc $attribs[$i] eq "submission_date" || lc $attribs[$i] eq "embargo_end_date" ){
                if ($values[$i] and $values[$i] =~ /^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z$/){
                }else{
                   $values[$i] = '';
                }
              }
	      next if ( !$values[$i] );

              if ($type =~ /human/ && lc $attribs[$i] eq "host_group"){
                $record->{host_group} = "Human";
		if ( $type eq "human-sh" ){
                   $record->{host_group} = "Human Southern Hemisphere";
                }
                $record->{host_species} = "Homo sapiens";
                $record->{host_common_name} = "Human";
              }elsif ($type eq 'animal' && lc $attribs[$i] eq "host_group"){
                $record->{host_group} = "Animal";
              }else{
                if ( $map{lc $attribs[$i]}->{'flag'} eq "Y" ){
                  my @arr;
                  if (lc $attribs[$i] eq "taxon_lineage_ids"){
                    @arr = split /\s/, $values[$i];
		  }elsif ( $solr[$i] eq "host_race" ){
		    @arr = $values[$i];
                  }else{
                    @arr = split /,/, $values[$i];
                  }
                  foreach my $av (@arr){
                    $av =~ s/^\s+//g;
                    $av =~ s/\s+$//g;
                    push @{$record->{$solr[$i]}},lookup($solr[$i], "$av");
                  }
                }else{
                  #Ignore collection_latitude and collection_longitude if NA, since SOLR accepts only number
                  if (not ((lc $attribs[$i] eq "collection_latitude" or lc $attribs[$i] eq "collection_longitude") and $values[$i] eq "NA")) {
                    $new_values = lookup($solr[$i], "$values[$i]");
                    $record->{$solr[$i]} = $new_values;

                    #Add geographic_group
                    if (lc $attribs[$i] eq "collection_country"){
                      $record->{"geographic_group"} = $country_mapping{group}->{$values[$i]};
                    }
                  }
                }
              }
            }

            #Assign values for accession id
            if (lc $attribs[$i] eq "contributing_institution"){
              $contributing_institution = $values[$i];
            }elsif (lc $attribs[$i] eq "sample_identifier"){
              $sampleid = $values[$i];
            }elsif (lc $attribs[$i] eq "virus_identifier"){
              $virusid = $values[$i];
            }
          }else{
            $row_number = $values[$i];
          }
        }

        # Assign BVBRC Acession ID
        my $bvbrc_accession_id = $accession_id_map{$sampleid};
        if ($bvbrc_accession_id eq ''){
          my $bvbrc_accession_info = "";
          $seq++;
          my $seq_ext = sprintf("%010d", $seq);
          my $dcode = $ds_by_value{$contributing_institution};
          $bvbrc_accession_id = "BVBRC_".$dcode.$seq_ext;
          if( $type =~ /human/ or $type eq "animal" ){
            $bvbrc_accession_info .= "$row_number,$sampleid,$bvbrc_accession_id\n";
          }else{
            $bvbrc_accession_info .= "$row_number,$sampleid,$bvbrc_accession_id,$virusid\n";
          }
          $para->write_to_file($bvbrc_accession_file,$bvbrc_accession_info);

          $accession_id_map{$sampleid} = $bvbrc_accession_id;
        }
        $record->{"sample_accession"} = $bvbrc_accession_id;
        push(@values, ($bvbrc_accession_id));
        $csv->print($temp_file, \@values);
        print STDERR "Assigning accession id $bvbrc_accession_id to the sample $sampleid for row $row_number\n";

        push @records, $record;
      }
    }else{
      my ($cde, $str, $pos) = $csv->error_diag();
      die "ERROR: Line could not be parsed: $entry\nREASON: $cde, $str, $pos\n";
    }
  }
  close FH;
  close $temp_file;
  move($processed_file_tmp, $processed_file) or die "Move failed from $processed_file_tmp to $processed_file: $!";

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
   my ($attr,$data) = @_;

   my %common = ( 'Y'   => 'Yes',
     'N'   => 'No',
     'NON' => 'None',
     'U'   => 'Unknown',
     'NA'  => 'Not Applicable',
     'P'   => 'Positive',
     'N'   => 'Negative',
   );

   my $new_v; 
   if($attr eq "collection_country" ){
      if( $country{$data} ){
         $new_v = $country{$data}.",";
      }else{
         $new_v = $data.",";
      }
   }elsif($attr eq "host_common_name"){
     $data =~ s/(\S+)/\u\L$1/g;
     $new_v = $data.",";
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
   $new_v =~ s/\"//g;
   return $new_v;
}
