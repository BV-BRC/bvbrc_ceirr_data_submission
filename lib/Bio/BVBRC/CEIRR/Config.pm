package Bio::BVBRC::CEIRR::Config;
use strict;
use File::Basename;
use Module::Metadata;
use vars qw (@ISA @EXPORT_OK %configParams);

@ISA = qw ( Exporter ) ;
@EXPORT_OK = qw ( %configParams );

my $module_path = dirname(Module::Metadata->find_module_by_name(__PACKAGE__));

%configParams = (
  'LIB_DIR'  =>  $module_path,
);

my $para = \%configParams;

sub new
{
    my $class = shift;
    my $self = {
        _type => shift,
    };
    bless $self, $class;
    return $self;
}

sub setType {
    my ($self, $type) = @_;
    $self->{_type} = $type if defined($type);
    return $self->{_type};
}

sub getType {
    my($self) = @_;
    return $self->{_type};
}

sub getLookUp {
   my ($self) = @_;
   my $ldir = $para->{'LIB_DIR'};
   my $type  = $self->{_type};

   my %lookup;
   my $lfile = "${type}_lookup.txt";
   if( -e "$ldir/$lfile" ){
      open(FH, '<', "$ldir/$lfile") or die $!;
      while( my $line = <FH> ) {
         chomp $line;
         my @lines = split /,/, $line;
         $lookup{$lines[0]}->{$lines[1]} = "$lines[2]";
      }
      close FH;
   }
   return %lookup;
}

sub write_to_file {
   my ($self, $fname, $str, $opt) = @_;
  
   if ( $opt ){
      open(FH, '>', $fname) or die $!;
   }else{
      open(FH, '>>', $fname) or die $!;
   }
   print FH $str;
   close(FH);
}

sub getFileHeader {
  my ($self) = @_;
  my $libdir = $para->{'LIB_DIR'};
  my $type  = $self->{_type};
  my $fheader;
  my %hh;
  open(FH, '<', "$libdir/header.txt") or die $!;
  while( my $line = <FH> ) {   
    my @lines = split /=/, $line;
    if ( $type eq $lines[0] ) {
       $fheader = lc $lines[1];
    }
  }
  close FH;
  chomp $fheader;
  my @arr = split(',', $fheader);
  foreach (@arr) {
     $hh{$_} = $_;
  }
  return %hh; 
}

sub getCountryMapping {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};

   my %country;

   open(FH, '<', "$libdir/country_mapping.txt") or die $!;
   #my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      $country{'name'}->{$lines[0]} = $lines[0];
      $country{'group'}->{$lines[0]} = $lines[1];
   }
   close FH;
   return %country;
}

sub getCountryCode {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};
 
   my %country;

   open(FH, '<', "$libdir/country_code.txt") or die $!;
   my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      $country{'code'}->{lc $lines[0]} = $lines[0];          
      $country{'name'}->{lc $lines[0]} = $lines[1];
   }
   close FH;  
   return %country;
}

sub getUSAStates {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};

   my %states;

   open(FH, '<', "$libdir/usa_states.txt") or die $!;
   my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      $states{'code'}->{lc $lines[0]} = $lines[0];
      $states{'name'}->{lc $lines[1]} = $lines[1];
   }
   close FH;
   return %states;
}

sub getDataSourceCode {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};
   
   my %source;

   open(FH, '<', "$libdir/data_source.txt") or die $!;
   my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      my $val = uc $lines[1];
      $source{$val} = "$lines[3]";
   }
   close FH;
   return %source;
}

sub getDataVocab {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};
   my $type = $self->{_type};
   my %vocab;

   my $fname = "$type\.vocab.txt";
   open(FH, '<', "$libdir/$fname") or die $!; 
   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /=/, $line;
      my $field = lc $lines[0];
      my @fv = split /,/, $lines[1];
      
      foreach my $v (@fv) {
         $vocab{$field}->{lc $v} = $v;
      }
   }
   close FH;
   return %vocab;
}

sub getDataCode {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};
   
   my %code;

   open(FH, '<', "$libdir/data_code.txt") or die $!;
   my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      my $typ = uc $lines[0];
      $code{$typ}->{lc $lines[1]} = $lines[1];
   }
   close FH;
   return %code; 
}

sub getSpecies {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};

   my %spec;

   open(FH, '<', "$libdir/species_dic.txt") or die $!;
   my $header = <FH>;  # skip header

   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /\|/, $line;
      my $cn = lc $lines[0];
      my $hs = lc $lines[1];
      $spec{'host_name'}->{$cn} = $lines[0];
      $spec{'host_species'}->{$hs} = $lines[1];
   }
   close FH;
   return %spec;
}

sub getHeaderMap () {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};

   my %map;
   open(FH, '<', "$libdir/surv_header_map.txt") or die $!;
   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /=/, $line;
      my $old = lc $lines[0];
      my $new = lc $lines[1];
      $map{$old} = $new;
   }
   close FH;

   return %map;
}

sub getMapToSolr () {
   my ($self) = @_;
   my $libdir = $para->{'LIB_DIR'};
   my $type = $self->{_type};

   my %map;
   my $file;
   if ( $type eq "human" || $type eq "animal" || $type eq "human-sh" ) {
      $file = "surv_to_solr.txt";
   }else{
      $file = "serology_to_solr.txt";
   }
   open(FH, '<', "$libdir/$file") or die $!;
   my $header = <FH>;  # skip header
   while( my $line = <FH> ) {
      chomp $line;
      my @lines = split /,/, $line;
      my $ceir = lc $lines[1];
      my $solr = lc $lines[2];
      my $role = $lines[3];
      my $flag = $lines[4];
      if ( $ceir ne "" ){
         $map{$ceir}->{'solr'} = $solr;
         $map{$ceir}->{'role'} = $role if ( $role ne "" );
         $map{$ceir}->{'flag'} = $flag;
      }
   }
   close FH;

   return %map;
}

1;
