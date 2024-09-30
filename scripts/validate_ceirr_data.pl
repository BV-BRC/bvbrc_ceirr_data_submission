#!/usr/bin/perl
#======================================================================
# validate_data.pl
#======================================================================
use strict;
use Data::Dumper;
use Benchmark;
use Getopt::Std;
use POSIX qw/strftime/;
use File::Basename;
use Text::CSV;
use P3DataAPI;

use Bio::BVBRC::CEIRR::Config; 

sub usage{
  my $str = <<STR;
Usage: validate_data.pl -f file_path -t type[human,animal,serology,viral]
STR

  print STDERR "$str\n";
}

my $para;
my %options;
my $file = '';
my $type = '';
my $delimiter = ',';

getopts("f:t:i:",\%options);
if ( !defined $options{f} && !defined $options{t}){
  usage();
  exit(1);
} else{
  $file = $options{f};
  $type = $options{t};
}

$para = new Bio::BVBRC::CEIRR::Config; 
$para->setType($type);
my $jobDir = dirname($file);

my $valid_file   = "sample_valid.csv";
my $invalid_file = "sample_invalid.csv";
my $validation_report_file = "validation_report.csv";

my %header  = $para->getFileHeader($type);
my %headermap = $para->getHeaderMap();
my $flag = 0;

my $header_check = "P";
my $value_check  = "P";

#
# validate headers
#
print STDERR "Validating $file header...\n";

# define csv parser object
my $csv = Text::CSV->new({ binary => 1, sep_char => $delimiter }); 

# Read the header of data file for validation
open(FH, '<', "$file") or die "Failed to open $file: $!";
my $dummy_line = <FH>; # data template type
my $header_line = <FH>; # actual header
close FH;

chomp $header_line;
$header_line =~ s/["\s]+//gi; #remove double quotes and possible spaces from header

# Validate header columns
my @fh;
print STDERR "Header is $header_line\n"; 
if ($csv->parse($header_line)) {
  @fh = $csv->fields();
} else {
  my ($cde, $str, $pos) = $csv->error_diag ();
  print STDERR "$cde, $str, $pos\n";
  die "Cannot parse the header. Make sure it is a valid CSV object.\n";
}

foreach my $h (@fh){
  my $nh = lc $h;
  if( $header{$nh} ){
    $header_check="P";
  }else{
    if( $headermap{$nh} ){
      if( $header{$headermap{$nh}} ){
        $header_check="F";
        die "Header $nh is deprecated, replace it with new one: $header{$headermap{$nh}}";
      }else{
        $header_check="F";
        die "Invalid file header: $nh($header{$nh})";
      }
    }else{
      $header_check="F";
      print STDERR "Invalid file header column: $nh\n";
      die "Invalid file header: $nh($header{$nh})";
    }
  }
}

print STDERR "Header validation passed.\n";

# Validate field values 
my @data;
open(FH, '<', "$file") or die "Failed to open $file: $!";
while( my $line = <FH> ){
  chomp $line;

  $line =~ s/\r//g; #remove ^M characters 
  if ($csv->parse($line)) {
    my @lines = $csv->fields();
    s{^\s+|\s+$}{}g foreach @lines;
    push (@data, \@lines);
  } else {
    my ($cde, $str, $pos) = $csv->error_diag ();
    print STDERR "WARNING: Line could not be parsed: $line\nREASON: $cde, $str, $pos\n";
  }
}
close FH;

my %country    = $para->getCountryCode();
my %states     = $para->getUSAStates();
my %datasource = $para->getDataSourceCode();
my %vocab      = $para->getDataVocab();
my %datacode   = $para->getDataCode();
my %species    = $para->getSpecies();


my @headers  = @{$data[1]};
my $row_size = scalar @data;
my $col_size = scalar @headers;
my $sample_size = $row_size - 2;

print STDERR "Sample count: $sample_size, Column size: $col_size\n";

my %hash_data; 
for( my $j = 0; $j < $row_size; $j++ ){
   for( my $i = 0; $i < $col_size; $i++ ){
      $hash_data{lc $headers[$i]}->{$j} = $data[$j][$i];
   }
}

# Define api object
my $api = P3DataAPI->new; 

# Write headers
$para->write_to_file($invalid_file,"$dummy_line","new");
$para->write_to_file($invalid_file,"$header_line\n");
$para->write_to_file($valid_file,"Row,$header_line\n", "new");
$para->write_to_file($validation_report_file, "Row,Field,Value,Message\n", "new");

my $valid_sample_count = 0;
my $invalid_sample_count = 0;
for( my $j = 2; $j < $row_size; $j++ ){  # Ignore first two rows (DataTemplate and header)
   my $row_number = $j - 1;
   my $invalid_values = "";
   my $data_row = "";

   for ( my $i = 0; $i < $col_size; $i++ ) {
      my $header = lc $headers[$i];
      my $dval = "$data[$j][$i]";
      #print STDERR "ROW: $row_number, HEADER: $header, VALUE: $dval\n";

      my $check;
      if ( $header =~ m/state_province/ ) {
         $check = validate_value($header,"$dval","$hash_data{'collection_country'}->{$j}");
      } elsif ( $header =~ m/influenza_test_result/ or $header =~ m/influenza_test_interpretation/ ) {
         my @ttype = split (',', $hash_data{'influenza_test_type'}->{$j});
         my $num_ttype = @ttype;
         my @tresult = split (',', $dval);
         my $num_tresult = @tresult;
         if ( $num_tresult == $num_ttype ){
            $check = validate_value($header,"$dval","$hash_data{'influenza_test_type'}->{$j}");
         } else {
            $check = "Invalid: Number of Influenza_Test_Type, Influenza_Test_Antigen, Influenza_Test_Result, Influenza_Test_Interpretation should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      }elsif ( $header =~ m/sars-cov-2_test_result/ or $header =~ m/sars-cov-2_test_interpretation/ ){
         my @ttype = split (',', $hash_data{'sars-cov-2_test_type'}->{$j});
         my $num_ttype = @ttype;
         my @tresult = split (',', $dval);
         my $num_tresult = @tresult;
         if ( $num_tresult == $num_ttype ){
            $check = validate_value($header,"$dval","$hash_data{'sars-cov-2_test_type'}->{$j}");
         }else{
            $check = "Invalid: Number of SARS-CoV-2_Test_Type, SARS-CoV-2_Test_Antigen, SARS-CoV-2_Test_Result, SARS-CoV-2_Test_Interpretation should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      }elsif ( $header =~ m/definition_of_positive_sample/ ){
        my @ttype = split (',', $hash_data{'influenza_test_type'}->{$j});
        my $num_ttype = @ttype;

        my @dops = split (',', $dval);
        my $num_dops = @dops;
        if ( $num_dops == $num_ttype ){
          $check = validate_value($header,"$dval","$hash_data{'influenza_test_interpretation'}->{$j}");
        }else{
          $check = "Invalid: Number of Influenza_Test_Type and Definition_of_Positive_Sample should match. (Error_63_INVALID_NUMBER_ENTRIES)";
        }
      } elsif ( $header =~ m/collection_date/ or $header =~ m/sample_receipt_date/ ) {
         my $cdate = $hash_data{'collection_date'}->{$j};
         my $rdate = $hash_data{'sample_receipt_date'}->{$j};
         if ( !$cdate and !$rdate ){
            $check = "Invalid: At least one date must be entered for receipt or collection date. Fields cannot both be blank. (Error_118_DATE_MISSING)";
         } else {
            $check = validate_value($header,$dval);
         }
      } elsif ( $header =~ m/other_pathogen_test_result/ ) {
         my @otype = split (',', $hash_data{'other_pathogens_tested'}->{$j});
         my $num_otype = @otype;
         my @tresult = split (',', $dval);
         my $num_tresult = @tresult;
         if ( $num_tresult == $num_otype ) {
            $check = validate_value($header,"$dval","$hash_data{'other_pathogens_tested'}->{$j}");
         } else {
            $check = "Invalid: Number of Other_Pathogens_Tested and Other_Pathogen_Test_Result should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      } elsif ( $header =~ /duration_of_(poultry|wild_bird|swine|human|covid_human)_exposure/ ) {
         my ($dur, $exp_type) = $header =~ m/^(duration_of_)(.*)/;
         $check = validate_value($header,"$dval","$hash_data{$exp_type}->{$j}");
      } elsif( $header =~ /type_exposure/ ) {
         my $exp_value = "$hash_data{poultry_exposure}->{$j},$hash_data{wild_bird_exposure}->{$j},$hash_data{swine_exposure}->{$j}";
         $check = validate_value($header,"$dval","$exp_value");
      } elsif ( $header =~ /use_of_personal_protective_equipment/ ) {
         my $exp_value;
         if ( $type eq "human_sh" ){
            $exp_value = "$hash_data{covid_human_exposure}->{$j}";
         }else{
            $exp_value = "$hash_data{poultry_exposure}->{$j},$hash_data{wild_bird_exposure}->{$j},$hash_data{swine_exposure}->{$j},$hash_data{human_exposure}->{$j}";
         }
         $check = validate_value($header,"$dval","$exp_value");
      } elsif ( $header =~ /hospitalization_duration|intensive_care_unit|ventilation|oxygen_saturation|ecmo|dialysis/ ) {
         my $h_value = $hash_data{hospitalized}->{$j};
         $check = validate_value($header,"$dval","$h_value");
      } elsif ( $header =~ /initiation_of_vasoactive_treatment|vasoactive_treatment_dosage|duration_of_vasoactive_treatment/ ){
         my $at = $hash_data{vasoactive_treatment}->{$j};
         my @atval = split (',', $at);
         my $num_at = @atval;
         my @inval = split (',', $dval);
         my $num_in = @inval;
         if ( $num_at == $num_in ){
           $check = validate_value($header,"$dval","$at");
         }else{
           $check = "Invalid: Number of Vasoactive_Treatment; Initiation_of_Vasoactive_Treatment; Vasoactive_Treatment_Dosage and Duration_of_Vasoactive_Treatment should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      } elsif ( $header =~ /initiation_of_antiviral_treatment|treatment_dosage|duration_of_antiviral_treatment/ ) {
         my $at = $hash_data{antiviral_treatment}->{$j};
         my @atval = split (',', $at);
         my $num_at = @atval;
         my @inval = split (',', $dval);
         my $num_in = @inval;
         if ( $num_at == $num_in ) {
            $check = validate_value($header,"$dval","$at");
         } else {
            $check = "Invalid: Number of Antiviral_Treatment; Initiation_of_Antiviral_Treatment; Treatment_Dosage and Duration_of_Antiviral_Treatment should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      } elsif ( $header =~ /days_elapsed_to_influenza_vaccination|source_of_above_vaccine_information|vaccine_lot_number|vaccine_manufacturer|vaccine_dosage/ ){
         my $ivt = $hash_data{influenza_vaccination_type}->{$j};
         my @ivval = split (',', $ivt);
         my $num_ivt = @ivval;
         my @deval = split (',', $dval);
         my $num_de = @deval;
         if ( $num_ivt == $num_de ) {
            $check = validate_value($header,"$dval","$ivt");
         } else {
            $check = "Invalid: Number of Influenza_Vaccination_Type; Days_Elapsed_to_Influenza_Vaccination; Source_of_Above_Vaccine_Information; Vaccine_Lot_Number; Vaccine_Manufacturer; and Vaccine_Dosage should match. (Error_63_INVALID_NUMBER_ENTRIES)";
         }
      } elsif ( $header =~ m/types_of_allergies/ ) {
         $check = validate_value($header,"$dval",$hash_data{chronic_conditions}->{$j});
      } elsif ( $header =~ m/packs_per_day_for_how_many_years/ ) {
         $check = validate_value($header,"$dval",$hash_data{tobacco_use}->{$j}); 
      } elsif ( $header =~ m/host_habitat/ ) {
         $check = validate_value($header,"$dval",$hash_data{host_natural_state}->{$j});
      } elsif ( $header =~ m/influenza_subtype/ ) {
         my $itr = $hash_data{'influenza_test_interpretation'}->{$j};
         $check = validate_value($header,"$dval","$itr");
      } else {
         $check = validate_value($header,$dval);
      }

      #$data_row .= '"'.$dval.'",';
      $data_row .= '"'.$check.'",';

      if ( $check =~ m/(Invalid:)/ ) {
         $invalid_values .= "$row_number,$header,\"$dval\",$check\n";  
         #$invalid_values .= "$header: ".$dval." - $check\n"; 
      }
   }
   
   chop $data_row;

   if ( $invalid_values =~ m/(Invalid)/ ) {
      #print STDERR "DATA_ROW: $data_row\n";
      $para->write_to_file($invalid_file,"$data_row\n");
      $para->write_to_file($validation_report_file, "$invalid_values");
      $invalid_sample_count++;
   } else {
      $para->write_to_file($valid_file,"$row_number,$data_row\n");
      $valid_sample_count++;
   }
}   

print STDERR "There are $valid_sample_count valid samples in the submission.\n";
print STDERR "There are $invalid_sample_count invalid samples in the submission.\n";

print $valid_sample_count;

1;

sub validate_value{
   my ($header,$dval,$parse_v) = @_;

   #$dval =~ s/^\s+$//;
   $dval =~ s/^\s+//;
   $dval =~ s/\s+$//;
   $parse_v =~ s/^\s+//;
   $parse_v =~ s/\s+$//;    
   my $lcval = lc $dval;
   my $val_len = length($dval);
   my $check = $dval;

   ## study_identifier
   ## A unique Study Identifier generated by the iDPCC by combining the Center-generated Study Code and a random 5-digit number
   ## only check data format 
   ## Maximum allowed length: 50 characters
   ## will not check study_identifier unique
   if ( $header =~ m/study_identifier/ ){
      if ( $val_len > 50 ){
         #SERVICE#print STDERR "  invalid study_identifier:  $dval (maximum length 50 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $dval =~ /^(\w+)\_([0-9]{5})/i ){
            #SERVICE#print STDERR " passed\n";
         }else{
            #SERVICE#print STDERR "  invalid study_identifier: $dval(incorrect format)\n";
            #$check = "Invalid: Study_Identifier should be a valid Study Identifier. (Error_161_INVALID_STUDY_CODE)";
         }
      }
   }
   ## Contributing_Institution
   ## an Institution Code value registered with the iDPCC
   if ( $header =~ m/institution/ ){
      if ( !$datasource{uc $dval} ){
         #SERVICE#print STDERR "  invalid institution code: $dval\n";
         $check = "Invalid: Value must be an Institution Code value registered with the iDPCC (Error_1_INVALID_VALUE)";
      } else {
         $check = uc $dval;
      }
   }
   ## Embargo_End_Date, Collection_Date, Sample_Receipt_Date
   ## D-Mon-YYYY
   ## DD-Mon-YYYY
   ## D-Mon-YY
   ## DD-Mon-YY
   ## Mon-YYYY
   ## Mon-YY
   ## YYY
   ## YY
   ## Month is NOT case-sensitive
   ## Date cannot be future date
   ## NA for Embargo_End_Date
   ## U for Collection_Date,Sample_Receipt_Date 
   foreach my $field ("embargo_end_date","collection_date","sample_receipt_date"){ 
      if ( $header =~ m/$field/ ){
         if ( $vocab{$field}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{$field}->{lc $dval};
         }else{
            if ( $dval =~ /^([0-9]{1,2})\-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\-([0-9]{2,4})/i ){
               #SERVICE#print STDERR "  passed\n";
            }elsif( $dval =~ /^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\-([0-9]{2,4})/i ){
               #SERVICE#print STDERR "  passed\n";
            }elsif( $dval =~ /^([0-9]{2,4})$/ ){
               #SERVICE#print STDERR "  passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry:  $dval\n";
               $check = "Invalid: Value has invalid date format.";
            }
         }
      }
   }
   ## Collection_POI
   ## Maximum allowed length: 100 characters
   if ( $header =~ m/collection_poi/i ){
      if ( $val_len > 100 ){
         #SERVICE#print STDERR "  invalid data entry:  $dval (maximum length: 100 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $dval =~ m/\w+\,?\w+/){
            #SERVICE#print STDERR " passed\n";
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Invalid data entry. (Error_1_INVALID_VALUE)";
         }
      }
   }   
   ## Definition_of_Positive_Sample
   ## Maximum length: 100 characters
   ## If multiple comma-separated tests are listed under Influenza_Test_Type, Definition_of_Positive_Sample must list the same number of comma-separated test results.
   ## Enter NA if the Definition_of_Positive_Sample is unavailable or Influenza_Test_Interpretation is negative N or unknown U.
   if ( $header =~ m/definition_of_positive_sample/ ){
      my @dps = split (',', $dval);
      my @iti = split (',', $parse_v);
      $check = "";
      my $i = 0;
      foreach my $rt (@dps){
         $rt =~ s/^\s+//g;
         $rt =~ s/\s+$//g;
         $val_len = length($rt);
         my $itiv = $iti[$i];
         $itiv =~ s/^\s+//g;
         $itiv =~ s/\s+$//g;
         if ( $val_len > 100 ){
            #SERVICE#print STDERR "  invalid data entry:  $dval (maximum length: 100 characters)\n";
            $check .= "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{'definition_of_positive_sample'}->{lc $rt} ){
               #SERVICE#print STDERR " passed\n";
               $check .= "$vocab{'definition_of_positive_sample'}->{lc $rt},";
            }else{
               if ( $type eq "serology" ){
                  if ( $rt =~ m/\w+/){
                     #SERVICE#print STDERR " passed\n";
                     $check .= "$rt,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $rt\n";
                     $check .= "Invalid: Invalid data entry. (Error_1_INVALID_VALUE),";
                  }
               }else{
                  if ($itiv eq "N" || $itiv eq "U" ){
                     #SERVICE#print STDERR "  invalid data entry: $rt (test type entry: NA)\n";
                     $check .= "Invalid: Valume must be NA instead of N or U. (Error_1_INVALID_VALUE),";
                  }else{
                     if ( $rt =~ m/\w+/){
                        #SERVICE#print STDERR " passed\n";
                        $check .= "$rt,";
                     }else{
                        #SERVICE#print STDERR "  invalid data entry: $rt\n";
                        $check .= "Invalid: Invalid data entry. (Error_1_INVALID_VALUE),";
                     }
                  }
               }
            }
         }
         $i++;
      }
      chop $check;
   }
   ## Sample_Identifier
   ## Value must be unique and not match any previously submitted samples 
   ## Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _
   ## Maximum allowed length: 50 characters
   if ( $header =~ m/sample_identifier/ ){ 
      my @res = $api->query('surveillance', ['eq', 'sample_identifier', $dval], ['select', 'sample_identifier']);
      my $res_size = scalar @res;
      if ( $res_size > 0 ){
         if( $type eq "human" or $type eq "animal" ){
            #SERVICE#print STDERR "  invalid sample_identifier: existing $dval\n";
            $check = "Invalid: Sample_Identifier should be unique across all iDPCC data. (Error_66_NON_UNIQ_SAMPL_ID)";
         }else{
            #SERVICE#print STDERR " passed\n";
         }
      }else{
         if ( $type eq "viral" and $type eq "serology" ){
            #SERVICE#print STDERR "  invalid sample_identifier: not matching $dval\n";
            $check = "Invalid: The Sample_Identifier initially assigned to the surveillance sample must be provided. (Error_1_INVALID_VALUE)";
         }else{
            if ( $val_len > 50 ){
               #SERVICE#print STDERR "  invalid data entry:  $dval (maximum length: 50 characters)\n";
               $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
            }else{
               if ( $dval =~ /[a-zA-Z0-9\-\_]/ and $dval !~ /[\'\!\@\#\$\%\^\&\*\(\)\=\+\>\<\?\"\/\\\:\.\,\;]/ ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid sample identifier:  $dval, invalid characters\n";
                  $check = "Invalid: Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _ (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }
   ## Virus_Identifier
   ## Value must be unique and not match any previously submitted one, and its Sample_Identifier must match one of either human or animal 
   ## Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _
   ## Maximum allowed length: 50 characters
   if ( $header =~ m/virus_identifier/ ){
      if ( $vocab{virus_identifier}->{lc $dval} ){
         #SERVICE#print STDERR " passed\n";
         $check =  $vocab{virus_identifier}->{lc $dval};
      }else{
         my @res = $api->query('serology', ['eq', 'virus_identifier', $dval], ['select', 'virus_identifier']);
         my $res_size = scalar @res;
         if ( $res_size > 0 ){
            #SERVICE#print STDERR "  invalid virus_identifier: existing $dval\n";
            $check = "Invalid: Virus_Identifier should be unique. (Error_66_NON_UNIQ_SAMPL_ID)";
         }else{
            if ( $val_len > 50 ){
               #SERVICE#print STDERR "  invalid data entry:  $dval (maximum length: 50 characters)\n";
               $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
            }else{
               if (  $dval =~ /[a-zA-Z0-9\-\_]/ and $dval !~ /[\'\!\@\#\$\%\^\&\*\(\)\=\+\>\<\?\"\/\\\:\.\,\;]/ ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid virus identifier:  $dval, invalid characters\n";
                  $check = "Invalid: Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _ (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }

   ## Sample_Material,Education,Host_Capture_status,Host_Natural_State,Sample_Transport_Medium
   ## Maximum allowed length: 30 characters
   ## The entry must be one and only one member of the Value List, case-sensitive and all-caps. NOTE: User can enter other value by prefixing 'OTH-'
   foreach my $field ("sample_material","education","host_capture_status","host_natural_state","sample_transport_medium"){
      if ( $header =~ m/$field/ ){
         if ( $val_len > 30 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 30 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 30 characters (Error_75_INVALID_FIELD_LENGTH_OTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ m/^OTH-[\s+]?\w+/i ){
                  #SERVICE#print STDERR " passed\n";
                  $check = $dval;
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Field value should be one of valid values as in list. NOTE: User can enter other value by prefixing 'OTH-' (Error_1_INVALID_VALUE)";
               }
            }
         } 
      }
   }
   ## Subject_Identifier
   ## Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _
   ## Maximum allowed length: 50 characters
   if ( $header =~ m/subject_identifier/ ){
      if ( $val_len > 50 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 30 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $dval =~ /[a-zA-Z0-9\-\_]/ and $dval !~ /[\'\!\@\#\$\%\^\&\*\(\)\=\+\>\<\?\"\/\\\:\.\,\;]/ ){
            #SERVICE#print STDERR " passed\n";
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _ (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Subject_Age
   ## Enter Subject_Age as a whole number in years if three or older. 
   ## Values equal to or less than 3 can be an integer or number.
   ## If Subject_Age is less than three, submit in the format: 1 month/12=0.083 years or 35 months/12=2.91 years. 
   ## All ages greater than 90 must be entered as >90 to de-identify individual health information.
   ## Any age greater than three can be prefixed with > to de-identify individual health information.
   ## Maximum length: 17 characters
   ##
   ## Human: Any age greater than three can be prefixed with > to de-identify individual health information.
   ## Human-SH: For ages between 3 and 90, prefix with < if there is a need to mask the individual health information.
   if ( $header =~ m/subject_age/ ){
      if ( $val_len > 17 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 17 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'subject_age'}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{'subject_age'}->{lc $dval};
         }else{
            my $val = $dval;
            $val =~ s/^>//g;
            $val =~ s/^<//g;
            #if ( $val =~ /[a-zA-z\*\!\@\#\$\%\^\&\(\)=\+\<\?\"\/\:\']/g ){ # checking non-numeric characters
            if ( $val =~ /^\d+(\.\d+)?$/ ){    # check numeric
               if ( $val >= 3 ){
                  if ( $val =~ /(\.\d+)$/ ){
                     #SERVICE#print STDERR "  invalid subject_age entry: $dval (must be a whole number)\n";
                     $check = "Invalid: Value must be a whole number (Error_90_INVALID_AGE)";
                  }else{
                     if ( $val > 90 ){
                        if ( $dval =~ m/^>/ ){
                           #SERVICE#print STDERR " passed\n";
                        }else{
                           #SERVICE#print STDERR "  invalid subject_age entry: $dval (must be prefixed with > to de-identify)\n";
                           $check = "Invalid: All values greater than 90 must be entered as >90 (Error_90_INVALID_AGE)";
                        }
                     }else{
                        #SERVICE#print STDERR " passed\n";
                     }
                  } 
               }else{
                  #SERVICE#print STDERR " passed\n";
               }
            }else{
               #SERVICE#print STDERR "  invalid subject_age entry: $dval\n";
               $check = "Invalid: Value must not have non-numeric characters (Error_90_INVALID_AGE)";
            }
         }
      }
   }
   ## Collector_Name
   ## Maximum allowed length: 100 characters
   if ( $header =~ m/collector_name/ ){
      if ( $val_len > 100 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 100 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'collector_name'}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{'collector_name'}->{lc $dval};
         }else{
            if ( $dval =~ /[\*\!\@\#\$\%\^\&\(\)=\+\<\?\"\/\:\']/g ){ # checking special characters
               #SERVICE#print STDERR "  invalid subject_age entry: $dval (must not contain a special character)\n";
               $check = "Invalid: Value must not contain special characters. (Error_1_INVALID_VALUE)";
            }else{
               #SERVICE#print STDERR " passed\n";
            }
         }
      }
   }
   ## Days_Elapsed_to_Sample_Collection, Days_Elapsed_to_Disease_Status
   ## The value should be number, Not Collected, Not Provided or Restricted Access.
   ## 0 if the sample was collected on the same day of the subject's enrollment.
   ## a negative number using a hyphen if the sample was collected prior to the subject's enrollment in the study.
   foreach my $field ("days_elapsed_to_sample_collection","days_elapsed_to_disease_status"){
      if ( $header =~ m/$field/ ){
         if ( $val_len > 17 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 17 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ /^[-]?\d+/ ) {
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval (must not contain any characters except -)\n";
                  $check = "Invalid: Value should be number, Not Collected, Not Provided or Restricted Access. (Error_143_INVALID_NUMBER_INSDC)";
               }
            }
         }
      }
   }
   ## Collection_Season
   ## The value should be YYYY,YYYY-YYYY, Not Collected, Not Provided or Restricted Access.
   if ( $header =~ m/collection_season/ ){
      if ( $vocab{'collection_season'}->{lc $dval} ){
         #SERVICE#print STDERR " passed\n";
         $check = $vocab{'collection_season'}->{lc $dval};
      }else{
         if ( $dval =~ /[a-zA-Z\*\!\@\#\$\%\^\&\(\)=\+\<\?\"\/\:\']/g ){ # checking any characters except -
            #SERVICE#print STDERR "  invalid data entry: $dval (must not contain any characters except -)\n";
            $check = "Invalid: Value should be number, Not Collected, Not Provided or Restricted Access. (Error_143_INVALID_NUMBER_INSDC)";
         }else{
            if ( $dval =~ m/^\d{4}$/ ){ 
               #SERVICE#print STDERR " passed\n";
            }elsif ( $dval =~ m/^\d{4}\-\d{4}$/){
               my ($v1,$v2) = split ('-', $dval);
               if ($v1 < $v2 ){
                  #SERVICE#print STDERR " passed\n"; 
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval (must be sequential)\n";
                  $check = "Invalid: Year ranges must be sequential (Error_145_SEQUENTIAL_YEAR_RANGE)";
               }
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval (must be 4 digits for year)\n";
               $check = "Invalid: Value must be 4 digits for year (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Collection_Country
   ## Field value should be a valid ISO three-letter country code or U
   if ( $header =~ m/country/ ){
      if ( $vocab{'collection_country'}->{lc $dval} ){
         #SERVICE#print STDERR " passed\n";
         $check = $vocab{'collection_country'}->{lc $dval};
      }else{
         if ( $country{'code'}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $country{'code'}->{lc $dval};
         }else{
            #SERVICE#print STDERR "  invalid country name: $dval (must be a three letter code)\n";
            $check = "Invalid: Field value should be a valid ISO three-letter country code or U (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Collection_State_Province
   ## Validate US state names only 
   if ( $header =~ m/state_province/ ){
      if ( $val_len > 50 ){
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'collection_state_province'}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{'collection_state_province'}->{lc $dval};
         }else{
            if ( uc $parse_v eq "USA" ){
               if ( $states{'name'}->{lc $dval} ){
                  #SERVICE#print STDERR " passed\n";
                  $check = $states{'name'}->{lc $dval};
               }else{
                  #SERVICE#print STDERR "  invalid USA state name: $dval\n";
                  $check = "Invalid: Field value should be a correct USA state name (Error_1_INVALID_VALUE)";
               }
            }else{
               if ( $states{'name'}->{lc $dval} ){
                  #SERVICE#print STDERR "  invalid state/province name for $parse_v\n";
                  $check = "Invalid: Field value should be a correct state/province name (Error_1_INVALID_VALUE)";
               }else{
                  #SERVICE#print STDERR " passed\n";
               }
            }
         }
      }
   }
   ## Influenza_Test_Type, SARS-CoV-2_Test_Type
   ## iDPCC Data Dictionary
   ## The entry must be one or more comma-separated members of the Value List.
   foreach my $field ("influenza_test_type","sars-cov-2_test_type"){
      if ( $header =~ m/$field/ ){
         my @ttype = split (',', $dval);
         $check ="";
         foreach my $t (@ttype){
            $t =~ s/^\s+//g;
            $t =~ s/\s+$//g;
            if ( $vocab{$field}->{lc $t} ){
               #SERVICE#print STDERR " passed\n";
               $check .= "$vocab{$field}->{lc $t},";
            }else{
               if ( $datacode{'TEST_TYPE'}->{lc $t} ){
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$datacode{'TEST_TYPE'}->{lc $t},";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $t\n";
                  $check .= "Invalid: Values must reference iDPCC Data Dictionary. (Error_13_REFER_TO_DATA_DICTIONARY),";
               }
            }
         }
         chop $check;
      }
   }
   ## Influenza_Test_Result, SARS-CoV-2_Test_Result
   ## The numerical result(s) of Influenza_Test_Type
   ## If multiple comma-separated tests are listed under Influenza_Test_Type, Influenza_Test_Result must list the same number of comma-separated test results
   ## If multiple comma-separated tests are listed under SARS-CoV-2_Test_Type, SARS-CoV-2_Test_Result must list the same number of comma-separated test results
   ## Enter NA if the value under Influenza_Test_Type/SARS-CoV-2_Test_Type is NA.
   foreach my $field ("influenza_test_result","sars-cov-2_test_result"){
      if ( $header =~ m/$field/ ){
         my @tresult = split (',', $dval);
         my @ttype   = split (',', $parse_v);
         $check ="";
         my $i = 0;
         foreach my $re (@tresult){
            $re =~ s/^\s+//;
            $re =~ s/\s+$//;
            $ttype[$i] =~ s/^\s+//;
            $ttype[$i] =~ s/\s+$//;
            if ( $vocab{$field}->{lc $re} ){
               if (uc $ttype[$i] eq "NA"){
                  if ( uc $re eq "NA"){
                     #SERVICE#print STDERR " passed\n";
                     $check .= "NA,";
                  }elsif ( uc $re eq "U"){
                     #SERVICE#print STDERR " passed\n";
                     $check .= "U,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $rt (test type entry: NA)\n";
                     $check .= "Invalid: $field value must be NA or U if test type is NA (Error_1_INVALID_VALUE),";
                  }
               }else{
                  #SERVICE#print STDERR " passed\n";
                  $check .= "$vocab{$field}->{lc $re},";
               }
            }else{
               #if ( $re =~ /^[0-9]\d*(\.\d+)?$/ || $re =~ /^<|>|(<=)|(>=)\d+\.?\d+$/ ){
               if ( $re =~ /^[<|>|(<=)|(>=)]?\d+(\.\d+)?$/ ){
                  #SERVICE#print STDERR " passed\n";
                  $check .= "$re,";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $re (must be numeric)\n";
                  $check .= "Invalid: Field $field value must be number, U, or NA (Error_86_INVALID_NUM_U_NA),";
               }
            }
         }
         chop $check;
      }
   }
   ## Influenza_Test_Interpretation, SARS-CoV-2_Test_Interpretation
   ## The entry must be one or more comma-separated members of the Value List, all-caps
   ## Enter NA if the value under Influenza_Test_Type/SARS-CoV-2_Test_Interpretation is NA.
   foreach my $field ("influenza_test_interpretation","sars-cov-2_test_interpretation"){
      if ( $header =~ m/$field/ ){
         my @tinter = split (',', $dval);
         my @ttype  = split (',', $parse_v);
         $check ="";
         my $i = 0;
         foreach my $rt (@tinter){
            $rt =~ s/\s+//g;
            if ( $vocab{$field}->{lc $rt} ){
               if (uc $ttype[$i] eq "NA"){
                  if ( uc $rt eq "NA"){
                     #SERVICE#print STDERR " passed\n";
                     $check .= "NA,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $rt (test type entry: NA)\n";
                     $check .= "Invalid: Field value should be NA if test type is NA. (Error_1_INVALID_VALUE),";
                  }
               }else{
                  #SERVICE#print STDERR " passed\n";
                  $check .= "$vocab{$field}->{lc $rt},";
               }
            }else{
               #SERVICE#print STDERR "  invalid data entry: $rt\n";
               $check .= "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE),";
            }
            $i++;
         }
         chop $check;
      }
   }
   ## Other_Pathogens_Tested
   ## A list for organisms tested for other than influenza virus
   ## Maximum allowed length: 200 characters
   if ( $header =~ m/other_pathogens_tested/ ){
      $check ="";
      my @otype = split (',', $dval);
      my $num_otype = @otype;
      foreach my $ot (@otype){
         $ot =~ s/^\s+//g;
         $ot =~ s/\s+$//g;  
         $val_len = length($ot);
         if ( $val_len > 200 ){
            #SERVICE#print STDERR "  invalid data entry: $ot (maximum length: 200 characters)\n";
            $check .= "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH),";
         }else{
            if ( $vocab{'other_pathogens_tested'}->{lc $ot} ){
               #SERVICE#print STDERR "  passed\n";
               $check .= "$vocab{'other_pathogens_tested'}->{lc $ot},";
            }else{
               if ( $ot =~ /\w/ ){
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$ot,";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $ot\n";
                  $check .= "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),";
               }
            }
         }
      }
      chop $check;
   }
   ## Other_Pathogen_Test_Result
   ## 
   if ( $header =~ m/other_pathogen_test_result/ ){
      $check ="";
      my $i = 0;
      my @oresult = split (',', $dval);
      my @ttype  = split (',', $parse_v);
      foreach my $or (@oresult){
         $or =~ s/^\s+//g;
         $or =~ s/\s+$//g;  
         if ( $vocab{'other_pathogen_test_result'}->{lc $or} ){
            if (uc $ttype[$i] eq "NA"){
               if ( uc $or eq "NA"){
                   #SERVICE#print STDERR " passed\n";
                  $check .= "NA,";
               }else{
                   #SERVICE#print STDERR "  invalid data entry: $or (test type entry: NA)\n";
                   $check .= "Invalid: Field value should be NA if test type is NA. (Error_1_INVALID_VALUE),";
               }
            }else{
               #SERVICE#print STDERR " passed\n";
               $check .= "$or,";
            }
         }else{
            #SERVICE#print STDERR "  invalid data entry: $or\n";
            $check .= "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE),"; 
         }
         $i++;
      }
      chop $check;
   }
   ## Poultry_Exposure, Wild_Bird_Exposure, Swine_Exposure, Human_Exposure, COVID_Human_Exposure
   ## The entry must be one and only one member of the Value List. Values are case-sensitive and must be entered in all-caps.
   foreach my $tp ("poultry","wild_bird","swine","human","covid_human"){
      my $expo = $tp."_exposure";
      if ( $header =~ m/^$expo/ ){
         if ( $vocab{$expo}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{$expo}->{lc $dval};
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Duration_of_Poultry_Exposure, Duration_of_Wild_Bird_Exposure, Duration_of_Swine_Exposure, Duration_of_Human_Exposure
   ## If Y is selected entered for Poultry_Exposure, Duration_of_Poultry_Exposure must be a number or Not Collected, Not Provided, or Restricted Access. 
   ## If N is entered for Poultry_Exposure, Duration_of_Poultry_Exposure must be NA.
   ## If Not Collected is entered for Poultry_Exposure, Duration_of_Poultry_Exposure must be Not Collected. 
   ## If Not Provided is entered for Poultry_Exposure, Duration_of_Poultry_Exposure must be Not Provided. 
   ## If Restricted Access is entered for Poultry_Exposure, Duration_of_Poultry_Exposure must be Restricted Access. 
   foreach my $tp ("poultry","wild_bird","swine","human"){
      my $dura = "duration_of_".$tp."_exposure";
      if ( $header =~ m/^$dura/ ){
         if ( $val_len > 17 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 17 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$dura}->{lc $dval} ){
               if ( uc $parse_v eq "Y"){
                  if ( uc $dval eq "NA" ){
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Field value should NOT be NA if exposure is Y. (Error_1_INVALID_VALUE)";
                  }else{
                     $check =  $vocab{$dura}->{lc $dval};
                  }
               }elsif ( uc $parse_v eq "N" ){
                  if ( uc $dval eq "NA" ){
                     #SERVICE#print STDERR " passed\n";
                     $check =  $vocab{$dura}->{lc $dval};
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Field value should be NA if exposure is N. (Error_1_INVALID_VALUE)";
                  }
               }else{
                  if ( uc "$parse_v" eq uc "$dval" ){
                     #SERVICE#print STDERR " passed\n";
                     $check = $vocab{$dura}->{lc $dval};
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: If exposure is Y; duration must be a number or Not Collected; Not Provided; Restricted Access. If exposure is N; duration must be NA. (Error_146_INVALID_EXP_DURATION)";
                  }
               } 
            }else{
               if ( $dval =~ /^\d+(\.\d+)?$/ and uc $parse_v eq "Y" ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Value must be a number. (Error_144_INVALID_NUMBER_NA_INSDC)";
               }
            }
         }
      }
   } 
   ## Type_Exposure
   ## The entry must be three comma-separated members of the Value List. NOTE: User can enter other value by prefixing 'OTH-'
   ## If Y is entered for Poultry_Exposure, Wild_Bird_Exposure, or Swine_Exposure, Type_Exposure cannot be NA.
   ## If N is entered for Poultry_Exposure, Wild_Bird_Exposure, or Swine_Exposure, Type_Exposure must be NA for the specific exposure within the comma-separated list.
   ## If Not Collected is entered for Poultry_Exposure, Wild_Bird_Exposure, or Swine_Exposure, Type_Exposure must be Not Collected for the specific exposure within the comma-separated list.
   ## If Not Provided is entered for Poultry_Exposure, Wild_Bird_Exposure, or Swine_Exposure, Type_Exposure must be Not Provided for the specific exposure within the comma-separated list.
   ## If Restricted Access is entered for Poultry_Exposure, Wild_Bird_Exposure, or Swine_Exposure, Type_Exposure must be Restricted Access for the specific exposure within the comma-separated list.
   if ( $header =~ /type_exposure/ ){
      my @tval = split (',', $dval);
      my @expv = split (',', $parse_v);
      my $num_val = @tval;
      $check = "";
      if ( $num_val == 3 ){
         #SERVICE#print STDERR "\n";
         my @type =  ("poultry","wild_bird","swine");
         my $p = 0;
         foreach my $t (@tval){
            $t =~ s/^\s+//g;
            $t =~ s/\s+$//g;
            $expv[$p] =~ s/^\s+//g;
            $expv[$p] =~ s/\s+$//g;
            my $expo = $type[$p]."_exposure";
            if ( $vocab{'type_exposure'}->{lc $t} ){
               if ( uc $expv[$p] eq "Y" ){
                  if ( uc $t eq "NA" ){
                     #SERVICE#print STDERR "  invalid data entry: $t\n";
                     $check .= "Invalid: Field value must not be NA if exposure is Y (Error_147_TYPE_EXPOSURE_NA),";
                  }else{
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "$vocab{'type_exposure'}->{lc $t},";
                  }
               }elsif ( uc $expv[$p] eq "N" ){
                  if ( uc $t eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "NA,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $t\n";
                     $check .= "Invalid: Field value must be NA if exposure is N (Error_147_TYPE_EXPOSURE_NA),";
                  }
               }else{
                  if ( uc "$expv[$p]" eq uc "$t" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "$vocab{'type_exposure'}->{lc $t},";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $t, for $p\n";
                     $check .= "Invalid: The comma-separated values should match with the values in Poultry_Exposure; Wild_Bird_Exposure; or Swine_Exposure. (Error_150_INVALID_TYPE_EXPOSURE),";
                  }
               }
            }else{
               if ( $t =~ m/^OTH-[\s+]?\w+/i ){
                  if ( $val_len > 60 ){
                     #SERVICE#print STDERR "  invalid data entry: $t\n";
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 60 characters (Error_70_INVALID_FIELD_LENGTH),";
                  }else{
                     #SERVICE#print STDERR "  passed\n";
                     $t =~ s/oth/OTH/i;
                     $check .= "$t,";
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $t\n";
                  $check .= "Invalid: Field value should be of valid values as in list. (Error_1_INVALID_VALUE),";
               }
            }
            $p++;
         }
      }else{
         #SERVICE#print STDERR "  invalid data entry: $dval (must be three members)\n";
         $check .= "Invalid: Field value should be three valid values of the list. (E_151_INVALID_TYPE_EXPOSURE_COUNT),";
      }
      chop $check;
   }
   ## Use_of_Personal_Protective_Equipment
   ## Maximum allowed length: 50 characters
   ## If N is entered for Poultry_Exposure, Wild_Bird_Exposure, Swine_Exposure, and Human_Exposure, Use_of_Personal_Protective_Equipment must be NA.
   if ( $header =~ /use_of_personal_protective_equipment/ ){
      my @expv = split (',', $parse_v);
      if ( $val_len > 50 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 50 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'use_of_personal_protective_equipment'}->{lc $dval} ){
            my $flag = 'NN';
            foreach my $v (@expv){
               $v =~ s/\s+//g;
               $flag = 'YY' if ( uc $v ne "N" );
            }
            if ( $flag eq "NN" ){
               if ( uc $dval eq "NA" ){
                  #SERVICE#print STDERR "  passed\n";
                  $check = $vocab{'use_of_personal_protective_equipment'}->{lc $dval};
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval (at least one exposure is not N)\n";
                  $check = "Invalid: If N is entered for Poultry_Exposure; Wild_Bird_Exposure; Swine_Exposure; and Human_Exposure; Use_of_Personal_Protective_Equipment must be NA. (Error_149_ONLY_NA_ALLOWED)";
               }
            }else{
               #SERVICE#print STDERR "  passed\n";
               $check = $vocab{'use_of_personal_protective_equipment'}->{lc $dval};
            }
         }else{
            if ( $dval =~ m/\w+/ ){
               #SERVICE#print STDERR "  passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Hospitalized
   ## Field value should be one of valid values as in list.
   ## Maximum length: 17 characters
   foreach my $field ("hospitalized"){
      if ( $header =~ m/^$field$/ ){
         if ( $vocab{$field}->{lc $dval} ){
            #SERVICE#print STDERR "  passed\n";
            $check = $vocab{$field}->{lc $dval};
         }else{
            if ( $type eq "human" ){
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }else{
               my ($dvala,$dvalb) = split ('/',$dval);
               if ( uc $dvala eq "N" ){
                  if ( uc $dvalb eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check = uc $dval;
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
                  }
               }elsif(uc $dvala eq "Y" ){
                  if ( $dvalb =~ /^\d+/ or uc $dvalb eq "U" ){
                     $check = uc $dval;
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }
   ## Hospitalization_Duration
   ## Length of time subject was hospitalized in days
   ## If N is entered for Hospitalized, Hospitalization_Duration must be NA.
   if ( $header =~ m/^hospitalization_duration/ ){
      if ( $vocab{'hospitalization_duration'}->{lc $dval} ){
         if ( uc $parse_v eq "N" ){
            if ( uc $dval eq "NA" ){
               #SERVICE#print STDERR "  passed\n";
               $check = $vocab{'hospitalization_duration'}->{lc $dval};
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: If N is selected entered for Hospitalized; Hospitalization_Duration must be NA. (Error_149_ONLY_NA_ALLOWED)";
            }
         }else{
            #SERVICE#print STDERR "  passed\n";
            $check = $vocab{'hospitalization_duration'}->{lc $dval};
         }
      }else{
         if ( $dval =~ /^\d+[\.\d+]?$/ and uc $parse_v eq "Y" ){
            #SERVICE#print STDERR "  passed\n";
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Value must be a number; NA; or an INSDC value. (Error_144_INVALID_NUMBER_NA_INSDC)";
         }
      }
   }
   ## Intensive_Care_Unit
   ## Field value should be one of valid values as in list. 
   ## If N is entered for Hospitalized, Intensive_Care_Unit must be NA.
   ## Maximum length: 17 characters
   foreach my $field ("intensive_care_unit"){
      if ( $header =~ m/^$field/ ){
         if ( $vocab{$field}->{lc $dval} ){
            if ( uc $parse_v eq "N" and $type eq "human" ){
               if ( uc $dval eq "NA" ){
                  #SERVICE#print STDERR "  passed\n";
                  $check = $vocab{$field}->{lc $dval};
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: If N is selected entered for Hospitalized, Intensive_Care_Unit must be NA. (Error_149_ONLY_NA_ALLOWED)";
               }
            }else{
               #SERVICE#print STDERR "  passed\n";
               $check = $vocab{$field}->{lc $dval};
            }
         }else{
            if ( $type eq "human" ){
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE)";
            }else{
               my ($dvala,$dvalb) = split ('/',$dval);
               if ( uc $dvala eq "N" ){
                  if ( uc $dvalb eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check = uc $dval;
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
                  }
               }elsif(uc $dvala eq "Y" ){
                  if ( $dvalb =~ /^\d+/ or uc $dvalb eq "U" ){
                     $check = uc $dval;
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }

   ## Ventilation
   ## Field value should be one of valid values as in list.
   ## Multiple ventilation types can be entered as comma-separated values.
   ## If N is entered for Hospitalized, Ventilation must be NA
   foreach my $field ("ventilation"){
      if ( $header =~ m/^$field$/ ){
         $check = "";
         my @vval = split (',', $dval);
         foreach my $v (@vval){
            $v =~ s/^\s+//g;
            $v =~ s/\s+$//g; 
            if ( $vocab{'ventilation'}->{lc $v} ){
               if ( uc $parse_v eq "N" ){
                  if ( uc $dval eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "NA,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $dval\n";
                     $check .= "Invalid: If N is selected entered for Hospitalized; Ventilation must be NA. (Error_149_ONLY_NA_ALLOWED),";
                  }
               }else{
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$vocab{'ventilation'}->{lc $v},";
               }
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check .= "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE),";
            }
         }
         chop $check;
      }
   }
   ## Oxygen_Saturation, ECMO, Dialysis
   ## Saturation levels reported as a percentage (%) and as whole number.
   ## If N is entered for Hospitalized, Oxygen_Saturation must be NA.
   foreach my $field ("oxygen_saturation","ecmo","dialysis"){
      if ( $header =~ m/^$field/ ){
         if ( $vocab{$field}->{lc $dval} ){
            if ( uc $parse_v eq "N" ){
               if ( uc $dval eq "NA" ){
                  #SERVICE#print STDERR "  passed\n";
                  $check = "NA";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: If N is entered for Hospitalized, value must be NA (Error_149_ONLY_NA_ALLOWED)";
               }
            }else{
               #SERVICE#print STDERR "  passed\n";
               $check = "$vocab{$field}->{lc $dval}";
            }
         }else{
            if ( $dval =~ /^\d+/ and $header eq "oxygen_saturation" ) {
               #SERVICE#print STDERR "  passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Field value should be one of valid values as in list. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Antiviral_Treatment, Vasoactive_Treatment, Other_Treatments
   ## Field value should be one of valid values as in list. NOTE: User can enter other value by prefixing 'OTH-'
   ## Maximum allowed length: 100 characters
   ## If multiple treatments were administered then this field must contain information for each administration.
   foreach my $field ("antiviral_treatment","vasoactive_treatment","other_treatments"){
      if ( $header =~ m/^$field/ ){
         $check = "";
         my @atval = split (',', $dval);
         foreach my $v (@atval){
            $v =~ s/^\s+//g; #trim if multiple values
            $v =~ s/\s+$//g;
            if ( $vocab{$field}->{lc $v} ){
               #SERVICE#print STDERR "  passed\n";
               $check .= "$vocab{$field}->{lc $v},";
            }else{
               if ( $v =~ m/^OTH-[\s+]?\w+/i ){
                  $val_len = length($v);
                  if ( $val_len > 100 ){
                     #SERVICE#print STDERR "  invalid data entry: $v(over 100 characters)\n";
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_75_INVALID_FIELD_LENGTH_OTH),";
                  }else{
                     #SERVICE#print STDERR "  passed\n";
                     $v =~ s/oth/OTH/i;
                     $check .= "$v,";
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $v\n";
                  $check .= "Invalid: Field value must be OTH-text (Error_1_INVALID_VALUE),";
               }
            }
         }
         chop $check;
      }
   }
   ## Initiation_of_Antiviral_Treatment, Treatment_Dosage, Duration_of_Antiviral_Treatment
   ## Initiation_of_Vasoactive_Treatment, Vasoactive_Treatment_Dosage, Duration_of_Vasoactive_Treatment
   ## Number of days after onset of clinical symptoms antiviral treatment was initiated in days
   ## If multiple treatments were administered then this field must contain information for each administration provided in Antiviral_Treatment.
   ## If NON is entered for Antiviral_Treatment, Initiation_of_Antiviral_Treatment must be NA.
   ## Number could be positive, negative, or 0.
   foreach my $field ("initiation_of_antiviral_treatment","treatment_dosage","duration_of_antiviral_treatment","initiation_of_vasoactive_treatment","vasoactive_treatment_dosage","duration_of_vasoactive_treatment"){
      if ( $header =~ m/^$field$/ ){
         $check = "";
         my @atval = split (',', $parse_v);
         my @inval = split (',', $dval);
                
         my $n = 0;
         foreach my $v (@inval){
            $v =~ s/^\s+//g;
            $v =~ s/\s+$//g; 
            if ( $vocab{$field}->{lc $v} ){
               my $at = $atval[$n];
               $at =~ s/^\s+//g;
               $at =~ s/\s+$//g;   
               if ( uc $at eq "NON" ){
                  if ( uc $v eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "NA,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $v\n";
                     $check .= "Invalid: If NON is selected entered for Antiviral_Treatment; value must be NA (Error_149_ONLY_NA_ALLOWED),";
                  }
               }else{
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$vocab{$field}->{lc $v},";
               }
            }else{
               if ( $field eq "treatment_dosage" ){
                  if ( $v =~ /^[(\d+]|[(\d+\.?\d+)]\s?(mg|ml)?/i ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "$v,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $v\n";
                     $check .= "Invalid: Value must be a number; NA; or an INSDC value. (Error_144_INVALID_NUMBER_NA_INSDC),";
                  }
               }else{
                  $v =~ s/\s+//g;
                  if ( $v =~ /^[+-]?\d+(\.\d+)?/ ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "$v,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $v\n";
                     $check .= "Invalid: Value must be a number; NA; or an INSDC value. (Error_144_INVALID_NUMBER_NA_INSDC),";
                  }
               }
            }
            $n++;
         }
         chop $check;
      }
   }
   ## Influenza_Vaccination_Type, Other_Vaccinations
   ## If multiple vaccinations were administered then this field must contain information for each vaccination.
   ## Maximum allowed length: 50 characters
   foreach my $field ("influenza_vaccination_type","other_vaccinations"){
      if ( $header =~ m/$field/ ){
         $check = "";
         my @ivval = split (',', $dval);
         foreach my $v (@ivval){
            $v =~ s/^\s+//g;
            $v =~ s/\s+$//g;
            $val_len = length($v);
            if ( $val_len > 50 ){
               #SERVICE#print STDERR "  invalid data entry: $dval(over 50 characters)\n";
               $check .= "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH),";
            }else{
               if ( $vocab{$field}->{lc $v} ){
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$vocab{$field}->{lc $v},";
               }else{
                  if ( $v =~ m/^\w+/ ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "$v,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $v\n";
                     $check .= "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),";
                  }
               }
            }
         }
         chop $check;
      }
   }
   ## Days_Elapsed_to_Influenza_Vaccination,Source_of_Above_Vaccine_Information,Vaccine_Lot_Number,Vaccine_Manufacturer,Vaccine_Dosage
   ## If NON is entered for Influenza_Vaccination_Type, Days_Elapsed_to_Influenza_Vaccination must be NA.
   foreach my $field ("days_elapsed_to_influenza_vaccination","source_of_above_vaccine_information","vaccine_lot_number","vaccine_manufacturer","vaccine_dosage"){
      if ( $header =~ m/^$field$/ ){
         my @ivval = split (',', $parse_v);
         my @deval = split (',', $dval);

         $check = "";
         my $n = 0;
         foreach my $v (@deval){
            $v =~ s/^\s+//g;
            $v =~ s/\s+$//g;
            $ivval[$n] =~ s/^\s+//g;
            $ivval[$n] =~ s/\s+$//g;
            $val_len = length ($v); 
            if ( $vocab{$field}->{lc $v} ){
               if ( uc $ivval[$n] eq "NON" ){
                  if ( uc $v eq "NA" ){
                     #SERVICE#print STDERR "  passed\n";
                     $check .= "NA,";
                  }else{
                     #SERVICE#print STDERR "  invalid data entry: $v\n";
                     $check .= "Invalid: If NON is entered for Influenza_Vaccination_Type, value must be NA (Error_149_ONLY_NA_ALLOWED),";
                  }
               }else{
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$vocab{$field}->{lc $v},";
               }
            }else{
               if ( $field eq "days_elapsed_to_influenza_vaccination" ){
                  if ( $val_len > 17 ){
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH),";
                  }else{
                     if ( $v =~ /^[-]?(\d+)\;[-]?(\d+)$/ ){
                        #SERVICE#print STDERR "  passed\n";
                        $check .= "$v,";
                     }else{
                        $check .= "Invalid: Value must be two numbers separated by a semicolon; NA; or an INSDC value. (Error_152_INVALID_DURATION_NA_INSDC),";
                     }
                  }
               }elsif ( $field eq "vaccine_lot_number" or $field eq "vaccine_manufacturer" ){
                  if ( $val_len > 50 ){
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH),";
                  }else{
                     if ( $v =~ /^\w+/ ){
                        #SERVICE#print STDERR "  passed\n";
                        $check .= "$v,";
                     }else{
                        #SERVICE#print STDERR "  invalid data entry: $v\n";
                        $check .= "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),";
                     }
                  }
               }elsif ($field eq "vaccine_dosage"){
                  if ( $val_len > 50 ){
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH),";
                  }else{ 
                     if ( $v =~ /^[(\d+\s?]|[(\d+\.?\d+)\s?][ml|mg]/i ){
                        #SERVICE#print STDERR "  passed\n";
                        $check .= "$v,";
                     }elsif ( $v =~ /^(\d+\.?\d*)[ml|mg]/ ){
                        substr $v, -2, 0, ' '; #add space before mL
                        $check .= "$v,"
                     }else{
                        #SERVICE#print STDERR "  invalid data entry: $v\n";
                        $check .= "Invalid: Measurement and its unit must be valid. (Error_91_INVALID_MEASUREMENT_W_UNIT),";
                     }
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $v\n";
                  $check .= "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),";
               }
            }
            $n++;
         }
         chop $check;
      }
   }
   ## Influenza_Like_Illness,Chills,Conjunctivitis,Cough,Diarrhea,Fever,Headache,Loss_of_Appetite,Malaise,Myalgia,Nausea,Runny_Nose,Shortness_of_Breath,Sore_Throat,Vomiting,Wheezing
   ## Maximum length: 17 characters
   ## Y,Number,Number first number could be negative
   ## Y,Number,U      first number could be negative
   ## Y,U,Number
   ## Y,U,U
   foreach my $field ("influenza_like_illness","chills","conjunctivitis","cough","diarrhea","fever","headache","loss_of_appetite","malaise","myalgia","nausea","runny_nose","shortness_of_breath","sore_throat","vomiting","wheezing","nasal_congestion","odynophagia","cyanosis","who_covid_disease_severity","abdominal_pain","anosmia"){
      if ( $header =~ /^$field$/ ){
         if ( $val_len > 17 ){
            #SERVICE#print STDERR "  invalid data entry: $dval(over 17 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR "  passed\n";
               $check = "$vocab{$field}->{lc $dval}";
            }else{
               $dval =~ s/\s+//g;
               if ( uc($dval) =~ /^Y,([-]?\d+,\d+)|([-]?\d+,U)|(U,\d+)|(U,U)/ ){
                  #SERVICE#print STDERR "  passed\n";
                  $check = uc $dval;
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Symptom data must be corrected. (Error_92_INVALID_SYMPTOM)";
               }
            }
         }
      }
   }
   ## Other_Symptoms
   ## Maximum allowed length: 100 characters
   ## If symptom is selected, two values must be entered following a comma.
   ## If multiple additional symptoms were recorded, separate the individual symptoms with a semicolon.
   ## If N, Not Collected, Not Provided, or Restricted Access is entered, no additional values are needed.
   if ( $header =~ m/other_symptoms/ ){
      $check = "";
      if ( $val_len > 100 ){
         #SERVICE#print STDERR "  invalid data entry: $dval(over 100 characters)\n";
         $check .= "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH),";
      }else{
         if ( $vocab{'other_symptoms'}->{lc $dval} ){
            #SERVICE#print STDERR "  passed\n";
            $check .= "$vocab{'other_symptoms'}->{lc $dval},";
         }else{
            my @oval = split (';', $dval);
            my $val_size = @oval;
            foreach my $v (@oval){
               $v =~ s/^\s+//g; 
               $v =~ s/\s+$//g;   
               if ( "$v" =~ /^(([\w\s]+?),[-]?(\d+),(\d+))|(([\w\s]+?),[-]?(\d+),U)|(([\w\s]+?),U,[-]?(\d+))|(([\w\s]+?),U,U)/ ){
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$v;";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $v\n";
                  $check .= "Invalid: Symptoms data must be corrected. (Error_93_INVALID_SYMPTOM_OTHER),";
               }
            }
         }
      }
      chop $check;
   }
   ## Chronic_Conditions, Maintenance_Medication, Infections_Within_Five_Years
   ## Maximum allowed length: 250 characters
   ## The entry must be one or more comma-separated members of the Value List. NOTE: User can enter other value by prefixing 'OTH-'
   foreach my $field ("chronic_conditions","maintenance_medication","infections_within_five_years"){
      if ( $header =~ m/$field/ ){
         $check = "";
         my @cval = split (',', $dval);
         foreach my $v (@cval){
            $v =~ s/^\s+//g;
            $v =~ s/\s+$//g;
            if ( $vocab{$field}->{lc $v} ){
               #SERVICE#print STDERR " passed\n";
               $check .= "$vocab{$field}->{lc $v},";
            }else{
               if ( $v =~ m/^OTH-[\s+]?\w+/i ){
                  $val_len = length($v);
                  if ( $val_len > 250 ){
                     #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 250 characters)\n";
                     $check .= "Invalid: Value length is $val_len. Maximum allowed length is 250 characters (Error_75_INVALID_FIELD_LENGTH_OTH),";
                  }else{
                     #SERVICE#print STDERR " passed\n";
                     $v =~ s/oth/OTH/i;
                     $check .= "$v,";
                  }
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $v\n";
                  $check .= "Invalid: Field value must be OTH-text (Error_1_INVALID_VALUE),";
               }
            }
         }
         chop $check;
      }
   }
   ## Types_of_Allergies
   ## Maximum allowed length: 100 characters
   ## NA is not allowed if ALL is selected for Chronic_Conditions.
   if ( $header =~ m/types_of_allergies/ ){
      if ( $val_len > 100 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 100 characters)\n";
         $check = "Invalid:  Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'types_of_allergies'}->{lc $dval} ){
            if ( uc $dval eq "NA" ){
               if ( uc $parse_v eq "ALL" ){
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: NA is not allowed if ALL is selected for Chronic_Conditions. (Error_108_INVALID_ALLERGY_TYPE)";
               }else{
                  #SERVICE#print STDERR " passed\n";
                  $check = "NA";
               }
            }else{
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{'types_of_allergies'}->{lc $dval};
            }
         }else{
            if ( $dval =~ m/^([\w\s\,]+?)/ ){
               #SERVICE#print STDERR " passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   } 
   ## Human_Leukocyte_Antigens,Genbank_Accession_Numbers,HPAI_H5N1
   ## Maximum allowed length: 50 characters
   foreach my $field ("human_leukocyte_antigens","genbank_accession_numbers","hpai_h5n1"){
      if ( $header =~ m/$field/ ){
         if ( $val_len > 50 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 50 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ m/^([\w\s]+?)/ ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }
   ## Packs_Per_Day_for_How_Many_Years
   ## Maximum allowed length: 50 characters
   ## IF Y is selected for tobacco_Use, NA is not allowed.
   ## If N is selected for Tobacco_Use, must be NA.
   if ( $header =~ m/packs_per_day_for_how_many_years/ ){
      if ( $val_len > 50 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 50 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'packs_per_day_for_how_many_years'}->{lc $dval} ){
            if ( uc $dval eq "NA" ){
               if ( uc $parse_v eq "N" ){
                  #SERVICE#print STDERR " passed\n";
                  $check = "NA";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: If Tobacco_Use is 'Y'; then 'NA' is not allowed. (Error_148_NA_NOT_ALLOWED)";
               }
            }else{
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{'packs_per_day_for_how_many_years'}->{lc $dval};
            }
         }else{
            if ( $dval =~ m/^(\d+|(\d+\.?\d+)) (pack|pk)\/day for (\d+|(\d+\.?\d+)) year/i ){
               #SERVICE#print STDERR " passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Travel_History
   ## List all subject's travel destinations over the past month, as listed in the ISO 3166 Standard Country Codes
   ## Maximum allowed length: 200 characters
   ## If multiple values were recorded, separate the individual one with a comma. 
   if ( $header =~ m/travel_history/ ){
      $check = "";
      my @cval = split (',', $dval);
      foreach my $v (@cval){
         $v =~ s/^\s+//g;
         $v =~ s/\s+$//g;
         $val_len = length $v;
         if ( $val_len > 200 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 200 characters)\n";
            $check .= "Invalid: Value length is $val_len. Maximum allowed length is 200 characters (Error_70_INVALID_FIELD_LENGTH),";
         }else{
            if ( $vocab{'travel_history'}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check .= "$vocab{'travel_history'}->{lc $v},";
            }else{
               if ( $country{'code'}->{lc $v} or $country{'name'}->{lc $v} ){
                  #SERVICE#print STDERR " passed\n";
                  $check .= "$country{'code'}->{lc $v},";
               }else{
                  #SERVICE#print STDERR "  invalid country code/name: $v\n";
                  $check .= "Invalid: Field value should be a valid country code. (Error_1_INVALID_VALUE),";
               }
            }
         }
      }
      chop $check;
   }
   ## Primary_Living_Situation,Profession,Host_ID_Type
   ## Maximum allowed length: 50 characters
   foreach my $field ("primary_living_situation","profession","host_id_type"){
      if ( $header =~ m/^$field$/ ){
         if ( $val_len > 50 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 50 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ /^\w+/ ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }
   ## Comments
   ##  Value length is $val_len. Maximum allowed length is 2000 characters
   if ( $header =~ m/comments/ ){
      if ( $val_len > 2000 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 2000 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 2000 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         #SERVICE#print STDERR " passed\n";
      }
   }
   ## Collection_Latitude, Collection_Longitude
   ## Maximum allowed length: 20 characters
   ## The entry must be in WGS84 decimal degree format
   ## latitude <DD> values range from -90 to 90; longitude <DDD> values range from -180 to 180.
   foreach my $field ("collection_latitude","collection_longitude"){
      if ( $header =~ m/^$field$/ ){
         if ( $vocab{$field}->{lc $dval} ){
            #SERVICE#print STDERR "  passed\n";
            $check = $vocab{$field}->{lc $dval};
         }else{
            if ( $dval =~ /[-]?(\d+)$/ ){
               if ( $field =~ /collection_latitude/ ){
                  if ( $dval > 90 or $dval < -90 ){
                     #SERVICE#print STDERR "  invalid data entry: $dval (range from -90 to 90)\n";
                     $check = "Invalid: Field value range should be from -90 to 90 (Error_79_INVALID_WGS84_FORMAT)"; 
                  }else{
                     #SERVICE#print STDERR " passed\n";
                  }
               }
               if ( $field =~ /collection_longitude/ ){
                  if ( $dval > 180 or $dval < -180 ){
                     #SERVICE#print STDERR "  invalid data entry: $dval (range from -180 to 180)\n";
                     $check = "Invalid: Field value range should be from -180 to 180 (Error_79_INVALID_WGS84_FORMAT)";
                  }else{
                     #SERVICE#print STDERR " passed\n";
                  }
               }
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Host_Species
   ## If the entry is not ENV or U, the host species name is validated against the iDPCC Species Dictionary.
   if ( $header =~ m/host_species/ ){
      if ( $vocab{'host_species'}->{lc $dval} ){
         #SERVICE#print STDERR " passed\n";
         $check = $vocab{'host_species'}->{lc $dval};
      }else{
         if ( $species{'host_species'}->{$lcval} ){
            #SERVICE#print STDERR " passed\n";
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Field value doesn't exist in the iDPCC Species Dictionary (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Host_Common_Name
   ## If the entry is U, the host common name is validated against the iDPCC Species Dictionary.
   if ( $header =~ m/host_common_name/ ){
      if ( $vocab{'host_common_name'}->{lc $dval} ){
         #SERVICE#print STDERR " passed\n";
         $check = $vocab{'host_common_name'}->{lc $dval};
      }else{
         if ( $species{'host_name'}->{$lcval} ){
            #SERVICE#print STDERR " passed\n";
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Field value doesn't exist in the iDPCC Species Dictionary (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Host_Identifier
   ## Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _
   ## Maximum allowed length: 50 characters
   if ( $header =~ m/host_identifier/ ){
      if ( $val_len > 50 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 50 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 50 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'host_identifier'}->{lc $dval} ){
            #SERVICE#print STDERR " passed\n";
            $check = $vocab{'host_identifier'}->{lc $dval};
         }else{
            if ( $dval =~ /[a-zA-Z0-9\-\_]/ and $dval !~ /[\'\!\@\#\$\%\^\&\*\(\)\=\+\>\<\?\"\/\\\:\.\,\;]/ ){
               #SERVICE#print STDERR " passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Allowed characters include alphanumeric, hyphen, and underscore: a-z, A-Z, 0-9, -, _ (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Host_Habitat
   ## The entry must be one and only one member of the Value List.
   ## Values are case-sensitive and must be entered in all-caps.
   ## If Host_Natural_State is WLD, then Host_Habitat must be MIG/RES/U/OTH-text.
   ## If Host_Natural_State is DOM, then Host_Habitat must be FRF/FRI/SEH/MAR/U/OTH-text.
   ## If Host_Natural_State is U, then Host_Habitat should be U.
   ## If Host_Natural_State is OTH-text, then Host_Habitat should be Other and provided as OTH-text.
   ## Maximum allowed length: 30 characters
   if ( $header =~ m/host_habitat/ ){
      if ( $val_len > 30 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 30 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 30 characters (Error_70_INVALID_FIELD_LENGTH)";
      }else{
         if ( $vocab{'host_habitat'}->{lc $dval} ){
            if ( uc $parse_v eq "DOM" and uc($dval) =~ /FRF|FRI|SEH|MAR|U/ ){ 
               #SERVICE#print STDERR " passed\n";
               $check = uc $dval;
            }elsif ( uc $parse_v eq "U" and uc $dval eq "U" ){
               #SERVICE#print STDERR " passed\n";
               $check = uc $dval;
            }elsif ( uc $parse_v eq "WLD" and uc($dval) =~ /MIG|RES|U/ ){
               #SERVICE#print STDERR " passed\n";
               $check = uc $dval;
            }elsif ( uc $parse_v eq "CPW" and uc($dval) =~ /MAR|RES|U/ ){
               $check = uc $dval;
            }elsif ( uc $parse_v eq "DMF" and uc $dval eq "U" ){
               $check = uc $dval;
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: If Host_Natural_State is WLD; then Host_Habitat must be MIG/RES/U/OTH-text. If Host_Natural_State is DOM; then Host_Habitat must be FRF/FRI/SEH/MAR/U/OTH-text. If Host_Natural_State is U; then Host_Habitat should be U. (Error_71_HOST_NAT_STATE_WLD/Error_72_HOST_NAT_STATE_DOM/Error_73_HOST_NAT_STATE_U)";
            }
         }else{
            if ( $dval =~ m/^OTH-[\s+]?\w+/i ){
               #SERVICE#print STDERR " passed\n";
               $dval =~ s/oth/OTH/i;
               $check = $dval;
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: If Host_Natural_State is OTH-text; then Host_Habitat should be Other and provided as OTH-text. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Strain_Name,Influenza_Test_Antigen
   ## Maximum allowed length: 150 characters
   ## Influenza A virus: Antigenic Type/Host of Origin/Geographical Origin/Strain Number/Year of Isolation (Subtype)
   ## Influenza B, C, or D virus: Antigenic Type/Host of Origin/Geographical Origin/Strain Number/Year of Isolation
   ## SARS-CoV-2 and other viruses: Virus Name/Host of Origin/Geographical Origin/Strain Number/Year of IsolationAntigenic Type/Host of Origin/Geographical Origin/Strain Number/Year of Isolation
   foreach my $field ("strain_name","influenza_test_antigen"){ 
      if ( $header =~ m/$field/ ){
         if ( $val_len > 150 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 150 characters)\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 150 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ m/^\w+[\/]?/g ){
                  #SERVICE#print STDERR " passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   } 
   ## Influenza_Subtype
   ## Maximum allowed length: 6 characters
   ## The entry must be one and only one member of the Value List.
   ## H and N are case-sensitive and must be entered in all-caps.
   ## If only one of the subtypes have been tested, then use the format H5Nx or HxN1.
   ## Enter HxNx if Influenza_Test_Result is positive but inconclusive.
   ## Enter NA if Influenza_Test_Interpretation is negative."
   if ( $header =~ m/influenza_subtype/ ){
      if ( $val_len > 6 ){
         #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 150 characters)\n";
         $check = "Invalid: Value length is $val_len. Maximum allowed length is 6 character (Error_70_INVALID_FIELD_LENGTH),";
      }else{
         if ( $vocab{'influenza_subtype'}->{lc $dval} ){
            if ( uc $parse_v eq "N" ){
               if ( uc $dval eq "NA" ){
                  #SERVICE#print STDERR " passed\n";
                  $check = "NA";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $v for $dval\n";
                  $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }else{
               #SERVICE#print STDERR " passed\n";
               $check = $vocab{'influenza_subtype'}->{lc $dval};
            }
         }else{
            if ( $dval =~ m/^(H[x|\d]N[x|\d])/ ){
               #SERVICE#print STDERR " passed\n";
            }else{
               #SERVICE#print STDERR "  invalid data entry: $dval\n";
               $check = "Invalid: Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Subject Race
   ## The entry must be one or more comma-separated members of the Value List.
   ## Maximum allowed length: 150 characters
   if ( $header =~ m/subject_race/) {
      $check = "";
      my @deval = split (',', $dval);
      foreach my $t (@deval){
         $t =~ s/^\s+//;
         $t =~ s/\s+$//;
         $val_len = length($t);
         if ( $val_len > 150 ){
            #SERVICE#print STDERR "  invalid data entry: $dval (maximum length: 150 characters)\n";  
            $check .= "Invalid: Value length is $val_len. Maximum allowed length is 150 characters (Error_70_INVALID_FIELD_LENGTH),";
         }else{
            if ( $vocab{subject_race}->{lc $t} ){
               #SERVICE#print STDERR "  passed\n";
               $check .= "$vocab{subject_race}->{lc $t},"; 
            }else{
               #SERVICE#print STDERR "  invalid data entry: ($t)\n";
               $check .= "Invalid: Field value should be one of valid values as in the list. Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),"; 
            } 
         }
      }
      chop $check;
   }   
   ## Longitudinal_Study,Pregnancy,Trimester_of_Pregnancy,Disease_Outcome,Influenza_Like_Illness_Over_the_Past_Year,Breastfeeding,Alcohol_or_Other_Drug_Dependence,Tobacco_Use,Host_Health,Host_Age,Host_Sex
   ## subject_gender,subject_ethnicity
   ## The entry must be one and only one member of the Value List.  Values are case-sensitive and must be entered in all-caps.
   foreach my $field ("longitudinal_study","pregnancy","trimester_of_pregnancy","disease_outcome","influenza_like_illness_over_the_past_year","breastfeeding","alcohol_or_other_drug_dependence","tobacco_use",
                      "host_health","host_age","host_sex","subject_gender","subject_ethnicity"){
      if ( $header =~ m/^$field$/ ){
         if ( $vocab{$field}->{lc $dval} ){
            #SERVICE#print STDERR "  passed\n";
            $check = $vocab{$field}->{lc $dval};
         }else{
            #SERVICE#print STDERR "  invalid data entry: $dval\n";
            $check = "Invalid: Field value should be one of valid values as in the list. Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
         }
      }
   }
   ## Subject_Height,Subject_Weight
   ## Maximum length: 17 characters
   foreach my $field ("subject_height","subject_weight"){
      if ( $header =~ m/^$field$/ ){
         if ( $val_len > 17 ) {
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR "  passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               if ( $dval =~ /^\d+(\.\d+)?$/ ){
                  #SERVICE#print STDERR "  passed\n";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: ($dval)\n";
                  $check = "Invalid: Field value should be one of valid values as in the list. Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
               }
            }
         }
      }
   }
   ## Nursing_Home_Residence,Daycare_Attendance,Disease_Status,WHO_COVID_Disease_Severity
   ## Maximum length: 17 characters
   foreach my $field ("nursing_home_residence","daycare_attendance","disease_status","who_covid_disease_severity"){
      if ( $header =~ m/^$field$/ ){
         if ( $val_len > 17 ) {
            #SERVICE#print STDERR "  invalid data entry: maximum length: 17 characters\n";
            $check = "Invalid: Value length is $val_len. Maximum allowed length is 17 characters (Error_70_INVALID_FIELD_LENGTH)";
         }else{
            if ( $vocab{$field}->{lc $dval} ){
               #SERVICE#print STDERR "  passed\n";
               $check = $vocab{$field}->{lc $dval};
            }else{
               #SERVICE#print STDERR "  invalid data entry: ($dval)\n";
               $check = "Invalid: Field value should be one of valid values as in the list. Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE)";
            }
         }
      }
   }
   ## Chest_Imaging_Interpretation
   ## Maximum length: 100 characters
   ## The entry must be one or more comma-separated members of the Value List. NRM should be entered as s single value only.
   ## NA should be entered if the subject was not administered a chest x-ray.
   if ( $header =~ m/chest_imaging_interpretation/ ) {
      $check = "";
      my @deval = split (',', $dval);
      my $c = @deval;
      foreach my $t (@deval){
         $t =~ s/^\s+//;
         $t =~ s/\s+$//;
         $val_len = length($t);
         if ( $val_len > 100 ){
            #SERVICE#print STDERR "  invalid data entry: maximum length: 100 characters\n";
            $check .= "Invalid: Value length is $val_len. Maximum allowed length is 100 characters (Error_70_INVALID_FIELD_LENGTH),";
         }else{
            if ( $vocab{chest_imaging_interpretation}->{lc $t} ){
               if ( uc $t eq "NRM" and $c > 1 ) {
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: NRM should be entered as s single value only. (Error_1_INVALID_VALUE),";
               }else{
                  #SERVICE#print STDERR "  passed\n";
                  $check .= "$vocab{chest_imaging_interpretation}->{lc $t},";
               }
            }else{
               if ( $t =~ m/^OTH-[\s+]?\w+/i ){
                  #SERVICE#print STDERR " passed\n";
                  $dval =~ s/oth/OTH/i;
                  $check .= "$t,";
               }else{
                  #SERVICE#print STDERR "  invalid data entry: $dval\n";
                  $check = "Invalid: Field value should be one of valid values as in the list. Please correct the value based on iDPCC data standards. (Error_1_INVALID_VALUE),";
               }
            }
         }
      }
      chop $check;
   }

   return $check;
}
