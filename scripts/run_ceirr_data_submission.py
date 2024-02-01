#!/usr/bin/env python

import argparse
import json
import os
import subprocess
import sys

VALIDATE_DATA_SCRIPT = "validate_ceirr_data"
PROCESS_DATA_SCRIPT = "process_ceirr_data"
CONVERT_DATA_SCRIPT = "convert_ceirr_data"

INPUT_VALID_FILE_NAME = "sample_valid.csv"
INPUT_PROCESSED_FILE_NAME = "sample_processed.csv"
DATA_TEMPLATE_HEADER = "#DataTemplate:"

if "P3_BASE_URL" in os.environ:
  SOLR_DB_URL = "{0}/services/SolrProxy".format(os.environ["P3_BASE_URL"])
else:
  SOLR_DB_URL = "http://ash.cels.anl.gov:7099"
SOLR_CORE_SURVEILLANCE_NAME = "surveillance_new"
SOLR_CORE_SEROLOGY_NAME = "serology_new"

def determineType(input_file):
  with open(input_file) as f:
    first_line = f.readline()
    template_type = first_line.replace(DATA_TEMPLATE_HEADER, '')
    if template_type[0] == "\"":
      template_type = template_type[1:]

    if template_type.startswith('Animal_Surveillance'):
      return 'animal'
    elif template_type.startswith('Human_Surveillance-Southern_Hemisphere'):
      return 'human-sh'
    elif template_type.startswith('Human_Surveillance'):
      return 'human'
    elif template_type.startswith('Serological_Data'):
      return 'serology'
    elif template_type.startswith('Viral_Isolate'):
      return 'viral'
    else:
      print("Invalid data template definition. Example: #DataTemplate:Animal_Surveillance_v2.3, #DataTemplate:Human_Surveillance_v2.5, #DataTemplate:Serological_Data_v2.2, #DataTemplate:Viral_Isolate_Data_v2.2")
      sys.exit(-1) 
 
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Surveillance Data Submission Script")
  parser.add_argument("-j", "--jfile", help="json file for job", required=True)
  parser.add_argument("-o", "--output", help="Output directory. defaults to current directory", required=False, default=".")

  args = parser.parse_args()

  #Load job data
  job_data = None
  try:
    with open(args.jfile, "r") as j:
      job_data = json.load(j)
  except Exception as e:
    print("Error in opening job file:\n %s" %(e))
    sys.exit(-1)

  if not job_data:
    print("job_data is null")
    sys.exit(-1)
  print(job_data)

  #Setup output directory
  output_dir = args.output
  output_dir = os.path.abspath(output_dir)
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  os.chdir(output_dir)

  output_file = os.path.join(output_dir, job_data["output_file"] + ".txt")

  #Create input file 
  input_file_name = job_data["ceirr_data"].rsplit("/", 1)[1] 
  input_file = os.path.join(output_dir, input_file_name)
  #Fetch CEIRR data file from workspace
  try:
    fetch_ceirr_data_cmd = ["p3-cp", "ws:%s" %(job_data["ceirr_data"]), input_file]
    subprocess.check_call(fetch_ceirr_data_cmd, shell=False)
  except Exception as e:
    print("Error copying ceirr data file from workspace:\n %s" %(e))
    sys.exit(-1)

  if os.path.getsize(input_file) == 0:
    print("Surveillance data file is empty.")
    sys.exit(-1)

  submission_type = determineType(input_file)
  print("Submission type is %s for %s" %(submission_type, job_data["ceirr_data"]))
  if submission_type == "human" or submission_type == "animal" or submission_type == "human-sh":
    solr_core_name = SOLR_CORE_SURVEILLANCE_NAME 
  else:
    solr_core_name = SOLR_CORE_SEROLOGY_NAME

  #Validate data
  try:
    validate_data_cmd = [VALIDATE_DATA_SCRIPT, "-f", input_file, "-t", submission_type]
    valid_sample_count = subprocess.check_output(validate_data_cmd, shell=False)
  except Exception as e:
    print("Error running validate_ceirr_data for %s:\n %s" %(input_file, e))
    sys.exit(-1)

  if valid_sample_count == '0':
    print("There is no valid sample submitted. Please fix the validation errors and resubmit your data.\n")
  else:
    input_file_valid = os.path.join(output_dir, INPUT_VALID_FILE_NAME) 

    if os.path.getsize(input_file_valid) == 0:
      print("Valid input file is empty.")
      sys.exit(-1)

    #Process data
    try:
      process_data_cmd = [PROCESS_DATA_SCRIPT, "-f", input_file_valid, "-t", submission_type, "-u", input_file]
      processed_sample_count = subprocess.check_output(process_data_cmd, shell=False)
    except Exception as e:
      print("Error running process_ceirr_data for %s:\n %s" %(input_file_valid, e))
      sys.exit(-1)

    #Convert data
    try:
      input_file_processed = os.path.join(output_dir, INPUT_PROCESSED_FILE_NAME)

      if os.path.getsize(input_file_processed) == 0:
        print("Processed input file is empty.")

      convert_data_cmd = [CONVERT_DATA_SCRIPT, "-f", input_file_processed, "-t", submission_type, '-n', processed_sample_count]
      json_file_name = subprocess.check_output(convert_data_cmd, shell=False)
    except Exception as e:
      print("Error running convert_ceirr_data for %s:\n %s" %(input_file_valid, e))
      sys.exit(-1)

    #Post data to SOLR
    try:
      json_file_path = os.path.join(output_dir, os.fsdecode(json_file_name))

      if os.path.getsize(json_file_path) == 0:
        print("JSON file is empty: %s" %(json_file_path))
        sys.exit(-1)
      else:
        print("Inserting new data from %s" %(json_file_path))
        submit_data_cmd = ["p3-solr-insert", "--insert", solr_core_name, json_file_path] 
        if SOLR_DB_URL:
          submit_data_cmd.insert(1, "--url")
          submit_data_cmd.insert(2, SOLR_DB_URL)
        output = subprocess.check_output(submit_data_cmd, shell=False)
        print("Inserted new data from %s" %(json_file_path))
    except Exception as e:
      print("Error posting %s to SOLR:\n %s" %(json_file_path, e))
      sys.exit(-1)
