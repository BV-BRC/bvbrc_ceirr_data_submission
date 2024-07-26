# CEIRR Data Submission Service
The CEIRR Data Submission Service allows users to validate and submit virus surveillance data in human and animal subjects, as well as serological and viral isolate data. This service follows [CEIRR Data Standards](https://www.ceirr-network.org/resources/data-standards).

> Current pipeline supports Animal Surveillance, Human Surveillance, Human Surveillance-Southern Hemisphere, Serological Data, and Viral Isolate Data.

### Data Template Headers

For the CEIRR Data Submission Service to identify the type of submission, the data template header must be included in the CEIRR data file. The following headers should be used:

- **#DataTemplate:Animal_Surveillance_v2.3**
- **#DataTemplate:Human_Surveillance_v2.6**
- **#DataTemplate:Serological_Data_v2.2**
- **#DataTemplate:Viral_Isolate_Data_v2.2**

Example JSON input for the service:

```
{
    "ceirr_data": "/olson@patricbrc.org/PATRIC-QA/DataSets/CEIRRDataSubmission/animal_surveillance_v2_3_success.csv",
    "output_file": "animal_success",
    "output_path": "/olson@patricbrc.org/home/testn"
}
```
