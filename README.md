# iSeq_EC_LM_mlv
Analysis files for E. coli and Listeria iSeq multi-laboratory project

## How to run
Create a directory with the fastq files from the run you want to analyze and run the following:

`iSeq_ECLM_analysis.py *_R1_001.fastq.gz`

Among other files, this should produce one called iSeq_report_\<Date-Time\>.docx. Download this file and update the relevant sections to complete the report
  
## How it works
When the user calls iSeq_ECLM_analysis.py, it will generate two bash files: LaunchJobs_\<Date-Time\>.sh and JobList_\<Date-Time\>. JobList is the set of commands to be run on each sample. After ensuring that the environment is set properly, LaunchJobs will run these commands in parallel, produce summary output tables, and generate the report.
 
 In the JobList file, you can see that the following are done for each sample:
 
 * get_read_metrics.py: Runs FastQC, parses the output, and feeds it to a temporary file called read_metrics_nohead.tsv.  
 * SKESA: Creates a *de novo* assembly the reads, saved as \<Sample ID\>_\<Machine ID\>.fasta.
 * Quast: Computes assembly quality metrics and feeds them to a temporary file called assembly_stats_nohead.tsv.  
 * AMRFinder: Identifies AMR genes present in each sample. A temporary file is generated containing the genes detected for comparison to those dectected in the reference geneome.
 * ABRicate: Searches *E. coli* assemblies for virulence factors in [VPECdb](https://github.com/pkmitchell/VPECdb). A similar temporary file is generate to that from the previous step.
 * check_genes.py: Compares the AMR and virulence genes from the previous two steps against those found in the reference sequence, feeding summary output to a temporary file called tmp_AnnotationTable_\<Date-Time\>.tsv. Other temporary files associated with this step are then deleted.
 * CFSAN SNP Pipeline: After a few setup commands, the CFSAN SNP pipeline is launched to map reads to the reference sequence.
 
 After these steps are complete for all samples, LaunchJobs will convert the temporary output files for read QC, assembly QC, and gene annotation to their final forms and generate an output table of alignment metrics from the SNP pipeline output. It will also run MLST and E. coli serotype predction and create a tables of their output. These are then merged to create the subtype table. Mixed in with these steps are some commands to clean up output and to record software versions. 
 
 Finally, LaunchJobs copies the R markdown file, updates it to include the \<Date-Time\> identifiers used for this run, and generates the output word file.
 


## Setup of this repo
Note that the this repo currently serves as an accessible secondary location for files stored on cbsuahdc01. As such, some hard-coded file paths do not correspond to the directory structure of this repository.
### scripts
This directory contains the main script used to run the pipeline, iSeq_ECLM_analysis.py, as well as two helper scripts that it calls.
### R_templates
This directory contains the R markdown template used to generate the output report and a word file used as a template for formatting
### conda_environment
This directory contains yml and spec files for the conda environment built for this analysis



