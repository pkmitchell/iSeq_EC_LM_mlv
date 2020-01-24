#!/workdir/miniconda3/envs/iSeq_EC_LM_mlv/bin/python

import sys
import string
import time
from subprocess import call
from subprocess import check_output
import argparse
import re
import multiprocessing as mp


#Argument Parser
parser = argparse.ArgumentParser(description="Assemble genomes using SKESA (or, optionally, SPAdes), or generate scripts to do so but do not submit immediately")
parser.add_argument("-w", "--wait", help="Launch job later (Default = Off)", action="store_true")
parser.add_argument("-t", "--trim", help="Turn on read trimming (trimmomatic ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20)", action="store_true")
parser.add_argument("-q", "--no_qc", help="Turn off default mapping/assembly QC output", action="store_true")
parser.add_argument("-p", "--threads", help="Maximum number of processors to use (Default = 8)", type = int, default=8)
#parser.add_argument("-n", "--jobs", help="Maximum number of jobs to run simultaneously, equivalent to perl_fork_univ.pl argument (Default = 6)", type = int, default=6)
parser.add_argument("-m", "--max_mem", help="Maximum total memory allocation to simultaneous processes in GB (Default = 24)", type = int, default=24)
parser.add_argument("-d", "--spades", help="Use SPAdes instead of SKESA for assembly (Default = Off)", action="store_true")
parser.add_argument("Fastq", help="Input forward fastqs to be assembled", nargs='+')
args = parser.parse_args()

#Creating dictionary to match file stems to reference sequences
RefDict = {"804-002-iSeq-EC-11":"ECP19-198", "804-002-iSeq-EC-12":"ECP19-598", "804-002-iSeq-EC-13":"ECP19-798", "804-002-iSeq-EC-14":"ECP19-2498","804-002-iSeq-LM-19":"LMP18-H8393","804-002-iSeq-LM-20":"LMP18-H2446", "804.002-iSeq-EC-15":"N17EC0320", "804.002-iSeq-EC-16":"N17EC0326","804.002-iSeq-EC-17":"N17ECO616","804.002-iSeq-EC-18":"N17EC1164"}
RefDir = "/workdir/iSeq_Ecoli/ref_output/"

#Determining the number of jobs and threads per job
ncpus = args.threads #The total number of CPUs to be used simultaneously across all jobs
tot_cpus = mp.cpu_count()
if ncpus > tot_cpus:
	ncpus = tot_cpus

#njobs = args.jobs #The number of jobs to be run simultaneously
totjobs = len(args.Fastq) #The total number of input sequences, which is also the total number of jobs to run

#Determining the number of threads per job and number of jobs to run simultaneously using the bounds provided by the user
#if njobs > totjobs: 
#	njobs = totjobs
	
memper = int(args.max_mem)#/njobs)

#This is a bit dumb and clunky, might be worth improving at some point to handle different input scenarios
if memper <=10:
	print("Seems like not enough memory for the number of jobs requested.\n")
	print("Reducing number of simultaneous processes")
	memper= 16
	njobs= int(args.max_mem/memper)
elif memper < 16:
	print("Might not be enough memory, trying anyway")

#if njobs > ncpus:
#	njobs = ncpus
	
threadsper = int(ncpus)#/njobs)

strT = time.strftime("%Y%m%d_%H%M%S")
#This creates an input file for perl_fork_univ.pl with one line per input sequence
#It is then called by the launch script generated in the next step

fn1="JobList_" + strT
ns=open(fn1,'w')

prefs = []

for i in range(0, totjobs):
	fs=args.Fastq[i]
	rs=None
	stem=re.split(r'_S[0-9][0-9]*_L001_R1_001.fastq', fs)[0]
	rs=fs.replace('_R1_', '_R2_')
	RefStem = RefDict[stem]
	
	mID = check_output(["zgrep", "-m1", "-o", "FS[0-9]\+", fs])
	mID = mID.decode().strip()
	
	prefix=stem.split("/")[-1] + "_" + mID
	prefs.append(prefix)
	
	if "EC" in prefix:
		species = "EC"
	elif "LM" in prefix:
		species = "LM"
		
	ns.write("get_read_metrics.py -p 2 -s " + species + " -i " + prefix + " " + fs + " " + rs + " >>read_metrics_nohead.tsv && ")
	if args.trim == True:
		ns.write("trimmomatic PE -threads " + str(threadsper) + " -phred33 " + fs + " " + rs + " trim_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_1.fastq.gz trim_" + prefix + "_2.fastq.gz unpaired_" + prefix + "_2.fastq.gz ILLUMINACLIP:/workdir/miniconda3/envs/bacWGS/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 &>trim_" + prefix + ".log &&")
		ns.write(" rm unpaired_" + prefix + "_1.fastq.gz unpaired_" + prefix + "_2.fastq.gz &&")
		if args.spades == True:
			ns.write(" spades.py --pe1-1 " + "trim_" + prefix + "_1.fastq.gz --pe1-2 trim_" + prefix + "_2.fastq.gz -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&")
			ns.write(" quast.py -r " + RefDir + RefStem + ".fasta --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + "/contigs.fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv &&")
			ns.write(" mv " + prefix + "/contigs.fasta ./" + prefix + ".fasta && mv " + prefix + "/spades.log ./" + prefix + "_spades.log")
		else:
			ns.write(" skesa --cores " + str(threadsper) + " --memory " + str(memper) + " --fastq " + "trim_" + prefix + "_1.fastq.gz trim_" + prefix + "_2.fastq.gz >" + prefix + ".fasta 2>" + prefix + "_skesa.log &&")
			ns.write(" quast.py -r " + RefDir + RefStem + ".fasta --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + ".fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv")
	else:
		if args.spades == True:
			ns.write(" spades.py --pe1-1 " + fs + " --pe1-2 " + rs + " -o ./" + prefix + "/ --careful -t " + str(threadsper) + " -m " + str(memper) + " &&")
			ns.write(" quast.py -r " + RefDir + RefStem + ".fasta --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + "/contigs.fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv &&")
			ns.write(" mv " + prefix + "/contigs.fasta ./" + prefix + ".fasta && mv " + prefix + "/spades.log ./" + prefix + "_spades.log")
		else:
			ns.write(" skesa --cores " + str(threadsper) + " --memory " + str(memper) + " --fastq " + fs + " " + rs + " >" + prefix + ".fasta 2>" + prefix + "_skesa.log &&")
			ns.write(" quast.py -r " + RefDir + RefStem + ".fasta --threads " + str(threadsper) + " -o " + prefix + " --min-contig 1 -l \"" + prefix + "\" --contig-thresholds 500,1000 " + prefix + ".fasta && tail -n 1 " + prefix + "/transposed_report.tsv >>assembly_stats_nohead.tsv")
	###if args.no_qc == False:
	###	ns.write(" && bbmap.sh -Xmx"+ str(memper) + "g nodisk ref=" + RefStem + ".fasta in1=" + fs + " in2=" + rs + " covstats=" + prefix + "_covstats.txt covhist=" + prefix + "_covhist.txt t=" + str(threadsper) + " &>" + prefix + "_bbmap.log")

	ns.write(" && amrfinder --threads " + str(threadsper))
	if species == "EC":
		ns.write(" -O Escherichia")
	ns.write(" -n " + prefix + ".fasta >" + prefix + "_amrfinder.tsv 2>" + prefix + "_amrfinder.log")
	ns.write(" && cut -d$'\t' -f 6 " + prefix + "_amrfinder.tsv|sort >tmp" + prefix + "_ARGs")
	if species == "EC":
		ns.write(" && abricate --db VPECdb --mincov 50 --minid 90 --threads " + str(threadsper) + " " + prefix + ".fasta" + ">" + prefix + "_abricate_VPECdb.tsv 2>" + prefix +"_abricate_VPECdb.log")
		ns.write(" && cut -d$'\t' -f 6 " + prefix + "_abricate_VPECdb.tsv|sort >tmp" + prefix + "_VFs")
	ns.write(" && check_genes.py tmp" + prefix + "_ARGs " + RefDir + RefStem + "_ARGs ")
	if species == "EC":
		ns.write("tmp" + prefix + "_VFs " + RefDir + RefStem + "_VFs ")
	ns.write(">>tmp_AnnotationTable_" + strT + ".tsv")
	ns.write(" && rm tmp" + prefix + "*")
	#Setup and run cfsan snp pipeline
	ns.write(" && mkdir " + prefix + "_snp && mkdir " + prefix + "_snp/fastqs && cp -l " + stem + "*.fastq.gz " + prefix + "_snp/fastqs && run_snp_pipeline.sh -c snppipeline.conf -o " + prefix + "_snp -s " + prefix + "_snp " + RefDir + RefStem + ".fasta")
	ns.write(";\n")

ns.close()
 
#This step generates the bash script to set the environment and launch perl_fork_univ.pl, calling the previously generated input file 

fn2="./LaunchJobs_" + strT + ".sh"
ss=open(fn2, 'w')

ss.write("#!/bin/bash \n")
ss.write("\n")
ss.write("export PATH=/workdir/miniconda3/bin:$PATH\n\n")
ss.write("source activate iSeq_EC_LM_mlv\n\n")
ss.write("cfsan_snp_pipeline data configurationFile\n")
ss.write("export CLASSPATH=/local/workdir/miniconda3/envs/iSeq_EC_LM_mlv/share/varscan-2.4.4-0/VarScan.jar:$CLASSPATH\n")
ss.write("export CLASSPATH=/local/workdir/miniconda3/envs/iSeq_EC_LM_mlv/share/picard-2.21.6-0/picard.jar:$CLASSPATH\n")
ss.write("export CLASSPATH=/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar:$CLASSPATH\n")
ss.write("sed -i 's/MaxCpuCores=/MaxCpuCores=" + str(threadsper) + "/' snppipeline.conf\n")
ss.write("\nmkdir FastQC_output\n")
ss.write("\nsh ./" + fn1 + " &>" + fn1 + ".log \n")
ss.write("(head -n 1 " + prefs[0] + "/transposed_report.tsv && cat assembly_stats_nohead.tsv) >assembly_stats.tsv && rm assembly_stats_nohead.tsv\n")

ss.write("\nfastqc --version >>software_versions_" + strT + ".txt\n")
ss.write("echo \"Isolate\tReads\tF Length\tR Length\tF Q\tR Q\tEst. Coverage\" >ReadMetrics_" + strT + ".tsv\n")
ss.write("sort read_metrics_nohead.tsv >>ReadMetrics_" + strT + ".tsv && rm read_metrics_nohead.tsv\n")
if args.trim == True:
	#ss.write("\necho \"trimmomatic $(trimmomatic -version)\" >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir Trimmomatic_output && mv trim* Trimmomatic_output\n")
if args.spades == True:
	#ss.write("\nspades.py -v >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir SPAdes_logs && mv *_spades.log SPAdes_logs\n")
else:
	#ss.write("\nskesa -v >>software_versions_" + strT + ".txt\n")
	ss.write("mkdir SKESA_logs && mv *_skesa.log SKESA_logs\n")
#ss.write("\nquast.py -v >>software_versions_" + strT + ".txt\n")
if args.no_qc == False:
	ss.write("echo -e \"Isolate\tTotal Length\tGenome Fraction (%)\tUnaligned Length\tContigs\tN50\tNG50\" >AssemblyMetrics_" + strT + ".tsv\n")
	ss.write("for s in " + ' '.join(prefs) + ";\n do i=$(grep $s assembly_stats.tsv |awk -F$'\t' '{print $1,\"\t\",$8,\"\t\",$29,\"\t\",$28,\"\t\",$6,\"\t\",$12,\"\t\",$13}');\n echo -e \"$i\";done >>AssemblyMetrics_" + strT + ".tsv\n")

#Finish annotation table and Clean up AMRFinder and ABRicate output
ss.write("\necho -e \"Isolate\tMissed ARGs\tUnexpected ARGs\tMissed VFs\tUnexpected VFs\" >AnnotationTable_" + strT + ".tsv\n")
ss.write("sort tmp_AnnotationTable_" + strT + ".tsv >>AnnotationTable_" + strT + ".tsv && rm tmp_AnnotationTable_" + strT + ".tsv\n")
#ss.write("\nconda list|grep \"amrfinder\"|tr -s ' ' |cut -d \" \" -f 1,2 >>software_versions_" + strT + ".txt\n")
#ss.write("echo \"AMRFinder database$(ls -l $(which amrfinder|sed 's/amrfinder/data/')|grep \"latest\"|cut -d \">\" -f 2)\" >>software_versions_" + strT + ".txt\n")
#amrfiles = [p + "_amrfinder.tsv" for p in prefs]
#ss.write("make_amr_table.py " +  ' '.join(amrfiles) + " >ARG_table_" + strT + ".tsv\n")
ss.write("mkdir AMRFinder_output && mv *amrfinder* AMRFinder_output\n")
ss.write("abricate -V >>software_versions_" + strT + ".txt\n")
ss.write("mkdir ABRicate_output && mv *_abricate_VPECdb.* ABRicate_output\n")

#Get alignment metrics
ss.write("\necho -e \"Isolate\tReads Mapped (%)\tInsert Size\tRead Depth\tSNPs\" >AlignmentMetrics_" + strT + ".tsv\n")
ss.write("for s in " + ' '.join(prefs) + ";\n do tail -n 1 ${s}_snp/metrics.tsv|awk -v s=\"$s\" -F$'\t' '{print s,\"\t\",$8,\"\t\",$10,\"\t\",$11,\"\t\",$15}';done >>AlignmentMetrics_" + strT + ".tsv\n")


#Run MLST on all assemblies
mlst_in = [p + ".fasta" for p in prefs]
ss.write("\necho -e \"Isolate\tST\tAlleles\" >MLST_" + strT + ".tsv\n")
ss.write("mlst --threads " + str(ncpus) + " " + ' '.join(mlst_in) + " 2>MLST_" + strT + ".log|cut -d$'\t' -f 1,3-|sed 's/.fasta//g' >>MLST_" + strT + ".tsv\n")
ss.write("sed -i 's/)\t/) /g' MLST_" + strT + ".tsv\n")
#ss.write("mlst --version >>software_versions_" + strT + ".txt\n")

#Run ecytper on E. coli assemblies
EC_re = re.compile(".*-EC-.*")
ec_prefs = list(filter(EC_re.match, prefs))
#ss.write("\nectyper -V >>software_versions_" + strT + ".txt\n")
fastas = [p + ".fasta" for p in ec_prefs]
ss.write("ectyper -i " + ','.join(fastas) + " --cores " + str(ncpus) + " --verify -o ectyper_out_" + strT + " &>ectyper_" + strT + ".log\n")
ss.write("cut -d$'\t' -f 1,2,3 ectyper_out_" + strT + "/output.tsv >ECOL_sero_" + strT + ".tsv\n")
ss.write("mv ectyper_" + strT + ".log ectyper_out_" + strT + "\n")

#Make subtyping table
ss.write("\njoin -j 1 --header -e NA -a2 -o 2.1,1.2,1.3,2.2,2.3  -t$'\t' ECOL_sero_" + strT + ".tsv MLST_" + strT + ".tsv >SubtypeTable_" + strT + ".tsv\n")

#Make software version table
ss.write("\necho -e \"Read QC\t$(fastqc -v)\" >Software_table_" + strT + ".tsv\n")
ss.write("echo -e \"Assembly and QC\t$(skesa -v), $(quast.py -v)\" >>Software_table_" + strT + ".tsv\n")
ss.write("echo -e \"Alignment\t$(run_snp_pipeline.sh --version)\" >>Software_table_" + strT + ".tsv\n")
ss.write("echo -d \"Subtyping\t$(mlst -v), $(ectyper -V)\" >>Software_table_" + strT + ".tsv\n")
ss.write("echo -e \"Gene Annotation\t$(conda list|grep \"amrfinder\"|tr -s ' ' |cut -d \" \" -f 1,2)")
ss.write(" db:$(ls -l $(which amrfinder|sed 's/amrfinder/data/')|grep \"latest\"|cut -d \">\" -f 2|cut -d \" \" -f 2)")
ss.write(", $(abricate -V)\" >>Software_table_" + strT + ".tsv\n")


#Generate report
ss.write("\nsed 's/DATETIME/" + strT + "/g' /workdir/iSeq_Ecoli/report_templates/testreport.Rmd >iSeq_report_" + strT + ".Rmd\n")
ss.write("/programs/R-3.5.0/bin/Rscript -e \"rmarkdown::render('iSeq_report_" + strT + ".Rmd', clean=TRUE)\"\n")

ss.write("\nconda deactivate")
ss.close()
call(["chmod", "u+x", fn2])

#Optionally run the job
if not args.wait:
	call([fn2])
else:
	print("To run, call " + fn2)