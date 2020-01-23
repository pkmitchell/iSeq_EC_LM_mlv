#!/workdir/miniconda3/envs/iSeq_EC_LM_mlv/bin/python

import sys
import subprocess
import csv

stem=sys.argv[1].replace("_ARGs", "")
stem=stem.replace("tmp", "")
oname = "tmp_comm_" + stem
subprocess.call(["comm", sys.argv[1], sys.argv[2]], stdout=open(oname,'w'))
csvr=csv.reader(open(oname,'r'),delimiter="\t")
missing=[]
unexpected=[]

for line in csvr:
	unexpected.append(line[0].strip())
	if len(line) > 1:
		missing.append(line[1].strip())

subprocess.call(["rm", oname])
	
UG = ', '.join(list(filter(None, unexpected)))
if not UG:
	UG = "None"
	
MG = ', '.join(list(filter(None, missing)))
if not MG:
	MG = "None"

if len(sys.argv)==5:
	vf_oname = "tmp_comm_" + stem + "_VFs" 
	subprocess.call(["comm", sys.argv[3], sys.argv[4]], stdout=open(vf_oname,'w'))
	vf_csvr=csv.reader(open(vf_oname,'r'),delimiter="\t")
	vf_missing=[]
	vf_unexpected=[]

	for line in vf_csvr:
		vf_unexpected.append(line[0].strip())
		if len(line) > 1:
			vf_missing.append(line[1].strip())
		
	subprocess.call(["rm", vf_oname])
	vf_UG = ', '.join(list(filter(None, vf_unexpected)))
	if not vf_UG:
		vf_UG = "None"
		
	vf_MG = ', '.join(list(filter(None, vf_missing)))
	if not vf_MG:
		vf_MG = "None"
else:
	vf_UG = "NA"
	vf_MG = "NA"
	
sys.stdout.write(stem + "\t" + MG + "\t" + UG + "\t" + vf_MG + "\t" + vf_UG + "\n")
sys.stdout.close()