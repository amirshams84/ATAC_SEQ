# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Project: ENCODE ATAC_SEQ
# Aim: Snakemake workflow for pre_process
# snakemake --snakefile ./snakefile/pre_process.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile ./snakefile/pre_process.py --configfile Encode.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# ################################### IMPORT ##################################### #


import os
import sys
import re
import glob
from snakemake.utils import read_job_properties
library_path = os.path.abspath(workflow.basedir + "/../library/")
sys.path.append(library_path)
import utility
# ################################### FUNCTIONS ################################## #


def get_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + SAMPLE_DELIMITER + "*" + "." + SAMPLE_SUFFIX, recursive=True)


def get_forward_fastq(wildcards):
	"""
	"""
	forward_fastq_List = []
	for sample, sample_Dict in sample_treatment_Dict.items():
		#
		if sample_Dict[TREATMENT_COLUMN] == wildcards.design and sample == wildcards.sample:
			#
			forward_fastq_List.extend(glob.glob(DATADIR + "/**/" + sample + SAMPLE_DELIMITER + "*" + FORWARD_DELIMITER + "*" + "." + SAMPLE_SUFFIX, recursive=True))
		else:
			pass
	return forward_fastq_List
# ################################### CONFIGURATION ############################## #


# +++++++++++++++++++++++++++++++++++++
#MAIN PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../bash_script")
R_Script_path = os.path.abspath(workflow.basedir + "/../R_script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../python_script")
Script_Path = os.path.abspath(workflow.basedir + "/../template")
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
INFOLINK = config_general_Dict["INFOLINK"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#DIRECTORY
config_directory_Dict = config["DIRECTORY"]
WORKDIR = utility.fix_path(config_directory_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_directory_Dict["DATADIR"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
FORMAT = config_data_Dict["FORMAT"].lower()
LAYOUT = config_data_Dict["LAYOUT"].lower()
SAMPLE_DELIMITER = config_data_Dict["SAMPLE_DELIMITER"]
SAMPLE_SUFFIX = config_data_Dict["SAMPLE_SUFFIX"].lower()
GENOME = config_data_Dict["GENOME"].lower()
####
if LAYOUT == "paired":
	#
	if PLATFORM == "illumina":
		FORWARD_DELIMITER = "R1"
		REVERSE_DELIMITER = "R2"
	elif PLATFORM == "sra":
		FORWARD_DELIMITER = ".1"
		REVERSE_DELIMITER = ".2"
	else:
		print("This pipeline only can be used with illumina-sra paired-end fastq format")
		print("Aborting!!!!")
		sys.exit(2)
elif LAYOUT == "single":
	FORWARD_DELIMITER = ""
	REVERSE_DELIMITER = ""
else:
	print("This pipeline only can be used with single/paired-end fastq format")
	print("Aborting!!!!")
	sys.exit(2)
####
if SAMPLE_SUFFIX[0] == ".":
	SAMPLE_SUFFIX = SAMPLE_SUFFIX[1:]
if SAMPLE_SUFFIX not in ["fastq", "fastq.gz"]:
	print("This pipeline only can bed used with fastq or fastq.gz data")
	print("Aborting!!!!")
	sys.exit(2)
else:
	pass
####
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#SUBSAMPLE
config_subsample_Dict = config["SUBSAMPLE"]
SUBSAMPLE_PCT = config_subsample_Dict["SUBSAMPLE_PCT"]
USE_SUBSAMPLE = config_subsample_Dict["USE_SUBSAMPLE"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#UTILITIES
config_utilities_Dict = config["UTILITIES"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
PROCESSORS, MEMORY = utility.get_cluster_info(sys.argv)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
SAMPLE_COLUMN = config_metadata_Dict["SAMPLE_COLUMN"]
TREATMENT_COLUMN = config_metadata_Dict["TREATMENT_COLUMN"]
TREATMENT_LIST = list(config_metadata_Dict["TREATMENT_LIST"])
sample_treatment_Dict = utility.build_sample_treatment_dict(METADATA_FILE, SAMPLE_COLUMN)
metadata_Dict = utility.build_metadata_dict(sample_treatment_Dict, TREATMENT_COLUMN, TREATMENT_LIST)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PRE_PROCESS
config_pre_process_Dict = config["PRE_PROCESS"]
# ------------------------------------
# ################################### WILDCARDS ################################ #


pre_process_List = []
for sample, sample_Dict in sample_treatment_Dict.items():
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		pre_process_List
# ################################### PIPELINE RULES ########################## #


if LAYOUT == "paired":
	#
	rule pre_process_paired:
		"""
		"""
		input:
			fastq_List = get_forward_fastq
		output:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq",
		threads: PROCESSORS
		message: "pre_process_paired: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_reverse = each_fastq.replace(FORWARD_DELIMITER, REVERSE_DELIMITER)
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + SAMPLE_SUFFIX, "", each_fastq_basename)
				shell("""
					#
					##
					total_start_time="$(date -u +%s)"

					RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/pre_process
					mkdir -p $RESULT_PATH
					
					REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $REPORT_PATH
					##
					#
					printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
					printf "%s\\n" "pre_process_paired: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
					printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
					printf "INPUT1: %s\\n" "{each_fastq}"  | tee >(cat >&2)
					printf "INPUT2: %s\\n" "{each_fastq_reverse}"  | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "OUTPUT1: %s\\n" "{output.processed_fwd_fastq}.tmp"  | tee >(cat >&2)
					printf "OUTPUT2: %s\\n" "{output.processed_rev_fastq}.tmp"  | tee >(cat >&2)
					printf "OUTPUT3: %s\\n" "$REPORT_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "OUTPUT4: %s\\n" "{output.processed_fwd_fastq}"  | tee >(cat >&2)
					printf "OUTPUT5: %s\\n" "{output.processed_rev_fastq}"  | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "###################################- EXECUTION -##############################" | tee >(cat >&2)
					printf "%s\\n" "module load cutadapt/2.3 || exit 1"  | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "" | tee >(cat >&2)
					#
					##
					start_time="$(date -u +%s)"

					module load cutadapt/2.3 || exit 1

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "cutadapt {config_pre_process_Dict[CUTADAPT_PAIRED]} --cores={threads} --output={output.processed_fwd_fastq}.tmp \\
					--paired-output={output.processed_rev_fastq}.tmp {each_fastq} {each_fastq_reverse} > $REPORT_PATH/{each_fastq_begining}.cutadapt.txt" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "" | tee >(cat >&2)
					#
					##
					start_time="$(date -u +%s)"

					cutadapt {config_pre_process_Dict[CUTADAPT_PAIRED]} --cores={threads} --output={output.processed_fwd_fastq}.tmp --paired-output={output.processed_rev_fastq}.tmp \\
					{each_fastq} {each_fastq_reverse} > $REPORT_PATH/{each_fastq_begining}.cutadapt.txt

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "cat {output.processed_fwd_fastq}.tmp >> {output.processed_fwd_fastq}" | tee >(cat >&2)
					printf "%s\\n" "cat {output.processed_rev_fastq}.tmp >> {output.processed_rev_fastq}" | tee >(cat >&2)
					printf "%s\\n" "rm {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}.tmp" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "" | tee >(cat >&2)
					#
					##
					start_time="$(date -u +%s)"

					cat {output.processed_fwd_fastq}.tmp >> {output.processed_fwd_fastq}
					cat {output.processed_rev_fastq}.tmp >> {output.processed_rev_fastq}
					rm {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}.tmp

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "###################################- FINALIZING -#############################" | tee >(cat >&2)
					
					total_end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "%s\\n" "#" | tee >(cat >&2)
					printf "TOTAL ELAPSED TIME: %s seconds\\n" "$(($total_end_time-$total_start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				""")

			else:
				pass
			shell("""
				#
				##

				total_start_time="$(date -u +%s)"

				RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/pre_process
				mkdir -p $RESULT_PATH

				REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/pre_process
				mkdir -p $REPORT_PATH
				##
				#
				printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
				printf "%s\\n" "merging_preprocessed_paired: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
				printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{output.processed_fwd_fastq}"  | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{output.processed_rev_fastq}"  | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.processed_fwd_fastq}"  | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.processed_rev_fastq}"  | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$REPORT_PATH/{wildcards.sample}.R1_fastqc.html"  | tee >(cat >&2)
				printf "OUTPUT4: %s\\n" "$REPORT_PATH/{wildcards.sample}.R2_fastqc.html"  | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "module load usearch/11.0.667 || exit 1"  | tee >(cat >&2)
				printf "%s\\n" "module load fastqc/0.11.8 || exit 1"  | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				module load usearch/11.0.667 || exit 1
				module load fastqc/0.11.8 || exit 1

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "\\t%s\\n" "" | tee >(cat >&2)
				printf "%s\\n" "if [ {USE_SUBSAMPLE} = "TRUE" ]; then" | tee >(cat >&2)
				printf "\\t%s\\n" "usearch -fastx_subsample {output.processed_fwd_fastq} -reverse {output.processed_rev_fastq} \\
				-fastqout {output.processed_fwd_fastq}.tmp -output2 {output.processed_rev_fastq}.tmp -sample_pct '{SUBSAMPLE_PCT}' " | tee >(cat >&2)
				printf "\\t%s\\n" "mv {output.processed_rev_fastq}.tmp {output.processed_fwd_fastq}" | tee >(cat >&2)
				printf "\\t%s\\n" "mv {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}" | tee >(cat >&2)
				printf "%s\\n" "fi" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				if [ {USE_SUBSAMPLE} = "TRUE" ]; then

					usearch -fastx_subsample {output.processed_fwd_fastq} -reverse {output.processed_rev_fastq} -fastqout {output.processed_fwd_fastq}.tmp -output2 {output.processed_rev_fastq}.tmp -sample_pct '{SUBSAMPLE_PCT}'
					mv {output.processed_rev_fastq}.tmp {output.processed_fwd_fastq}
					mv {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}
					
				fi

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $REPORT_PATH -f fastq --threads {threads} {output.processed_fwd_fastq}" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $REPORT_PATH -f fastq --threads {threads} {output.processed_rev_fastq}" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				fastqc -o $REPORT_PATH -f fastq --threads {threads} {output.processed_fwd_fastq}
				fastqc -o $REPORT_PATH -f fastq --threads {threads} {output.processed_rev_fastq}

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "DONE!!!!"  | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "###################################- FINALIZING -#############################" | tee >(cat >&2)
				
				total_end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "TOTAL ELAPSED TIME: %s seconds\\n" "$(($total_end_time-$total_start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			""")

elif LAYOUT == "single":
	#
	rule Pre_Process_Single:
		input:
			fastq_List = get_fastq
		output:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
		threads: PROCESSORS
		message: "PreProcess_Single: {wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + SAMPLE_SUFFIX, "", each_fastq_basename)
				shell("""
					module load cutadapt/2.1
					module load fastqc/0.11.8
					QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $QC_PATH
					#
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "cutadapt/2.1"  | tee >(cat >&2)
					printf "%s\\n" "fastqc/0.11.8"  | tee >(cat >&2)
					printf "%s\\n" "Removing low complexity, remnant adapter and short reads."  | tee >(cat >&2)
					printf "INPUT1: %s\\n" "{each_fastq}"  | tee >(cat >&2)
					printf "OUTPUT1: %s\\n" "{output.processed_fastq}"  | tee >(cat >&2)
					printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "cutadapt {config_pre_process_Dict[CUTADAPT_SINGLE]} --cores={threads} {each_fastq} >> {output.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "%s\\n" "fastqc -o $QC_PATH -f fastq --threads {threads} {each_fastq}" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					cutadapt {config_pre_process_Dict[CUTADAPT_SINGLE]} --cores={threads} {each_fastq} >> {output.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.cutadapt.txt
					fastqc -o $QC_PATH -f fastq --threads {threads} {each_fastq}
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)

				""")
