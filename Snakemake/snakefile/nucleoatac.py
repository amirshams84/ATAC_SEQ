# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Annotate
# snakemake --snakefile nucleoatac.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile nucleoatac.py --configfile Encode.json --rulegraph | dot -Tsvg > nucleoatac.svg
# ################################### IMPORT ##################################### #


import os
import sys
import pandas
library_path = os.path.abspath(workflow.basedir + "/../library/")
sys.path.append(library_path)
import utility
# ################################### FUNCTIONS ################################## #


def build_design_Dict(metadata_Dict):
	#
	design_Dict = {}
	for sample, sample_Dict in metadata_Dict.items():
		#
		if sample_Dict["Design"] not in design_Dict:
			#
			design_Dict[sample_Dict["Design"]] = {}
			design_Dict[sample_Dict["Design"]]["Case"] = []
			design_Dict[sample_Dict["Design"]]["Control"] = []
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
		elif sample_Dict["Design"] in design_Dict:
			#
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
	return design_Dict
 ################################### CONFIGURATION ############################## #


# ++++++++++++++++++++++++++++++++++++
#PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../bash_script")
R_Script_path = os.path.abspath(workflow.basedir + "/../R_script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../python_script")
Script_Path = os.path.abspath(workflow.basedir + "/../template")
# -----------------------------------
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
LAYOUT = config_data_Dict["LAYOUT"].lower()
GENOME = config_data_Dict["GENOME"].lower()
# -----------------------------------
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
#UTILITIES
config_utilities_Dict = config["UTILITIES"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PEAK_CALLING
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
MACS2_NARROW_PARAMETERS = config_peak_calling_Dict["MACS2_NARROW"]
MACS2_BROAD_PARAMETERS = config_peak_calling_Dict["MACS2_BROAD"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]
EFFECTIVE_GENOME_SIZE = config_reference_Dict["EFFECTIVE_GENOME_SIZE"]
# ------------------------------------
# ################################### WILDCARDS ################################ #


post_alignment_List = []
broadpeak_List = []
nucleoatac_List = []
for sample, sample_Dict in sample_treatment_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	#
	broadpeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	nucleoatac_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucleoatac_signal.bedgraph.gz".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
for design in metadata_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	#
	broadpeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	nucleoatac_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{pooled_case}.nucleoatac_signal.bedgraph.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		nucleoatac_List
# ################################### PIPELINE RULES ########################## #


rule nucleoatac:
	input:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
	output:
		nucleoatac_signal_bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucleoatac_signal.bedgraph.gz",
		nucleoatac_occ_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.occpeaks.bed.gz",
		nucleoatac_nucpos_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucpos.bed.gz",
		nucleoatac_nfrpos_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nfrpos.bed.gz",
	priority: 996
	message: "nucleoatac: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = 200000
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/nucleoatac
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/nucleoatac
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "nucleoatac: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.broadpeak_bed}" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.nucleoatac_signal_bedgraph}" | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.nucleoatac_occ_bed}" | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.nucleoatac_nucpos_bed}" | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.nucleoatac_nfrpos_bed}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load nucleoatac/0.3.4 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load ucsc/373 || exit 1
			module load nucleoatac/0.3.4 || exit 1

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
			printf "%s\\n" "cd $SCRATCH_PATH" | tee >(cat >&2)
			printf "%s\\n" "nucleoatac run --bed {input.broadpeak_bed} --cores {threads} --bam {input.processed_bam} \\
			--fasta {config_reference_Dict[WG_FASTA]} --out {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/nucleoatac/{wildcards.sample}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			cd $SCRATCH_PATH
			
			nucleoatac run --bed {input.broadpeak_bed} --cores {threads} --bam {input.processed_bam} \\
			--fasta {config_reference_Dict[WG_FASTA]} --out {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/nucleoatac/{wildcards.sample}

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
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "TOTAL ELAPSED TIME: %s seconds\\n" "$(($total_end_time-$total_start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")
