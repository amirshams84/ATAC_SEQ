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
R_Script_path = os.path.abspath(workflow.basedir + "/../R_Script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../Python_Script")
Script_Path = os.path.abspath(workflow.basedir + "/../Script")
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
LAYOUT = config_data_Dict["LAYOUT"].lower()
GENOME = config_data_Dict["GENOME"].lower()
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
config_cluster_Dict = config["CLSUTER_CONFIG"]
PROCESSORS = config_cluster_Dict["PROCESSORS"]
MEMORY = config_cluster_Dict["MEMORY"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
metadata_Dict = utility.build_metadata_dict(METADATA_FILE)
design_Dict = build_design_Dict(metadata_Dict)
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
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	#
	broadpeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	nucleoatac_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucleoatac_signal.bedgraph.gz".format(design=sample_Dict["Design"], sample=sample))
for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=sample_Dict["Design"], pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#
	broadpeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=sample_Dict["Design"], pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	nucleoatac_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{pooled_case}.nucleoatac_signal.bedgraph.gz".format(design=sample_Dict["Design"], pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		nucleoatac_List
# ################################### PIPELINE RULES ########################## #


rule nucleoATAC:
	input:
		broadPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
	output:
		nucleoatac_Bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucleoatac_signal.bedgraph.gz",
		occ_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.occpeaks.bed.gz",
		nucpos_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nucpos.bed.gz",
		nfrpos_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/nucleoatac/{sample}.nfrpos.bed.gz",
	priority: 996
	message: "nucleoATAC: {wildcards.design}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load nucleoatac || exit 1
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load nucleoatac || exit 1" | tee >(cat >&2)
			printf "%s\\n" "nucleoatac run --bed {input.broadPeak_Bed} --cores {threads} --bam {input.processed_bam} --fasta {config_reference_Dict[WG_FASTA]} --out {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/nucleoatac/{wildcards.sample}"  | tee >(cat >&2)

			nucleoatac run --bed {input.broadPeak_Bed} --cores {threads} --bam {input.processed_bam} --fasta {config_reference_Dict[WG_FASTA]} --out {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/nucleoatac/{wildcards.sample}

		""")
