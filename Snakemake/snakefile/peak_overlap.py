# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Project: ENCODE ATAC_SEQ
# Aim: Snakemake workflow for Peak Analysis
# snakemake --snakefile peak_overlap.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile peak_overlap.py --configfile Yoko.json --rulegraph | dot -Tsvg > peak_analysis.svg
# ################################### IMPORT ##################################### #


import os
import sys
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


def get_narrowpeak(wildcards):
	"""
	"""
	narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}.narrowPeak.gz".format(design=wildcards.design, case=sample))
	return narrowPeak_List


def get_broadpeak(wildcards):
	"""
	"""
	broadPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}.broadPeak.gz".format(design=wildcards.design, case=sample))
	return broadPeak_List


def get_pooled_narrowpeak(wildcards):
	"""
	"""
	pooled_narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
		else:
			pass
	pooled_narrowPeak_List = list(set(pooled_narrowPeak_List))
	return pooled_narrowPeak_List


def get_pooled_bed(wildcards):
	"""
	"""
	pooled_bed_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_bed_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bed.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
		else:
			pass
	pooled_bed_List = list(set(pooled_bed_List))
	return pooled_bed_List


def get_pooled_broadpeak(wildcards):
	"""
	"""
	pooled_broadPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
		else:
			pass
	pooled_broadPeak_List = list(set(pooled_broadPeak_List))
	return pooled_broadPeak_List


def get_narrowpeak_controlled(wildcards):
	"""
	"""
	narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=wildcards.design, case=sample, control=wildcards.control))
	return narrowPeak_List


def get_broadpeak_controlled(wildcards):
	"""
	"""
	broadPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=wildcards.design, case=sample, control=wildcards.control))
	return broadPeak_List


def get_pooled_controlled_narrowpeak(wildcards):
	"""
	"""
	pooled_narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"]), control=wildcards.control))
	pooled_narrowPeak_List = list(set(pooled_narrowPeak_List))
	return pooled_narrowPeak_List


def get_pooled_controlled_broadpeak(wildcards):
	"""
	"""
	pooled_broadPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"]), control=wildcards.control))
	pooled_broadPeak_List = list(set(pooled_broadPeak_List))
	return pooled_broadPeak_List

# ################################### CONFIGURATION ############################## #

# +++++++++++++++++++++++++++++++++++++
#PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../Bash_Script")
R_Script_path = os.path.abspath(workflow.basedir + "/../R_Script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../Python_Script")
Script_Path = os.path.abspath(workflow.basedir + "/../Script")
# ------------------------------------
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
# ################################### WILDCARDS ################################ #


post_alignment_List = []
peak_calling_List = []
peak_overlap_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))


for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	##
	#
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
			#pass
	##		
	for control in design_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}_VS_{control}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}_VS_{control}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}_VS_{control}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}_VS_{control}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		#pass
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_overlap_List
# ################################### PIPELINE RULES ########################## #


rule narrowpeak_overlap:
	"""
	"""
	input:
		narrowpeak_List = get_narrowpeak,
		pooled_narrowpeak = get_pooled_narrowpeak,
		pooled_bed = get_pooled_bed,
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "narrowpeak_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.overlapped}.macs2_narrowPeak_"NR}} {{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf "%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n", \$1,\$2,\$3,\$5}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "narrowpeak_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.narrowpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/{wildcards.overlapped}.narrowPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1

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
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowpeak}) -b <(zcat -f {input.narrowpeak_List}) | \\
			cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowpeak}) -b <(zcat -f {input.narrowpeak_List}) | \\
			cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.sorted
			bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}
			tabix -f -p bed {output.narrowpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.narrowpeak_bed}.sorted > {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bdg}.tmp
			cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq
			slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}
			bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh {wildcards.overlapped} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/{wildcards.overlapped}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh {wildcards.overlapped} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/{wildcards.overlapped}.narrowPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				rm -rf {output.narrowpeak_bed}.tmp
				rm -rf {output.narrowpeak_bed}.edited
				rm -rf {output.narrowpeak_bed}.sorted
				rm -rf {output.narrowpeak_bed}.fix
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.edited
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

rule narrowpeak_IDR:
	"""
	"""
	input:
		narrowpeak_List = get_narrowpeak,
		pooled_narrowpeak = get_pooled_narrowpeak,
		pooled_bed = get_pooled_bed,
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "narrowpeak_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}"
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")

			declare -a input_List=({input.narrowpeak_List})


			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if (\$2<0) \$2=0; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf \\\"%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n\\\", \$1,\$2,\$3,\$5}}\'"
			AWK_COMMAND5="awk -v p=0.05 \'BEGIN{{print -log(p)/log(10)}}\'"
			AWK_COMMAND6="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} \$12>=\$idr_thresh_transformed {{print \$0}}\'"

			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "narrowpeak_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.narrowpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/{wildcards.IDR}.narrowPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1

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
			printf "%s\\n" "declare -a input_List=({input.narrowpeak_List})" | tee >(cat >&2)
			printf "%s\\n" "idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_narrowpeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}.narrowPeak.txt" | tee >(cat >&2)
			printf "%s\\n" "mv {output.narrowpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}.narrowPeak.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			declare -a input_List=({input.narrowpeak_List})
			idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_narrowpeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}.narrowPeaks.txt
			mv {output.narrowpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}.narrowPeaks.png

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
			printf "%s\\n" "idr_thresh_transformed=\$($AWK_COMMAND5)" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND6 {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.filt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s {output.narrowpeak_bed}.filt ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")
			
			awk 'BEGIN{{OFS="\\t"}} $12>='"${{idr_thresh_transformed}}"' {{print $0}}' {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.filt

			if [ -s {output.narrowpeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.filt | \\
				sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.tmp | \\
				sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited
			fi

			bgzip -c {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}
			tabix -f -p bed {output.narrowpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.narrowpeak_bed}.edited > {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowpeak_bed}.edited > {output.narrowpeak_bdg}.tmp
			cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq
			slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}
			bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh {wildcards.IDR} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/{wildcards.IDR}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh {wildcards.IDR} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/{wildcards.IDR}.narrowPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "mv {output.narrowpeak_bed}.tmp {output.narrowpeak_bed}.raw" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.filt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				mv {output.narrowpeak_bed}.tmp {output.narrowpeak_bed}.raw
				rm -rf {output.narrowpeak_bed}.filt
				rm -rf {output.narrowpeak_bed}.edited
				rm -rf {output.narrowpeak_bed}.fix
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.edited
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


rule narrowpeak_controlled_overlap:
	"""
	"""
	input:
		narrowpeak_List = get_narrowpeak_controlled,
		pooled_narrowpeak = get_pooled_controlled_narrowpeak,
		pooled_bed = get_pooled_bed,
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "narrowpeak_controlled_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}_VS_{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.overlapped}.macs2_narrowPeak_"NR}} {{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf "%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n", \$1,\$2,\$3,\$5}}\'"

			sample_Name=$(basename {output.narrowpeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "narrowpeak_controlled_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|${{sample_Name}}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.narrowpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1

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
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowpeak}) -b <(zcat -f {input.narrowpeak_List}) | \\
			cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			sample_Name=$(basename {output.narrowpeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}

			intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowpeak}) -b <(zcat -f {input.narrowpeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.sorted
			bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}
			tabix -f -p bed {output.narrowpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.narrowpeak_bed}.sorted > {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bdg}.tmp
			cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq
			slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}
			bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} \\
			{output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				rm -rf {output.narrowpeak_bed}.tmp
				rm -rf {output.narrowpeak_bed}.edited
				rm -rf {output.narrowpeak_bed}.sorted
				rm -rf {output.narrowpeak_bed}.fix
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.edited
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


rule narrowpeak_controlled_IDR:
	"""
	"""
	input:
		narrowpeak_List = get_narrowpeak_controlled,
		pooled_narrowpeak = get_pooled_controlled_narrowpeak,
		pooled_bed = get_pooled_bed,
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "narrowpeak_controlled_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}_VS_{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")

			sample_Name=$(basename {output.narrowpeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if (\$2<0) \$2=0; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf \\\"%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n\\\", \$1,\$2,\$3,\$5}}\'"
			AWK_COMMAND5="awk -v p=0.05 \'BEGIN{{print -log(p)/log(10)}}\'"
			AWK_COMMAND6="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} \$12>=\$idr_thresh_transformed {{print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "narrowpeak_controlled_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|${{sample_Name}}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.narrowpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1

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
			printf "%s\\n" "declare -a input_List=({input.narrowpeak_List})" | tee >(cat >&2)
			printf "%s\\n" "idr --samples \${{input_List[@]:0:2}} --peak-list {input.pooled_narrowpeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}.narrowPeak.txt" | tee >(cat >&2)
			printf "%s\\n" "mv {output.narrowpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}.narrowPeak.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			declare -a input_List=({input.narrowpeak_List})
			idr --samples ${{input_List[@]:0:2}}  --peak-list {input.pooled_narrowpeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.txt
			mv {output.narrowpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png

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
			printf "%s\\n" "idr_thresh_transformed=\$($AWK_COMMAND5)" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND6 {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.filt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s {output.narrowpeak_bed}.filt ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.narrowpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")
			
			awk 'BEGIN{{OFS="\\t"}} $12>='"${{idr_thresh_transformed}}"' {{print $0}}' {output.narrowpeak_bed}.tmp > {output.narrowpeak_bed}.filt

			if [ -s {output.narrowpeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowpeak_bed}.edited
			fi

			bgzip -c {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}
			tabix -f -p bed {output.narrowpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.edited > {output.narrowpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.narrowpeak_bed}.edited > {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowpeak_bed}.edited > {output.narrowpeak_bdg}.tmp
			cat {output.narrowpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowpeak_bdg}.uniq
			slopBed -i {output.narrowpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.edited > {output.narrowpeak_bdg}
			bedGraphToBigWig {output.narrowpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/{wildcards.IDR}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "mv {output.narrowpeak_bed}.tmp {output.narrowpeak_bed}.raw" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.filt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				mv {output.narrowpeak_bed}.tmp {output.narrowpeak_bed}.raw
				rm -rf {output.narrowpeak_bed}.filt
				rm -rf {output.narrowpeak_bed}.edited
				rm -rf {output.narrowpeak_bed}.fix
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.edited
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


rule broadpeak_overlap:
	"""
	"""
	input:
		broadpeak_List = get_broadpeak,
		pooled_broadpeak = get_pooled_broadpeak,
		pooled_bed = get_pooled_bed,
	output:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "broadpeak_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.overlapped}.macs2_broadPeak_"NR}} {{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf "%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n", \$1,\$2,\$3,\$5}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "broadpeak_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.broadpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/{wildcards.overlapped}.broadPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1

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
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadpeak}) -b <(zcat -f {input.broadpeak_List}) | \\
			cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadpeak_bed}.sorted > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadpeak}) -b <(zcat -f {input.broadpeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}.sorted
			bgzip -c {output.broadpeak_bed}.sorted > {output.broadpeak_bed}
			tabix -f -p bed {output.broadpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.broadpeak_bed}.sorted > {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bdg}.tmp
			cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq
			slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}
			bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh {wildcards.overlapped} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/{wildcards.overlapped}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh {wildcards.overlapped} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/{wildcards.overlapped}.broadPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				rm -rf {output.broadpeak_bed}.tmp
				rm -rf {output.broadpeak_bed}.edited
				rm -rf {output.broadpeak_bed}.sorted
				rm -rf {output.broadpeak_bed}.fix
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.edited
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


rule broadpeak_IDR:
	"""
	"""
	input:
		broadpeak_List = get_broadpeak,
		pooled_broadpeak = get_pooled_broadpeak,
		pooled_bed = get_pooled_bed,
	output:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "broadpeak_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}"
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"
			
			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")

			declare -a input_List=({input.broadpeak_List})

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if (\$2<0) \$2=0; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf \\\"%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n\\\", \$1,\$2,\$3,\$5}}\'"
			AWK_COMMAND5="awk -v p=0.05 \'BEGIN{{print -log(p)/log(10)}}\'"
			AWK_COMMAND6="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} \$11>=\$idr_thresh_transformed {{print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "broadpeak_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.broadpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/{wildcards.IDR}.broadPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1

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
			printf "%s\\n" "declare -a input_List=({input.broadpeak_List})" | tee >(cat >&2)
			printf "%s\\n" "idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_broadpeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}.broadPeak.txt" | tee >(cat >&2)
			printf "%s\\n" "mv {output.broadpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}.broadPeak.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			declare -a input_List=({input.broadpeak_List})
			idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_broadpeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}_IDR.txt
			mv {output.broadpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}_IDR.png

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
			printf "%s\\n" "idr_thresh_transformed=\$($AWK_COMMAND5)" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND6 {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.filt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s {output.broadpeak_bed}.filt ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadpeak_bed}.edited > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")
			awk 'BEGIN{{OFS="\\t"}} $11>='"${{idr_thresh_transformed}}"' {{print $0}}' {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.filt

			if [ -s {output.broadpeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited
			fi

			bgzip -c {output.broadpeak_bed}.edited > {output.broadpeak_bed}
			tabix -f -p bed {output.broadpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.broadpeak_bed}.edited > {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.edited > {output.broadpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.broadpeak_bed}.edited > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadpeak_bed}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bed}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadpeak_bed}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadpeak_bed} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadpeak_bed}.edited > {output.broadpeak_bdg}.tmp
			cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq
			slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}
			bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh {wildcards.IDR} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/{wildcards.IDR}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh {wildcards.IDR} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/{wildcards.IDR}.broadPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "mv {output.broadpeak_bed}.tmp {output.broadpeak_bed}.raw" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.filt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				mv {output.broadpeak_bed}.tmp {output.broadpeak_bed}.raw
				rm -rf {output.broadpeak_bed}.filt
				rm -rf {output.broadpeak_bed}.edited
				rm -rf {output.broadpeak_bed}.sorted
				rm -rf {output.broadpeak_bed}.fix
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.edited
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


rule broadpeak_controlled_overlap:
	"""
	"""
	input:
		broadpeak_List = get_broadpeak_controlled,
		pooled_broadpeak = get_pooled_controlled_broadpeak,
		pooled_bed = get_pooled_bed,
	output:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "broadpeak_controlled_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.overlapped}_VS_{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi

			sample_Name=$(basename {output.broadpeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="{wildcards.overlapped}.macs2_broadPeak_"NR}} {{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf "%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n", \$1,\$2,\$3,\$5}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "broadpeak_controlled_overlap: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|${{sample_Name}}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.broadpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1

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
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadpeak}) -b <(zcat -f {input.broadpeak_List}) | \\
			cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadpeak_bed}.sorted > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			sample_Name=$(basename {output.broadpeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}

			intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadpeak}) -b <(zcat -f {input.broadpeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}.sorted
			bgzip -c {output.broadpeak_bed}.sorted > {output.broadpeak_bed}
			tabix -f -p bed {output.broadpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.broadpeak_bed}.sorted > {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bdg}.tmp
			cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq
			slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}
			bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				rm -rf {output.broadpeak_bed}.tmp
				rm -rf {output.broadpeak_bed}.edited
				rm -rf {output.broadpeak_bed}.sorted
				rm -rf {output.broadpeak_bed}.fix
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.edited
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


rule broadpeak_controlled_IDR:
	"""
	"""
	input:
		broadpeak_List = get_broadpeak_controlled,
		pooled_broadpeak = get_pooled_controlled_broadpeak,
		pooled_bed = get_pooled_bed,
	output:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "broadpeak_controlled_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.IDR}_VS_{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")

			sample_Name=$(basename {output.broadpeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}

			declare -a input_List=({input.broadpeak_List})

			AWK_COMMAND1="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= 0.5) || (\$21/s2 >= 0.5)) {{print \$0}}}}\'"
			AWK_COMMAND2="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{\$4="${{sample_Name}}.macs2_broadPeak_"NR}} {{if (\$2<0) \$2=0; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9}}\'"
			AWK_COMMAND3="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{if (\$5>1000) \$5=1000; print \$0}}\'"
			AWK_COMMAND4="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} {{printf \\\"%s\\\\t%d\\\\t%d\\\\t%2.3f\\\\n\\\", \$1,\$2,\$3,\$5}}\'"
			AWK_COMMAND5="awk -v p=0.05 \'BEGIN{{print -log(p)/log(10)}}\'"
			AWK_COMMAND6="awk \'BEGIN{{FS=\\\"\\\\t\\\";OFS=\\\"\\\\t\\\"}} \$11>=\$idr_thresh_transformed {{print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "broadpeak_controlled_IDR: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|${{sample_Name}}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.broadpeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadpeak}"  | tee >(cat >&2)
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_bed}"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadpeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadpeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadpeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadpeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadpeak_bigwig}"  | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt"  | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load idr/2.0.3 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			

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
			printf "%s\\n" "declare -a input_List=({input.broadpeak_List})" | tee >(cat >&2)
			printf "%s\\n" "idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_broadpeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \\
			--use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}.broadPeak.txt" | tee >(cat >&2)
			printf "%s\\n" "mv {output.broadpeak_bed}.tmp.png $REPORT_PATH/${{sample_Name}}.broadPeak.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			declare -a input_List=({input.broadpeak_List})
			idr --samples ${{input_List[@]:0:2}}  --peak-list {input.pooled_broadpeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --peak-merge-method sum --output-file {output.broadpeak_bed}.tmp --max-iter 10000 --log-output-file $REPORT_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.txt
			mv {output.broadpeak_bed}.tmp.png $REPORT_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png

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
			printf "%s\\n" "idr_thresh_transformed=\$($AWK_COMMAND5)" | tee >(cat >&2)
			printf "%s\\n" "$AWK_COMMAND6 {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.filt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s {output.broadpeak_bed}.filt ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "$AWK_COMMAND2 {output.broadpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadpeak_bed}.edited > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			idr_thresh_transformed=$(awk -v p=0.05 "BEGIN{{print -log(p)/log(10)}}")
			
			awk 'BEGIN{{OFS="\\t"}} $11>='"${{idr_thresh_transformed}}"' {{print $0}}' {output.broadpeak_bed}.tmp > {output.broadpeak_bed}.filt

			if [ -s {output.broadpeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="${{sample_Name}}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="${{sample_Name}}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadpeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadpeak_bed}.edited
			fi
			
			bgzip -c {output.broadpeak_bed}.edited > {output.broadpeak_bed}
			tabix -f -p bed {output.broadpeak_bed}

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
			printf "%s\\n" "$AWK_COMMAND3 {output.broadpeak_bed}.edited > {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.edited > {output.broadpeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "$AWK_COMMAND4 {output.broadpeak_bed}.edited > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadpeak_bed}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bed}.uniq" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadpeak_bed}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bed}.edited > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadpeak_bed} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadpeak_bed}.edited > {output.broadpeak_bdg}.tmp
			cat {output.broadpeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadpeak_bdg}.uniq
			slopBed -i {output.broadpeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.edited > {output.broadpeak_bdg}
			bedGraphToBigWig {output.broadpeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}

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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.pooled_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt

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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "mv {output.broadpeak_bed}.tmp {output.broadpeak_bed}.raw" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.filt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bigbed}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				mv {output.broadpeak_bed}.tmp {output.broadpeak_bed}.raw
				rm -rf {output.broadpeak_bed}.filt
				rm -rf {output.broadpeak_bed}.edited
				rm -rf {output.broadpeak_bed}.fix
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.edited
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
