# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
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
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_overlap_List
# ################################### PIPELINE RULES ########################## #


rule NarrowPeak_Overlap:
	"""
	"""
	input:
		narrowPeak_List = get_narrowpeak,
		pooled_narrowPeak = get_pooled_narrowpeak,
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "NarrowPeak_Overlap: {wildcards.design}|{wildcards.overlapped}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{FS="\\t";OFS="\\t"}} {{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}'
			AWK_COMMAND2='BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Overlap Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.narrowPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowPeak}) -b <(zcat -f {input.narrowPeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND2' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			intersectBed -wo -f 0.5 -a <(zcat -f {input.pooled_narrowPeak}) -b <(zcat -f {input.narrowPeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bdg}.tmp
			cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq
			slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.edited
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf {output.narrowPeak_bed}.fix
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")


rule NarrowPeak_Controlled_Overlap:
	"""
	"""
	input:
		narrowPeak_List = get_narrowpeak_controlled,
		pooled_narrowPeak = get_pooled_controlled_narrowpeak,
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "NarrowPeak_Controlled_Overlap: {wildcards.design}|{wildcards.overlapped}|{wildcards.control}"
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{FS="\\t";OFS="\\t"}} {{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.1) || ($21/s2 >= 0.1)) {{print $0}}}}'
			AWK_COMMAND2='BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Overlap Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.narrowPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_narrowPeak}) -b <(zcat -f {input.narrowPeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND2' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			intersectBed -wo -f 0.5 -a <(zcat -f {input.pooled_narrowPeak}) -b <(zcat -f {input.narrowPeak_List}) | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bdg}.tmp
			cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq
			slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.edited
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf {output.narrowPeak_bed}.fix
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")

rule NarrowPeak_IDR:
	"""
	"""
	input:
		narrowPeak_List = get_narrowpeak,
		pooled_narrowPeak = get_pooled_narrowpeak,
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "NarrowPeak_IDR: {wildcards.design}|{wildcards.IDR}"
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak IDR Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.narrowPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.1 --plot --use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowPeak_bed}.tmp --max-iter 10000 --verbose" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "mv {output.narrowPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_IDR.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND1' {output.narrowPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.narrowPeak_bed}.edited > {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			declare -a input_List=({input.narrowPeak_List})
			idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_narrowPeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowPeak_bed}.tmp --max-iter 10000 --log-output-file $QC_PATH/{wildcards.IDR}_IDR.txt
			mv {output.narrowPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_IDR.png
			#
			awk 'BEGIN{{OFS="\\t"}} $12>=1.30 {{print $0}}' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.filt
			if [ -s {output.narrowPeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited
			fi
			
			bgzip -c {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowPeak_bed}.edited > {output.narrowPeak_bdg}.tmp
			cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq
			slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.edited
				rm -rf {output.narrowPeak_bed}.fix
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
			
		""")


rule NarrowPeak_Controlled_IDR:
	"""
	"""
	input:
		narrowPeak_List = get_narrowpeak_controlled,
		pooled_narrowPeak = get_pooled_controlled_narrowpeak,
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "NarrowPeak_Controlled_IDR: {wildcards.design}|{wildcards.IDR}|{wildcards.control}"
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak IDR Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.narrowPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_narrowPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.1 --plot --use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowPeak_bed}.tmp --max-iter 10000 --verbose" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "mv {output.narrowPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND1' {output.narrowPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.narrowPeak_bed}.edited > {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			declare -a input_List=({input.narrowPeak_List})
			idr --samples ${{input_List[@]:0:2}}  --peak-list {input.pooled_narrowPeak} --input-file-type narrowPeak --output-file-type narrowPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --input-file-type narrowPeak --peak-merge-method sum --output-file {output.narrowPeak_bed}.tmp --max-iter 10000 --log-output-file $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.txt
			mv {output.narrowPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png

			awk 'BEGIN{{OFS="\\t"}} $12>=1.30 {{print $0}}' {output.narrowPeak_bed}.tmp > {output.narrowPeak_bed}.filt
			if [ -s {output.narrowPeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_narrowPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.narrowPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.narrowPeak_bed}.edited
			fi
			bgzip -c {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.edited > {output.narrowPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.narrowPeak_bed}.edited > {output.narrowPeak_bdg}.tmp
			cat {output.narrowPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.narrowPeak_bdg}.uniq
			slopBed -i {output.narrowPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.edited > {output.narrowPeak_bdg}
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.edited
				rm -rf {output.narrowPeak_bed}.fix
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
			
		""")


rule BroadPeak_Overlap:
	"""
	"""
	input:
		broadPeak_List = get_broadpeak,
		pooled_broadPeak = get_pooled_broadpeak,
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_VS_.*|.*_IDR_.*|\\.).)*}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "BroadPeak_Overlap: {wildcards.design}|{wildcards.overlapped}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{FS="\\t";OFS="\\t"}} {{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}'
			AWK_COMMAND2='BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Overlap Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.broadPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadPeak}) -b <(zcat -f {input.broadPeak_List}) | cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND2' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bed}.edited > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.broadPeak_bed}.sorted > {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			intersectBed -wo -f 0.5 -a <(zcat -f {input.pooled_broadPeak}) -b <(zcat -f {input.broadPeak_List}) | cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bed}.edited > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bdg}.tmp
			cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq
			slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.edited
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf {output.broadPeak_bed}.fix
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")


rule BroadPeak_Controlled_Overlap:
	"""
	"""
	input:
		broadPeak_List = get_broadpeak_controlled,
		pooled_broadPeak = get_pooled_controlled_broadpeak,
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped, ((?!.*_IDR_.*|\\.).)*}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "BroadPeak_Controlled_Overlap: {wildcards.design}|{wildcards.overlapped}|{wildcards.control}"
	resources:
		mem_mb = MEMORY
	
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/narrowpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{FS="\\t";OFS="\\t"}} {{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}'
			AWK_COMMAND2='BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Overlap Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.broadPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "intersectBed -wo -f 1E-9 -a <(zcat -f {input.pooled_broadPeak}) -b <(zcat -f {input.broadPeak_List}) | cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND2' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bed}.edited > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.broadPeak_bed}.sorted > {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			intersectBed -wo -f 0.5 -a <(zcat -f {input.pooled_broadPeak}) -b <(zcat -f {input.broadPeak_List}) | cut -f 1-9 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.tmp
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{$4="{wildcards.overlapped}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bed}.edited > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bdg}.tmp
			cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq
			slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.edited
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf {output.broadPeak_bed}.fix
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")


rule BroadPeak_IDR:
	"""
	"""
	input:
		broadPeak_List = get_broadpeak,
		pooled_broadPeak = get_pooled_broadpeak,
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "BroadPeak_IDR: {wildcards.design}|{wildcards.IDR}"
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak IDR Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.broadPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "idr --samples {input.broadPeak_List} --peak-list {input.pooled_broadPeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.1 --plot --use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadPeak_bed}.tmp --max-iter 10000 --verbose" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "mv {output.broadPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_IDR.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND1' {output.broadPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.edited > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.broadPeak_bed}.edited > {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.broadPeak_bed}.edited > {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			declare -a input_List=({input.broadPeak_List})
			idr --samples ${{input_List[@]:0:2}} --peak-list {input.pooled_broadPeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadPeak_bed}.tmp --max-iter 10000 --log-output-file $QC_PATH/{wildcards.IDR}_IDR.txt
			mv {output.broadPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_IDR.png

			awk 'BEGIN{{OFS="\\t"}} $12>=1.30 {{print $0}}' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.filt
			if [ -s {output.broadPeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited
			fi

			bgzip -c {output.broadPeak_bed}.edited > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.edited > {output.broadPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadPeak_bed}.edited > {output.broadPeak_bdg}.tmp
			cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq
			slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.edited
				rm -rf {output.broadPeak_bed}.fix
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")


rule BroadPeak_Controlled_IDR:
	"""
	"""
	input:
		broadPeak_List = get_broadpeak_controlled,
		pooled_broadPeak = get_pooled_controlled_broadpeak,
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR, ((?!.*_OVERLAPPED_.*|\\.).)*}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "BroadPeak_IDR: {wildcards.design}|{wildcards.IDR}|{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load idr/2.0.3 || exit 1
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling/broadpeak
			mkdir -p $QC_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			AWK_COMMAND1='BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'
			AWK_COMMAND3='BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			AWK_COMMAND4='BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load idr/2.0.3" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak IDR Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.broadPeak_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			index=$(($index+1))
			printf "INPUT%s: %s\\n" "${{index}}" "{input.pooled_broadPeak}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "idr --samples {input.broadPeak_List} --peak-list {input.pooled_broadPeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.1 --plot --use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadPeak_bed}.tmp --max-iter 10000 --verbose" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "mv {output.broadPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND1' {output.broadPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.edited > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND3' {output.broadPeak_bed}.edited > {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND4' {output.broadPeak_bed}.edited > {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.edited" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bdg}.edited" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			declare -a input_List=({input.broadPeak_List})
			idr --samples ${{input_List[@]:0:2}}  --peak-list {input.pooled_broadPeak} --input-file-type broadPeak --output-file-type broadPeak --rank signal.value --soft-idr-threshold 0.05 --plot \
			--use-best-multisummit-IDR --input-file-type broadPeak --peak-merge-method sum --output-file {output.broadPeak_bed}.tmp --max-iter 10000 --log-output-file $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.txt
			mv {output.broadPeak_bed}.tmp.png $QC_PATH/{wildcards.IDR}_VS_{wildcards.control}_IDR.png

			awk 'BEGIN{{OFS="\\t"}} $12>=1.30 {{print $0}}' {output.broadPeak_bed}.tmp > {output.broadPeak_bed}.filt
			if [ -s {output.broadPeak_bed}.filt ]; then
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.filt | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited
			else
				awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.IDR}_VS_{wildcards.control}.macs2_broadPeak_"NR}} {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' {output.broadPeak_bed}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.broadPeak_bed}.edited
			fi

			bgzip -c {output.broadPeak_bed}.edited > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{FS="\\t";OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.edited > {output.broadPeak_bed}.fix
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.fix {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n", $1,$2,$3,$5}}' {output.broadPeak_bed}.edited > {output.broadPeak_bdg}.tmp
			cat {output.broadPeak_bdg}.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.broadPeak_bdg}.uniq
			slopBed -i {output.broadPeak_bdg}.uniq -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.edited
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.edited > {output.broadPeak_bdg}
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.edited
				rm -rf {output.broadPeak_bed}.fix
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.edited
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
			
		""")
