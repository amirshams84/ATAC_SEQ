# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Calling
# snakemake --snakefile peak_calling.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile peak_calling.py --configfile Encode.json --rulegraph | dot -Tsvg > peak_calling.svg
# ################################### IMPORT ##################################### #


import os
import sys
import re
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
# ################################### CONFIGURATION ############################## #


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
overlap_peak_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bed.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##
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
		

# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_calling_List

# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++

rule Peak_Calling_Narrow:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "Peak_Calling_Narrow: {wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load samtools/1.9 || exit 1
			module load macs/2.1.2 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND_1='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Calling Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND_1' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.narrowPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.narrowPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.narrowPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then" | tee >(cat >&2)
			line="\\$line"
			printf "\\t%s\\n" "while IFS= read line; do grep -nr -m 1 --perl-regex '$line' {output.narrowPeak_bdg} | cut -d':' -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "while IFS= read line; do sed '${{line}}d' {output.narrowPeak_bdg} > {output.narrowPeak_bdg}.tmp; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bdg.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}
			#
			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			#
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE --pseudocount 0.00001
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}
			#
			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowPeak_bdg} > {output.narrowPeak_bdg}.uniq
			cut -f1,2,3 {output.narrowPeak_bdg}.uniq | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowPeak_bdg}.uniq | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowPeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			#
			if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				cat $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowPeak_bdg}.uniq | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done
				cat $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "${{line}}d" {output.narrowPeak_bdg}.uniq > {output.narrowPeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.uniq.fix > {output.narrowPeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.narrowPeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			else
				bedGraphToBigWig {output.narrowPeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			fi
			
			#
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.uniq.fix
				rm -rf {output.narrowPeak_bdg}.uniq.fix.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")

rule Peak_Calling_Narrow_Controlled:
	input:
		processed_case_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
		#
		processed_case_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz",
		processed_case_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz.tbi",
		processed_control_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz",
		processed_control_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz.tbi",
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "Peak_Calling_Narrow_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load samtools/1.9 || exit 1
			module load macs/2.1.2 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f {Script_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Script_Path}/bigNarrowPeak.as
			fi
			#
			AWK_COMMAND_1='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			sample_Name=$(basename {output.narrowPeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Calling Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND_1' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.narrowPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.narrowPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.narrowPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then" | tee >(cat >&2)
			line="\\$line"
			printf "\\t%s\\n" "while IFS= read line; do grep -nr -m 1 --perl-regex '$line' {output.narrowPeak_bdg} | cut -d':' -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "while IFS= read line; do sed '${{line}}d' {output.narrowPeak_bdg} > {output.narrowPeak_bdg}.tmp; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.narrowPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bdg.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			sample_Name=$(basename {output.narrowPeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}

			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bedToBigBed -as={Script_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}

			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE --pseudocount 0.00001
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}
			#
			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowPeak_bdg} > {output.narrowPeak_bdg}.uniq
			cut -f1,2,3 {output.narrowPeak_bdg}.uniq | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowPeak_bdg}.uniq | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowPeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			#
			if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				cat $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowPeak_bdg}.uniq | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done
				cat $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "${{line}}d" {output.narrowPeak_bdg}.uniq > {output.narrowPeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.uniq.fix > {output.narrowPeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.narrowPeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			else
				bedGraphToBigWig {output.narrowPeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			fi
			
			#
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf {output.narrowPeak_bdg}.tmp
				rm -rf {output.narrowPeak_bdg}.uniq
				rm -rf {output.narrowPeak_bdg}.uniq.fix
				rm -rf {output.narrowPeak_bdg}.uniq.fix.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")

rule Peak_Calling_Broad:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			#
			module load samtools/1.9 || exit 1
			module load macs/2.1.2 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			#
			AWK_COMMAND_1='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Calling Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND_1' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.broadPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.broadPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.broadPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then" | tee >(cat >&2)
			line="\\$line"
			printf "\\t%s\\n" "while IFS= read line; do grep -nr -m 1 --perl-regex '$line' {output.broadPeak_bdg} | cut -d':' -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "while IFS= read line; do sed '${{line}}d' {output.broadPeak_bdg} > {output.broadPeak_bdg}.tmp; done< $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bdg.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}

			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE --pseudocount 0.00001
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}
			#
			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadPeak_bdg} > {output.broadPeak_bdg}.uniq
			cut -f1,2,3 {output.broadPeak_bdg}.uniq | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadPeak_bdg}.uniq | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadPeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			#
			if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				cat $OUT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.broadPeak_bdg}.uniq | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done
				cat $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "${{line}}d" {output.broadPeak_bdg}.uniq > {output.broadPeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.uniq.fix > {output.broadPeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.broadPeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			else
				bedGraphToBigWig {output.broadPeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			fi
			
			#
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.uniq.fix
				rm -rf {output.broadPeak_bdg}.uniq.fix.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")

rule Peak_Calling_Broad_Controlled:
	input:
		processed_case_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
		#
		processed_case_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz",
		processed_case_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz.tbi",
		processed_control_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz",
		processed_control_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz.tbi",
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			#
			module load samtools/1.9 || exit 1
			module load macs/2.1.2 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			#
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f {Script_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Script_Path}/bigBroadPeak.as
			fi
			#
			#
			AWK_COMMAND_1='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			sample_Name=$(basename {output.broadPeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Calling Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND_1' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.broadPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.broadPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.broadPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then" | tee >(cat >&2)
			line="\\$line"
			printf "\\t%s\\n" "while IFS= read line; do grep -nr -m 1 --perl-regex '$line' {output.broadPeak_bdg} | cut -d':' -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "while IFS= read line; do sed '${{line}}d' {output.broadPeak_bdg} > {output.broadPeak_bdg}.tmp; done< $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadPeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bdg.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			sample_Name=$(basename {output.broadPeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}

			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bedToBigBed -as={Script_Path}/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}

			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE --pseudocount 0.00001
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}
			#
			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadPeak_bdg} > {output.broadPeak_bdg}.uniq
			cut -f1,2,3 {output.broadPeak_bdg}.uniq | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadPeak_bdg}.uniq | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadPeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			#
			if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				cat $OUT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.broadPeak_bdg}.uniq | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done
				cat $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "${{line}}d" {output.broadPeak_bdg}.uniq > {output.broadPeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.uniq.fix > {output.broadPeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.broadPeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			else
				bedGraphToBigWig {output.broadPeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			fi
			
			#
			if [ -f {output.broadPeak_bigwig} ]; then
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf {output.broadPeak_bdg}.tmp
				rm -rf {output.broadPeak_bdg}.uniq
				rm -rf {output.broadPeak_bdg}.uniq.fix
				rm -rf {output.broadPeak_bdg}.uniq.fix.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt
			fi
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")

