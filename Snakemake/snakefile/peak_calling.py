# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Project: ENCODE ATAC_SEQ
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

# ################################### CONFIGURATION ############################## #


# ++++++++++++++++++++++++++++++++++++
#PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../bash_script")
R_Script_path = os.path.abspath(workflow.basedir + "/../R_script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../python_script")
Template_Path = os.path.abspath(workflow.basedir + "/../template")
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
# ################################### WILDCARDS ################################ #


post_alignment_List = []
peak_calling_List = []
overlap_peak_List = []
for sample, sample_Dict in sample_treatment_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))

for design in metadata_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bed.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	##
	for case in metadata_Dict[design]["Case"]:
		for control in metadata_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
	##
	for control in metadata_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"]), control=control))
		

# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_calling_List

# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++

rule peak_calling_narrow:
	"""
	"""
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "peak_calling_narrow: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
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

			if [ ! -f {Template_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Template_Path}/bigNarrowPeak.as
			fi
			
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			AWK_BED_MAX_BED_SIGNAL_FILTER="awk \'BEGIN{{OFS=\\\"\\\\t\\\"}}{{if (\$5>1000) \$5=1000; print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "peak_calling_narrow: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
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
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
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
			module load macs/2.1.2 || exit 1
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
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $RESULT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $RESULT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowpeak_bed}.sorted
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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt

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
			printf "%s\\n" "$AWK_BED_MAX_BED_SIGNAL_FILTER {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Template_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.tmp
			bedToBigBed -as={Template_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_narrow --outdir $RESULT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.tmp > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_narrow --outdir $RESULT_PATH --method FE --pseudocount 0.00001
			slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.tmp > {output.narrowpeak_bdg}

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
			printf "%s\\n" "LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowpeak_bdg} > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.narrowpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d | sed 's/\\\\t/\\\\t.*\\\\t/' >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			line="\$line"
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "\$line" {output.narrowpeak_bdg}.uniq | \\
			cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done" | tee >(cat >&2)
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "\${{line}}d" {output.narrowpeak_bdg}.uniq > {output.narrowpeak_bdg}.uniq.fix; done" | tee >(cat >&2)
			printf "\\t%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.uniq.fix > {output.narrowpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowpeak_bdg} > {output.narrowpeak_bdg}.uniq
			cut -f1,2,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			#
			if [ -s $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowpeak_bdg}.uniq | cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done
				cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "${{line}}d" {output.narrowpeak_bdg}.uniq > {output.narrowpeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.uniq.fix > {output.narrowpeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.narrowpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}
			else
				bedGraphToBigWig {output.narrowpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}
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
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowpeak_bed}.tmp
				rm -rf {output.narrowpeak_bed}.sorted
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.uniq.fix
				rm -rf {output.narrowpeak_bdg}.uniq.fix.sorted
				rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
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

rule peak_calling_narrow_controlled:
	"""
	"""
	input:
		processed_case_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_case_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz",
		processed_case_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz.tbi",
		#
		processed_control_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
		processed_control_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz",
		processed_control_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz.tbi",
	output:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz",
		narrowpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz.tbi",
		narrowpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bb",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bdg",
		narrowpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	message: "peak_calling_narrow_controlled: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
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

			if [ ! -f {Template_Path}/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O {Template_Path}/bigNarrowPeak.as
			fi

			sample_Name=$(basename {output.narrowpeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}

			AWK_BED_MAX_BED_SIGNAL_FILTER="awk \'BEGIN{{OFS=\\\"\\\\t\\\"}}{{if (\$5>1000) \$5=1000; print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "peak_calling_narrow_controlled: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.case}_VS_{wildcards.control}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_case_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_case_bed_index}"  | tee >(cat >&2)
			printf "INPUT5: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT6: %s\\n" "{input.processed_control_bam_index}"  | tee >(cat >&2)
			printf "INPUT7: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "INPUT8: %s\\n" "{input.processed_control_bed_index}"  | tee >(cat >&2)
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
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
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
			module load macs/2.1.2 || exit 1
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
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $RESULT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $RESULT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowpeak_bed}.sorted
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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_case_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_case_bed} {output.narrowpeak_bed} $REPORT_PATH/${{sample_Name}}.narrowPeak.frip.txt

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
			printf "%s\\n" "$AWK_BED_MAX_BED_SIGNAL_FILTER {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Template_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowpeak_bed}.sorted > {output.narrowpeak_bed}.tmp
			bedToBigBed -as={Template_Path}/bigNarrowPeak.as -type=bed6+4 {output.narrowpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigbed}

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
			printf "%s\\n" "macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_narrow --outdir $RESULT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.tmp > {output.narrowpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $RESULT_PATH --method FE --pseudocount 0.00001
			slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.tmp > {output.narrowpeak_bdg}

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
			printf "%s\\n" "LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowpeak_bdg} > {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.narrowpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d | sed 's/\\\\t/\\\\t.*\\\\t/' >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			line="\$line"
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "\$line" {output.narrowpeak_bdg}.uniq | \\
			cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done" | tee >(cat >&2)
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "\${{line}}d" {output.narrowpeak_bdg}.uniq > {output.narrowpeak_bdg}.uniq.fix; done" | tee >(cat >&2)
			printf "\\t%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.uniq.fix > {output.narrowpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.narrowpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.narrowpeak_bdg} > {output.narrowpeak_bdg}.uniq
			cut -f1,2,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowpeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
			#

			if [ -s $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowpeak_bdg}.uniq | cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done
				cat $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt | while read line; do sed "${{line}}d" {output.narrowpeak_bdg}.uniq > {output.narrowpeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowpeak_bdg}.uniq.fix > {output.narrowpeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.narrowpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}
			else
				bedGraphToBigWig {output.narrowpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.narrowpeak_bigwig}
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
			printf "%s\\n" "if [ -f {output.narrowpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.narrowpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.narrowpeak_bigwig} ]; then
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowpeak_bed}.tmp
				rm -rf {output.narrowpeak_bed}.sorted
				rm -rf {output.narrowpeak_bdg}.tmp
				rm -rf {output.narrowpeak_bdg}.uniq
				rm -rf {output.narrowpeak_bdg}.uniq.fix
				rm -rf {output.narrowpeak_bdg}.uniq.fix.sorted
				rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $RESULT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
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

rule peak_calling_broad:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "peak_calling_broad: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
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

			if [ ! -f {Template_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Template_Path}/bigBroadPeak.as
			fi

			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			AWK_BED_MAX_BED_SIGNAL_FILTER="awk \'BEGIN{{OFS=\\\"\\\\t\\\"}}{{if (\$5>1000) \$5=1000; print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "peak_calling_broad: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
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
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
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
			module load macs/2.1.2 || exit 1
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
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $RESULT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadpeak_bed}.sorted > {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadpeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $RESULT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadpeak_bed}.sorted
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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt

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
			printf "%s\\n" "$AWK_BED_MAX_BED_SIGNAL_FILTER {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Template_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.tmp
			bedToBigBed -as={Template_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_broad --outdir $RESULT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.tmp > {output.broadpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_broad --outdir $RESULT_PATH --method FE --pseudocount 0.00001
			slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.tmp > {output.broadpeak_bdg}

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
			printf "%s\\n" "LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadpeak_bdg} > {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.broadpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.broadpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.broadpeak_bdg}.uniq | sort | uniq -d | sed 's/\\\\t/\\\\t.*\\\\t/' >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then" | tee >(cat >&2)
			line="\$line"
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "\$line" {output.broadpeak_bdg}.uniq | \\
			cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done" | tee >(cat >&2)
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "\${{line}}d" {output.broadpeak_bdg}.uniq > {output.broadpeak_bdg}.uniq.fix; done" | tee >(cat >&2)
			printf "\\t%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.uniq.fix > {output.broadpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.broadpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.broadpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadpeak_bdg} > {output.broadpeak_bdg}.uniq
			cut -f1,2,3 {output.broadpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadpeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			#
			if [ -s $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				cat $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.broadpeak_bdg}.uniq | cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done
				cat $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "${{line}}d" {output.broadpeak_bdg}.uniq > {output.broadpeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.uniq.fix > {output.broadpeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.broadpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}
			else
				bedGraphToBigWig {output.broadpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}
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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg
				rm -rf {output.broadpeak_bed}.tmp
				rm -rf {output.broadpeak_bed}.sorted
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.uniq.fix
				rm -rf {output.broadpeak_bdg}.uniq.fix.sorted
				rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt
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

rule peak_calling_broad_controlled:
	"""
	"""
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
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz",
		broadpeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz.tbi",
		broadpeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bb",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bdg",
		broadpeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "peak_calling_broad_controlled: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
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

			if [ ! -f {Template_Path}/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O {Template_Path}/bigBroadPeak.as
			fi

			sample_Name=$(basename {output.broadpeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}

			AWK_BED_MAX_BED_SIGNAL_FILTER="awk \'BEGIN{{OFS=\\\"\\\\t\\\"}}{{if (\$5>1000) \$5=1000; print \$0}}\'"
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "peak_calling_broad_controlled: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.case}_VS_{wildcards.control}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_case_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_case_bed_index}"  | tee >(cat >&2)
			printf "INPUT5: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT6: %s\\n" "{input.processed_control_bam_index}"  | tee >(cat >&2)
			printf "INPUT7: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "INPUT8: %s\\n" "{input.processed_control_bed_index}"  | tee >(cat >&2)
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
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
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
			module load macs/2.1.2 || exit 1
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
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $RESULT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadpeak_bed}.sorted" | tee >(cat >&2)
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

			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $RESULT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadpeak_bed}.sorted
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
			printf "%s\\n" "bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_case_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bash {Bash_Script}/frip_score.sh ${{sample_Name}} {input.processed_case_bed} {output.broadpeak_bed} $REPORT_PATH/${{sample_Name}}.broadPeak.frip.txt

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
			printf "%s\\n" "$AWK_BED_MAX_BED_SIGNAL_FILTER {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as={Template_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			awk 'BEGIN{{OFS=FS}} {{if ($5>1000) $5=1000; print $0}}' {output.broadpeak_bed}.sorted > {output.broadpeak_bed}.tmp
			bedToBigBed -as={Template_Path}/bigBroadPeak.as -type=bed6+3 {output.broadpeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigbed}

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
			printf "%s\\n" "macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg \\
			--o-prefix ${{sample_Name}}.macs2_broad --outdir $RESULT_PATH --method FE --pseudocount 0.00001" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.tmp > {output.broadpeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			macs2 bdgcmp -t $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $RESULT_PATH --method FE --pseudocount 0.00001
			slopBed -i $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.tmp > {output.broadpeak_bdg}

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
			printf "%s\\n" "LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadpeak_bdg} > {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2,3 {output.broadpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,2 {output.broadpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "cut -f1,3 {output.broadpeak_bdg}.uniq | sort | uniq -d | sed 's/\\\\t/\\\\t.*\\\\t/' >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -s $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then" | tee >(cat >&2)
			line="\$line"
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "\$line" {output.broadpeak_bdg}.uniq | \\
			cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done" | tee >(cat >&2)
			printf "\\t%s\\n" "cat $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "\${{line}}d" {output.broadpeak_bdg}.uniq > {output.broadpeak_bdg}.uniq.fix; done" | tee >(cat >&2)
			printf "\\t%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.uniq.fix > {output.broadpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.broadpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "else" | tee >(cat >&2)
			printf "\\t%s\\n" "bedGraphToBigWig {output.broadpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			LC_COLLATE=C sort -u -k1,1 -k2,2n -k3,3n -s {output.broadpeak_bdg} > {output.broadpeak_bdg}.uniq
			cut -f1,2,3 {output.broadpeak_bdg}.uniq | sort | uniq -d > $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadpeak_bdg}.uniq | sort | uniq -d >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadpeak_bdg}.uniq | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
			#

			if [ -s $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				cat $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt | while read line; do grep -nr -m 1 --perl-regex "$line" {output.broadpeak_bdg}.uniq | cut -d":" -f1 >> $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done
				cat $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt | while read line; do sed "${{line}}d" {output.broadpeak_bdg}.uniq > {output.broadpeak_bdg}.uniq.fix; done
				LC_COLLATE=C sort -k1,1 -k2,2n {output.broadpeak_bdg}.uniq.fix > {output.broadpeak_bdg}.uniq.fix.sorted
				bedGraphToBigWig {output.broadpeak_bdg}.uniq.fix.sorted {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}
			else
				bedGraphToBigWig {output.broadpeak_bdg}.uniq {config_reference_Dict[CHROM_SIZE]} {output.broadpeak_bigwig}
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
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "if [ -f {output.broadpeak_bigwig} ]; then" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bed}.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.tmp" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq.fix" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf {output.broadpeak_bdg}.uniq.fix.sorted" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt" | tee >(cat >&2)
			printf "\\t%s\\n" "rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt" | tee >(cat >&2)
			printf "%s\\n" "fi" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			if [ -f {output.broadpeak_bigwig} ]; then
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $RESULT_PATH/${{sample_Name}}.macs2_broad_FE.bdg
				rm -rf {output.broadpeak_bed}.tmp
				rm -rf {output.broadpeak_bed}.sorted
				rm -rf {output.broadpeak_bdg}.tmp
				rm -rf {output.broadpeak_bdg}.uniq
				rm -rf {output.broadpeak_bdg}.uniq.fix
				rm -rf {output.broadpeak_bdg}.uniq.fix.sorted
				rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $RESULT_PATH/${{sample_Name}}.broad.duplicate_line.txt
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
