# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Project: ENCODE ATAC_SEQ
# Aim: Snakemake workflow for alignment
# snakemake --snakefile alignment.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile alignment.py --configfile Encode.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# ################################### IMPORT ##################################### #

import os
import sys
import re
library_path = os.path.abspath(workflow.basedir + "/../library/")
sys.path.append(library_path)
import utility
# ################################### FUNCTIONS ################################## #


def get_case_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in sample_treatment_Dict.items():
		if sample_Dict[TREATMENT_COLUMN] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{case}.bam".format(design=wildcards.design, case=sample))
	return bam_List
# ################################### CONFIGURATION ############################## #


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
PLATFORM = config_data_Dict["PLATFORM"].lower()
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
config_picard_Dict = config["PICARD"][PLATFORM]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ################################### WILDCARDS ################################ #


pre_process_List = []
alignment_List = []
for sample, sample_Dict in sample_treatment_Dict.items():
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))

for design in metadata_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))

# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		alignment_List
# ################################### PIPELINE RULES ########################## #


if LAYOUT == "paired":
	#
	rule alignment_paired:
		"""
		"""
		input:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq",
		output:
			bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai"
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "alignment_paired: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fwd_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq", "", each_fastq_basename)
			shell("""
				#
				##
				total_start_time="$(date -u +%s)"

				RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment
				mkdir -p $RESULT_PATH
				
				REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $REPORT_PATH
				
				SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
				mkdir -p $SCRATCH_PATH

				AWK_ALIGNMENT_FILTER="awk \'BEGIN{{OFS=FS}}{{if ( \$3 != \\\"chrUn\\\" && \$3 !~ /chrUn/ && \$3 !~ /random/ ) print \$0}}\'"
				##
				#
				printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
				printf "%s\\n" "alignment_paired: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
				printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{input.processed_fwd_fastq}" | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{input.processed_rev_fastq}" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.bam}" | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.bam_index}" | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$RESULT_PATH/{each_fastq_begining}.R1.mapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT4: %s\\n" "$RESULT_PATH/{each_fastq_begining}.R2.mapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT5: %s\\n" "$RESULT_PATH/{each_fastq_begining}.R1.unmapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT6: %s\\n" "$RESULT_PATH/{each_fastq_begining}.R2.unmapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT7: %s\\n" "$REPORT_PATH/{each_fastq_begining}.alignment.txt" | tee >(cat >&2)
				printf "OUTPUT8: %s\\n" "$REPORT_PATH/{each_fastq_begining}.picard.txt" | tee >(cat >&2)
				printf "OUTPUT9: %s\\n" "$REPORT_PATH/{each_fastq_begining}.samtools.flagstat.txt" | tee >(cat >&2)
				printf "OUTPUT10: %s\\n" "$REPORT_PATH/{each_fastq_begining}.samtools.idxstats.txt" | tee >(cat >&2)
				printf "OUTPUT11: %s\\n" "$REPORT_PATH/{each_fastq_begining}_fastqc.html" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "###################################- EXECUTION -##############################" | tee >(cat >&2)
				printf "%s\\n" "module load bowtie/2-2.3.5 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load picard/2.18.27 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load fastqc/0.11.8 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				module load bowtie/2-2.3.5 || exit 1
				module load samtools/1.9 || exit 1
				module load picard/2.18.27 || exit 1
				module load bedtools/2.27.1 || exit 1
				module load deeptools/3.1.3 || exit 1
				module load ucsc/373 || exit 1
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
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} \\
				--un-conc-gz $RESULT_PATH/{wildcards.sample}_R%.unmapped.fastq.gz --al-conc-gz $RESULT_PATH/{wildcards.sample}_R%.mapped.fastq.gz 2> $REPORT_PATH/{each_fastq_begining}.alignment.txt | \\
				samtools view --threads {threads} -Sh -f 2 -F 256 /dev/stdin | $AWK_ALIGNMENT_FILTER | \\
				samtools sort --threads {threads} -O bam -T $RESULT_PATH/{wildcards.sample} -o {output.bam}.tmp -" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				cd $SCRATCH_PATH
				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} \\
				--un-conc-gz $RESULT_PATH/{wildcards.sample}_R%.unmapped.fastq.gz --al-conc-gz $RESULT_PATH/{wildcards.sample}_R%.mapped.fastq.gz 2> $REPORT_PATH/{each_fastq_begining}.alignment.txt | \\
				samtools view --threads {threads} -Sh -f 2 -F 256 /dev/stdin | awk 'BEGIN{{OFS=FS}}{{if ( $3 != \"chrUn\" && $3 !~ /chrUn/ && $3 !~ /random/ ) print $0}}' | \\
				samtools sort --threads {threads} -O bam -T $RESULT_PATH/{wildcards.sample} -o {output.bam}.tmp -

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
				printf "%s\\n" "java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar \$PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp \\
				OUTPUT={output.bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$REPORT_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp \\
				OUTPUT={output.bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$REPORT_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}
				rm -rf {output.bam}.tmp
				samtools index -@ {threads} -b {output.bam}

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
				printf "%s\\n" "samtools flagstat --threads {threads} {output.bam} > $REPORT_PATH/{each_fastq_begining}.samtools.flagstat.txt" | tee >(cat >&2)
				printf "%s\\n" "samtools idxstats --threads {threads} {output.bam} > $REPORT_PATH/{each_fastq_begining}.samtools.idxstats.txt" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $REPORT_PATH -f bam --threads {threads} {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "" | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"

				samtools flagstat --threads {threads} {output.bam} > $REPORT_PATH/{each_fastq_begining}.samtools.flagstat.txt
				samtools idxstats --threads {threads} {output.bam} > $REPORT_PATH/{each_fastq_begining}.samtools.idxstats.txt
				fastqc -o $REPORT_PATH -f bam --threads {threads} {output.bam}

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

elif LAYOUT == "single":
	rule Alignment_Single:
		"""
		"""
		input:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
		output:
			bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Single: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq", "", each_fastq_basename)
			shell("""
				#
				module load bowtie/2-2.3.5 || exit 1
				module load samtools/1.9 || exit 1
				module load picard/2.18.27 || exit 1
				module load fastqc/0.11.8 || exit 1
				module load qualimap/2.2.1 || exit 1
				#
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				#
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bowtie/2-2.3.5" | tee >(cat >&2)
				printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
				printf "%s\\n" "picard/2.18.27" | tee >(cat >&2)
				printf "%s\\n" "fastqc/0.11.8" | tee >(cat >&2)
				printf "%s\\n" "qualimap/2.2.1" | tee >(cat >&2)
				printf "%s\\n" "Align processed reads to specified genome and sort it, mark it with duplicate, index it"  | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{input.processed_fwd_fastq}"  | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{input.processed_rev_fastq}"  | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.bam}"  | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.bam_index}"  | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.alignment.txt"  | tee >(cat >&2)
				printf "OUTPUT4: %s\\n" "$QC_PATH/{each_fastq_begining}.picard.txt"  | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -U {input.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads}  -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads}  -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp -
				java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}
				samtools index -@ {threads} -b {output.bam}
				rm -rf {output.bam}.tmp
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			""")

rule pooling_case_replicates:
	"""
	"""
	input:
		bam_List = get_case_bam
	output:
		pooled_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam",
		pooled_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam.bai"
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "pooling_case_replicates: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"
			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment
			mkdir -p $RESULT_PATH
				
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
			mkdir -p $REPORT_PATH

			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "pooling_case_replicates: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.pooled_case}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			declare -a bam_List=({input.bam_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.pooled_bam}" | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.pooled_bam_index}" | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "$REPORT_PATH/{wildcards.pooled_case}.picard.txt" | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "$REPORT_PATH/{wildcards.pooled_case}.samtools.flagstat.txt" | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "$REPORT_PATH/{wildcards.pooled_case}.samtools.idxstats.txt" | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$REPORT_PATH/{wildcards.pooled_case}_fastqc.html" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load picard/2.18.27 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load fastqc/0.11.8 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load picard/2.18.27 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
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
			printf "%s\\n" "samtools merge --threads {threads} {output.pooled_bam}.unsorted {input.bam_List}" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			samtools merge --threads {threads} {output.pooled_bam}.unsorted {input.bam_List}
			samtools sort --threads {threads} -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted

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
			printf "%s\\n" "java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar \$PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted \\
			OUTPUT={output.pooled_bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$REPORT_PATH/{wildcards.pooled_case}.picard.txt {config_picard_Dict}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} {output.pooled_bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.pooled_bam}.unsorted {output.pooled_bam}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted \\
			OUTPUT={output.pooled_bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$REPORT_PATH/{wildcards.pooled_case}.picard.txt {config_picard_Dict}
			samtools index -@ {threads} {output.pooled_bam}
			rm -rf {output.pooled_bam}.unsorted {output.pooled_bam}.sorted

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
			printf "%s\\n" "samtools flagstat --threads {threads} {output.pooled_bam} > $REPORT_PATH/{wildcards.pooled_case}.samtools.flagstat.txt" | tee >(cat >&2)
			printf "%s\\n" "samtools idxstats --threads {threads} {output.pooled_bam} > $REPORT_PATH/{wildcards.pooled_case}.samtools.idxstats.txt" | tee >(cat >&2)
			printf "%s\\n" "fastqc -o $REPORT_PATH -f bam --threads {threads} {output.pooled_bam}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			samtools flagstat --threads {threads} {output.pooled_bam} > $REPORT_PATH/{wildcards.pooled_case}.samtools.flagstat.txt
			samtools idxstats --threads {threads} {output.pooled_bam} > $REPORT_PATH/{wildcards.pooled_case}.samtools.idxstats.txt
			fastqc -o $REPORT_PATH -f bam --threads {threads} {output.pooled_bam}

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
