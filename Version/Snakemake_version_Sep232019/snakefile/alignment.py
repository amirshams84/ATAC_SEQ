# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
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


def get_case_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{case}.bam".format(design=wildcards.design, case=sample))
	return bam_List
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
EXECUTION_MODE = config_general_Dict["EXECUTION_MODE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
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
config_picard_Dict = config["PICARD"][PLATFORM]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ################################### WILDCARDS ################################ #


pre_process_List = []
alignment_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq".format(design=sample_Dict["Design"], sample=sample))
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		alignment_List
# ################################### PIPELINE RULES ########################## #


if LAYOUT == "paired":
	rule Alignment_Paired:
		"""
		"""
		input:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq",
		output:
			bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
			bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bed.gz",
			bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bigwig",
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = 200000
		message: "Alignment_Paired: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fwd_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq", "", each_fastq_basename)
			shell("""
				#
				##
				DATA_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment
				mkdir -p $DATA_PATH
				
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				
				SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
				mkdir -p $SCRATCH_PATH

				AWK_ALIGNMENT_FILTER="awk \'BEGIN{{OFS=FS}}{{if ( \$3 != \\\"chrUn\\\" && \$3 !~ /chrUn/ && \$3 !~ /random/ ) print \$0}}\'"
				##
				#
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "Name: %s\\n" "ALIGNMENT" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{input.processed_fwd_fastq}" | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{input.processed_rev_fastq}" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.bam}" | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.bam_index}" | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$DATA_PATH/{each_fastq_begining}.R1.mapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT4: %s\\n" "$DATA_PATH/{each_fastq_begining}.R2.mapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT5: %s\\n" "$DATA_PATH/{each_fastq_begining}.R1.unmapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT6: %s\\n" "$DATA_PATH/{each_fastq_begining}.R2.unmapped.fastq.gz" | tee >(cat >&2)
				printf "OUTPUT7: %s\\n" "$QC_PATH/{each_fastq_begining}.alignment.txt" | tee >(cat >&2)
				printf "OUTPUT8: %s\\n" "$QC_PATH/{each_fastq_begining}.picard.txt" | tee >(cat >&2)
				printf "OUTPUT9: %s\\n" "$QC_PATH/{each_fastq_begining}.samtools.flagstat.txt" | tee >(cat >&2)
				printf "OUTPUT10: %s\\n" "$QC_PATH/{each_fastq_begining}.samtools.idxstats.txt" | tee >(cat >&2)
				printf "OUTPUT11: %s\\n" "$QC_PATH/{each_fastq_begining}.R1.mapped_fastqc.html" | tee >(cat >&2)
				printf "OUTPUT12: %s\\n" "$QC_PATH/{each_fastq_begining}.R2.mapped_fastqc.html" | tee >(cat >&2)
				printf "OUTPUT13: %s\\n" "$QC_PATH/{each_fastq_begining}.R1.unmapped_fastqc.html" | tee >(cat >&2)
				printf "OUTPUT14: %s\\n" "$QC_PATH/{each_fastq_begining}.R2.unmapped_fastqc.html" | tee >(cat >&2)
				printf "OUTPUT15: %s\\n" "$QC_PATH/{each_fastq_begining}.PBC.txt" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "module load bowtie/2-2.3.5 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load picard/2.18.27 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "module load fastqc/0.11.8 || exit 1" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "cd $SCRATCH_PATH" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} \\
				--un-conc-gz $DATA_PATH/{wildcards.sample}_R%.unmapped.fastq.gz --al-conc-gz $DATA_PATH/{wildcards.sample}_R%.mapped.fastq.gz 2> $QC_PATH/{each_fastq_begining}.alignment.txt | \\
				samtools view --threads {threads} -Sh -f 2 -F 256 /dev/stdin | $AWK_ALIGNMENT_FILTER | \\
				samtools sort --threads {threads} -O bam -T $DATA_PATH/{wildcards.sample} -o {output.bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar \$PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp \\
				OUTPUT={output.bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "samtools flagstat --threads {threads} {output.bam} > $QC_PATH/{each_fastq_begining}.samtools.flagstat.txt" | tee >(cat >&2)
				printf "%s\\n" "samtools idxstats --threads {threads} {output.bam} > $QC_PATH/{each_fastq_begining}.samtools.idxstats.txt" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $QC_PATH -f fastq --threads {threads} $DATA_PATH/{wildcards.sample}_R1.mapped.fastq.gz" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $QC_PATH -f fastq --threads {threads} $DATA_PATH/{wildcards.sample}_R2.mapped.fastq.gz" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $QC_PATH -f fastq --threads {threads} $DATA_PATH/{wildcards.sample}_R1.unmapped.fastq.gz" | tee >(cat >&2)
				printf "%s\\n" "fastqc -o $QC_PATH -f fastq --threads {threads} $DATA_PATH/{wildcards.sample}_R2.unmapped.fastq.gz" | tee >(cat >&2)
				printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.bam} --outfile $QC_PATH/{wildcards.sample}.PBC.txt --cores {threads}" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				#
				##
				start_time="$(date -u +%s)"
				module load bowtie/2-2.3.5 || exit 1
				module load samtools/1.9 || exit 1
				module load picard/2.18.27 || exit 1
				module load fastqc/0.11.8 || exit 1
				module load bedtools/2.27.1 || exit 1
				module load ucsc/373 || exit 1
				module load deeptools/3.1.3 || exit 1

				cd $SCRATCH_PATH

				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} \\
				--un-conc-gz $DATA_PATH/{wildcards.sample}_R%.unmapped.fastq.gz --al-conc-gz $DATA_PATH/{wildcards.sample}_R%.mapped.fastq.gz 2> $QC_PATH/{each_fastq_begining}.alignment.txt | \\
				samtools view --threads {threads} -Sh -f 2 -F 256 /dev/stdin | awk 'BEGIN{{OFS=FS}}{{if ( $3 != \"chrUn\" && $3 !~ /chrUn/ && $3 !~ /random/ ) print $0}}' | \\
				samtools sort --threads {threads} -O bam -T $DATA_PATH/{wildcards.sample} -o {output.bam}.tmp -

				java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp \\
				OUTPUT={output.bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_picard_Dict}

				rm -rf {output.bam}.tmp
				samtools index -@ {threads} -b {output.bam}

				bedtools bamtobed -i {output.bam} > {output.bed}.tmp
				LC_COLLATE=C sort -k1,1 -k2,2n {output.bed}.tmp > {output.bed}.sorted
				bgzip -c {output.bed}.sorted > {output.bed}
				tabix -f -p bed {output.bed}
				rm -rf {output.bed}.tmp {output.bed}.sorted

				bamCoverage --bam {output.bam} --outFileName {output.bigwig} --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
				--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]}


				samtools flagstat --threads {threads} {output.bam} > $QC_PATH/{each_fastq_begining}.samtools.flagstat.txt
				samtools idxstats --threads {threads} {output.bam} > $QC_PATH/{each_fastq_begining}.samtools.idxstats.txt

				fastqc -o $QC_PATH -f bam --threads {threads} {output.bam}
				
				printf "%s\\n" "#!/bin/bash" > $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "module load python/2.7" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.bam} --outfile $QC_PATH/{wildcards.sample}.PBC.txt --cores {threads}" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh
				printf "%s\\n" "echo 'Pipeline execution successfully finished at: '\$(date)" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC.sh

				cd $QC_PATH
				sbatch --mem=300G --cpus-per-task=32 --partition=largemem --time=1-00:00:00 ./{wildcards.design}_{wildcards.sample}_PBC.sh

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "%s\\n" "#" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
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


rule Pooling_Case_Replicates:
	input:
		bam_List = get_case_bam
	output:
		pooled_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam",
		pooled_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam.bai",
		pooled_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bed.gz",
		pooled_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bigwig",
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = 200000
	message: "Pooling_Case_Replicates: {wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			#
			##
			DATA_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment
			mkdir -p $DATA_PATH
			
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
			mkdir -p $QC_PATH

			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH
			##
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "Name: %s\\n" "POOLING_CASE_REPLICATE" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			declare -a bam_List=({input.bam_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.pooled_bam}" | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.pooled_bam_index}" | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "$QC_PATH/{wildcards.pooled_case}.picard.txt" | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "$QC_PATH/{wildcards.pooled_case}.samtools.flagstat.txt" | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "$QC_PATH/{wildcards.pooled_case}.samtools.idxstats.txt" | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "$QC_PATH/{wildcards.pooled_case}.PBC.txt" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load picard/2.18.27 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load fastqc/0.11.8 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "samtools merge --threads {threads} {output.pooled_bam}.unsorted {input.bam_List}" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted" | tee >(cat >&2)
			printf "%s\\n" "java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar \$PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted \\
			OUTPUT={output.pooled_bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$QC_PATH/{wildcards.pooled_case}.picard.txt {config_picard_Dict}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} {output.pooled_bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.pooled_bam}.unsorted {output.pooled_bam}.sorted" | tee >(cat >&2)
			printf "%s\\n" "samtools flagstat --threads {threads} {output.pooled_bam} > $QC_PATH/{wildcards.pooled_case}.samtools.flagstat.txt" | tee >(cat >&2)
			printf "%s\\n" "samtools idxstats --threads {threads} {output.pooled_bam} > $QC_PATH/{wildcards.pooled_case}.samtools.idxstats.txt" | tee >(cat >&2)
			printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.pooled_bam} --outfile $QC_PATH/{wildcards.pooled_case}.PBC.txt --cores {threads}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"
			module load bowtie/2-2.3.5 || exit 1
			module load samtools/1.9 || exit 1
			module load picard/2.18.27 || exit 1
			module load fastqc/0.11.8 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load ucsc/373 || exit 1
			module load deeptools/3.1.3 || exit 1

			samtools merge --threads {threads} {output.pooled_bam}.unsorted {input.bam_List}
			samtools sort --threads {threads} -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted
			
			java -Xms50000M -Xmx50000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted \\
			OUTPUT={output.pooled_bam} TMP_DIR=$SCRATCH_PATH METRICS_FILE=$QC_PATH/{wildcards.pooled_case}.picard.txt {config_picard_Dict}

			samtools index -@ {threads} {output.pooled_bam}
			rm -rf {output.pooled_bam}.unsorted {output.pooled_bam}.sorted
			
			bedtools bamtobed -i {output.pooled_bam} > {output.pooled_bed}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.pooled_bed}.tmp > {output.pooled_bed}.sorted
			bgzip -c {output.pooled_bed}.sorted > {output.pooled_bed}
			tabix -f -p bed {output.pooled_bed}
			rm -rf {output.pooled_bed}.tmp {output.pooled_bed}.sorted

			bamCoverage --bam {output.pooled_bam} --outFileName {output.pooled_bigwig} --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
			--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]}

			samtools flagstat --threads {threads} {output.pooled_bam} > $QC_PATH/{wildcards.pooled_case}.samtools.flagstat.txt
			samtools idxstats --threads {threads} {output.pooled_bam} > $QC_PATH/{wildcards.pooled_case}.samtools.idxstats.txt

			fastqc -o $QC_PATH -f bam --threads {threads} {output.pooled_bam}

			printf "%s\\n" "#!/bin/bash" > $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "module load python/2.7" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.pooled_bam} --outfile $QC_PATH/{wildcards.pooled_case}.PBC.txt --cores {threads}" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh
			printf "%s\\n" "echo 'Pipeline execution successfully finished at: '\$(date)" >> $QC_PATH/{wildcards.design}_{wildcards.pooled_case}_PBC.sh

			cd $QC_PATH
			sbatch --mem=300G --cpus-per-task=32 --partition=largemem --time=1-00:00:00 ./{wildcards.design}_{wildcards.pooled_case}_PBC.sh

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)

		""")
