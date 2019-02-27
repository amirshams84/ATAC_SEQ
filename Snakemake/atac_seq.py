shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Dec-28-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for NGS DATA Quality Control
# source /data/shamsaddinisha/conda/etc/profile.d/conda.sh
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --debug-dag --cores=50
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --cluster-config cluster.yml --debug-dag --cores=50 --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} --time={cluster.time} {cluster.extra}"
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --rulegraph | dot -Tsvg > atac_seq.svg
#https://hpc.nih.gov/~RTB/Amir/trial-project/trial-experiment/hg38/genome_browser/ucsc/genomes.txt
# ################################### TO DO LIST ##################################### #
"""
dynamically assign cpu and memory core
better conda environment
igv session
work with illumina and sra(R1,R2,r1,r2)
- controlled
- IDR
- bigNarrowPeak

Centipede
IGV track hub
ucsc track hub
"""
# ################################### IMPORT ##################################### #


import os
import re
from os.path import join
import sys
import glob
import logging as log
import itertools
import collections
import multiprocessing
import pandas
sys.path.append(os.path.abspath("./library/"))
import utility
import execution_report
# ################################### WILDCARDS FUNCTIONS ############################# #


def get_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + FASTQ_DELIMITER + "*.fastq*", recursive=True)


def get_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{case}.picard.filter_mapq.filter_chromosome.bam".format(design=wildcards.design, case=sample))
	return bam_List


def get_narrowpeak(wildcards):
	"""
	"""
	narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{case}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=wildcards.design, case=sample))
	return narrowPeak_List


def get_pooled_narrowpeak(wildcards):
	"""
	"""
	pooled_narrowPeak_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			pooled_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
		else:
			pass
	pooled_narrowPeak_List = list(set(pooled_narrowPeak_List))
	return pooled_narrowPeak_List

# ################################### CLUSTER CONFIGURATION ########################### #


PROCESSORS = 20
MEMORY = 50000
# ################################### PIPELINE CONFIGURATION ######################### #


configfile: "atac_seq.json"
localrules: End_Point
# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
INFOLINK = config_general_Dict["INFOLINK"]
TITLE = config_general_Dict["TITLE"]
GENOME = config_general_Dict["GENOME"]
GENOME = GENOME.lower()
FASTQ_DELIMITER = config_general_Dict["FASTQ_DELIMITER"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CONDA
config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CONDA_PY2 = config_conda_Dict["ATAC_Seq_py2"]
CONDA_PY3 = config_conda_Dict["ATAC_Seq_py3"]
ACTIVATE_CONDA_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_CONDA_PY3 = config_conda_Dict["ACTIVATE_PY3"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
metadata_Dict = utility.build_metadata_dict(METADATA_FILE)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PRE_PROCESS
config_pre_process_Dict = config["PRE_PROCESS"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#SUBSAMPLE
config_subsample_Dict = config["SUBSAMPLE"]
SUBSAMPLE_PCT = config_subsample_Dict["SUBSAMPLE_PCT"]
USE_SUBSAMPLE = config_subsample_Dict["USE_SUBSAMPLE"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_post_alignment_Dict = config["POST_ALIGNMENT"]
CHROMOSOME_FILTER_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["FILTER_CHROMOSOME"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_peak_analysis_Dict = config["PEAK_ANALYSIS"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_track_hub_Dict = config["TRACK_HUB"]
IGV_server_List = [WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/igv_trackhub.xml"]
# ################################### WILDCARDS ###################################### #

design_Dict = {}
#
pre_process_processed_fastq_List = []
#
alignment_bowtie2_List = []
alignment_picard_List = []
alignment_bamqc_List = []
#
post_alignment_filter_mapq_List = []
post_alignment_filter_chromosome_List = []
#
peak_calling_macs2_narrowPeak_List = []
#
peak_analysis_homer_peak_annotate_List = []

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
	#
	#PRE_PROCESS
	##PROCESSED_FASTQ
	pre_process_processed_fastq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample}.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	#ALIGNMENT
	##BOWTIE2
	alignment_bowtie2_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	alignment_bowtie2_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample}.bam.bai".format(design=sample_Dict["Design"], sample=sample))
	##PICARD
	alignment_picard_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample}.picard.bam".format(design=sample_Dict["Design"], sample=sample))
	alignment_picard_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample}.picard.bam.bai".format(design=sample_Dict["Design"], sample=sample))
	##QUALIMAP
	alignment_bamqc_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample}_qualimap/qualimapReport.html".format(design=sample_Dict["Design"], sample=sample))
	#POST-ALIGNMENT
	##FILTER_MAPQ
	post_alignment_filter_mapq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample}.picard.filter_mapq.bam".format(design=sample_Dict["Design"], sample=sample))
	##FILTER_CHROMOSOME
	post_alignment_filter_chromosome_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample}.picard.filter_mapq.filter_chromosome.bam".format(design=sample_Dict["Design"], sample=sample))
	
	#PEAK-CALLING
	##MACS2_NARROWPEAK
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.bigwig".format(design=sample_Dict["Design"], sample=sample))
	#PEAK_ANALYSIS
	##HOMER_PEAK_ANNOTATE
	peak_analysis_homer_peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{sample}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate.txt".format(design=sample_Dict["Design"], sample=sample))
else:
	pass


post_alignment_pooling_List = []
peak_analysis_naive_overlap_List = []
peak_analysis_idr_overlap_List = []
for design in design_Dict:
	#
	##POOLING
	post_alignment_pooling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case}.picard.filter_mapq.filter_chromosome.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#PEAK-CALLING
	##MACS2_NARROWPEAK
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.bigwig".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#POST-ANALYSIS
	##NAIVE_OVERLAP
	peak_analysis_naive_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{naive_overlapped}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=design, naive_overlapped="_NAIVE_OVERLAPPED_".join(design_Dict[design]["Case"])))
	##IDR
	peak_analysis_idr_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz".format(design=design, idr_overlapped="_IDR_OVERLAPPED_".join(design_Dict[design]["Case"])))
	##HOMER_PEAK_ANNOATE
	peak_analysis_homer_peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{pooled_case}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_analysis_homer_peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{naive_overlapped}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate.txt".format(design=design, naive_overlapped="_NAIVE_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_analysis_homer_peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{idr_overlapped}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate.txt".format(design=design, idr_overlapped="_IDR_OVERLAPPED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################## #


pre_process_List = []
pre_process_List.extend(pre_process_processed_fastq_List)
#
alignment_List = []
alignment_List.extend(alignment_bowtie2_List)
alignment_List.extend(alignment_picard_List)
alignment_List.extend(alignment_bamqc_List)
#
post_alignment_List = []
post_alignment_List.extend(post_alignment_filter_mapq_List)
post_alignment_List.extend(post_alignment_filter_chromosome_List)
post_alignment_List.extend(post_alignment_pooling_List)
#
peak_calling_List = []
peak_calling_List.extend(peak_calling_macs2_narrowPeak_List)
#
peak_analysis_List = []
peak_analysis_List.extend(peak_analysis_naive_overlap_List)
peak_analysis_List.extend(peak_analysis_idr_overlap_List)
peak_analysis_List.extend(peak_analysis_homer_peak_annotate_List)
#
genome_browser_List = []
genome_browser_List.extend(IGV_server_List)
rule End_Point:
	input:
		pre_process_List
		+ alignment_List
		+ post_alignment_List
		+ peak_calling_List
		+ peak_analysis_List
		+ genome_browser_List

# ################################### PIPELINE RULES ########################### #

ruleorder:
	Cutadapt > Bowtie2 > Picard > BamQC > Filter_mapQ > Filter_Chromosome > Pooling > Macs2_NarrowPeak > Macs2_NarrowPeak_Naive_Overlap > Macs2_NarrowPeak_IDR_Overlap > Homer_Peak_Annotate > IGV_Trackhub

#+++++++++++++++++++++++++++++
##REALM1: PRE_PROCESS
#+++++++++++++++++++++++++++++


rule Cutadapt:
	input:
		fastq_List = get_fastq
	output:
		processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.fastq.gz",
		discarded_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.discarded.fastq.gz",
	log:
		cutadapt_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.cutadapt.stdout",
		cutadapt_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.cutadapt.stderr",
		usearch_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/usearch/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.usearch.stdout",
		usearch_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/usearch/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.usearch.stderr",
	priority: 999
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Cutadapt: {wildcards.design}|{wildcards.sample}"
	run:
		for each_fastq in input.fastq_List:
			#
			each_fastq_basename = os.path.basename(each_fastq)
			each_fastq_begining = re.sub(".fastq.*", "", each_fastq_basename)
			#
			if each_fastq.split(".")[-1] == "gz":
				#
				each_processed_fastq = each_fastq_begining + ".processed.fastq"
				each_discarded_fastq = each_fastq_begining + ".discarded.fastq"
				shell("""
					{ACTIVATE_CONDA_PY2}
					##
					##EXECUTION
					start_time="$(date -u +%s)"
					mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt
					cutadapt_Path={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt
					cutadapt {config_pre_process_Dict[CUTADAPT]} --too-short-output=$cutadapt_Path/{each_discarded_fastq}.gz --output=$cutadapt_Path/{each_processed_fastq}.gz {each_fastq} 1> {log.cutadapt_stdout} 2> {log.cutadapt_stderr}

					gunzip < $cutadapt_Path/{each_processed_fastq}.gz > $cutadapt_Path/{each_processed_fastq}.tmp
					mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/usearch
					usearch_Path={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/usearch
					module load usearch 1> {log.usearch_stdout} 2> {log.usearch_stderr}

					if [ {USE_SUBSAMPLE} = "TRUE" ]; then
						usearch -fastx_subsample $cutadapt_Path/{each_processed_fastq}.tmp -fastqout $usearch_Path/{each_processed_fastq}.subsample -sample_pct {SUBSAMPLE_PCT} >> {log.usearch_stdout} 2>> {log.usearch_stderr}
						usearch -fastq_filter $usearch_Path/{each_processed_fastq}.subsample {config_pre_process_Dict[USEARCH]} -relabel "{each_fastq_begining};Read_" -fastqout $usearch_Path/{each_processed_fastq} 1>> {log.usearch_stdout} 2>> {log.usearch_stderr}
						rm -rf $usearch_Path/{each_processed_fastq}.subsample
					else
						usearch -fastq_filter $cutadapt_Path/{each_processed_fastq}.tmp {config_pre_process_Dict[USEARCH]} -relabel "{each_fastq_begining};Read_" -fastqout $usearch_Path/{each_processed_fastq} 1>> {log.usearch_stdout} 2>> {log.usearch_stderr}
					fi

					cat $usearch_Path/{each_processed_fastq} >> {output.processed_fastq}.tmp
					cat $cutadapt_Path/{each_discarded_fastq}.gz >> {output.discarded_fastq}
					rm -rf $cutadapt_Path/{each_processed_fastq}.tmp $cutadapt_Path/{each_processed_fastq}.gz $usearch_Path/{each_processed_fastq}
					rm -rf $cutadapt_Path/{each_discarded_fastq}.gz
					end_time="$(date -u +%s)"
					elapsed_time=$(($end_time-$start_time))
					##
					##
				""")
			else:
				pass
		else:
			pass
		shell("""
			gzip < {output.processed_fastq}.tmp > {output.processed_fastq}
			rm -rf {output.processed_fastq}.tmp
		""")


#+++++++++++++++++++++++++++++
##REALM2: ALIGNMENT
#+++++++++++++++++++++++++++++


rule Bowtie2:
	input:
		processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.fastq.gz",
	output:
		bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
	log:
		bowtie2_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.report",
		bowtie2_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stdout",
		bowtie2_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stderr",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Bowtie2: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""

			{ACTIVATE_CONDA_PY2}
			##
			##EXECUTION
			start_time="$(date -u +%s)"
			bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} {input.processed_fastq} 2> {log.bowtie2_stdout} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample} -o {output.bowtie2_Bam} - 1> {log.bowtie2_stdout} 2> {log.bowtie2_stderr}
			samtools index -@ {threads} -b {output.bowtie2_Bam} >> {log.bowtie2_stdout} 2>> {log.bowtie2_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

		""")


rule Picard:
	input:
		bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
	output:
		picard_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam",
		picard_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam.bai",
	log:
		picard_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.picard.report",
		picard_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.picard.stdout",
		picard_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.picard.stderr",

	priority: 995
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Picard: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			export _JAVA_OPTIONS="-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}" 1> {log.picard_stderr} 2> {log.picard_stderr}
			picard MarkDuplicates INPUT={input.bowtie2_Bam} OUTPUT={output.picard_Bam} METRICS_FILE={log.picard_report} {config_alignment_Dict[PICARD]} 1>> {log.picard_stderr} 2>> {log.picard_stderr}
			samtools index -@ {threads} -b {output.picard_Bam} 2>> {log.picard_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

		""")


rule BamQC:
	input:
		picard_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam",
	output:
		qualimap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}_qualimap/qualimapReport.html",
		
	log:
		qualimap_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.qualimap.report",
		qualimap_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.qualimap.stdout",
		qualimap_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.qualimap.stderr",
		#
		fastqc_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.fastqc.report",
		fastqc_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.fastqc.stdout",
		fastqc_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bamqc/{sample, ((?!.*_VS_.*|\\.).)*}.alignment.bamqc.fastqc.stderr",
	priority: 994
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "BamQC: {wildcards.design}|{wildcards.sample}"
	run:

		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##EXECUTION
			start_time="$(date -u +%s)"

			unset DISPLAY
			qualimap bamqc -bam {input.picard_Bam} -nt {threads} -outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bamqc/{wildcards.sample}_qualimap {config_alignment_Dict[QUALIMAP]} 1> {log.qualimap_stdout} 2> {log.qualimap_stderr}
			mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bamqc/{wildcards.sample}_fastqc
			fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bamqc/{wildcards.sample}_fastqc -f bam_mapped --threads {threads} {input.picard_Bam} 1> {log.fastqc_stdout} 2> {log.fastqc_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

		""")

#+++++++++++++++++++++++++++++
##REALM3: POST-ALIGNMENT
#+++++++++++++++++++++++++++++


rule Filter_mapQ:
	input:
		picard_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam",
		picard_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam.bai",
	output:
		filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.bam",
		filter_mapq_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.bam.report",
	log:
		filter_mapq_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_mapq.report",
		filter_mapq_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_mapq.stdout",
		filter_mapq_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_mapq.stderr",
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Filter_mapQ: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_MAPQ]} {input.picard_Bam} -b -O BAM -o {output.filter_mapq_Bam} 1> {log.filter_mapq_stdout} 2> {log.filter_mapq_stderr}
			samtools index -@ {threads} -b {output.filter_mapq_Bam} 2>> {log.filter_mapq_stderr}
			samtools flagstat --threads {threads} {output.filter_mapq_Bam} > {output.filter_mapq_Bam_Report} 2>> {log.filter_mapq_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")


rule Filter_Chromosome:
	input:
		filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.bam",
	output:
		filter_chromosome_filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam",
		filter_chromosome_filter_mapq_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam.report",
	log:
		filter_chromosome_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_chromosome.report",
		filter_chromosome_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_chromosome.stdout",
		filter_chromosome_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_chromosome.stderr",
	priority: 991
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Filter_Chromosome: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			#EXECUTION
			start_time="$(date -u +%s)"
			samtools view --threads {threads} -h {input.filter_mapq_Bam} | awk -F'\\t' '{CHROMOSOME_FILTER_PROCESS}' | samtools view --threads {threads} -Shb - > {output.filter_chromosome_filter_mapq_Bam} 2> {log.filter_chromosome_stderr}
			samtools index -@ {threads} -b {output.filter_chromosome_filter_mapq_Bam} 2>> {log.filter_chromosome_stderr}
			samtools flagstat --threads {threads} {output.filter_chromosome_filter_mapq_Bam} > {output.filter_chromosome_filter_mapq_Bam_Report} 2>> {log.filter_chromosome_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")

rule Pooling:
	input:
		bam_List = get_bam
	output:
		pooled_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case, .*_POOLED_.*}.picard.filter_mapq.filter_chromosome.bam",
		pooled_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case, .*_POOLED_.*}.picard.filter_mapq.filter_chromosome.bam.report",
	log:
		pooling_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.report",
		pooling_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stdout",
		pooling_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stderr",
	priority: 990
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Pooling: {wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			#EXECUTION
			start_time="$(date -u +%s)"
			samtools merge --threads {threads} {output.pooled_Bam} {input.bam_List} 1> {log.pooling_stdout} 2> {log.pooling_stderr}
			samtools index -@ {threads} {output.pooled_Bam} 2>> {log.pooling_stderr}
			samtools flagstat --threads {threads} {output.pooled_Bam} > {output.pooled_Bam_Report} 2>> {log.pooling_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")


#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++


rule Macs2_NarrowPeak:
	input:
		processed_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam",
	output:
		macs2_narrowPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz",
		macs2_narrowPeak_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.bigwig",
		
	log:
		macs2_narrowPeak_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_narrowPeak.report",
		macs2_narrowPeak_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_narrowPeak.stdout",
		macs2_narrowPeak_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_narrowPeak.stderr",
	priority: 989
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Macs2_NarrowPeak: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			##
			#EXECUTION
			start_time="$(date -u +%s)"
			sample_Name=$(basename {input.processed_Bam})
			sample_Name=${{sample_Name%.picard.filter_mapq.filter_chromosome.bam}}
			#NARROWPEAK
			macs2 callpeak --treatment {input.processed_Bam} --name ${{sample_Name}}.macs2_narrow {config_peak_calling_Dict[MACS2_NARROWPEAK]} --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ 1> {log.macs2_narrowPeak_stdout} 2> {log.macs2_narrowPeak_stderr}
			
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.macs2_narrowPeak_Bed}.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bedtools intersect -v -a {output.macs2_narrowPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.macs2_narrowPeak_Bed}.sorted.filtered.tmp 2>> {log.macs2_narrowPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.macs2_narrowPeak_Bed}.sorted.filtered.tmp > {output.macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bgzip -c {output.macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp > {output.macs2_narrowPeak_Bed} 2>> {log.macs2_narrowPeak_stderr}
			tabix -f -p bed {output.macs2_narrowPeak_Bed} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			rm -f {output.macs2_narrowPeak_Bed}.sorted.tmp {output.macs2_narrowPeak_Bed}.sorted.filtered.tmp {output.macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}

			#SIGNAL COVERAGE:FE
			macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_treat_pileup.bdg \
			-c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ {config_peak_calling_Dict[MACS2_FE_SIGNAL]} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}


			slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{output.macs2_narrowPeak_Bigwig}.tmp.bdg 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}

			LC_COLLATE=C sort -k1,1 -k2,2n {output.macs2_narrowPeak_Bigwig}.tmp.bdg   > {output.macs2_narrowPeak_Bigwig}.sorted.tmp.bdg  2>> {log.macs2_narrowPeak_stderr}
			#cat {output.macs2_narrowPeak_Bigwig}.sorted.tmp | sort -u -k1,1 -k2,2 -k3,3 -s > {output.macs2_narrowPeak_Bigwig}.uniq.bdg
			#LC_COLLATE=C sort -k1,1 -k2,2n {output.macs2_narrowPeak_Bigwig}.uniq.bdg > {output.macs2_narrowPeak_Bigwig}.sorted.uniq.bdg
			
			bedGraphToBigWig {output.macs2_narrowPeak_Bigwig}.sorted.tmp.bdg  {config_reference_Dict[CHROM_SIZE]} {output.macs2_narrowPeak_Bigwig} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}


			rm -rf {output.macs2_narrowPeak_Bed}*.tmp
			rm -rf {output.macs2_narrowPeak_Bigwig}*.bdg
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")


rule Macs2_NarrowPeak_Naive_Overlap:
	input:
		narrowPeak_List = get_narrowpeak,
		pooled_narrowPeak_List = get_pooled_narrowpeak
	output:
		naive_overlapped_narrowPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{naive_overlapped, .*_NAIVE_OVERLAPPED_.*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz",
		naive_overlapped_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{naive_overlapped, .*_NAIVE_OVERLAPPED_.*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.bigwig",
	log:
		Naive_Overlap_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{naive_overlapped, .*_NAIVE_OVERLAPPED_.*}.peak_calling.naive_overlap.stdout",
		Naive_Overlap_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{naive_overlapped, .*_NAIVE_OVERLAPPED_.*}.peak_calling.naive_overlap.stderr",
	priority: 988
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Naive_Overlap: {wildcards.design}|{wildcards.naive_overlapped}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			##
			#EXECUTION
			start_time="$(date -u +%s)"
			intersectBed -wo {config_peak_analysis_Dict[NAIVE_OVERLAP]} -a <(zcat -f {input.pooled_narrowPeak_List}) -b <(zcat -f {input.narrowPeak_List}) 2> {log.Naive_Overlap_stderr} \
			| awk 'BEGIN{{FS="\\t";OFS="\\t"}}  {{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.naive_overlapped_narrowPeak}.tmp 2>> {log.Naive_Overlap_stderr}
			
			
			cat {output.naive_overlapped_narrowPeak}.tmp | awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.naive_overlapped}.macs2_narrow_peak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.naive_overlapped_narrowPeak}.edited.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.naive_overlapped_narrowPeak}.edited.tmp > {output.naive_overlapped_narrowPeak}.sorted.tmp 2>> {log.Naive_Overlap_stderr}
			bgzip -c {output.naive_overlapped_narrowPeak}.sorted.tmp > {output.naive_overlapped_narrowPeak} 2>> {log.Naive_Overlap_stderr}
			tabix -f -p bed {output.naive_overlapped_narrowPeak} 1>> {log.Naive_Overlap_stdout} 2>> {log.Naive_Overlap_stderr}

			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n" , $1,$2,$3,$5}}' {output.naive_overlapped_narrowPeak}.sorted.tmp > {output.naive_overlapped_Bigwig}.bdg
			cat {output.naive_overlapped_Bigwig}.bdg | sort -u -k1,1 -k2,2 -k3,3 -s > {output.naive_overlapped_Bigwig}.uniq.bdg
			
			slopBed -i {output.naive_overlapped_Bigwig}.uniq.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.naive_overlapped_Bigwig}.uniq.edited.bdg 1>> {log.Naive_Overlap_stdout} 2>> {log.Naive_Overlap_stderr}

			LC_COLLATE=C sort -k1,1 -k2,2n {output.naive_overlapped_Bigwig}.uniq.edited.bdg > {output.naive_overlapped_Bigwig}.sorted.uniq.edited.bdg 2>> {log.Naive_Overlap_stderr}
			bedGraphToBigWig {output.naive_overlapped_Bigwig}.sorted.uniq.edited.bdg {config_reference_Dict[CHROM_SIZE]} {output.naive_overlapped_Bigwig} 1>> {log.Naive_Overlap_stdout} 2>> {log.Naive_Overlap_stderr}

			
			rm -rf {output.naive_overlapped_narrowPeak}.tmp {output.naive_overlapped_narrowPeak}.sorted.tmp {output.naive_overlapped_narrowPeak}.edited.tmp
			rm -rf {output.naive_overlapped_Bigwig}.bdg {output.naive_overlapped_Bigwig}.uniq.bdg {output.naive_overlapped_Bigwig}.uniq.edited.bdg {output.naive_overlapped_Bigwig}.sorted.uniq.edited.bdg

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

		""")


rule Macs2_NarrowPeak_IDR_Overlap:
	input:
		narrowPeak_List = get_narrowpeak,
		pooled_narrowPeak_List = get_pooled_narrowpeak
	output:
		idr_overlapped_narrowPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped, .*_IDR_OVERLAPPED_.*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz",
		idr_overlapped_narrowPeak_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped, .*_IDR_OVERLAPPED_.*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz.report",
		idr_overlapped_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped, .*_IDR_OVERLAPPED_.*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.bigwig",
	log:
		IDR_Overlap_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped, .*_IDR_OVERLAPPED_.*}.peak_calling.IDR_overlap.stdout",
		IDR_Overlap_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{idr_overlapped, .*_IDR_OVERLAPPED_.*}.peak_calling.IDR_overlap.stderr",
	priority: 987
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "IDR_Overlap: {wildcards.design}|{wildcards.idr_overlapped}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY3}
			##
			#EXECUTION
			start_time="$(date -u +%s)"
			idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak_List} --output-file {output.idr_overlapped_narrowPeak}.tmp {config_peak_analysis_Dict[IDR]} --log-output-file {output.idr_overlapped_narrowPeak_Report} > {log.IDR_Overlap_stdout} 2> {log.IDR_Overlap_stderr}
			idr_thresh_transformed=$(awk 'BEGIN{{print -log({config_peak_analysis_Dict[PVALUE]})/log(10)}}')
			awk 'BEGIN{{OFS="\\t"}} $12>='"${{idr_thresh_transformed}}"' {{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0"}}' {output.idr_overlapped_narrowPeak}.tmp | sort | uniq | LC_COLLATE=C sort -k1,1 -k2,2n > {output.idr_overlapped_narrowPeak}.filtered.tmp
			cat {output.idr_overlapped_narrowPeak}.filtered.tmp | awk 'BEGIN{{OFS="\\t"}} {{$4="{wildcards.idr_overlapped}.macs2_narrow_peak_"NR}} {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.idr_overlapped_narrowPeak}.edited.filtered.tmp
			
			{ACTIVATE_CONDA_PY2}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.idr_overlapped_narrowPeak}.edited.filtered.tmp > {output.idr_overlapped_narrowPeak}.sorted.tmp 2>> {log.IDR_Overlap_stderr}
			bgzip -c {output.idr_overlapped_narrowPeak}.sorted.tmp > {output.idr_overlapped_narrowPeak} 2>> {log.IDR_Overlap_stderr}
			tabix -f -p bed {output.idr_overlapped_narrowPeak} 1>> {log.IDR_Overlap_stdout} 2>> {log.IDR_Overlap_stderr}
			
			awk 'BEGIN{{OFS="\\t"}} {{printf "%s\\t%d\\t%d\\t%2.3f\\n" , $1,$2,$3,$5}}' {output.idr_overlapped_narrowPeak}.sorted.tmp > {output.idr_overlapped_Bigwig}.bdg
			cat {output.idr_overlapped_Bigwig}.bdg | sort -u -k1,1 -k2,2 -k3,3 -s > {output.idr_overlapped_Bigwig}.uniq.bdg
			
			slopBed -i {output.idr_overlapped_Bigwig}.uniq.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.idr_overlapped_Bigwig}.uniq.edited.bdg 1>> {log.IDR_Overlap_stdout} 2>> {log.IDR_Overlap_stderr}

			LC_COLLATE=C sort -k1,1 -k2,2n {output.idr_overlapped_Bigwig}.uniq.edited.bdg > {output.idr_overlapped_Bigwig}.sorted.uniq.edited.bdg 2>> {log.IDR_Overlap_stderr}
			bedGraphToBigWig {output.idr_overlapped_Bigwig}.sorted.uniq.edited.bdg {config_reference_Dict[CHROM_SIZE]} {output.idr_overlapped_Bigwig} 1>> {log.IDR_Overlap_stdout} 2>> {log.IDR_Overlap_stderr}

			rm -rf {output.idr_overlapped_narrowPeak}.tmp {output.idr_overlapped_narrowPeak}.filtered.tmp {output.idr_overlapped_narrowPeak}.edited.filtered.tmp {output.idr_overlapped_narrowPeak}.sorted.tmp
			rm -rf {output.idr_overlapped_Bigwig}.bdg  {output.idr_overlapped_Bigwig}.uniq.bdg {output.idr_overlapped_Bigwig}.uniq.edited.bdg {output.idr_overlapped_Bigwig}.sorted.uniq.edited.bdg

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")

#+++++++++++++++++++++++++++++
##REALM5: PEAK-ANALYSIS
#+++++++++++++++++++++++++++++


rule Homer_Peak_Annotate:
	input:
		narrowPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.gz",
	output:
		narrowPeak_homer_annotate = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate.txt",
		narrowPeak_homer_annotate_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_analysis/homer_peak_annotate/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.macs2_narrowPeak.narrowPeak.homer_annotate_report.txt",
	priority: 986
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Homer_Peak_Annotate: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY3}
			##
			#EXECUTION
			start_time="$(date -u +%s)"

			gunzip < {input.narrowPeak} > {output.narrowPeak_homer_annotate}.tmp
			annotatePeaks.pl {output.narrowPeak_homer_annotate}.tmp hg38 > {output.narrowPeak_homer_annotate} 2> {output.narrowPeak_homer_annotate_report}
			rm -f {output.narrowPeak_homer_annotate}.tmp

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
		""")




#+++++++++++++++++++++++++++++
##REALM6: GENOME-BROWSER
#+++++++++++++++++++++++++++++

rule IGV_Trackhub:
	input:
		alignment_List + post_alignment_List + peak_calling_List
	output:
		igv_trackhub = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/igv_trackhub.xml",
	priority: 900
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		IGV_String = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
		<Global name="''' + EXPERIMENT + '''"  infolink="''' + INFOLINK + '''" version="1">
		'''
		for design in design_Dict:
			#
			##
			IGV_String += '''
				<Category name="''' + design + '''">
			'''
			##BAM
			IGV_String += '''
					<Category name="Bam">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bam", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			##PEAK
			IGV_String += '''
					<Category name="Peak">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.narrowPeak.gz", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.broadPeak.gz", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			##BIGWIG
			IGV_String += '''
					<Category name="Bigwig">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bigwig", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			IGV_String += '''
				</Category>
			'''
		IGV_String += '''
		</Global>
		'''
		f = open(output.igv_trackhub, "w")
		f.write(IGV_String)
		f.close()


