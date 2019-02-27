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


def get_forward_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + "*R1*" + ".fastq.gz", recursive=True)


def get_reverse_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + "*R2*" + ".fastq.gz", recursive=True)


def get_processed_forward_fastq(wildcards):
	"""
	"""
	if USE_SUBSAMPLE is False:
		#
		return WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample}.R1.processed.fastq.gz".format(design=wildcards.design, sample=wildcards.sample)
	else:
		return WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample}.R1.processed.subsample.fastq.gz".format(design=wildcards.design, sample=wildcards.sample)


def get_processed_reverse_fastq(wildcards):
	"""
	"""
	if USE_SUBSAMPLE is False:
		#
		return WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample}.R2.processed.fastq.gz".format(design=wildcards.design, sample=wildcards.sample)
	else:
		return WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample}.R2.processed.subsample.fastq.gz".format(design=wildcards.design, sample=wildcards.sample)


def get_noMito_bam_list(wildcards):
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/samtools/{case}.noMito.bam".format(design=wildcards.design, case=sample))
	return bam_List


def get_filter_tn5_filter_chromosome_Bed(wildcards):
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{case}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz".format(design=wildcards.design, case=sample))
	return bam_List


def get_filter_tn5_filter_chromosome_noduplicate_Bed(wildcards):
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{case}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz".format(design=wildcards.design, case=sample))
	return bam_List
# ################################### CLUSTER CONFIGURATION ###################### #


PROCESSORS = 20
MEMORY = 50000
# ################################### PIPELINE CONFIGURATION ######################## #


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
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#SUBSAMPLE
config_subsample_Dict = config["SUBSAMPLE"]
SUBSAMPLE_PCT = config_subsample_Dict["SUBSAMPLE_PCT"]
USE_SUBSAMPLE = config_subsample_Dict["USE_SUBSAMPLE"]
if USE_SUBSAMPLE == "TRUE":
	USE_SUBSAMPLE = True
else:
	USE_SUBSAMPLE = False
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
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_post_alignment_Dict = config["POST_ALIGNMENT"]
CHROMOSOME_FILTER_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["FILTER_CHROMOSOME"])
AWK_TN5_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["FILTER_TN5"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_track_hub_Dict = config["TRACK_HUB"]
IGV_server_List = [WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/igv_trackhub.xml"]
UCSC_server_List = [WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/ucsc_trackhub.txt"]
#data_share_Link_prefix = config_track_hub_Dict["DATA_SHARE_LINK"]
# ------------------------------------


# ################################### WILDCARDS ###################################### #

design_Dict = {}
#
pre_process_processed_fastq_List = []
pre_process_subsample_fastq_List = []
pre_process_fastqc_List = []
#
alignment_bowtie2_List = []
alignment_picard_List = []
alignment_qualimap_List = []
#
post_alignment_filter_mapq_List = []
post_alignment_filter_duplicate_List = []
post_alignment_filter_chromosome_List = []
post_alignment_filter_tn5_List = []
#
peak_calling_macs2_narrowPeak_List = []
peak_calling_macs2_broadPeak_List = []
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
	else:
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
	pre_process_processed_fastq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	pre_process_processed_fastq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample}.R2.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	##FASTQC
	pre_process_fastqc_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample}.R1.processed_fastqc.zip".format(design=sample_Dict["Design"], sample=sample))
	pre_process_fastqc_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample}.R2.processed_fastqc.zip".format(design=sample_Dict["Design"], sample=sample))
	##SUBSAMPLE_FASTQ
	pre_process_subsample_fastq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample}.R1.processed.subsample.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	pre_process_subsample_fastq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample}.R2.processed.subsample.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	#ALIGNMENT
	##BOWTIE2
	alignment_bowtie2_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	alignment_bowtie2_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/bowtie2/{sample}.bam.bai".format(design=sample_Dict["Design"], sample=sample))
	##PICARD
	alignment_picard_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample}.picard.bam".format(design=sample_Dict["Design"], sample=sample))
	alignment_picard_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample}.picard.bam.bai".format(design=sample_Dict["Design"], sample=sample))
	##QUALIMAP
	alignment_qualimap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/qualimap/{sample}/qualimapReport.html".format(design=sample_Dict["Design"], sample=sample))
	#POST-SLIGNMENT
	##FILTER_MAPQ
	post_alignment_filter_mapq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample}.picard.filter_mapq.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_filter_mapq_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample}.picard.filter_mapq_discarded.bam".format(design=sample_Dict["Design"], sample=sample))
	##FILTER_DUPLICATE
	post_alignment_filter_duplicate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample}.picard.filter_mapq.noduplicate.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_filter_duplicate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample}.picard.filter_mapq.duplicate.bam".format(design=sample_Dict["Design"], sample=sample))
	##FILTER_CHROMOSOME
	post_alignment_filter_chromosome_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample}.picard.filter_mapq.filter_chromosome.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_filter_chromosome_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample}.picard.filter_mapq.noduplicate.filter_chromosome.bam".format(design=sample_Dict["Design"], sample=sample))
	##FILTER_TN5
	post_alignment_filter_tn5_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_filter_tn5_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz".format(design=sample_Dict["Design"], sample=sample))
	#PEAK-CALLING
	##MACS2_NARROWPEAK
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bedgraph.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bigwig".format(design=sample_Dict["Design"], sample=sample))
	##MACS2_BROADPEAK
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bedgraph.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bigwig".format(design=sample_Dict["Design"], sample=sample))
else:
	pass

# #######
post_alignment_pooling_List = []
for design in design_Dict:
	#
	##POOLING
	post_alignment_pooling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_pooling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#PEAK-CALLING
	##MACS2_NARROWPEAK
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bedgraph.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_macs2_narrowPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bigwig".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##MACS2_BROADPEAK
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bedgraph.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_macs2_broadPeak_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{pooled_case}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bigwig".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ###################################### #


pre_process_List = []
pre_process_List.extend(pre_process_processed_fastq_List)
pre_process_List.extend(pre_process_fastqc_List)
pre_process_List.extend(pre_process_subsample_fastq_List)
#
alignment_List = []
alignment_List.extend(alignment_bowtie2_List)
alignment_List.extend(alignment_picard_List)
alignment_List.extend(alignment_qualimap_List)
#
post_alignment_List = []
post_alignment_List.extend(post_alignment_filter_mapq_List)
post_alignment_List.extend(post_alignment_filter_duplicate_List)
post_alignment_List.extend(post_alignment_filter_chromosome_List)
post_alignment_List.extend(post_alignment_filter_tn5_List)
post_alignment_List.extend(post_alignment_pooling_List)
#
peak_calling_List = []
peak_calling_List.extend(peak_calling_macs2_narrowPeak_List)
peak_calling_List.extend(peak_calling_macs2_broadPeak_List)
#
genome_browser_List = []
genome_browser_List.extend(IGV_server_List)
genome_browser_List.extend(UCSC_server_List)
rule End_Point:
	input:
		pre_process_List
		+ alignment_List
		+ post_alignment_List
		+ peak_calling_List
		+ genome_browser_List


# ################################### PIPELINE RULES ###################################### #


#+++++++++++++++++++++++++++++
##REALM1: PRE_PROCESS
#+++++++++++++++++++++++++++++
rule Cutadapt:
	input:
		forward_fastq_List = get_forward_fastq,
		reverse_fastq_List = get_reverse_fastq,
	output:
		processed_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.discarded.fastq.gz",
	log:
		cutadapt_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.cutadapt.report",
	priority: 999
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Cutadapt: {wildcards.design}|{wildcards.sample}"
	run:

		for each_forward_fastq, each_reverse_fastq in zip(input.forward_fastq_List, input.reverse_fastq_List):
			#

			forward_fastq = os.path.basename(each_forward_fastq)
			reverse_fastq = os.path.basename(each_reverse_fastq)

			forward_ending = re.sub(".*" + "R1", "", forward_fastq)
			reverse_ending = re.sub(".*" + "R2", "", reverse_fastq)

			#
			processed_forward_fastq = forward_fastq.replace(forward_ending, ".processed.fastq.gz")
			processed_reverse_fastq = reverse_fastq.replace(reverse_ending, ".processed.fastq.gz")
			discarded_forward_fastq = forward_fastq.replace(forward_ending, ".discarded.fastq.gz")
			discarded_reverse_fastq = reverse_fastq.replace(reverse_ending, ".discarded.fastq.gz")
			cutadapt_stdout = forward_fastq.replace(forward_ending, ".cutadapt.stdout")
			cutadapt_stderr = forward_fastq.replace(forward_ending, ".cutadapt.stderr")
			#
			shell("""
				{ACTIVATE_CONDA_PY2}
				mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/
				processed_forward_fastq={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{processed_forward_fastq}
				processed_reverse_fastq={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{processed_reverse_fastq}
				discarded_forward_fastq={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{discarded_forward_fastq}
				discarded_reverse_fastq={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{discarded_reverse_fastq}
				cutadapt_stdout={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{cutadapt_stdout}
				cutadapt_stderr={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{cutadapt_stderr}

				##
				##EXECUTION
				start_time="$(date -u +%s)"
				cutadapt {config_pre_process_Dict[CUTADAPT_PE]} --too-short-output=$discarded_forward_fastq --too-short-paired-output=$discarded_reverse_fastq \
				--output=$processed_forward_fastq --paired-output=$processed_reverse_fastq {each_forward_fastq} {each_reverse_fastq} 1> $cutadapt_stdout 2> $cutadapt_stderr
				cat $discarded_forward_fastq >> {output.discarded_R1}
				cat $discarded_reverse_fastq >> {output.discarded_R2}
				cat $processed_forward_fastq >> {output.processed_R1}
				cat $processed_reverse_fastq >> {output.processed_R2}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#Report
				printf "%s\\t" "input:" > {log.cutadapt_report}
				printf "%s\\t" "{each_forward_fastq}" >> {log.cutadapt_report}
				printf "%s\\n" "{each_reverse_fastq}" >> {log.cutadapt_report}
				
				printf "%s\\t" "output:" >> {log.cutadapt_report}
				printf "%s\\t" "$processed_forward_fastq" >> {log.cutadapt_report}
				printf "%s\\t" "$processed_reverse_fastq" >> {log.cutadapt_report}
				printf "%s\\t" "$discarded_forward_fastq" >> {log.cutadapt_report}
				printf "%s\\n" "$discarded_reverse_fastq" >> {log.cutadapt_report}
				
				printf "%s\\t" "commandline:" >> {log.cutadapt_report}
				printf "%s\\n" "cutadapt {config_pre_process_Dict[CUTADAPT_PE]} --too-short-output=$discarded_forward_fastq --too-short-paired-output=$discarded_reverse_fastq --output=$processed_forward_fastq --paired-output=$processed_reverse_fastq {each_forward_fastq} {each_reverse_fastq}" >> {log.cutadapt_report}
				
				printf "%s\\t" "stdout:" >> {log.cutadapt_report}
				printf "%s\\n" "$cutadapt_stdout" >> {log.cutadapt_report}
				
				printf "%s\\t" "stderr:" >> {log.cutadapt_report}
				printf "%s\\n" "$cutadapt_stderr" >> {log.cutadapt_report}
				
				printf "%s\\t" "elapsed_time:" >> {log.cutadapt_report}
				printf "%s\\n" "$elapsed_time" >> {log.cutadapt_report}
				#
				#


			""")
		else:
			pass

rule FastQC:
	input:
		processed_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.discarded.fastq.gz",
	output:
		fastqc_processed_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed_fastqc.zip",
		fastqc_processed_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed_fastqc.zip",
		fastqc_discarded_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.discarded_fastqc.zip",
		fastqc_discarded_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.discarded_fastqc.zip",
	log:
		fastqc_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.fastqc.report",
	priority: 998
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "FastQC: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/
			fastq_file=$(basename {input.processed_R1})
		
			fastqc_html={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.html}}
			fastqc_zip={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.zip}}
			fastqc_stdout={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R1*/R1.processed.fastqc.stdout}}
			fastqc_stderr={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R1*/R1.processed.fastqc.stderr}}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.processed_R1} 1> $fastqc_stdout 2> $fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.fastqc_report}
			printf "%s\\n" "{input.processed_R1}" >> {log.fastqc_report}
			
			printf "%s\\t" "output:" >> {log.fastqc_report}
			printf "%s\\t" "$fastqc_html" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_zip" >> {log.fastqc_report}
			
			printf "%s\\t" "commandline:" >> {log.fastqc_report}
			printf "%s\\n" "fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.processed_R1}" >> {log.fastqc_report}
			
			printf "%s\\t" "stdout:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stdout" >> {log.fastqc_report}
			
			printf "%s\\t" "stderr:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stderr" >> {log.fastqc_report}
			
			printf "%s\\t" "elapsed_time:" >> {log.fastqc_report}
			printf "%s\\n" "$elapsed_time" >> {log.fastqc_report}
			#
			#


			fastq_file=$(basename {input.processed_R2})
		
			fastqc_html={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.html}}
			fastqc_zip={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.zip}}
			fastqc_stdout={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R2*/R2.processed.fastqc.stdout}}
			fastqc_stderr={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R2*/R2.processed.fastqc.stderr}}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.processed_R2} 1> $fastqc_stdout 2> $fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" >> {log.fastqc_report}
			printf "%s\\n" "{input.processed_R2}" >> {log.fastqc_report}
			
			printf "%s\\t" "output:" >> {log.fastqc_report}
			printf "%s\\t" "$fastqc_html" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_zip" >> {log.fastqc_report}
			
			printf "%s\\t" "commandline:" >> {log.fastqc_report}
			printf "%s\\n" "fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.processed_R1}" >> {log.fastqc_report}
			
			printf "%s\\t" "stdout:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stdout" >> {log.fastqc_report}
			
			printf "%s\\t" "stderr:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stderr" >> {log.fastqc_report}
			
			printf "%s\\t" "elapsed_time:" >> {log.fastqc_report}
			printf "%s\\n" "$elapsed_time" >> {log.fastqc_report}
			#
			#

			fastq_file=$(basename {input.discarded_R1})
		
			fastqc_html={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.html}}
			fastqc_zip={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.zip}}
			fastqc_stdout={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R1*/R1.discarded.fastqc.stdout}}
			fastqc_stderr={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R1*/R1.discarded.fastqc.stderr}}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.discarded_R1} 1> $fastqc_stdout 2> $fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" >> {log.fastqc_report}
			printf "%s\\n" "{input.discarded_R1}" >> {log.fastqc_report}
			
			printf "%s\\t" "output:" >> {log.fastqc_report}
			printf "%s\\t" "$fastqc_html" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_zip" >> {log.fastqc_report}
			
			printf "%s\\t" "commandline:" >> {log.fastqc_report}
			printf "%s\\n" "fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.discarded_R1}" >> {log.fastqc_report}
			
			printf "%s\\t" "stdout:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stdout" >> {log.fastqc_report}
			
			printf "%s\\t" "stderr:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stderr" >> {log.fastqc_report}
			
			printf "%s\\t" "elapsed_time:" >> {log.fastqc_report}
			printf "%s\\n" "$elapsed_time" >> {log.fastqc_report}
			#
			#

			fastq_file=$(basename {input.discarded_R2})
		
			fastqc_html={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.html}}
			fastqc_zip={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/.fastq*/_fastqc.zip}}
			fastqc_stdout={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R2*/R2.discarded.fastqc.stdout}}
			fastqc_stderr={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{fastq_file/R2*/R2.discarded.fastqc.stderr}}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.discarded_R2} 1> $fastqc_stdout 2> $fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" >> {log.fastqc_report}
			printf "%s\\n" "{input.discarded_R2}" >> {log.fastqc_report}
			

			printf "%s\\t" "output:" >> {log.fastqc_report}
			printf "%s\\t" "$fastqc_html" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_zip" >> {log.fastqc_report}
			
			printf "%s\\t" "commandline:" >> {log.fastqc_report}
			printf "%s\\n" "fastqc -o {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {input.discarded_R2}" >> {log.fastqc_report}
			
			printf "%s\\t" "stdout:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stdout" >> {log.fastqc_report}
			
			printf "%s\\t" "stderr:" >> {log.fastqc_report}
			printf "%s\\n" "$fastqc_stderr" >> {log.fastqc_report}
			
			printf "%s\\t" "elapsed_time:" >> {log.fastqc_report}
			printf "%s\\n" "$elapsed_time" >> {log.fastqc_report}
			#
			#
		""")

rule Subsample:
	input:
		processed_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/processed_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
	output:
		subsample_R1 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.subsample.fastq.gz",
		subsample_R2 = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.subsample.fastq.gz",
	log:
		#SUBSAMPLE
		subsample_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.subsample.report",
		subsample_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.subsample.stdout",
		subsample_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/subsample_fastq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.pre_process.subsample.stderr",
	priority: 997
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Subsample: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			mkdir -p {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/subsample_fastq/
			
			##
			##EXECUTION
			start_time="$(date -u +%s)"
			fastq_R1={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/subsample_fastq/{wildcards.sample}.R1.processed.fastq
			fastq_R2={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/subsample_fastq/{wildcards.sample}.R2.processed.fastq
			subsample_R1={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/subsample_fastq/{wildcards.sample}.R1.processed.subsample.fastq
			subsample_R2={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/pre_process/subsample_fastq/{wildcards.sample}.R2.processed.subsample.fastq
			module load usearch > {log.subsample_stdout} 2> {log.subsample_stderr}
			gunzip < {input.processed_R1} > $fastq_R1 2>> {log.subsample_stderr}
			gunzip < {input.processed_R2} > $fastq_R2 2>> {log.subsample_stderr}
			usearch -fastx_subsample $fastq_R1 -reverse $fastq_R2 -fastqout $subsample_R1 -output2 $subsample_R2 -sample_pct {SUBSAMPLE_PCT} >> {log.subsample_stdout} 2>> {log.subsample_stderr}
			rm -f $fastq_R1 $fastq_R2
			gzip $subsample_R1 2>> {log.subsample_stderr}
			gzip $subsample_R2 2>> {log.subsample_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.subsample_report}
			printf "%s\\t" "{input.processed_R1}" >> {log.subsample_report}
			printf "%s\\n" "{input.processed_R2}" >> {log.subsample_report}
			
			printf "%s\\t" "output:" >> {log.subsample_report}
			printf "%s\\t" "{output.subsample_R1}" >> {log.subsample_report}
			printf "%s\\n" "{output.subsample_R2}" >> {log.subsample_report}
			
			printf "%s\\t" "commandline:" >> {log.subsample_report}
			printf "%s\\t" "gunzip < {input.processed_R1} > $fastq_R1" >> {log.subsample_report}
			printf "%s\\t" "gunzip < {input.processed_R2} > $fastq_R2" >> {log.subsample_report}
			printf "%s\\t" "usearch -fastx_subsample $fastq_R1 -reverse $fastq_R2 -fastqout $subsample_R1 -output2 $subsample_R2 -sample_pct {SUBSAMPLE_PCT}" >> {log.subsample_report}
			printf "%s\\t" "rm -f $fastq_R1 $fastq_R2" >> {log.subsample_report}
			printf "%s\\t" "gzip $subsample_R1" >> {log.subsample_report}
			printf "%s\\n" "gzip $subsample_R2" >> {log.subsample_report}
			
			printf "%s\\t" "stdout:" >> {log.subsample_report}
			printf "%s\\n" "{log.subsample_stdout}" >> {log.subsample_report}
			
			printf "%s\\t" "stderr:" >> {log.subsample_report}
			printf "%s\\n" "{log.subsample_stderr}" >> {log.subsample_report}
			
			printf "%s\\t" "elapsed_time:" >> {log.subsample_report}
			printf "%s\\n" "$elapsed_time" >> {log.subsample_report}
			#
			#

			""")

#+++++++++++++++++++++++++++++
##REALM2: ALIGNMENT
#+++++++++++++++++++++++++++++

rule Bowtie2:
	input:
		processed_R1 = get_processed_forward_fastq,
		processed_R2 = get_processed_reverse_fastq
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
			
			mapped_R1={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.mapped.fastq.gz
			mapped_R2={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.mapped.fastq.gz
			unmapped_R1={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.unmapped.fastq.gz
			unmapped_R2={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.unmapped.fastq.gz
			
			##
			##Execution
			start_time="$(date -u +%s)"
			bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_R1} -2 {input.processed_R2} \
			--un-conc-gz {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.unmapped.fastq.gz \
			--al-conc-gz {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.mapped.fastq.gz \
			2> {log.bowtie2_stdout} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample} -o {output.bowtie2_Bam} - 1> {log.bowtie2_stdout} 2> {log.bowtie2_stderr}
			samtools index -@ {threads} -b {output.bowtie2_Bam} >> {log.bowtie2_stdout} 2>> {log.bowtie2_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.bowtie2_report}
			printf "%s\\t" "{input.processed_R1}" >> {log.bowtie2_report}
			printf "%s\\n" "{input.processed_R2}" >> {log.bowtie2_report}
			
			printf "%s\\t" "output:" >> {log.bowtie2_report}
			printf "%s\\t" "$mapped_R1" >> {log.bowtie2_report}
			printf "%s\\t" "$mapped_R2" >> {log.bowtie2_report}
			printf "%s\\t" "$unmapped_R1" >> {log.bowtie2_report}
			printf "%s\\t" "$unmapped_R2" >> {log.bowtie2_report}
			printf "%s\\t" "{output.bowtie2_Bam}" >> {log.bowtie2_report}
			printf "%s\\n" "{output.bowtie2_index_Bam}" >> {log.bowtie2_report}
			
			printf "%s\\t" "commandline:" >> {log.bowtie2_report}
			printf "%s\\t" "bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_R1} -2 {input.processed_R2} --un-conc-gz {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.unmapped.fastq.gz --al-conc-gz {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.mapped.fastq.gz" >> {log.bowtie2_report}
			printf "%s\\t" "samtools sort --threads {threads} -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample} -o {output.bowtie2_Bam}" >> {log.bowtie2_report}
			printf "%s\\n" "samtools index -@ {threads} -b {output.bowtie2_Bam}" >> {log.bowtie2_report}
			
			printf "%s\\t" "stdout:" >> {log.bowtie2_report}
			printf "%s\\n" "{log.bowtie2_stdout}" >> {log.bowtie2_report}
			
			printf "%s\\t" "stderr:" >> {log.bowtie2_report}
			printf "%s\\n" "{log.bowtie2_stderr}" >> {log.bowtie2_report}

			printf "%s\\t" "elapsed_time:" >> {log.bowtie2_report}
			printf "%s\\n" "$elapsed_time" >> {log.bowtie2_report}
			#
			#
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
			export _JAVA_OPTIONS="-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}"
			picard MarkDuplicates INPUT={input.bowtie2_Bam} OUTPUT={output.picard_Bam} METRICS_FILE={log.picard_stdout} {config_alignment_Dict[PICARD]} > {log.picard_stderr}
			samtools index -@ {threads} -b {output.picard_Bam} &>> {log.picard_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.picard_report}
			printf "%s\\n" "{input.bowtie2_Bam}" >> {log.picard_report}
			
			printf "%s\\t" "output:" >> {log.picard_report}
			printf "%s\\t" "{output.picard_Bam}" >> {log.picard_report}
			printf "%s\\n" "{output.picard_index_Bam}" >> {log.picard_report}

			printf "%s\\t" "commandline:" >> {log.picard_report}
			printf "%s\\t" "export _JAVA_OPTIONS='-Xms{MEMORY}M -Xmx{MEMORY}M -XX:ParallelGCThreads={threads}'" >> {log.picard_report}
			printf "%s\\t" "picard MarkDuplicates INPUT={input.bowtie2_Bam} OUTPUT={output.picard_Bam} METRICS_FILE={log.picard_stdout} {config_alignment_Dict[PICARD]}" >> {log.picard_report}
			printf "%s\\n" "samtools index -@ {threads} -b {output.picard_Bam}" >> {log.picard_report}
			
			printf "%s\\t" "stdout:" >> {log.picard_report}
			printf "%s\\n" "{log.picard_stdout}" >> {log.picard_report}
			
			printf "%s\\t" "stderr:" >> {log.picard_report}
			printf "%s\\n" "{log.picard_stderr}" >> {log.picard_report}

			printf "%s\\t" "elapsed_time:" >> {log.picard_report}
			printf "%s\\n" "$elapsed_time" >> {log.picard_report}
			#
			#
		""")

rule Qualimap:
	input:
		picard_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam",
		picard_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.bam.bai",
	output:
		qualimap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}/qualimapReport.html",
	log:
		qualimap_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.qualimap.report",
		qualimap_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.qualimap.stdout",
		qualimap_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.qualimap.stderr",
	priority: 994
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Qualimap: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			qualimap bamqc -bam {input.picard_Bam} -nt {threads} -outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/qualimap/{wildcards.sample} {config_alignment_Dict[QUALIMAP]} 1> {log.qualimap_stdout} 2> {log.qualimap_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.qualimap_report}
			printf "%s\\n" "{input.picard_Bam}" >> {log.qualimap_report}
			
			printf "%s\\t" "output:" >> {log.qualimap_report}
			printf "%s\\n" "{output.qualimap}" >> {log.qualimap_report}
			
			printf "%s\\t" "commandline:" >> {log.qualimap_report}
			printf "%s\\n" "qualimap bamqc -bam {input.picard_Bam} -nt {threads} -outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/qualimap/{wildcards.sample} {config_alignment_Dict[QUALIMAP]}" >> {log.qualimap_report}
			
			printf "%s\\t" "stdout:" >> {log.qualimap_report}
			printf "%s\\n" "{log.qualimap_stdout}" >> {log.qualimap_report}
			
			printf "%s\\t" "stderr:" >> {log.qualimap_report}
			printf "%s\\n" "{log.qualimap_stderr}" >> {log.qualimap_report}

			printf "%s\\t" "elapsed_time:" >> {log.qualimap_report}
			printf "%s\\n" "$elapsed_time" >> {log.qualimap_report}
			#
			#
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
		filter_mapq_discarded_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq_discarded.bam",
		filter_mapq_discarded_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq_discarded.bam.report",
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
			samtools sort -n --threads {threads} {output.filter_mapq_Bam} -o {output.filter_mapq_Bam}.nmsort.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools fixmate --threads {threads} {output.filter_mapq_Bam}.nmsort.tmp {output.filter_mapq_Bam}.nmsort.fixmate.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools sort --threads {threads} {output.filter_mapq_Bam}.nmsort.fixmate.tmp -o {output.filter_mapq_Bam} 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools index -@ {threads} -b {output.filter_mapq_Bam} &>> {log.filter_mapq_stderr}
			samtools flagstat --threads {threads} {output.filter_mapq_Bam} > {output.filter_mapq_Bam_Report} 2>> {log.filter_mapq_stderr}
			rm -f {output.filter_mapq_Bam}.nmsort.tmp {output.filter_mapq_Bam}.nmsort.fixmate.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			
			samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_MAPQ_DISCARDED]} {input.picard_Bam} -b -O BAM -o {output.filter_mapq_discarded_Bam} 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools sort -n --threads {threads} {output.filter_mapq_discarded_Bam} -o {output.filter_mapq_discarded_Bam}.nmsort.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools fixmate --threads {threads} {output.filter_mapq_discarded_Bam}.nmsort.tmp {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools sort --threads {threads} {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp -o {output.filter_mapq_discarded_Bam} 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			samtools index -@ {threads} -b {output.filter_mapq_discarded_Bam} &>> {log.filter_mapq_stderr}
			samtools flagstat --threads {threads} {output.filter_mapq_discarded_Bam} > {output.filter_mapq_discarded_Bam_Report} 2>> {log.filter_mapq_stderr}
			rm -f {output.filter_mapq_discarded_Bam}.nmsort.tmp {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp 1>> {log.filter_mapq_stdout} 2>> {log.filter_mapq_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.filter_mapq_report}
			printf "%s\\n" "{input.picard_Bam}" >> {log.filter_mapq_report}
			
			printf "%s\\t" "output:" >> {log.filter_mapq_report}
			printf "%s\\t" "{output.filter_mapq_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "{output.filter_mapq_Bam_Report}" >> {log.filter_mapq_report}
			printf "%s\\t" "{output.filter_mapq_discarded_Bam}" >> {log.filter_mapq_report}
			printf "%s\\n" "{output.filter_mapq_discarded_Bam_Report}" >> {log.filter_mapq_report}

			printf "%s\\t" "commandline:" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_MAPQ]} {input.picard_Bam} -b -O BAM -o {output.filter_mapq_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.filter_mapq_Bam} -o {output.filter_mapq_Bam}.nmsort" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_mapq_Bam}.nmsort -O {output.filter_mapq_Bam}.nmsort.fixmate.tmp" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.filter_mapq_Bam}.nmsort.fixmate.tmp -o {output.filter_mapq_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.filter_mapq_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.filter_mapq_Bam} > {output.filter_mapq_Bam_Report}" >> {log.filter_mapq_report}
			printf "%s\\t" "rm -f {output.filter_mapq_Bam}.nmsort.tmp {output.filter_mapq_Bam}.nmsort.fixmate.tmp" >> {log.filter_mapq_report}
			
			printf "%s\\t" "samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_MAPQ_DISCARDED]} {input.picard_Bam} -b -O BAM -o {output.filter_mapq_discarded_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.filter_mapq_discarded_Bam} -o {output.filter_mapq_discarded_Bam}.nmsort" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_mapq_discarded_Bam}.nmsort -O {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp -o {output.filter_mapq_discarded_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.filter_mapq_discarded_Bam}" >> {log.filter_mapq_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.filter_mapq_discarded_Bam} > {output.filter_mapq_discarded_Bam_Report}" >> {log.filter_mapq_report}
			printf "%s\\n" "rm -f {output.filter_mapq_discarded_Bam}.nmsort.tmp {output.filter_mapq_discarded_Bam}.nmsort.fixmate.tmp" >> {log.filter_mapq_report}

			printf "%s\\t" "stdout:" >> {log.filter_mapq_report}
			printf "%s\\n" "{log.filter_mapq_stdout}" >> {log.filter_mapq_report}
			
			printf "%s\\t" "stderr:" >> {log.filter_mapq_report}
			printf "%s\\n" "{log.filter_mapq_stderr}" >> {log.filter_mapq_report}

			printf "%s\\t" "elapsed_time:" >> {log.filter_mapq_report}
			printf "%s\\n" "$elapsed_time" >> {log.filter_mapq_report}
			#
			#
		""")

rule Filter_Duplicate:
	input:
		filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.bam",
	output:
		noduplicate_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.bam",
		noduplicate_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.bam.report",
		duplicate_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.duplicate.bam",
		duplicate_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.duplicate.bam.report",
	log:
		filter_duplicate_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_duplicate.report",
		filter_duplicate_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_duplicate.stdout",
		filter_duplicate_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_duplicate.stderr",
	priority: 992
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Filter_Duplicate: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##EXECUTION
			start_time="$(date -u +%s)"
			samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_OUT_DUPLICATE]} {input.filter_mapq_Bam} -b -O BAM -o {output.noduplicate_Bam} 1> {log.filter_duplicate_stdout} 2> {log.filter_duplicate_stderr}
			samtools sort -n --threads {threads} {output.noduplicate_Bam} -o {output.noduplicate_Bam}.nmsort.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools fixmate --threads {threads} {output.noduplicate_Bam}.nmsort.tmp {output.noduplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools sort --threads {threads} {output.noduplicate_Bam}.nmsort.fixmate.tmp -o {output.noduplicate_Bam} 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools flagstat --threads {threads} {output.noduplicate_Bam} > {output.noduplicate_Bam_Report} 2>> {log.filter_duplicate_stderr}
			samtools index -@ {threads} -b {output.noduplicate_Bam} &>> {log.filter_duplicate_stderr}
			rm -f {output.noduplicate_Bam}.nmsort.tmp {output.noduplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			
			samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_IN_DUPLICATE]} {input.filter_mapq_Bam} -b -O BAM -o {output.duplicate_Bam} 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools sort -n --threads {threads} {output.duplicate_Bam} -o {output.duplicate_Bam}.nmsort.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools fixmate --threads {threads} {output.duplicate_Bam}.nmsort.tmp {output.duplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools sort --threads {threads} {output.duplicate_Bam}.nmsort.fixmate.tmp -o {output.duplicate_Bam} 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			samtools flagstat --threads {threads} {output.duplicate_Bam} > {output.duplicate_Bam_Report} 2>> {log.filter_duplicate_stderr}
			samtools index -@ {threads} -b {output.duplicate_Bam} &>> {log.filter_duplicate_stderr}
			rm -f {output.duplicate_Bam}.nmsort.tmp {output.duplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_duplicate_stdout} 2>> {log.filter_duplicate_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.filter_duplicate_report}
			printf "%s\\n" "{input.filter_mapq_Bam}" >> {log.filter_duplicate_report}
			
			printf "%s\\t" "output:" >> {log.filter_duplicate_report}
			printf "%s\\t" "{output.noduplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "{output.noduplicate_Bam_Report}" >> {log.filter_duplicate_report}
			printf "%s\\t" "{output.duplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\n" "{output.duplicate_Bam_Report}" >> {log.filter_duplicate_report}

			printf "%s\\t" "commandline:" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_OUT_DUPLICATE]} {input.filter_mapq_Bam} -b -O BAM -o {output.noduplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.noduplicate_Bam} -o {output.noduplicate_Bam}.nmsort.tmp" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.noduplicate_Bam}.nmsort.tmp {output.noduplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.noduplicate_Bam}.nmsort.fixmate.tmp -o {output.noduplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.noduplicate_Bam} > {output.noduplicate_Bam_Report}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.noduplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "rm -f {output.noduplicate_Bam}.nmsort.tmp {output.noduplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_duplicate_report}
			
			printf "%s\\t" "samtools view --threads {threads} -h {config_post_alignment_Dict[FILTER_IN_DUPLICATE]} {input.filter_mapq_Bam} -b -O BAM -o {output.duplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.duplicate_Bam} -o {output.duplicate_Bam}.nmsort.tmp" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.duplicate_Bam}.nmsort.tmp {output.duplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.duplicate_Bam}.nmsort.fixmate.tmp -o {output.duplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.duplicate_Bam} > {output.duplicate_Bam_Report}" >> {log.filter_duplicate_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.duplicate_Bam}" >> {log.filter_duplicate_report}
			printf "%s\\n" "rm -f {output.duplicate_Bam}.nmsort.tmp {output.duplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_duplicate_report}

			printf "%s\\t" "stdout:" >> {log.filter_duplicate_report}
			printf "%s\\n" "{log.filter_duplicate_stdout}" >> {log.filter_duplicate_report}
			
			printf "%s\\t" "stderr:" >> {log.filter_duplicate_report}
			printf "%s\\n" "{log.filter_duplicate_stderr}" >> {log.filter_duplicate_report}

			printf "%s\\t" "elapsed_time:" >> {log.filter_duplicate_report}
			printf "%s\\n" "$elapsed_time" >> {log.filter_duplicate_report}
			#
			#
		""")

rule Filter_Chromosome:
	input:
		filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_mapq/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.bam",
		noduplicate_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_duplicate/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.bam",
	output:
		filter_chromosome_filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam",
		filter_chromosome_filter_mapq_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam.report",
		filter_chromosome_noduplicate_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.bam",
		filter_chromosome_noduplicate_Bam_Report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.bam.report",
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
			samtools sort -n --threads {threads} {output.filter_chromosome_filter_mapq_Bam} -o {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp 1> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools fixmate --threads {threads} {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools sort --threads {threads} {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp -o {output.filter_chromosome_filter_mapq_Bam} 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools index -@ {threads} -b {output.filter_chromosome_filter_mapq_Bam} &>> {log.filter_chromosome_stderr}
			samtools flagstat --threads {threads} {output.filter_chromosome_filter_mapq_Bam} > {output.filter_chromosome_filter_mapq_Bam_Report} 2>> {log.filter_chromosome_stderr}
			rm -f {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}

			samtools view --threads {threads} -h {input.noduplicate_Bam} | awk -F'\\t' '{CHROMOSOME_FILTER_PROCESS}' | samtools view --threads {threads} -Shb - > {output.filter_chromosome_noduplicate_Bam} 2>> {log.filter_chromosome_stderr}
			samtools sort -n --threads {threads} {output.filter_chromosome_noduplicate_Bam} -o {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools fixmate --threads {threads} {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools sort --threads {threads} {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp -o {output.filter_chromosome_noduplicate_Bam} 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			samtools index -@ {threads} -b {output.filter_chromosome_noduplicate_Bam} &>> {log.filter_chromosome_stderr}
			samtools flagstat --threads {threads} {output.filter_chromosome_noduplicate_Bam} > {output.filter_chromosome_noduplicate_Bam_Report} 2>> {log.filter_chromosome_stderr}
			rm -f {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp 1>> {log.filter_chromosome_stdout} 2>> {log.filter_chromosome_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			CHROMOSOME_FILTER_COMMAND='{config_post_alignment_Dict[FILTER_CHROMOSOME]}'

			#
			#Report
			printf "%s\\t" "input:" > {log.filter_chromosome_report}
			printf "%s\\t" "{input.filter_mapq_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\n" "{input.noduplicate_Bam}" >> {log.filter_chromosome_report}
			
			printf "%s\\t" "output:" >> {log.filter_chromosome_report}
			printf "%s\\t" "{output.filter_chromosome_filter_mapq_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "{output.filter_chromosome_filter_mapq_Bam_Report}" >> {log.filter_chromosome_report}
			printf "%s\\t" "{output.filter_chromosome_noduplicate_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\n" "{output.filter_chromosome_noduplicate_Bam_Report}" >> {log.filter_chromosome_report}

			printf "%s\\t" "commandline:" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools view --threads {threads} -h {input.filter_mapq_Bam} | awk -F'\\t' '$CHROMOSOME_FILTER_COMMAND' | samtools view --threads {threads} -Shb - > {output.filter_chromosome_filter_mapq_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.filter_chromosome_filter_mapq_Bam} -o {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp -o {output.filter_chromosome_filter_mapq_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.filter_chromosome_filter_mapq_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.filter_chromosome_filter_mapq_Bam} > {output.filter_chromosome_filter_mapq_Bam_Report}"  >> {log.filter_chromosome_report}
			printf "%s\\t" "rm -f {output.filter_chromosome_filter_mapq_Bam}.nmsort.tmp {output.filter_chromosome_filter_mapq_Bam}.nmsort.fixmate.tmp" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools view --threads {threads} -h {input.noduplicate_Bam} | awk -F'\\t' '$CHROMOSOME_FILTER_COMMAND' | samtools view --threads {threads} -Shb - > {output.filter_chromosome_noduplicate_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {output.filter_chromosome_noduplicate_Bam} -o {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools sort --threads {threads} {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp -o {output.filter_chromosome_noduplicate_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools index -@ {threads} -b {output.filter_chromosome_noduplicate_Bam}" >> {log.filter_chromosome_report}
			printf "%s\\t" "samtools flagstat --threads {threads} {output.filter_chromosome_noduplicate_Bam} > {output.filter_chromosome_noduplicate_Bam_Report}" >> {log.filter_chromosome_report}
			printf "%s\\n" "rm -f {output.filter_chromosome_noduplicate_Bam}.nmsort.tmp {output.filter_chromosome_noduplicate_Bam}.nmsort.fixmate.tmp" >> {log.filter_chromosome_report}

			printf "%s\\t" "stdout:" >> {log.filter_chromosome_report}
			printf "%s\\n" "{log.filter_chromosome_stdout}" >> {log.filter_chromosome_report}
			
			printf "%s\\t" "stderr:" >> {log.filter_chromosome_report}
			printf "%s\\n" "{log.filter_chromosome_stderr}" >> {log.filter_chromosome_report}

			printf "%s\\t" "elapsed_time:" >> {log.filter_chromosome_report}
			printf "%s\\n" "$elapsed_time" >> {log.filter_chromosome_report}
			#
			#

		""")

rule Filter_Tn5:
	input:
		filter_chromosome_filter_mapq_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.bam",
		filter_chromosome_noduplicate_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_chromosome/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.bam",
	output:
		filter_tn5_filter_chromosome_filter_mapq_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz",
		filter_tn5_filter_chromosome_noduplicate_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz",
	log:
		filter_tn5_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_tn5.report",
		filter_tn5_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_tn5.stdout",
		filter_tn5_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.filter_tn5.stderr",

	priority: 990
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Filter_Tn5: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			#EXECUTION
			start_time="$(date -u +%s)"
			samtools sort -n --threads {threads} {input.filter_chromosome_filter_mapq_Bam} -o {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp 1> {log.filter_tn5_stdout} 2> {log.filter_tn5_stderr}
			samtools fixmate --threads {threads} {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			bedtools bamtobed -cigar -i {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp 2>> {log.filter_tn5_stderr}
			cat {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp | awk -F'\\t' '{AWK_TN5_PROCESS}' > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp 2>> {log.filter_tn5_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted 2>> {log.filter_tn5_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted > {output.filter_tn5_filter_chromosome_filter_mapq_Bed} 2>> {log.filter_tn5_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_filter_mapq_Bed} 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}

			samtools sort -n --threads {threads} {input.filter_chromosome_noduplicate_Bam} -o {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp 1> {log.filter_tn5_stdout} 2> {log.filter_tn5_stderr}
			samtools fixmate --threads {threads} {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			bedtools bamtobed -cigar -i {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp 2>> {log.filter_tn5_stderr}
			cat {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp | awk -F'\\t' '{AWK_TN5_PROCESS}' > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp 2>> {log.filter_tn5_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted 2>> {log.filter_tn5_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted > {output.filter_tn5_filter_chromosome_noduplicate_Bed} 2>> {log.filter_tn5_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_Bed} 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted 1>> {log.filter_tn5_stdout} 2>> {log.filter_tn5_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			TN5_AWK_COMMAND='{config_post_alignment_Dict[FILTER_TN5]}'

			#
			#Report
			printf "%s\\t" "input:" > {log.filter_tn5_report}
			printf "%s\\t" "{input.filter_chromosome_filter_mapq_Bam}" >> {log.filter_tn5_report}
			printf "%s\\n" "{input.filter_chromosome_noduplicate_Bam}" >> {log.filter_tn5_report}

			printf "%s\\t" "output:" >> {log.filter_tn5_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_filter_mapq_Bed}" >> {log.filter_tn5_report}
			printf "%s\\n" "{output.filter_tn5_filter_chromosome_noduplicate_Bed}" >> {log.filter_tn5_report}

			printf "%s\\t" "commandline:" >> {log.filter_tn5_report}
			printf "%s\\t" "samtools sort -n --threads {threads} {input.filter_chromosome_filter_mapq_Bam} -o {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "bedtools bamtobed -cigar -i {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "cat {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp | awk -F'\\t' '$TN5_AWK_COMMAND' > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted" >> {log.filter_tn5_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted > {output.filter_tn5_filter_chromosome_filter_mapq_Bed}" >> {log.filter_tn5_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_filter_mapq_Bed}" >> {log.filter_tn5_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.bed.tmp {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_filter_mapq_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted" >> {log.filter_tn5_report}

			printf "%s\\t" "samtools sort -n --threads {threads} {input.filter_chromosome_noduplicate_Bam} -o {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "samtools fixmate --threads {threads} {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "bedtools bamtobed -cigar -i {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "cat {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp | awk -F'\\t' '$TN5_AWK_COMMAND' > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp > {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted" >> {log.filter_tn5_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted > {output.filter_tn5_filter_chromosome_noduplicate_Bed}" >> {log.filter_tn5_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_Bed}" >> {log.filter_tn5_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tmp" >> {log.filter_tn5_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.bed.tmp {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp" >> {log.filter_tn5_report}
			printf "%s\\n" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_Bed}.nmsort.fixmate.tn5.bed.tmp.sorted" >> {log.filter_tn5_report}

			printf "%s\\t" "stdout:" >> {log.filter_tn5_report}
			printf "%s\\n" "{log.filter_tn5_stdout}" >> {log.filter_tn5_report}
			
			printf "%s\\t" "stderr:" >> {log.filter_tn5_report}
			printf "%s\\n" "{log.filter_tn5_stderr}" >> {log.filter_tn5_report}

			printf "%s\\t" "elapsed_time:" >> {log.filter_tn5_report}
			printf "%s\\n" "$elapsed_time" >> {log.filter_tn5_report}
			#
			#

		""")

rule Pooling:
	input:
		filter_tn5_filter_chromosome_Bed_List = get_filter_tn5_filter_chromosome_Bed,
		filter_tn5_filter_chromosome_noduplicate_Bed_List = get_filter_tn5_filter_chromosome_noduplicate_Bed,
	output:
		pooling_filter_tn5_filter_chromosome_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case, .*_POOLED_.*}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz",
		pooling_filter_tn5_noduplicate_filter_chromosome_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case, .*_POOLED_.*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz",
	log:
		pooling_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.report",
		pooling_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stdout",
		pooling_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stderr",
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
			cat {input.filter_tn5_filter_chromosome_Bed_List} > {output.pooling_filter_tn5_filter_chromosome_Bed} 2>> {log.pooling_stderr}
			gunzip < {output.pooling_filter_tn5_filter_chromosome_Bed} > {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp 2>> {log.pooling_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp > {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp 2>> {log.pooling_stderr}
			bgzip -c {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp > {output.pooling_filter_tn5_filter_chromosome_Bed} 2>> {log.pooling_stderr}
			tabix -f -p bed {output.pooling_filter_tn5_filter_chromosome_Bed} 1> {log.pooling_stdout} 2>> {log.pooling_stderr}
			rm -f {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp 1>> {log.pooling_stdout} 2>> {log.pooling_stderr}

			cat {input.filter_tn5_filter_chromosome_noduplicate_Bed_List} > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed} 2>> {log.pooling_stderr}
			gunzip < {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed} > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp 2>> {log.pooling_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp 2>> {log.pooling_stderr}
			bgzip -c {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed} 2>> {log.pooling_stderr}
			tabix -f -p bed {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed} 1>> {log.pooling_stdout} 2>> {log.pooling_stderr}
			rm -f {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp 1>> {log.pooling_stdout} 2>> {log.pooling_stderr}

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s" "input:" > {log.pooling_report}
			declare -a bam_List=({input.filter_tn5_filter_chromosome_Bed_List})

			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "\\t%s" "${{bam_List[$case_index]}}" >> {log.pooling_report}
			done
			declare -a bam_List=({input.filter_tn5_filter_chromosome_noduplicate_Bed_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($index+1))
				printf "\\t%s" "${{bam_List[$case_index]}}" >> {log.pooling_report}
			done

			printf "\\n" >> {log.pooling_report}

			printf "%s\\t" "output:" >> {log.pooling_report}
			printf "%s\\t" "{output.pooling_filter_tn5_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\n" "{output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}" >> {log.pooling_report}

			printf "%s\\t" "commandline:" >> {log.pooling_report}
			printf "%s\\t" "cat {input.filter_tn5_filter_chromosome_Bed_List} > {output.pooling_filter_tn5_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "gunzip < {output.pooling_filter_tn5_filter_chromosome_Bed} > {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp" >> {log.pooling_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp > {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp" >> {log.pooling_report}
			printf "%s\\t" "bgzip -c {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp > {output.pooling_filter_tn5_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "tabix -f -p bed {output.pooling_filter_tn5_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "rm -f {output.pooling_filter_tn5_filter_chromosome_Bed}.tmp {output.pooling_filter_tn5_filter_chromosome_Bed}.nmsort.tmp" >> {log.pooling_report}
			
			printf "%s\\t" "cat {input.filter_tn5_filter_chromosome_noduplicate_Bed_List} > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "gunzip < {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed} > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp" >> {log.pooling_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp" >> {log.pooling_report}
			printf "%s\\t" "bgzip -c {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp > {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "tabix -f -p bed {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}" >> {log.pooling_report}
			printf "%s\\t" "rm -f {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.tmp {output.pooling_filter_tn5_noduplicate_filter_chromosome_Bed}.nmsort.tmp" >> {log.pooling_report}
			
			printf "%s\\t" "stdout:" >> {log.pooling_report}
			printf "%s\\n" "{log.pooling_stdout}" >> {log.pooling_report}
			
			printf "%s\\t" "stderr:" >> {log.pooling_report}
			printf "%s\\n" "{log.pooling_stderr}" >> {log.pooling_report}

			printf "%s\\t" "elapsed_time:" >> {log.pooling_report}
			printf "%s\\n" "$elapsed_time" >> {log.pooling_report}
			#
			#
		""")


#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++


rule Macs2_NarrowPeak:
	input:
		filter_tn5_filter_chromosome_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz",
		filter_tn5_filter_chromosome_noduplicate_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz",
	output:
		filter_tn5_filter_chromosome_macs2_narrowPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.narrowPeak.gz",
		filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bedgraph.gz",
		filter_tn5_filter_chromosome_macs2_narrowPeak_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_narrowPeak.bigwig",
		#
		filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_narrowPeak.narrowPeak.gz",
		filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_narrowPeak.bedgraph.gz",
		filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_narrowPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_narrowPeak.bigwig",
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
			sample_Name=$(basename {input.filter_tn5_filter_chromosome_Bed})
			sample_Name=${{sample_Name%.bed.gz}}
			
			macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_Bed} --name ${{sample_Name}}.macs2_narrow {config_peak_calling_Dict[MACS2_NARROWPEAK]} \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ 1> {log.macs2_narrowPeak_stdout} 2> {log.macs2_narrowPeak_stderr}
			raw_peak_count=$(wc -l {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak | awk '{{printf "%f", $1/1000000}}')
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bedtools intersect -v -a {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp 2>> {log.macs2_narrowPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed} 2>> {log.macs2_narrowPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}

			macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_treat_pileup.bdg \
			-c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ -m FE 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph} 2>> {log.macs2_narrowPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			bedGraphToBigWig {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bigwig} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}

			sample_Name=$(basename {input.filter_tn5_filter_chromosome_noduplicate_Bed})
			sample_Name=${{sample_Name%.bed.gz}}
			
			macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_noduplicate_Bed} --name ${{sample_Name}}.macs2_narrow {config_peak_calling_Dict[MACS2_NARROWPEAK]} \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			raw_peak_count=$(wc -l {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak | awk '{{printf "%f", $1/1000000}}')
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bedtools intersect -v -a {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp 2>> {log.macs2_narrowPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed} 2>> {log.macs2_narrowPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}

			macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_treat_pileup.bdg \
			-c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ -m FE 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp 2>> {log.macs2_narrowPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph} 2>> {log.macs2_narrowPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			bedGraphToBigWig {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bigwig} 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp 1>> {log.macs2_narrowPeak_stdout} 2>> {log.macs2_narrowPeak_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.macs2_narrowPeak_report}
			printf "%s\\t" "{input.filter_tn5_filter_chromosome_Bed}" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "{input.filter_tn5_filter_chromosome_noduplicate_Bed}" >> {log.macs2_narrowPeak_report}
			
			printf "%s\\t" "output:" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}">> {log.macs2_narrowPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bigwig}">> {log.macs2_narrowPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bigwig}" >> {log.macs2_narrowPeak_report}
			
			printf "%s\\t" "commandline:" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_Bed} --name ${{sample_Name}}.macs2_narrow {config_peak_calling_Dict[MACS2_NARROWPEAK]} --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bedtools intersect -v -a {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ -m FE " >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bedGraphToBigWig {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bigwig}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_macs2_narrowPeak_Bedgraph}.sorted.tmp" >> {log.macs2_narrowPeak_report}

			printf "%s\\t" "macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_noduplicate_Bed} --name ${{sample_Name}}.macs2_narrow {config_peak_calling_Dict[MACS2_NARROWPEAK]} --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bedtools intersect -v -a {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/ -m FE " >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_narrowPeak/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}" >> {log.macs2_narrowPeak_report}
			printf "%s\\t" "bedGraphToBigWig {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bigwig}" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_narrowPeak_Bedgraph}.sorted.tmp" >> {log.macs2_narrowPeak_report}


			printf "%s\\t" "stdout:" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "{log.macs2_narrowPeak_stdout}" >> {log.macs2_narrowPeak_report}
			
			printf "%s\\t" "stderr:" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "{log.macs2_narrowPeak_stderr}" >> {log.macs2_narrowPeak_report}

			printf "%s\\t" "elapsed_time:" >> {log.macs2_narrowPeak_report}
			printf "%s\\n" "$elapsed_time" >> {log.macs2_narrowPeak_report}
			#
			#
		""")


rule Macs2_BroadPeak:
	input:
		filter_tn5_filter_chromosome_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.bed.gz",
		filter_tn5_filter_chromosome_noduplicate_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/filter_tn5/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.bed.gz",
	output:
		filter_tn5_filter_chromosome_macs2_broadPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.broadPeak.gz",
		filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bedgraph.gz",
		filter_tn5_filter_chromosome_macs2_broadPeak_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.filter_chromosome.filter_tn5.macs2_broadPeak.bigwig",
		#
		filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_broadPeak.broadPeak.gz",
		filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_broadPeak.bedgraph.gz",
		filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.picard.filter_mapq.noduplicate.filter_chromosome.filter_tn5.macs2_broadPeak.bigwig",
	log:
		macs2_broadPeak_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_broadPeak.report",
		macs2_broadPeak_stdout = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_broadPeak.stdout",
		macs2_broadPeak_stderr = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/macs2_broadPeak/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.macs2_broadPeak.stderr",
	priority: 989
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Macs2_BroadPeak: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			##
			#EXECUTION
			start_time="$(date -u +%s)"
			sample_Name=$(basename {input.filter_tn5_filter_chromosome_Bed})
			sample_Name=${{sample_Name%.bed.gz}}
			
			macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_Bed} --name ${{sample_Name}}.macs2_broad {config_peak_calling_Dict[MACS2_BROADPEAK]} \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ 1> {log.macs2_broadPeak_stdout} 2> {log.macs2_broadPeak_stderr}
			raw_peak_count=$(wc -l {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak | awk '{{printf "%f", $1/1000000}}')
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bedtools intersect -v -a {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp 2>> {log.macs2_broadPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed} 2>> {log.macs2_broadPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}

			macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_treat_pileup.bdg \
			-c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ -m ppois -S $raw_peak_count 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph} 2>> {log.macs2_broadPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			bedGraphToBigWig {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bigwig} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}

			sample_Name=$(basename {input.filter_tn5_filter_chromosome_noduplicate_Bed})
			sample_Name=${{sample_Name%.bed.gz}}
			
			macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_noduplicate_Bed} --name ${{sample_Name}}.macs2_broad {config_peak_calling_Dict[MACS2_BROADPEAK]} \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			raw_peak_count=$(wc -l {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak | awk '{{printf "%f", $1/1000000}}')
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bedtools intersect -v -a {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp 2>> {log.macs2_broadPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed} 2>> {log.macs2_broadPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}

			macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_treat_pileup.bdg \
			-c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad \
			--outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ -m ppois -S $raw_peak_count 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp 2>> {log.macs2_broadPeak_stderr}
			bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph} 2>> {log.macs2_broadPeak_stderr}
			tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			bedGraphToBigWig {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bigwig} 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp 1>> {log.macs2_broadPeak_stdout} 2>> {log.macs2_broadPeak_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#Report
			printf "%s\\t" "input:" > {log.macs2_broadPeak_report}
			printf "%s\\t" "{input.filter_tn5_filter_chromosome_Bed}" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "{input.filter_tn5_filter_chromosome_noduplicate_Bed}" >> {log.macs2_broadPeak_report}
			
			printf "%s\\t" "output:" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}">> {log.macs2_broadPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_macs2_broadPeak_Bigwig}">> {log.macs2_broadPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "{output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bigwig}" >> {log.macs2_broadPeak_report}
			
			printf "%s\\t" "commandline:" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_Bed} --name ${{sample_Name}}.macs2_broad {config_peak_calling_Dict[MACS2_BROADPEAK]} --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bedtools intersect -v -a {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ -m ppois -S $raw_peak_count" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bedGraphToBigWig {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bigwig}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_macs2_broadPeak_Bedgraph}.sorted.tmp" >> {log.macs2_broadPeak_report}

			printf "%s\\t" "macs2 callpeak --treatment {input.filter_tn5_filter_chromosome_noduplicate_Bed} --name ${{sample_Name}}.macs2_broad {config_peak_calling_Dict[MACS2_BROADPEAK]} --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bedtools intersect -v -a {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bed}.sorted.filtered.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "macs2 bdgcmp -t {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/ -m ppois -S $raw_peak_count" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "slopBed -i {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling/macs2_broadPeak/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "LC_COLLATE=C sort -k1,1 -k2,2n {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp  > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bgzip -c {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp > {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "tabix -f -p bed {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}" >> {log.macs2_broadPeak_report}
			printf "%s\\t" "bedGraphToBigWig {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp {config_reference_Dict[CHROM_SIZE]} {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bigwig}" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "rm -f {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.tmp {output.filter_tn5_filter_chromosome_noduplicate_macs2_broadPeak_Bedgraph}.sorted.tmp" >> {log.macs2_broadPeak_report}


			printf "%s\\t" "stdout:" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "{log.macs2_broadPeak_stdout}" >> {log.macs2_broadPeak_report}
			
			printf "%s\\t" "stderr:" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "{log.macs2_broadPeak_stderr}" >> {log.macs2_broadPeak_report}

			printf "%s\\t" "elapsed_time:" >> {log.macs2_broadPeak_report}
			printf "%s\\n" "$elapsed_time" >> {log.macs2_broadPeak_report}
			#
			#
		""")


#+++++++++++++++++++++++++++++
##REALM5: PEAK-PROCESS
#+++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++
##REALM6: GENOME-BROWSER
#+++++++++++++++++++++++++++++

rule IGV_Trackhub:
	input:
		alignment_List + post_alignment_List + peak_calling_List
	output:
		igv_trackhub = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/igv_trackhub.xml",
	priority: 988
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
			##BED
			IGV_String += '''
					<Category name="Bed">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bed.gz", recursive=True):
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
			##BEDGRAPH
			IGV_String += '''
					<Category name="Bedgraph">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bedgraph.gz", recursive=True):
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


rule UCSC_Trackhub:
	input:
		alignment_List + post_alignment_List + peak_calling_List
	output:
		ucsc_trackhub = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/genome_browser/ucsc_trackhub.txt"
	priority: 988
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		ucsc_String = ''
		for design in design_Dict:
			#
			##BAM
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bam", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				ucsc_String += '''
				track type=bam name="''' + design + '_' + sample_name + '''" description="" bamColorMode=gray maxWindowToDraw=200000 db=hg38 visibility=pack bigDataUrl=''' + accessible_each_track + '''
				'''
			##BIGWIG
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/" + design + "/**/*.bigwig", recursive=True):
				#
				accessible_each_track = config_track_hub_Dict["HPC_DATASHARE"] + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				ucsc_String += '''
				track type=bigWig name="''' + design + '_' + sample_name + '''" description="" db=hg38 bigDataUrl=''' + accessible_each_track + '''
				'''

			
		f = open(output.ucsc_trackhub, "w")
		f.write(ucsc_String)
		f.close()






