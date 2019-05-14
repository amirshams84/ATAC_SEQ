shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Nov-20-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for NGS DATA Quality Control
# source /data/shamsaddinisha/conda/etc/profile.d/conda.sh
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --debug-dag --cores=50
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --rulegraph | dot -Tsvg > atac_seq.svg
# ################################### TO DO LIST ##################################### #
"""
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
# ################################### CLUSTER CONFIGURATION ###################### #


PROCESSORS = 10
MEMORY = 4000
# ################################### DATA CONFIGURATION ######################## #


configfile: "atac_seq.json"

# ++++++++++++++++++++++++++++++++++++
config_general_Dict = config["GENERAL"]
TITLE = config_general_Dict["TITLE"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
GENOME = config_general_Dict["GENOME"]
GENOME = GENOME.lower()
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CONDA_PY2 = config_conda_Dict["ATAC_Seq_py2"]
CONDA_PY3 = config_conda_Dict["ATAC_Seq_py3"]
ACTIVATE_CONDA_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_CONDA_PY3 = config_conda_Dict["ACTIVATE_PY3"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
metadata_Dict = utility.build_metadata_dict(METADATA_FILE)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_pre_process_Dict = config["PRE_PROCESS"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_alignment_Dict = config["ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_post_alignment_Dict = config["POST_ALIGNMENT"]
CHROMOSOME_FILTER_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["CHROMOSOME_FILTER"])
AWK_TN5_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["TN5_PROCESS"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
config_track_hub_Dict = config["TRACK_HUB"]
#data_share_Link_prefix = config_track_hub_Dict["DATA_SHARE_LINK"]

# ------------------------------------

# ################################### WILDCARDS FUNCTIONS ############################# #


def get_forward_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + "*R1*" + "fastq.gz", recursive=True)


def get_reverse_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + "*R2*" + "fastq.gz", recursive=True)


def get_bed_list(wildcards):
	bed_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bed_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{case}.processed.bed.gz".format(design=wildcards.design, case=sample))
	return bed_List


def get_bam_list(wildcards):
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{case}.processed.bam".format(design=wildcards.design, case=sample))
	return bam_List


def get_bam_tracks(wildcards):
	bam_List = []
	for design in design_Dict:
		if design == wildcards.design:
			bam_List.extend(glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/*.bam"))
	return bam_List

# ################################### WILDCARDS ################################## #
localrules: End_Point
design_Dict = {}
rule_pre_process_List = []
rule_bowtie2_List = []
rule_post_alignment_List = []
rule_macs2_List = []
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
	rule_pre_process_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	rule_bowtie2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	rule_post_alignment_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))
else:
	pass

pooled_bam_List = []
overlapped_bam_List = []
pooled_bed_List = []
overlapped_bed_List = []

for design in design_Dict:
	#
	pooled_bam_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	overlapped_bam_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{overlapped_case}.processed.bam".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	pooled_bed_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{pooled_case}.processed.bed.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	overlapped_bed_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{overlapped_case}.processed.bed.gz".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	#
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{overlapped_case}.narrowPeak.gz".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{overlapped_case}.broadPeak.gz".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	#
	for case in design_Dict[design]["Case"]:
		#
		for control in design_Dict[design]["Control"]:
			#
			rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
			rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
			rule_macs2_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{overlapped_case}_VS_{control}.narrowPeak.gz".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))

REPORTDIR = WORKDIR + "/" + TITLE + "/" + GENOME + "/execution_report"
report_Dict = collections.OrderedDict()
report_Dict = config["EXECUTION_REPORT"]
rule_execution_report_List = [REPORTDIR + "/provenance_overview_diagram.svg"]

IGV_SERVER_DIR = WORKDIR + "/" + TITLE + "/" + GENOME + "/IGV_Server"
IGV_server_file_List = [IGV_SERVER_DIR + "/" + EXPERIMENT + ".xml"]
UCSC_SERVER_DIR = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/UCSC_Server"
IGV_server_file_List = [UCSC_SERVER_DIR + "/" + PROJECT + "_" + EXPERIMENT + "_" + "UCSC_SERVER.txt"]

# ################################### RULES ###################################### #


rule End_Point:
	input: rule_pre_process_List + rule_bowtie2_List + rule_post_alignment_List + pooled_bam_List + pooled_bed_List + overlapped_bam_List + overlapped_bed_List + rule_macs2_List + rule_execution_report_List + IGV_server_file_List

rule Execution_Report:
	input: rule_pre_process_List + rule_bowtie2_List + rule_post_alignment_List + pooled_bam_List + pooled_bed_List + overlapped_bam_List + overlapped_bed_List + rule_macs2_List
	output:
		REPORTDIR + "/provenance_overview_diagram.svg",
	run:
		shell("""
			module load graphviz
			snakemake --snakefile atac_seq.py --configfile atac_seq.json --rulegraph | dot -Tsvg > {output[0]}
			
			unzip /data/RTB/datashare/Amir/DROPBOX/assets.zip -d {REPORTDIR} > {REPORTDIR}/stdout 2>{REPORTDIR}/stderr
			rm -rf {REPORTDIR}/__MACOSX
			chmod -R 2770 {REPORTDIR}/
		
			
			""")
		for each_title in report_Dict:
			#
			for each_subtitle in report_Dict[each_title]:
				#
				each_subtitle_String = ""
				report_List = glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/*/report/" + "/*." + each_title + "." + each_subtitle + ".report")
				for each_report in report_List:
					#
					with open(each_report, 'r') as report_file:
						#
						report_file_String = report_file.read()
						each_subtitle_String += report_file_String
				else:
					pass
				sidebar_String = execution_report.build_sidebar_html_string(report_Dict, "./", REPORTDIR, "ATAC_Seq", each_title, each_subtitle, "execution_log")
				if each_subtitle == "provenance":
					#
					body_String = execution_report.build_body_html_string("./provenance_overview_diagram.svg", each_title, each_subtitle, "provenance")
				else:
					#
					body_String = execution_report.build_body_html_string(each_subtitle_String, each_title, each_subtitle, "execution_log")
				javascript_String = ""
				main_html_String = execution_report.build_main_html_string("./", REPORTDIR, "ATAC_Seq", sidebar_String, body_String, javascript_String)
				utility.write_string_down(main_html_String, REPORTDIR + "/atac_seq_" + each_title + "_" + each_subtitle + "_" + "execution_log.html")
			else:
				pass
		else:
			pass


rule Pre_Process_PE:
	input:
		forward_fastq_List = get_forward_fastq,
		reverse_fastq_List = get_reverse_fastq,
	output:
		processed_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.discarded.fastq.gz",
	log:
		cutadapt_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.preprocess.cutadapt.report",
		fastqc_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.preprocess.fastqc.report",
	priority: 999
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "pre_process: {wildcards.design}|{wildcards.sample}"
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
				processed_forward_fastq={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{processed_forward_fastq}
				processed_reverse_fastq={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{processed_reverse_fastq}
				discarded_forward_fastq={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{discarded_forward_fastq}
				discarded_reverse_fastq={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{discarded_reverse_fastq}
				cutadapt_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{cutadapt_stdout}
				cutadapt_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/cutadapt/{cutadapt_stderr}

				##
				##
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
				#
				printf "<tr>\\n" >> {log.cutadapt_report}
				printf "<td>%s</td>" "cutadapt" >> {log.cutadapt_report}
				printf "<td>%s</td>" "{wildcards.design}" >> {log.cutadapt_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.cutadapt_report}
				
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "INPUT1:" "{each_forward_fastq}" "INPUT2:" "{each_reverse_fastq}" >> {log.cutadapt_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p>" "OUTPUT1:" "$processed_forward_fastq" "OUTPUT2:" "$processed_reverse_fastq" "OUTPUT3:" "$discarded_forward_fastq" "OUTPUT4:" "$discarded_reverse_fastq" >> {log.cutadapt_report}
				printf "<p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT5:" "{output.processed_R1}" "OUTPUT6:" "{output.processed_R2}" "OUTPUT7:" "{output.discarded_R1}" "OUTPUT8:" "{output.discarded_R2}" >> {log.cutadapt_report}
				printf "<td><p><code>%s</code></p>" "cutadapt {config_pre_process_Dict[CUTADAPT_PE]} --too-short-output=$discarded_forward_fastq --too-short-paired-output=$discarded_reverse_fastq --output=$processed_forward_fastq --paired-output=$processed_reverse_fastq {each_forward_fastq} {each_reverse_fastq}" >> {log.cutadapt_report}
				printf "<p><code>%s</code></p>" "cat $processed_forward_fastq >> {output.processed_R1}" >> {log.cutadapt_report}
				printf "<p><code>%s</code></p>" "cat $processed_reverse_fastq >> {output.processed_R2}" >> {log.cutadapt_report}
				printf "<p><code>%s</code></p>" "cat $discarded_forward_fastq >> {output.discarded_R1}" >> {log.cutadapt_report}
				printf "<p><code>%s</code></p></td>" "cat $discarded_reverse_fastq >> {output.discarded_R2}" >> {log.cutadapt_report}
				printf "<td><p><small>%s</small></p></td>" "$cutadapt_stdout" >> {log.cutadapt_report}
				printf "<td><p><small>%s</small></p></td>" "$cutadapt_stderr" >> {log.cutadapt_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.cutadapt_report}
				printf "</tr>\\n" >> {log.cutadapt_report}
				#
				#

			""")
		else:
			pass

		shell("""
			{ACTIVATE_CONDA_PY2}
			mkdir -p {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/
			processed_fq1=$(basename {output.processed_R1})
		
			processed_fq1_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq1/.fastq*/_fastqc.html}}
			processed_fq1_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq1/.fastq*/_fastqc.zip}}
			processed_fq1_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq1/R1*/R1.processed.fastqc.stdout}}
			processed_fq1_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq1/R1*/R1.processed.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.processed_R1} 1> $processed_fq1_fastqc_stdout 2> $processed_fq1_fastqc_stderr
			
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "fastqc" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "{output.processed_R1}" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$processed_fq1_fastqc_html" "OUTPUT2:" "$processed_fq1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.processed_R1}" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$processed_fq1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$processed_fq1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			processed_fq2=$(basename {output.processed_R2})

			processed_fq2_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq2/.fastq*/_fastqc.html}}
			processed_fq2_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq2/.fastq*/_fastqc.zip}}
			processed_fq2_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq2/R2*/R2.processed.fastqc.stdout}}
			processed_fq2_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{processed_fq2/R2*/R2.processed.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.processed_R2} 1> $processed_fq2_fastqc_stdout 2> $processed_fq2_fastqc_stderr
			
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "fastqc" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "{output.processed_R2}" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$processed_fq2_fastqc_html" "OUTPUT2:" "$processed_fq2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.processed_R2}" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$processed_fq2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$processed_fq2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			discarded_fq1=$(basename {output.discarded_R1})
			
			discarded_fq1_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq1/.fastq*/_fastqc.html}}
			discarded_fq1_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq1/.fastq*/_fastqc.zip}}
			discarded_fq1_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq1/R1*/R1.discarded.fastqc.stdout}}
			discarded_fq1_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq1/R1*/R1.discarded.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.discarded_R1} 1> $discarded_fq1_fastqc_stdout 2> $discarded_fq1_fastqc_stderr
			
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "fastqc" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "{output.discarded_R1}" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq1_fastqc_html" "OUTPUT2:" "$discarded_fq1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.discarded_R1}" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			discarded_fq2=$(basename {output.discarded_R2})
			
			discarded_fq2_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq2/.fastq*/_fastqc.html}}
			discarded_fq2_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq2/.fastq*/_fastqc.zip}}
			discarded_fq2_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq2/R2*/R2.discarded.fastqc.stdout}}
			discarded_fq2_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/${{discarded_fq2/R2*/R2.discarded.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.discarded_R2} 1> $discarded_fq2_fastqc_stdout 2> $discarded_fq2_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "fastqc" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "{output.discarded_R2}" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq2_fastqc_html" "OUTPUT2:" "$discarded_fq2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/pre_process/fastqc/ -f fastq --threads {threads} {output.discarded_R2}" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#
		""")


rule Bowtie2_PE:
	input:
		processed_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/cutadapt/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz"
	output:
		alignment_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		alignment_Report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment_report.txt",
	log:
		#BOWTIE2
		bowtie2_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.report",
		bowtie2_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stdout",
		bowtie2_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stderr",
		#FASTQC
		fastqc_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.report",
		fastqc_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.stdout",
		fastqc_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/fastqc/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.stderr",
	priority: 998
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Bowtie2: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			
			{ACTIVATE_CONDA_PY2}

			mapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.mapped.fastq.gz
			mapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.mapped.fastq.gz
			unmapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.unmapped.fastq.gz
			unmapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.unmapped.fastq.gz
			
			##
			##
			start_time="$(date -u +%s)"
			bowtie2 {config_alignment_Dict[BOWTIE2_PE]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_R1} -2 {input.processed_R2} \
			--un-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.unmapped.fastq.gz \
			--al-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.mapped.fastq.gz \
			2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample} -o {output.alignment_Bam} - 1> {log.bowtie2_stdout} 2> {log.bowtie2_stderr}
			samtools index -@ {threads} -b {output.alignment_Bam} >> {log.bowtie2_stdout} 2>> {log.bowtie2_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.bowtie2_report}
			printf "<td>%s</td>" "BOWTIE2" >> {log.bowtie2_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.bowtie2_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.bowtie2_report}
			
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "INPUT1:" "{input.processed_R1}" "INPUT2:" "{input.processed_R2}" >> {log.bowtie2_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.alignment_Bam}" "OUTPUT2:" "{output.alignment_Report}" "OUTPUT3:" "$mapped_R1" "OUTPUT4:" "$mapped_R2" "OUTPUT5:" "$unmapped_R1" "OUTPUT6:" "$unmapped_R2" >> {log.bowtie2_report}
			
			printf "<td><code>%s</code><br>" "bowtie2 {config_alignment_Dict[BOWTIE2_PE]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]}  -1 {input.processed_R1} -2 {input.processed_R2} --un-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.unmapped.fastq.gz --al-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R%.processed.mapped.fastq.gz 2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample} -o {output.alignment_Bam} -" >> {log.bowtie2_report}
			printf "<code>%s</code></td>" "samtools index -@ {threads} -b {output.alignment_Bam}" >> {log.bowtie2_report}
			
			printf "<td><p><small>%s</small></p></td>" "{log.bowtie2_stdout}" >> {log.bowtie2_report}
			printf "<td><p><small>%s</small></p></td>" "{log.bowtie2_stderr}" >> {log.bowtie2_report}

			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.bowtie2_report}
			printf "</tr>\\n" >> {log.bowtie2_report}
			#
			#

			mapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.mapped.fastq.gz
			mapped_R1_Name=$(basename $mapped_R1)
			mapped_R1_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R1_Name/.fastq*/_fastqc.html}}
			mapped_R1_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R1_Name/.fastq*/_fastqc.zip}}
			mapped_R1_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R1_Name/R1*/R1.processed.mapped.fastqc.stdout}}
			mapped_R1_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R1_Name/R1*/R1.processed.mapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $mapped_R1 1> $mapped_R1_fastqc_stdout 2> $mapped_R1_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$mapped_R1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$mapped_R1_fastqc_html" "OUTPUT2:" "$mapped_R1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $mapped_R1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$mapped_R1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$mapped_R1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			mapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.mapped.fastq.gz
			mapped_R2_Name=$(basename $mapped_R2)
			mapped_R2_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R2_Name/.fastq*/_fastqc.html}}
			mapped_R2_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R2_Name/.fastq*/_fastqc.zip}}
			mapped_R2_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R2_Name/R2*/R2.processed.mapped.fastqc.stdout}}
			mapped_R2_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{mapped_R2_Name/R2*/R2.processed.mapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $mapped_R2 1> $mapped_R2_fastqc_stdout 2> $mapped_R2_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$mapped_R2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$mapped_R2_fastqc_html" "OUTPUT2:" "$mapped_R2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $mapped_R2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$mapped_R2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$mapped_R2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#


			unmapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R1.processed.unmapped.fastq.gz
			unmapped_R1_Name=$(basename $unmapped_R1)
			unmapped_R1_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R1_Name/.fastq*/_fastqc.html}}
			unmapped_R1_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R1_Name/.fastq*/_fastqc.zip}}
			unmapped_R1_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R1_Name/R1*/R1.processed.unmapped.fastqc.stdout}}
			unmapped_R1_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R1_Name/R1*/R1.processed.unmappedd.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $unmapped_R1 1> $unmapped_R1_fastqc_stdout 2> $unmapped_R1_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$unmapped_R1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$unmapped_R1_fastqc_html" "OUTPUT2:" "$unmapped_R1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $unmapped_R1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$unmapped_R1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$unmapped_R1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			unmapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bowtie2/{wildcards.sample}_R2.processed.unmapped.fastq.gz
			unmapped_R2_Name=$(basename $unmapped_R2)
			unmapped_R2_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R2_Name/.fastq*/_fastqc.html}}
			unmapped_R2_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R2_Name/.fastq*/_fastqc.zip}}
			unmapped_R2_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R2_Name/R2*/R2.processed.unmapped.fastqc.stdout}}
			unmapped_R2_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/${{unmapped_R2_Name/R2*/R2.processed.unmapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $unmapped_R2 1> $unmapped_R2_fastqc_stdout 2> $unmapped_R2_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$unmapped_R2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$unmapped_R2_fastqc_html" "OUTPUT2:" "$unmapped_R2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/fastqc/ -f fastq --threads {threads} $unmapped_R2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$unmapped_R2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$unmapped_R2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#
		""")


rule Post_Alignment:
	input:
		alignment_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bowtie2/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
	output:
		#PICARD
		marked_duplicate_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.marked_duplicate.bam",
		picard_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard_report.txt",
		#QUALIMAP
		qualimap_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}/qualimapReport.html",
		#SAMTOOLS
		processed_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		duplicate_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.duplicate.bam",
		discarded_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.discarded.bam",
		#BEDTOOLS
		Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bed.gz",
		filtered_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.filtered.bed.gz",
		processed_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bed.gz",

	log:
		#PICARD
		picard_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.report",
		picard_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.stdout",
		picard_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/picard/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.stderr",
		#QUALIMAP
		qualimap_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.report",
		qualimap_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.stdout",
		qualimap_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/qualimap/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.stderr",
		#SAMTOOLS
		samtools_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.report",
		samtools_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.stdout",
		samtools_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.stderr",
		#TAGALIGN
		bedtools_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.bedtools.report",
		bedtools_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.bedtools.stdout",
		bedtools_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.bedtools.stderr",

	priority: 997
	threads: PROCESSORS
	message: "Post_Alignment: {wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			Sample_Name=$(basename {input.alignment_Bam})
			Sample_Name=${{Sample_Name%.bam}}

			##
			##
			start_time="$(date -u +%s)"
			export _JAVA_OPTIONS="-Xms{MEMORY}M -Xmx{MEMORY}M -XX:ParallelGCThreads={threads}"
			picard MarkDuplicates INPUT={input.alignment_Bam} OUTPUT={output.marked_duplicate_Bam} METRICS_FILE={output.picard_report} {config_post_alignment_Dict[PICARD]} 1> {log.picard_stdout} 2> {log.picard_stderr}
			samtools index -@ {threads} -b {output.marked_duplicate_Bam}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.picard_report}
			printf "<td>%s</td>" "PICARD" >> {log.picard_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.picard_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.picard_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.alignment_Bam}" >> {log.picard_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.marked_duplicate_Bam}" "OUTPUT2:" "{output.picard_report}" >> {log.picard_report}
			printf "<td><p><code>%s</code></p><p><code>%s</code></p>" "export _JAVA_OPTIONS='-Xms{MEMORY}M -Xmx{MEMORY}M -XX:ParallelGCThreads={threads}'" "picard MarkDuplicates INPUT={input.alignment_Bam} OUTPUT={output.marked_duplicate_Bam} METRICS_FILE={output.picard_report} {config_post_alignment_Dict[PICARD]}" >> {log.picard_report}
			printf "<p><code>%s</code></p></td>" "samtools index -@ {threads} -b {output.marked_duplicate_Bam}" >> {log.picard_report}
			printf "<td><p><small>%s</small></p></td>" "{log.picard_stdout}" >> {log.picard_report}
			printf "<td><p><small>%s</small></p></td>" "{log.picard_stderr}" >> {log.picard_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.picard_report}
			#
			#

			##
			##
			start_time="$(date -u +%s)"
			qualimap bamqc -bam {output.marked_duplicate_Bam} -nt {threads} -outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/qualimap/{wildcards.sample} {config_post_alignment_Dict[QUALIMAP]} 1> {log.qualimap_stdout} 2> {log.qualimap_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.qualimap_report}
			printf "<td>%s</td>" "QUALIMAP" >> {log.qualimap_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.qualimap_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.qualimap_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{output.marked_duplicate_Bam}" >> {log.qualimap_report}
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.qualimap_report}" >> {log.qualimap_report}
			printf "<td><p><code>%s</code></p></td>" "qualimap bamqc -bam {output.marked_duplicate_Bam} -nt {threads} -outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/qualimap/{wildcards.sample}  {config_post_alignment_Dict[QUALIMAP]}" >> {log.qualimap_report}

			printf "<td><p><small>%s</small></p></td>" "{log.picard_stdout}" >> {log.qualimap_report}
			printf "<td><p><small>%s</small></p></td>" "{log.picard_stderr}" >> {log.qualimap_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.qualimap_report}
			#
			#


			##
			##
			start_time="$(date -u +%s)"
			samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_FILTER]} {output.marked_duplicate_Bam} -b -O BAM -o {output.processed_Bam} 1> {log.samtools_stdout} 2> {log.samtools_stderr}
			samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_DUPLICATE]} {output.marked_duplicate_Bam} -b -O BAM -o {output.duplicate_Bam} 1>> {log.samtools_stdout} 2>> {log.samtools_stderr}
			samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_DISCARDED]} {output.marked_duplicate_Bam} -b -O BAM -o {output.discarded_Bam} 1>> {log.samtools_stdout} 2>> {log.samtools_stderr}
			samtools index -@ {threads} -b {output.processed_Bam}
			samtools index -@ {threads} -b {output.duplicate_Bam}
			samtools index -@ {threads} -b {output.discarded_Bam}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##


			#
			#
			printf "<tr>\\n" > {log.samtools_report}
			printf "<td>%s</td>" "SAMTOOLS" >> {log.samtools_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.samtools_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.samtools_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{output.marked_duplicate_Bam}" >> {log.samtools_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.processed_Bam}" "OUTPUT2:" "{output.duplicate_Bam}" "OUTPUT3:" "{output.discarded_Bam}" >> {log.samtools_report}
			printf "<td><p><code>%s</code></p><p><code>%s</code></p>" "samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_FILTER]} {output.marked_duplicate_Bam} -b -O BAM -o {output.processed_Bam}" "samtools index -@ {threads} -b {output.processed_Bam}" >> {log.samtools_report}
			printf "<p><code>%s</code></p><p><code>%s</code></p>" "samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_DUPLICATE]} {output.marked_duplicate_Bam} -b -O BAM -o {output.duplicate_Bam}" "samtools index -@ {threads} -b {output.duplicate_Bam}" >> {log.samtools_report}
			printf "<p><code>%s</code></p><p><code>%s</code></p></td>" "samtools view --threads {threads} -h {config_post_alignment_Dict[SAMTOOLS_DISCARDED]} {output.marked_duplicate_Bam} -b -O BAM -o {output.discarded_Bam}" "samtools index -@ {threads} -b {output.discarded_Bam}" >> {log.samtools_report}

			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stdout}" >> {log.samtools_report}
			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stderr}" >> {log.samtools_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.samtools_report}
			#
			#

			##
			##
			start_time="$(date -u +%s)"

			samtools sort -n {output.processed_Bam} -O BAM -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.bam 1> {log.bedtools_stdout} 2> {log.bedtools_stderr}
			samtools fixmate --threads {threads} -r {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.bam {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.fixed.bam 1>> {log.bedtools_stdout} 2>> {log.bedtools_stderr}
			bedtools bamtobed -cigar -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.fixed.bam | gzip -nc > {output.Bed} 2>> {log.bedtools_stderr}
			bedtools intersect -v -a <(zcat -f {output.Bed} ) -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) | gzip -nc > {output.filtered_Bed} 2>> {log.bedtools_stderr}
			zcat {output.filtered_Bed} | awk -F'\\t' '{CHROMOSOME_FILTER_PROCESS}' | awk -F'\\t' '{AWK_TN5_PROCESS}' | gzip -nc > {output.processed_Bed} 2>> {log.bedtools_stderr}

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##
			CHROMOSOME_FILTER_COMMAND='{config_post_alignment_Dict[CHROMOSOME_FILTER]}'
			TN5_AWK_COMMAND='{config_post_alignment_Dict[TN5_PROCESS]}'
			#
			#
			printf "<tr>\\n" > {log.bedtools_report}
			printf "<td>%s</td>" "BEDTOOLS" >> {log.bedtools_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.bedtools_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.bedtools_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{output.processed_Bam}" >> {log.bedtools_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.Bed}" "OUTPUT2:" "{output.filtered_Bed}" "OUTPUT3:" "{output.processed_Bed}" >> {log.bedtools_report}

			printf "<td><p><code>%s</code></p>" "samtools sort -n {output.processed_Bam} -O BAM -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.bam" >> {log.bedtools_report}
			printf "<p><code>%s</code></p>" "samtools fixmate --threads {threads} -r {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.bam {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.fixed.bam" >> {log.bedtools_report}
			printf "<p><code>%s</code></p>" "bedtools bamtobed -cigar -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/alignment/bedtools/{wildcards.sample}.nmsort.fixed.bam | gzip -nc > {output.Bed}" >> {log.bedtools_report}
			printf "<p><code>%s</code></p>" "bedtools intersect -v -a <(zcat -f {output.Bed} ) -b <(zcat -f {config_reference_Dict[BLACK_LIST]} ) | gzip -nc > {output.filtered_Bed}" >> {log.bedtools_report}
			printf "<p><code>%s</code></p></td>" "zcat {output.filtered_Bed} | awk -F'\\t' '$CHROMOSOME_FILTER_COMMAND' | awk -F'\\t' '$TN5_AWK_COMMAND' | gzip -nc > {output.processed_Bed}" >> {log.bedtools_report}

			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stdout}" >> {log.bedtools_report}
			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stderr}" >> {log.bedtools_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.bedtools_report}
			#
			#
		""")


rule Pooling:
	input:
		bed_List = get_bed_list,
		bam_List = get_bam_list,
	output:
		pooled_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{pooled_case, .*_POOLED_.*}.processed.bed.gz",
		pooled_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{pooled_case, .*_POOLED_.*}.processed.bam",
	log:
		pooling_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.report",
		pooling_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stdout",
		pooling_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stderr",
	priority: 996
	message: "Pooling: {wildcards.pooled_case}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##
			start_time="$(date -u +%s)"
			cat {input.bed_List} > {output.pooled_Bed}
			samtools merge --threads {threads} {output.pooled_Bam} {input.bam_List} 1> {log.pooling_stdout} 2> {log.pooling_stderr}
			samtools index -@ {threads} {output.pooled_Bam} 1> {log.pooling_stdout} 2> {log.pooling_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.pooling_report}
			printf "<td>%s</td>" "POOLING" >> {log.pooling_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.pooling_report}
			printf "<td>%s</td>" "{wildcards.pooled_case}" >> {log.pooling_report}

			printf "<td>" >> {log.pooling_report}
			declare -a bam_List=({input.bam_List})

			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{bam_List[$case_index]}}" >> {log.pooling_report}
			done

			declare -a bed_List=({input.bed_List})

			for case_index in "${{!bed_List[@]}}"
			do
				index=$(($index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{bed_List[$case_index]}}" >> {log.pooling_report}
			done

			printf "</td>" >> {log.pooling_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.pooled_Bam}" "OUTPUT2:" "{output.pooled_Bed}" >> {log.pooling_report}

			printf "<td><p><code>%s</code></p><p><code>%s</code></p><p><code>%s</code></p></td>" "cat {input.bed_List} > {output.pooled_Bed}" "samtools merge --threads {threads} {output.pooled_Bam} {input.bam_List}" "samtools index -@ {threads} {output.pooled_Bam}" >> {log.pooling_report}

			printf "<td><p><small>%s</small></p></td>" "{log.pooling_stdout}" >> {log.pooling_report}
			printf "<td><p><small>%s</small></p></td>" "{log.pooling_stderr}" >> {log.pooling_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.pooling_report}
			#
			#

		""")


rule Overlapping:
	input:
		bed_List = get_bed_list,
		bam_List = get_bam_list,
	output:
		overlapped_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{overlapped_case, .*_OVERLAPPED_.*}.processed.bed.gz",
		overlapped_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/samtools/{overlapped_case, .*_OVERLAPPED_.*}.processed.bam",
	log:
		overlapping_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.report",
		overlapping_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.stdout",
		overlapping_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.stderr",
	priority: 996
	message: "Overlapping: {wildcards.overlapped_case}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			##
			##
			start_time="$(date -u +%s)"
			declare -a bam_List=({input.bam_List})
			bam_group_A=""
			bam_group_B=""
			for case_index in "${{!bam_List[@]}}"
			do
				if [ "$case_index" -eq 0 ] ; then
					bam_group_A+="${{bam_List[$case_index]}} "
				else
					bam_group_B+="${{bam_List[$case_index]}} "
				fi
			done
			bedtools intersect -wa -a $bam_group_A -b $bam_group_B {config_post_alignment_Dict[BEDTOOLS_OVERLAPPING]} > {output.overlapped_Bam}  2> {log.overlapping_stderr}
			samtools index -@ {threads} {output.overlapped_Bam} 1>> {log.overlapping_stdout} 2>> {log.overlapping_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			##
			##
			start_time="$(date -u +%s)"
			declare -a bed_List=({input.bed_List})
			bed_group_A=""
			bed_group_B=""
			for case_index in "${{!bed_List[@]}}"
			do
				if [ "$case_index" -eq 0 ] ; then
					bed_group_A+="${{bed_List[$case_index]}} "
				else
					bed_group_B+="${{bed_List[$case_index]}} "
				fi
			done
			bedtools intersect -wa -a $bed_group_A -b $bed_group_B {config_post_alignment_Dict[BEDTOOLS_OVERLAPPING]} | gzip -nc > {output.overlapped_Bed}  2> {log.overlapping_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##



			#
			#
			printf "<tr>\\n" > {log.overlapping_report}
			printf "<td>%s</td>" "OVERLAPPING" >> {log.overlapping_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.overlapping_report}
			printf "<td>%s</td>" "{wildcards.overlapped_case}" >> {log.overlapping_report}

			printf "<td>" >> {log.overlapping_report}
			declare -a bam_List=({input.bam_List})

			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{bam_List[$case_index]}}" >> {log.overlapping_report}
			done

			declare -a bed_List=({input.bed_List})

			for case_index in "${{!bed_List[@]}}"
			do
				index=$(($index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{bed_List[$case_index]}}" >> {log.overlapping_report}
			done

			printf "</td>" >> {log.overlapping_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.overlapped_Bam}" "OUTPUT2:" "{output.overlapped_Bed}" >> {log.overlapping_report}

			printf "<td><p><code>%s</code></p><p><code>%s</code></p><p><code>%s</code></p></td>" "bedtools intersect -wa -a $bed_group_A -b $bed_group_B {config_post_alignment_Dict[BEDTOOLS_OVERLAPPING]} | gzip -nc > {output.overlapped_Bed}" "bedtools intersect -wa -a $bam_group_A -b $bam_group_B {config_post_alignment_Dict[BEDTOOLS_OVERLAPPING]} > {output.overlapped_Bam}" "samtools index -@ {threads} {output.overlapped_Bam}" >> {log.overlapping_report}

			printf "<td><p><small>%s</small></p></td>" "{log.overlapping_stdout}" >> {log.overlapping_report}
			printf "<td><p><small>%s</small></p></td>" "{log.overlapping_stderr}" >> {log.overlapping_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.overlapping_report}
			#
			#

		""")


rule Macs2:
	input:
		processed_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{sample, ((?!.*_VS_.*|\\.).)*}.processed.bed.gz",
	output:
		#Macs2_Broad
		narrowPeak_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.gz",
		narrowPeak_Bigwig = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.bigwig",
		#Macs2_Broad
		broadPeak_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.gz",
		broadPeak_Bigwig = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.bigwig",
	log:
		#Macs2_narrow
		macs2_narrow_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.report",
		macs2_narrow_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.stdout",
		macs2_narrow_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.stderr",
		#Macs2_broad
		macs2_broad_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.report",
		macs2_broad_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.stdout",
		macs2_broad_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{sample, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.stderr",
	priority: 996
	message: "Macs2: {wildcards.design}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {input.processed_Bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			##
			##
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.processed_Bed} --name ${{sample_Name}}_narrow {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ 1> {log.macs2_narrow_stdout} 2> {log.macs2_narrow_stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz 2>> {log.macs2_narrow_stderr}
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}_narrow \
			--outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph  2>> {log.macs2_narrow_stderr}
			mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##


			#
			#
			printf "<tr>\\n" > {log.macs2_narrow_report}
			printf "<td>%s</td>" "NARROWPEAK" >> {log.macs2_narrow_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.macs2_narrow_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.macs2_narrow_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.processed_Bed}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s %s</small></p>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT2:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT3:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT4:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT5:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT6:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p></td>" "OUTPUT7:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig" >> {log.macs2_narrow_report}


			printf "<td><p><code>%s</code></p>" "macs2 callpeak --treatment {input.processed_Bed} --name ${{sample_Name}}_narrow {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}_narrow --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p></td>" "bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig" >> {log.macs2_narrow_report}

			printf "<td><p><small>%s</small></p></td>" "{log.macs2_narrow_stdout}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s</small></p></td>" "{log.macs2_narrow_stderr}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.macs2_narrow_report}
			#
			#

			sample_Name=$(basename {input.processed_Bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			##
			##
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.processed_Bed} --name ${{sample_Name}}_broad {config_peak_calling_Dict[MACS2_BROAD]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ 1> {log.macs2_broad_stdout} 2> {log.macs2_broad_stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz 2>> {log.macs2_broad_stderr}
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg --o-prefix ${{sample_Name}}_broad \
			--outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph  2>> {log.macs2_broad_stdout}
			mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.macs2_broad_report}
			printf "<td>%s</td>" "BROADPEAK" >> {log.macs2_broad_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.macs2_broad_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.macs2_broad_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.processed_Bed}" >> {log.macs2_broad_report}
			printf "<td><p><small>%s %s</small></p>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT2:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT3:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT4:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT5:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT6:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p></td>" "OUTPUT7:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig" >> {log.macs2_broad_report}


			printf "<td><p><code>%s</code></p>" "macs2 callpeak --treatment {input.processed_Bed} --name ${{sample_Name}}_broad {config_peak_calling_Dict[MACS2_BROAD]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg --o-prefix ${{sample_Name}}_broad --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p></td>" "bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig" >> {log.macs2_broad_report}

			printf "<td><p><small>%s</small></p></td>" "{log.macs2_broad_stdout}" >> {log.macs2_broad_report}
			printf "<td><p><small>%s</small></p></td>" "{log.macs2_broad_stderr}" >> {log.macs2_broad_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.macs2_broad_report}
			#
			#

		""")


rule Macs2_Controlled:
	input:
		case_processed_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{case, ((?!.*_VS_.*|\\.).)*}.processed.bed.gz",
		control_processed_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/alignment/bedtools/{control, ((?!.*_VS_.*|\\.).)*}.processed.bed.gz",
	output:
		#Macs2_Broad
		narrowPeak_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.narrowPeak.gz",
		narrowPeak_Bigwig = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.narrowPeak.bigwig",
		#Macs2_Broad
		broadPeak_Bed = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.broadPeak.gz",
		broadPeak_Bigwig = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.broadPeak.bigwig",
	log:
		#Macs2_narrow
		macs2_narrow_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.report",
		macs2_narrow_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.stdout",
		macs2_narrow_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.narrow_peak.stderr",
		#Macs2_broad
		macs2_broad_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/report/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.report",
		macs2_broad_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.stdout",
		macs2_broad_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/macs2/{case, ((?!.*_VS_.*|\\.).)*}_VS_{control, ((?!.*_VS_.*|\\.).)*}.peak_calling.broad_peak.stderr",
	priority: 996
	message: "Macs2: {wildcards.design}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {output.narrowPeak_Bed})
			
			sample_Name=${{sample_Name%.processed.bed.gz}}

			##
			##
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.case_processed_Bed} --control {input.control_processed_Bed} --name ${{sample_Name}}_narrow {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ 1> {log.macs2_narrow_stdout} 2> {log.macs2_narrow_stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz 2>> {log.macs2_narrow_stderr}
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}_narrow \
			--outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph  2>> {log.macs2_narrow_stderr}
			mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig 1>> {log.macs2_narrow_stdout} 2>> {log.macs2_narrow_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##


			#
			#
			printf "<tr>\\n" > {log.macs2_narrow_report}
			printf "<td>%s</td>" "NARROWPEAK" >> {log.macs2_narrow_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.macs2_narrow_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.macs2_narrow_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.case_processed_Bed}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s %s</small></p></td>" "INPUT2:" "{input.control_processed_Bed}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s %s</small></p>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT2:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT3:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT4:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT5:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT6:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><small>%s %s</small></p></td>" "OUTPUT7:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig" >> {log.macs2_narrow_report}


			printf "<td><p><code>%s</code></p>" "macs2 callpeak --treatment {input.case_processed_Bed} --control {input.control_processed_Bed} --name ${{sample_Name}}_narrow {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_peaks.narrowPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}_narrow --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p>" "tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph.gz" >> {log.macs2_narrow_report}
			printf "<p><code>%s</code></p></td>" "bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_narrow.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.narrowPeak.bigwig" >> {log.macs2_narrow_report}

			printf "<td><p><small>%s</small></p></td>" "{log.macs2_narrow_stdout}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s</small></p></td>" "{log.macs2_narrow_stderr}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.macs2_narrow_report}
			#
			#

			sample_Name=$(basename {output.narrowPeak_Bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}

			##
			##
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.case_processed_Bed} --control {input.control_processed_Bed} --name ${{sample_Name}}_broad {config_peak_calling_Dict[MACS2_BROAD]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ 1> {log.macs2_broad_stdout} 2> {log.macs2_broad_stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz 2>> {log.macs2_broad_stderr}
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg --o-prefix ${{sample_Name}}_broad \
			--outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} \
			{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph  2>> {log.macs2_broad_stdout}
			mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig 1>> {log.macs2_broad_stdout} 2>> {log.macs2_broad_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.macs2_broad_report}
			printf "<td>%s</td>" "BROADPEAK" >> {log.macs2_broad_report}
			printf "<td>%s</td>" "{wildcards.design}" >> {log.macs2_broad_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.macs2_broad_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.case_processed_Bed}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s %s</small></p></td>" "INPUT2:" "{input.control_processed_Bed}" >> {log.macs2_narrow_report}
			printf "<td><p><small>%s %s</small></p>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT2:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT3:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT4:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT5:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p>" "OUTPUT6:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><small>%s %s</small></p></td>" "OUTPUT7:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig" >> {log.macs2_broad_report}


			printf "<td><p><code>%s</code></p>" "macs2 callpeak --treatment {input.case_processed_Bed} --control {input.control_processed_Bed} --name ${{sample_Name}}_broad {config_peak_calling_Dict[MACS2_BROAD]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_peaks.broadPeak | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_control_lambda.bdg --o-prefix ${{sample_Name}}_broad --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/ -m ppois -S 1.0" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p>" "tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph.gz" >> {log.macs2_broad_report}
			printf "<p><code>%s</code></p></td>" "bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}_broad.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/macs2/${{sample_Name}}.broadPeak.bigwig" >> {log.macs2_broad_report}

			printf "<td><p><small>%s</small></p></td>" "{log.macs2_broad_stdout}" >> {log.macs2_broad_report}
			printf "<td><p><small>%s</small></p></td>" "{log.macs2_broad_stderr}" >> {log.macs2_broad_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.macs2_broad_report}
			#
			#

		""")


rule build_IGV_server:
	input:
		rule_post_alignment_List + pooled_bam_List + pooled_bed_List + overlapped_bam_List + overlapped_bed_List + rule_macs2_List
	output:
		IGV_server_xml = IGV_SERVER_DIR + "/" + EXPERIMENT + ".xml"
	run:
		IGV_String = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
		<Global name="''' + PROJECT + '''"  infolink="https://myrtb.nih.gov/Microarray/Lists/Projects" version="1">
			<Category name="''' + EXPERIMENT + '''">
		'''
		for design in design_Dict:
			#
			# alignment bam
			IGV_String += '''
				<Category name="''' + design + '''">
			'''
			####BAM
			IGV_String += '''
					<Category name="Bam">
			'''
			#raw bam
			IGV_String += '''
						<Category name="Raw Bam">
			'''
			for each_bam in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/picard/*.marked_duplicate.bam", recursive=True):
				#
				accessible_each_bam = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bam.split("RTB")[1]
				sample_name = each_bam.split("/alignment/picard/")[1].split(".marked_duplicate.bam")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bam + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			
			#processed bam
			IGV_String += '''
						<Category name="Processed Bam">
			'''
			for each_bam in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/samtools/*.processed.bam", recursive=True):
				#
				accessible_each_bam = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bam.split("RTB")[1]
				sample_name = each_bam.split("/alignment/samtools/")[1].split(".processed.bam")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bam + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#duplicate bam
			IGV_String += '''
						<Category name="Duplicate Bam">
			'''
			for each_bam in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/samtools/*.duplicate.bam", recursive=True):
				#
				accessible_each_bam = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bam.split("RTB")[1]
				sample_name = each_bam.split("/alignment/samtools/")[1].split(".duplicate.bam")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bam + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#discarded bam
			IGV_String += '''
						<Category name="Discarded Bam">
			'''
			for each_bam in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/samtools/*.discarded.bam", recursive=True):
				#
				accessible_each_bam = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bam.split("RTB")[1]
				sample_name = each_bam.split("/alignment/samtools/")[1].split(".discarded.bam")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bam + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			IGV_String += '''	</Category>''' + '\n'
			############
			####BEDTOOLS
			IGV_String += '''
					<Category name="Bed">
			'''
			#raw Bed
			IGV_String += '''
						<Category name="Raw Bed">
			'''
			for each_bed in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/bedtools/*[!processed].bed.gz", recursive=True):
				#
				accessible_each_bed = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bed.split("RTB")[1]
				sample_name = each_bed.split("/alignment/bedtools/")[1].split(".bed.gz")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bed + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#Processed Bed
			IGV_String += '''
						<Category name="Processed Bed">
			'''
			for each_bed in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/bedtools/*.processed.bed.gz", recursive=True):
				#
				accessible_each_bed = config_track_hub_Dict["DATA_SHARE_LINK"] + each_bed.split("RTB")[1]
				sample_name = each_bed.split("/alignment/bedtools/")[1].split(".processed.bed.gz")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_bed + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			IGV_String += '''	</Category>''' + '\n'
			###########
			####PEAK
			IGV_String += '''
					<Category name="PEAK">
			'''
			#narrow peak bed
			IGV_String += '''
						<Category name="Narrow Peak Bed">
			'''
			for each_peak in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/macs2/*.narrowPeak.gz", recursive=True):
				#
				accessible_each_peak = config_track_hub_Dict["DATA_SHARE_LINK"] + each_peak.split("RTB")[1]
				sample_name = each_peak.split("/peak_calling/macs2/")[1].split(".narrowPeak.gz")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_peak + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#
			#narrow peak signal
			IGV_String += '''
						<Category name="Narrow Peak Signal">
			'''
			for each_peak in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/macs2/*.narrowPeak.bigwig", recursive=True):
				#
				accessible_each_peak = config_track_hub_Dict["DATA_SHARE_LINK"] + each_peak.split("RTB")[1]
				sample_name = each_peak.split("/peak_calling/macs2/")[1].split(".narrowPeak.bigwig")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_peak + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#
			#broad peak bed
			IGV_String += '''
						<Category name="Broad Peak Bed">
			'''
			for each_peak in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/macs2/*.broadPeak.gz", recursive=True):
				#
				accessible_each_peak = config_track_hub_Dict["DATA_SHARE_LINK"] + each_peak.split("RTB")[1]
				sample_name = each_peak.split("/peak_calling/macs2/")[1].split(".broadPeak.gz")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_peak + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			#
			#broad peak signal
			IGV_String += '''
						<Category name="Broad Peak Signal">
			'''
			for each_peak in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/macs2/*.broadPeak.bigwig", recursive=True):
				#
				accessible_each_peak = config_track_hub_Dict["DATA_SHARE_LINK"] + each_peak.split("RTB")[1]
				sample_name = each_peak.split("/peak_calling/macs2/")[1].split(".broadPeak.bigwig")[0]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_peak + '''"/>
				'''
			IGV_String += '''		</Category>''' + '\n'
			IGV_String += '''	</Category>''' + '\n'
			IGV_String += '''	</Category>''' + '\n'

		IGV_String += '''	</Category>''' + '\n'
		IGV_String += '''</Global>'''
		f = open(output.IGV_server_xml, "w")
		f.write(IGV_String)
		f.close()


run build_UCSC_server:
	input:
		rule_post_alignment_List + pooled_bam_List + pooled_bed_List + overlapped_bam_List + overlapped_bed_List + rule_macs2_List
	output:
		UCSC_server_text = UCSC_SERVER_DIR + "/" + PROJECT + "_" + EXPERIMENT + "_" + "UCSC_SERVER.txt"
	run:
		#Build






