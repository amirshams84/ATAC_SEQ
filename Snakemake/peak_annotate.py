# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Annotate
# snakemake --snakefile peak_annotate.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile peak_annotate.py --configfile Encode.json --rulegraph | dot -Tsvg > Peak_Calling.svg
# ################################### IMPORT ##################################### #


import os
import sys
import pandas
sys.path.append(os.path.abspath("./library/"))
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
#CONDA
config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CONDA_PY2 = config_conda_Dict["ATAC_Seq_py2"]
CONDA_PY3 = config_conda_Dict["ATAC_Seq_py3"]
ACTIVATE_CONDA_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_CONDA_PY3 = config_conda_Dict["ACTIVATE_PY3"]
# ------------------------------------
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
config_qualimap_Dict = config["QUALIMAP"][GENOME]
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
config_qualimap_Dict = config["QUALIMAP"][GENOME]
config_reference_Dict = config["REFERENCE"][GENOME]
EFFECTIVE_GENOME_SIZE = config_reference_Dict["EFFECTIVE_GENOME_SIZE"]
# ------------------------------------
# ################################### WILDCARDS ################################ #


post_alignment_List = []
peak_calling_List = []

peak_annotate_List = []
peak_overlap_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	#
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{sample}/{sample}.narrowPeak.homer_annotate.txt".format(design=sample_Dict["Design"], sample=sample))
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{sample}/{sample}.broadPeak.homer_annotate.txt".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{pooled_case}/{pooled_case}.narrowPeak.homer_annotate.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{pooled_case}/{pooled_case}.broadPeak.homer_annotate.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	#
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{overlapped}/{overlapped}.narrowPeak.homer_annotate.txt".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{overlapped}/{overlapped}.broadPeak.homer_annotate.txt".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	##
	#
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	#
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{IDR}/{IDR}.narrowPeak.homer_annotate.txt".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{IDR}/{IDR}.broadPeak.homer_annotate.txt".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
			#
			peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{case}_VS_{control}/{case}_VS_{control}.narrowPeak.homer_annotate.txt".format(design=design, case=case, control=control))
			peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{case}_VS_{control}/{case}_VS_{control}.broadPeak.homer_annotate.txt".format(design=design, case=case, control=control))
	##
	for control in design_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{pooled_case}_VS_{control}/{pooled_case}_VS_{control}.narrowPeak.homer_annotate.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{pooled_case}_VS_{control}/{pooled_case}_VS_{control}.broadPeak.homer_annotate.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}_VS_{control}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}_VS_{control}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		#
		#
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{overlapped}_VS_{control}/{overlapped}_VS_{control}.narrowPeak.homer_annotate.txt".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{overlapped}_VS_{control}/{overlapped}_VS_{control}.broadPeak.homer_annotate.txt".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}_VS_{control}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		peak_overlap_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}_VS_{control}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		#
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{IDR}_VS_{control}/{IDR}_VS_{control}.narrowPeak.homer_annotate.txt".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		peak_annotate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{IDR}_VS_{control}/{IDR}_VS_{control}.broadPeak.homer_annotate.txt".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_annotate_List
# ################################### PIPELINE RULES ########################## #


rule Homer_NarrowPeak_Annotate:
	"""
	"""
	input:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bdg",
	output:
		narrowPeak_annotate = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{sample}/{sample}.narrowPeak.homer_annotate.txt",
	priority: 996
	threads: PROCESSORS
	message: "Homer_NarrowPeak_Annotate: {wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		bedgraph_List = []
		if "_VS_" in wildcards.sample:
			sample_name_List = wildcards.sample.split("_VS_")
			bedgraph_List.append(input.narrowPeak_bdg)
			for each_sample in sample_name_List:
				#
				bedgraph_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{each_sample}.narrowPeak.bdg".format(design=design, each_sample=each_sample))
		else:
			bedgraph_List.append(input.narrowPeak_bdg)

		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load homer/4.10.1 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_annotate/narrowpeak/{wildcards.sample}
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_annotate/narrowpeak/
			mkdir -p $QC_PATH
			#
			homer_annotate_file=$OUT_PATH/annotate.txt
			homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load homer/4.10.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Annotate Narrow"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.narrowPeak_bed}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_annotate}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "gunzip < {input.narrowPeak_bed} > $OUT_PATH/{wildcards.sample}.narrowPeak.bed" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_file=$OUT_PATH/annotate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "annotatePeaks.pl $OUT_PATH/{wildcards.sample}.narrowPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} -go $OUT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "for file in $OUT_PATH/*.txt" | tee >(cat >&2)
			printf "%s\\n" "do" | tee >(cat >&2)
			printf "\\t%s\\n" "name=$(basename \\$file)" | tee >(cat >&2)
			printf "\\t%s\\n" "name={wildcards.sample}.narrowPeak.homer_\\$name" | tee >(cat >&2)
			printf "\\t%s\\n" "mv \\$file $OUT_PATH/\\$name" | tee >(cat >&2)
			printf "%s\\n" "done" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			gunzip < {input.narrowPeak_bed} > $OUT_PATH/{wildcards.sample}.narrowPeak.bed
			homer_annotate_file=$OUT_PATH/annotate.txt
			homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt
			annotatePeaks.pl $OUT_PATH/{wildcards.sample}.narrowPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} -go $OUT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file
			for file in $OUT_PATH/*.txt
			do
				name=$(basename $file)
				name={wildcards.sample}.narrowPeak.homer_$name
				mv $file $OUT_PATH/$name
			done
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")
		processed_column_List = []
		target_Path = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/".format(design=design)
		homer_DF = pandas.read_csv(output.narrowPeak_annotate, sep="\t", low_memory=False, index_col=None)
		column_List = homer_DF.columns.values.tolist()
		column_List[0] = column_List[0].split(" (")[0]
		for each_column in column_List:
			#
			if target_Path in each_column:
				string = each_column.split(target_Path)[1]
				string = string.split(".bdg")[0] + "_read_coverage_counts"

				processed_column_List.append(string)
			else:
				processed_column_List.append(each_column)
		#
		homer_DF.columns = processed_column_List
		homer_DF.to_csv(output.narrowPeak_annotate, sep='\t', index=False, header=True)


rule Homer_BroadPeak_Annotate:
	"""
	"""
	input:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bdg",
	output:
		broadPeak_annotate = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{sample}/{sample}.broadPeak.homer_annotate.txt",
	priority: 996
	threads: PROCESSORS
	message: "Homer_BroadPeak_Annotate: {wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		bedgraph_List = []
		if "_VS_" in wildcards.sample:
			sample_name_List = wildcards.sample.split("_VS_")
			bedgraph_List.append(input.broadPeak_bdg)
			for each_sample in sample_name_List:
				#
				bedgraph_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{each_sample}.broadPeak.bdg".format(design=design, each_sample=each_sample))
		else:
			bedgraph_List.append(input.broadPeak_bdg)

		shell("""
			#
			module load macs/2.1.2 || exit 1
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load ucsc/373 || exit 1
			module load homer/4.10.1 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_annotate/broadpeak/{wildcards.sample}
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_annotate/broadpeak/
			mkdir -p $QC_PATH
			#
			homer_annotate_file=$OUT_PATH/annotate.txt
			homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load macs/2.1.2 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load deeptools/3.1.3 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load homer/4.10.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "Description:"  | tee >(cat >&2)
			printf "%s\\n" "Peak Annotate Broad"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.broadPeak_bed}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.broadPeak_bdg}"  | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_annotate}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "gunzip < {input.broadPeak_bed} > $OUT_PATH/{wildcards.sample}.broadPeak.bed" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_file=$OUT_PATH/annotate.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "annotatePeaks.pl $OUT_PATH/{wildcards.sample}.broadPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} -go $OUT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "for file in $OUT_PATH/*.txt" | tee >(cat >&2)
			printf "%s\\n" "do" | tee >(cat >&2)
			printf "\\t%s\\n" "name=$(basename \\$file)" | tee >(cat >&2)
			printf "\\t%s\\n" "name={wildcards.sample}.broadPeak.homer_\\$name" | tee >(cat >&2)
			printf "\\t%s\\n" "mv \\$file $OUT_PATH/\\$name" | tee >(cat >&2)
			printf "%s\\n" "done" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			gunzip < {input.broadPeak_bed} > $OUT_PATH/{wildcards.sample}.broadPeak.bed
			homer_annotate_file=$OUT_PATH/annotate.txt
			homer_annotate_statistics_file=$OUT_PATH/annotate_statistics.txt
			annotatePeaks.pl $OUT_PATH/{wildcards.sample}.broadPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} -go $OUT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file
			for file in $OUT_PATH/*.txt
			do
				name=$(basename $file)
				name={wildcards.sample}.broadPeak.homer_$name
				mv $file $OUT_PATH/$name
			done
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")
		processed_column_List = []
		target_Path = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/".format(design=design)
		homer_DF = pandas.read_csv(output.broadPeak_annotate, sep="\t", low_memory=False, index_col=None)
		column_List = homer_DF.columns.values.tolist()
		column_List[0] = column_List[0].split(" (")[0]
		for each_column in column_List:
			#
			if target_Path in each_column:
				string = each_column.split(target_Path)[1]
				string = string.split(".bdg")[0] + "_read_coverage_counts"

				processed_column_List.append(string)
			else:
				processed_column_List.append(each_column)
		#
		homer_DF.columns = processed_column_List
		homer_DF.to_csv(output.broadPeak_annotate, sep='\t', index=False, header=True)







