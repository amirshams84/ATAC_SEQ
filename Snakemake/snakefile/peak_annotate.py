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
import glob
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


def get_design(file_path, genome):
	"""
	"""
	each_path = os.path.abspath(file_path)
	each_sample = each_path.split(genome + "/")[1].split("/")[0]
	return each_sample
 ################################### CONFIGURATION ############################## #


# ++++++++++++++++++++++++++++++++++++
#PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../Bash_Script")
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
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
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

rule narrowpeak_annotate:
	"""
	"""
	input:
		narrowpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz",
		narrowpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bdg",
	output:
		narrowpeak_annotate = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/narrowpeak/{sample}/{sample}.narrowPeak.homer_annotate.txt",
	priority: 996
	threads: PROCESSORS
	message: "narrowpeak_annotate: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		bedgraph_List = []
		peak_basename = os.path.basename(input.narrowpeak_bed)
		peak_begining = re.sub(".narrowPeak.gz", "", peak_basename)
		peak_design = get_design(input.narrowpeak_bed, GENOME)
		bedgraph_List.append(input.narrowpeak_bdg)
		target_List = ["_POOLED_", "_OVERLAPPED_", "_IDR_", "_VS_"]
		for each_target in target_List:
			if each_target in peak_begining:
				#
				sample_name_List = peak_begining.split(each_target)
				for each_sample in sample_name_List:
					#
					bedgraph_List.extend(glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + peak_design + "/peak_calling/narrowpeak/*" + each_sample + "*.narrowPeak.bdg"))
				else:
					pass
			else:
				pass
		else:
			pass
		bedgraph_List = list(set(bedgraph_List))
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_annotate/narrowpeak/{wildcards.sample}
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_annotate/narrowpeak/{wildcards.sample}
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "narrowpeak_annotate: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.narrowpeak_bed}" | tee >(cat >&2)
			declare -a bam_List=({bedgraph_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				index=$(($index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowpeak_annotate}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load homer/4.10.1 || exit 1" | tee >(cat >&2)
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
			module load homer/4.10.1 || exit 1

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
			printf "%s\\n" "gunzip < {input.narrowpeak_bed} > $RESULT_PATH/{wildcards.sample}.narrowPeak.bed" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			gunzip < {input.narrowpeak_bed} > $RESULT_PATH/{wildcards.sample}.narrowPeak.bed

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
			printf "%s\\n" "homer_annotate_file=$RESULT_PATH/annotate.txt" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_statistics_file=$RESULT_PATH/annotate_statistics.txt" | tee >(cat >&2)
			printf "%s\\n" "annotatePeaks.pl $RESULT_PATH/{wildcards.sample}.narrowPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} \\
			-go $RESULT_PATH -annStats \$homer_annotate_statistics_file -cpu {threads} > \$homer_annotate_file" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			homer_annotate_file=$RESULT_PATH/annotate.txt
			homer_annotate_statistics_file=$RESULT_PATH/annotate_statistics.txt
			annotatePeaks.pl $RESULT_PATH/{wildcards.sample}.narrowPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} \\
			-go $RESULT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file

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
			printf "%s\\n" "for file in $RESULT_PATH/*.txt" | tee >(cat >&2)
			printf "%s\\n" "do" | tee >(cat >&2)
			printf "\\t%s\\n" "name=\$(basename \$file)" | tee >(cat >&2)
			printf "\\t%s\\n" "name={wildcards.sample}.narrowPeak.homer_\$name" | tee >(cat >&2)
			printf "\\t%s\\n" "mv \$file $RESULT_PATH/\$name" | tee >(cat >&2)
			printf "%s\\n" "done" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			for file in $RESULT_PATH/*.txt
			do
				name=$(basename $file)
				name={wildcards.sample}.narrowPeak.homer_$name
				mv $file $RESULT_PATH/$name
			done

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
		processed_column_List = []
		target_Path = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + peak_design + "/peak_calling/narrowpeak/"
		homer_DF = pandas.read_csv(output.narrowpeak_annotate, sep="\t", low_memory=False, index_col=False)
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
		homer_DF.to_csv(output.narrowpeak_annotate, sep='\t', index=False, header=True)


rule broadpeak_annotate:
	input:
		broadpeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		broadpeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bdg",
	output:
		broadpeak_annotate = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_annotate/broadpeak/{sample}/{sample}.broadPeak.homer_annotate.txt",
	priority: 996
	threads: PROCESSORS
	message: "broadpeak_annotate: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		bedgraph_List = []
		peak_basename = os.path.basename(input.broadpeak_bed)
		peak_begining = re.sub(".broadPeak.gz", "", peak_basename)
		peak_design = get_design(input.broadpeak_bed, GENOME)
		bedgraph_List.append(input.broadpeak_bdg)
		target_List = ["_POOLED_", "_OVERLAPPED_", "_IDR_", "_VS_"]
		for each_target in target_List:
			if each_target in peak_begining:
				sample_name_List = peak_begining.split(each_target)
				for each_sample in sample_name_List:
					#
					bedgraph_List.extend(glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + peak_design + "/peak_calling/broadpeak/*" + each_sample + "*.broadPeak.bdg"))
				else:
					pass
			else:
				pass
		else:
			pass
		bedgraph_List = list(set(bedgraph_List))
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_annotate/broadpeak/{wildcards.sample}
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_annotate/broadpeak/{wildcards.sample}
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH
			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "broadpeak_annotate: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.broadpeak_bed}" | tee >(cat >&2)
			declare -a bam_List=({bedgraph_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				index=$(($index+1))
				printf "INPUT%s: %s\\n" "${{index}}" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadpeak_annotate}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load bedtools/2.27.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load ucsc/373 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "module load homer/4.10.1 || exit 1" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load ucsc/373 || exit 1
			module load homer/4.10.1 || exit 1

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
			printf "%s\\n" "gunzip < {input.broadpeak_bed} > $RESULT_PATH/{wildcards.sample}.broadPeak.bed" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			gunzip < {input.broadpeak_bed} > $RESULT_PATH/{wildcards.sample}.broadPeak.bed

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
			printf "%s\\n" "homer_annotate_file=${{RESULT_PATH}}/annotate.txt" | tee >(cat >&2)
			printf "%s\\n" "homer_annotate_statistics_file=${{RESULT_PATH}}/annotate_statistics.txt" | tee >(cat >&2)
			printf "%s\\n" "annotatePeaks.pl $RESULT_PATH/{wildcards.sample}.broadPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} \\
			-go $RESULT_PATH -annStats \$homer_annotate_statistics_file -cpu {threads} > \$homer_annotate_file" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			printf "%s\\n" "for file in $RESULT_PATH/*.txt" | tee >(cat >&2)
			printf "%s\\n" "do" | tee >(cat >&2)
			printf "\\t%s\\n" "name=\$(basename \$file)" | tee >(cat >&2)
			printf "\\t%s\\n" "name={wildcards.sample}.broadPeak.homer_\$name" | tee >(cat >&2)
			printf "\\t%s\\n" "mv \$file $RESULT_PATH/\$name" | tee >(cat >&2)
			printf "%s\\n" "done" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			homer_annotate_file=$RESULT_PATH/annotate.txt
			homer_annotate_statistics_file=$RESULT_PATH/annotate_statistics.txt
			annotatePeaks.pl $RESULT_PATH/{wildcards.sample}.broadPeak.bed {GENOME} -bedGraph {bedgraph_List} -gsize {EFFECTIVE_GENOME_SIZE} \\
			-go $RESULT_PATH -annStats $homer_annotate_statistics_file -cpu {threads} > $homer_annotate_file
			for file in $RESULT_PATH/*.txt
			do
				name=$(basename $file)
				name={wildcards.sample}.broadPeak.homer_$name
				mv $file $RESULT_PATH/$name
			done

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
		processed_column_List = []
		target_Path = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + peak_design + "/peak_calling/broadpeak/"
		homer_DF = pandas.read_csv(output.broadpeak_annotate, sep="\t", low_memory=False, index_col=None)
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
		homer_DF.to_csv(output.broadpeak_annotate, sep='\t', index=False, header=True)
