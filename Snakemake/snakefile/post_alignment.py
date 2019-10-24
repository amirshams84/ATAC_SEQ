# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Project: ENCODE ATAC_SEQ
# Aim: Snakemake workflow for post_alignment
# snakemake --snakefile post_alignment.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile post_alignment.py --configfile Encode.json --rulegraph | dot -Tsvg > post_alignment.svg
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
# ------------------------------------
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
# ++++++++++++++++++++++++++++++++++++
#POST_ALIGNMENT
config_post_alignment_Dict = config["POST_ALIGNMENT"]
# ------------------------------------
# ################################### WILDCARDS ################################ #

alignment_List = []
post_alignment_List = []
for sample, sample_Dict in sample_treatment_Dict.items():
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict[TREATMENT_COLUMN], sample=sample))

for design in metadata_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(metadata_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		post_alignment_List
# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM3: POST-ALIGNMENT
#+++++++++++++++++++++++++++++
rule post_alignment:
	"""
	"""
	input:
		bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
		bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
	output:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
		processed_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bb",
		processed_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bigwig",
	priority: 993
	threads: PROCESSORS
	message: "post_alignment: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/post_alignment
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/post_alignment
			mkdir -p $REPORT_PATH
			
			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			if [ ! -f {Template_Path}/{GENOME}.blacklist.bed.gz ]; then
				wget {config_reference_Dict[BLACK_LIST]} -O {Template_Path}/{GENOME}.blacklist.bed.gz
			fi

			AWK_ALIGNMENT_FILTER="awk \'BEGIN{{OFS=FS}}{{if ( \$3 != \\\"chrM\\\" && \$3 != \\\"chrUn\\\" && \$3 !~ /chrUn/ && \$3 !~ /random/ ) print \$0}}\'"
			AWK_ATAC_SEQ_FILTER="awk \'BEGIN {{FS=\\\"\\\\t\\\"; OFS=\\\"\\\\t\\\"}} {{if (\$6 == \\\"+\\\") {{\$2 = \$2 + 4}} else if (\$6 == \\\"-\\\") {{\$3 = \$3 - 5}} print \$0}}\'"

			##
			#
			printf "%s\\n" "###################################- JOB INFO -################################" | tee >(cat >&2)
			printf "%s\\n" "post_alignment: {EXPERIMENT}|{TITLE}|{GENOME}|{wildcards.design}|{wildcards.sample}" | tee >(cat >&2)
			printf "%s\\n" "###################################- INPUT/OUTPUT -############################" | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.bam}" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.bam_index}" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.processed_bam}" | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.processed_bam_index}" | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.processed_bed}" | tee >(cat >&2)
			printf "OUTPUT4: %s\\n" "{output.processed_bed_index}" | tee >(cat >&2)
			printf "OUTPUT5: %s\\n" "{output.processed_bigbed}" | tee >(cat >&2)
			printf "OUTPUT6: %s\\n" "{output.processed_bigwig}" | tee >(cat >&2)
			printf "OUTPUT7: %s\\n" "$REPORT_PATH/{wildcards.sample}.processed.samtools.flagstat.txt" | tee >(cat >&2)
			printf "OUTPUT8: %s\\n" "$REPORT_PATH/{wildcards.sample}.processed.samtools.idxstats.txt" | tee >(cat >&2)
			printf "OUTPUT9: %s\\n" "$REPORT_PATH/{wildcards.sample}.processed.PBC.txt" | tee >(cat >&2)
			printf "OUTPUT10: %s\\n" "$REPORT_PATH/{wildcards.sample}.processed.SPP.txt" | tee >(cat >&2)
			printf "OUTPUT11: %s\\n" "$REPORT_PATH/{wildcards.sample}.fragment_length_distribution.txt" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "module load samtools/1.9 || exit 1" | tee >(cat >&2)
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
			module load bedtools/2.27.1 || exit 1
			module load ucsc/373 || exit 1
			module load fastqc/0.11.8 || exit 1
			module load deeptools/3.1.3 || exit 1

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
			printf "%s\\n" "samtools view --threads {threads} -h -F 780 -f 2 -q 30 {input.bam} | $AWK_ALIGNMENT_FILTER | \\
			samtools view --threads {threads} -Shb - > $RESULT_PATH/{wildcards.sample}.nomito.bam" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			samtools view --threads {threads} -h -F 780 -f 2 -q 30 {input.bam} | \\
			awk 'BEGIN{{OFS=FS}}{{if ( $3 != \"chrM\" && $3 != \"chrUn\" && $3 !~ /chrUn/ && $3 !~ /random/ ) print $0}}' | \\
			samtools view --threads {threads} -Shb - > $RESULT_PATH/{wildcards.sample}.nomito.bam

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
			printf "%s\\n" "bedtools intersect -v -abam $RESULT_PATH/{wildcards.sample}.nomito.bam \\
			-b <(zcat -f {Template_Path}/{GENOME}.blacklist.bed.gz ) > $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -O bam $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam -o {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf $RESULT_PATH/{wildcards.sample}.nomito.bam $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bedtools intersect -v -abam $RESULT_PATH/{wildcards.sample}.nomito.bam \\
			-b <(zcat -f {Template_Path}/{GENOME}.blacklist.bed.gz ) > $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam
			samtools sort --threads {threads} -O bam $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam -o {output.processed_bam}
			samtools index -@ {threads} -b {output.processed_bam}
			rm -rf $RESULT_PATH/{wildcards.sample}.nomito.bam $RESULT_PATH/{wildcards.sample}.nomito.blkfilt.bam

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
			printf "%s\\n" "bedtools bamtobed -i {output.processed_bam} | $AWK_ATAC_SEQ_FILTER > {output.processed_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.processed_bed}.tmp > {output.processed_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.processed_bed}.sorted > {output.processed_bed}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bedtools bamtobed -i {output.processed_bam} | \\
			awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} {{if ($6 == \"+\") {{$2 = $2 + 4}} else if ($6 == \"-\") {{$3 = $3 - 5}} print $0}}' > {output.processed_bed}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.processed_bed}.tmp > {output.processed_bed}.sorted
			bgzip -c {output.processed_bed}.sorted > {output.processed_bed}

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
			printf "%s\\n" "bedToBigBed {output.processed_bed}.sorted {config_reference_Dict[CHROM_SIZE]} {output.processed_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.processed_bed}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.processed_bed}.tmp {output.processed_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bedToBigBed {output.processed_bed}.sorted {config_reference_Dict[CHROM_SIZE]} {output.processed_bigbed}
			tabix -f -p bed {output.processed_bed}
			rm -rf {output.processed_bed}.tmp {output.processed_bed}.sorted

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
			printf "%s\\n" "bamCoverage --bam {output.processed_bam} --outFileName {output.processed_bigwig} --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
			--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]}" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			bamCoverage --bam {output.processed_bam} --outFileName {output.processed_bigwig} --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
			--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]}

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
			printf "%s\\n" "samtools flagstat --threads {threads} {output.processed_bam} > $REPORT_PATH/{wildcards.sample}.processed.samtools.flagstat.txt" | tee >(cat >&2)
			printf "%s\\n" "samtools idxstats --threads {threads} {output.processed_bam} > $REPORT_PATH/{wildcards.sample}.processed.samtools.idxstats.txt" | tee >(cat >&2)
			printf "%s\\n" "fastqc -o $REPORT_PATH -f bam --threads {threads} {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.processed_bam} --outfile $REPORT_PATH/{wildcards.sample}.processed.PBC.txt --cores {threads}" | tee >(cat >&2)
			printf "%s\\n" "Rscript {R_Script_path}/run_spp_nodups.R -c={output.processed_bam} -odir=$REPORT_PATH/ -rf -savp > $REPORT_PATH/{wildcards.sample}.processed.SPP.txt" | tee >(cat >&2)
			printf "%s\\n" "bash {Bash_Script}/fragment_length_distribution.sh {output.processed_bam} $REPORT_PATH/{wildcards.sample}.fragment_length_distribution.txt" | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "" | tee >(cat >&2)
			#
			##
			start_time="$(date -u +%s)"

			samtools flagstat --threads {threads} {output.processed_bam} > $REPORT_PATH/{wildcards.sample}.processed.samtools.flagstat.txt
			samtools idxstats --threads {threads} {output.processed_bam} > $REPORT_PATH/{wildcards.sample}.processed.samtools.idxstats.txt
			fastqc -o $REPORT_PATH -f bam --threads {threads} {output.processed_bam}

			printf "%s\\n" "#!/bin/bash" > $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "module load python/2.7" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "module load R/3.5" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.processed_bam} --outfile $REPORT_PATH/{wildcards.sample}.processed.PBC.txt --cores {threads}" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "Rscript {R_Script_path}/run_spp_nodups.R -c={output.processed_bam} -odir=$REPORT_PATH/ -rf -savp > $REPORT_PATH/{wildcards.sample}.processed.SPP.txt" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "bash {Bash_Script}/fragment_length_distribution.sh {output.processed_bam} $REPORT_PATH/{wildcards.sample}.fragment_length_distribution.txt" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			printf "%s\\n" "echo 'Pipeline execution successfully finished at: '\$(date)" >> $REPORT_PATH/{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh
			cd $REPORT_PATH
			sbatch --mem=300G --cpus-per-task=52 --partition=largemem --time=1-00:00:00 ./{wildcards.design}_{wildcards.sample}_processed_PBC_SPP.sh

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


