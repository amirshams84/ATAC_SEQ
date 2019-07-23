# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for POST post_alignment
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
PLATFORM = config_data_Dict["PLATFORM"].lower()
LAYOUT = config_data_Dict["LAYOUT"].lower()
GENOME = config_data_Dict["GENOME"].lower()
READ_LENGTH = config_data_Dict["READ_LENGTH"]
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
# ++++++++++++++++++++++++++++++++++++
#POST_ALIGNMENT
config_post_alignment_Dict = config["POST_ALIGNMENT"]

#AWK_CHROMOSOME_FILTER = utility.build_snakemake_awk(config_post_alignment_Dict["FILTER_CHROMOSOME"])
#print(AWK_CHROMOSOME_FILTER)
#AWK_ATAC_STRAND_SHIFT = utility.build_snakemake_awk(config_post_alignment_Dict["ATAC_STRAND_SHIFT"])
#print(AWK_ATAC_STRAND_SHIFT)
# ------------------------------------
# ################################### WILDCARDS ################################ #

alignment_List = []
post_alignment_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		post_alignment_List
# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM3: POST-ALIGNMENT
#+++++++++++++++++++++++++++++

rule Post_Alignment:
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
	message: "Post_Alignment: {wildcards.design}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			module load samtools/1.9 || exit 1
			module load bedtools/2.27.1 || exit 1
			module load deeptools/3.1.3 || exit 1
			module load fastqc/0.11.8 || exit 1
			module load qualimap/2.2.1 || exit 1
			module load ucsc/373 || exit 1
			module load R/3.5 || exit 1
			#
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/post_alignment
			mkdir -p $OUT_PATH
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/post_alignment
			mkdir -p $QC_PATH
			#
			#
			if [ ! -f {Script_Path}/{GENOME}.blacklist.bed.gz ]; then
				wget {config_reference_Dict[BLACK_LIST]} -O {Script_Path}/{GENOME}.blacklist.bed.gz
			fi
			#
			
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			
			samtools view --threads {threads} -h -F 1804 -f 2 -q 30 {input.bam} | awk 'BEGIN{{OFS=FS}}{{if ( $3 != \"chrM\" && $3 != \"chrUn\" && $3 != \"chrEBV\" && $3 !~ /chrUn/ && $3 !~ /random/ ) print $0}}' | samtools view --threads {threads} -Shb - > $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.bam
			bedtools intersect -v -abam $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.bam -b <(zcat -f {Script_Path}/{GENOME}.blacklist.bed.gz ) > $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.blkfilt.bam
			samtools sort --threads {threads} -O bam $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.blkfilt.bam -o {output.processed_bam}
			samtools index -@ {threads} -b {output.processed_bam}
			rm -rf $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.bam $OUT_PATH/{wildcards.sample}.mapped.dupfilt.chrfilt.blkfilt.bam
			#
			#ATAC_SEQ
			bedtools bamtobed -i {output.processed_bam} | awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} {{if ($6 == \"+\") {{$2 = $2 + 4}} else if ($6 == \"-\") {{$3 = $3 - 5}} print $0}}' > {output.processed_bed}.tmp
			#
			LC_COLLATE=C sort -k1,1 -k2,2n {output.processed_bed}.tmp > {output.processed_bed}.sorted
			bgzip -c {output.processed_bed}.sorted > {output.processed_bed}
			bedToBigBed {output.processed_bed}.sorted {config_reference_Dict[CHROM_SIZE]} {output.processed_bigbed}

			tabix -f -p bed {output.processed_bed}
			rm -rf {output.processed_bed}.tmp {output.processed_bed}.sorted

			bamCoverage --bam {output.processed_bam} --outFileName {output.processed_bigwig} --binSize 5 --normalizeUsing RPGC --extendReads 150 --outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]}

			samtools flagstat --threads {threads} {output.processed_bam} > $QC_PATH/{wildcards.sample}.processed.samtools.txt
			fastqc -o $QC_PATH -f bam --threads {threads} {output.processed_bam}

			printf "%s\\n" "#!/bin/bash" > $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "module load python/2.7" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "module load R/3.5" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "python {Python_Script_path}/bamQC.py --infile {output.processed_bam} --outfile $QC_PATH/{wildcards.sample}_processed_PBC.txt --cores {threads}" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "Rscript {R_Script_path}/run_spp_nodups.R -c={output.processed_bam} -odir=$QC_PATH/ -rf -savp > $QC_PATH/{wildcards.sample}_processed_SPP.txt" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			printf "%s\\n" "echo 'Pipeline execution successfully finished at: '$(date)" >> $QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			cd $QC_PATH
			sbatch --mem=300G --cpus-per-task=52 --partition=largemem --time=3-00:00:00 ./$QC_PATH/{wildcards.design}_{wildcards.sample}_PBC_SPP.sh
			
		""")

