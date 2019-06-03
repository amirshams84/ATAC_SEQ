# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Analysis
# snakemake --snakefile peak_analysis.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile peak_analysis.py --configfile Yoko.json --rulegraph | dot -Tsvg > peak_analysis.svg
# ################################### IMPORT ##################################### #


import os
import sys
import glob
import itertools
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


def get_narrowpeak_signal(wildcards):
	"""
	"""
	narrowpeak_signal_List = []
	narrowpeak_signal_List = glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/*.narrowPeak.bigwig".format(design=wildcards.design), recursive=True)
	narrowpeak_signal_List = sorted(narrowpeak_signal_List, key=os.path.getmtime)
	return narrowpeak_signal_List


def get_narrowpeak_bed(wildcards):
	"""
	"""
	narrowpeak_bed_List = []
	narrowpeak_bed_List = glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/*.narrowPeak.gz".format(design=wildcards.design), recursive=True)
	narrowpeak_bed_List = sorted(narrowpeak_bed_List, key=os.path.getmtime)
	return narrowpeak_bed_List


def get_bam(wildcards):
	"""
	"""
	bam_List = []
	bam_List = glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/*.bam".format(design=wildcards.design), recursive=True)
	bam_List = sorted(bam_List, key=os.path.getmtime)
	return bam_List
# ################################### CONFIGURATION ############################## #

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
# ################################### WILDCARDS ################################ #



peak_analysis_List = []

for design in design_Dict:
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.matrix.npz".format(design=design))
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.fingerprint.pdf".format(design=design))
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.signal_enrichment.flag".format(design=design))
else:
	pass

# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_analysis_List
# ################################### PIPELINE RULES ########################## #


rule Signal_Correlations:
	"""
	"""
	input:
		narrowpeak_signal_List = get_narrowpeak_signal,
	output:
		signal_correlation_matrix = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.matrix.npz",
		signal_correlation_pearson_heatmap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.pearson.heatmap.pdf",
		signal_correlation_spearman_heatmap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.spearman.heatmap.pdf",
		signal_correlation_pca = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.pca.pdf",
	priority: 996
	threads: PROCESSORS
	message: "Signal_Correlations: {wildcards.design}"
	resources:
		mem_mb = MEMORY
	run:
		narrowpeak_signal_label_String = ""
		narrowpeak_signal_label_List = []
		for each_signal in input.narrowpeak_signal_List:
			narrowpeak_signal_label_List.append(os.path.basename(each_signal).split(".")[0])
		narrowpeak_signal_label_String = " ".join(narrowpeak_signal_label_List)
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
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			multiBigwigSummary bins --bwfiles {input.narrowpeak_signal_List} --labels {narrowpeak_signal_label_String} --binSize 10000 --distanceBetweenBins 0 --numberOfProcessors {threads} --outFileName {output.signal_correlation_matrix}
			plotCorrelation -in {output.signal_correlation_matrix} --corMethod pearson --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --skipZeros --plotFileFormat pdf --plotHeight 12 --plotWidth 12 --removeOutliers -o {output.signal_correlation_pearson_heatmap}
			plotCorrelation -in {output.signal_correlation_matrix} --corMethod spearman --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --skipZeros --plotFileFormat pdf --plotHeight 12 --plotWidth 12 --removeOutliers -o {output.signal_correlation_spearman_heatmap}
			plotPCA --corData {output.signal_correlation_matrix} --plotFileFormat pdf --plotFile {output.signal_correlation_pca}

			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")


rule Fingerprint:
	"""
	"""
	input:
		bam_List = get_bam
	output:
		fingerprint_pdf = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.fingerprint.pdf",
		histogram_pdf = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.histogram.pdf",
		coverage_pdf = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.coverage.pdf",
	priority: 996
	threads: PROCESSORS
	message: "Fingerprint: {wildcards.design}"
	resources:
		mem_mb = MEMORY
	run:
		bam_label_String = ""
		bam_label_List = []
		for each_bam in input.bam_List:
			bam_label_List.append(os.path.basename(each_bam).split(".")[0])
		bam_label_String = " ".join(bam_label_List)
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
			printf "%s\\n" "#" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			printf "%s\\n" "#" | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			plotFingerprint --bamfiles {input.bam_List} --labels {bam_label_String} --binSize 500 --plotFileFormat pdf --skipZeros --numberOfProcessors {threads} --plotFile {output.fingerprint_pdf}
			bamPEFragmentSize --bamfiles {input.bam_List} --binSize 1000 --samplesLabel {bam_label_String} --plotFileFormat pdf --numberOfProcessors {threads} --histogram {output.histogram_pdf}
			plotCoverage --bamfiles {input.bam_List} --labels {bam_label_String} --skipZeros --plotFileFormat pdf --numberOfProcessors {threads} --plotFile {output.coverage_pdf}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "----------------------------------------------------------------------------" | tee >(cat >&2)
		""")

rule Signal_Enrichment:
	"""
	"""
	input:
		narrowpeak_signal_List = get_narrowpeak_signal,
		narrowpeak_bed_List = get_narrowpeak_bed,
	output:
		signal_enrichment_flag = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}.narrowPeak.signal_enrichment.flag",
		
	priority: 996
	threads: PROCESSORS
	message: "Signal_Enrichment: {wildcards.design}"
	resources:
		mem_mb = MEMORY
	run:
		for each_signal, each_bed in zip(input.narrowpeak_signal_List, input.narrowpeak_bed_List):
			each_label = os.path.basename(each_signal).split(".")[0]
			each_signal_enrichment_matrix = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}".format(design=design) + each_label + ".narrowPeak.signal_enrichment.matrix.npz",
			each_signal_enrichment_reference_heatmap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}".format(design=design) + each_label + ".narrowPeak.signal_enrichment_reference.heatmap.pdf",
			each_signal_enrichment_region_heatmap = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_analysis/narrowpeak/{design}".format(design=design) + each_label + ".narrowPeak.signal_enrichment_region.heatmap.pdf",
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
				OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_analysis/narrowpeak
				mkdir -p $OUT_PATH
				gunzip < {each_bed} > $OUT_PATH/narrowpeak.bed
				echo "computeMatrix reference-point --referencePoint TSS -S {each_signal} -R $OUT_PATH/narrowpeak.bed --outFileName {each_signal_enrichment_matrix} --numberOfProcessors {threads} -a 1500 -b 1500 -bs 10"
				computeMatrix reference-point --referencePoint TSS -S {each_signal} -R $OUT_PATH/narrowpeak.bed --outFileName {each_signal_enrichment_matrix} --numberOfProcessors {threads} -a 1500 -b 1500 -bs 10
				plotHeatmap -m {each_signal_enrichment_matrix} -out {each_signal_enrichment_reference_heatmap} --sortUsing max --plotFileFormat pdf
				computeMatrix scale-regions -S {each_signal} -R $OUT_PATH/narrowpeak.bed --outFileName {each_signal_enrichment_matrix} --numberOfProcessors {threads} -a 1500 -b 1500 -bs 10
				plotHeatmap -m {each_signal_enrichment_matrix} -out {each_signal_enrichment_region_heatmap} --sortUsing max --plotFileFormat pdf
			""")
		else:
			pass
		shell("""
			touch {output.signal_enrichment_flag}
		""")






