import os,sys
import pandas
import glob
import zipfile
import math
import gzip
basic_statistic_Dict = {}


def get_chromosome(genome):
	chromosome_List = []
	if genome in ["hg38", "hg19"]:
		chr_List = list(range(1, 23, 1))

		for each_number in chr_List:
			chromosome_List.append("chr" + str(each_number))
		chromosome_List.extend(["chrX", "chrY", "chrM", "chrEBV"])
		return chromosome_List
	elif genome in ["mm9", "mm10"]:
		chr_List = list(range(1, 20, 1))
		for each_number in chr_List:
			chromosome_List.append("chr" + str(each_number))
		chromosome_List.extend(["chrX", "chrY", "chrM"])
		return chromosome_List
	else:
		pass


def get_design(file_path, genome):
	"""
	"""
	each_path = os.path.abspath(file_path)
	each_sample = each_path.split(genome + "/")[1].split("/")[0]
	return each_sample

###################OLD
def get_cutadapt_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	cutadapt_report_List = glob.glob(main_Path + "/" + genome + "/**/report/pre_process/*.cutadapt.txt", recursive=True)
	for each_file in cutadapt_report_List:
		#
		each_design = get_design(each_file, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		#
		each_sample = os.path.basename(each_file).split("_")[0]
		if each_sample not in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample] = {}
		else:
			pass
		#
		for line in open(each_file):
			if "Total read pairs processed:" in line:
				read_count = line.split("Total read pairs processed:")[1].replace(" ", "").replace(",", "").rstrip()
				break
			else:
				pass
		else:
			pass
		basic_statistic_Dict[each_design][each_sample]["Raw_Read_count"] = int(read_count)
	else:
		pass

	return True


def get_alignment_pbc_report_old(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	pbc_report_List = glob.glob(main_Path + "/**/report/alignment/*.PBC.txt", recursive=True)
	for each_file in pbc_report_List:
		each_design = get_design(each_file, genome)
		each_sample = os.path.basename(each_file).split(".PBC.txt")[0]
		each_pbc_DF = pandas.read_csv(each_file, sep="\t", index_col=False)
		basic_statistic_Dict[each_design][each_sample]["NRF"] = each_pbc_DF.loc[0]["NRF"]
		basic_statistic_Dict[each_design][each_sample]["PBC1"] = each_pbc_DF.loc[0]["PBC1"]
		basic_statistic_Dict[each_design][each_sample]["PBC2"] = each_pbc_DF.loc[0]["PBC2"]
	return True
##########################


def get_flowcell_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	cutadapt_report_List = glob.glob(main_Path + "/" + genome + "/**/report/pre_process/*.cutadapt.txt", recursive=True)
	for each_file in cutadapt_report_List:
		#
		each_design = get_design(each_file, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		#
		each_sample = os.path.basename(each_file).split("_")[0]
		each_flowcell = os.path.basename(each_file).split("_")[1].split("_S")[0]
		if "Flowcell_ID" not in basic_statistic_Dict[each_design][each_sample]:
			basic_statistic_Dict[each_design][each_sample]["Flowcell_ID"] = [each_flowcell]
		else:
			basic_statistic_Dict[each_design][each_sample]["Flowcell_ID"].append(each_flowcell)
	else:
		pass
	return True


def get_fastqc_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	fastqc_report_List = glob.glob(main_Path + "/" + genome + "/**/report/pre_process/*.R1.processed_fastqc.zip", recursive=True)
	for each_fastqc_zip in fastqc_report_List:
		each_design = get_design(each_fastqc_zip, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		#
		each_sample = os.path.basename(each_fastqc_zip).split(".R1.processed_fastqc.zip")[0]
		if each_sample not in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample] = {}
		else:
			pass
		#
		each_fastqc_directory = zipfile.ZipFile(each_fastqc_zip)
		for each_file in each_fastqc_directory.namelist():
			if not os.path.isdir(each_file):
				if each_file.endswith('fastqc_data.txt'):
					for line in each_fastqc_directory.open(each_file):
						if "Total Sequences".encode() in line:
							read_count = line.split("Total Sequences".encode())[1].replace("\t".encode(), "".encode()).rstrip().decode()
							break
						else:
							pass
		basic_statistic_Dict[each_design][each_sample]["Raw_Read_count"] = int(read_count)
		each_fastqc_directory.close()
	print("get_fastqc_report_DONE!!!")
	return True


def get_alignment_flagstat_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	flagstat_report_List = glob.glob(main_Path + "/" + genome + "/**/report/alignment/*.samtools.flagstat.txt", recursive=True)
	for each_file in flagstat_report_List:
		#
		each_design = get_design(each_file, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		each_sample = os.path.basename(each_file).split(".samtools.flagstat.txt")[0]
		if each_sample not in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample] = {}
		else:
			pass
		for line in open(each_file):
			if "total" in line:
				mapped_count = line.split("+")[0].rstrip()
			elif "duplicates" in line:
				duplicates_count = line.split("+")[0].rstrip()
			else:
				pass
		else:
			pass
		
		basic_statistic_Dict[each_design][each_sample]["Mapped_count"] = math.ceil(int(mapped_count) / 2)
		basic_statistic_Dict[each_design][each_sample]["Duplicate_count"] = math.ceil(int(duplicates_count) / 2)
		basic_statistic_Dict[each_design][each_sample]["Unique_count"] = math.ceil((int(mapped_count) - int(duplicates_count)) / 2)
	else:
		pass
	target_name = ""
	for each_design in basic_statistic_Dict:
		total_read = 0
		flowcell_List = []
		for each_sample in basic_statistic_Dict[each_design]:
			if "POOLED" in each_sample:
				target_name = each_sample
				continue
			else:
				total_read += basic_statistic_Dict[each_design][each_sample]["Raw_Read_count"]
				flowcell_List.extend(basic_statistic_Dict[each_design][each_sample]["Flowcell_ID"])
		else:
			pass
		if "POOLED" in target_name:
			#
			basic_statistic_Dict[each_design][target_name]['Raw_Read_count'] = total_read
			basic_statistic_Dict[each_design][target_name]['Flowcell_ID'] = list(set(flowcell_List))
		else:
			pass
	else:
		pass
	print("get_alignment_flagstat_report_Done!!!")

	return True


def get_post_alignment_flagstat_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	flagstat_report_List = glob.glob(main_Path + "/" + genome + "/**/report/post_alignment/*.processed.samtools.flagstat.txt", recursive=True)
	for each_file in flagstat_report_List:
		#
		each_design = get_design(each_file, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		each_sample = os.path.basename(each_file).split(".processed.samtools.flagstat.txt")[0]
		if each_sample not in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample] = {}
		else:
			pass
		for line in open(each_file):
			if "total" in line:
				mapped_count = line.split("+")[0].rstrip()
			elif "duplicates" in line:
				duplicates_count = line.split("+")[0].rstrip()
			else:
				pass
		else:
			pass
		
		basic_statistic_Dict[each_design][each_sample]["Filtered_Mapped_count"] = math.ceil(int(mapped_count) / 2)
		basic_statistic_Dict[each_design][each_sample]["Filtered_Mapped_Duplicate_count"] = math.ceil(int(duplicates_count) / 2)
	print("get_post_alignment_flagstat_report_DONE!!!")
	return True


def get_alignment_idxstats_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""

	chromosome_List = get_chromosome(genome)
	idxstats_report_List = glob.glob(main_Path + "/**/report/alignment/*.samtools.idxstats.txt", recursive=True)
	for each_file in idxstats_report_List:
		each_design = get_design(each_file, genome)
		each_sample = os.path.basename(each_file).split(".samtools.idxstats.txt")[0]
		each_idxstats_DF = pandas.read_csv(each_file, sep="\t", header=None, index_col=0)
		each_idxstats_DF = each_idxstats_DF.filter(items=chromosome_List, axis='index')
		for each_chromosme in chromosome_List:
			basic_statistic_Dict[each_design][each_sample][each_chromosme] = math.ceil(each_idxstats_DF.loc[each_chromosme, 2] / 2)
		basic_statistic_Dict[each_design][each_sample]["Mitochondrial_count"] = basic_statistic_Dict[each_design][each_sample]["chrM"]
	print("get_alignment_idxstats_report_DONE!!!")
	return True


def get_alignment_pbc_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	pbc_report_List = glob.glob(main_Path + "/**/report/post_alignment/*.processed.PBC.txt", recursive=True)
	for each_file in pbc_report_List:
		each_design = get_design(each_file, genome)
		each_sample = os.path.basename(each_file).split(".processed.PBC.txt")[0]
		each_pbc_DF = pandas.read_csv(each_file, sep="\t", index_col=False)
		basic_statistic_Dict[each_design][each_sample]["NRF"] = each_pbc_DF.loc[0]["NRF"]
		basic_statistic_Dict[each_design][each_sample]["PBC1"] = each_pbc_DF.loc[0]["PBC1"]
		basic_statistic_Dict[each_design][each_sample]["PBC2"] = each_pbc_DF.loc[0]["PBC2"]
	print("get_alignment_pbc_report_Done!!!")
	return True


def get_duplication_ratio(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	fastqc_report_List = glob.glob(main_Path + "/" + genome + "/**/report/post_alignment/*.processed_fastqc.zip", recursive=True)
	for each_fastqc_zip in fastqc_report_List:
		each_design = get_design(each_fastqc_zip, genome)
		if each_design not in basic_statistic_Dict:
			basic_statistic_Dict[each_design] = {}
		else:
			pass
		#
		each_sample = os.path.basename(each_fastqc_zip).split(".processed_fastqc.zip")[0]
		if each_sample not in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample] = {}
		else:
			pass
		#
		each_fastqc_directory = zipfile.ZipFile(each_fastqc_zip)
		table_string = ''
		for each_file in each_fastqc_directory.namelist():
			if not os.path.isdir(each_file):
				if each_file.endswith('fastqc_data.txt'):
					file_handle = each_fastqc_directory.open(each_file)
					for line in file_handle:
						if line.startswith(b'>>Sequence Duplication Levels'):
							for line in file_handle:
								table_string += line.decode()
								if line.startswith(b'>>END_MODULE'):
									break
				
		each_fastqc_directory.close()
		duplication_ratio_DF = pandas.DataFrame([x.split('\t') for x in table_string.split('\n')])
		duplication_ratio_DF.drop(duplication_ratio_DF.head(1).index, inplace=True)
		duplication_ratio_DF.drop(duplication_ratio_DF.tail(2).index, inplace=True)
		new_header = duplication_ratio_DF.iloc[0]
		duplication_ratio_DF = duplication_ratio_DF.iloc[1:]
		duplication_ratio_DF.columns = new_header
		for index, row in duplication_ratio_DF.iterrows():
			if row["Percentage of deduplicated"] == "0.0":
				basic_statistic_Dict[each_design][each_sample][row["#Duplication Level"] + "_Duplicates"] = float(0.0)
			else:

				basic_statistic_Dict[each_design][each_sample][row["#Duplication Level"] + "_Duplicates"] = math.ceil((float(row["Percentage of deduplicated"]) * basic_statistic_Dict[each_design][each_sample]["Duplicate_count"]) / 100)
			
	return True


def get_peak_count(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	peak_count_List = glob.glob(main_Path + "/**/peak_calling/narrowpeak/*.narrowPeak.gz", recursive=True)
	peak_header_List = ["chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit"]
	for each_file in peak_count_List:

		if "IDR" in each_file:
			continue
		elif "OVERLAPPED" in each_file:
			continue
		
		each_design = get_design(each_file, genome)
		each_sample = os.path.basename(each_file).split(".narrowPeak.gz")[0]
		
		each_peak_file_DF = pandas.read_csv(each_file, sep="\t", header=None, index_col=False, compression='gzip', error_bad_lines=False)
		each_peak_file_DF.columns = peak_header_List
		
		
		basic_statistic_Dict[each_design][each_sample]["Peak_count"] = each_peak_file_DF.shape[0]
		basic_statistic_Dict[each_design][each_sample]["Peak_average_signal"] = each_peak_file_DF["signalValue"].mean()
		basic_statistic_Dict[each_design][each_sample]["Peak_average_pValue"] = each_peak_file_DF["pValue"].mean()
		basic_statistic_Dict[each_design][each_sample]["Peak_average_qValue"] = each_peak_file_DF["qValue"].mean()

	return True


def get_frip_score(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	frip_score_List = glob.glob(main_Path + "/**/report/peak_calling/**/*.frip.txt", recursive=True)
	for each_file in frip_score_List:
		#
		each_frip_DF = pandas.read_csv(each_file, sep="\t", index_col=False)
		frip_ratio = float(each_frip_DF.loc[0]["read_in_peak_count"]) / float(each_frip_DF.loc[0]["total_read"])
		total_peak = each_frip_DF.loc[0]["total_peak"]

		each_design = get_design(each_file, genome)
		
		if "narrowPeak" in each_file:
			the_sample = os.path.basename(each_file).split(".narrowPeak.frip.txt")[0]
			if "IDR" in the_sample:
				peak_type = "narrowPeak_IDR"
			elif "OVERLAPPED" in the_sample:
				peak_type = "narrowPeak_OVERLAPPED"
			else:
				peak_type = "narrowPeak"
		elif "broadPeak" in each_file:
			the_sample = os.path.basename(each_file).split(".broadPeak.frip.txt")[0]
			if "IDR" in the_sample:
				peak_type = "broadPeak_IDR"
			elif "OVERLAPPED" in the_sample:
				peak_type = "broadPeak_OVERLAPPED"
			else:
				peak_type = "broadPeak"
		else:
			pass

		if "IDR" in peak_type:
			#
			sample_List = []
			for each_sample in the_sample.split("_IDR_"):
				sample_List.append(each_sample)
				basic_statistic_Dict[each_design][each_sample][peak_type + "_count(pValue<0.05)"] = total_peak
				basic_statistic_Dict[each_design][each_sample][peak_type + "_FRIP(pValue<0.05)"] = frip_ratio
			else:
				pass
			pooled_sample = "_POOLED_".join(sample_List)
			
			basic_statistic_Dict[each_design][pooled_sample][peak_type + "_count(pValue<0.05)"] = total_peak
			basic_statistic_Dict[each_design][pooled_sample][peak_type + "_FRIP(pValue<0.05)"] = frip_ratio
		elif "OVERLAPPED" in peak_type:
			#
			sample_List = []
			for each_sample in the_sample.split("_OVERLAPPED_"):
				sample_List.append(each_sample)
				basic_statistic_Dict[each_design][each_sample][peak_type + "_count"] = total_peak
				basic_statistic_Dict[each_design][each_sample][peak_type + "_FRIP"] = frip_ratio
			else:
				pass
			pooled_sample = "_POOLED_".join(sample_List)
			
			basic_statistic_Dict[each_design][pooled_sample][peak_type + "_count"] = total_peak
			basic_statistic_Dict[each_design][pooled_sample][peak_type + "_FRIP"] = frip_ratio
		else:
			basic_statistic_Dict[each_design][the_sample][peak_type + "_count(qValue<0.05)"] = total_peak
			basic_statistic_Dict[each_design][the_sample][peak_type + "_FRIP(qValue<0.05)"] = frip_ratio

	return True


def get_homer_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	homer_report_List = glob.glob(main_Path + "/**/peak_annotate/**/*.homer_annotate_statistics.txt", recursive=True)
	for each_homer_report in homer_report_List:
		#
		each_design = get_design(each_homer_report, genome)
		#
		if "narrowPeak" in each_homer_report:
			the_sample = os.path.basename(each_homer_report).split(".narrowPeak.homer_annotate_statistics.txt")[0]
			if "IDR" in the_sample:
				peak_type = "narrowPeak_IDR"
			elif "OVERLAPPED" in the_sample:
				peak_type = "narrowPeak_OVERLAPPED"
			else:
				peak_type = "narrowPeak"
		elif "broadPeak" in each_homer_report:
			the_sample = os.path.basename(each_homer_report).split(".broadPeak.homer_annotate_statistics.txt")[0]
			if "IDR" in the_sample:
				peak_type = "broadPeak_IDR"
			elif "OVERLAPPED" in the_sample:
				peak_type = "broadPeak_OVERLAPPED"
			else:
				peak_type = "broadPeak"
		else:
			pass

		each_homer_DF = pandas.read_csv(each_homer_report, sep="\t", index_col=False)
		for row_index in range(len(each_homer_DF)):
			#
			if each_homer_DF.loc[row_index]["Annotation"] == "Annotation":
				continue
			Annotation_peak_count = each_homer_DF.loc[row_index]["Annotation"] + "_Peak_count"
			Annotation_log = each_homer_DF.loc[row_index]["Annotation"] + "_Log2_Enrichment"

			if "IDR" in peak_type:
				#
				sample_List = []
				for each_sample in the_sample.split("_IDR_"):
					sample_List.append(each_sample)
					basic_statistic_Dict[each_design][each_sample][peak_type + "_" + Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
					basic_statistic_Dict[each_design][each_sample][peak_type + "_" + Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
				else:
					pass
				pooled_sample = "_POOLED_".join(sample_List)
				
				basic_statistic_Dict[each_design][pooled_sample][peak_type + "_" + Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
				basic_statistic_Dict[each_design][pooled_sample][peak_type + "_" + Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
			elif "OVERLAPPED" in peak_type:
				#
				sample_List = []
				for each_sample in the_sample.split("_OVERLAPPED_"):
					sample_List.append(each_sample)
					basic_statistic_Dict[each_design][each_sample][peak_type + "_" + Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
					basic_statistic_Dict[each_design][each_sample][peak_type + "_" + Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
				else:
					pass
				pooled_sample = "_POOLED_".join(sample_List)
				
				basic_statistic_Dict[each_design][pooled_sample][peak_type + "_" + Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
				basic_statistic_Dict[each_design][pooled_sample][peak_type + "_" + Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
			else:
				basic_statistic_Dict[each_design][the_sample][peak_type + "_" + Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
				basic_statistic_Dict[each_design][the_sample][peak_type + "_" + Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
	return True


def get_idr_report(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	idr_report = glob.glob(main_Path + "/**/report/peak_calling/narrowpeak/*_IDR.txt", recursive=True)
	for each_file in idr_report:
		#
		each_design = get_design(each_file, genome)
		for line in open(each_file):
			if "Number of peaks passing IDR cutoff" in line:
				idr_value = line.split(" - ")[1].split(" (")[0].split("/")[0]
				break
			else:
				pass
		else:
			pass
		for each_sample in basic_statistic_Dict[each_design]:
			basic_statistic_Dict[each_design][each_sample]["IDR_Peak_count(<0.05)"] = idr_value
	return True


def get_homer_report_old(main_Path, genome, basic_statistic_Dict):
	"""
	"""
	homer_report = glob.glob(main_Path + "/**/peak_annotate/narrowpeak/**/*.homer_annotate_statistics.txt", recursive=True)
	for each_file in homer_report:
		if "IDR" in each_file:
			continue
		elif "OVERLAPPED" in each_file:
			continue
		else:
			each_design = get_design(each_file, genome)
			each_sample = os.path.basename(each_file).split(".narrowPeak.homer_annotate_statistics.txt")[0]
			
			each_homer_DF = pandas.read_csv(each_file, sep="\t", header=0, index_col=False)
			
			for row_index in range(len(each_homer_DF)):
				#
				if each_homer_DF.loc[row_index]["Annotation"] == "Annotation":
					continue
				Annotation_peak_count = each_homer_DF.loc[row_index]["Annotation"] + "_Peak_count"
				Annotation_log = each_homer_DF.loc[row_index]["Annotation"] + "_Log2_Enrichment"
				if Annotation_peak_count not in basic_statistic_Dict[each_design][each_sample]:
					basic_statistic_Dict[each_design][each_sample][Annotation_peak_count] = float(each_homer_DF.loc[row_index]['Number of peaks'])
				if Annotation_log not in basic_statistic_Dict[each_design][each_sample]:
					basic_statistic_Dict[each_design][each_sample][Annotation_log] = float(each_homer_DF.loc[row_index]['Log2 Enrichment'])
			
	return True


get_fastqc_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)
get_flowcell_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)

get_alignment_flagstat_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)

get_post_alignment_flagstat_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)

get_alignment_pbc_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)


get_frip_score(sys.argv[1], sys.argv[2], basic_statistic_Dict)


#get_peak_count(sys.argv[1], sys.argv[2], basic_statistic_Dict)
#get_idr_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)

get_homer_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)
print(basic_statistic_Dict)
get_alignment_idxstats_report(sys.argv[1], sys.argv[2], basic_statistic_Dict)
get_duplication_ratio(sys.argv[1], sys.argv[2], basic_statistic_Dict)

###########################################
column_List = [
	"Flowcell_ID", "Raw_Read_count", "Mapped_count", "Mitochondrial_count", "Duplicate_count", "Filtered_Mapped_count", "Filtered_Mapped_Duplicate_count", "NRF", "PBC1", "PBC2",
	"narrowPeak_count(qValue<0.05)", "narrowPeak_FRIP(qValue<0.05)",
	"narrowPeak_OVERLAPPED_count", "narrowPeak_OVERLAPPED_FRIP",
	"narrowPeak_IDR_count(pValue<0.05)", "narrowPeak_IDR_FRIP(pValue<0.05)",
	"broadPeak_count(qValue<0.05)", "broadPeak_FRIP(qValue<0.05)",
	"broadPeak_OVERLAPPED_count", "broadPeak_OVERLAPPED_FRIP",
	"broadPeak_IDR_count(pValue<0.05)", "broadPeak_IDR_FRIP(pValue<0.05)"
]

Duplicate_List = ["Unique_count", "1_Duplicates", "2_Duplicates", "3_Duplicates", "4_Duplicates", "5_Duplicates", ">10_Duplicates", ">100_Duplicates", ">500_Duplicates", ">1k_Duplicates", ">10k+_Duplicates"]
chromosome_List = get_chromosome(sys.argv[2])

#HomerList = ["3UTR", "5UTR", "Promoter", "Intergenic", "CpG-Island", "miRNA", "ncRNA", "TTS", "pseudo", "Exon", "Intron", "snoRNA", "scRNA", "rRNA", "Retroposon", "RC?", "RNA", "LINE", "srpRNA"]
#HomerList.extend(["SINE", "tRNA", "DNA?", "DNA", "LTR?", "Low_complexity", "LTR", "Simple_repeat", "snRNA", "SINE?", "Unknown", "Satellite"])
HomerList = ["narrowPeak_5UTR", "narrowPeak_Promoter", "narrowPeak_3UTR", "narrowPeak_Exon", "narrowPeak_Intron", "narrowPeak_Intergenic"]
HomerList.extend(["broadPeak_5UTR", "broadPeak_Promoter", "narrowPeak_3UTR", "broadPeak_Exon", "broadPeak_Intron", "broadPeak_Intergenic"])

HomerList.extend(["narrowPeak_IDR_5UTR", "narrowPeak_IDR_Promoter", "narrowPeak_IDR_3UTR", "narrowPeak_IDR_Exon", "narrowPeak_IDR_Intron", "narrowPeak_IDR_Intergenic"])
HomerList.extend(["broadPeak_IDR_5UTR", "broadPeak_IDR_Promoter", "broadPeak_IDR_3UTR", "broadPeak_IDR_Exon", "broadPeak_IDR_Intron", "broadPeak_IDR_Intergenic"])

HomerList.extend(["narrowPeak_OVERLAPPED_5UTR", "narrowPeak_OVERLAPPED_Promoter", "narrowPeak_OVERLAPPED_3UTR", "narrowPeak_OVERLAPPED_Exon", "narrowPeak_OVERLAPPED_Intron", "narrowPeak_OVERLAPPED_Intergenic"])
HomerList.extend(["broadPeak_OVERLAPPED_5UTR", "broadPeak_OVERLAPPED_Promoter", "broadPeak_OVERLAPPED_3UTR", "broadPeak_OVERLAPPED_Exon", "broadPeak_OVERLAPPED_Intron", "broadPeak_OVERLAPPED_Intergenic"])

#HomerList.extend(["3UTR", "5UTR", "Promoter", "Exon", "Intron", "Intergenic", "CpG-Island"])



#HomerList = list(set(HomerList))
HomerList_new = []
for each_element in HomerList:
	HomerList_new.append(each_element + "_Log2_Enrichment")
	HomerList_new.append(each_element + "_Peak_count")



column_List.extend(chromosome_List)
column_List.extend(Duplicate_List)
column_List.extend(HomerList_new)

dataframe_List = []
design_List = []
for each_design in basic_statistic_Dict:
	#
	design_List.append(each_design)
	each_design_DF = pandas.DataFrame(basic_statistic_Dict[each_design]).T
	
	#each_design_DF = each_design_DF[column_List]
	each_design_DF = each_design_DF.reindex(columns=column_List)
	each_design_DF.index.name = "Sample_name"

	dataframe_List.append(each_design_DF)
basic_statistic_DF = pandas.concat(dataframe_List, ignore_index=False, keys=design_List, names=["Design"])

basic_statistic_DF.reset_index(level=['Design'], inplace=True, drop=False)
basic_statistic_DF = basic_statistic_DF.T
basic_statistic_DF.index.name = "Sample_name"
basic_statistic_DF.to_csv("./ATAC_QC_stats.txt", sep="\t", header=True, index=True, float_format='%.3f')
