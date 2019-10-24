import os
import sys
import pandas
import re
# ################################### UTILITY FUNCTIONS ########################## #


def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		raise Exception("The file Path: " + file_Path + " is not accessible or the premission to read is not granted!!!")
		raise Exception("ABORTIG!!!")
		sys.exit(2)
		return False


def is_path_exist(the_Path):
	"""
	"""
	if os.path.isdir(the_Path) and os.access(the_Path, os.W_OK) and os.access(the_Path, os.R_OK):
		return True
	else:
		raise Exception("The Path: " + the_Path + " is not accessible or the premission to read/Write is not granted!!!")
		raise Exception("ABORTIG!!!")
		sys.exit(2)
		return False


def fix_path(the_Path):
	"""
	"""
	if the_Path[-1] == "/":
		the_Path = the_Path[:-1]
	return the_Path


def write_string_down(the_String, the_file_Path):
	#
	f = open(the_file_Path, "w")
	f.write(the_String)
	f.close()
	return True


def build_sample_treatment_dict(treatment_file_Path, sample_column):
	"""
	"""
	if is_file_exist(treatment_file_Path) is True:
		pass
	sample_treatment_Dict = {}

	sample_treatment_DF = pandas.read_csv(
		treatment_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col=sample_column
	)
	transposed_sample_treatment_DF = sample_treatment_DF.transpose()
	sample_treatment_Dict = transposed_sample_treatment_DF.to_dict()
	return sample_treatment_Dict


def build_snakemake_awk(awk_String):
	"""
	"""
	snakemake_awk_String = ""
	dupliucate_string_List = ['\\', '{', '}']
	escape_character_List = ['$', '"']
	for each_letter in list(awk_String):
		if each_letter not in dupliucate_string_List:
			#
			snakemake_awk_String += each_letter
		elif each_letter in escape_character_List:
			#
			snakemake_awk_String += '\\' + each_letter

		else:
			snakemake_awk_String += each_letter * 2
	return snakemake_awk_String


def slugify(target_String):
	#replacing non-AlphNumeric with underscore
	slugified_String = re.sub('[^0-9a-zA-Z_-]+', '_', target_String)
	#slugified_String = slugified_String.replace("___", "__")
	return slugified_String


def get_cluster_info(parameters_List):
	"""
	"""
	for each_element in parameters_List:
		#
		if "--cluster=" not in each_element:
			#
			processors = 10
			memory = 10000
			return(processors, memory)
		if "--cluster=" in each_element:
			cluster_String = each_element
			cluster_List = cluster_String.split(" ")
			for each_cluster_element in cluster_List:
				if "--cpus-per-task=" in each_cluster_element:
					processors = int(each_cluster_element.split("=")[1])
				elif "--mem=" in each_cluster_element:
					memory = int(each_cluster_element.split("=")[1])
				else:
					pass
			else:
				pass
	else:
		pass
	return(processors, memory)


def build_metadata_dict(sample_treatment_Dict, treatment_column, treatment_List):
	"""
	"""
	metadata_Dict = {}
	for sample, sample_Dict in sample_treatment_Dict.items():
		#
		if sample_Dict[treatment_column] not in treatment_List:
			#Skip unwanted treatment
			continue
		if sample_Dict[treatment_column] not in metadata_Dict:
			#
			metadata_Dict[sample_Dict[treatment_column]] = {}
			metadata_Dict[sample_Dict[treatment_column]]["Case"] = []
			metadata_Dict[sample_Dict[treatment_column]]["Control"] = []
			#
			if sample_Dict["Type"] == "CASE":
				#
				metadata_Dict[sample_Dict[treatment_column]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				metadata_Dict[sample_Dict[treatment_column]]["Control"].append(sample)
			else:
				pass
		elif sample_Dict[treatment_column] in metadata_Dict:
			#
			if sample_Dict["Type"] == "CASE":
				#
				metadata_Dict[sample_Dict[treatment_column]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				metadata_Dict[sample_Dict[treatment_column]]["Control"].append(sample)
			else:
				pass
		else:
			pass
	else:
		pass
	return metadata_Dict

