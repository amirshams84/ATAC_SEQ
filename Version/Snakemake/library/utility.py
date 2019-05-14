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


def build_metadata_dict(metadata_file_Path):
	if is_file_exist(metadata_file_Path) is True:
		pass
	metadata_Dict = {}

	metadata_DF = pandas.read_table(
		metadata_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col="Name"
	)
	transposed_metadata_DF = metadata_DF.transpose()
	metadata_Dict = transposed_metadata_DF.to_dict()
	return metadata_Dict


def build_snakemake_awk(awk_String):
	"""
	"""
	snakemake_awk_String = ""
	dupliucate_string_List = ['\\', '{', '}']
	for each_letter in list(awk_String):
		if each_letter not in dupliucate_string_List:
			#
			snakemake_awk_String += each_letter
		else:
			snakemake_awk_String += each_letter * 2
	return snakemake_awk_String


def slugify(target_String):
	#replacing non-AlphNumeric with underscore
	slugified_String = re.sub('[^0-9a-zA-Z_-]+', '_', target_String)
	#slugified_String = slugified_String.replace("___", "__")
	return slugified_String

