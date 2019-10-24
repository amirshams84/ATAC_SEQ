import os
import sys



def convert_narrowpeak_to_hammock(narrowpeak_file_Path, hammock_file_path):
	"""
	"""
	hammock_String = ""
	f = open(narrowpeak_file_Path, "rU")
	iterator = 1
	for i in f:
		i = i.rstrip()
		line = i.split("\t")

		hammock_String += line[0] + '\t' + line[1] + '\t' + line[2] + '\tscorelst:[' + line[6] + ',' + line[7] + ',' + line[8] + '],id:' + iterator
		iterator += 1
		if len(line[3]) > 1:
			hammock_String += 'name:"' + line[3] + '",'
		if line[5] != '.':
			hammock_String += 'strand:"' + line[5] + '",'
		if line[9] != '-1':
			hammock_String += 'sbstroke:[' + line[9] + ']'
		hammock_String += '\n'
	else:
		pass
	o = open(hammock_file_path, 'w')
	o.write(hammock_String)
	o.close()
	print("Conversion was successfull")
	return True
