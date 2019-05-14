shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-17-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for NGS DATA Quality Control
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --debug-dag --cores=50
# snakemake --snakefile atac_seq.py --configfile atac_seq.json --rulegraph | dot -Tsvg > atac_seq.svg
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
	if is_path_exist(the_Path) is True:
		pass

	if the_Path[-1] == "/":
		the_Path = the_Path[:-1]
	return the_Path


def write_string_down(the_String, the_file_Path):
	#
	f = open(the_file_Path, "w")
	f.write(the_String)
	f.close()
	return True
# ################################### MAIN FUNCTIONS ############################# #


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


def get_forward_fastq(wildcards):
	return glob.glob(DATADIR + "/" + wildcards.sample + "*R1*" + metadata_Dict[wildcards.sample]["Extension"])


def get_case_bam(wildcards):
	processed_bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				processed_bam_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{case}.processed.bam".format(design=wildcards.design, case=sample))
	return processed_bam_List


def build_sidebar_html_string(pipeline_Dict, reference_Path, file_Path, pipeline_title, rule_title, rule_subtitle, report_type):
	#
	sidebar_String = ''
	for each_rule_title in pipeline_Dict:
		if each_rule_title == rule_title:
			#active_mode
			sidebar_String += '''
				<li class="active">
					<a data-toggle="collapse" href="#''' + each_rule_title + '''" aria-expanded="true">
						<i class="ti-panel"></i>
							<p>''' + each_rule_title + '''<b class="caret"></b></p>
					</a>
					<div class="collapse in" id="''' + each_rule_title + '''">
						<ul class="nav">
			'''
			for each_rule_subtitle in pipeline_Dict[each_rule_title]:
				#
				if each_rule_subtitle == rule_subtitle:
					#active_mode
					sidebar_String += '''
							<li class="active">
								<a href="''' + reference_Path + '''/''' + pipeline_title + '''_''' + each_rule_title + '''_''' + each_rule_subtitle + '''_''' + report_type + '''.html">
									<span class="sidebar-mini">''' + each_rule_subtitle[0].upper() + '''</span>
									<span class="sidebar-normal">''' + each_rule_subtitle.upper() + '''</span>
								</a>
							</li>
					'''
				elif each_rule_subtitle != rule_subtitle:
					#deactive_mode
					sidebar_String += '''
							<li>
								<a href="''' + reference_Path + '''/''' + pipeline_title + '''_''' + each_rule_title + '''_''' + each_rule_subtitle + '''_''' + report_type + '''.html">
									<span class="sidebar-mini">''' + each_rule_subtitle[0].upper() + '''</span>
									<span class="sidebar-normal">''' + each_rule_subtitle.upper() + '''</span>
								</a>
							</li>
					'''
				else:
					pass
		elif each_rule_title != rule_title:
			#deactive_mode
			sidebar_String += '''
				<li>
					<a data-toggle="collapse" href="#''' + each_rule_title + '''" aria-expanded="true">
						<i class="ti-panel"></i>
							<p>''' + each_rule_title + '''<b class="caret"></b></p>
					</a>
					<div class="collapse" id="''' + each_rule_title + '''">
						<ul class="nav">
			'''
			for each_rule_subtitle in pipeline_Dict[each_rule_title]:
				#deactive_mode
					sidebar_String += '''
							<li>
								<a href="''' + reference_Path + '''/''' + pipeline_title + '''_''' + each_rule_title + '''_''' + each_rule_subtitle + '''_''' + report_type + '''.html">
									<span class="sidebar-mini">''' + each_rule_subtitle[0].upper() + '''</span>
									<span class="sidebar-normal">''' + each_rule_subtitle.upper() + '''</span>
								</a>
							</li>
					'''
		sidebar_String += '''
						</ul>
					</div>
				</li>
		<!--#####################################################################################-->
		'''
	return sidebar_String


def build_body_html_string(body_Content, rule_title, rule_subtitle, report_type):
	#
	body_String = ''
	body_String += '''
						<div class="navbar-minimize">
							<button id="minimizeSidebar" class="btn btn-fill btn-icon"><i class="ti-more-alt"></i></button>
						</div>
						<div class="navbar-header">
							<button type="button" class="navbar-toggle">
								<span class="sr-only">Toggle navigation</span>
								<span class="icon-bar bar1"></span>
								<span class="icon-bar bar2"></span>
								<span class="icon-bar bar3"></span>
							</button>
							<a class="navbar-brand" href="#Dashboard"><small>''' + rule_title.upper() + "_" + rule_subtitle.upper() + '''</small></a>
						</div>
				<div class="collapse navbar-collapse">

					<!--form class="navbar-form navbar-left navbar-search-form" role="search">
						<div class="input-group">
							<span class="input-group-addon"><i class="fa fa-search"></i></span>
							<input type="text" value="" class="form-control" placeholder="Search...">
						</div>
					</form-->

					<ul class="nav navbar-nav navbar-right">
						<li>
							<a href="#stats" class="dropdown-toggle btn-magnify" data-toggle="dropdown">
								<i class="ti-panel"></i>
								<p>Stats</p>
							</a>
						</li>
						<li class="dropdown">
							<a href="#notifications" class="dropdown-toggle btn-rotate" data-toggle="dropdown">
								<i class="ti-bell"></i>
								<span class="notification">5</span>
								<p class="hidden-md hidden-lg">
									Notifications
									<b class="caret"></b>
								</p>
							</a>
							<ul class="dropdown-menu">
								<li><a href="#not1">Notification 1</a></li>
								<li><a href="#not2">Notification 2</a></li>
								<li><a href="#not3">Notification 3</a></li>
								<li><a href="#not4">Notification 4</a></li>
								<li><a href="#another">Another notification</a></li>
							</ul>
						</li>
						<li>
							<a href="#settings" class="btn-rotate">
								<i class="ti-settings"></i>
								<p class="hidden-md hidden-lg">
									Settings
								</p>
							</a>
						</li>
					</ul>
				</div>
			</div>
		</nav>

		<div class="content">
			<div class="container-fluid">
	'''

	if report_type == "execution_log":
		#
		body_String += """
				<div class="row">
					<div class="col-md-12">
						<div class="card">
							<div class="card-content">
								<div class="fresh-datatables">
									<table id="datatables" class="table table-striped table-no-bordered table-hover" cellspacing="0" width="100%" style="width:100%">
										<thead>
												<tr>
													<td>Task</td><td>Wildcard</td><td>Input</td><td>Output</td><td>Script</td><td>StdOut</td><td>StdErr</td><td>Elapsed Time</td>
												</tr>
										</thead>
										<tfoot>
												<tr>
													<td>Task</td><td>Wildcard</td><td>Input</td><td>Output</td><td>Script</td><td>StdOut</td><td>StdErr</td><td>Elapsed Time</td>
												</tr>
										</tfoot>
										<tbody>

		"""
		body_String += body_Content
		body_String += """
										</tbody>
									</table>
								</div>
							</div>
							<div class="card-footer">
							<hr>
							<div class="footer-title">table</div>
								<div class="pull-right">
									<button class="btn btn-info btn-fill btn-icon btn-sm">
										<i class="ti-plus"></i>
									</button>
								</div>
							</div>
						</div>
				</div>
			</div>
		"""
	elif report_type == "provenance":
		#
		body_String += """
				<div class="row">
					<div class="col-md-12">
						<div class="card">
							<div class="card-content">
		"""
		body_String += '''<img style="margin:0px auto;display:block" src="''' + body_Content + '''" alt="provenance_tag_goes_here">'''
		body_String += """
								</div>
							</div>
							<div class="card-footer">
							<hr>
							<div class="footer-title">table</div>
								<div class="pull-right">
									<button class="btn btn-info btn-fill btn-icon btn-sm">
										<i class="ti-plus"></i>
									</button>
								</div>
							</div>
						</div>
				</div>
			</div>
		"""

	return body_String


def build_main_html_string(reference_Path, file_Path, pipeline_title, sidebar_String, body_String, javascript_String):
	#
	main_html_String = ''
	main_html_String += '''
	<!doctype html>
		<html lang="en">
			<head>
				<meta charset="utf-8" />
				<link rel="apple-touch-icon" sizes="76x76" href="''' + reference_Path + '''/assets/img/apple-icon.png">
				<link rel="icon" type="image/png" sizes="96x96" href="''' + reference_Path + '''/assets/img/favicon.png">
				<meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1" />

				<title>''' + pipeline_title + '''</title>

				<meta content='width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=0' name='viewport' />
				<meta name="viewport" content="width=device-width" />

				<!-- Bootstrap core CSS     -->
				<link href="''' + reference_Path + '''/assets/css/bootstrap.min.css" rel="stylesheet" />

				<!--  Paper Dashboard core CSS    -->
				<link href="''' + reference_Path + '''/assets/css/paper-dashboard.css" rel="stylesheet"/>


				<!--  CSS for Demo Purpose, don't include it in your project     -->
				<link href="''' + reference_Path + '''/assets/css/demo.css" rel="stylesheet" />


				<!--  Fonts and icons     -->
				<link href="http://maxcdn.bootstrapcdn.com/font-awesome/latest/css/font-awesome.min.css" rel="stylesheet">
				<link href='https://fonts.googleapis.com/css?family=Muli:400,300' rel='stylesheet' type='text/css'>
				<link href="''' + reference_Path + '''/assets/css/themify-icons.css" rel="stylesheet">
			</head>
			<!--#####################################################################################-->
			<body>
				<!--#####################################################################################-->
												<!-- SIDE BAR -->
				<!--#####################################################################################-->
				<div class="wrapper">
					<div class="sidebar" data-background-color="white" data-active-color="danger">
					<!--
						Tip 1: you can change the color of the sidebar's background using: data-background-color="white | brown"
						Tip 2: you can change the color of the active button using the data-active-color="primary | info | success | warning | danger"
					-->
					<div class="logo">
						<a href="#" class="simple-text logo-mini">''' + pipeline_title[0:2].upper() + '''</a>
						<a href="#" class="simple-text logo-normal">''' + pipeline_title.upper() + '''</a>
					</div>
					<div class="sidebar-wrapper">
						<ul class="nav">
	'''
	main_html_String += sidebar_String
	main_html_String += '''
						</ul>
					</div>
				</div>
				<!--#####################################################################################-->
													<!--  END OF SIDE BAR -->
				<!--#####################################################################################-->
				<div class="main-panel">
					<nav class="navbar navbar-default">
							<div class="container-fluid">
	'''
	main_html_String += body_String
	main_html_String += '''
						</div>
					</div>
					<footer class="footer">
						<div class="container-fluid">
							<nav class="pull-left">
								<!--ul>
									<li>
										<a href="http://www.creative-tim.com">
											Creative Tim
										</a>
									</li>
									<li>
										<a href="http://blog.creative-tim.com">
											Blog
										</a>
									</li>
									<li>
										<a href="http://www.creative-tim.com/license">
											Licenses
										</a>
									</li>
								</ul-->
							</nav>
							<!--div class="copyright pull-right">
								&copy; <script>document.write(new Date().getFullYear())</script>, made with <i class="fa fa-heart heart"></i> by <a href="http://www.creative-tim.com">Creative Tim</a>
							</div-->
						</div>
					</footer>
					</div>
				</div>
			</body>
			<!--#####################################################################################-->
											<!--  END OF BODY -->
			<!--#####################################################################################-->
			<!--#####################################################################################-->
											<!--  JAVASCRIPT -->
			<!--#####################################################################################-->
			<!--   Core JS Files. Extra: TouchPunch for touch library inside jquery-ui.min.js   -->
			<script src="''' + reference_Path + '''/assets/js/jquery-3.1.1.min.js" type="text/javascript"></script>
			<script src="''' + reference_Path + '''/assets/js/jquery-ui.min.js" type="text/javascript"></script>
			<script src="''' + reference_Path + '''/assets/js/perfect-scrollbar.min.js" type="text/javascript"></script>
			<script src="''' + reference_Path + '''/assets/js/bootstrap.min.js" type="text/javascript"></script>

			<!--  Forms Validations Plugin -->
			<script src="''' + reference_Path + '''/assets/js/jquery.validate.min.js"></script>

			<!-- Promise Library for SweetAlert2 working on IE -->
			<script src="''' + reference_Path + '''/assets/js/es6-promise-auto.min.js"></script>

			<!--  Plugin for Date Time Picker and Full Calendar Plugin-->
			<script src="''' + reference_Path + '''/assets/js/moment.min.js"></script>

			<!--  Date Time Picker Plugin is included in this js file -->
			<script src="''' + reference_Path + '''/assets/js/bootstrap-datetimepicker.js"></script>

			<!--  Select Picker Plugin -->
			<script src="''' + reference_Path + '''/assets/js/bootstrap-selectpicker.js"></script>

			<!--  Switch and Tags Input Plugins -->
			<script src="''' + reference_Path + '''/assets/js/bootstrap-switch-tags.js"></script>

			<!-- Circle Percentage-chart -->
			<script src="''' + reference_Path + '''/assets/js/jquery.easypiechart.min.js"></script>

			<!--  Charts Plugin -->
			<script src="''' + reference_Path + '''/assets/js/chartist.min.js"></script>

			<!--  Notifications Plugin    -->
			<script src="''' + reference_Path + '''/assets/js/bootstrap-notify.js"></script>

			<!-- Sweet Alert 2 plugin -->
			<script src="''' + reference_Path + '''/assets/js/sweetalert2.js"></script>

			<!-- Vector Map plugin -->
			<script src="''' + reference_Path + '''/assets/js/jquery-jvectormap.js"></script>

			<!--  Google Maps Plugin    -->
			<script src="https://maps.googleapis.com/maps/api/js?key=YOUR_KEY_HERE"></script>

			<!-- Wizard Plugin    -->
			<script src="''' + reference_Path + '''/assets/js/jquery.bootstrap.wizard.min.js"></script>

			<!--  Bootstrap Table Plugin    -->
			<script src="''' + reference_Path + '''/assets/js/bootstrap-table.js"></script>

			<!--  Plugin for DataTables.net  -->
			<script src="''' + reference_Path + '''/assets/js/jquery.datatables.js"></script>

			<!--  Full Calendar Plugin    -->
			<script src="''' + reference_Path + '''/assets/js/fullcalendar.min.js"></script>

			<!-- Paper Dashboard PRO Core javascript and methods for Demo purpose -->
			<script src="''' + reference_Path + '''/assets/js/paper-dashboard.js"></script>

			<!-- Paper Dashboard PRO DEMO methods, don't include it in your project! -->
			<script src="''' + reference_Path + '''/assets/js/demo.js"></script>
			<!-- PLOTLY -->
			<script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>

			<script type="text/javascript">
				$(document).ready(function() {

				$('#datatables').DataTable({
				"pagingType": "full_numbers",
				"lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
				responsive: true,
				language: {
				search: "_INPUT_",
				searchPlaceholder: "Search records",
				}
				});


				var table = $('#datatables').DataTable();
				// Edit record
				table.on( 'click', '.edit', function () {
				$tr = $(this).closest('tr');

				var data = table.row($tr).data();
				alert( 'You press on Row: ' + data[0] + ' ' + data[1] + ' ' + data[2] + '\\'s row.' );
				} );

				// Delete a record
				table.on( 'click', '.remove', function (e) {
				$tr = $(this).closest('tr');
				table.row($tr).remove().draw();
				e.preventDefault();
				} );

				//Like record
				table.on( 'click', '.like', function () {
				alert('You clicked on Like button');
				});

				});
			</script>
	
	'''
	main_html_String += javascript_String
	main_html_String += '''
			<!--#####################################################################################-->
											<!--  END OF JAVASCRIPT -->
			<!--#####################################################################################-->
	</html>
	'''
	return main_html_String
# ################################### CONFIGURATION ############################## #


localrules: Execution_Report
configfile: "atac_seq.json"
PROCESSORS = 10
MEMORY = 4000

config_general_Dict = config["GENERAL"]
TITLE = config_general_Dict["TITLE"]
GENOME = config_general_Dict["GENOME"]
WORKDIR = fix_path(config_general_Dict["WORKDIR"])
DATADIR = fix_path(config_general_Dict["DATADIR"])
EXECDIR = fix_path(config_general_Dict["EXECDIR"])
TEMPDIR = fix_path(config_general_Dict["TEMPDIR"])

METADATA = config["METADATA"]
metadata_Dict = build_metadata_dict(METADATA)

config_pre_processing_Dict = config["PRE_PROCESSING"]
config_alignment_Dict = config["ALIGNMENT"]
config_post_alignment_Dict = config["POST_ALIGNMENT"]
config_reference_Dict = config["REFERENCE"][GENOME]

config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]

config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CHIP_Seq_py2 = config_conda_Dict["CHIP_Seq_py2"]
CHIP_Seq_py3 = config_conda_Dict["CHIP_Seq_py3"]
ACTIVATE_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_PY3 = config_conda_Dict["ACTIVATE_PY3"]


report_Dict = collections.OrderedDict()
report_Dict = config["REPORT"]
# ################################### WILDCARDS ################################## #


design_Dict = {}
PRE_PROCESSING_List = []
ALIGNMENT_List = []
POST_ALIGNMENT_List = []
NARROW_PEAK_CALLING_List = []
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
		#
		if sample_Dict["Type"] == "CASE":
			#
			design_Dict[sample_Dict["Design"]]["Case"].append(sample)
		elif sample_Dict["Type"] == "CONTROL":
			#
			design_Dict[sample_Dict["Design"]]["Control"].append(sample)
	#
	PRE_PROCESSING_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	POST_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PEAK_CALLING/MACS2/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))


POOLED_CASE_List = []
POOLED_NARROW_PEAK_CALLING_List = []
OVERLAPPED_CASE_List = []
OVERLAPPED_NARROW_PEAK_CALLING_List = []

for design in design_Dict:
	#
	POOLED_CASE_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	POOLED_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PEAK_CALLING/MACS2/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	OVERLAPPED_CASE_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{overlapped_case}.processed.bam".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	OVERLAPPED_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PEAK_CALLING/MACS2/{overlapped_case}.narrowPeak.gz".format(design=design, overlapped_case="_OVERLAPPED_".join(design_Dict[design]["Case"])))

REPORTDIR = WORKDIR + "/" + TITLE + "/" + GENOME + "/Execution_Report"
# ################################### RULES ###################################### #


rule Execution_Report:
	input: PRE_PROCESSING_List + ALIGNMENT_List + POST_ALIGNMENT_List + POOLED_CASE_List + OVERLAPPED_CASE_List + NARROW_PEAK_CALLING_List + POOLED_NARROW_PEAK_CALLING_List + OVERLAPPED_NARROW_PEAK_CALLING_List
	output:
		REPORTDIR + "/provenance_overview_diagram.svg",
	run:
		shell("""
			module load graphviz
			snakemake --snakefile atac_seq.py --configfile atac_seq.json --rulegraph | dot -Tsvg > {output[0]}
			unzip /data/RTB/datashare/Amir/DROPBOX/assets.zip -d {REPORTDIR}
			rm -rf {REPORTDIR}/__MACOSX
			""")
		for each_rule_title in report_Dict:
			for each_rule_subtitle in report_Dict[each_rule_title]:
				each_subtitle_String = ""
				for each_file in glob.glob(WORKDIR + "/" + TITLE + "/" + GENOME + "/A/REPORT/" + "/*." + each_rule_title + "." + each_rule_subtitle + ".report"):
					with open(each_file, 'r') as myfile:

						each_subtitle_String += myfile.read()
				sidebar_String = build_sidebar_html_string(report_Dict, "./", REPORTDIR, "CHIP_Seq", each_rule_title, each_rule_subtitle, "execution_log")
				if each_rule_subtitle == "provenance":
					#
					body_String = build_body_html_string("./provenance_overview_diagram.svg", each_rule_title, each_rule_subtitle, "provenance")
				else:
					#
					body_String = build_body_html_string(each_subtitle_String, each_rule_title, each_rule_subtitle, "execution_log")
				javascript_String = ""
				main_html_String = build_main_html_string("./", REPORTDIR, "CHIP_Seq", sidebar_String, body_String, javascript_String)
				write_string_down(main_html_String, REPORTDIR + "/chip_seq_" + each_rule_title + "_" + each_rule_subtitle + "_" + "execution_log.html")


rule Pre_Process:
	input:
		fq1_List = get_forward_fastq,
	output:
		processed_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.discarded.fastq.gz",
	log:
		fastqc_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.preprocess.fastqc.report",
		cutatap_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.preprocess.cutadapt.report",
	priority: 99
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "PRE_PROCESSING"
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}
			declare -a fq1_List=({input.fq1_List})

			for i in "${{fq1_List[@]}}"
			do
				fq1=$i
				fq2=${{fq1/R1/R2}}
				fq1_Name=$(basename $fq1)
				fq2_Name=$(basename $fq2)

				fq1_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/.fastq*/_fastqc.html}}
				fq1_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/.fastq*/_fastqc.zip}}
				fq1_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.fastqc.stdout}}
				fq1_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.fastqc.stderr}}


				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $fq1 1> $fq1_fastqc_stdout 2> $fq1_fastqc_stderr
				#mv $fq1_fastqc_html ${{fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $fq1_fastqc_zip ${{fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#fq1_fastqc_html=${{fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#fq1_fastqc_zip=${{fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
				
				printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$fq1" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$fq1_fastqc_html" "OUTPUT2:" "$fq1_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $fq1" >> {log.fastqc_report}

				printf "<td><p><small>%s</small></p></td>" "$fq1_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$fq1_fastqc_stderr" >> {log.fastqc_report}
				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#

				fq2_fastqc_html={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/.fastq*/_fastqc.html}}
				fq2_fastqc_zip={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/.fastq*/_fastqc.zip}}
				fq2_fastqc_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/R2*/R2.fastqc.stdout}}
				fq2_fastqc_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/R2*/R2.fastqc.stderr}}


				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $fq2 1> $fq2_fastqc_stdout 2> $fq2_fastqc_stderr
				#mv $fq2_fastqc_html ${{fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $fq2_fastqc_zip ${{fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#fq2_fastqc_html=${{fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#fq2_fastqc_zip=${{fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
				printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "$fq2" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$fq2_fastqc_html" "OUTPUT2:" "$fq2_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $fq2" >> {log.fastqc_report}

				printf "<td><p><small>%s</small></p></td>" "$fq2_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$fq2_fastqc_stderr" >> {log.fastqc_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#


				trimmed_fq1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.processed.fastq.gz}}
				trimmed_fq2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/R2*/R2.processed.fastq.gz}}
				discarded_fq1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.discarded.fastq.gz}}
				discarded_fq2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq2_Name/R2*/R2.discarded.fastq.gz}}
				cutadapt_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.cutadapt.stdout}}
				cutadapt_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/${{fq1_Name/R1*/R1.cutadapt.stderr}}

				##
				##
				start_time="$(date -u +%s)"
				cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output=$discarded_fq1 --too-short-paired-output=$discarded_fq2 --output=$trimmed_fq1 --paired-output=$trimmed_fq2 $fq1 $fq2 1> $cutadapt_stdout 2> $cutadapt_stderr
				cat $discarded_fq1 >> {output.discarded_R1}
				cat $discarded_fq2 >> {output.discarded_R2}
				cat $trimmed_fq1 >> {output.processed_R1}
				cat $trimmed_fq2 >> {output.processed_R2}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.cutatap_report}
				printf "<td>%s</td>" "CUTADAPT" >> {log.cutatap_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.cutatap_report}
				
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "INPUT1:" "$fq1" "INPUT2:" "$fq2" >> {log.cutatap_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$trimmed_fq1" "OUTPUT2:" "$trimmed_fq2" "OUTPUT3:" "$discarded_fq1" "OUTPUT4:" "$discarded_fq2" >> {log.cutatap_report}
				
				printf "<td><code>%s</code></td>" "cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output=$discarded_fq1 --too-short-paired-output=$discarded_fq2 --output=$trimmed_fq1 --paired-output=$trimmed_fq2 $fq1 $fq2" >> {log.cutatap_report}
				
				printf "<td><p><small>%s</small></p></td>" "$cutadapt_stdout" >> {log.cutatap_report}
				printf "<td><p><small>%s</small></p></td>" "$cutadapt_stderr" >> {log.cutatap_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.cutatap_report}
				printf "</tr>\\n" >> {log.cutatap_report}
				#
				#


				trimmed_fq1_fastqc_html=${{trimmed_fq1/.fastq*/_fastqc.html}}
				trimmed_fq1_fastqc_zip=${{trimmed_fq1/.fastq*/_fastqc.zip}}
				trimmed_fq1_fastqc_stdout=${{trimmed_fq1/R1*/R1.processed.fastqc.stdout}}
				trimmed_fq1_fastqc_stderr=${{trimmed_fq1/R1*/R1.processed.fastqc.stderr}}

				
				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $trimmed_fq1 1> $trimmed_fq1_fastqc_stdout 2> $trimmed_fq1_fastqc_stderr
				#mv $trimmed_fq1_fastqc_html ${{trimmed_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $trimmed_fq1_fastqc_zip ${{trimmed_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#trimmed_fq1_fastqc_html=${{trimmed_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#trimmed_fq1_fastqc_zip=${{trimmed_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
				printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$trimmed_fq1" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$trimmed_fq1_fastqc_html" "OUTPUT2:" "$trimmed_fq1_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $trimmed_fq1"  >> {log.fastqc_report}
				
				printf "<td><p><small>%s</small></p></td>" "$trimmed_fq1_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$trimmed_fq1_fastqc_stderr" >> {log.fastqc_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#


				trimmed_fq2_fastqc_html=${{trimmed_fq2/.fastq*/_fastqc.html}}
				trimmed_fq2_fastqc_zip=${{trimmed_fq2/.fastq*/_fastqc.zip}}
				trimmed_fq2_fastqc_stdout=${{trimmed_fq2/R2*/R2.processed.fastqc.stdout}}
				trimmed_fq2_fastqc_stderr=${{trimmed_fq2/R2*/R2.processed.fastqc.stderr}}

				
				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $trimmed_fq2 1> $trimmed_fq2_fastqc_stdout 2> $trimmed_fq2_fastqc_stderr
				#mv $trimmed_fq2_fastqc_html ${{trimmed_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $trimmed_fq2_fastqc_zip ${{trimmed_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#trimmed_fq2_fastqc_html=${{trimmed_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#trimmed_fq2_fastqc_zip=${{trimmed_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
				printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$trimmed_fq2" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$trimmed_fq2_fastqc_html" "OUTPUT2:" "$trimmed_fq2_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $trimmed_fq2"  >> {log.fastqc_report}
				
				printf "<td><p><small>%s</small></p></td>" "$trimmed_fq2_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$trimmed_fq2_fastqc_stderr" >> {log.fastqc_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#

				discarded_fq1_fastqc_html=${{discarded_fq1/.fastq*/_fastqc.html}}
				discarded_fq1_fastqc_zip=${{discarded_fq1/.fastq*/_fastqc.zip}}
				discarded_fq1_fastqc_stdout=${{discarded_fq1/R1*/R1.discarded.fastqc.stdout}}
				discarded_fq1_fastqc_stderr=${{discarded_fq1/R1*/R1.discarded.fastqc.stderr}}

				
				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq1 1> $discarded_fq1_fastqc_stdout 2> $discarded_fq1_fastqc_stderr
				#mv $discarded_fq1_fastqc_html ${{discarded_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $discarded_fq1_fastqc_zip ${{discarded_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#discarded_fq1_fastqc_html=${{discarded_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
				#discarded_fq1_fastqc_zip=${{discarded_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
				printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$discarded_fq1" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq1_fastqc_html" "OUTPUT2:" "$discarded_fq1_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq1"  >> {log.fastqc_report}
				
				printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stderr" >> {log.fastqc_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#


				discarded_fq2_fastqc_html=${{discarded_fq2/.fastq*/_fastqc.html}}
				discarded_fq2_fastqc_zip=${{discarded_fq2/.fastq*/_fastqc.zip}}
				discarded_fq2_fastqc_stdout=${{discarded_fq2/R2*/R2.discarded.fastqc.stdout}}
				discarded_fq2_fastqc_stderr=${{discarded_fq2/R2*/R2.discarded.fastqc.stderr}}

				
				##
				##
				start_time="$(date -u +%s)"
				fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq2 1> $discarded_fq2_fastqc_stdout 2> $discarded_fq2_fastqc_stderr
				#mv $discarded_fq2_fastqc_html ${{discarded_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#mv $discarded_fq2_fastqc_zip ${{discarded_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				#discarded_fq2_fastqc_html=${{discarded_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
				#discarded_fq2_fastqc_zip=${{discarded_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
				end_time="$(date -u +%s)"
				elapsed_time=$(($end_time-$start_time))
				##
				##

				#
				#
				printf "<tr>\\n" >> {log.fastqc_report}
				printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
				printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
				printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$discarded_fq2" >> {log.fastqc_report}
				printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq2_fastqc_html" "OUTPUT2:" "$discarded_fq2_fastqc_zip" >> {log.fastqc_report}
				printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq2"  >> {log.fastqc_report}
				
				printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stdout" >> {log.fastqc_report}
				printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stderr" >> {log.fastqc_report}

				printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
				printf "</tr>\\n" >> {log.fastqc_report}
				#
				#
			done

			processed_fq1={output.processed_R1}
			
			processed_fq1_fastqc_html=${{processed_fq1/.fastq*/_fastqc.html}}
			processed_fq1_fastqc_zip=${{processed_fq1/.fastq*/_fastqc.zip}}
			processed_fq1_fastqc_stdout=${{processed_fq1/R1*/R1.processed.fastqc.stdout}}
			processed_fq1_fastqc_stderr=${{processed_fq1/R1*/R1.processed.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $processed_fq1 1> $processed_fq1_fastqc_stdout 2> $processed_fq1_fastqc_stderr
			#mv $processed_fq1_fastqc_html ${{processed_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
			#mv $processed_fq1_fastqc_zip ${{processed_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			#processed_fq1_fastqc_html=${{processed_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
			#processed_fq1_fastqc_zip=${{processed_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$processed_fq1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$processed_fq1_fastqc_html" "OUTPUT2:" "$processed_fq1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $processed_fq1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$processed_fq1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$processed_fq1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			processed_fq2={output.processed_R2}

			processed_fq2_fastqc_html=${{processed_fq2/.fastq*/_fastqc.html}}
			processed_fq2_fastqc_zip=${{processed_fq2/.fastq*/_fastqc.zip}}
			processed_fq2_fastqc_stdout=${{processed_fq2/R2*/R2.processed.fastqc.stdout}}
			processed_fq2_fastqc_stderr=${{processed_fq2/R2*/R2.processed.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $processed_fq2 1> $processed_fq2_fastqc_stdout 2> $processed_fq2_fastqc_stderr
			#mv $processed_fq2_fastqc_html ${{processed_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
			#mv $processed_fq2_fastqc_zip ${{processed_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			#processed_fq2_fastqc_html=${{processed_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
			#processed_fq2_fastqc_zip=${{processed_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$processed_fq2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$processed_fq2_fastqc_html" "OUTPUT2:" "$processed_fq2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $processed_fq2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$processed_fq2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$processed_fq2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			discarded_fq1={output.discarded_R1}
			
			discarded_fq1_fastqc_html=${{discarded_fq1/.fastq*/_fastqc.html}}
			discarded_fq1_fastqc_zip=${{discarded_fq1/.fastq*/_fastqc.zip}}
			discarded_fq1_fastqc_stdout=${{discarded_fq1/R1*/R1.discarded.fastqc.stdout}}
			discarded_fq1_fastqc_stderr=${{discarded_fq1/R1*/R1.discarded.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq1 1> $discarded_fq1_fastqc_stdout 2> $discarded_fq1_fastqc_stderr
			#mv $discarded_fq1_fastqc_html ${{discarded_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
			#mv $discarded_fq1_fastqc_zip ${{discarded_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			#discarded_fq1_fastqc_html=${{discarded_fq1_fastqc_html/_fastqc.html/.fastqc.html}}
			#discarded_fq1_fastqc_zip=${{discarded_fq1_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$discarded_fq1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq1_fastqc_html" "OUTPUT2:" "$discarded_fq1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$discarded_fq1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			discarded_fq2={output.discarded_R2}
			
			discarded_fq2_fastqc_html=${{discarded_fq2/.fastq*/_fastqc.html}}
			discarded_fq2_fastqc_zip=${{discarded_fq2/.fastq*/_fastqc.zip}}
			discarded_fq2_fastqc_stdout=${{discarded_fq2/R2*/R2.discarded.fastqc.stdout}}
			discarded_fq2_fastqc_stderr=${{discarded_fq2/R2*/R2.discarded.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq2 1> $discarded_fq2_fastqc_stdout 2> $discarded_fq2_fastqc_stderr
			#mv $discarded_fq2_fastqc_html ${{discarded_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
			#mv $discarded_fq2_fastqc_zip ${{discarded_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			#discarded_fq2_fastqc_html=${{discarded_fq2_fastqc_html/_fastqc.html/.fastqc.html}}
			#discarded_fq2_fastqc_zip=${{discarded_fq2_fastqc_zip/_fastqc.zip/.fastqc.zip}}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$discarded_fq2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$discarded_fq2_fastqc_html" "OUTPUT2:" "$discarded_fq2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $discarded_fq2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$discarded_fq2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

		""")


rule Alignment:
	input:
		processed_R1 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PRE_PROCESSING/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz"
	output:
		alignment_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		alignment_Report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment_report.txt",
	log:
		#BOWTIE2
		bowtie2_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.report",
		bowtie2_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stdout",
		bowtie2_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.bowtie2.stderr",
		#FASTQC
		fastqc_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.report",
		fastqc_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.stdout",
		fastqc_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.alignment.fastqc.stderr",

	priority: 98
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "ALIGNMENT"
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}

			mapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R1.processed.mapped.fastq.gz
			mapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R2.processed.mapped.fastq.gz
			unmapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R1.processed.unmapped.fastq.gz
			unmapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R2.processed.unmapped.fastq.gz
			
			##
			##
			start_time="$(date -u +%s)"
			bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -1 {input.processed_R1} -2 {input.processed_R2} \
			--un-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R%.processed.unmapped.fastq.gz \
			--al-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R%.processed.mapped.fastq.gz \
			2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/ALIGNMENT/{wildcards.sample} -o {output.alignment_Bam} - 1> {log.bowtie2_stdout} 2> {log.bowtie2_stderr}
			samtools index -@ {threads} -b {output.alignment_Bam} >> {log.bowtie2_stdout} 2>> {log.bowtie2_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.bowtie2_report}
			printf "<td>%s</td>" "BOWTIE2" >> {log.bowtie2_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.bowtie2_report}
			
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "INPUT1:" "{input.processed_R1}" "INPUT2:" "{input.processed_R2}" >> {log.bowtie2_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.alignment_Bam}" "OUTPUT2:" "{output.alignment_Report}" "OUTPUT3:" "$mapped_R1" "OUTPUT4:" "$mapped_R2" "OUTPUT5:" "$unmapped_R1" "OUTPUT6:" "$unmapped_R2" >> {log.bowtie2_report}
			
			printf "<td><code>%s</code><br>" "bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 {input.processed_R1} -2 {input.processed_R2} --un-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R%.processed.unmapped.fastq.gz --al-conc-gz {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R%.processed.mapped.fastq.gz 2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/ALIGNMENT/BOWTIE2/{wildcards.sample} -o {output.alignment_Bam} -" >> {log.bowtie2_report}
			printf "<code>%s</code></td>" "samtools index -@ {threads} -b {output.alignment_Bam}" >> {log.bowtie2_report}
			
			printf "<td><p><small>%s</small></p></td>" "{log.bowtie2_stdout}" >> {log.bowtie2_report}
			printf "<td><p><small>%s</small></p></td>" "{log.bowtie2_stderr}" >> {log.bowtie2_report}

			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.bowtie2_report}
			printf "</tr>\\n" >> {log.bowtie2_report}
			#
			#

			mapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R1.processed.mapped.fastq.gz
			mapped_R1_fastqc_html=${{mapped_R1/.fastq*/_fastqc.html}}
			mapped_R1_fastqc_zip=${{mapped_R1/.fastq*/_fastqc.zip}}
			mapped_R1_fastqc_stdout=${{mapped_R1/R1*/R1.processed.mapped.fastqc.stdout}}
			mapped_R1_fastqc_stderr=${{mapped_R1/R1*/R1.processed.mapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $mapped_R1 1> $mapped_R1_fastqc_stdout 2> $mapped_R1_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$mapped_R1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$mapped_R1_fastqc_html" "OUTPUT2:" "$mapped_R1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $mapped_R1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$mapped_R1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$mapped_R1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			mapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R2.processed.mapped.fastq.gz
			mapped_R2_fastqc_html=${{mapped_R2/.fastq*/_fastqc.html}}
			mapped_R2_fastqc_zip=${{mapped_R2/.fastq*/_fastqc.zip}}
			mapped_R2_fastqc_stdout=${{mapped_R2/R2*/R2.processed.mapped.fastqc.stdout}}
			mapped_R2_fastqc_stderr=${{mapped_R2/R2*/R2.processed.mapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $mapped_R2 1> $mapped_R2_fastqc_stdout 2> $mapped_R2_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$mapped_R2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$mapped_R2_fastqc_html" "OUTPUT2:" "$mapped_R2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $mapped_R2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$mapped_R2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$mapped_R2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#


			unmapped_R1={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R1.processed.unmapped.fastq.gz
			unmapped_R1_fastqc_html=${{unmapped_R1/.fastq*/_fastqc.html}}
			unmapped_R1_fastqc_zip=${{unmapped_R1/.fastq*/_fastqc.zip}}
			unmapped_R1_fastqc_stdout=${{unmapped_R1/R1*/R1.processed.unmapped.fastqc.stdout}}
			unmapped_R1_fastqc_stderr=${{unmapped_R1/R1*/R1.processed.unmappedd.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $unmapped_R1 1> $unmapped_R1_fastqc_stdout 2> $unmapped_R1_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$unmapped_R1" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$unmapped_R1_fastqc_html" "OUTPUT2:" "$unmapped_R1_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $unmapped_R1" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$unmapped_R1_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$unmapped_R1_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#

			unmapped_R2={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/{wildcards.sample}_R2.processed.unmapped.fastq.gz
			unmapped_R2_fastqc_html=${{unmapped_R2/.fastq*/_fastqc.html}}
			unmapped_R2_fastqc_zip=${{unmapped_R2/.fastq*/_fastqc.zip}}
			unmapped_R2_fastqc_stdout=${{unmapped_R2/R2*/R2.processed.unmapped.fastqc.stdout}}
			unmapped_R2_fastqc_stderr=${{unmapped_R2/R2*/R2.processed.unmapped.fastqc.stderr}}

			##
			##
			start_time="$(date -u +%s)"
			fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $unmapped_R2 1> $unmapped_R2_fastqc_stdout 2> $unmapped_R2_fastqc_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" >> {log.fastqc_report}
			printf "<td>%s</td>" "FASTQC" >> {log.fastqc_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.fastqc_report}
			
			printf "<td><p><small>%s %s</small></p></td>" "INPUT:" "$unmapped_R2" >> {log.fastqc_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "$unmapped_R2_fastqc_html" "OUTPUT2:" "$unmapped_R2_fastqc_zip" >> {log.fastqc_report}
			printf "<td><code>%s</code></td>" "fastqc -o {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PRE_PROCESSING/ -f fastq --threads {threads} $unmapped_R2" >> {log.fastqc_report}

			printf "<td><p><small>%s</small></p></td>" "$unmapped_R2_fastqc_stdout" >> {log.fastqc_report}
			printf "<td><p><small>%s</small></p></td>" "$unmapped_R2_fastqc_stderr" >> {log.fastqc_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.fastqc_report}
			printf "</tr>\\n" >> {log.fastqc_report}
			#
			#
			""")


rule Post_Alignment:
	input:
		alignment_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
	output:
		#PICARD
		marked_duplicate_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.marked_duplicate.bam",
		picard_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.picard_report.txt",
		#QUALIMAP
		qualimap_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}_qualimap/qualimapReport.html",
		#SAMTOOLS
		processed_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		duplicate_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.duplicate.bam",
		discarded_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.discarded.bam",

	log:
		#PICARD
		picard_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.report",
		picard_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.stdout",
		picard_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.picard.stderr",
		#QUALIMAP
		qualimap_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.report",
		qualimap_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.stdout",
		qualimap_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.qualimap.stderr",
		#SAMTOOLS
		samtools_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.report",
		samtools_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.stdout",
		samtools_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.post_alignment.samtools.stderr",
	priority: 97
	threads: PROCESSORS
	message: "POST_ALIGNMENT"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}
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
			qualimap bamqc -bam {output.marked_duplicate_Bam} -nt {threads} -outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/ALIGNMENT/{wildcards.sample}_qualimap {config_post_alignment_Dict[QUALIMAP]} 1> {log.qualimap_stdout} 2> {log.qualimap_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.qualimap_report}
			printf "<td>%s</td>" "QUALIMAP" >> {log.qualimap_report}
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.qualimap_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{output.marked_duplicate_Bam}" >> {log.qualimap_report}
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.qualimap_report}" >> {log.qualimap_report}
			printf "<td><p><code>%s</code></p></td>" "qualimap bamqc -bam {output.marked_duplicate_Bam} -nt {threads} -outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/ALIGNMENT/{wildcards.sample}_qualimap  {config_post_alignment_Dict[QUALIMAP]}" >> {log.qualimap_report}

			printf "<td><p><small>%s</small></p></td>" "{log.picard_stdout}" >> {log.qualimap_report}
			printf "<td><p><small>%s</small></p></td>" "{log.picard_stderr}" >> {log.qualimap_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.qualimap_report}
			#
			#

			##
			##
			AWK_COMMAND='{{if({config_post_alignment_Dict[FILTER]}){{print $0}}}}'
			start_time="$(date -u +%s)"

			samtools view --threads {threads} -h {output.marked_duplicate_Bam} -b | bedtools intersect -v -abam stdin -b {config_reference_Dict[BLACK_LIST]} | awk '$AWK_COMMAND' |  samtools view --threads {threads} -F 1804 -f 2 -q 30 -b - -o {output.processed_Bam} 1> {log.samtools_stdout} 2> {log.samtools_stderr}
			samtools view --threads {threads} -h {output.marked_duplicate_Bam} -b | bedtools intersect -v -abam stdin -b {config_reference_Dict[BLACK_LIST]} | awk '$AWK_COMMAND' |  samtools view --threads {threads} -F 780 -f 1026 -b - -o {output.duplicate_Bam} 1>> {log.samtools_stdout} 2>> {log.samtools_stderr}
			samtools view --threads {threads} -h {output.marked_duplicate_Bam} | samtools view --threads {threads} -F 2 -b - -o {output.discarded_Bam} 1>> {log.samtools_stdout} 2>> {log.samtools_stderr}
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
			printf "<td>%s</td>" "{wildcards.sample}" >> {log.samtools_report}

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{output.marked_duplicate_Bam}" >> {log.samtools_report}
			printf "<td><p><small>%s %s</small></p><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.processed_Bam}" "OUTPUT2:" "{output.duplicate_Bam}" >> {log.samtools_report}
			printf "<td><p><code>%s</code></p><p><code>%s</code></p>" "samtools view --threads {threads} -h {output.marked_duplicate_Bam} -b | bedtools intersect -v -abam stdin -b {config_reference_Dict[BLACK_LIST]} | awk '$AWK_COMMAND' | samtools view --threads {threads} -F 1804 -f 2 -q 30 -b - -o {output.processed_Bam}" "samtools index -@ {threads} -b {output.processed_Bam}" >> {log.samtools_report}
			printf "<p><code>%s</code></p><p><code>%s</code></p>" "samtools view --threads {threads} -h {output.marked_duplicate_Bam} -b | bedtools intersect -v -abam stdin -b {config_reference_Dict[BLACK_LIST]} | awk '$AWK_COMMAND' | samtools view --threads {threads} -F 780 -f 1026 -b - -o {output.duplicate_Bam}" "samtools index -@ {threads} -b {output.duplicate_Bam}" >> {log.samtools_report}
			printf "<p><code>%s</code></p><p><code>%s</code></p></td>" "samtools view --threads {threads} -h {output.marked_duplicate_Bam} | samtools view --threads {threads} -F 2 -b - -o {output.discarded_Bam}" "samtools index -@ {threads} -b {output.discarded_Bam}" >> {log.samtools_report}

			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stdout}" >> {log.samtools_report}
			printf "<td><p><small>%s</small></p></td>" "{log.samtools_stderr}" >> {log.samtools_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.samtools_report}
			#
			#


		""")


rule Pooling_Replicates:
	input:
		case_bam_List = get_case_bam
	output:
		pooled_bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{pooled_case, .*_POOLED_.*}.processed.bam",
	log:
		pooling_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.report",
		pooling_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stdout",
		pooling_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{pooled_case, .*_POOLED_.*}.post_alignment.pooling.stderr",
	threads: PROCESSORS
	message: "POOLING"
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}

			##
			##
			start_time="$(date -u +%s)"
			samtools merge --threads {threads} {output.pooled_bam} {input.case_bam_List} 1> {log.pooling_stdout} 2> {log.pooling_stderr}
			samtools index -@ {threads} {output.pooled_bam} 1> {log.pooling_stdout} 2> {log.pooling_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##


			#
			#
			printf "<tr>\\n" > {log.pooling_report}
			printf "<td>%s</td>" "POOLING" >> {log.pooling_report}
			printf "<td>%s</td>" "{wildcards.pooled_case}" >> {log.pooling_report}

			printf "<td>" >> {log.pooling_report}
			declare -a case_bam_List=({input.case_bam_List})

			for case_index in "${{!case_bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{case_bam_List[$case_index]}}" >> {log.pooling_report}
			done

			printf "</td>" >> {log.pooling_report}
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.pooled_bam}" >> {log.pooling_report}

			printf "<td><p><code>%s</code></p><p><code>%s</code></p></td>" "samtools merge --threads {threads} {output.pooled_bam} {input.case_bam_List}" "samtools index -@ {threads} {output.pooled_bam}" >> {log.pooling_report}

			printf "<td><p><small>%s</small></p></td>" "{log.pooling_stdout}" >> {log.pooling_report}
			printf "<td><p><small>%s</small></p></td>" "{log.pooling_stderr}" >> {log.pooling_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.pooling_report}
			#
			#

			""")


rule Overlapping_Replicates:
	input:
		case_bam_List = get_case_bam
	output:
		overlapped_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{overlapped_case, .*_OVERLAPPED_.*}.processed.bam",
	log:
		overlapping_report = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/REPORT/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.report",
		overlapping_stdout = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.stdout",
		overlapping_stderr = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{overlapped_case, .*_OVERLAPPED_.*}.post_alignment.overlapping.stderr",
	threads: PROCESSORS
	message: "OVERLAP"
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}
			##
			##
			start_time="$(date -u +%s)"
			declare -a case_bam_List=({input.case_bam_List})
			group_A=""
			group_B=""
			for case_index in "${{!case_bam_List[@]}}"
			do
				if [ "$case_index" -eq 0 ] ; then
					group_A+="${{case_bam_List[$case_index]}} "
				else
					group_B+="${{case_bam_List[$case_index]}} "
				fi
			done
			bedtools intersect -wa -a $group_A -b $group_B -f 0.9 -r > {output.overlapped_Bam}  2> {log.overlapping_stderr}
			samtools index -@ {threads} {output.overlapped_Bam} 1>> {log.overlapping_stdout} 2>> {log.overlapping_stderr}
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > {log.overlapping_report}
			printf "<td>%s</td>" "OVERLAP" >> {log.overlapping_report}
			printf "<td>%s</td>" "{wildcards.overlapped_case}" >> {log.overlapping_report}

			printf "<td>" >> {log.overlapping_report}
			declare -a case_bam_List=({input.case_bam_List})

			for case_index in "${{!case_bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "<p><small>%s %s</small></p>" "INPUT${{index}}:" "${{case_bam_List[$case_index]}}" >> {log.overlapping_report}
			done

			printf "</td>" >> {log.overlapping_report}
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{output.overlapped_Bam}" >> {log.overlapping_report}

			printf "<td><p><code>%s</code></p><p><code>%s</code></p></td>" "bedtools intersect -wa -a $group_A -b $group_B -f 0.9 -r > {output.overlapped_Bam}" "samtools index -@ {threads} {output.overlapped_Bam}" >> {log.overlapping_report}

			printf "<td><p><small>%s</small></p></td>" "{log.overlapping_stdout}" >> {log.overlapping_report}
			printf "<td><p><small>%s</small></p></td>" "{log.overlapping_stderr}" >> {log.overlapping_report}
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> {log.overlapping_report}
			#
			#
			""")


rule PEAK_CALLING:
	input:
		processed_Bam = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/ALIGNMENT/{sample, ((?!.*_VS_.*|\\.).)*}.processed.bam",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/" + GENOME + "/{design}/PEAK_CALLING/MACS2/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_PY2}
			sample_Name=$(basename {input.processed_Bam})
			sample_Name=${{sample_Name%.processed.bam}}
			peak_calling_narrow_report={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/REPORT/${{sample_Name}}.peak_calling.narrow_peak.report
			peak_calling_narrow_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.peak_calling.narrow_peak.stdout
			peak_calling_narrow_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.peak_calling.narrow_peak.stderr

			##
			##
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.processed_Bam} -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/ 1> $peak_calling_narrow_stdout 2> $peak_calling_narrow_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > $peak_calling_narrow_report
			printf "<td>%s</td>" "NARROW PEAK CALLING" >> $peak_calling_narrow_report
			printf "<td>%s</td>" "$sample_Name" >> $peak_calling_narrow_report

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{input.processed_Bam}" >> $peak_calling_narrow_report
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak" >> $peak_calling_narrow_report
			printf "<td><p><code>%s</code></p></td>" "macs2 callpeak --treatment {input.processed_Bam} -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/" >> $peak_calling_narrow_report

			printf "<td><p><small>%s</small></p></td>" "$peak_calling_narrow_stdout" >> $peak_calling_narrow_report
			printf "<td><p><small>%s</small></p></td>" "$peak_calling_narrow_stderr" >> $peak_calling_narrow_report
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> $peak_calling_narrow_report
			#
			#

			black_list_filtering_report={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/REPORT/${{sample_Name}}.post_peak_calling.black_list_filtering.report
			black_list_filtering_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.post_peak_calling.black_list_filtering.stdout
			black_list_filtering_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.post_peak_calling.black_list_filtering.stderr

			##
			##
			start_time="$(date -u +%s)"
			AWK_COMMAND1='BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}'
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk '$AWK_COMMAND1' | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz 2> $black_list_filtering_stderr
			AWK_COMMAND2='BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz ) -b <(zcat -f {config_reference_Dict[BLACK_LIST]} )  | grep -P 'chr[\\dXY]+[\\t]' | awk '$AWK_COMMAND2' | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz 2> $black_list_filtering_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > $black_list_filtering_report
			printf "<td>%s</td>" "BLACK LIST FILTERING" >> $black_list_filtering_report
			printf "<td>%s</td>" "$sample_Name" >> $black_list_filtering_report

			printf "<td><p><small>%s %s</small></p>" "INPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak" >> $black_list_filtering_report
			printf "<p><small>%s %s</small></p></td>" "INPUT2:" "{config_reference_Dict[BLACK_LIST]}" >> $black_list_filtering_report
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz" >> $black_list_filtering_report
			printf "<td><p><code>%s</code></p><p><code>%s</code></p></td>" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk '$AWK_COMMAND1' | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz" "bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz ) -b <(zcat -f {config_reference_Dict[BLACK_LIST]} )  | grep -P 'chr[\\dXY]+[\\t]' | awk '$AWK_COMMAND2' | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz" >> $black_list_filtering_report
			
			printf "<td><p><small>%s</small></p></td>" "$black_list_filtering_stdout" >> $black_list_filtering_report
			printf "<td><p><small>%s</small></p></td>" "$black_list_filtering_stderr" >> $black_list_filtering_report
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> $black_list_filtering_report
			#
			#

			pvalue_filtering_report={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/REPORT/${{sample_Name}}.post_peak_calling.pvalue_filtering.report
			pvalue_filtering_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.post_peak_calling.pvalue_filtering.stdout
			pvalue_filtering_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.post_peak_calling.pvalue_filtering.stderr

			##
			##
			start_time="$(date -u +%s)"
			AWK_COMMAND3='BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}'
			zcat {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz | awk '$AWK_COMMAND3' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz 2> $pvalue_filtering_stderr
			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > $pvalue_filtering_report
			printf "<td>%s</td>" "PVALUE FILTERING" >> $pvalue_filtering_report
			printf "<td>%s</td>" "$sample_Name" >> $pvalue_filtering_report

			printf "<td><p><small>%s %s</small></p></td>" "INPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz" >> $pvalue_filtering_report
			printf "<td><p><small>%s %s</small></p></td>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz" >> $pvalue_filtering_report
			printf "<td><p><code>%s</code></p></td>" "zcat {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz | awk '$AWK_COMMAND3' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz" >> $pvalue_filtering_report
			
			printf "<td><p><small>%s</small></p></td>" "$pvalue_filtering_stdout" >> $pvalue_filtering_report
			printf "<td><p><small>%s</small></p></td>" "$pvalue_filtering_stderr" >> $pvalue_filtering_report
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> $pvalue_filtering_report
			#
			#

			signal_processing_report={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/REPORT/${{sample_Name}}.peak_calling.signal_p-value.report
			signal_processing_stdout={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.peak_calling.signal_p-value.stdout
			signal_processing_stderr={WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.peak_calling.signal_p-value.stderr

			##
			##
			start_time="$(date -u +%s)"
			raw_peaks=$(wc -l {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk '{{printf "%f", $1/1000000}}')
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name \
			--outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/ -m ppois -S ${{raw_peaks}} 1> $signal_processing_stdout 2> $signal_processing_stderr

			slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph 1>> $signal_processing_stdout 2>> $signal_processing_stderr
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.sorted.bedgraph 2>> $signal_processing_stderr
			mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph.gz 1>> $signal_processing_stdout 2>> $signal_processing_stderr
			bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bigwig 1>> $signal_processing_stdout 2>> $signal_processing_stderr

			end_time="$(date -u +%s)"
			elapsed_time=$(($end_time-$start_time))
			##
			##

			#
			#
			printf "<tr>\\n" > $signal_processing_report
			printf "<td>%s</td>" "SIGNAL PROCESSING" >> $signal_processing_report
			printf "<td>%s</td>" "$sample_Name" >> $signal_processing_report

			printf "<td><p><small>%s %s</small></p>" "INPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg" >> $signal_processing_report
			printf "<p><small>%s %s</small></p></td>" "INPUT2:" "{config_reference_Dict[CHROM_SIZE]}" >> $signal_processing_report
			printf "<td><p><small>%s %s</small></p>" "OUTPUT1:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz" >> $signal_processing_report
			printf "<p><small>%s %s</small></p>" "OUTPUT2:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_ppois.bdg" >> $signal_processing_report
			printf "<p><small>%s %s</small></p>" "OUTPUT3:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph.gz" >> $signal_processing_report
			printf "<p><small>%s %s</small></td></p>" "OUTPUT4:" "{WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bigwig" >> $signal_processing_report
			printf "<td><p><code>%s</code></p>" "raw_peaks=$(wc -l {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk '{{printf "%f", $1/1000000}}')" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name --outdir {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/ -m ppois -S ${{raw_peaks}}" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "slopBed -i {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.sorted.bedgraph" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "mv {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.sorted.bedgraph {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "bgzip -c {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph.gz" >> $signal_processing_report
			printf "<p><code>%s</code></p>" "tabix -f -p bed {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph.gz" >> $signal_processing_report
			printf "<p><code>%s</code></td></p>" "bedGraphToBigWig {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bedgraph {config_reference_Dict[CHROM_SIZE]} {WORKDIR}/{TITLE}/{GENOME}/{wildcards.design}/PEAK_CALLING/MACS2/${{sample_Name}}.bigwig" >> $signal_processing_report

			printf "<td><p><small>%s</small></p></td>" "$signal_processing_stdout" >> $signal_processing_report
			printf "<td><p><small>%s</small></p></td>" "$signal_processing_stderr" >> $signal_processing_report
			printf "<td><p><small>%s Seconds</small></p></td>" "$elapsed_time" >> $signal_processing_report
			#
			#
			""")


