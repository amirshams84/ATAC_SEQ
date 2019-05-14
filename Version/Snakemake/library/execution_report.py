import os
import sys
import pandas
# ################################### UTILITY FUNCTIONS ########################## #


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
													<td>Task</td><td>Design</td><td>Wildcard</td><td>Input</td><td>Output</td><td>Script</td><td>StdOut</td><td>StdErr</td><td>Elapsed Time</td>
												</tr>
										</thead>
										<tfoot>
												<tr>
													<td>Task</td><td>Design</td><td>Wildcard</td><td>Input</td><td>Output</td><td>Script</td><td>StdOut</td><td>StdErr</td><td>Elapsed Time</td>
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
