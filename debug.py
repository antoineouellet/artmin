# ------------------------------------------------------------------------------------------------
# ARTmin: Mineral classification
# version May 2018
# ------------------------------------------------------------------------------------------------

# INPUT FILES #
analyses_file = 'artmin1154.xlsx'
purephases_database_file = 'DBartmin_v11.xlsx'
solidsolutions_database_file = 'DBsolidsolutions_v11.xlsx'
mixed_analyses_file = ''									
guesses_file = ''


# FORMAT OF INPUT #
col_id_analyses = 0					# Column for unique id, analyses (0-indexed) (int)
col_chem_analyses = 3				# First column for chemical analysis (int)
col_id_purephases = 0				# Column for unique id, pure phases (0-indexed) (int)
col_generic_pp = 1					# Column for generic name, pure phases (0-indexed) (int)
col_specific_pp = 0					# Column for specific name, pure phases (0-indexed) (int)
col_chem_purephases = 8				# First column for chemical composition, pure phases (int)
col_id_solisolutions = None			# Column for unique id, solid solutions (0-indexed) (int or 
										# None, def = None). If None, excel sheet names are used
col_generic_ss = None 				# Column for generic name (0-indexed) (int or None). If None
										# excel sheet names are used.
col_specific_ss = 0					# Column for generic name, solid solutions (0-indexed) (int)
col_chem_solidsolutions = 8			# First column for chemical composition, solid solutions (int)


# MODIFY DATABASE #
remove_pp = []						# Pure phases names to be removed from database (list of str)
remove_ss = []						# Solid solutions to be removed from database (list of str)
add_pp = []							# Add phases to database (chemical formula) (list of str)


# CLASSIFICATION SCHEME #
rigorous_ss = False				# Compute rigorously the solid solutions (time consuming)? 
									# (bool, def = False) 
mixed_analyses_run = False		# Do a run for mixed analyses? (bool, def = False)
omit_traces = 0					# Value under which traces are omitted (float, [0,1], def = 0)
major_elements = 1				# Value over which presence of element is considered essential 
									# (float, [0,1], def = 1)
success_cutoff = 0.1			# Normalized distance cutoff for success classification 
									# (float, [0,1], def = 0.1)
distinct_cutoff = 0.02 			# Normalized distance cutoff for unicity of classification 
									# (float, [0,1], def = 0.02)
show_only_equivalent = False	# Print out only phases within distinct_cutoff distance from best 
									# (bool, def = False). False always yields three best matches

# Relevant elements to show when present in analysis, with cut-off values (molar) from which they 
# are shown. (Dict = {key: val}) (key format is str) (val format is float, [0,1])
anomal_elem = {'Cr':0.01, 'V':0.01, 'La':0.01, 'U':0.01}


# OUTPUTS #
export_results = True						# Export results? (bool)
format_output = 'csv'						# Format of export, xlsx or csv (str)
print_progress = True						# Print progress to command window? (bool)
print_all_comments = True					# Print all comments to command window? (bool)
print_log = False							# Print .txt logfile of command window output? (bool)


# ------------------------------------------ ATTENTION -------------------------------------------
# Ne pas descendre plus bas si vous ne vous appelez pas Alexandre ou Antoine ou si vous ne passez
# pas vos longues soirées solitaires devant un écran noir plein de lignes de codes avec un
# Mountain Dew.
# ------------------------------------------------------------------------------------------------















































# ------------------------------------------------------------------------------------------------
# Bonjour à toi apprenti-Sith,

# Bienvenue du côté obscur de la Force. Ton nouveau maître: Stack Overflow.

# Hors des libraires Python disponibles en téléchargement libre (numpy, scipy, pandas, math, re, 
# collections, fastnumbers, sklearn, time et datetime), les méthodes supplémentaires proviennent 
# des modules maison artmin.py et artmin_database.py. hdbscan nécessite un interprétateur C++
# (sous Windows, Visual Studio).

# Bon chance.
# ------------------------------------------------------------------------------------------------


# Library imports
import artmin as am 
import numpy as np
import pandas as pd
import time
import datetime

t_i = time.clock()

# Name files for exports
start_time = datetime.datetime.now()
index_point = analyses_file.find('.')
name = analyses_file[:index_point]
if export_results:
	filename_output = 'classified_' + str(name) + '.' + str(format_output)
if print_log:
	str_start_time = str(start_time)
	str_start_time = int(''.join(c for c in str_start_time if c.isdigit()))
	str_start_time = str(str_start_time)
	str_start_time = str_start_time[:12]
	logfilename = 'log_' + str(str_start_time) + '.txt'


# Initialize output logger to .txt if 
if print_log:
	class Tee:
		def write(self, *args, **kwargs):
			self.out1.write(*args, **kwargs)
			self.out2.write(*args, **kwargs)
		def flush(self):
			self.out1.flush()
			self.out2.flush()
		def __init__(self, out1, out2):
			self.out1 = out1
			self.out2 = out2
	import sys
	sys.stdout = Tee(open(logfilename, 'w'), sys.stdout)


# Print initializing comments
if print_all_comments:
	print('ARTmin: Mineral classification')
	print('version May 2018')
	print('_____________________________________________________________________________________\n')
	print('Classification run on:'.ljust(32), start_time,'\n')
	
	print('Analyses to be classified:'.ljust(32), analyses_file)
	print('Pure phases database used:'.ljust(32), purephases_database_file)
	print('Solid solutions database used:'.ljust(32), solidsolutions_database_file)
	if remove_pp:
		print('Pure phases removed from db:'.ljust(32), remove_pp)
	if remove_ss:
		print('Solid solutions removed from db:'.ljust(32), remove_ss)
	if add_pp:
		print('Phases added to db:'.ljust(32), add_pp)
	print('')

	print('Type of classification: '.ljust(32), end=' ')
	if mixed_analyses_run:
		print('Mixed analyses separation')
		print('Mixed analyses separation is not implemented yet.')
	else:
		print('Regular')
	if rigorous_ss:
		print('Treatment of solid solutions:'.ljust(32), 'Rigorous')
	else:
		print('Treatment of solid solutions:'.ljust(32), 'Simplified')
	print('Classify according to:'.ljust(32), end = ' ')
	print('Non-equivalent phases\n') if show_only_equivalent else \
	print('Three best matches')
	print('Unicity of solution at:'.ljust(32), str(100*distinct_cutoff) + '%')
	print('Success of classification at:'.ljust(32), str(100*success_cutoff) + '%' , '\n')
	if omit_traces > 0:
		print('Traces omitted at: '.ljust(33) + str(100*omit_traces) + '%')
	if major_elements < 1:
		print('Elements are essential from: '.ljust(33) +str(100*major_elements) + '%')

	if export_results:
		print('Results exported to:'.ljust(32), filename_output)
	else:
		print('Results not exported')
	if print_log:
		print('Log report exported to:'.ljust(32), logfilename)
	else:
		print('No log report printed')
	print('_____________________________________________________________________________________\n')


# ------------------------------------------------------------------------------------------------
# Load data from files to ARTmin
ANALYSES = am.load_analyses(analyses_file, print_comm=print_progress)
DB_pp = am.load_pp(purephases_database_file, print_comm=print_progress)
DB_ss = am.load_ss(solidsolutions_database_file, print_comm=print_progress)

# Prepare data from raw to usable
(analyses, values_analyses, id_analyses, db_pp, values_pp, id_pp, generic_pp, specific_pp, \
db_ss, values_ss, id_ss, generic_ss, specific_ss) = am.prepare_data(ANALYSES, DB_pp, DB_ss, \
col_ida=col_id_analyses, col_chema=col_chem_analyses, col_idpp=col_id_purephases, \
col_gen_pp=col_generic_pp, col_spec_pp=col_specific_pp, col_chempp=col_chem_purephases, \
col_idss=col_id_solisolutions, col_gen_ss=col_generic_ss, col_spec_ss=col_specific_ss, \
col_chemss=col_chem_solidsolutions, omit_traces=omit_traces, print_comm=print_progress)
# ------------------------------------------------------------------------------------------------


# Print comments on imported data
if print_all_comments:
	print('Number of analyses:'.ljust(32), len(values_analyses))
	print('Number of pure phases in db:'.ljust(32), len(values_pp))
	print('Number of solid solutions in db:'.ljust(32), len(values_ss),'\n')


# ------------------------------------------------------------------------------------------------
# Compute distances from all analyses to all database entries
(DIST, SS_OPTI) = am.anal_to_all(values_analyses, values_pp, values_ss, id_analyses=id_analyses,
id_pp=id_pp, id_ss=id_ss, rigorous=rigorous_ss, print_comm=print_progress)

# Classify all computed distances to retrieve most likely phases for each analysis
classified = am.classify1(DIST, SS_OPTI, generic_pp, specific_pp, generic_ss, specific_ss,
success_cutoff=success_cutoff, distinct_cutoff=distinct_cutoff, print_comm=print_progress)
# ------------------------------------------------------------------------------------------------


# Print final comments, diagnostics and simplified report
if print_all_comments:

	# Values for success rate
	n = len(classified)
	uca = [None]*n
	aca = [None]*n
	nca = [None]*n
	j = 0
	for i in classified.index:
		uca[j] = classified['success'][i] and classified['unique'][i]
		aca[j] = classified['success'][i] and not classified['unique'][i]
		nca[j] = not classified['success'][i]
		j = j+1
	uca = np.sum(uca)
	aca = np.sum(aca)
	nca = np.sum(nca)

	uca_r = '{0:.1f}'.format(100*uca/n) + '%'
	aca_r = '{0:.1f}'.format(100*aca/n) + '%'
	nca_r = '{0:.1f}'.format(100*nca/n) + '%'
	uca_str = str(uca)
	aca_str = str(aca)
	nca_str = str(nca)
	n = str(n)
	l = len(n)

	finish_time = datetime.datetime.now()
	print('_____________________________________________________________________________________\n')
	print('Classification finished on:'.ljust(32), finish_time, '\n')
	print('Unique analyses:'.ljust(33) + uca_str.rjust(l) + '/' + n + ' = ' + uca_r.rjust(5))
	print('Ambiguous analyses:'.ljust(33) + aca_str.rjust(l) + '/' + n + ' = ' + aca_r.rjust(5))
	print('Failed classification:'.ljust(33) + nca_str.rjust(l) + '/' + n + ' = ' + nca_r.rjust(5))

	# Values for simplified report
	pd.set_option('display.width', 1000)
	pd.set_option('display.max_rows', 10000)
	# Group, sum and sort values according to generic mineral name and succes of classification
	classif = classified.groupby(['success', 'unique', 'generic1']).count()
	ind_set = set(classif.index.get_level_values(2))
	ind_list = list(ind_set)
	
	generic_report = pd.DataFrame(columns=['Classified', 'C(%)', ' Ambiguous', 'A(%)', \
	'    Failed', 'F(%)', '     Total', 'T(%)'], index=ind_list)
	if uca != 0:
		bin_uca = classif['d1(%)'].loc[True][True]
		generic_report['Classified'] = bin_uca
	if aca != 0:
		bin_aca = classif['d1(%)'].loc[True][False]
		generic_report[' Ambiguous'] = bin_aca
	if nca != 0:
		bin_nca1 = classif['d1(%)'].loc[False][True]
		bin_nca2 = classif['d1(%)'].loc[False][False]
		bin_nca = bin_nca1.add(bin_nca2, fill_value=0)
		generic_report['    Failed'] = bin_nca
	del generic_report.index.name
	generic_report = generic_report.fillna(value = 0)
	generic_report = pd.DataFrame(generic_report, dtype=int)
	generic_report = generic_report.sort_values('Classified', ascending=False)


	generic_report['     Total'] = generic_report.sum(axis=1)
	tot = generic_report.sum(axis = 0)
	generic_report['C(%)'] = 100*generic_report['Classified']/tot['Classified']
	generic_report['C(%)'] = generic_report['C(%)'].round(decimals=1)
	generic_report['A(%)'] = 100*generic_report[' Ambiguous']/tot[' Ambiguous']
	generic_report['A(%)'] = generic_report['A(%)'].round(decimals=1)
	generic_report['F(%)'] = 100*generic_report['    Failed']/tot['    Failed']
	generic_report['F(%)'] = generic_report['F(%)'].round(decimals=1)
	generic_report['T(%)'] = 100*generic_report['     Total']/tot['     Total']
	generic_report['T(%)'] = generic_report['T(%)'].round(decimals=1)
	generic_report = generic_report.fillna(value = 0)
	print('\n',generic_report,'\n')

	if export_results:
		print('Detailed classification in:'.ljust(32), filename_output)
	
	t_f = time.clock()
	elaps = t_f - t_i
	str_elaps = am.sec_to_hms(elaps)
	print('Classification completed in:'.ljust(32), str_elaps, '\n')








col = ANALYSES.columns.values[col_chem_analyses:]
ind = analyses[:, col_id_analyses]
df = pd.DataFrame(values_analyses, columns=col, index=ind)

from sklearn import datasets
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

clusterer = am.cluster_analyses(ANALYSES, col_ida=col_id_analyses, col_chema=col_chem_analyses)

classes = list(generic_report.index)[0:10]
print(classes)
class_colors = sns.color_palette('hls', len(classes))
print(class_colors)
i = 0
colors_of_points = [None]*len(classified)
for x in classified['generic1']:
	if x in classes:
		colors_of_points[i] = class_colors[classes.index(x)]
	else:
		colors_of_points[i] = (0.5, 0.5, 0.5)
	i = i + 1
#cluster_member_colors = [sns.desaturate(x,p) for x, p in zip(cluster_colors, clusterer.probabilities_)]

pca = PCA(n_components=2)
X_r = pca.fit(values_analyses).transform(values_analyses)
import matplotlib.patches as mpatches

plt.scatter(X_r[:,0], X_r[:,1], c=colors_of_points, alpha=0.5)

#legend constr
top25classes = classes
recs = []
for i in range(0,10):
	recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colors[i]))
plt.legend(recs,top25classes)
plt.show()

thefile = open('clusters.txt', 'w')
for clust in clusterer.labels_:
	thefile.write('%s\n' % clust)