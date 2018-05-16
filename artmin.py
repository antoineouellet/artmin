# artmin.py module by Antoine R. Ouellet

# A collection of functions to be used through ARTmin mineral identification
# routine as developped by IOS Services Geoscientifiques
# First version (Apr. 2018)

# Functions may be codependent and modifying one could potentially harm a second one in a
# non-trivial way. For further custimization of functions in this library, it is therefore
# recommended to only modify copies of functions with different names.

import numpy as np
import scipy as sp
import pandas as pd
import math
import time
import hdbscan
from scipy.optimize import minimize
from scipy.optimize import nnls
from scipy.optimize import lsq_linear


def sec_to_hms(laps):
#	Convert seconds to hour:minute:second.

#	Parameters
#	----------
#	laps: (int or float) time in seconds.

#	Returns
#	-------
#	Returns string of format 24h:60m:60s.

	h = np.floor(laps/3600)
	m = np.floor((laps - h*3600)/60)
	s = np.floor(laps - h*3600 - m*60)
	h = str(int(h))
	m = str(int(m))
	s = str(int(s))
	str_time = h + 'h:' + m + 'm:' + s + 's'

	return str_time


def load_analyses(analyses_file, sheet_name=0, print_comm=True):
#	Load analyses from Excel file.

#	Parameters
#	----------
#	analyses_file:	str, name of Excel file + extension

#	Returns
#	-------
#	DataFrame of analyses data
	
	if print_comm:
		print('Loading analytical data...', end=' ', flush=True)
	ANALYSES = pd.read_excel(analyses_file, sheet_name=sheet_name)
	if print_comm:
		print('Done.')

	return ANALYSES


def load_pp(purephases_database, sheet_name=0, print_comm=True):
#	Load pure phases database from Excel file.

#	Parameters
#	----------
#	purephases_database: str, name of Excel file + extension

#	Returns
#	-------
#	DataFrame of pure phases database data

	if print_comm:
		print('Loading pure phases database...', end=' ', flush=True)
	DB_pp = pd.read_excel(purephases_database, sheet_name=sheet_name)
	if print_comm:
		print('Done.')

	return DB_pp


def load_ss(solidsolutions_database, sheet_name=None, print_comm=True):
#	Load solid solutions database from Excel file.

#	Parameters
#	----------
#	solidsolutions_database: 	str, name of Excel file + extension
#	sheet_name: 				list (of strings), name of sheets (i.e. phase solid solutions) 
#								to import. Def = None, imports all sheets.

#	Returns
#	-------
#	Ordered Dictionary of DataFrames. Key is phase solid solution and DataFrames are chemical
#	composition of poles.

	if print_comm:
		print('Loading solid solutions database...', end=' ', flush=True)
	DB_ss = pd.read_excel(solidsolutions_database, sheet_name=sheet_name)
	if print_comm:
		print('Done.')

	return DB_ss


def prepare_data(ANALYSES, DB_pp, DB_ss, col_ida=0, col_chema=5, col_idpp=0, col_gen_pp=1,
	col_spec_pp=0, col_chempp=8, col_idss=None, col_gen_ss=None, col_spec_ss=0, col_chemss=8, 
	omit_traces=0, print_comm=True):
#	Prepare data from DataFrame imports to format usable in further calculations.

#	Parameters
#	----------
#	ANALYSES: 	DataFrame, data for all analyses.
#	DB_pp:		DataFrame, data for all pure phases.
#	DB_ss:		Ordered Dictionary of DataFrames, data for solid solutions.

#	Optional
#	--------
#	col_[...]:	int, column (0-indexed) in which data is to be found. id = unique id of analysis,
#				chem = first column of chemistry, gen = generic name, spec =  specific name,
#				a = analysis, pp = pure phases, ss = solid solutions.
#	omit_traces: cut off value under which an element's value is set to zero (analysis renormalized)

#	Returns
#	-------
#	Simple arrays and lists (or dictionary of arrays and lists) of data to be used.

	analyses = ANALYSES.values										# array version
	headers_analyses = (list(ANALYSES.columns.values))				# list of headers
	values_analyses = analyses[:, col_chema:]						# array version of values
	values_analyses = np.float64(values_analyses)					# converted to float64
	(n_analyses, n_elements) = np.shape(values_analyses)			# number of analyses and elem.
	id_analyses = list(ANALYSES.ix[:,col_ida])						# unique id of analysis

	db_pp = DB_pp.values											# array version of database
	headers_pp = (list(DB_pp.columns.values))						# list of database headers
	values_pp = db_pp[:, col_chempp:]								# array of database values
	values_pp = np.float64(values_pp)								# converted to float64
	n_pp = len(values_pp)											# number of pure phases in db
	id_pp = list(DB_pp.ix[:,col_idpp])								# id of pure phase
	generic_pp = list(DB_pp.ix[:,col_gen_pp])						# first name of pure phase
	specific_pp = list(DB_pp.ix[:,col_spec_pp])						# second name of pure phase
	for i in range(len(generic_pp)):								# if no more generic name,
		if type(generic_pp[i]) != str:									# use specific
			generic_pp[i] = specific_pp[i]

	db_ss = {k: v.values for k, v in DB_ss.items()}					# dictionary of arrays of db
	values_ss = {k: v[:,col_chemss:] for k, v in db_ss.items()}		# dictionary of array of values
	values_ss = {k: np.float64(v) for k, v in values_ss.items()}	# converted to float64
	n_ss = len(values_ss)											# number of solid solutions
	id_ss = list(DB_ss.keys())										# list of solid solutions name
	if type(col_gen_ss) == int:										# label1, if in a column
		generic_ss = {k: v[:,col_gen_ss:] for k, v in db_ss.items()}
	elif not col_gen_ss:
		generic_ss = id_ss.copy()
	else:
		raise ValueError('col_gen_ss must be integer or None')
	specific_ss = {k: v[:,col_spec_ss] for k, v in db_ss.items()}	# label2 from column

	if 0 < omit_traces < 1:
		values_analyses[values_analyses < omit_traces] = 0
		tot_a = np.sum(values_analyses, axis=1)
		tot_a = np.transpose(np.tile(tot_a, (n_elements,1)))
		values_analyses = values_analyses/tot_a
		values_pp[values_pp < omit_traces] = 0
		tot_pp = np.sum(values_pp, axis=1)
		tot_pp = np.transpose(np.tile(tot_pp, (n_elements,1)))
		values_pp = values_pp/tot_pp
	elif omit_traces == 0:
		pass
	else:
		raise ValueError('omit_traces is bound to [0,1)')

	if print_comm:
		print('Data succesfully built.\n')

	return analyses, values_analyses, id_analyses, db_pp, values_pp, id_pp, generic_pp, \
	specific_pp, db_ss, values_ss, id_ss, generic_ss, specific_ss


def euclid_distance(compo1, compo2):
#	Compute distance (euclidian) between two minerals (or two groups) of fixed composition.

#	Parameters
#	----------
#	compo1 and compo2: 1-d array of floats (or integers), or 2-d arrays of floats, with different
#		minerals in different lines.

#	Returns
#	-------
#	Eulidian distance (float) between two n-dimensional compositions. 

#	Notes
#	-----
#	-	Euclidian distance defined as d = sqrt(sum_i((x2_i-x1_i)^2)))
#	-	If compositions are properly normalized, distance d is bound in {0 <= d <= sqrt(2)}
	
	sh1 = np.shape(compo1)
	sh2 = np.shape(compo2)
	if sh1 != sh2:
		raise ValueError('Shapes of inputs must match')
	if len(sh1) == 1:
		N = 0
	elif len(sh2) == 2:
		N = 1
	else:
		raise ValueError('Unexpected shape of input')

	coeff = compo1 - compo2
	dist = np.linalg.norm(coeff, axis=N)
	
	return(dist)


def euclid_dist_to_ss(compo, solid_solution, method = 'SLSQP'):
#	Compute distance (euclidian) between a mineral of fixed composition and a solid solution 
#	mineral.

#	Parameters
#	----------
#	compo: 1-d array of floats (or integers)
#	solid_solution: 2-d array of floats (or integers). Different lines are different end-members
#		and columns different chemical elements

#	Optional
#	--------
#	method:	String or callable. Type of numerical solver, see scipy.optimize.minimize 
#	documentation. Default is 'SLSQP'.

#	Returns
#	-------
#	- Eulidian distance (float) between a given 'compo' and an optimal solid solution 
#	- Array of optimal solid solution (fractions of its constituents)

#	Notes
#	-----
#	-	Euclidian distance defined as d = sqrt(sum_i((x2_i-x1_i)^2)))
#	-	If compositions are properly normalized, distance d is bound in {0 <= d <= sqrt(2)}
#	-	Answer with corresponding distance are optimized in a least-square sense. That is,
#		resulting solid solution is the "closest" (i.e. euclidian distance) to fixed phase
#		composition 'compo'.
#	-	Constrained optimization from https://www.youtube.com/watch?v=cXHvC_FGx24

	n_poles = len(solid_solution)
	A = np.transpose(solid_solution)			# Matrix of solid solution composition
	b = np.transpose(compo)						# Composition vector of 'compo'

	# Define bounds on x
	minmax = (0.0, 1.0)	
	bound = [None]*n_poles
	bound = [minmax for i in range(n_poles)]
	# Define normalization constraint
	def constraint(x):
		return np.sum(x) - 1
	con = {'type': 'eq', 'fun': constraint}
	# Set of constraints
	#cons = [con1, con2]
	# Define function to be minimized, that is |Ax - b|
	def f(x):
		y = np.dot(A, x) - b
		return np.linalg.norm(y)
	# Initial guess, random
	x0 = sp.zeros(len(solid_solution))

	# Minimize distance function 'f' by finding optimal 'x'
	sol = minimize(f, x0, method = method, constraints = con, bounds = bound)
	dist = sol.fun
	ss_opti = sol.x

	return (dist, ss_opti)


def simpeuclid_dist_to_ss(compo, solid_solution):
#	Compute distance (euclidian) between a mineral of fixed composition and a simplified solid 
#	solution mineral.

#	Parameters
#	----------
#	compo: 1-d array of floats (or integers)
#	solid_solution: 2-d array of floats (or integers). Different lines are different end-members
#		and columns different chemical elements

#	Returns
#	-------
#	- Eulidian distance (float) between a given 'compo' and an optimal simplified solid solution 
#	- Array of optimal simplified solid solution (fractions of its constituents)

#	Notes
#	-----
#	This version omits the normalization constraint, HOWEVER it is much faster. This
#	simplification may yield aberrations, artificially low distances and unrealistic solid
#	solutions. Still, in almost all cases when it is used in a classification scheme, final
#	results are comparable so it may be used when calculation times are an issue.
	
	# Input data
	A = np.transpose(solid_solution)
	b = np.transpose(compo)

	# Non-negative linear optimization
	(ss_opti, dist) = nnls(A,b)				# Fortran wrapper for non-negative linear solver

	return (dist, ss_opti)


def anal_to_all(values_anal, values_pp, values_ss, id_analyses=[None], id_pp=[None], id_ss=[None], 
	rigorous=True, toDF=True, print_comm=True):
#	Compute the distance of an analysis or multiple analyses to all pure phases and solid solutions
#	in databases.

#	Parameters
#	----------
#	values_anal:  Array of chemical composition of analyses (different analyses on different lines)
#	values_pp:	  Array of chemical composition of pure phases
#	values_ss:    Dictionary (keys = solid solutions) of arrays of chemical composition (lines
#				  = end-members)

#	Optional
#	--------
#	rigorous:		Boolean, compute optimal solid solution rigorously (euclidian_dist_to_ss) or in
#					simplified form (simpeuclidian_dist_to_ss)
#	toDF:       	Boolean, write to DataFrame. If false, written to array, default is True
#	analyses_id:	List, id-labels for analysis, if written to DataFrame. Default is sequential
#					number, 0-indexed.
#	pp_id:			List, id-labels for pure phases, if written to DataFrame. Default is 
#					sequential number, 0-indexed.
#	ss_id:			List, id-labels for solid solutions, if written to DataFrame. Default is 
#					sequential number, 0-indexed.

#	Returns
#	-------
#	-	pd.DataFrame (or array) of distances (index = analyses, columns = solid solutions)
#	-	pd.DataFrame (or array) of optimal solid solutions constituents

	n_anal = len(values_anal)							# Number of analyses
	n_pp = len(values_pp)								# Number of pure phases in DB
	n_ss = len(values_ss)								# Number of solid solutions in DB
	dist_pp = np.zeros((n_anal, n_pp))					# pre-allocate distance array
	dist_ss = np.zeros((n_anal, n_ss))					# pre-allocate distance array
	ss_opti = np.zeros((n_anal, n_ss))					# pre-allocate optimal solid solution array
	ss_opti = np.array(ss_opti, dtype=np.ndarray)

	# Comments
	if print_comm:
		print('Calculating...', end='\n', flush=True)
		t = -1
		start = time.clock()
		if rigorous:
			est_time = n_anal*0.16
			str_est_time = sec_to_hms(est_time)
			print('[Estimated time for calculations: ' + str_est_time + ']')
		if not rigorous:
			est_time = n_anal*0.0035
			str_est_time = sec_to_hms(est_time)
			print('[Estimated time for calculations: ' + str_est_time + ']')

	# Loop over all analyses
	for i in range(n_anal):
		
		# Print progress
		if print_comm:
			t_prev = t
			t = math.floor(20*(i+1)/n_anal)
			if t != t_prev:
				print('/' + str(5*t) + '%', end='', flush=True)

		analysis = values_anal[i]
		# Loop over all solid solutions
		j = 0
		for key in values_ss:
			ssolution = values_ss[key]
			if rigorous: 
				(dist_ss[i,j], ss_opti[i,j]) = euclid_dist_to_ss(analysis, ssolution)
			else:
				(dist_ss[i,j], ss_opti[i,j]) = simpeuclid_dist_to_ss(analysis, ssolution)
			j += 1
		# Simultaneously compute distance with all pure phases
		analysis = np.tile(analysis, (n_pp, 1))			# Broadcast for parallel calculations
		coeff = analysis - values_pp
		dist_pp[i][:] = np.linalg.norm(coeff, axis = 1)	# Assign distance i to final array

	# Index/headers as numbers if not well behaved label-list
	if len(id_analyses) != n_anal:
		id_analyses = range(n_anal)
	if (len(id_pp) + len(id_ss)) != (n_pp + n_ss):
		mineral_id = range((n_pp + n_ss))
	else:
		mineral_id = id_pp
		mineral_id.extend(id_ss)

	# Prepare output as DataFrame or array depending on input
	dist = np.concatenate((dist_pp, dist_ss), axis = 1)
	if toDF:
		DIST = pd.DataFrame(dist, index = id_analyses, columns = mineral_id)
		SS_OPTI = pd.DataFrame(ss_opti, index = id_analyses, columns = id_ss)
	else:
		DIST = dist
		SS_OPTI = ss_opti

	if print_comm:
		laps = time.clock() - start
		str_laps = sec_to_hms(laps)
		print('\n' + '[Finished calculating in ' + str_laps + ']', '\n')

	return DIST, SS_OPTI


def classify1(DIST, SS_OPTI, generic_pp, specific_pp, generic_ss, specific_ss, success_cutoff=0.1,
	distinct_cutoff=0.02, print_comm=True):
#	Deprecated.
#	Classify analysis to most likely phase from total distance DataFrame.

#	Parameters
#	----------
#	DIST:			DataFrame of distances. Index = analysis and Columns = database phases

#	Optional
#	--------
#	print_comm:		Boolean, print commentaries and progression, default is True

#	Returns
#	-------
#	DataFrame indexed for analysis id. Contains maximum 3 cases of fitting minerals, each with
	
	# Print comments, estimated time
	if print_comm:
		start = time.clock()
		print('Classifying...')
		est_time = len(DIST)*0.0002
		str_est_time = sec_to_hms(est_time)
		print('[Estimated time for classifying: ' + str_est_time + ']')
	# End of comments on time

	# Values for classification
	dist = DIST.values 								# array of all distances
	dist = dist/np.sqrt(2)							# normalized array of all distances
	generic = generic_pp + generic_ss 				# list of generic names
	id_anal = (list(DIST.index.values))				# list of IDs of analyses
	ind = np.argsort(dist, axis=1)					# indices of sorted distances for a given anal.
	(n_anal, n_ps) = np.shape(dist)					# shape of the data
	n_ss = len(SS_OPTI.columns)						# number of solid solutions
	n_pp = n_ps - n_ss 								# number of pure phases

	# Columns for success and unicity of solution
	unique = [(dist[i][ind[i][1]]-dist[i][ind[i][0]]) >= distinct_cutoff for i in range(n_anal)]
	success = [dist[i][ind[i][0]] < success_cutoff for i in range(n_anal)]

	# Build the three best anwers (distance, type, )
	dist1 = ['{0:.2f}'.format(dist[i][ind[i][0]]) for i in range(n_anal)]
	type1 = [('solid solution' if ind[i][0] >= n_pp else 'pure phase') for i in range(n_anal)]
	specific1 = np.array(np.zeros(n_anal), dtype=str)
	generic1 = [generic[ind[i][0]] for i in range(n_anal)]
	for i in range(n_anal):
		if type1[i] == 'pure phase':
			specific1[i] = specific_pp[ind[i][0]]
		else:
			maxx = max(SS_OPTI[generic1[i]][id_anal[i]])
			maxx = '{0:.1f}'.format(maxx)
			imax = np.argmax(SS_OPTI[generic1[i]][id_anal[i]])
			specific1[i] = specific_ss[generic1[i]][imax] + '(' + str(maxx) + ')'

	dist2 = ['{0:.2f}'.format(dist[i][ind[i][1]]) for i in range(n_anal)]
	type2 = [('solid solution' if ind[i][1] >= n_pp else 'pure phase') for i in range(n_anal)]
	specific2 = np.array(np.zeros(n_anal), dtype=str)
	generic2 = [generic[ind[i][1]] for i in range(n_anal)]
	if generic2 != generic2:
		generic2 = specific2
	for i in range(n_anal):
		if type2[i] == 'pure phase':
			specific2[i] = specific_pp[ind[i][1]]
		else:
			maxx = max(SS_OPTI[generic2[i]][id_anal[i]])
			maxx = '{0:.1f}'.format(maxx)
			imax = np.argmax(SS_OPTI[generic2[i]][id_anal[i]])
			specific2[i] = specific_ss[generic2[i]][imax] + '(' + str(maxx) + ')'

	dist3 = ['{0:.2f}'.format(dist[i][ind[i][2]]) for i in range(n_anal)]
	type3 = [('solid solution' if ind[i][2] >= n_pp else 'pure phase') for i in range(n_anal)]
	specific3 = np.array(np.zeros(n_anal), dtype=str)
	generic3 = [generic[ind[i][2]] for i in range(n_anal)]
	if generic3 != generic3:
		generic3 = specific3	
	for i in range(n_anal):
		if type3[i] == 'pure phase':
			specific3[i] = specific_pp[ind[i][2]]
		else:
			maxx = max(SS_OPTI[generic3[i]][id_anal[i]])
			maxx = '{0:.1f}'.format(maxx)
			imax = np.argmax(SS_OPTI[generic3[i]][id_anal[i]])
			specific3[i] = specific_ss[generic3[i]][imax] + '(' + str(maxx) + ')'

	classified = pd.DataFrame(columns=['unique', 'success', 'd1(%)', 'type1', 'generic1', \
		'specific1', 'd2(%)', 'type2', 'generic2', 'specific2', 'd3(%)', 'type3', 'generic3', \
		'specific3'], index=id_anal)
	
	classified['unique'] = unique
	classified['success'] = success
	classified['d1(%)'] = dist1
	classified['type1'] = type1
	classified['generic1'] = generic1
	classified['specific1'] = specific1
	classified['d2(%)'] = dist2
	classified['type2'] = type2
	classified['generic2'] = generic2
	classified['specific2'] = specific2
	classified['d3(%)'] = dist3
	classified['type3'] = type3
	classified['generic3'] = generic3
	classified['specific3'] = specific3

	classified.to_csv('classified.csv')

	# Print elapsed time
	if print_comm:
		laps = time.clock() - start
		str_laps = sec_to_hms(laps)
		print('[Finished classifying in ' + str_laps + ']')

	return classified


def classify2(DIST, SS_OPTI, DB_pp, DB_ss, print_comm=True):
	pass





































def build_ss(values_ss, print_comm=True):
#	Compute the coefficient values for function representing the solid solution.

#	Parameters
#	----------
#	values_ss:	Array or dictionary of arrays (float). Solid solutions composition. Different poles
#				on different lines, chemical elements in columns.

#	Optional
#	--------
#	print_comm:	Boolean, print commentaries to command windown. Def = True.

#	Returns
#	-------
#	Vector array or dictionary of vector arrays (float). Each component ai is coefficient for
#	corresponding chemical element xi (as input in values_ss).

#	a1*x1 + a2*x2 + ... + an*xn - 1 = 0

	if type(values_ss) == dict:
		built_ss = dict(keys = values_ss.keys)
		for key in values_ss:
			S = values_ss[key]							# Solid solution matrix
			if type(S) != np.ndarray:
				raise ValueError('Dictionary items must be arrays')
			n_poles, n_elements = np.shape(S)
			a0 = np.ones((n_poles,1))					# Arbitrary constant
			S_inv = np.linalg.pinv(S)					# Moore-Penrose inverse
			a = np.dot(S_inv, a0)						# Resulting coefficient vector
			built_ss[key] = a
	elif type(values_ss) == np.ndarray:
		S = values_ss
		n_poles, n_elements = np.shape(S)
		a0 = np.ones((n_poles,1))					# Arbitrary constant
		S_inv = np.linalg.pinv(S)					# Moore-Penrose inverse
		a = np.dot(S_inv, a0)						# Resulting coefficient vector
		built_ss = a
	else:
		raise ValueError('input must be array or dictionary of arrays')

	return built_ss


def euclid_phase_to_ss(compo, solid_solution):
#	Deprecated
#	Compute distance (euclidian) between a mineral of fixed composition and a solid solution 
#	mineral.

#	Parameters
#	----------
#	compo: 			1-d array of floats (or integers), chemical composition x = [xi]
#	solid_solution: 1-d array of floats (or integers), coefficient of hyperplane, a = [ai(xi)],
#					representing a solid solution.

#	Returns
#	-------
#	- Eulidian distance (float) between a given 'compo' and an optimal solid solution 
#	- Array of optimal solid solution (fractions of its constituents)

#	Notes
#	-----
#	-	Euclidian distance defined as d = sqrt(sum_i((x2_i-x1_i)^2)))
	
	m1, n1 = np.shape(compo)
	m2, n2 = np.shape(solid_solution)
	print(m1, n1, m2, n2)
	print(np.dot(compo,solid_solution))

	if (m1, n1) == (m2, n2):
		compo = np.transpose(compo)
		m1, n1 = np.shape(compo)
		m2, n2 = np.shape(solid_solution)
	if (m1, n1) == (n2, m2):
		a = solid_solution
		dist = abs(np.dot(compo, a) - 1)/np.linalg.norm(a)
	else:
		raise ValueError('Could not perform dot product on arguments because of shape')

	return dist


def anal_to_purephases(values_anal, values_pp, toDF = True, analyses_id = [None], 
	pp_id = [None]):
#	Compute the distance of an analysis or multiple analyses to all pure phases in database.

#	Parameters
#	----------
#	values_anal:  Array of chemical composition of analyses (different analyses on different lines)
#	values_pp:    Array of chemical composition of pure phases

#	Optional
#	--------
#	toDF:			Boolean, write to DataFrame. If false, written to array, default is True
#	analyses_id:	List, id-labels for analysis, if written to DataFrame. Default is sequential
#					number, 0-indexed.
#	pp_id:			List, id-labels for pure phases, if written to DataFrame. Default is 
#					sequential number, 0-indexed.

#	Returns
#	-------
#	pd.DataFrame (or array) of distances (index = analyses, columns = pure phases)

	n_anal = len(values_anal)							# Number of analyses
	n_pp = len(values_pp)								# Number of phases in DB
	dist = np.zeros((n_anal, n_pp))						# pre-allocate distance array

	# Loop over all analyses
	for i in range(n_anal):
		analysis = values_anal[i]
		analysis = np.tile(analysis, (n_pp, 1))			# Broadcast for parallel calculations
		coeff = analysis - values_pp
		dist[i][:] = np.linalg.norm(coeff, axis = 1)	# Assign distance i to final array

	# Index/headers as numbers if not well behaved label-list
	if len(analyses_id) != n_anal:
		analyses_id = range(n_anal)
	if len(pp_id) != n_pp:
		pp_id = range(n_pp)
	# Prepare output as DataFrame or array depending on input
	if toDF:
		DIST = pd.DataFrame(dist, index = analyses_id, columns = pp_id)
	else:
		DIST = dist

	return DIST


def anal_to_solidsolutions(values_anal, values_ss, rigorous = True, toDF = True, analyses_id
	= [None], ss_id = [None]):
#	Compute the distance of an analysis or multiple analyses to all solid solutions in database.

#	Parameters
#	----------
#	values_anal:  Array of chemical composition of analyses (different analyses on different lines)
#	values_ss:    Dictionary (keys = solid solutions) of arrays of chemical composition (lines
#				  = end-members)

#	Optional
#	--------
#	rigorous:		Boolean, compute optimal solid solution rigorously (euclidian_dist_to_ss) or in
#					simplified form (simpeuclidian_dist_to_ss)
#	toDF:       	Boolean, write to DataFrame. If false, written to array, default is True
#	analyses_id:	List, id-labels for analysis, if written to DataFrame. Default is sequential
#					number, 0-indexed.
#	ss_id:			List, id-labels for solid solutions, if written to DataFrame. Default is 
#					sequential number, 0-indexed.

#	Returns
#	-------
#	pd.DataFrame (or array) of distances (index = analyses, columns = solid solutions)

	n_anal = len(values_anal)							# Number of analyses
	n_ss = len(values_ss)								# Number of solid solutions in DB
	dist = np.zeros((n_anal, n_ss))						# pre-allocate distance array
	sol = np.zeros((n_anal, n_ss))
	sol = np.array(sol, dtype = list)

	# Loop over all analyses
	for i in range(n_anal):
		analysis = values_anal[i]
		# Loop over all solid solutions
		j = 0
		for key in values_ss:
			ssolution = values_ss[key]
			if rigorous: 
				(dist[i,j], b) = euclid_dist_to_ss(analysis, ssolution)
			else:
				(dist[i,j], b) = simpeuclid_dist_to_ss(analysis, ssolution)
			j += 1

	# Index/headers as numbers if not well behaved label-list
	if len(analyses_id) != n_anal:
		analyses_id = range(n_anal)
	if len(ss_id) != n_ss:
		ss_id = range(n_ss)
	# Prepare output as DataFrame or array depending on input
	if toDF:
		DIST = pd.DataFrame(dist, index = analyses_id, columns = ss_id)
	else:
		DIST = dist

	return DIST


def separate_mixes(mixed_analyses, possible_componts, toDF=True):
#	Find the optimal mix of two minerals for a given analysis.

#	Parameters
#	----------
#	mixed_analyses:		2-d array (or list of lists) of analyses compositions suspected to be a mix
#						of two other minerals
#	possible_componts:	2-d array (or list of lists) of mineral compositions suspected to be the
#						components of the mixed analysis

#	Returns
#	-------
#	min_dist:			float, minimal distance from analysis yield by optimal mix
#	phases_mix:			list of strings, the two phases name of the optimal mix
#	prop:				proportion

#	Notes
#	-----
#	Very slow, should be parallelized in some way.

	n_mix = len(mixed_analyses)
	n_componts = len(possible_componts)

	boun = [(0.0, 1.0)]								# Bounds on x
	results = np.zeros((n_mix, 4))					# Pre-allocate results array
	results = np.array(results, dtype = str)		# Make it a string array
	headers_results = ['d', 'x', 'Mineral1(x)', 'Mineral2(1-x)']

	# Loop on all mixed analyses to separate
	for i in range(n_mix):
		dist = np.ones((n_componts, n_componts))		# Pre-allocate distance matrix
		X = np.zeros((n_componts, n_componts))			# Pre-allocate proportion matrix
		mix_anal = mixed_analyses[i]					# Select i^th mixed analysis
		for j in range(1, n_componts):					# Loop on first mineral of pair
			mineral1 = possible_componts[j]					# Select first mineral
			for k in range(0, j):							# Loop on second mineral of pair
				mineral2 = possible_componts[k]					# Select second mineral
				def f(x):										# Define function to optimize
					compo_mix = x*mineral1 + (1-x)*mineral2
					y = euclid_distance(mix_anal, compo_mix)
					return y
				x0 = 0.5										# Initial guess
				sol = minimize(f, x0, method = 'SLSQP', bounds = boun)
				dist[j, k] = sol.fun 							# Allocate optimal distance
				X[j, k] = sol.x[0]
		ij = np.argmin(dist)							# flattened indice of min distance 
		ij = np.unravel_index(ij, dist.shape)			# unraveled indices
		results[i][0] = '{0:.2f}'.format(np.amin(dist))
		results[i][1] = '{0:.2f}'.format(X[ij])
		results[i][2] = ij[0]
		results[i][3] = ij[1]

	# Prepare output as DataFrame or array depending on input
	if toDF:
		RES = pd.DataFrame(results, columns = headers_results)
	else:
		RES = dist

	return RES


def separate_mixes_fromguess(mixed_analyses, main_componts, secondary_componts):
	pass

def cluster_analyses(ANALYSES, col_ida=0, col_chema=5):
# http://hdbscan.readthedocs.io/en/latest/

	analyses = ANALYSES.values
	values_analyses = analyses[:, col_chema:]
	col = ANALYSES.columns.values[col_chema:]
	ind = analyses[:, col_ida]
	df = pd.DataFrame(values_analyses, columns=col, index=ind)
	n_anal = len(df)
	print(n_anal)
	
	clusterer = hdbscan.HDBSCAN(min_cluster_size=int(n_anal*0.005), min_samples=1).fit(df)
	print(clusterer.labels_)
	print(clusterer.labels_.max())
	print(clusterer)

	return clusterer

def plot_PCA_classified(classified):
	pass