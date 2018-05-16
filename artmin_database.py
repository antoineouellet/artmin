# artmin_database.py module by Antoine R. Ouellet

# A collection of functions to be used in the creation of databases for ARTmin
# a mineral identification routine as developped by IOS Services 
# Geoscientifiques.
# First version (Mar. 2018)

# Throughout library and databases Z and Ln refers to indifferentiated REEs and lanthanides.

import numpy as np
import pandas as pd
import re
from fastnumbers import fast_real
from collections import Counter


# Create ordered list of all elements
elements = 'H	He	Li	Be	B	C	N	O	F	Ne	Na	Mg	Al	Si	P	S	Cl	Ar	K	Ca	Sc	Ti\
V	Cr	Mn	Fe	Co	Ni 	Cu	Zn	Ga	Ge	As	Se	Br	Kr	Rb	Sr	Y	Zr	Nb	Mo	Tc	Ru	Rh	Pd	Ag\
Cd	In	Sn	Sb	Te	I	Xe	Cs	Ba	La	Ce	Pr 	Nd	Pm	Sm	Eu	Gd	Tb	Dy	Ho	Er	Tm	Yb	Lu	Hf\
Ta	W	Re	Os	Ir	Pt	Au	Hg	Tl	Pb	Bi	Po	At	Rn	Fr	Ra	Ac	Th 	Pa	U	Np	Pu	Am	Cm	Bk\
Cf	Es	Fm	Md	No	Lr	Rf	Db	Sg	Bh	Hs	Mt	Ds	Rg	Cn	Nh	Fl	Mc	Lv	Ts	Og  Z Ln'
elements = re.findall(r'([A-Z][a-z]*)', elements)
REEs = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',\
'Z', 'Ln']

# Create ordered list of ARTmin elements
am_elements = 'F	Na	Mg	Al	Si	P	S	Cl	K	Ca	Sc	Ti	V	Cr	Mn	Fe	Co	Ni	Cu	Zn	Ga\
Ge	As	Se	Br	Rb	Sr	Y	Zr	Nb	Mo	Ru	Rh	Pd	Ag	Cd	In	Sn	Sb	Te	I	Cs	Ba	La	Hf	Ta\
W	Re	Os	Ir	Pt	Au	Hg	Tl	Pb	Bi	Th	U'
am_elements = re.findall(r'([A-Z][a-z]*)', am_elements)


def add_coeff(formula):
#	Adds coefficient (1) next to elements that were not attached to a number (implicit 1).

#	Parameters
#	----------
#	formula: String for a chemical formula in clean empirical format.
#		Clean empirical is defined to be (only): chemical elements (as letter symbols), numbers(as 
#		int or float) and parentheses. All numbers should be placed to the right of the element or 
#		group of elements it affects.

#	Returns
#	-------
#	String formula with all elements attached to a coefficient.

	# Loop over all characters of formula and add 1 when needed
	i = 0
	while i < len(formula) - 1:
		if re.match(r'[A-Z]', formula[i]):					# if letter is capital
			if re.match(r'[A-Z]', formula[i+1]):				# if followed by capital
				formula = formula[:i + 1] + '1' + formula[i+1:]		# add 1
				i = i + 2											# move forward
			elif re.match(r'[a-z]', formula[i+1]):				# if followed by lower case
				i = i + 1											# next character
			elif re.match(r'[()]', formula[i+1]):				# if followed by parenthesis
				formula = formula[:i+1] + '1' + formula[i+1:]		# add 1
				i = i + 2											# move forward
			elif re.match(r'[\d.]', formula[i+1]):				# if followed by number
				i = i + 1											# next character
		elif re.match(r'[a-z]', formula[i]):				# if letter is lower case
			if re.match(r'[A-Z]', formula[i+1]):				# if followed by capital
				formula = formula[:i+1] + '1' + formula[i+1:]		# add 1
				i = i + 2											# move forward
			elif re.match(r'[a-z]', formula[i+1]):				# if followed by lower case			
				print('error: 2 lower cases')						# error
			elif re.match(r'[()]', formula[i+1]):				# if followed by parenthesis
				formula = formula[:i+1] + '1' + formula[i+1:]		# add 1
				i = i + 2											# move forward
			elif re.match(r'[\d.]', formula[i+1]):				# if followed by number		
					i = i + 1										# next character
		elif re.match(r'[\d.]', formula[i]):				# if character is number
			if re.match(r'[A-Z]', formula[i+1]):				# if followed by capital
				i = i + 1											# next character
			elif re.match(r'[a-z]', formula[i+1]):				# if followed by lower case
				print('error: number followed by lower case')		# error
			elif re.match(r'[()]', formula[i+1]):				# if followed by parenthesis
				i = i + 1											# next character
			elif re.match(r'[\d.]', formula[i+1]):				# if followed by number
					i = i + 1										# next character
		elif re.match(r'[()]', formula[i]):					# if character is parenthesis
			i = i + 1											# next character
		else:												# if any other type of character
			print('error: not valid character')					# error
	# Scan last character and add 1 only if letter
	if re.match(r'[A-Za-z]', formula[i]):				 	# if last character is a letter
		formula = formula[:] + '1'								# add 1

	return formula


def distribute_coeff(formula):
#	Distribute all factors to elements, free of parentheses. Elements will repeat in the output
#	if they repeated in the input. This is dealt with later in 'distribute_elements()'.

#	Parameters
#	----------
#	formula: String for a chemical formula in clean empirical format with added coefficients.
#		Clean empirical with added coefficients is defined to be (only): chemical elements (as 
#		letter symbols), numbers(as int or float) and parentheses. All numbers should be placed to 
#		the right of the element or group of elements it affects. No element should be free of 
#		coefficient (be it 1).

#	Returns
#	-------
#	String formula with all elements attached to a coefficient, and no external factors left (i.e. 
#	free of parentheses).

	fformula = formula 							# temporary working formula
	test = re.search(r'([()]+)', formula)		# test for presence of parenthesis

	# Iterations while there are still parentheses left
	while test:

		# Find all innermost parentheses without a factor
		z = re.findall(r'([(][^()]*[)])(?![\d])', fformula)			# identify useless groups of ()
		ind_zz = [m.start() for m in re.finditer(r'([(][^()]*[)])(?![\d])', fformula)]	# locate
		zz = z.copy()												# copy to build new groups
		# Loop on the different groups of () to delete its parentheses
		for i in range(0, len(z)):
			zz = re.findall(r'([(][^()]*[)])(?![\d])', fformula)	# identify next group of ()
			ind_zz = [m.start() for m in re.finditer(r'([(][^()]*[)])(?![\d])', fformula)]	# locate
			zz[0] = zz[0][1:-1]										# crop parentheses
			fformula = fformula[:ind_zz[0]] + zz[0] + fformula[ind_zz[0] + len(zz[0]) + 2:]	# update

		# Find all innermost parentheses with a factor
		x = re.findall(r'([(][^()]*[)][\d][\.]*[\d]*)', fformula)	# identify useful groups of ()
		ind_x = [m.start() for m in re.finditer(r'([(][^()]*[)][\d][\.]*[\d]*)', fformula)]	# locate
		xx=x.copy()													# copy to build new groups
		# Loop on the different groups of () to distribute coefficient
		for i in range(0, len(x)):
			y = re.findall(r'([\d.]+)', x[i])			# find coefficients inside ()
			yy = y.copy()								#copy to build new coefficient list
			# Loop on numbers within the () group, do multiplication
			for j in range(0, len(y) - 1):			
				yy = [fast_real(i) for i in yy]			# convert coefficient from str to number
				yy[j] = yy[j]*yy[-1]					# apply multiplication
				yy = [str(i) for i in yy]				# convert coefficient from number to str
				ind_yy = [m.start() for m in re.finditer(r'([\d.]+)', xx[i])]	# locate all coeff
				index = ind_yy[j]												# locate j^th coeff
				xx[i] = xx[i][:index] + yy[j] + xx[i][index+len(y[j]):] 		# change j^th coeff
			# End of changing numbers within one parenthesis group
			xx[i] = re.sub(r'([)][\d.]+)','', xx[i])	# remove outside coefficient and )
			xx[i] = re.sub(r'([(])', '', xx[i])			# remove (
			fformula = fformula.replace(x[i], xx[i])	# update formula
		# End of loop within one group
		
		test = re.search(r'([()]+)', fformula)		# test for presence of parenthesis
	
	return fformula


def parse_elements(formula):
#	Calculate number of each element for a formula unit.

#	Parameters
#	----------
#	formula: String for a chemical formula in clean empirical format.
#		Clean empirical is defined to be (only): chemical elements (as letter symbols), numbers(as
#		int or float) and parentheses. All numbers should be placed to the right of the element or 
#		group of elements it affects.

#	Returns
#	-------
#	Dictionary of chemical composition of given formula. Keys are chemical elements and values are 
#	number of atoms of given element in the formula. 
#	e.g.: Fe3Al2(SiO4)3 -> {'Fe': 3, 'Al': 2, 'Si': 3, 'O': 12} 
	
	formula = add_coeff(formula)
	formula = distribute_coeff(formula)

	# From the given formula, retrieve number associated to element
	x = re.findall(r'([A-Z][a-z]*)([\d]+[.]?[\d]*)', formula)	# list of tuples for el. and num.
	L = len(x)
	y = [None]*L												# create an empty list of lists
	y = [list(x[i]) for i in range(0, L)]						# copy values from tuples to lists
	for i in range(0, L):										# convert number strings to real
		y[i][1] = fast_real(y[i][1])							# numbers

	cnt = Counter()												# create empty counter
	for key, value in y:										# count values to element								
		cnt[key] += value
	chemical_compo = dict(cnt)									# convert to regular dictionary
	
	return chemical_compo


def nparse_elements(formula):
#	From a given chemical formula, calculate normalized molar composition

#	Parameters
#	----------
#	formula: String for a chemical formula in clean empirical format.
#		Clean empirical is defined to be (only): chemical elements (as letter symbols), numbers(as
#		int or float) and parentheses. All numbers should be placed to the right of the element or
#		group of elements it affects.

#	Returns
#	-------
#	Dictionary of chemical composition of given formula. Keys are chemical elements and values are 
#	normalized proportions of element in the formula. 
#	e.g.: Fe3Al2(SiO4)3 -> {'Fe': 0.15, 'Al': 0.1, 'Si': 0.15, 'O': 0.6} 

	chemical_compo = parse_elements(formula)

	# Normalization
	tot = sum(chemical_compo.values())
	chemical_compo = {key: chemical_compo[key]/tot for key in chemical_compo}

	return chemical_compo

# !!!!!!!ADD NPARSE FOR CUSTOM SET OF ELEMENTS, LIKE ARTMIN ELEMENTS


def read_database(input_file, mineral_name_col = 0, mineral_formula_col = 6, compo_col = 8, 
	start_line = 1):
#	Read a database from .xlsx/xls file.

#	Parameters
#	----------
#	input_file: 			.xls/.xlsx file of pure phases database

#	Optional
#	--------
#	mineral_name_col:		column (0-indexed) of file containing mineral name (default ARTmin
#							format, col = 0)
#	mineral_formula_col:	column (0-indexed) of file containing chemical formula in clean 
#							empirical (default ARTmin format, col = 6)
#	compo_col:		 		first column (0-indexed) of file from which chemical composition should
#							written/overwritten (default ARTmin format, col = 8)
#	start_line:				first line (0-indexed) of data

#	Returns
#	-------
#	Data for ARTmin calculations.

#	Notes
#	-----
#	-	Identical version in artmin.py library. Please be careful when modifying one without
#		modifying the other.

	# Read input file for database
	DB = pd.read_excel(input_file)					# complete database in pd.DataFrame
	headers = (list(DB.columns.values))				# list of database headers
	db = DB.values									# db is array version of database
	DB_infos = DB[headers[:compo_col]]				# pd.DataFrame without possible numericals
	db_infos = DB_infos.values
	DB_values = DB[headers[compo_col:]]
	db_values = DB_values.values					# array of database values
	db_values = np.float64(db_values)				# converted to float64 np.Array for calc
	
	(n_phases, n_elements) = np.shape(db_values)	# size of database
	names = db[:, mineral_name_col]					# mineral_names
	formulas = db[start_line-1:, mineral_formula_col]	# list of chemical formulas

	return DB


def create_database(input_file, mineral_formula_col = 6, compo_col = 8, start_line = 1, 
	export_type = None, export_name = 'database_allelements.xlsx'):
#	Create a numerical database from file of formulas (all chemical elements considered).

#	Parameters
#	----------
#	input_file: 			.xls/.xlsx file of pure phases database

#	Optional
#	--------
#	mineral_formula_col:	column (0-indexed) of file containing chemical formula in clean 
#							empirical (default ARTmin format, col = 6)
#	compo_col:		 		first column (0-indexed) of file from which chemical composition should
#							written/overwritten (default ARTmin format, col = 8)
#	start_line:				first line (0-indexed) of data
#	export_type:			format of export output. Excel('Excel', 'excel', 'XLS' or 'xls') or
#							CSV('CSV' or 'csv'). Default is None.
#	export_name:			string for name of output file. File extension and 'export_type' have
#							to agree. Default is 'database_allelements.xlsx'

#	Returns
#	-------
#	DataFrame of new complete database AND writes a .xlsx file IF specified.

#	Notes
#	-----
#	-	All chemical elements are considered along with 'Z'(last element) for all undifferentiated
#		REEs and 'Ln'. Otherwise elements are ordered for atomic number.
#	-	Excel export is much slower than CSV.
#	-	CSV Import/Export does not deal well with special characters.
	
	# Read input file for database
	DB = pd.read_excel(input_file)					# complete database in pd.DataFrame
	headers = (list(DB.columns.values))				# list of database headers
	DB_infos = DB[headers[:compo_col]]				# pd.DataFrame without possible numericals
	db = DB.values									# db is array version of database
	db_formulas = db[start_line - 1:, mineral_formula_col]	# list of chemical formulas
	
	# Calculate composition for each formula
	db_compo = [nparse_elements(db_formulas[i]) for i in range(0, len(db_formulas))]
	# Build numerical part of dataframe
	DF = pd.DataFrame(columns = elements, data = db_compo)		# read composition from dict
	DF = pd.DataFrame.fillna(DF, value = 0)						# replace NaN with 0
	# Assemble dataframe
	df = pd.concat([DB_infos, DF], axis = 1)
	
	# Export final database
	if export_type in ('Excel', 'excel', 'XLS', 'xls'):
		df.to_excel(export_name, index = False)					# export to excel
	elif export_type in ('CSV', 'csv'):
		df.to_csv(export_name, index = False)					# export to csv
	else:
		pass 													# do not export

	return df


def create_artmindatabase(input_file, mineral_formula_col = 6, compo_col = 8, start_line = 1, 
	export_type = None, export_name = 'database_artmin.xlsx'):
#	Create a numerical database from formulas (ARTmin chemical elements considered).

#	Parameters
#	----------
#	input_file: 			.xls/.xlsx file of pure phases database

#	Optional
#	--------
#	mineral_formula_col:	column (0-indexed) of file containing chemical formula in clean 
#							empirical (default ARTmin format, col = 6)
#	compo_col:		 		first column (0-indexed) of file from which chemical composition should
#							written/overwritten (default ARTmin format, col = 8)
#	start_line:				first line (0-indexed) of data
#	export_type:			format of export output. Excel('Excel', 'excel', 'XLS' or 'xls') or
#							CSV('CSV' or 'csv'). Default is None.
#	export_name:			string for name of output file. File extension and 'export_type' have
#							to agree. Default is 'database_artmin.xlsx'
#	overw:					Boolean. Overwrite input database, by default False. If set to True,
#							export_name is bypassed and input_file is used.

#	Returns
#	-------
#	DataFrame of new ARTmin database AND writes a .xlsx file

#	Notes
#	-----
#	-	Only 58 elements used in ARTmin are considered, all REEs grouped under 'La' (including
#		'La'-'Lu', 'Z' and 'Ln'.
#	-	Excel export is much slower than CSV.
#	-	CSV Import/Export does not deal well with special characters.
	
	# Read input file for database
	DB = pd.read_excel(input_file)					# complete database in pd.DataFrame
	headers = (list(DB.columns.values))				# list of database headers
	DB_infos = DB[headers[:compo_col]]				# pd.DataFrame without possible numericals
	db = DB.values									# db is array version of database
	db_formulas = db[start_line - 1:, mineral_formula_col]	# list of chemical formulas
	
	# Calculate composition for each formula
	db_compo = [nparse_elements(db_formulas[i]) for i in range(0, len(db_formulas))]
	# Group REEs
	for i in range(len(db_compo)):
		db_compo[i]['La'] = sum([db_compo[i].get(REEs[j], 0) for j in range(len(REEs))])
	# Build numerical part of dataframe
	DF = pd.DataFrame(columns = am_elements, data = db_compo)	# read composition from dict
	DF = pd.DataFrame.fillna(DF, value = 0)						# replace NaN with 0
	DF = DF.div(DF.sum(axis=1), axis=0)							# normalize to elements
	DF = pd.DataFrame.fillna(DF, value = 0)						# replace NaN with 0
	# Assemble dataframe
	df = pd.concat([DB_infos, DF], axis = 1)
	
	# Export final database
	if export_type in ('Excel', 'excel', 'XLS', 'xls'):
			df.to_excel(export_name, index = False)				# export to excel
	elif export_type in ('CSV', 'csv'):
		df.to_csv(export_name, index = False)					# export to csv
	else:
		pass 													# do not export
	
	return df


def create_customdatabase(input_file, db_elements, mineral_formula_col = 6, compo_col = 8,
	start_line = 1, export_type = None, export_name = 'database_customelements.xlsx', 
	groupREE = False):
#	Create a numerical database from formulas (custom chemical elements considered).

#	Parameters
#	----------
#	input_file: 			.xls/.xlsx file of pure phases database
#	db_elements:			list of elements (as string of element symbol, e.g.: Ca)

#	Optional
#	--------
#	mineral_formula_col:	column (0-indexed) of file containing chemical formula in clean 
#							empirical (default ARTmin format, col = 6)
#	compo_col:		 		first column (0-indexed) of file from which chemical composition should
#							written/overwritten (default ARTmin format, col = 8)
#	start_line:				first line (0-indexed) of data
#	export_type:			format of export output. Excel('Excel', 'excel', 'XLS' or 'xls') or
#							CSV('CSV' or 'csv'). Default is None.
#	export_name:			string for name of output file. File extension and 'export_type' have
#							to agree. Default is 'database_customelements.xlsx'
#	groupREE:				Boolean. Group all REEs under the label 'La'. Default is False

#	Returns
#	-------
#	DataFrame of new ARTmin database AND writes a .xlsx file

#	Notes
#	-----
#	-	Only 58 elements used in ARTmin are considered, all REEs grouped under 'La'.
#	-	Excel export is much slower than CSV.
#	-	CSV Import/Export does not deal well with special characters.
	
	# Read input file for database
	DB = pd.read_excel(input_file)					# complete database in pd.DataFrame
	headers = (list(DB.columns.values))				# list of database headers
	DB_infos = DB[headers[:compo_col]]				# pd.DataFrame without possible numericals
	db = DB.values									# db is array version of database
	db_formulas = db[start_line - 1:, mineral_formula_col]	# list of chemical formulas
	
	# Calculate composition for each formula
	db_compo = [nparse_elements(db_formulas[i]) for i in range(0, len(db_formulas))]
	# Group REEs
	if groupREE:
		for i in range(len(db_compo)):
			db_compo[i]['La'] = sum([db_compo[i].get(REEs[j], 0) for j in range(len(REEs))])
	# Build numerical part of dataframe
	DF = pd.DataFrame(columns = db_elements, data = db_compo)	# read composition from dict
	DF = pd.DataFrame.fillna(DF, value = 0)						# replace NaN with 0
	DF = DF.div(DF.sum(axis=1), axis=0)							# normalize to elements
	DF = pd.DataFrame.fillna(DF, value = 0)						# replace NaN with 0
	# Assemble dataframe
	df = pd.concat([DB_infos, DF], axis = 1)
	
	# Export final database
	if export_type in ('Excel', 'excel', 'XLS', 'xls'):
		df.to_excel(export_name, index = False)					# export to excel
	elif export_type in ('CSV', 'csv'):
		df.to_csv(export_name, index = False)					# export to csv
	else:
		pass 													# do not export
	
	return df


def remove_purephase(phase_name, database):
#	Remove phase from pure phases database

#	Parameters
#	----------
#	phase: 	If string or list of strings: Phase names to be removed, as they appear in database
#			If integer or list of integers: n^th line (0-indexed)
#	database: pd.DataFrame, pure phase database

#	Returns
#	-------
#	pd.DataFrame of new updates

#	Notes
#	-----
#	Phases are removed from database in the session only, not actually written to external file.
#	Name is supposed to be in first column.

	db = database.values
	names = db[:, 0]

	if type(phase_name) == list:
		n_tot = len(phase_name)
		for i in range(n_tot):
			if type(phase_name[i]) == str:
				phase_index = np.argwhere(names == phase_name[i])
				if not phase_index:
					raise ValueError('Phase name not found')
			elif type(phase_name[i]) == int:
				phase_index = phase_name[i]
			else:
				raise ValueError('Invalid type for phase name') 
			database = database.drop(phase_index)
		database = database.reset_index(drop = True)
	elif type(phase_name) == str:
		phase_index = np.argwhere(names == phase_name)
		if not phase_index:
			raise ValueError('Phase name not found')
		database = database.drop(phase_index)
		database = database.reset_index(drop = True)
	elif type(phase_name) == int:
		database = database.drop(phase_index)
		database = database.reset_index(drop = True)
	else:
		raise ValueError('Invalid type for phase name')
	
	return database


def remove_solidsolution(ssolution, ssolution_database):
#	Remove

#	Parameters
#	----------

#	Optional
#	--------

#	Returns
#	-------

#	Notes
#	-----
	
	return database


def add_purephase(phase_formula, database):
#	Add phase from pure phases database

#	Parameters
#	----------
#	phase: 	String or list of strings: Phase (as chemical formulas) to be added
#	database: pd.DataFrame, pure phase database

#	Returns
#	-------
#	pd.DataFrame of new updates

#	Notes
#	-----
#	Phases are added in database in the session only, not actually written to external file.


	return database


def add_solidsolution(phase, database):
#	Add

#	Parameters
#	----------

#	Optional
#	--------

#	Returns
#	-------

#	Notes
#	-----
	
	return database


def distancedb_id(database, id_p, d_sup = 0.0001, d_inf = 0):
#	Check if mineral 'id_p' from database is within distance d: [d_inf, d_sup) from other minerals.

#	Parameters
#	----------
#	database: pd.DataFrame pure phases database
#	id_p: integer, corresponding to 'Database ID' column in database

#	Optional
#	--------
#	d_sup: integer/float, superior value of distance interval (default 0.0001)
#	d_inf: integer/float, inferior value of distance interval (default 0)

#	Returns
#	-------
#	Single key dictionary. Key is phase from 'id_p' and value is list of minerals within
#	specified distance interval. Will include itself if d_inf = 0.

#	Notes
#	-----
#	Main interest of this function is to check database for other minerals equivalent to 'id_p'.
#	To do so, leave d_inf = 0 and set d_sup to an arbitrarily low value (such as default 0.0001).
#	d_sup may be raised for not-equivalent but dangerously close pairs (say d = 0.01).
#	For most applications d_inf should be 0.

	# Read database
	DB = database 									# complete database in pd.DataFrame
	db = DB.values									# db is array version of database
	headers = (list(DB.columns.values))				# list of database headers
	db_values = db[:,8:]							# array of database values
	db_values = np.float64(db_values)				# converted to float64 np.Array for calc
	(n_phases, n_elements) = np.shape(db_values)
	names = db[:,0]

	# Function for distance between two endmembers
	def distance1(compo1, compo2):					# Euclidian distance between two phases
		coeff = compo1 - compo2
		dist = np.linalg.norm(coeff)
		return(dist)

	if d_sup < d_inf:
		raise ValueError('d_sup should not be smaller than d_inf')

	D = np.ones(n_phases)							# Pre-allocate distance vector
	for j in range(0, n_phases):					# Loop over all database phases
		D[j] = distance1(db_values[id_p, :], db_values[j, :])# Distance between phase and j^th
	x = np.argwhere(np.logical_and(D >= d_inf, D < d_sup))		# Locate all phases whithin range	
	x = x.flatten()
	
	res_phases = [None]*len(x)						# Pre-allocate resulting phases list
	for i in range(0, len(x)):						# Create strings in resulting list
		indice = x[i]
		res_phases[i] = names[indice]

	return {names[id_p] : res_phases}


def distancedb_seqcheck(database, id_p, d_sup = 0.0001, d_inf = 0):
#	Check if database minerals from 'id_p' and forward are within distance d: [d_inf, d_sup) from 
#	other minerals. This code starts at 'id_p' and moves forward (i = i + 1) as long as resulting
#	ensemble is empty or contains only the mineral itself.

#	Parameters
#	----------
#	database: pd.DataFrame of pure phases database
#	id_p: integer, corresponding to 'Database ID' column in database

#	Optional
#	--------
#	d_sup: integer/float, superior value of distance interval (default 0.0001)
#	d_inf: integer/float, inferior value of distance interval (default 0)

#	Returns
#	-------
#	Single key dictionary. Key is phase from 'id_p' and value is list of minerals within
#	specified range. Also, number of this mineral following the id_p

#	Notes
#	-----
#	Main interest of this function is to check database for other minerals equivalent to 'id'.
#	To do so, leave d_inf = 0 and set d_sup to an arbitrarily low value (such as default 0.0001).
#	d_sup may be raised for not-equivalent but dangerously close pairs (say d = 0.01).
#	No clear interest in raising d_inf.

	# Read database
	DB = database 									# complete database in pd.DataFrame
	db = DB.values									# db is array version of database
	headers = (list(DB.columns.values))				# list of database headers
	db_values = db[:,8:]							# array of database values
	db_values = np.float64(db_values)				# converted to float64 np.Array for calc
	(n_phases, n_elements) = np.shape(db_values)
	names = db[:,0]

	# Function for distance between two endmembers
	def distance1(phase1, phase2):					# Euclidian distance between two phases
		coeff = phase1 - phase2
		dist_sq = np.sum(np.power(coeff, 2))
		dist = np.sqrt(dist_sq)
		return(dist)

	if d_sup < d_inf:
		raise ValueError('d_sup should not be smaller than d_inf')

	test = True
	criterion = 1 if d_inf == 0 else 0
	# Loop while no groups of phases has met specified interval
	while test:
		D = np.ones(n_phases)						# Pre-allocate distance vector
		for j in range(0, n_phases):				# Loop over all database phases
			D[j] = distance1(db_values[id_p, :], db_values[j, :])# Distance between phase and j^th
		x = np.argwhere(np.logical_and(D >= d_inf, D < d_sup))	# Locate all phases whithin range	
		x = x.flatten()	
		if len(x) > criterion:						# If result is not empty and more than itself							
			test = False								# stop iterating
			res_phases = [None]*len(x)					# pre-allocate list of phases
			# Build resulting list
			for k in range(0, len(x)):
				indice = x[k]						
				res_phases[k] = names[indice]
			result = {names[id_p] : res_phases}
		else:										# Keep going to next phase in database
			test = True
			id_p = id_p + 1

	return {names[id_p] : res_phases}, id_p


def distancedb_totalcheck(database, d_sup = 0.0001, d_inf = 0):
#	Check database for groups of minerals within distance d: [d_inf, d_sup) from themselves.

#	Parameters
#	----------
#	database: pd.DataFrame of pure phases database

#	Optional
#	--------
#	d_sup: integer/float, superior value of distance interval (default 0.0001)
#	d_inf: integer/float, inferior value of distance interval (default 0)

#	Returns
#	-------
#	Dictionary. Keys are phases from database and value is list of minerals within
#	specified range.

#	Notes
#	-----
#	Careful as calculation times quickly ramp up with size of database ({n=100: t=1s}, 
#	{n=1000: t=10s}, {n=5300: t=220s}, including database loading at around 1s)
#
#	Values may repeat. That is for example, Andalusite-Sillimanite-Kyanite polymorphs
#	(all chemically equivalent) will yield three different dictionary entries (for d being very 
#	small) that are:
#	'Andalusite': 	['Andalusite', 'Sillimanite', 'Kyanite']
#	'Sillimanite': 	['Andalusite', 'Sillimanite', 'Kyanite']
#	'Kyanite': 		['Andalusite', 'Sillimanite', 'Kyanite'].
#
#	This is done on purpose so that if, given phase1, phase2, phase3, and d_12 = 0.1 and
#	d_23 = 0.1, but d_13 = 0.2 and function is run for d = 0.11, complete answer set is obtained:
#	'phase1': ['phase1', 'phase2']
#	'phase2': ['phase1', 'phase2', 'phase3']
#	'phase3': ['phase2', 'phase3'],
#	and not a subset of it.

	# Read database
	DB = database  									# complete database in pd.DataFrame
	db = DB.values									# db is array version of database
	headers = (list(DB.columns.values))				# list of database headers
	db_values = db[:,8:]							# array of database values
	db_values = np.float64(db_values)				# converted to float64 np.Array for calc
	(n_phases, n_elements) = np.shape(db_values)
	names = db[:,0]

	# Function for distance between two endmembers
	def distance1(phase1, phase2):					# Euclidian distance between two phases
		coeff = phase1 - phase2
		dist_sq = np.sum(np.power(coeff, 2))
		dist = np.sqrt(dist_sq)
		return(dist)

	if d_sup < d_inf:
		raise ValueError('d_sup should not be smaller than d_inf')

	criterion = 1 if d_inf == 0 else 0
	result = {}									# Pre-allocate result dictionary
	# Loop on all phases of database
	for i in range (n_phases):
		D = np.ones(n_phases)						# Pre-allocate distance vector
		for j in range(0, n_phases):				# Loop over all database phases for distance
			D[j] = distance1(db_values[i, :], db_values[j, :])# Distance between phase and j^th
		x = np.argwhere(np.logical_and(D >= d_inf, D < d_sup))# Locate all phases whithin range	
		x = x.flatten()	
		if len(x) > criterion:						# If result is not empty and more than itself
			res_phases = [None]*len(x)					# pre-allocate list of phases
			for k in range(0, len(x)):					# build resulting list
				indice = x[k]								# where k^th phase is
				res_phases[k] = names[indice]				# name of k^th phase
			result[names[i]] = res_phases				# append dictionary with new group

	return result


def formula_check(database, formula_col = 6):
#	Check if database empirical formulas respect some basic rules

#	Parameters
#	----------
#	database: pd.DataFrame file of pure phases database

#	Optional
#	--------
#	formula_col: column of database in which formula is located (default ARTmin col = 6).

#	Returns
#	-------
#	Dictionary. 'Key' is problem and 'value' is list of problematic formulas.
#
#	Keys:
#		'num_start': 		Formula starts with a number
#		'invalid_char':		Invalid (not accepted) character in formula
#		'unmatch_par':		Parentheses (open-close) do not match
#		'unmatch_float':	Number is associated with two points

#	Notes
#	-----
#	Eventually should be coded as a regex
#	\(?([a-zA-Z]{1,2})([1-9]*)([+-]?)([1-9]*)\)?([1-9]*)

	# Read database
	DB = database									# complete database in pd.DataFrame
	db = DB.values									# db is array version of database
	headers = (list(DB.columns.values))				# list of database headers
	db_formulas = db[:,6]							# array of database values
	(n_phases, n_elements) = np.shape(db)

	result = {}

	# Check for formulas that start with a number
	result['num_start'] = []
	for i in range(n_phases):						# loop on formulas
		formula = db_formulas[i]
		prob = formula[0].isdigit()						# problem if starts with a digit
		# add to problem dictionary
		if prob:
			result['num_start'] = result['num_start'] + [formula]

	# Check if only accepted characters are present
	A = elements + ['1','2','3','4','5','6','7','8','9','0','(',')','.']
	result['invalid_char'] = []					
	for i in range(n_phases):						# loop on formulas
		prob = False
		formula = db_formulas[i]
		j = 0
		while j < len(formula) and not prob:			# loop on characters
			char = formula[j]								# single character
			dual_char = formula[j: j+2]						# dual character
			if dual_char in A:							# check if two characters are allowed
				j = j + 2
			elif char in A:								# if not, check single character
				j = j + 1
			else:										# if not, problem
				prob = True
				result['invalid_char'] = result['invalid_char'] + [formula]
				j = j + 2

	# Check if parentheses match
	result['unmatch_par'] = []
	for i in range(n_phases):						# loop on formulas
	    balanced = True
	    formula = db_formulas[i]
	    cnt = 0
	    j = 0
	    while j < len(formula) and balanced:			# loop on characters
	        char = formula[j]								
	        if char == "(":									# if opening parenthesis, + 1
	            cnt = cnt + 1
	           
	        elif char == ")":								# if closing parenthesis
	            if cnt == 0:									# unbalanced if no corresponding (
	                balanced = False
	            else:											# else, - 1
	            	cnt = cnt - 1
	        else:											# if any other character
	        	pass 											# do nothing
	        j = j + 1										# next character
	    # Result for given formula
	    if cnt > 0 or not balanced:
	    	result['unmatch_par'] = result['unmatch_par'] + [formula]

	# Check if numbers are valid float
	result['unmatch_float'] = []
	for i in range(n_phases):
		formula = db_formulas[i]
		prob_float = re.findall(r'([\d][\.]+[\d]+[\.]+)', formula)
		if prob_float:
			result['unmatch_float'] = result['unmatch_float'] + [formula]

	return result