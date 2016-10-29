import sys, re, math, os, numpy
# ============================================ Student info methods================================================
def get_student_name():
	# @TO_STUDENT: Write your name here
	student_name = "Faizy Ahsan"
	if not student_name:
		raise ValueError("Error: you forgot to add your name in get_student_name method.")
	return student_name

def get_student_id():
	# @TO_STUDENT: Write your student id here
	student_id = "260517720"
	if not student_id:
		raise ValueError("Error: you forgot to add your student id in get_student_id method.")
	return student_id
# =================================================================================================================

# =================================== Validate input/output methods================================================
def validate_Q3_1_input_format(subopt_result):
	if not isinstance(subopt_result, list) or [sr for sr in subopt_result if re.search("[^\(\)\.]", sr)]:
		raise ValueError("Error: your input should be list of strings (secondary structures, ie alphabet is '().')")

def validate_Q3_1_output_format(result):
	if not isinstance(result, list) or [sr for sr in result if not isinstance(sr, list)]:
		raise ValueError("Error: your output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ].")
	if [sr for sr in result if sr[0] >= sr[1]]:
		raise ValueError("Error: check your i and j values. They should satisfy: i > j.")
# =================================================================================================================
# ================================================== Helper methods================================================
def parse_subopt_result_file(filepath):
	'''
	Parsing of a standard txt file that contains result of
		"RNAsubopt -p __k__ < myFasta.fasta > subopt_result_file.txt"
		where __k__ is parameter. (Filename chosen randomly. Please, feel free to use your own names.)
	@args:
	filepath: (full or relative) path to the subopt_result_file.txt.
	(Filename chosen randomly. Please, feel free to use your own names.)

	@return: subopt_result: list of the strings (assumed to be secondary structures)
	'''
	subopt_result = []
	with open(filepath, 'r') as f:
		for i, line in enumerate(f):
			if i < 2:
				continue
			subopt_result.append(line.strip())
	return subopt_result

def parse_dot_ps_file(filepath):
	'''
	Parsing of a dot.ps file that contains result of RNAfold program
	@args:
	filepath: (full or relative) path to the dot.ps.
	@return:
	dot_ps_result: list f lists with i, j, freq_i_j
	'''
	dot_ps_result = []
	with open(filepath, 'r') as f:
		is_data = False
		for line in f:
			if not is_data and line.startswith('%start of base pair probability data'):
				is_data = True
				continue
			elif is_data and line.startswith('showpage'):
				break
			elif is_data:
				# take only first 3 numbers
				data_line = line.split()[:3]
				dot_ps_result.append(
					[int(data_line[0]), int(data_line[1]), float(data_line[2])]
				)
	return dot_ps_result

# =================================================================================================================
def get_answer_Q3_1(subopt_result):
	'''
	This method should be implemented by student.
	@args:
	subopt_result: a list of the secondary structures produced by RNAsubopt -p <k> for particular input

	@return:
	result: list of lists (example is below) with indices and relevant frequency.
	example [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16], ....... ]

	@note:
	check input/output as advised in code. Question will be marked as 0 in case of not following the formats.
	'''
	# basic check for the proper input
	validate_Q3_1_input_format(subopt_result)
	# @TO_STUDENT: Write your code here

	#subopt_result = ['((.))', '(..).']
	#print (subopt_result)
	# declare result
	result = []
	
	N = len( subopt_result )

	# store pair counts
	pc = dict()

	for seq in subopt_result:
		tmpstack = []

		for s in range(0, len(seq) ):
			currpair = ()
			if seq[s] == '(':
				tmpstack.append(s)
			elif seq[s] == '.':
				continue
			else:
				currpair = (tmpstack.pop(), s) 

				if not currpair in pc.keys():
					pc[ currpair ] = 1
				else:
					pc[ currpair ] = pc[ currpair ] + 1
	
	for k in pc.keys():
		prob = pc[k] / N

		res = list(k)
		res.append( prob )

		result.append( res )



	# @TO_STUDENT: use result variable for results. below is an example of an expected format for result.
	#result = [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16] ]

	# @TO_STUDENT: output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ]
	# use validate_Q3_output_format(result) to validate the output
	validate_Q3_1_output_format(result)
	#print ( result )
	#sys.exit()
	return result

def get_answer_Q3_2(q3_1_result, dot_ps_result):
	'''
	This method should be implemented by student.
	Compare output from RNAfold and result of question3_1 for the same sequence and return an error (see text assignment)
	result_error is expected to be numeric
	'''
	result_error = 0
	# @TO_STUDENT: Write your code here (trust me, answer is not 0 :-) )

	#print ( q3_1_result)
	#print ( dot_ps_result )

	#dot_ps_result = [ [1,2,0.5], [2,4,0.9]]
	#q3_1_result = [ [0,1,0.25], [1,3,0.81] ]

	for rfold in dot_ps_result:
		curre = 0.0
		#print ( rfold )

		for q3 in q3_1_result:
			
			if (q3[0]+1) == rfold[0] and (q3[1]+1) == rfold[1]:
				#print ( q3 )

				curre = ( ( rfold[2] * rfold[2] ) - q3[2] )
				#print (curre)
				curre = curre * curre

				result_error = result_error + curre

				#print (rfold, q3, curre, result_error )
				break

	
	result_error = math.sqrt( result_error )

	#print ( result_error )
	#sys.exit()
	return result_error

# @TO_STUDENT: You can test your methods below by calling methods. Workflow is given already (can be changed).
# @TO_STUDENT: Everything below this point will not be considered as a solution and will be deleted for grading.
# @TO_STUDENT: Advise: be careful with idents. Use only tabs, or only FOUR spaces. NEVER mix them.

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

for k in ( 10, 200000):
	e = []
	for j in range(0, 1):
		ce = 0.0
		cmd = ' RNAsubopt -p ' + str(k) + ' < seq.fasta > t '
		os.system( cmd )

		subopt_result_filepath = "/Users/faizy/Desktop/comp598/a1/t"
		dot_ps_filepath = "/Users/faizy/Desktop/comp598/a1/dot.ps"

		# parsing RNAsubopt result file
		subopt_result = parse_subopt_result_file(subopt_result_filepath)

		#print (subopt_result) 
		#validate_Q3_1_input_format(subopt_result) 


		## solving quesion Q3_1
		q3_1_result = get_answer_Q3_1(subopt_result)

		#print ( q3_1_result )
		#print ( "\n \n \n \n \n " )
		#
		## parsing dot.ps file
		dot_ps_result = parse_dot_ps_file(dot_ps_filepath)
		#print ( dot_ps_result )
		#
		## solving question Q3_2
		q3_2_result = get_answer_Q3_2(q3_1_result, dot_ps_result)

		ce =  q3_2_result 

		e.append( ce )

		#sys.exit()
	
	print ("k: ", k, " mean: ", numpy.mean(e), " std: ", numpy.std(e) )



