
import sys, os, numpy as np, random #, collections

class AnalyseAncestor:

	def __init__(self):

		return


	''' Objective 1 '''
	''' set tree as false for normal stockholm and true for modified stockholm'''

	def stockholmParser( self, filename, tree ):

		alignment = dict(); consensus = []

		for line in open(filename, 'r'):
			#print (line)
			line = line.replace('\n', '')

			if len(line) < 1:
				continue

			elif line[0].isalpha():
				line = line.split()
				if line[0] not in alignment.keys():
					if not tree:
						alignment[ line[0] ] = line[1]
					else:
						alignment[ line[1] ] = line[2]
				else:
					alignment[ line[0] ] = alignment[ line[0] ] + line[1]

			elif line[:12] == '#=GC SS_cons':
				line = line.split()
				consensus.append( line[2] )

		consensus = ''.join(consensus)

		return [ alignment, consensus ]
	
	''' Objective 2, 3 and 4'''


	''' takes tree list for eg. ((A,B),C)
	    returns list of nodes
	'''

	def parseTree( self, treelist):

		tfid = open( treelist, 'r')
		treelist = tfid.readlines(); treelist = ''.join(treelist)
		tfid.close()

		nodestack = []
		parent = dict() #collections.OrderedDict()
		leaves = []

		for t in treelist:
			if t.isalpha():
				nodestack.append( t )
				leaves.append( t )

			if t == ')':
				#parent = dict()
				right_child = nodestack.pop(); 
				left_child = nodestack.pop()
				parentid = left_child + right_child
				parent[ parentid ] = [left_child, right_child]
				nodestack.append( parentid )



		return parent, leaves

	''' parent is a dict that points to its left and right child name
	    leaves are the list of leave names
	    cost_matrix is the penalty matrix
	    alignment is a dictionary with species as keys and values as sequences
	'''

	def sankoff(self, parent, leaves, cost_matrix, alignment): 
		
		inf = float( 'inf' )
		#cost_matrix = [[0,2.5,1,2],[2,0,2,1],[1,2,0,2],[2,1,2,0]]
		len_cost_matrix = len(cost_matrix)

		if len_cost_matrix == 4:
			bpindex = {'A':0,'C':1,'G':2,'U':3}
			bpval = {0:'A', 1:'C', 2:'G', 3:'U' }
		else:
			bpindex = {'A':0,'C':1,'G':2,'U':3,'-':4}
			bpval   = {0:'A', 1:'C', 2:'G', 3:'U',4:'-' }


		# alignment length
		for l in leaves:
			lenalign = len(alignment[l])
			break

		for pos in range( lenalign ):
			cost = dict()
			
			for l in leaves:
       				cost[l] = np.asarray( [inf]* len_cost_matrix, dtype = float)
       				cost[l][ bpindex[ alignment[l][pos] ] ] = 0
			
			for p in sorted( parent, key = len):
				left_cost = cost[ parent[p][0] ]
				right_cost = cost[ parent[p][1] ]
				#parentcost = np.zeros(4, dtype = float)
				currcost = [0] * len_cost_matrix
				
				for ci, c in enumerate(cost_matrix):
					#print (c)
					c = np.asarray(c, dtype = float)
					min_left = min( left_cost + c )
					min_right = min( right_cost + c)
					#print ( c, left_cost, min_left, right_cost, min_right )

					#cost_sum = left_cost + 2*c + right_cost

					currcost[ci] = min_left + min_right #min( cost_sum )

				cost[p] = np.asarray( currcost, dtype = float )
				bp = bpval[ currcost.index( min(currcost) ) ]
				
				if p not in alignment.keys():
       					alignment[p] = 	bp
				else:
       					alignment[p] = alignment[p] + bp 

			#for c,v in cost.items():
			#	print (c, v)

		return alignment
			
			
		

	''' Objective 5 '''

	# compute structure using RNAfold
	def computeStruc( self, sequence):
		cmd = "echo "+ sequence + " > tmpq6; RNAfold < tmpq6 | tail -n1 | cut -d' ' -f1 > tmpq6out"
		os.system( cmd )
		for line in open( 'tmpq6out', 'r'):
			res = line
			res = res.replace('\n', '')
			break
		return res

	# compute distance using RNAdistance
	def bp_distance( self, fold, target ):
		target = target.replace('<','('); target = target.replace('>',')'); 
		target = target.replace('_','.'); target = target.replace('-','.'); target = target.replace(':','.'); target = target.replace(',','.');
		fout = open('tmp2q6', 'w')
		fout.write(fold + '\n' + target)
		fout.close()
		cmd = "RNAdistance < tmp2q6 | cut -d' ' -f2  > tmp2q6out"
		os.system( cmd )
		res = ''
		for line in open( 'tmp2q6out', 'r'):
			res = line
			res = res.replace('\n', '')
			break
		#print ( fold, target, res )
		return int(res)
	

	''' Objective 6 and 7 '''

	''' Takes a sequence and a structure. Forces the bp dependency of 
	structure onto sequence '''

	def validpair( self, a, b ):
		valid = [ 'au', 'ua', 'cg', 'gc', 'gu', 'ug' ]

		if [a,b] in valid:
			return True
		else:
			return False

	def getvalidpair( self, a, b):

		glist = [ ['G', 'C'], ['G', 'U'] ]
		ulist = [ ['U', 'A'], ['U', 'G'] ]
		gaplist = [ ['A', 'U'], ['U', 'A'], ['C', 'G'], ['G', 'C'], ['G', 'U'], ['U', 'G'] ]

		if a == 'A':
			return [a, 'U']
		elif a == 'C':
			return [ a, 'G' ]
		elif a == 'G':
			return glist[ random.randint(0,1) ]
		elif a == 'U':
			return ulist[ random.randint(0,1) ]
		elif a == '_':
			if b == '_':
				return gaplist[ random.randint(0,5) ]
			else:
				return self.getvalidpair( b, a )



	def extendedSankoffwithBpDependency(self, sequence, target ):
		target = target.replace('<','('); target = target.replace('>',')'); 
		target = target.replace('_','.'); target = target.replace('-','.'); target = target.replace(':','.');target = target.replace(',','.');
		#print ( target )
		sequence = list( sequence )
		slist = []
		spair = []
		
		for si, s in enumerate( target ):

			if s == ')':
				left = slist.pop()
				right = si
				spair.append( [left,right] )

			elif s == '(':
				slist.append( si )

		for pair in spair:
			left = pair[0]; right = pair[1]

			if not self.validpair( sequence[left], sequence[right] ):
				[sequence[left], sequence[right]] = self.getvalidpair( sequence[left], sequence[right] )

		return ''.join(sequence)

	''' Objective 8 '''
	''' Identify structure stability '''

	def computeMFE( self, sequence):
		cmd = "echo "+ sequence + " > rnafoldtmp; RNAfold -p < tmpq6 | tail -n1 | cut -d' ' -f8 > rnafoldtmpout"
		os.system( cmd )
		for line in open( 'rnafoldtmpout', 'r'):
			res = line
			res = res.replace('\n', ''); res = res.replace(';', '')
			break
		return float(res)


		return


if __name__ == '__main__':

	if sys.argv[1] == '-h':
		print ( 'python rnaAncestor.py stockholmfilename treefilename normal_stockholmfile' )
		print (' set normal_stockholmfile as True for just parsing of stockholm file' )
		print (' For other analyses set normal_stockholmfile as False and use modified stockholm file' )
		sys.exit()

	filename = sys.argv[1] #'comp561-master/RF00128_seed.stockholm.txt' #sys.argv[1]
	if sys.argv[3] == 'True':
		normal_stockholm = True
	else:
		normal_stockholm = False

	ancestorobject = AnalyseAncestor()

	[alignment, consensus] = ancestorobject.stockholmParser( filename, normal_stockholm )

	gc_count = 0

	# Objective 1

	print('your file: ', filename )
	print('consensus structure: ', consensus )
	print('Alignments: ')
	for k,v in alignment.items():
		print (k +' ' + v)
		#print ( len(v) )
		#sys.exit()
		gc_count = gc_count + v.count('G')
		gc_count = gc_count + v.count('C')
	print ( 'GC count: ', gc_count )

	if not normal_stockholm:
		sys.exit()
	#sys.exit()

	#[alignment, consensus] = ancestorobject.stockholmParser( filename, True )

	treefile = sys.argv[2]

	[parent, leaves ] = ancestorobject.parseTree( treefile )

	for p in sorted( parent, key = len ) :
		print (p)

	for l in leaves:
		print ( l )

	# Objective 2, 3 and 4
	# Execute sankoff
	cost_matrix = [[0,2,1,2,2],[2,0,2,1,2],[1,2,0,2,2],[2,1,2,0,2],[2,2,2,2,0]]

	allalignment = ancestorobject.sankoff( parent, leaves, cost_matrix, alignment)

	print ('All Alignments')
	for k,v in allalignment.items():
		print (k,v)

	newallalignment = dict()

	mfelist = []
	dist1 = []
	dist2 = []

	# Objective 5,6,7,8

	# compute structure
	for p in sorted( parent, key = len):
		
		parent_struc = ancestorobject.computeStruc( allalignment[p] ) 
		print (p, parent_struc)
		print( ancestorobject.bp_distance( parent_struc, consensus ) )
		dist1.append( ancestorobject.bp_distance( parent_struc, consensus ) )
		#print ( 'ACACGAC' ) #allalignment[p] )
		#print ( consensus )

		newallalignment[p] = ancestorobject.extendedSankoffwithBpDependency( allalignment[p], consensus )

		parent_struc = ancestorobject.computeStruc( newallalignment[p] ) 
		print (p, parent_struc)
		dist2.append( ancestorobject.bp_distance( parent_struc, consensus ) )
		print( ancestorobject.bp_distance( parent_struc, consensus ) )
		mfe = ancestorobject.computeMFE( newallalignment[p] ) 
		print ('MFE: ', mfe )

		mfelist.append(mfe)
		#break

	print ( 'Mean MFE:', np.mean(mfelist), 'Mean dist1: ', np.mean(dist1), 'Mean dist2: ', np.mean(dist2) )

	
#
#	# test sankoff
#	cost_matrix = [[0,2.5,1,2.5],[2.5,0,2.5,1],[1,2.5,0,2.5],[2.5,1,2.5,0]]
#	alignment = {'A':'C' , 'B':'A', 'C':'C', 'D':'A', 'E':'G' }
#	allalignment = ancestorobject.sankoff( parent, leaves, cost_matrix, alignment)
#
#	print (allalignment)





#	treelist = '((A,B),C)'; alignment = { 'A':'ACGU', 'B':'ACGU', 'C':'AAGC'}
#
#	slist = ancestorobject.parseTree( treelist, alignment )
#
#	for k,v in slist.items():
#		print (k, v['id'], v['parent'])
		#for i,j in v.items():
		#	print ('val', i, j)

#	ancestorobject.sankoff()

	




