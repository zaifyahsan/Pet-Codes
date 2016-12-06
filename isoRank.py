import networkx as nx, numpy as np, sys, pickle, matplotlib.pyplot as plt

from scipy import linalg as LA


def createGraph( networkfilename ):

	g = nx.Graph()
	
	#ctr = 0
	for line in open( networkfilename, 'r' ):
		#if ctr == 0:
		#	ctr = ctr + 1
		#	continue
		line = line.replace('\n', '')
		line = line.split()
		#print (line)
		g.add_edge( line[0], line[1] )

	#print (g['3'])

	return g

if __name__ == '__main__':

	n1filename = sys.argv[1]; n2filename = sys.argv[2]

	n1 = nx.Graph()

	n1 = createGraph( n1filename ) 
	#sys.exit()
	n2 = createGraph( n2filename )

	n1nodes = np.sort(n1.nodes())
	n2nodes = np.sort(n2.nodes())

	#print ( n1nodes, n2nodes )

	# get the neighbors

	n1dict = dict(); n2dict = dict()

	for item in n1nodes:
		n1dict[ item ] = n1.neighbors( item )

	for item in n2nodes:
		n2dict[ item ] = n2.neighbors( item )

	m = len(n1nodes)
	n = len(n2nodes)
	mn =  m * n #len(n1nodes) * len(n2nodes)

	#print (m, n, mn )
	#sys.exit()

	A = np.zeros( ( mn, mn ), dtype = float )


	for i in range(len(n1nodes)) :
		for j in range(len(n2nodes)):
	
			arow = np.zeros( (m, n), dtype = float )
	
			for ni in range( len(n1nodes) ):
	
				for nj in range( len(n2nodes) ):
	
					if n1nodes[ni] in n1dict[ n1nodes[i] ] and n2nodes[nj] in n2dict[ n2nodes[j] ]:
						#print ( ni, n1nodes[ni], n1dict[ n1nodes(n )
						#arow[ni, nj] = 
						a = 1.0 / ( len( n1dict[ n1nodes[ni]] ) * ( len( n2dict[ n2nodes[nj]] ) ) ) 
						#print (a)
						arow[ni, nj] = a
	
			arow = np.reshape( arow, (1, mn) )
	
			#print ( arow.shape )
	
			A[ (n)*i + j : ] = arow

	print ('Matrix A created' )

	#print ('A')
	#print ( A )

	np.save( n1filename.split('.')[0] + n2filename.split('.')[0] + 'mat', A )

	#np.load(n1filename.split('.')[0] + n2filename.split('.')[0] + 'mat')
	
		
	
	Reigval, Reigvec = LA.eig( A )

	np.save(n1filename.split('.')[0] + n2filename.split('.')[0] + 'eigvec', Reigvec[0] )

	#np.load( n1filename.split('.')[0] + n2filename.split('.')[0] + 'eigvec' )

	reqvec = Reigvec[0]

	threshold = 0.01

	reqvec = np.extract( reqvec > threshold, reqvec )

	reqveclist = []

	for rv in reqvec:
		reqveclist.append( abs(rv) )

	## plot reqveclist
	#plt.figure()

	#plt.hist( reqveclist, bins= 100 )

	#plt.savefig( n1filename.split('.')[0] + n2filename.split('.')[0] + '.png' )

	#print ( 'sorted', np.sort(reqveclist) )
	#sys.exit()

	reqcomblist = []

	for n1 in range( len(n1nodes) ):
		for n2 in range( len(n2nodes) ):
			reqcomblist.append( ( n1nodes[n1], n2nodes[n2] ) )

	
	lamda = min( reqveclist )

	#print ( len(reqveclist), len(reqcomblist), reqveclist, lamda )


	count = min( m, n ) #sum( i >= lamda for i in reqveclist )
	#print ('count: ', count, lamda)

	subgraph = nx.Graph()
	chk = 0

	print ('Required Alignment:' )

	while count > 0 and chk < len(reqveclist):

		currmaxindex = reqveclist.index( max( reqveclist ) )
		reqveclist[ currmaxindex ] = lamda - 1

		curredge = reqcomblist[ currmaxindex ]

		#print ( curredge[0], curredge[1] )
		#sys.exit()

		if not subgraph.has_node( curredge[0] ) and not subgraph.has_node( curredge[1] ):
			print ( curredge[0], curredge[1] )
			subgraph.add_edge( *curredge )
			count = count - 1
		chk = chk + 1
	





