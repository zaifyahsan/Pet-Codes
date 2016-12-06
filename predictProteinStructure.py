import re, numpy as np, matplotlib.pyplot as plt

def prepareData():

	alpha = {"A":0.,
	              "R":0.21,
	              "N":0.65,
	              "D":0.69,
	              "C":0.68,
	              "E":0.4,
	              "Q":0.39,
	              "G":1.,
	              "H":0.61,
	              "I":0.41,
	              "L":0.21,
	              "K":0.26,
	              "M":0.24,
	              "F":0.54,
	              "P":3.16,
	              "S":0.5,
	              "T":0.66,
	              "W":0.49,
	              "Y":0.53,
	              "V":0.61}
	
	#Benchmark 
	bench_seq =   'VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG'
	bench_str = '   HHHHHHHHHHHHHHTTSHHHHHHHHHHHHHHH HHHHHT HHHHT  SHHHHHH HHHHHHHHHHHHHHHHHHTTTT  HHHHHHHHHHHHHTT   HHHHHHHHHHHHHHHHHH TTTTSHHHHHHHHHHHHHHHHHHHHHHHHHT   '
	
	#print ( len(bench_seq), len(bench_str) )
	
	true_str = re.sub( r'[^H]', 'C', bench_str )
	
	#print ( true_str )

	return alpha, bench_seq, true_str

class ProteinStructure:
	
	def __init__(self, alpha):
		self.alpha = alpha

	'''
	Takes sequence and lambda
	Returns predicted structure based on lambda as a threshold
	'''
	def predictStr( self, seq, l ):
		predicted_str = ''
	
		for res in seq:
			if alpha[res] < l:
				predicted_str = predicted_str + 'H'
			else:
				predicted_str = predicted_str + 'C'
	
		return predicted_str

	'''
	Takes sequence and lambda
	Returns predicted structure with improved algo
	'''
	def improvedpredictStr( self, seq, l ):
		predicted_str = ''

		for i in range(0, len(seq) ):
			before = '' ; after = ''

			if i >= 4:
				before = seq[i-4:i]
			else:
				before = seq[:i]

			if i <= len(seq) - 5 :
				after = seq[i: i+5]
			else:
				after = seq[i:]

			subseq = before + after

			prop = [ self.alpha[x] for x in subseq ]

			avg_prop = sum( prop ) / float(len(prop))

			if avg_prop < l:
				predicted_str = predicted_str + 'H'
			else:
				predicted_str = predicted_str + 'C'
	
		return predicted_str

	def accuracy( self, pred_str, true_str, l ):
	
		#pred_str = predictStr( bench_seq, l )
	
		tp = 0; fp = 0; tn = 0; fn = 0;
	
		for i in range(0, len(true_str) ):
		
			if pred_str[i] == 'H':
				if true_str[i] == 'H':
					tp = tp + 1
				else:
					fp = fp + 1
		
			else: 
				if true_str[i] == 'C':
					tn = tn + 1
				else:
					fn = fn + 1
		
		tpr = ( tp / float( tp + fn ) ) 
		
		fpr = ( fp / float( fp + tn ) ) 
	
		return [ tpr, fpr ]

if __name__ == "__main__":

	# collect data
	alpha, bench_seq, true_str = prepareData()

	''' q2 part 1, implement the simple algorithm '''

	# take lambda value
	l = 0.66
	
	protein = ProteinStructure( alpha )
	
	# predicted structure

	pred_str = protein.predictStr( bench_seq, l )
	
	print ( 'Predicted Structure for lambda :', l, '\n', pred_str )

	''' q2 part 2, compute TPR and FPR '''
	
	tpr, fpr = protein.accuracy( pred_str, true_str, l )

	print ( 'TPR and FPR for lambda: ', l, ' are ', tpr, ' ', fpr )

	''' q2 part 3, AUC plot for different values of lambda '''

	tpr = []; fpr = []
	
	for l in np.arange( 0, 4, 0.005 ):
		pred_str = protein.predictStr( bench_seq, l )
	
		[t, f] = protein.accuracy( pred_str, true_str, l)
		tpr.append(t)
		fpr.append(f)
	
	print ( 'AUC for primary algorithm is: ', np.trapz(tpr, fpr) )

	# plot
	plt.figure()
	plt.plot( fpr, tpr )
	plt.xlabel( ' False Positive Rate ' )
	plt.ylabel( ' True Positive Rate ' )
	plt.title( ' Primary Algorithm ' )

	plt.savefig('primary.png')

	''' q2 part 4, predict structure with improved algorithm '''
	# take lambda value
	l = 0.66
	print ( 'Improved predicted structure for lambda : ', l, '\n', protein.improvedpredictStr( bench_seq, l ) )

	''' q2 part 5, compare AUC of both algorithms '''

	oldtpr = tpr; oldfpr = fpr

	tpr = []; fpr = []
	
	for l in np.arange( 0, 4, 0.005 ):
		pred_str = protein.improvedpredictStr( bench_seq, l )
	
		[t, f] = protein.accuracy( pred_str, true_str, l)
		tpr.append(t)
		fpr.append(f)
	
	print ( 'AUC for improved algorithm is: ', np.trapz(tpr, fpr) )

	# compare the plots

	plt.figure()

	plt.plot( oldfpr, oldtpr, 'r', label = 'primary algo.' )
	plt.plot( fpr, tpr, 'g', label = 'improved algo.' )

	plt.xlabel(' False Positive Rate ')
	plt.ylabel(' True Positive Rate ')
	plt.title( ' AUC comparison ' )
	plt.legend( loc = 4 )
	plt.savefig( 'auccomparison.png' )

	
	#l = 0.66
	#
	#print ( predictStr( bench_seq, l ) )
	#
	#pred_str = predictStr( bench_seq, l) 
	#
	## Compute TPR FPR
	##def computeAUC( pred_str, true_str ):
	#
	#tpr = []
	#fpr = []
	#
	##for l in np.arange( 0, 4, 0.005 ):
	#
	#
	#l = 0.66
	#
	#print ( predictStr( bench_seq, l ) )
	#
	#print ( accuracy( bench_seq, true_str, l ) )
	#
	




