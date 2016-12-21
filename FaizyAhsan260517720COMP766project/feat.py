# Feature Importance
#from sklearn import datasets
#from sklearn import metrics

import numpy as np, sys
from sklearn.ensemble import ExtraTreesClassifier

if sys.argv[1] == '-h':
	print ( 'python feat.py traindatafilename ' )
	print ( 'compatible with python 2.7.9')
	print (' traindatafile should have 8 columns. 1 is row id, 2-7 features, and 8 is difficulty level' )
	print ('look at newworkdata.txt')
	print (' train and test files should have more than one line' )
	sys.exit()



# load the datasets
data = np.loadtxt( sys.argv[1] ) 

x = data[:, 1:7]
y = data[:, 7]

# fit an Extra Trees model to the data
model = ExtraTreesClassifier()

model.fit( x, y)

# display the relative importance of each attribute
print(model.feature_importances_)
