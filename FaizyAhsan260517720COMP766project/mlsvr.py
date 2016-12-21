import sys, numpy as np, random, pickle
from sklearn.grid_search import GridSearchCV
from sklearn.svm import SVR
from sklearn import preprocessing

if sys.argv[1] == '-h':
	print ( 'python mlsvr.py traindatafilename c gamma epsilon testfilename' )
	print (' traindatafile should have 8 columns. 1 is row id, 2-7 features, and 8 is difficulty level' )
	print ('look at newworkdata.txt')
	print (' train and test files should have more than one line' )
	sys.exit()

traindatafilename = sys.argv[1]
testdatafilename = sys.argv[5]

data = np.loadtxt( traindatafilename)
testdata = np.loadtxt( testdatafilename )


xtrain = data[:,1:7]
ytrain = data[:,7]

xtest = testdata[:,1:7]
ytest = testdata[:,7]


#print ( x.shape, y )
print ( 'data loaded.' )

print ( 'Train size: ', xtrain.shape, ytrain.shape, ' Test size: ', xtest.shape, ytest.shape )

# learn model
#svrmodel = SVR( C = clf.best_params_['C'], gamma = clf.best_params_['gamma'], epsilon = clf.best_params_['epsilon'] )
c = float(sys.argv[2])
g = float( sys.argv[3])
e = float( sys.argv[4]) 

svrmodel = SVR( C = c, gamma = g, epsilon = e ) 

svrmodel = svrmodel.fit( xtrain, ytrain )

# test model

ypred = svrmodel.predict( xtest )

# report accuracy
print ( 'R-sqrd error on Train set: ', c, g, e, round(svrmodel.score( xtrain, ytrain ), 3) )

print ( 'R-sqrd error on Test set: ', c, g, e, round(svrmodel.score( xtest, ytest ), 3) )

print ('Difficulty level of test data')
[ print(y) for y in ypred ]

# do feature selection

