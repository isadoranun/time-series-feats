from Feature import FeatureSpace
import pandas
import numpy as np

def Table(a):

	all_features = np.random.uniform(size=(len(a.result(method='array')),1), low=1, high = 2)

	all_features[:,0] = a.result(method= 'array')
	df = pandas.DataFrame(all_features)
	df.index = a.result(method= 'features')
	df.reset_index(level=0, inplace=True)
	df.columns =["Feature", "Value"]
	pandas.set_option('display.float_format', lambda x: '%.3f' % x)
	return df