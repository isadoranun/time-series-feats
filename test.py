

from Feature import FeatureSpace
import numpy as np



data = np.random.randint(0,10000, 100000000)
a = FeatureSpace(category='all', automean=[0,0])
print a.featureList
a=a.calculateFeature(data)
print a.result(method='')
print a.result(method='dict')
