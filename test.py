

from Feature import FeatureSpace
import numpy as np

#"test"

data = np.random.uniform(-5,-3, 1000)
second_data = np.random.uniform(-5,-3, 1000)
error= np.random.uniform(0.000001,1, 1000)
mjd= np.random.uniform(40000,50000, 1000)



a = FeatureSpace(category='all',featureList=None, automean=[0,0], StetsonL=second_data ,  B_R=second_data, Beyond1Std=error, StetsonJ=second_data, MaxSlope=mjd, LinearTrend=mjd, EtaB_R=second_data)
# a = FeatureSpace(category='basic', automean=[0,0])
#print a.featureList
a=a.calculateFeature(data)
#print a.result(method='')

print a.result(method='dict')
print data
print mjd