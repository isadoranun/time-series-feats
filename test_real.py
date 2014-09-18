
from Feature import FeatureSpace
import numpy as np
from import_lightcurve import LeerLC_MACHO

#Opening the light curve
lc = LeerLC_MACHO('1.3444.614')
[data, mjd, error,second_data] = lc.leerLC()
 

#Calculating the features
a = FeatureSpace(category='all',featureList=None, automean=[0,0], StetsonL=second_data ,  B_R=second_data, Beyond1Std=error, StetsonJ=second_data, MaxSlope=mjd, LinearTrend=mjd, Eta_B_R=second_data, Eta_e=mjd, Q31B_R=second_data, PeriodLS=mjd)
a=a.calculateFeature(data)

print a.result(method='dict')

# nombres = a.result(method='features')
# guardar = np.vstack((nombres,a.result(method='array')))
# # a=np.vstack((previous_data,a))
# np.savetxt('test_real.csv', guardar, delimiter="," ,fmt="%s")