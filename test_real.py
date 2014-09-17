
from Feature import FeatureSpace
import numpy as np

def LeerLC_MACHO(fid):
    saltos_linea = 3
    delimiter = ' '
    for i in range(0,saltos_linea):
        fid.next()
    LC = []

    for lines in fid:
        str_line = lines.strip().split()
        floats = map(float, str_line)
        #numbers = (number for number in str_line.split())
        LC.append(floats)

    LC = np.asarray(LC)
    return LC

#Opening the blue band
fid=open('lc_1.3568.288.B.mjd','r')

lc = LeerLC_MACHO(fid)

data  = lc[:,1]
error = lc[:,2]
mjd = lc[:,0]

#Opening the red band

fid2=open('lc_1.3568.288.R.mjd','r')

lc2 = LeerLC_MACHO(fid2)

second_data  = lc2[:,1]

#Calculating the features
a = FeatureSpace(category='all',featureList=None, automean=[0,0], StetsonL=second_data ,  B_R=second_data, Beyond1Std=error, StetsonJ=second_data, MaxSlope=mjd, LinearTrend=mjd, Eta_B_R=second_data, Eta_e=mjd, Q31B_R=second_data, PeriodLS=mjd)
a=a.calculateFeature(data)


print a.result(method='dict')

nombres = a.result(method='features')
guardar = np.vstack((nombres,a.result(method='array')))
# a=np.vstack((previous_data,a))
np.savetxt('test_real.csv', guardar, delimiter="," ,fmt="%s")