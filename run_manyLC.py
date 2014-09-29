
from Feature import FeatureSpace
import numpy as np
from import_lightcurve import LeerLC_MACHO
from PreprocessLC import Preprocess_LC
import os.path

count = 0

for i in os.listdir(os.getcwd()):
# for i in os.listdir('/Users/isadoranun/Desktop/Features/time-series-feats/'):
    if i.endswith("B.mjd"): 

        count = count + 1

        lc_B = LeerLC_MACHO(i)
        lc_R = LeerLC_MACHO(i[:-5] + 'R.mjd')

#Opening the light curve

        [data, mjd, error] = lc_B.leerLC()
        [data2, mjd2, error2] = lc_R.leerLC()

        preproccesed_data = Preprocess_LC(data, mjd, error)
        [data, mjd, error] = preproccesed_data.Preprocess()

        preproccesed_data = Preprocess_LC(data2, mjd2, error2)
        [second_data, mjd2, error2] = preproccesed_data.Preprocess()

        a = FeatureSpace(category='all',featureList=None, automean=[0,0], StetsonL=second_data ,  B_R=second_data, Beyond1Std=error, StetsonJ=second_data, MaxSlope=mjd, LinearTrend=mjd, Eta_B_R=second_data, Eta_e=mjd, Q31B_R=second_data, PeriodLS=mjd, CAR_sigma=[mjd, error])
        a=a.calculateFeature(data)

        if count == 1:
        	nombres = a.result(method='features')
        	guardar = np.vstack((nombres, a.result(method='array')))
        	np.savetxt('test_real.csv', guardar, delimiter="," ,fmt="%s")

        else:
        	my_data = np.genfromtxt('test_real.csv', delimiter=',')
        	guardar = np.vstack((nombres,my_data[1:], a.result(method='array')))
        	np.savetxt('test_real.csv', guardar, delimiter="," ,fmt="%s")

#B_R = second_data, Eta_B_R = second_data, Eta_e = mjd, MaxSlope = mjd, PeriodLS = mjd, Q31B_R = second_data, StetsonJ = second_data, StetsonL = second_data)