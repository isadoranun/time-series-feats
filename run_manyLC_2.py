
from Feature import FeatureSpace
import numpy as np
from import_lightcurve import ReadLC_MACHO
from PreprocessLC import Preprocess_LC
from alignLC import Align_LC
import os.path
import pandas as pd

contador = 0
folder = 1
df = []
#guardar = np.zeros(shape=(1,42))

path = '/Users/isadoranun/Google Drive/lightcurves/'

ex_ts = pd.read_csv('/Users/isadoranun/Documents/MACHO_ts2.csv',index_col=0) 
ids = ex_ts.index.values


for j in os.listdir(path)[2:]:
    
    if os.path.isdir(path + j):

        print j

        for i in os.listdir(path + j):

            if i.endswith("B.mjd") and not i.startswith('.') and os.path.isfile(path + j +'/'+ i[:-5] + 'R.mjd'):


                a =[num for num, s in enumerate(ids) if i[3:-6] in s]

                print i[3:-5], a

                if a == []:

                    print 'nuevo'

                    lc_B = ReadLC_MACHO(path + j +'/'+ i[:])
                    lc_R = ReadLC_MACHO(path + j +'/'+ i[:-5] + 'R.mjd')

            #Opening the light curve

                    [mag, time, error] = lc_B.ReadLC()
                    [mag2, time2, error2] = lc_R.ReadLC()

                    preproccesed_data = Preprocess_LC(mag, time, error)
                    [mag, time, error] = preproccesed_data.Preprocess()

                    preproccesed_data = Preprocess_LC(mag2, time2, error2)
                    [mag2, time2, error2] = preproccesed_data.Preprocess()

                    if len(mag) != len(mag2):
                        [aligned_mag, aligned_mag2, aligned_time, aligned_error, aligned_error2] = Align_LC(time, time2, mag, mag2, error, error2)
                    else:
                        aligned_mag = mag
                        aligned_mag2 = mag2
                        aligned_time = time
                        aligned_error = error
                        aligned_error2 = error2

                    lc = np.array([mag, time, error, mag2, aligned_mag, aligned_mag2, aligned_time, aligned_error, aligned_error2])

                    a = FeatureSpace(Data='all', featureList=None)

                    try:
                        a=a.calculateFeature(lc)
                        idx = i[3:-6]
                        # print 'feat', (a.result(method='array')).shape
                        contador = contador + 1


                        if contador == 1:
                            df = pd.DataFrame(np.asarray(a.result(method='array')).reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index=[idx])
                        else:
                            print contador
                            df2 = pd.DataFrame(np.asarray(a.result(method='array')).reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index =[idx])
                            df = pd.concat([df, df2])

                    except:
                        pass
                
                 #   a = FeatureSpace(category='all',featureList=None, automean=[0,0], StetsonL=[aligned_mag2, aligned_mag] ,  B_R=mag2, Beyond1Std=error, StetsonJ=[aligned_mag2, aligned_mag], MaxSlope=time, LinearTrend=time, Eta_B_R=[aligned_mag2, aligned_mag], Eta_e=time, Q31B_R=[aligned_mag2, aligned_mag], PeriodLS=time, CAR_sigma=[time, error], SlottedA = time)

df.to_csv('/Users/isadoranun/Documents/MACHO_ts3.csv')     
                


        #B_R = mag2, Eta_B_R = mag2, Eta_e = time, MaxSlope = time, PeriodLS = time, Q31B_R = mag2, StetsonJ = mag2, StetsonL = mag2)