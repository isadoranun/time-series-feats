#!/usr/bin/env python

from Feature import FeatureSpace
import numpy as np
from import_lc_cluster import ReadLC_MACHO
from PreprocessLC import Preprocess_LC
from alignLC import Align_LC
import os.path
import tarfile
import sys
import pandas as pd

import warnings
warnings.filterwarnings('ignore')


path = '/Volumes/LaCie/MACHO_LMC/F_1/'
#path = '/Users/isadoranun/Dropbox/lightcurves/'
#path2 = '/Volumes/LaCie/MACHO_LMC/'
#path = os.getcwd()

count = 0
check = False

for j in os.listdir(path):
    
    if tarfile.is_tarfile(path + j):
    # .endswith(".tar"):
    	df = []
    	contador = 0
    	count = count + 1

        tar = tarfile.open(path + j, 'r')

        for member in tar.getmembers():

        	if member.name.endswith("B.mjd"):

        		id = member.name.split('lc_')[1]


        		for member2 in tar.getmembers():

        			if member2.name == (member.name[:-5] + 'R.mjd'):	

	    				f = tar.extractfile(member)
	    				g = tar.extractfile(member2)

	    				
	    				content1 = f.read().split('\n')
	    				content2 = g.read().split('\n')
        				
	    				lc_B = ReadLC_MACHO(content1)
	    				lc_R = ReadLC_MACHO(content2)
	    				[mag, time, error] = lc_B.ReadLC()
	    				[mag2, time2, error2] = lc_R.ReadLC()

	    				preproccesed_mag = Preprocess_LC(mag, time, error)
	    				[mag, time, error] = preproccesed_mag.Preprocess()

	    				preproccesed_mag = Preprocess_LC(mag2, time2, error2)
	    				[mag2, time2, error2] = preproccesed_mag.Preprocess()

	    				if len(mag) != len(mag2):
	    					[aligned_mag, aligned_mag2, aligned_time] = Align_LC(time, time2, mag, mag2, error, error2)
	    				else:
	    					aligned_mag = mag
	    					aligned_mag2 = mag2
	    					aligned_time = time

	    				lc = np.array([mag,time,error,mag2,aligned_mag, aligned_mag2, aligned_time])

	    				a = FeatureSpace(mag='all',featureList=None)

	    				try:
	    					a=a.calculateFeature(lc)
	    					idx = [id[:-6]]
	    					contador = contador + 1
	    					check = True


	    					if contador == 1:
	    						print "contador1"
	    						df = pd.magFrame(a.result(method='array').reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index =[idx])
	    						print "hice mi primer mag frame"	
	    					else:
	    						df2 = pd.magFrame(a.result(method='array').reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index =[idx])
	    						df = pd.concat([df, df2])

	    				except:
	    					pass
	if check:
		folder = (member.name.split('lc')[0]).split('/')[0]
		field = (member.name.split('lc')[0]).split('/')[1]
		file_name = folder + '_' + field + '.csv'
		df.to_csv(file_name)
		check = False


	# if count == 1:
	# 	folder = (member.name.split('lc')[0]).split('/')[0]
	# 	field = (member.name.split('lc')[0]).split('/')[1]
	# 	file_name = folder + '_' + field + '.csv'
	# 	nombres = np.hstack(("MACHO_Id" , a.result(method='features')))
	# 	guardar = np.vstack((nombres, guardar[1:]))
	# 	np.savetxt(file_name, guardar, delimiter="," ,fmt="%s")
	# 	guardar = np.zeros(shape=(1,2))
	# else:
	# 	nombres = np.hstack(("MACHO_Id" , a.result(method='features')))
	# 	my_mag = np.genfromtxt(file_name, delimiter=',', dtype=None)
	# 	guardar = np.vstack((nombres, my_mag[1:], guardar[1:] ))
	# 	np.savetxt(file_name, guardar, delimiter="," ,fmt="%s")
	# 	guardar = np.zeros(shape=(1,2))