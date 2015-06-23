#!/usr/bin/env python

from Feature import FeatureSpace
import numpy as np
import ReadLC_Catalina as R
import os.path
import tarfile
import sys
import pandas as pd

import warnings
warnings.filterwarnings('ignore')
#for subdir in /n/seasfs03/IACS/TSC/MACHO/MACHO_LMC/F*; do cd "$subdir"; for f in ./*tar; do readlink -f "$f" ;done; cd ..; done; 
#for f in /n/seasfs03/IACS/TSC/MACHO/MACHO_LMC/F_1/*tar;do sbatch /n/home10/inun/Extract_features/run_extractF2.sh  "$readlink -f "$f"" ; done;

def main(argv):
	path = argv + '/'

	count = 0
	check = False

	for filename in os.listdir(path):
	
		[mag, time, error] =  R.ReadLC_Catalina(path+filename)

		a = FeatureSpace(Data=['magnitude', 'time', 'error'], featureList=None)
		lc = np.array([mag,time,error])

		try:
			a=a.calculateFeature(lc)
			idx = filename.split('.')[0]
			count = count + 1
			

			if count == 1:
					
				df = pd.DataFrame(np.asarray(a.result(method='array')).reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index =[idx])
			
			else:
				
				df2 = pd.DataFrame(np.asarray(a.result(method='array')).reshape((1,len(a.result(method='array')))), columns = a.result(method='features'), index =[idx])
				df = pd.concat([df, df2])
			
			check = True	
				
		except:
			pass

	if check:
		file_name = path.split('/')[4] + '.csv'
		#df.to_csv('/n/home10/inun/Extract_features/'+file_name)
		df.to_csv('/n/regal/TSC/Catalina_features/'+file_name)

if __name__ == "__main__":
	main(sys.argv[1:][1]) 

