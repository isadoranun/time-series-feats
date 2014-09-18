
from Feature import FeatureSpace
import numpy as np

class LeerLC_MACHO:


    def __init__(self,id):

        self.id=id

    def leerLC(self):

        # Opening the blue band
        fid = open('lc_'+self.id+'.B.mjd','r')

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

        data  = LC[:,1]
        error = LC[:,2]
        mjd = LC[:,0]

        # Opening the red band
        fid2=open('lc_'+self.id+'.R.mjd','r')

        saltos_linea = 3
        delimiter = ' '
        for i in range(0,saltos_linea):
            fid2.next()
        LC2 = []

        for lines in fid2:
            str_line = lines.strip().split()
            floats = map(float, str_line)
            #numbers = (number for number in str_line.split())
            LC2.append(floats)

        LC2 = np.asarray(LC2)

        second_data  = LC2[:,1]

        return [data, mjd, error,second_data]

    