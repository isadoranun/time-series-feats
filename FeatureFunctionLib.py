import os,sys,time
import numpy as np
import pandas as pd
from scipy import stats

from Base import Base

class Rcs(Base):
    def __init__(self):
        self.category='timeSeries'
    def fit(self, data):
        sigma = np.std(data)
        N = len(data)
        m = np.mean(data)
        s = (np.cumsum(data)-m)*1.0/(N*sigma)
        R = np.max(s) - np.min(s)
        return R
   
class StestonK(Base):
    def __init__(self):
        self.category='timeSeries'
    def fit(self, data):
        N = len(data)
        sigmap = np.sqrt(N*1.0/(N-1)) * (data-np.mean(data))/np.std(data)
    
        K = 1/np.sqrt(N*1.0) * np.sum(np.abs(sigmap)) / np.sqrt(np.sum(sigmap**2))

        return K


class automean(Base):
    '''
    This is just a prototype, not a real feature
    '''
    def __init__(self, length): 
        self.category='basic'
        if len(length)!=2:
            print "need 2 parameters for feature automean"
            sys.exit(1)
        self.length = length[0]
        self.length2 = length[1]
    def fit(self, data):
        return np.mean(data)+self.length+self.length2

class meanvariance(Base):
    def __init__(self): 
        self.category='basic'
      
    def fit(self, data):
        return np.std(data)/np.mean(data)



    
class autocor(Base):
    def __init__(self, lag=1):
        self.category='timeSeries'
    def fit(self, data):
        return np.correlate(data, data)[0]

class StetsonL(Base):
    def __init__(self, second_data):
        self.category='timeSeries'
        if second_data is None:
            print "please provide another data series to compute StetsonL"
            sys.exit(1)
        self.data2 = second_data

    def fit(self, data):
        if len(data) != len(self.data2) :
            print " the lengh of 2 data series are not the same"
            sys.exit(1)

        N = len(data)
        sigmap = np.sqrt(N*1.0/(N-1)) * (data-np.mean(data))/np.std(data)
        sigmaq = np.sqrt(N*1.0/(N-1)) * (self.data2-np.mean(self.data2))/np.std(self.data2)
        sigma_i = sigmap * sigmaq
        J= 1.0/len(sigma_i) * np.sum(np.sign(sigma_i) * np.sqrt(np.abs(sigma_i)))
        K = 1/np.sqrt(N*1.0) * np.sum(np.abs(sigma_i)) / np.sqrt(np.sum(sigma_i**2))

        print J, K
        return J*K/0.798        

class Con(Base):
    '''
    Index introduced for selection of variable starts from OGLE database. 
    To calculate Con, we counted the number of three consecutive starts that are out of 2sigma range, and normalized by N-2
    '''
    def __init__(self, consecutiveStar=3):
        self.category='timeSeries'
        self.consecutiveStar = consecutiveStar

    def fit(self, data):

        N = len(data)
        if N < self.consecutiveStar:
            return 0
        sigma = np.std(data)
        m = np.mean(data)
        count=0
        
        for i in xrange(N-self.consecutiveStar+1):
            flag = 0
            for j in xrange(self.consecutiveStar):
                if (data[i+j] > m+2*sigma or data[i+j] < m-2*sigma) :
                    flag = 1
                else:
                    flag=0
                    break
            if flag:
                count = count+1
        return count*1.0/(N-self.consecutiveStar+1)


class VariabilityIndex(Base):
    '''
    The index is the ratio of mean of the square of successive difference to the variance of data points
    '''
    def __init__(self):
        self.category='timeSeries'
        

    def fit(self, data):

        N = len(data)
        sigma2 = np.var(data)
        
        return 1.0/((N-1)*sigma2) * np.sum(np.power(data[1:] - data[:-1] , 2))


class B_R(Base):
    '''
    average color for each MACHO lightcurve 
    mean(B1) - mean(B2)
    '''
    def __init__(self, second_data):
        self.category='basic'
        if second_data is None:
            print "please provide another data series to compute B_R"
            sys.exit(1)
        self.data2 = second_data
      
    def fit(self, data):
        return np.mean(data) - np.mean(self.data2)


# The categories of the following featurs should be revised

class Amplitude(Base):
    '''
    Half the difference between the maximum and the minimum magnitude
    '''

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        return (np.max(data) - np.min(data)) / 2

class Beyond1Std(Base):
    '''
    Percentage of points beyond one st. dev. from the weighted (by photometric errors) mean
    '''

    def __init__(self, error):
        self.category='basic'
        if error is None:
            print "please provide the measurement erros to compute Beyond1Std"
            sys.exit(1)
        self.error = error

    def fit(self, data):
        n = len(data)

        weighted_mean = np.average(data, weights= 1 / self.error**2)

        # Standard deviation with respect to the weighted mean
        var = 0
        for i in xrange(n):
            var += ((data[i]) - weighted_mean)**2
        std = np.sqrt( (1.0/(n-1)) * var )

        fraction = 0.0
        for i in xrange(n):
            if data[i] > weighted_mean + std or data[i] < weighted_mean - std:
                fraction += 1

        return fraction / n

class SmallKurtosis(Base):
    '''
    small sample kurtosis of the magnitudes. See http://www.xycoon.com/peakedness_small_sample_test_1.htm
    '''

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        n = len(data)
        mean = np.mean(data)
        std = np.std(data)

        suma = 0
        for i in xrange(n):
            suma += ((data[i] - mean) / std)**4

        c1 = float(n*(n + 1)) / ((n - 1)*(n - 2)*(n - 3))
        c2 = float(3 * (n - 1)**2) / ((n-2)*(n-3))

        return c1 * suma - c2

class Std(Base):
    '''
    standard deviation of the magnitudes.
    '''

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        return np.std(data)

class Skew(Base):
    '''
    skewness of the magnitudes
    '''

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        return stats.skew(data)

class StetsonJ(Base):
    '''
    Stetson (1996) variability index, a robust standard deviation
    '''
    def __init__(self, second_data):
        self.category='timeSeries'
        if second_data is None:
            print "please provide another data series to compute StetsonJ"
            sys.exit(1)
        self.data2 = second_data

    def fit(self, data):
        if len(data) != len(self.data2) :
            print " the lengh of 2 data series are not the same"
            sys.exit(1)
        
        N = len(data)

        sigmap = np.sqrt(N*1.0/(N-1)) * (data-np.mean(data))/np.std(data)
        sigmaq = np.sqrt(N*1.0/(N-1)) * (self.data2-np.mean(self.data2))/np.std(self.data2)
        
        sigma_i = sigmap * sigmaq
        
        J= 1.0/len(sigma_i) * np.sum(np.sign(sigma_i) * np.sqrt(np.abs(sigma_i)))

        return J

class MaxSlope(Base):
    '''
    Examining successive (time-sorted) magnitudes, the maximal first difference (value of delta magnitude over delta time)
    '''

    def __init__(self, mjd):
        self.category='timeSeries'
        if mjd is None:
            print "please provide the measurement times to compute MaxSlope"
            sys.exit(1)
        self.mjd = mjd

    def fit(self, data):
        max_slope = 0

        index = self.mjd

        for i in xrange(len(data) - 1):
            slope = float(data[i+1] - data[i]) / (index[i+1] - index[i])

            if slope > max_slope:
                max_slope = slope

        return max_slope

class MedianAbsDev(Base):

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        median = np.median(data)

        devs = []
        for i in xrange(len(data)):
            devs.append(abs(data[i] - median))

        return np.median(devs)

class MedianBRP(Base):
    '''
    Median buffer range percentage
    fraction (<= 1) of photometric points within amplitude/10 of the median magnitude
    '''

    def __init__(self):
        self.category='basic'

    def fit(self, data):
        median = np.median(data)
        amplitude = ( np.max(data) - np.min(data) ) / 10
        n = len(data)

        fraction = 0.0
        for i in xrange(n):
            if data[i] < median + amplitude and data[i] > median - amplitude:
                fraction += 1

        return fraction / n

class PairSlopeTrend(Base):
    '''
    considering the last 30 (time-sorted) measurements of source magnitude, 
    the fraction of increasing first differences minus the fraction of decreasing first differences.
    '''

    def __init__(self):
        self.category='timeSeries'

    def fit(self, data):
        data_last = data[-30:]

        inc = 0.0
        dec = 0.0

        for i in xrange(29):
            if data_last[i + 1] - data_last[i] > 0:
                inc += 1
            else:
                dec += 1

        return (inc - dec) / 30
