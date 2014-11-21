#!/usr/bin/env python

from Feature import FeatureSpace
import numpy as np
from import_lc_cluster import LeerLC_MACHO
from PreprocessLC import Preprocess_LC
from alignLC import Align_LC
import os.path
import tarfile
import sys
import pandas as pd
import pytest

@pytest.fixture
def fake_lc():
	data = np.random.normal(size=5000)
	mjd=np.arange(5000)
	error = np.random.normal(loc=0.01, scale =0.8, size=5000)
	second_data = np.random.normal(size=5000)
	mjd2=np.arange(5000)
	error2 = np.random.normal(loc=0.01, scale =0.8, size=5000)
	aligned_data = data
	aligned_second_data = second_data
	aligned_mjd = mjd
	return data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd

def test_Amplitude(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Amplitude'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Autocor(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Autocor'] )
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Automean(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Automean'] , Automean=[0,0])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_B_R(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['B_R'] , B_R=second_data)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Beyond1StdL(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Beyond1StdL'] , Beyond1StdL=error)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Bmean(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Bmean'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_CAR(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['CAR_sigma', 'CAR_tau', 'CAR_tmean'] , CAR_sigma=[mjd, error])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)


def test_Con(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Con'] , Con=1)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Eta_B_R(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Eta_B_R'] , Eta_B_R=[aligned_second_data, aligned_data, aligned_mjd])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Eta_e(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Eta_e'] )
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_LinearTrend(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['LinearTrend'] , LinearTrend = mjd)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_MaxSlope(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['MaxSlope'] , MaxSlope=mjd)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Meanvariance(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Meanvariance'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_MedianAbsDev(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['MedianAbsDev'] , MaxSlope=mjd)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_MedianBRP(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['MedianBRP'] , MaxSlope=mjd)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_PairSlopeTrend(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['PairSlopeTrend'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_PercentAmplitude(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['PercentAmplitude'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_PercentDifferenceFluxPercentile(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['PercentDifferenceFluxPercentile'])
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Period_Psi(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['PeriodLS', 'Period_fit','Psi_CS','Psi_eta'], PeriodLS = mjd, Psi_CS= mjd)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)

def test_Q31(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Q31'])
	a=a.calculateFeature(fake_lc[0])

def test_Q31B_R(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Q31B_R'], Q31B_R = [aligned_second_data, aligned_data])
	a=a.calculateFeature(fake_lc[0])

def test_Rcs(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Rcs'])
	a=a.calculateFeature(fake_lc[0])

def test_Skew(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Skew'])
	a=a.calculateFeature(fake_lc[0])

def test_SlottedA(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['SlottedA'], SlottedA = [mjd, 1])
	a=a.calculateFeature(fake_lc[0])

def test_SmallKurtosis(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['SmallKurtosis'])
	a=a.calculateFeature(fake_lc[0])

def test_Std(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Std'])
	a=a.calculateFeature(fake_lc[0])

