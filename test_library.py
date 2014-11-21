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
	data = np.random.normal(size=10000)
	mjd=np.arange(10000)
	error = np.random.normal(loc=0.01, scale =0.8, size=10000)
	second_data = np.random.normal(size=10000)
	mjd2=np.arange(10000)
	error2 = np.random.normal(loc=0.01, scale =0.8, size=10000)
	aligned_data = data
	aligned_second_data = second_data
	aligned_mjd = mjd
	return data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd

def test_con(fake_lc):
	# data, mjd, error, second_data, aligned_data, aligned_second_data, aligned_mjd = fake_lc()

	a = FeatureSpace(featureList=['Con'] , Con=1)
	a=a.calculateFeature(fake_lc[0])

	assert(a.result(method='array') >= 0.043 and a.result(method='array') <= 0.046)
