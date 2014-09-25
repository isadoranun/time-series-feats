import numpy as np

class Preprocess_LC:

	def __init__(self, mjd, data, error):

		self.N = len(mjd)
		self.m = np.mean(error)
		self.mjd = mjd
		self.data = data
		self.error = error

	def Preprocess(self):

	t_out = []
	data_out = []
	error_out = []	

		for i in xrange(N):
			if self.error[i] < (3 * self.m) & (np.absolute(self.data[i] - np.mean(data)) / np.std(x)) < 5 :
				t_out.append(self.t[i])
				data_out.append(self.data[i])
				error_out.append(self.error[i])



	return [data_out, t_out, error_out]
