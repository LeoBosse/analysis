#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

DtoR = np.pi/180.
RtoD = 1./ DtoR

class LightRay:
	def __init__(self):

		self.N = 1000

		self.InonPola = 0.5

		self.AoLP = np.linspace(-np.pi, np.pi, self.N)
		self.intensities = signal.gaussian(self.N, 10*DtoR)

		Ns = self.N
		rs_signal = np.zeros(Ns) + 0.5 * self.InonPola
		filter_orientation = np.linspace(0, 2 * np.pi, Ns, endpoint=False)
		for i_f, f in enumerate(filter_orientation):
			for ihist, hist in enumerate(self.intensities):
				rs_signal[i_f] += hist * np.cos(self.AoLP[ihist] - f) ** 2

		plt.plot(filter_orientation, rs_signal)
		plt.plot(filter_orientation, [0.5 * self.InonPola]*Ns)
		plt.show()

		self.V = np.average(rs_signal)
		self.Vcos = np.average(rs_signal * np.cos(2 * filter_orientation))
		self.Vsin = np.average(rs_signal * np.sin(2 * filter_orientation))

		self.I0 = 2 * self.V
		self.DoLP = 100 * 2 * np.sqrt(self.Vcos ** 2 + self.Vsin ** 2) / self.V
		self.AoRD = np.arctan2(self.Vsin, self.Vcos) / 2

		print(self.I0, self.DoLP, self.AoRD)

	def PlotHist(self, ):
		f3, ax = plt.subplots(1, figsize=(16, 8))

		ax = plt.subplot(111, projection='polar')
		ax.set_theta_zero_location("N")
		ax.set_theta_direction(-1)
		if double:
			ax.set_thetamin(-180)
			ax.set_thetamax(180)
		else:
			ax.set_thetamin(-90)
			ax.set_thetamax(90)

		if not double:
			bins = self.AoLP[:self.N]
			h = self.hst
		else:
			bins = np.append(self.AoLP[:self.N], 180*DtoR + self.AoLP[:self.N])
			h = np.append(self.hst, self.hst)
		bars = ax.bar(bins, h, width=self.width)

		ax.set_title("Weighted AoRD: I0 = " + str(np.format_float_scientific(I0, precision=3)) + " DoLP = " + str(np.round(DoLP, 1)) + " AoRD = " + str(np.round(AoRD*RtoD, 1)))

		if not double:
			ax.plot([AoRD, AoRD], [0, max(self.hst)], "r")
		else:
			ax.plot([-AoRD, AoRD], [max(self.hst), max(self.hst)], "r")

		if self.save_individual_plots:
			plt.savefig(self.path + self.save_name + '_AoRD_hist.png', bbox_inches='tight')




if __name__ == "__main__":
	ray = LightRay()
