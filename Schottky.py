""" 31.01.2017 KonstantinMurasov
This Module utilizes the Modules "DopingDrift.py" and "ChargeRedistribution.py"
to calculate the changes in doping density and depletio length due to migration
of doping atoms in a suffficient electric field, provided in "Units.py"

Please take note that the "Initial Conditions"-Graph is only valid without
applied external voltage.
"""


from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy

import Units as u
from SolverLU import *
from DopingDrift import *
from ChargeRedistribution import *


class Schottky(object) :
	def __init__(self, Bias, dt, nx): 
		#------------------------------------
		# 0. Initialize 
		self.X = linspace(0, u.DeviceLength, nx)
		self.dx = self.X[1] - self.X[0]
		self.nx = nx
		self.dt = dt
		self.Bias = Bias
		
		# 1. Calculate Initial Conditions
		self.Dieelek = ones(nx)*u.e0*u.Titanium_DielectricityFactor
		
		# Depletion Region Width
		self.W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight+self.Bias)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5
		
		# Initial Doping Density
		self.DopingDensity = zeros(self.nx)
		self.DopingDensity[:] = u.Titanium_DopingDensity
		
		# Initail Charge density
		self.ChargeDensity = zeros(nx)
		self.ChargeDensity[:int(self.W/u.DeviceLength*nx)] = u.Titanium_DopingCharge*u.Titanium_DopingDensity

		#Initial Potential Calculation
		self.Phi = zeros(nx)
		self.Phi[1] = self.Phi[0] = u.BarrierHeight+self.Bias; self.Phi[-1] = 0
		self.Phi = solverLU_backward(self.X, self.Phi, self.ChargeDensity, self.Dieelek)

		#Electrical Field
		E = E_Field(self.X, self.Phi)
		
		# Inital conditions Plot
		plt.figure("Initial Conditions")
		plt.plot(self.X, self.ChargeDensity)
		plt.plot(self.X, self.Phi,color = "G")
		plt.plot(self.X, E, color = "R")
		
		#Mirror Charge Correction
		#Phi = Phi - 1/((X+1e-1)*16*pi*u.e0*u.Titanium_DielectricityFactor)
		
		#Phi5 = solverLU_backward(X, Phi, ChargeDensity, Dieelek)
		#plt.plot(X, Phi5, color = "Red")	
		
	#------------------------------------
	def DopingDrift(self):
		# 2. Calculate Doping Drift
		DP = DriftRK( self.X , self.Phi, self.DopingDensity, self.dt)
		"""
		for i in range(Iterations-1) :
			print(self.Phi[-1])
			[self.ChargeDensity, self.W]  = ReloadDeltaPhi(self.X, DP, self.ChargeDensity, self.Phi, self.Dieelek, self.W, FermiLevel = self.Bias)
			self.Phi[-1] = 0.0; self.Phi[-2] = 0.0
			#print(i)
			self.Phi = solverLU(self.X, self.Phi, self.ChargeDensity, self.Dieelek)
			DP = DriftRK( self.X , self.Phi, DP, self.dt)
		"""
		[self.ChargeDensity, self.W]  = ReloadDeltaPhi(self.X, DP, self.ChargeDensity, self.Phi, self.Dieelek, self.W, FermiLevel = self.Bias)
		self.Phi[-1] = 0.0; self.Phi[-2] = 0.0
		self.Phi = solverLU_backward(self.X, self.Phi, self.ChargeDensity, self.Dieelek)
		self.DopingDensity = DP
		
	def PlotResults(self):
		plt.figure("Oxygen Vacancy Movement")
		
		plt.plot(self.X, self.DopingDensity, color = "Grey")
		plt.plot(self.X, self.ChargeDensity, color = "Red")
		#plt.plot(self.X, DP, color = "Blue")

		#plt.plot(X, E, color = "Purple")
		#plt.plot(X, DP-ChargeDensity)

		plt.figure("Potentials")
		E2 = E_Field(self.X, self.Phi)
		plt.plot(self.X, E2, color = "Blue")
		plt.plot(self.X, self.DopingDensity, color = "Grey")

		self.Phi[-1] = 0.0; self.Phi[-2] = 0.0
		Phi5 = solverLU_backward(self.X, self.Phi, self.ChargeDensity, self.Dieelek)
		plt.plot(self.X, Phi5, color = "Red")
		
		self.Phi[1] = self.Phi[0] = (u.BarrierHeight+self.Bias); self.Phi[-1] = 0.0
		self.Phi = solverLU(self.X, self.Phi, self.ChargeDensity, self.Dieelek)
		plt.plot(self.X, self.Phi, color = "Green")

		plt.show()
		
	def WriteToFile(self, File):
		savetxt(File, (self.Phi, self.DopingDensity, self.ChargeDensity), newline="\n")
	
	def Diffusion(self, DiffusionCoeffcient):
		self.DopingDensity = Diffusion(self.X, self.dt, self.DopingDensity, DiffusionCoeffcient)
