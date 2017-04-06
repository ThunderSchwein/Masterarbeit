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
				
		# Initial Doping Density
		self.DopingDensity = zeros(self.nx)
		self.DopingDensity[:] = u.Titanium_DopingDensity

		# Initail Charge density & Depletion Region Width
		self.ChargeDensity = zeros(nx)
		if(self.Bias > -u.BarrierHeight):
			self.W = 5
			self.ChargeDensity = 0
		else: 
			DepletionRegionCharge = u.Titanium_DopingCharge
			self.W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight+self.Bias)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5
			self.ChargeDensity[:int(self.W/u.DeviceLength*nx)] = u.Titanium_DopingCharge*u.Titanium_DopingDensity

		#Initial Potential Calculation
		self.Phi = zeros(nx)
		self.Phi[1] = self.Phi[0] = u.BarrierHeight+self.Bias; self.Phi[-1] = 0
		solverLU(self)
		
		#Electrical Field
		self.E = zeros(nx)
		self.E = E_Field(self.X, self.Phi, self.E)
		
		# Hydrogen Concentration
		self.Hydrogen = zeros(nx)
		self.H2Concentration = 0
		
		# Inital conditions Plot
		#plt.figure("Initial Conditions")
		#plt.plot(self.X, self.ChargeDensity)
		#plt.plot(self.X, self.Phi,color = "G")
		#plt.plot(self.X, self.E, color = "R")
		return
	#------------------------------------
	def DopingDrift(self):
		# 2. Calculate Doping Drift
		DriftRK(self)
		ReloadDeltaPhi(self)
		self.Phi[-1] = 0.0; self.Phi[-2] = 0.0
		solverLU(self)
		return
		
	def PlotResults(self):
		plt.figure("Charge & Doping Distribution")
		
		plt.xlabel("X [nm]")
		plt.ylabel("Concentration [1/nm^3]")
		
		P1, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Doping Density")
		P2, = plt.plot(self.X, self.ChargeDensity, color = "Red" , label = "Charge Density")
		P3, = plt.plot(self.X, self.Hydrogen	 , color = "Blue", label = "Incorporated Hydrogen Density")
		plt.legend()
		
		plt.figure("Potentials")
		
		plt.xlabel("X [nm]")
		plt.ylabel("Potential / Electrical Field / Doping Concentration")
		self.Phi[1] = self.Phi[0] = (u.BarrierHeight+self.Bias); self.Phi[-1] = 0.0
		solverLU(self)
		
		P4, = plt.plot(self.X, self.Phi, color = "Green", label = "Potential [V]")
		P5, = plt.plot(self.X, self.E  , color = "Blue" , label = "Electrical Field [V/nm]")
		
		P6, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Doping Density [1/nm^3]")
		plt.legend()
		
		plt.figure("Doping Density")
		
		plt.xlabel("X [nm]")
		plt.ylabel("Doping Concentration")
		P6, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Oxygen Vacancy Density [1/nm^3]")
		P7, = plt.plot(self.X, TotalResistance(self)/self.dx, color = "Green", label = "Resistance [Ohm/nm]")
		plt.legend()
		
		plt.show()
		
		return
		
	def WriteToFile(self, File):
		savetxt(File, (self.Phi, self.DopingDensity, self.ChargeDensity), newline="\n")
		return
	
	def Diffusion(self, DiffusionCoeffcient):
		self.DopingDensity = Diffusion(self.X, self.dt, self.DopingDensity, DiffusionCoeffcient)
		return
