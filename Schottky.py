# KonstantinMurasov
# Schottky.py
""" 
This Module utilizes the Modules "DopingDrift.py" and "ChargeRedistribution.py"
to calculate the changes in doping density and depletio length due to migration
of doping atoms in a suffficient electric field, provided in "Units.py"
"""
# The data, necessary for calculation, are saved in the Schottky-object.

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import pandas

import Units as u
from SolverLU import *
from DopingDrift import *
from ChargeRedistribution import *


class Schottky(object) :
	def __init__(self, Bias, dt, nx): 
		#------------------------------------
		# 0. Initialize 
			# X-values
		self.X = linspace(0, u.DeviceLength, nx)
			# Grid unit width /nm
		self.dx = self.X[1] - self.X[0]
			# Number of grid units
		self.nx = nx
			# Timestep /s
		self.dt = dt
			# Voltage applied to the metal side of the Schottky-junction
		self.Bias = Bias
		
		# 1. Calculate Initial Conditions
			# Dielectricity
		self.Dieelek = ones(nx)*u.e0*u.Titanium_DielectricityFactor
				
			# Initial (Uniform) Doping Density
		self.DopingDensity = zeros(self.nx)
		self.DopingDensity[:] = u.Titanium_DopingDensity
		
			# Initial Hydrogen Concentration
		self.Hydrogen = zeros(nx)
		self.H2Concentration = 0

		# Initail Charge density & Depletion Region Width
		self.ChargeDensity = zeros(nx)
		# Case of no depletion region
		# Never put depletion region width to 0 !
		if(self.Bias > -u.BarrierHeight):
			self.W = 1
			self.ChargeDensity = 0
		else: 
			DepletionRegionCharge = u.Titanium_DopingCharge
			self.W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight+self.Bias)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5
			self.ChargeDensity[:int(self.W/u.DeviceLength*nx)] = u.Titanium_DopingCharge*u.Titanium_DopingDensity
		print(sum(self.ChargeDensity))
		
		# 2. Initial Potential Calculation
			# Doping ions Potential
		self.Phi = zeros(nx)
		solverLU(self)
		
			#Electrical Field
		self.E = zeros(nx)
		self.E = E_Field(self)

		return
	#------------------------------------		
	def PlotResults(self):
	# Method for visualizing the state of the Schottky-Diode-Sample
	
		plt.figure("Charge & Doping Distribution")
		plt.xlabel("X [nm]")
		plt.ylabel("Concentration [1/nm^3]")
			# Doping Density plot (default color = black)
		P1, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Doping Density")
			# Charge Density plot (default color = red)
		P2, = plt.plot(self.X, self.ChargeDensity, color = "Red" , label = "Charge Density")
			# Hydrogen concentration plot (default color = blue)
		#P3, = plt.plot(self.X, self.Hydrogen	 , color = "Blue", label = "Incorporated Hydrogen Density")
		plt.legend()
		
		plt.figure("Potentials")
		plt.xlabel("X [nm]")
		plt.ylabel("Potential / Electrical Field / Doping Concentration")
			# Recalculating the potentials, before plotting them.
		solverLU(self)
		self.E = E_Field(self)
			# Doping ions potential plot (default color = green)
		P4, = plt.plot(self.X, self.Phi, color = "Green", label = "Potential [V]")
			# Doping ions potential plot (default color = blue)
		P5, = plt.plot(self.X, self.E  , color = "Blue" , label = "Electrical Field [V/nm]")
			# Plot of the doping density (deault color = black)
		P6, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Doping Density [1/nm^3]")
		plt.legend(loc = 4)
		
		plt.figure("Doping Density")
			# A presentation of the doping density with a higher line width (lw = 3)
		plt.xlabel("X [nm]")
		#df = pandas.DataFrame(self.DopingDensity)
		#df = pandas.DataFrame.rolling(df ,window = 20, center = False).mean()
		#P6, = plt.plot(self.X, df, color = "Grey", label = "Oxygen Vacancy Density [1/nm^3]", lw = 3)
		
		#plt.ylabel("Doping Concentration")
		
		#modes = ['valid']#'full', 'same' , 
		#D = zeros(self.nx)
		#for i in range(self.nx) : D[i] =  self.DopingDensity[-i]
		#for m in modes:
			#plt.plot(convolve(self.DopingDensity, D, mode=m));
				
			# Plotof the doping density (deault color = black)
		P6, = plt.plot(self.X, self.DopingDensity, color = "Grey", label = "Oxygen Vacancy Density [1/nm^3]", lw = 3)
		plt.legend()
		
		plt.show()
		
		return
		
	def DopingDrift(self, DriftAmount):
		# Calculation of Doping Drift
		# Drift proccess
		DriftRK(self, DriftAmount)	
		# Calculation of the corresponding charge density and depletion region width
		Reload(self)			
		# Calculation of the new doping ion potential
		solverLU(self)
		return
		
	def WriteToFile(self, File):
		# Saves the state of the Schottky Diode to an external file
		savetxt(File, (self.Phi, self.DopingDensity, self.ChargeDensity), newline="\n")
		return
	
	def Import(self, File):
		# Import Data of a Sample from a File.
		# The length of Data arrays must match.
		
		D =  loadtxt(File)
		self.Phi			= D[0]
		self.DopingDensity	= D[1]
		self.ChargeDensity	= D[2]
		self.Bias = self.Phi[0] - u.BarrierHeight 
		self.W = len(D[2][D[2] != 0])*self.dx
		self.E = E_Field(self)
		print("Import of Sample Data  from File '", File, "' has been successful.")
		return
	
	
	#def Diffusion(self, DiffusionCoeffcient):
		#self.DopingDensity = Diffusion(self.X, self.dt, self.DopingDensity, DiffusionCoeffcient)
		#return
