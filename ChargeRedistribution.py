""" 31.01.2017 Konstantin Murasov
This Module calculates the Distribution of charges, which
results from doping migration in a TiO2-Sample, with Data
provided by "Schottky.py"

Please take note that "ReloadDeltaPhi"-Method utilizes a
try-and-error approach.
"""

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *


def Reload(Sample):
	ReloadRK(Sample)
	return

"""
def Reload(X, ChargeDensity, DopingDensity, FermiLevel = 0):
	S = -sum(ChargeDensity)
	#print(S)
	i = 0
	D = DopingDensity
	while(S > D[i]):
		ChargeDensity[i] =  -2*D[i]
		S = S-2*D[i]
		i = i+1
	ChargeDensity[i] = 2*D[i] - S
	ChargeDensity[i+1 :] = 0
	
	return ChargeDensity

def ReloadFermi(Sample):
	#Sample.ChargeDensity = (1-exp(Sample.Phi/(u.Kb*u.T)))*Sample.DopingDensity*u.Titanium_DopingCharge
	
	#print("Redistributing the Electrons")
	nx = Sample.nx
	
	#Sample.ChargeDensity = zeros(nx)
	
	if(Sample.Bias < -u.BarrierHeight):		
		# Inital Guess
		#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
		#solverLU(Sample)
		if(Sample.W > u.DeviceLength):
			Sample.ChargeDensity = Sample.DopingDensity*u.Titanium_DopingCharge
			Sample.W = u.DeviceLength
		else:
			factor1 = 1.05
			factor2 = 1.01
			factor3 = 1.002
			factor4 = 1.0004
			# Improvement of Initial Guess
			if ((Sample.Phi[0]>(u.BarrierHeight+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				Sample.W = Sample.W*factor1
				#Sample.ChargeDensity[:] = 0.0
				#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				Sample.ChargeDensity = (exp(Sample.Phi))*Sample.DopingDensity*u.Titanium_DopingCharge
				solverLU(Sample)
				print(1, Sample.W, Sample.Phi[0])

			while (Sample.Phi[0]<(u.BarrierHeight+Sample.Bias)):
				Sample.W = Sample.W / factor2
				Sample.ChargeDensity[:] = 0.0
				#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				Sample.ChargeDensity = (1-exp(Sample.Phi/(u.Kb*u.T)))*Sample.DopingDensity*u.Titanium_DopingCharge
				solverLU(Sample)
				print(2, Sample.W, Sample.Phi[0])
				
			while ((Sample.Phi[0]>(u.BarrierHeight*0.99875+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				Sample.W = Sample.W*factor3
				#Sample.ChargeDensity[:] = 0.0
				#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				Sample.ChargeDensity = (1-exp(Sample.Phi/(u.Kb*u.T)))*Sample.DopingDensity*u.Titanium_DopingCharge
				solverLU(Sample)
				print(3, Sample.W, Sample.Phi[0])
			while (Sample.Phi[0]<(u.BarrierHeight*1.00125+Sample.Bias)):
				Sample.W = Sample.W /factor4
				Sample.ChargeDensity[:] = 0
				#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				Sample.ChargeDensity = (1-exp(Sample.Phi/(u.Kb*u.T)))*Sample.DopingDensity*u.Titanium_DopingCharge
				solverLU(Sample)
				print(4, Sample.W, Sample.Phi[0])
			#del DepletionRegionCharge
	else:
		Sample.W = 0#u.Titanium_FreePathLength
		Sample.Phi[0] = Sample.Bias + u.BarrierHeight
		Sample.ChargeDensity = zeros(nx)
		#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
		solverLU(Sample)

	#print(0, Sample.W, Sample.Phi[0])
	Sample.E = E_Field(Sample.X, Sample.Phi, Sample.E)
	return
	"""
	
def ReloadRK(Sample) :
	#print("Redistributing the Electrons")
	nx = Sample.nx
	
	Sample.ChargeDensity = zeros(nx)
	
	if(Sample.Bias < -u.BarrierHeight):
		factor1 = 1.05
		factor2 = 1.01
		factor3 = 1.002
		factor4 = 1.0004

		if((Sample.Bias > -u.BarrierHeight*0.99) & (Sample.Bias < -u.BarrierHeight*1.01)) :
					Sample.W = 1#u.Titanium_FreePathLength
					Sample.Phi = zeros(nx)
					Sample.ChargeDensity = zeros(nx)
		
		#Inital Guess
		#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
		#solverLU(Sample)
		
		if(Sample.W >= u.DeviceLength):
			Sample.ChargeDensity = Sample.DopingDensity*u.Titanium_DopingCharge
			Sample.W = u.DeviceLength
			solverLU(Sample)
		else:
			# Improvement of Initial Guess
			while ((Sample.Phi[0]>(u.BarrierHeight+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				Sample.W = Sample.W*factor1
				#Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				solverLU(Sample)
				#print(1, Sample.W, Sample.Phi[0])
			while (Sample.Phi[0]<(u.BarrierHeight+Sample.Bias)):
				Sample.W = Sample.W / factor2
				Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				solverLU(Sample)
				#print(2, Sample.W, Sample.Phi[0])
				
			while ((Sample.Phi[0]>(u.BarrierHeight*0.99875+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				Sample.W = Sample.W*factor3
				#Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				solverLU(Sample)
				#print(3, Sample.W, Sample.Phi[0])
			while (Sample.Phi[0]<(u.BarrierHeight*1.00125+Sample.Bias)):
				Sample.W = Sample.W /factor4
				Sample.ChargeDensity[:] = 0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				solverLU(Sample)
				#print(4, Sample.W, Sample.Phi[0])
			#del DepletionRegionCharge
	
	
	else:
		Sample.W = 1#u.Titanium_FreePathLength
		Sample.Phi[0] = Sample.Bias + u.BarrierHeight
		Sample.ChargeDensity = zeros(nx)
		#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
	solverLU(Sample)

	#print(0, Sample.W, Sample.Phi[0])
	Sample.E = E_Field(Sample.X, Sample.Phi, Sample.E)

	return
