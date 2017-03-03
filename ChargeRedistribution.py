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
	
def ReloadFermi(X, DopingDensity, Potential, ElectricField ,FermiLevel = 0):
	ChargeDensity = zeros(len(X))
	#ChargeDensity[1:] = (1/exp(-Potential[:-1]/u.Kb)) - 1 )*DopingDensity[1:]*u.Titanium_DopingCharge
	ChargeDensity[:-1] = (2/(1+exp(Potential[1:]/(u.Kb*0.001))) - 1)*DopingDensity[:-1]*u.Titanium_DopingCharge
	ChargeDensity[abs(ElectricField) <  2*abs(ElectricField[-1])] = 0
	
	return ChargeDensity

def ReloadDeltaPhi(Sample) :
	#print("Redistributing the Electrons")
	nx = Sample.nx
	Phi2 = Sample.Phi
	
	Sample.ChargeDensity = zeros(nx)
	
	if(Sample.Bias < abs(u.BarrierHeight)):
		DepletionRegionCharge = u.Titanium_DopingCharge
		factor1 = 1.003
		factor2 = 1.0015
		factor3 = 1.00075
		factor4 = 1.00015
		
		# Inital Guess
		#Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
		solverLU_backward(Sample)
		
		# Improvement of Initial Guess
		while (Phi2[0]>(u.BarrierHeight*0.995+Sample.Bias)):
			Sample.W = Sample.W*factor1
			Sample.ChargeDensity[:] = 0
			Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
			solverLU_backward(Sample)
			#print(1, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.005+Sample.Bias)):
			Sample.W = Sample.W / factor2
			Sample.ChargeDensity[:] = 0
			Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
			solverLU_backward(Sample)
			#print(2, W, Phi2[0])
			
		while (Phi2[0]>(u.BarrierHeight*0.999+Sample.Bias)):
			Sample.W = Sample.W*factor3
			Sample.ChargeDensity[:] = 0
			Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
			solverLU_backward(Sample)
			#print(3, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.001+Sample.Bias)):
			Sample.W = Sample.W /factor4
			Sample.ChargeDensity[:] = 0
			Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*DepletionRegionCharge
			solverLU_backward(Sample)
			#print(4, W, Phi2[0])
		del DepletionRegionCharge
	else:
		Sample.W = u.Titanium_FreePathLength
		Potential[0] = Sample.Bias - u.BarrierHeight
		Sample.ChargeDensity = zeros(nx)
		
	#print(0, W, Phi2[0])
	
	Sample.E = E_Field(Sample.X, Sample.Phi, Sample.E)
	return
