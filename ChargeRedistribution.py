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

def ReloadDeltaPhi(X, DopingDensity, ChargeDensity, Potential, Dieelek, W, FermiLevel = 0) :
	#print("Redistributing the Electrons")
	nx = len(X)
	# Detemining the ruling Potential
	Phi2 = zeros(len(X))
	Phi2[-1] = Phi2[-2] = 0
	Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
	
	ChargeDensity[:] = 0
	
	if(FermiLevel > -u.BarrierHeight):
		DepletionRegionCharge = u.Titanium_InversionCharge
		factor1 = 1/1.003
		factor2 = 1/1.0015
		factor3 = 1/1.00075
		factor4 = 1/1.00015
		
		ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
		Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
		
		while (Phi2[0]>(u.BarrierHeight*0.995+FermiLevel)):
			W = W*factor1
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = u.Titanium_InversionDensity*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(1, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.01+FermiLevel)):
			W = W / factor2
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = u.Titanium_InversionDensity*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(2, W, Phi2[0])
			
		while (Phi2[0]>(u.BarrierHeight*0.995+FermiLevel)):
			W = W*factor3
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = u.Titanium_InversionDensity*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(3, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.01+FermiLevel)):
			W = W /factor4
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = u.Titanium_InversionDensity*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(4, W, Phi2[0])		
	else: 
		DepletionRegionCharge = u.Titanium_DopingCharge
		factor1 = 1.003
		factor2 = 1.0015
		factor3 = 1.00075
		factor4 = 1.00015
		
		ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
		Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
		
		while (Phi2[0]>(u.BarrierHeight*0.995+FermiLevel)):
			W = W*factor1
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(1, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.01+FermiLevel)):
			W = W / factor2
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(2, W, Phi2[0])
			
		while (Phi2[0]>(u.BarrierHeight*0.995+FermiLevel)):
			W = W*factor3
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(3, W, Phi2[0])
		while (Phi2[0]<(u.BarrierHeight*1.01+FermiLevel)):
			W = W /factor4
			ChargeDensity[:] = 0
			ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*DepletionRegionCharge
			Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
			#print(4, W, Phi2[0])		
		
	#print(0, W, Phi2[0])
	
	return [ChargeDensity, W]
	
	"""
	W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight*2-Phi2[0]2)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5	
	ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
	Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
	print(1, W, Phi2[0], u.BarrierHeight*3-Phi2[0]*2)
	W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight*3-Phi2[0]*2)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5	
	ChargeDensity[:int(W/u.DeviceLength*nx)] = DopingDensity[:int(W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
	Phi2 = solverLU_backward(X, Phi2, ChargeDensity, Dieelek)
	#ChargeDensity[:int((W+34)/u.DeviceLength*nx)] = DopingDensity[:int((W+34)/u.DeviceLength*nx)]*u.Titanium_DopingCharge
	print(2, W, Phi2[0]) #,int(W/u.DeviceLength*nx), u.BarrierHeight*2-Phi2[0], solverLU_backward(X, Phi2, ChargeDensity, Dieelek)[0])	
	"""
	

"""
#-----------------------------------
# Depletion Region Width Correction

plt.figure()

plt.plot(X, Phi, color = "red")
DP = ones(nx)*u.Titanium_DopingDensity
#DP[:int(0.2*nx)] = 1
plt.plot(X, Phi, color= "Black")
Rho2 = zeros(nx)
W2 = W/0.962
n3 = zeros(nx)
#Rho2[:int(W2/u.DeviceLength*nx)] = -u.Titanium_DopingCharge*u.Titanium_DopingDensity
#Phi[0] = u.BarrierHeight; Phi[-1] = 0.0
n0 = (1-exp(-Phi/(u.Kb*300)))*u.Titanium_DopingDensity*(-u.Titanium_DopingCharge)
#n3[1:] = (1-exp(-Phi[:-1]/(u.Kb*1e-10)))*u.Titanium_DopingDensity*(-u.Titanium_DopingCharge)
n3 = ReloadFermi(X, DP, Phi)
Phi2 = solverLU_backward(X, Phi, n0, Dieelek)
#n0[DP == 0] = 0
plt.plot(X, n0, color = "Purple")
plt.plot(X, n3, color = "blue" )
plt.plot(X, Rho, color = "green")

#Phi2 = solverLU_backward(X, Phi, n0, Dieelek)
#plt.plot(X, Phi2)
"""
