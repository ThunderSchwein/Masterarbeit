from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *

def Drift ( X , Potential, DopingDensity, dt):
	Phi = Potential
	#Rho = ChargeDensity
	#Di = Dielectricity1D
	DP = zeros(len(DopingDensity))
	
	E = E_Field(X, Phi)
	dx = X[1]-X[0]
	
	for i in range(1, len(DP)):
		if (DopingDensity[i] != 0) & (abs(E[i]) > 4e-4) :
			DP[i] = DP[i] - DopingDensity[i]/2
			
			
			#DP[i+int(sign(E[i]))] = DP[i+int(sign(E[i]))] + DP[i]/2
			#DP[i+10] = DP[i+10] + DopingDensity[i]/2
			
			#k=1
			
			k = int(-E[i]*u.Titanium_DopingMobility/dx*dt)
			#f = k%1-0.5
			
			
			DP[i+k-1] = DP[i+k-1] + DopingDensity[i]/8
			DP[i+k]   = DP[i+k] + DopingDensity[i]/4
			DP[i+k+1] = DP[i+k+1] + DopingDensity[i]/8
			
			# include Gaussian Statistics
	return (DP+DopingDensity)

def DriftPoisson( X , Potential, DopingDensity, dt):
	Phi = Potential
	#Rho = ChargeDensity
	#Di = Dielectricity1D
	DP = zeros(len(DopingDensity))
	
	E = E_Field(X, Phi)
	dx = X[1]-X[0]
	
	for i in range(1, len(DP)):
		if (DopingDensity[i] != 0) & (abs(E[i]) > 4e-4) :
			DP[i] = DP[i] - DopingDensity[i]/2
			
			
			#DP[i+int(sign(E[i]))] = DP[i+int(sign(E[i]))] + DP[i]/2
			#DP[i+10] = DP[i+10] + DopingDensity[i]/2
			
			#k=1
			
			k = -E[i]*u.Titanium_DopingMobility/dx*dt			
			
			DP[i+k-1] = DP[i+k-1] + DopingDensity[i]/8
			DP[i+k]   = DP[i+k] + DopingDensity[i]/4
			DP[i+k+1] = DP[i+k+1] + DopingDensity[i]/8
			
			# include Gaussian Statistics
	return (DP+DopingDensity)
	

"""
	for i in range(1, len(DP)):
		if (DP[i] != 0) & (abs(E[i]) > 4e-4) :
			k = int(-E[i]*u.Titanium_mob/dx*dt)
			DP[i+k] = DP[i+k] + DopingDensity[i]/2
			DP[i]=0
	return DP
"""
