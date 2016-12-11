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
			
			k = i+int(-E[i]*u.Titanium_DopingMobility/dx*dt)
			#f = k%1-0.5
			
			
			DP[k-1] = DP[k-1] + DopingDensity[i]/8
			DP[k]   = DP[k] + DopingDensity[i]/4
			DP[k+1] = DP[k+1] + DopingDensity[i]/8
			
			# include Gaussian Statistics
	return (DP+DopingDensity)

# Runge-Kutta 4th Order Scheme
	
def DriftRK( X , Potential, DopingDensity, dt, DriftAmount = 0.2):
	Phi = Potential
	DP = zeros(len(DopingDensity))
	
	E = E_Field(X, Phi)
	dx = X[1]-X[0]
	
	for i in range(1,len(Phi)):
		if (DopingDensity[i] != 0) & (abs(E[i]) > 4e-4) :
			k1 = E[i]*u.Titanium_DopingMobility*dt
			k2 = E[int(i- k1/2)]*u.Titanium_DopingMobility*dt
			k3 = E[int(i- k2/2)]*u.Titanium_DopingMobility*dt
			k4 = E[int(i- k3)  ]*u.Titanium_DopingMobility*dt
			
			k0 = (k1/6 + k2/3 + k3/3 + k4/6)
			k = int(i-k0)
			if(k<10): k = 10
			
			#if (i < 1000) :	print(i, k0)
			DP[i] = DP[i] - DopingDensity[i]*DriftAmount
			
			DP[k]  = DP[k] + DopingDensity[i]*DriftAmount*(1-k0%1)
			DP[k-sign(k0)]  = DP[k-sign(k0)] + DopingDensity[i]*DriftAmount*(k0%1)
			
			#DP[k-1] = DP[k-1] + DopingDensity[i]/20
			#DP[k]   = DP[k] + DopingDensity[i]/10
			#DP[k+1] = DP[k+1] + DopingDensity[i]/20
			
	return DP+DopingDensity
	

"""
	for i in range(1, len(DP)):
		if (DP[i] != 0) & (abs(E[i]) > 4e-4) :
			k = int(-E[i]*u.Titanium_mob/dx*dt)
			DP[i+k] = DP[i+k] + DopingDensity[i]/2
			DP[i]=0
	return DP
"""
