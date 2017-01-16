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


	
def DriftRK( X , Potential, DopingDensity, dt, DriftAmount = 1.0):
	Phi = Potential
	DP = zeros(len(DopingDensity))
	
	E = E_Field(X, Phi)
	dx = X[1]-X[0]
	
	# Runge-Kutta 4th Order Scheme
	for i in range(1,len(Phi)):
		if (DopingDensity[i] != 0) & (abs(E[i]) > 4e-4) :
			k1 = E[i]*u.Titanium_DopingMobility*dt
			k2 = E[int(i- k1/2)]*u.Titanium_DopingMobility*dt
			k3 = E[int(i- k2/2)]*u.Titanium_DopingMobility*dt
			k4 = E[int(i- k3)  ]*u.Titanium_DopingMobility*dt
			
			k0 = (k1/6 + k2/3 + k3/3 + k4/6)
	# k0 is the change in the index		
	# k = new index of a given moving Doping partition
			k = int(i-k0)
	# Doping must not be allowed to leave the Sample
			if(k<10): k = k%10 + 2
			if(k>len(X)): k = len(X) - k*10 - 2
	
	# Removing the doping from its initial Position
			DP[i] = DP[i] - DopingDensity[i]*DriftAmount
	# Depositing the Doping at the new Position
			DP[k]  = DP[k] + DopingDensity[i]*DriftAmount*(1-k0%1)
			DP[k-sign(k0)]  = DP[k-sign(k0)] + DopingDensity[i]*DriftAmount*(k0%1)					
	return DP+DopingDensity
	
def Diffusion(X, dt, DopingDensity, DiffusionCoeffcient):
	DP = DopingDensity
	DP_neu = zeros(len(X))
	DP_final = zeros(len(X))
	#s = DiffusionCoeffcient*dt/(X[1]-X[0])**2
	s = 0.4
	
	#DP_neu[0] = DP[0]
	#DP_neu[-1] = DP[-1]
	DP_neu[0] = DP[0] +s*(DP[1]-DP[0])
	for i in range(len(X)-2):
		#DP_neu[i+1] = s*(DP[i+2] + DP[i])
		#DP_neu[i]   = DP_neu[i]	- s*(DP[i+2] - 2*DP[i+1] + DP[i]) #/2
		DP_neu[i+1]  = DP[i+1] 		+ s*(DP[i+2] - 2*DP[i+1] + DP[i]) # + DP_neu[i+1] 
		#DP_neu[i+2] = DP_neu[i+2]	- s*(DP[i+2] - 2*DP[i+1] + DP[i])/2
	DP_neu[-1] = DP[-1] +s*(DP[-1]-DP[-2])
	
	
	DP_final[0]	 = (DP_neu[0]	  + s*( DP_neu[0]+DP_neu[1] ))/(1+2*s)
	for i in range(len(X)-2):
		DP_final[i+1]  = (DP[i+1] + s*(DP_neu[i]+DP_neu[i+2]))/(1+2*s)
	DP_final[-1] = (DP_neu[-1]	  + s*(DP_neu[-2]+DP_neu[-1]))/(1+2*s)
	
	DP_final[0] = DP_final[0] + sum(DP) - sum(DP_final)
	
	#for i in range(len(X)-4):
		#DP_neu[i+2] = (1-s)*DP[i+2] - s*(DP[i]-2*DP[i+1]-2*DP[i+3]+DP[i+4])/2
		
	return DP_final
	
"""DP = DopingDensity
	DP0 = DopingDensity[0:2]
	
	s = DiffusionCoeffcient*dt*(X[1]-X[0])
	
	DP[1] = DP[1] + s*(DP0[2] - 2*DP0[1] + DP0[0])
	for i in range(len(X)-2) + 3 :
		DP0 = [DP0[1], DP0[2], DP[i]]
		DP[i] = - (DP0[0] + 2/s*(1-s)*DP0[1] + DP0[2]) + (2/s*(1+s)*DP[i-1] - DP[i-2] )
	
	return DP	
"""
"""
	for i in range(1, len(DP)):
		if (DP[i] != 0) & (abs(E[i]) > 4e-4) :
			k = int(-E[i]*u.Titanium_mob/dx*dt)
			DP[i+k] = DP[i+k] + DopingDensity[i]/2
			DP[i]=0
	return DP
"""
