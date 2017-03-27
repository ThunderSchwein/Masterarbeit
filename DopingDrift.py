""" 21.01.2017 Konstantin Murasov
This module uses the Information provided by "Schottky.py"
to calculate the migration of doping inside a TiO2-Sample.

Please take note that the Diffusion-Method is not correct.

"""

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *

def DriftRK(Sample, DriftAmount = 1.0):
	Phi = Sample.Phi
	E = Sample.E
	dx = Sample.dx
	nx = Sample.nx
	dt = Sample.dt
	
	DP = zeros(nx)
	
	# Runge-Kutta 4th Order Scheme
	for i in range(0,len(Phi)):
		if ((Sample.DopingDensity[i] != 0) &  (abs(E[i]) != 0)) :
			k1 = E[i]*u.Titanium_DopingMobility*dt/dx

			if((i+k1/2) < 0):	 		k2 = E[0]*u.Titanium_DopingMobility*dt/dx
			elif ((i+k1/2) > (nx-1)): 	k2 = E[nx -1]*u.Titanium_DopingMobility*dt/dx
			else : 						k2 = E[int(i + k1/2)]*u.Titanium_DopingMobility*dt/dx
			
			if(   (i+k2/2) < 0):	 	k3 = E[0]*u.Titanium_DopingMobility*dt/dx	
			elif ((i+k2/2) > (nx-1)):	k3 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
			else : 						k3 = E[int(i + k2/2)]*u.Titanium_DopingMobility*dt/dx
			
			if(   (i+k3) < 0):	 		k4 = E[0]*u.Titanium_DopingMobility*dt/dx
			elif ((i+k3) > (nx-1)):		k4 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
			else : 						k4 = E[int(i + k3)  ]*u.Titanium_DopingMobility*dt/dx
	
			k0 = (k1/6 + k2/3 + k3/3 + k4/6)
			#print(E[i],k0, sign(k0))
	# k0 is the change in the index		
	# k = new index of a given moving Doping partition
			k = int(i+k0)
	# Doping must not be allowed to leave the Sample
			if(k<1): k = 1
			if(k>(nx-2)): k = nx-2
	# Calculation the DriftAmount via Super-Duper-Formula
			DriftAmount =  ((1-exp(-(E[i]/u.E_min)**2))*dt)
			if(DriftAmount > 1) : DriftAmount = 1
	# Removing the doping from its initial Position	
			DP[i] = DP[i] - Sample.DopingDensity[i]*DriftAmount
	# Depositing the Doping at the new Position
			DP[k]  = DP[k] + Sample.DopingDensity[i]*DriftAmount*(1-k0%1)
			DP[k+1]  = DP[k+1] + Sample.DopingDensity[i]*DriftAmount*(k0%1)				
			del k0, k1, k2, k3, k4, k, DriftAmount
	#del dx, nx, E
	Sample.DopingDensity = DP+Sample.DopingDensity
	return
	
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
