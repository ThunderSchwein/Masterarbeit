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
from ChargeRedistribution import *

def Smoothing(Sample):
	DP = Sample.DopingDensity
	DP_neu = zeros(Sample.nx)
	DP_final = zeros(Sample.nx)
	s = 0.25
	
	DP_neu[0] = DP[0] +s*(DP[1]-DP[0])
	for i in range(Sample.nx-2):
		DP_neu[i+1]  = DP[i+1] 		+ s*(DP[i+2] - 2*DP[i+1] + DP[i])
	DP_neu[-1] = DP[-1] +s*(DP[-2]-DP[-1])

	Sample.DopingDensity = DP_neu
	return

def RunningMean(Sample, n):
	DP  = Sample.DopingDensity
	#Pivots
	P = []
	#for i in range(int(Sample.nx/n)):
		#
	#for i in range(int(Sample.nx/n)):
		#for j in range(n):
		#DP[i+int(.5*n)] = 0
	return
	
def DriftRK(Sample, DriftAmount_0 ):
	Phi = Sample.Phi
	E = Sample.E
	dx = Sample.dx
	nx = Sample.nx
	dt = Sample.dt
	DP = zeros(nx)
	DP0 = Sample.DopingDensity
	print("DriftRK")
	# Runge-Kutta 4th Order Scheme
	for h in range(int(1/dt)):
		for i in range(0,len(Phi)):
			if ((DP0[i] != 0) &  (abs(E[i]) > u.E_min)) :
				k1 = E[i]*u.Titanium_DopingMobility*dt/dx

				if((i+k1/2) < 0):	 		k2 = E[0]*u.Titanium_DopingMobility*dt/dx
				elif ((i+k1/2) > (nx-1)): 	k2 = E[nx -1]*u.Titanium_DopingMobility*dt/dx
				else : 						k2 = E[int(i + k1/2)]*u.Titanium_DopingMobility*dt/dx
				
				if(   (i+k2/2) < 0):	 	k3 = E[0]*u.Titanium_DopingMobility*dt/dx	
				elif ((i+k2/2) > (nx-1)):	k3 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
				else : 						k3 = E[int(i + k2/2)]*u.Titanium_DopingMobility*dt/dx
				
				if(   (i+k3) < 0):	 		k4 = E[0]*u.Titanium_DopingMobility*dt/dx
				elif ((i+k3) > (nx-1)):		k4 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
				else : 						k4 = E[int(i + k3)]*u.Titanium_DopingMobility*dt/dx
		
				k0 = (k1/6 + k2/3 + k3/3 + k4/6)
				#print(E[i],k0, sign(k0))
		# k0 is the displacement in terms of change in the index		
		# k = new index of a given moving Doping partition
				k = int(i+k0)
				#if(k0 != 0): print(k0)
		# Doping must not be allowed to leave the Sample
				if(k<1): k = 1
				if(k>(nx-2)): k = nx-2
		# Removing the doping from its initial Position	
				DriftAmount = DriftAmount_0*dt
				
				#DP[i] = DP[i] - (Sample.DopingDensity[i]+DP[i])*DriftAmount
				DP[i] = DP[i] - DP0[i]*DriftAmount
		# Depositing the Doping at the new Position
				DP[k]  = DP[k] + (DP0[i])*DriftAmount*(1-k0%1)
				DP[k+1]  = DP[k+1] + (DP0[i])*DriftAmount*(k0%1)
				del k0, k1, k2, k3, k4, k	
		Sample.DopingDensity = DP0 + DP
		Reload(Sample)
		#print("R", h, Sample.E[0], Sample.DopingDensity[0])
		Smoothing(Sample)
		Smoothing(Sample)
		Smoothing(Sample)
		
	return
