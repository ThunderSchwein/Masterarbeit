# Konstantin Murasov
# DopingDrift.py
""" 
This module uses the Information provided by "Schottky.py"
to calculate the migration of doping inside a TiO2-Sample.
"""

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from ChargeRedistribution import *

def Smoothing(Sample, s0 = 0.25):
	# As the Drift-RK method produces uneven doping profile
	# which can result in unnecessary spikes in electric field,
	# the resulting doping density must be smoothed.
	
	# Note that this smoothing-method has the form of a diffusion
	# equation and can be transformed into a diffusion-metod via
	# s = (Diffusion constant)*dt/dx^2
	
	DP = Sample.DopingDensity
	DP_neu = zeros(Sample.nx)
	# 2*s = part of the density moving to the neighbouring grid points
	s = s0
	
	# In the first and the last grid cell, the second derivative must be 
	# 
	DP_neu[0] = DP[0] +s*(DP[1]-DP[0])
	DP_neu[-1] = DP[-1] +s*(DP[-2]-DP[-1])
	
	for i in range(Sample.nx-2):
		DP_neu[i+1]  = DP[i+1] 		+ s*(DP[i+2] - 2*DP[i+1] + DP[i])
	
	Sample.DopingDensity = DP_neu
	return

#def RunningMean(Sample, n):
	DP  = Sample.DopingDensity
	#Pivots
	P = []
	#for i in range(int(Sample.nx/n)):
		#
	#for i in range(int(Sample.nx/n)):
		#for j in range(n):
		#DP[i+int(.5*n)] = 0
#	return
	
def DriftRK(Sample, DriftAmount_0 = 1 ):
# DriftRK-method calculates the drift of the doping atoms in a given sample,
# using the doping density and the electric field information saved in the
# arrays "Sample.DopingDensity" and "Sample.E"

# The variable DriftAmount_0 is 1 by default, so that the entirety of the 
# doping concentration in a grid cell is moved in every time step.
# Values <1 lead to a exponential drift effect, where the amount of the 
# dopnads move depends on the amount of the available dopands.
# Values >1 violate the condition that the amount doping atoms is always > 0. 

# Runge-Kutta cannot be implemented correctly, since then the electric field
# would have to be calculated after each drift process in every grid cell, 
# leading to greatly increased simulation time (factor Devicelength/dx).

	Phi = Sample.Phi	# Doping Ions Potential
	E = Sample.E		# Electric Field
	dx = Sample.dx		# Grid Cell Width
	nx = Sample.nx		# Number of grid cells
	dt = Sample.dt		# Duaration of a time step
	DP = zeros(nx)		# Intermediate change in the doping density
	DP0 = Sample.DopingDensity # Shortcut for the doping density of the sample

	# Runge-Kutta 4th Order Scheme
	for h in range(int(1/dt)):
	# To counteract the loss in precision due to the implementation of 
	# time-independent Runge-Kutta-Scheme, the change in the electric field
	# is calculated after a fraction of the doping ions, which corresponds
	# to (DriftAmount_0/dt), has been pushed.

		for i in range(0, nx):
		# for every grid cell, the amount and the width of displacement of the doping is calculated, ...
			if ((DP0[i] != 0) &  (abs(E[i]) > u.E_min)) :
			# if there are ions to be pushed and a sufficient electric field
			
				# Time-independent Runge-Kutta-Scheme
				k1 = E[i]*u.Titanium_DopingMobility*dt/dx

				# Catching ions before they escape the sample:
				if((i+k1/2) < 0):	 		k2 = E[0]*u.Titanium_DopingMobility*dt/dx
				elif ((i+k1/2) > (nx-1)): 	k2 = E[nx -1]*u.Titanium_DopingMobility*dt/dx
				# If they are not escaping :
				else : 						k2 = E[int(i + k1/2)]*u.Titanium_DopingMobility*dt/dx
				
				# Catching ions before they escape the sample:
				if(   (i+k2/2) < 0):	 	k3 = E[0]*u.Titanium_DopingMobility*dt/dx	
				elif ((i+k2/2) > (nx-1)):	k3 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
				# If they are not escaping :
				else : 						k3 = E[int(i + k2/2)]*u.Titanium_DopingMobility*dt/dx
				
				# Catching ions before they escape the sample:
				if(   (i+k3) < 0):	 		k4 = E[0]*u.Titanium_DopingMobility*dt/dx
				elif ((i+k3) > (nx-1)):		k4 = E[nx-1]*u.Titanium_DopingMobility*dt/dx
				# If they are not escaping :
				else : 						k4 = E[int(i + k3)]*u.Titanium_DopingMobility*dt/dx
		
				# k0 is the displacement in terms of change in the index
				k0 = (k1/6 + k2/3 + k3/3 + k4/6)
	
				# k = new index of a given moving Doping partition
				k = int(i+k0)
				# Doping must not be allowed to leave the Sample !
				if(k<1): k = 0
				if(k>(nx-2)): k = nx-2
				
				# Calculating the fraction of dopnads, which is moved per step	
				DriftAmount = DriftAmount_0*dt
				
				# Removing the doping from its initial Position	
				DP[i] = DP[i] - DP0[i]*DriftAmount
				# Depositing the Doping at the new Position
				for h in range(k+1):
					g = i + h
					DP[g] = DP[g] + DP0[i]*DriftAmount/k
				
				#DP[k]  = DP[k] + (DP0[i])*DriftAmount*(1-(k0%1)/2)
				#DP[k+1]  = DP[k+1] + (DP0[i])*DriftAmount*((k0%1)/2)
				
				# Setting the allocated memory free
				del k0, k1, k2, k3, k4, k
				
		# Calculating the resulting doping density		
		Sample.DopingDensity = DP0 + DP
		# Calcualating the corresponding electric field
		Reload(Sample)
		
		# Smoothing the data
		#Smoothing(Sample)
		#Smoothing(Sample)
		#Smoothing(Sample)
		
	return
