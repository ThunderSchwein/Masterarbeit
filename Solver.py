from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u

def E_Field(X, Phi, E):
	dx = X[1]-X[0]
	E = deepcopy(X)
	for i in range(len(X)-1) :
		E[i+1] = (Phi[i] - Phi[i+1])/dx
	E[0] = E[1]
	
	return E
	
def solverLU(Sample):
	dx = Sample.dx
	nx = Sample.nx
	Phi = Sample.Phi
	Rho = Sample.ChargeDensity
	
	#Phi_ende = deepcopy(Sample.Phi[-1])
	# Solution based only on present unbalanced Charge 
	Phi[-1] = Phi[-2] = 0	
	for j in range(nx-2):
		Phi[nx-3-j] = 2*Phi[nx-j-2] - Phi[nx-j-1] - Rho[nx-j-2]*dx**2/Sample.Dieelek[nx-j-2]
	
	if(Sample.Bias > abs(u.BarrierHeight)):
		#Phi_ende = deepcopy(Sample.Phi[-1])		
		
		#for j in range(len(Phi)-2):		
			#Phi[j+2] = 2*Phi[j+1] - Phi[j] - Rho[j+1]*dx**2/Sample.Dieelek[j+1]	
		# Second derivative is calculated by now
		
		# Calculation of the total resistance of the Sample
		N = TotalResistance(Sample)
		Ntotal = sum(N)
		
		Phi2 = zeros(nx)
		Phi2[0] = Sample.Bias + u.BarrierHeight - Phi[0]
		for i in range(len(Phi)-1):
			Phi2[i+1] = Phi2[i] - (N[i]/Ntotal)*Phi2[0]
		
		Sample.Phi = Phi2 + Phi
		
		#Sample.Phi = Sample.Phi + (Phi_ende-Sample.Phi[-1])*(Sample.X-Sample.X[0])/dx
	else:
		if(Sample.W >= u.DeviceLength):
			N = TotalResistance(Sample)
			Ntotal = sum(N)
			
			Phi2 = zeros(nx)
			Phi2[0] = Sample.Bias + u.BarrierHeight - Phi[0]
			for i in range(len(Phi)-1):
				Phi2[i+1] = Phi2[i] - (N[i]/Ntotal)*Phi2[0]
		
			Sample.Phi = Phi2 + Phi
			
			#Phi_start = deepcopy(Sample.Phi[0])
			#Sample.Phi = Sample.Phi + (Sample.Bias - u.BarrierHeight - Phi_start)*(Sample.X[-1]-Sample.X)/Sample.X[-1]
	return
	
def TotalResistance(Sample):
	dx = Sample.dx
	N = zeros(Sample.nx)
	for j in range(Sample.nx) :
		N[j] = u.R0/(exp((Sample.DopingDensity[j]-u.Titanium_DopingDensity)*1.8)+1)
		#if(  Sample.DopingDensity[j]  > 1.1*u.Titanium_DopingDensity) :  N[j] = u.R2*dx
		#elif(Sample.DopingDensity[j] <  0.9*u.Titanium_DopingDensity) :  N[j] = u.R3*dx
		#else: 															 N[j] = u.R1*dx
	return N
	
#-----------------------------------------------------------
