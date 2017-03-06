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
	
	if(Phi[0] > abs(u.BarrierHeight)):
		Sample.Phi[1] = Sample.Phi[0] = Sample.Bias + u.BarrierHeight
		#print(Sample.Phi[1])
		Phi_ende = deepcopy(Sample.Phi[-1])
		for j in range(len(Phi)-2):
			Sample.Phi[j+2] = 2*Sample.Phi[j+1] - Sample.Phi[j] - Rho[j+1]*dx**2/Sample.Dieelek[j+1]

		Sample.Phi = Sample.Phi + (Phi_ende-Sample.Phi[-1])*(Sample.X-Sample.X[0])/(Sample.X[-1]-Sample.X[0])
		#print(Sample.Phi[1])
		del Phi_ende
	else:
		Phi[-1] = Phi[-2] = 0	
		for j in range(len(Phi)-2):
			Phi[nx-3-j] = 2*Phi[nx-j-2] - Phi[nx-j-1] - Rho[nx-j-2]*dx**2/Sample.Dieelek[nx-j-2]

	return
	
#-----------------------------------------------------------
