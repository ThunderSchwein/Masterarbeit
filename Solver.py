from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u

"""
def solver(X, Phi, Rho, Dieelek):
	dx = X[1]-X[0]
	Phi2 = Phi[2:] + Phi[:-2]
	for i in range(len(Phi)-2) :
		Phi[i+1] = ( Phi2[i] - Rho[i+1]/Dieelek[i+1]*dx**2)/2
	return Phi
"""

def solverLU(X, Phi, Rho, Dieelek):
	dx = X[1]-X[0]
	Phi[1] = Phi[0]# = u.BarrierHeight
	
	Phi_ende = Phi[-1]
	
	Rho2 = zeros(len(X)-1)
	for i in range(len(X)-1) :
		Rho2[i] = (Rho[i+1] + Rho[i])/2

	for j in range(len(Phi)-2):
		Phi[j+2] = 2*Phi[j+1] - Phi[j] - Rho2[j+1]*dx**2/Dieelek[j+1]
	
	Phi = Phi + (Phi_ende-Phi[-1])*(X-X[0])/(X[-1]-X[0])
	return Phi

def E_Field(X, Phi, E):
	dx = X[1]-X[0]
	E = deepcopy(X)
	for i in range(len(X)-1) :
		E[i+1] = (Phi[i] - Phi[i+1])/dx
	E[0] = E[1]
	
	return E
	
def solverLU_backward(Sample):
	dx = Sample.dx
	nx = Sample.nx
	Phi = Sample.Phi
	Rho = Sample.ChargeDensity
	
	if(Phi[0] > abs(u.BarrierHeight)):
		for i in range(nx):
			Phi[i] = Phi[0] - (Phi[0]-Phi[-1]) * i/(nx-1)
	else :
		Rho2 = zeros(nx-1)
		Phi[-1] = Phi[-2] = 0
		for i in range(nx-1) :
			Rho2[i] = (Rho[i+1] + Rho[i])/2
			
		for j in range(len(Phi)-2):
			Phi[nx-3-j] = 2*Phi[nx-j-2] - Phi[nx-j-1] - Rho2[nx-j-2]*dx**2/Sample.Dieelek[nx-j-2]
			
		del Rho2
	return
	
#-----------------------------------------------------------


"""
#Test 2
n = 1e6
X = linspace(0,1000, n)
Rho = zeros(n)
Rho[:200000] = 2.5e-3
Rho[200000:500000] = 1e-3
Dieelek = ones(n)*u.Titanium_DielectricityFactor*u.e0
Phi = zeros(n)

Phi[0] = -40
Phi = solverLU(X, Phi, Rho, Dieelek)

#Phi = solverLU_backward(X, Phi, Rho, Dieelek)
E = E_Field(X,Phi)
#plt.plot(X, E)
plt.plot(Rho+.01)
#plt.plot(X, Phi)
plt.show()
"""
"""
nx = 10000
X = linspace(0, 2, nx)

#Test 1
plt.figure()

#A = loadtxt("Optimale Lösung")
#plt.plot(X,A, color = "r")
ES = 10
Dieelek = ones(nx)
Phi = zeros(nx)
Rho = zeros(nx)

#Doping density
DP = zeros(nx) 

# Rho neu definieren
#Rho[:int(.01*nx)] = 0.0009

# test 0
Phi[0] = 0; Phi[-1] = -8
Rho = 6*X
plt.plot(X, -X**3,color = "R")

# test 1
#Phi[0] = 5; Phi[-1] = -11
#Rho = 12*X**2
#plt.plot(X, -X**4 + 5,color = "R")

# test 2
#Phi[0] = 1/(4*pi**2); Phi[-1] = 1/(4*pi**2)
#Rho  = cos(2*pi*X) 
#plt.plot(X, cos(2*pi*X)/(4*pi**2), color = "R")

plt.plot(X,Rho)

Phi = solverLU_backward(X, Phi, Rho, Dieelek)
#E = E_Field(X, Phi)

plt.plot(X, Phi, color = "G")
#plt.plot(X,E)

#savetxt("Optimale Lösung", Phi)

plt.show()
"""
