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
	dx = X[1]-X[2]
	
	Phi[1] = Phi[0] #+ Rho[1]*dx
	#Phi[-2] = 2*( Phi[-1] - Rho[-1]*dx*dx )
	
	Phi_ende = Phi[-1]
	for j in range(len(Phi)-2):
		Phi[j+2] = 2*Phi[j+1] - Phi[j] - Rho[j+1]*dx**2/Dieelek[j+1]
	#Phi = (Phi + Phi[0])/2
	Phi = Phi + (Phi_ende-Phi[-1])*(X-X[0])/(X[-1]-X[0])
	return Phi

def E_Field(X, Phi):
	dx = X[1]-X[2]
	E = deepcopy(X)
	for i in range(len(Phi)-1) :
		E[i] = (Phi[i+1] - Phi[i])/dx
	E[-1] = 0
	return -E
	
def solverLU_backward(X, Phi, Rho, Dieelek):
	dx = X[1]-X[2]
	n = len(Phi)
	
	Phi[-2] = Phi[-1] #+ Rho[1]*dx
	#Phi[-2] = 2*( Phi[-1] - Rho[-1]*dx*dx )
	
	Phi_ende = Phi[0]
	for j in range(len(Phi)-2):
		Phi[n-3-j] = 2*Phi[n-j-2] - Phi[n-j-1] - Rho[n-j-2]*dx**2/Dieelek[n-j-2]
	#Phi = (Phi + Phi[0])/2
	#Phi = Phi + (Phi_ende-Phi[-1])*(X-X[0])/(X[-1]-X[0])
	return Phi
	
#-----------------------------------------------------------
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
