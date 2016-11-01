from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy

nx = 120000
X = linspace(-1,1, nx)


"""
def solver(X, Phi, Rho, Dieelek):
	Result = deepcopy(Phi)
	Phi2 = Phi[2:] + Phi[:-2]
	for i in range(len(Phi)-2) :
		Result[i+1] = ( Phi2[i] - Rho[i+1]/Dieelek[i+1]*dx**2)/2
	return Result

def solver(X, Phi, Rho, Dieelek):
	dx = X[1]-X[0]
	Phi2 = Phi[2:] + Phi[:-2]
	for i in range(len(Phi)-2) :
		Phi[i+1] = ( Phi2[i] - Rho[i+1]/Dieelek[i+1]*dx**2)/2
	return Phi


def solver(X, Phi, Rho, Dieelek):
	Phi2 = Phi[2:] + Phi[:-2]
	for i in range((len(Phi))//2) :
		j = 2*i+1
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
	Phi2 = Phi[2:] + Phi[:-2]
	for i in range((len(Phi))//2-1) :
		j = 2*(i+1)
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
		#print(i, 2*(i+1), 2*i+1) i, Gerade außer 0, ungerade außer letzten Eintrag
		#Phi[2*i+1] = ( Phi2[2*i] - Rho[2*i+1]/Dieelek[2*i+1]*dx**2)/2
	return Phi

def solver(X, Phi, Rho, Dieelek):
	dx = X[1]-X[0]
	#1
	Phi2 = Phi[2:] + Phi[:-2]
	for j in range(1,len(Phi)-1,3) :
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
	#2
	Phi2 = Phi[2:] + Phi[:-2]
	for j in range(2,len(Phi)-1,3) :
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
	
	#1
	Phi2 = Phi[2:] + Phi[:-2]
	for j in range(1,len(Phi)-1,3) :
		#j = 3*i+1
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
	#2
	Phi2 = Phi[2:] + Phi[:-2]
	for j in range(2,len(Phi)-1,3) :
		#j= 3*(i+1)-1
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
	#3
	for j in range(3,len(Phi)-1,3) :
		#j= 3*(i+1)
		Phi[j] = ( Phi2[j-1] - Rho[j]/Dieelek[j]*dx**2)/2
		
	return Phi
"""

def solverLU(X, Phi, Rho, Dieelek):
	dx = X[1]-X[2]
	
	Phi[1] = Phi[0] - Rho[-1]*dx
	#Phi[-2] = 2*( Phi[-1] - Rho[-1]*dx*dx )
	
	Phi_ende = Phi[-1]
	for j in range(len(Phi)-2):
		Phi[j+2] = 2*Phi[j+1] - Phi[j] + 2*Rho[j+1]*dx**2
	
	Phi = (Phi + Phi[0])/2
	Phi = Phi + (Phi_ende-Phi[-1])*(X-X[0])/2
	
	return Phi

#-----------------------------------------------------------

#Test 1
plt.figure()
A = loadtxt("Optimale Lösung")
#plt.plot(X,A, color = "r")

Rho = zeros(nx)
Dieelek = Rho + 1
Phi = zeros(nx) 
#for i in range(len(Phi)-1) :
	#Phi[i+1] = Phi[i] + 1

Rho = 12*X**2
#Rho  = cos(2*pi*X) 

plt.plot(X,Rho)

Dieelek = ones(nx)
Phi = zeros(nx)
Phi[0] = 1; Phi[-1] = 1


Phi = solverLU(X, Phi, Rho, Dieelek)
#for i in range(int(nx**1.5)) :
	#Phi = solver(X, Phi, Rho, Dieelek)

plt.plot(X, X**4,color = "R")
#plt.plot(X,-cos(2*pi*X), color = "R")
plt.plot(X, Phi, color = "G")

savetxt("Optimale Lösung", Phi)
plt.show()
