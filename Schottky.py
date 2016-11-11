from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *

nx = 10000
X = linspace(0, u.Length, nx)
dx = X[1] - X[0]

#Test 1
plt.figure()

Dieelek = ones(nx)*u.e0*u.Ti_es

Phi = zeros(nx)
Rho = zeros(nx)

#Doping density
DP = zeros(nx) 

# Depletion Region Width
W = (2*(u.Ti_es*u.e0)*u.barrier/(2*u.Ti_nd))**.5

# Rho neu definieren
Rho[:int(W/u.Length*nx)] = -2*u.Ti_nd
Phi[0] = u.barrier; Phi[-1] = 0.0

plt.plot(X,Rho)

Phi = solverLU_backward(X, Phi, Rho, Dieelek)
E = E_Field(X, Phi)
plt.plot(X, Phi,color = "G")
#plt.plot(X,E, color = "R")

#savetxt("Optimale Lösung", Phi)
plt.show()
