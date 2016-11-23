from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy

import Units as u
from SolverLU import *
from DopingDrift import *
from ChargeRedistribution import *

# 0. Define Grid
dt = 1
nx = 1000000
X = linspace(0, u.DeviceLength, nx)
dx = X[1] - X[0]

#------------------------------------
# 1. Calculate Initial Conditions
plt.figure("Initial Conditions")

Dieelek = ones(nx)*u.e0*u.Titanium_DielectricityFactor
# Depletion Region Width
W = (2*(u.Titanium_DielectricityFactor*u.e0)*u.BarrierHeight/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5

#Doping density
ChargeDensity = zeros(nx)
ChargeDensity[:int(W/u.DeviceLength*nx)] = -u.Titanium_DopingCharge*u.Titanium_DopingDensity
plt.plot(X,ChargeDensity)

Phi = zeros(nx)
#Phi[0] = u.barrier; Phi[-1] = 0.0
Phi = solverLU_backward(X, Phi, ChargeDensity, Dieelek)
E = E_Field(X, Phi)

plt.plot(X, Phi,color = "G")
plt.plot(X,E, color = "R")

#------------------------------------
# 2. Calculate Doping Drift
DP = zeros(nx)
DP[:] = u.Titanium_DopingDensity
DP1 = Drift( X , Phi, DP, dt)

for i in range(10) :
	#n1 = Reload(X, ChargeDensity, DP1)
	n1 = ReloadFermi(X, DP1, Phi, FermiLevel = 0)
	Phi3 = solverLU_backward(X, Phi, n1, Dieelek)
	DP1 = Drift( X , Phi, DP1, dt)
	
E2 = E_Field(X, Phi3)

plt.figure("Oxygen Vacancy Movement")
plt.plot(X, DP, color = "Grey")
plt.plot(X, n1, color = "Red")
plt.plot(X, DP1, color = "Blue")

plt.figure("Potentials")
plt.plot(X, E2, color = "Blue")
plt.plot(X, Phi3, color = "Green")
plt.plot(X, DP1, color = "Grey")

#solnge durchführen, bis es zum alten Wert konvergiert !!

#savetxt("Solution", (Phi3, DP1, n1))
plt.show()
