from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy

import Units as u
from SolverLU import *
from DopingDrift import *
from ChargeRedistribution import *

# 0. Define Grid
dt = 100
nx = 500000
X = linspace(0, u.DeviceLength, nx)
dx = X[1] - X[0]

#------------------------------------
# 1. Calculate Initial Conditions
plt.figure("Initial Conditions")

Dieelek = ones(nx)*u.e0*u.Titanium_DielectricityFactor
# Depletion Region Width
W = (2*(u.Titanium_DielectricityFactor*u.e0)*abs(u.BarrierHeight)/(u.Titanium_DopingCharge*u.Titanium_DopingDensity))**.5

#Doping density
ChargeDensity = zeros(nx)
ChargeDensity[:int(W/u.DeviceLength*nx)] = u.Titanium_DopingCharge*u.Titanium_DopingDensity
plt.plot(X,ChargeDensity)

Phi = zeros(nx)
Phi[1] = Phi[0] = u.BarrierHeight; Phi[-1] = 0.0
Phi = solverLU(X, Phi, ChargeDensity, Dieelek)
# Mirror Charge Correction
#Phi = Phi - 1/((X+1e-1)*16*pi*u.e0*u.Titanium_DielectricityFactor)

E = E_Field(X, Phi)

Phi5 = solverLU_backward(X, Phi, ChargeDensity, Dieelek)
plt.plot(X, Phi5, color = "Red")

plt.plot(X, Phi,color = "G")
plt.plot(X,E, color = "R")
#------------------------------------
# 2. Calculate Doping Drift
DP0 = zeros(nx)
DP0[:] = u.Titanium_DopingDensity
DP = DriftRK( X , Phi, DP0, dt)

for i in range(10) :
	#ChargeDensity = Reload(X, ChargeDensity, DP)
	#ChargeDensity = ReloadFermi(X, DP, Phi, E, FermiLevel = 0)
	[ChargeDensity,W]  = ReloadDeltaPhi(X, DP, ChargeDensity, Phi, Dieelek, W, FermiLevel = 0)
	
	Phi = solverLU(X, Phi, ChargeDensity, Dieelek)
	DP = DriftRK( X , Phi, DP, dt)

#ChargeDensity = ReloadFermi(X, DP, Phi, E, FermiLevel = 0)
[ChargeDensity,W]  = ReloadDeltaPhi(X, DP, ChargeDensity, Phi, Dieelek, W, FermiLevel = 0)

Phi = solverLU(X, Phi, ChargeDensity, Dieelek)
E2 = E_Field(X, Phi)

plt.figure("Oxygen Vacancy Movement")
plt.plot(X, DP0, color = "Grey")
plt.plot(X, ChargeDensity, color = "Red")
plt.plot(X, DP, color = "Blue")

#plt.plot(X, E, color = "Purple")
#plt.plot(X, DP-ChargeDensity)

plt.figure("Potentials")
plt.plot(X, E2, color = "Blue")
#plt.plot(X, Phi, color = "Green")
plt.plot(X, DP, color = "Grey")

Phi5 = solverLU_backward(X, Phi, ChargeDensity, Dieelek)
plt.plot(X, Phi5, color = "Red")
Phi[1] = Phi[0] = u.BarrierHeight;
Phi = solverLU(X, Phi, ChargeDensity, Dieelek)
plt.plot(X, Phi, color = "Green")

#plt.plot(X, -u.e*u.e/((X+1e-1)*16*pi*u.e0*u.Titanium_DielectricityFactor*u.ds)+Phi3, color = "Green")

#solnge durchführen, bis es zum alten Wert konvergiert !!

savetxt("Solution", (Phi, DP, ChargeDensity))
plt.show()
