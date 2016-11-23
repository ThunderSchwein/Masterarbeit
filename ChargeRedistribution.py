from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *

def Reload(X, ChargeDensity, DopingDensity, FermiLevel = 0):
	S = -sum(ChargeDensity)
	#print(S)
	i = 0
	D = DopingDensity
	while(S > D[i]):
		ChargeDensity[i] =  -2*D[i]
		S = S-2*D[i]
		i = i+1
	ChargeDensity[i] = 2*D[i] - S
	ChargeDensity[i+1 :] = 0
	
	return ChargeDensity
	
def ReloadFermi(X, DopingDensity, Potential, FermiLevel = 0):
	ChargeDensity = zeros(len(X))
	ChargeDensity[1:] = (exp(-Potential[:-1]/(u.Kb*1e-10)) - 1 )*DopingDensity[1:]*u.Titanium_DopingCharge
	return ChargeDensity
	
"""
#-----------------------------------
# Depletion Region Width Correction

plt.figure()

plt.plot(X, Phi, color = "red")
DP = ones(nx)*u.Titanium_DopingDensity
#DP[:int(0.2*nx)] = 1
plt.plot(X, Phi, color= "Black")
Rho2 = zeros(nx)
W2 = W/0.962
n3 = zeros(nx)
#Rho2[:int(W2/u.DeviceLength*nx)] = -u.Titanium_DopingCharge*u.Titanium_DopingDensity
#Phi[0] = u.BarrierHeight; Phi[-1] = 0.0
n0 = (1-exp(-Phi/(u.Kb*300)))*u.Titanium_DopingDensity*(-u.Titanium_DopingCharge)
#n3[1:] = (1-exp(-Phi[:-1]/(u.Kb*1e-10)))*u.Titanium_DopingDensity*(-u.Titanium_DopingCharge)
n3 = ReloadFermi(X, DP, Phi)
Phi2 = solverLU_backward(X, Phi, n0, Dieelek)
#n0[DP == 0] = 0
plt.plot(X, n0, color = "Purple")
plt.plot(X, n3, color = "blue" )
plt.plot(X, Rho, color = "green")

#Phi2 = solverLU_backward(X, Phi, n0, Dieelek)
#plt.plot(X, Phi2)
"""
