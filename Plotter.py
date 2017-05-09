from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
#from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S
from Run import Init_Sample

File = "V=10(05)T=201s"
D =  loadtxt(File)
nx = len(D[0,:])
Sample = Init_Sample(0, 1, nx)
Sample.Import(File)

#Sample.Bias = -30
#Sample.W = 100000000000
#solverLU(Sample)
#print(Sample.Phi[0])

Sample.PlotResults()

D =  loadtxt("Resistance10V(05)")

T	= D[0,:]
V   = D[1,:]
R   = D[2,:]
I   = D[3,:]
W 	= D[4,:]

plt.figure("Resistance")
R5, = plt.plot(V, R*1e-6, label = "Total R [MOhm]")
plt.xlabel("t[min]")
plt.legend()

plt.figure("Current")
R9, = plt.plot(V, I, label = "Current")
#R7, = plt.plot(V[:int(len(V)/2)], I[:int(len(V)/2)], label = "Voltage_rising")
#R6, = plt.plot(V[int(len(V)/2):], I[int(len(V)/2):], label = "Current_decreasing")
plt.xlabel("U[V]")
plt.legend(loc = 2)

plt.figure("Current2")
R8, = plt.plot(I*1e6, label = "Current")
plt.xlabel("U [V]")
plt.ylabel("I [µA]")
plt.legend(loc = 2)

plt.figure("Depletion")
R9, = plt.plot(T*5, W, label = "Breite der Verarmungsregion")
plt.legend(loc = 2)

plt.show()
