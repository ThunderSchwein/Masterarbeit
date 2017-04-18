from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
#from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S
from Run import Init_Sample, Import

Sample = Init_Sample(0, 1, 5001)
Import(Sample, "V=10(01)T=100s")
#Import(Sample, "V=0T=200min")
Sample.PlotResults()
plt.show()
D =  loadtxt("Resistance10V(01)")

T	= D[0,:]
V   = D[1,:]
R   = D[2,:]

plt.figure("Resistance")
R5, = plt.plot(V, R*1e-6, label = "Total R [MOhm]")
plt.xlabel("t[min]")
plt.legend()

plt.figure("Current")
R7, = plt.plot(V[int(len(V)/2):], V[int(len(V)/2):]/(R[int(len(V)/2):]*1e-6), label = "Current_rising")
R6, = plt.plot(V[:int(len(V)/2)], V[:int(len(V)/2)]/(R[:int(len(V)/2)]*1e-6), label = "Current_decreasing")
plt.xlabel("U[V]")
plt.legend()

#savetxt("Resistance10V(001)", (T, V, R), newline="\n")

#plt.figure("Capacity")
#plt.plot(T, C)
#plt.figure("Mott-Schottky")
#plt.plot(V, 1/(C*C))

plt.show()
