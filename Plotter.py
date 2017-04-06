from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S
from Run import Init_Sample, Import

Sample = Init_Sample(0, 1, 50001)
Import(Sample, "V=10(5)T=66000s")
#Import(Sample, "V=0T=200min")
Sample.PlotResults()

D =  loadtxt("Resistance10V(0)")

T	= D[0,:]
R   = D[1,:]
C   = D[2,:]

plt.figure("Resistance")
R5, = plt.plot(T/10, R*1e-6, label = "Total R [MOhm]")
R6, = plt.plot(T/10, C, label = "Capacity")
plt.xlabel("t[min]")
plt.legend()

#plt.figure("Capacity")
#plt.plot(T, C)
#plt.figure("Mott-Schottky")
#plt.plot(V, 1/(C*C))

plt.show()
