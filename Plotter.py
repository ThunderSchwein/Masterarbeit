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
#Import(Sample, "V=10-10(1)T=2000s")
Import(Sample, "V=0T=200min")
Sample.PlotResults()
"""
D =  loadtxt("Resistance0")
#[(T, R0, R_1, R1, C, Sample.W)]
#print(D)


T	= D[:,0]
R0	= D[:,1]
R_1	= D[:,2]
R1  = D[:,3]
C   = D[:,4]
W   = D[:,5]
V   = D[:,6]


plt.figure("Resistance")

R2, = plt.plot(T, R0/u.Sigma1, label = "Original R")
R3, = plt.plot(T, R_1/u.Sigma3, label = "High R")
R4, = plt.plot(T, R1/u.Sigma2, label = "Low R")

Rtotal = R1/u.Sigma2 + R_1/u.Sigma3 + R0/u.Sigma1
R5, = plt.plot(T, Rtotal, label = "Total R")
R6, = plt.plot(T, V/Rtotal*1000, label = "Curent/mA")
plt.legend()

plt.figure("Capacity")
plt.plot(T, C)
plt.plot(T, W)
plt.figure("Mott-Schottky")
plt.plot(V, 1/(C*C))
"""
plt.show()
