from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Diffusion
#import T2 as T2
import Schottky as S

Sample = S.Schottky(0, 100, 50000)
#Sample.DopingDrift(144)
#Sample.PlotResults()

D =  loadtxt("Solution")
#Sample.Phi				= D[0]
Sample.DopingDensity	= D[1]
#Sample.ChargeDensity	= D[2]

N0 = deepcopy(Sample.DopingDensity)

n = 100000
#Sample.Diffusion(3, 1)
for i in range(n):
	Sample.DopingDensity[:5000] = Diffusion(Sample.X[:5000], 10, Sample.DopingDensity[:5000], 1)
#def Diffusion(X, dt, DopingDensity, DiffusionCoeffcient):

print(amax(Sample.DopingDensity), mean(Sample.DopingDensity))
print(n)
print(sum(N0),sum(Sample.DopingDensity),sum(N0)-sum(Sample.DopingDensity))
print(sum(N0[:5]), sum(Sample.DopingDensity[:5]))
print(sum(N0[5:300]), sum(Sample.DopingDensity[5:300]))

plt.figure("Diffusion")
plt.plot(Sample.X, N0, color = "Black")
plt.plot(Sample.X, Sample.DopingDensity, color = "Red")

#savetxt("Solution", (Sample.Phi, Sample.DopingDensity, Sample.ChargeDensity))

plt.show()
