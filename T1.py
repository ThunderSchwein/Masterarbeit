from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK
#import T2 as T2
import Schottky as S

Sample = S.Schottky(0, 10, 50000)
Sample.DopingDrift(1)
Sample.PlotResults()

#N = deepcopy(Sample.DopingDensity)
#Sample.Diffusion(50, 0.1)
#plt.figure("Diffusion")
#print (sum(N)-sum(Sample.DopingDensity))
#plt.plot(Sample.X, N, color = "Black")
#plt.plot(Sample.X, Sample.DopingDensity, color = "Red")

savetxt("Solution", (Sample.Phi, Sample.DopingDensity, Sample.ChargeDensity))

plt.show()
