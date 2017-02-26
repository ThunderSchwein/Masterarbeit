from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S

def Init_Sample(Bias, dt, nx):
	return S.Schottky(Bias, dt, nx)

def OxygenDrift(Sample, Iterations):
	for i in range(Iterations):
		Sample.DopingDrift()
	
def OxygenDiff(Sample, Iterations, DiffCoefficient):
	for i in range(Iterations):
		Sample.Diffusion(DiffCoefficient)

def WriteToFile(Sample, File):
	Sample.WriteToFile(File)

def Sigma(Rho):
	if (Rho==u.Titanium_DopingDensity) : return u.Sigma1
	if (Rho <u.Titanium_DopingDensity) : return u.Sigma2
	if (Rho >u.Titanium_DopingDensity) : return u.Sigma3
	return 0

def Import(Sample, File):	
	D =  loadtxt(File)
	Sample.Phi				= D[0]
	Sample.DopingDensity	= D[1]
	Sample.ChargeDensity	= D[2]
	Sample.W = len(D[2][D[2] != 0])*Sample.dx
#--------------------------------
#Sample = S.Schottky(0, 1, 50001)
Bias = 0 # volt
dt = 0.1	# Seconds
nx = 50001	# 

Sample = Init_Sample(Bias,dt,nx)
"""
#Change the initial Charge Distribution
Sample.W = 550
Sample.DopingDensity[:int((Sample.W-67)/Sample.dx)] = 0.0
Sample.DopingDensity[1] = 50.0 - sum(Sample.DopingDensity)
#print(sum(Sample.DopingDensity))
"""
#ReadFile
Import(Sample, "InitialConditions")

#OxygenDrift(Sample, 500)
Sample.PlotResults()

# Voltage fct
V = []
for i in range(100):
	V = V+[i*0.05]
for i in range(200):
	V = V+[5-i*0.05]
for i in range(100):
	V = V + [-5 + i*0.05]
V = V + [0]

print(V)
#-------------------------------
T = 0
R = []

# Drift paired with Diffusion
DiffusionSteps = 1
DriftIterations = int(1/Sample.dt)

print("Start")
for i in range(len(V)) :#fix20 = 20mal Bilder Anzeigen
	for j in range(DriftIterations) :
		# Time actualization
		T = T+dt
		
		print("Drift")
		OxygenDrift(Sample, 1)
		
		print("Diffusing")
		OxygenDiff(Sample, DiffusionSteps, 1)
		#for h in range(DiffusionSteps) :
			#Sample.DopingDensity[:20000] = Diffusion(Sample.X[:20000], Sample.dt , Sample.DopingDensity[:20000], 1)
		
		[Sample.ChargeDensity, Sample.W]= ReloadDeltaPhi(Sample.X, Sample.DopingDensity, Sample.ChargeDensity, Sample.Phi, Sample.Dieelek, Sample.W, Sample.Bias)
		
		#Resistance Calculation
		R0 = R_1 = R1 = R3 = C = 0
		# R0 = original density, R_1 = lower, R1 = higher
		for g in range(nx) :
			Rho = Sample.DopingDensity[g]
			if (Rho==u.Titanium_DopingDensity) : R0  = R0  + Sample.dx
			if (Rho <u.Titanium_DopingDensity) : R_1 = R_1 + Sample.dx
			if (Rho >u.Titanium_DopingDensity) : R1  = R1  + Sample.dx
			#R0 = R0 + Sample.dx * Sigma(Sample.DopingDensity[g])
		C = sum(Sample.ChargeDensity)*Sample.dx
		R = R +[(T, R0, R_1, R1, C, Sample.W, Sample.Bias)]
		
	#print(Sample.Bias)
	print(i, Sample.W)
	#Sample.PlotResults()
	Sample.Bias = V[i]
	print("New Bias =", Sample.Bias)
	[Sample.ChargeDensity, Sample.W]= ReloadDeltaPhi(Sample.X, Sample.DopingDensity, Sample.ChargeDensity, Sample.Phi, Sample.Dieelek, Sample.W, Sample.Bias)
Sample.PlotResults()

savetxt("Resistance-1", R , newline="\n")
print(R)

WriteToFile(Sample, "RisingV0--5")

plt.show()

"""
N0 = deepcopy(Sample.DopingDensity)

n = 100000
#Sample.Diffusion(3, 1)
for i in range(n):
	Sample.DopingDensity[:5000] = Diffusion(Sample.X[:5000],Sample.dt, Sample.DopingDensity[:5000], 1)
#def Diffusion(X, dt, DopingDensity, DiffusionCoeffcient):

# Diffusion Info
print(amax(Sample.DopingDensity), mean(Sample.DopingDensity))
print(n)
print(sum(N0),sum(Sample.DopingDensity),sum(N0)-sum(Sample.DopingDensity))
print(sum(N0[:5]), sum(Sample.DopingDensity[:5]))
print(sum(N0[5:300]), sum(Sample.DopingDensity[5:300]))

plt.figure("Diffusion")
plt.plot(Sample.X, N0, color = "Black")
plt.plot(Sample.X, Sample.DopingDensity, color = "Red")
"""
#savetxt("Solution", (Sample.Phi, Sample.DopingDensity, Sample.ChargeDensity))
