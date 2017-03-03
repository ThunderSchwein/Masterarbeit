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
	Sample.E = E_Field(Sample.X, Sample.Phi, Sample.E)
#--------------------------------
def Run(Bias, dt, nx):
	Sample = Init_Sample(Bias,dt,nx)
	
	"""
	#Change the initial Charge Distribution
	Sample.W = 550
	Sample.DopingDensity[:int((Sample.W-67)/Sample.dx)] = 0.0
	Sample.DopingDensity[1] = 50.0 - sum(Sample.DopingDensity)
	#print(sum(Sample.DopingDensity))
	"""
	#ReadFile
	#Import(Sample, "InitialConditions")
	#Import(Sample, "V=0T=200min")

	#OxygenDrift(Sample, 500)
	Sample.PlotResults()

	"""
	# Voltage fct
	V = []
	for i in range(100):
		V = V+[i*0.05]
	for i in range(200):
		V = V+[5-i*0.05]
	for i in range(100):
		V = V + [-5 + i*0.05]
	V = V + [0]
	"""
	"""
	# Voltage fct
	V = []
	for i in range(100):
		V = V + [i*0.1]
	V = V + [0]
	#for i in range(100):
		#V = V + [5-i*0.1]
	#V = V + [0]
	"""
	V = ones(20)*0
	
	Sample.Bias = V[0]
	#print(V)
	
#-------------------------------
	T = 0
	R = []

	# Drift paired with Diffusion
	DiffusionSteps = 1
	DriftIterations = int(1/Sample.dt)*10

	print("Start")
	for i in range(len(V)) :#fix20 = 20mal Bilder Anzeigen
		for j in range(DriftIterations) :
			# Time actualization
			T = T+dt
			
			print("Drift")
			Sample.Phi[0] = Sample.Bias - u.BarrierHeight
			solverLU_backward(Sample)
			ReloadDeltaPhi(Sample)
			solverLU_backward(Sample)
			OxygenDrift(Sample, 1)
			
			#print("Diffusing")
			#OxygenDiff(Sample, DiffusionSteps, 1)
			#for h in range(DiffusionSteps) :
				#Sample.DopingDensity[:20000] = Diffusion(Sample.X[:20000], Sample.dt , Sample.DopingDensity[:20000], 1)
			
			#[Sample.ChargeDensity, Sample.W]= ReloadDeltaPhi(Sample.X, Sample.DopingDensity, Sample.ChargeDensity, Sample.Phi, Sample.Dieelek, Sample.W, Sample.Bias)
			
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
			del R0, R_1, R1, R3, C
		#print(Sample.Bias)
		print(sum(Sample.DopingDensity), T)
		print(i, Sample.W)
		Sample.PlotResults()
		Sample.Bias = V[i]
		print("New Bias =", Sample.Bias)
		ReloadDeltaPhi(Sample)
	Sample.PlotResults()

	savetxt("Resistance5-5V(0)", R , newline="\n")
	#print(R)

	WriteToFile(Sample, "V=5-5(0)")

	plt.show()
	return 0
	
if __name__ == '__main__':
	Run(0, 1, 50001)
