from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S

VoltageOnMetal = 0 # V
TimeStep = 1 # second
GridPointsNumber = 50001 # per u.DeviceLength = 5µm => 10 grid points per nm

#-------------------------------

def Init_Sample(Bias, dt, nx):
	print("A Sample has been created.")
	return S.Schottky(Bias, dt, nx)

def OxygenDrift(Sample, Iterations):
	for i in range(Iterations):
		Sample.DopingDrift()
	return
	
def OxygenDiff(Sample, Iterations, DiffCoefficient):
	for i in range(Iterations):
		Sample.Diffusion(DiffCoefficient)
	return
		
def WriteToFile(Sample, File):
	Sample.WriteToFile(File)
	print("Sucessfully saved the data in '", File, "' File.")
	return
	
def Sigma(Rho):
	if (Rho==u.Titanium_DopingDensity) : return u.Sigma1
	if (Rho <u.Titanium_DopingDensity) : return u.Sigma2
	if (Rho >u.Titanium_DopingDensity) : return u.Sigma3
	return 0

def Import(Sample, File):
	# Import Data of a Sample from a File.
	# The length of Data arrays must match.
	D =  loadtxt(File)
	Sample.Phi				= D[0]
	Sample.DopingDensity	= D[1]
	Sample.ChargeDensity	= D[2]
	Sample.Bias = Sample.Phi[0] - u.BarrierHeight 
	Sample.W = len(D[2][D[2] != 0])*Sample.dx
	Sample.E = E_Field(Sample.X, Sample.Phi, Sample.E)
	print("Import of Sample Data  from File '", File, "' has been successful.")
	return
	
def SetDepletionWidth(Sample, Width):
	#Changes the initial Charge Distribution
	Sample.W = Width
	Sample.ChargeDensity[:] = 0.0
	Sample.ChargeDensity[:int(Sample.W/Sample.dx)] = Sample.DopingDensity[:int(Sample.W/Sample.dx)]
	solverLU(Sample)
	print("Change of depletion Region length successful.")
	return
	
def AddVoltageTrace(Vstart, Vend, Steps = 0, V = []):
	for i in range(Steps):
		V = V+[Vstart + (Vend - Vstart)*i/Steps]
	print("Trace", Vstart, " ", Vend, "has been added to the Voltage trace.")
	return V

#--------------------------------
def Run(Bias, dt, nx):
	# Sample Initialization
	Sample = Init_Sample(Bias,dt,nx)
	
	# Read Sample Data from a File
	#Import(Sample, "V=0T=200min")
	Import(Sample, "V=10-10(2)T=3000s")
	
	# Plot Initial Conditions
	Sample.PlotResults()

	V = AddVoltageTrace(-10, 0, 100)
	#print(V)
	
#-------------------------------
	T = 0
	R = []

	# Drift paired with Diffusion
	DiffusionSteps = 1
	DriftIterations = int(1/Sample.dt)*10

	print("Start of Simulation Process")
	for i in range(len(V)) :#fix20 = 20mal Bilder Anzeigen
		Sample.Bias = V[i]
		print("New Bias =", Sample.Bias, "V")
	
		for j in range(DriftIterations) :
			# Time actualization
			T = T+dt
			
			print("Drift")
			Sample.Phi[0] = Sample.Bias + u.BarrierHeight
			ReloadDeltaPhi(Sample)
			OxygenDrift(Sample, 1)
			
			#print("Diffusing")
			#OxygenDiff(Sample, DiffusionSteps, 1)
			#for h in range(DiffusionSteps) :
				#Sample.DopingDensity[:20000] = Diffusion(Sample.X[:20000], Sample.dt , Sample.DopingDensity[:20000], 1)
			#ReloadDeltaPhi(Sample)
			
			#Resistance Calculation
			"""
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
			"""
		print(T, sum(Sample.DopingDensity))
		print("Time = ",T, "s, 'Drift' Method executed", (i+1)*DriftIterations, " times")
		print("Depletion region width = ", Sample.W, " nm")
		#Sample.PlotResults()
		ReloadDeltaPhi(Sample)
	print("End of Simulation Process")
	Sample.PlotResults()

	# Save Results
	#savetxt("Resistance5-5V(0)", R , newline="\n")
	WriteToFile(Sample, "V=10-10(3)T=4000s")

	plt.show()
	return 0
	
if __name__ == '__main__':
	Run(VoltageOnMetal, TimeStep, GridPointsNumber)
