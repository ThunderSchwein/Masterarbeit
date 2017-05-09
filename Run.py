from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Smoothing
from ChargeRedistribution import *
import Schottky as S

VoltageOnMetal = 0 # V
TimeStep = 0.1 # seconds
GridPointsNumber = 5001 # per u.DeviceLength = 5µm => 10 grid points p = 0.1
DriftAmount0 = 1

#-------------------------------

def Init_Sample(Bias, dt, nx):
	print("A Sample has been created.")
	return S.Schottky(Bias, dt, nx)

def OxygenDrift(Sample, Iterations, DriftAmount):
	for i in range(Iterations):
		Sample.DopingDrift(DriftAmount)
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
	
def SetDepletionWidth(Sample, Width):
	#Changes the initial Charge Distribution
	Sample.W = Width
	Sample.ChargeDensity[:] = 0.0
	Sample.ChargeDensity[:int(Sample.W/Sample.dx)] = Sample.DopingDensity[:int(Sample.W/Sample.dx)]
	solverLU(Sample)
	print("Change of depletion Region length successful.")
	return
	
def AddVoltageTrace(Vstart, Vend, Steps = 0, V = []):
	# If no V-array is given, a new V-array is created
	for i in range(Steps):
		V = V+[Vstart + (Vend - Vstart)*i/Steps]
	print("Trace", Vstart, " ", Vend, "has been added to the Voltage trace.")
	return V

#--------------------------------
def Run(Bias, dt, nx):
	# Sample Initialization
	Sample = Init_Sample(Bias,dt,nx)
	
	# Read Sample Data from a File
	Sample.Import("V=-30(02)T=43201s")

	print("Initial Depletion Region Width = ", Sample.W)
	
	
	#Sample.DopingDensity = zeros(Sample.nx)
	#Sample.DopingDensity[1:] = ones(Sample.nx-1)*0.000001
	#Sample.DopingDensity[0] = 0.5001
	#SetDepletionWidth(Sample, 5000)
	
	# Plot Initial Conditions
	Sample.PlotResults()
	
	
	# Create a Voltage Function
	
	#V = AddVoltageTrace(-30, -30, 1)
	
	V = AddVoltageTrace(0, 10, 100)
	V = AddVoltageTrace(10, 0, 100, V)
	V = AddVoltageTrace(0, -10, 100, V)
	V = AddVoltageTrace(-10, 0, 100, V)
	V = AddVoltageTrace(0,0,1,V)
	
	#print(V)
#-------------------------------
	T = [0]
	R = [sum(Resistivity(Sample)*Sample.dx)]
	I = [Current(Sample)]
	v = [V[0]]
	DepletionRegionWidth = [Sample.W]
	
	# Drift steps
	DriftIterations = 1
	SmoothingSteps = 0

	print("Start of Simulation Process")
	for i in range(len(V)) :#fix20 = 20mal Bilder Anzeigen
		Sample.Bias = V[i]
		print("Current Bias =", Sample.Bias, "V")
	
		for j in range(DriftIterations) :
			# Time actualization
			T = T +[T[-1]+1]
			v = v + [Sample.Bias]
			#print("Drift", j)
			Reload(Sample)
			OxygenDrift(Sample, 1, DriftAmount0)			
			for k in range(SmoothingSteps): Smoothing(Sample)
			Reload(Sample)
			
			R = R + [sum(Resistivity(Sample))*Sample.dx]
			print(sum(Resistivity(Sample))*Sample.dx,  min(abs(Sample.E)), max(abs(Sample.E)))#[int(Sample.W/Sample.dx):])*Sample.dx)#,Current(Sample)
			I = I + [Current(Sample)]
			DepletionRegionWidth = DepletionRegionWidth + [Sample.W]
			Smoothing(Sample, s0 = .5)
			
		print(int(T[-1]), sum(Sample.DopingDensity))
		print("Time = ",T[-1], "s, 'Drift' Method executed", (i+1)*DriftIterations, " times")
		print("Depletion region width = ", Sample.W, " nm")
		if(((len(T)-1)%100) < 0.01 ):Sample.PlotResults()
		#print(Sample.Phi[0])
	print("End of Simulation Process")
	Sample.PlotResults()

	# Save Results
	savetxt("Resistance10V(05)", (T, v, R, I, DepletionRegionWidth), newline="\n")
	WriteToFile(Sample, "V=10(05)T=201s")

	plt.show()
	return 0
	
if __name__ == '__main__':
	Run(VoltageOnMetal, TimeStep, GridPointsNumber)
