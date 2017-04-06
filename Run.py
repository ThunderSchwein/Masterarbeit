from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *
from DopingDrift import DriftRK, Diffusion
from ChargeRedistribution import *
import Schottky as S

VoltageOnMetal = 0 # V
TimeStep = 1 # seconds
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
	#Import(Sample, "V=0(1)T=43200s")
	#Import(Sample, "V=30(4)T=18000s")
	Import(Sample, "V=0(5)T=3600s")
	
	# Plot Initial Conditions
	Sample.PlotResults()

	V = AddVoltageTrace(10, 10, 100)
	#V = AddVoltageTrace(0, 10, 100)
	#V = AddVoltageTrace(10, -10, 200, V)
	#V = AddVoltageTrace(-10, 0, 100, V)
	#V = AddVoltageTrace(0,0,1,V)
	#print(V)
	Sample.W = u.DeviceLength/10
#-------------------------------
	T = [64800]
	R = [sum(TotalResistance(Sample))*Sample.dx]
	C = [sum(Sample.ChargeDensity)*Sample.dx]
	
	# Drift paired with Diffusion
	DiffusionSteps = 0
	DriftIterations = 10

	print("Start of Simulation Process")
	for i in range(len(V)) :#fix20 = 20mal Bilder Anzeigen
		Sample.Bias = V[i]
		print("Current Bias =", Sample.Bias, "V")
	
		for j in range(DriftIterations) :
			# Time actualization
			T = T +[T[-1]+dt]
			
			print("Drift")
			Sample.Phi[0] = Sample.Bias + u.BarrierHeight
			ReloadDeltaPhi(Sample)
			OxygenDrift(Sample, 1)
			
			#print("Diffusing")
			OxygenDiff(Sample, DiffusionSteps, 1)
			#for h in range(DiffusionSteps) :
				#Sample.DopingDensity[:20000] = Diffusion(Sample.X[:20000], Sample.dt , Sample.DopingDensity[:20000], 1)
			#ReloadDeltaPhi(Sample)
			
			#Resistance Calculation
			C = C + [sum(Sample.ChargeDensity)*Sample.dx]
			R = R + [sum(TotalResistance(Sample))*Sample.dx]
			
		print(T[-1], sum(Sample.DopingDensity))
		print("Time = ",T[-1], "s, 'Drift' Method executed", (i+1)*DriftIterations, " times")
		print("Depletion region width = ", Sample.W, " nm")
		if((V[i]%1)< 0.001) : Sample.PlotResults()
		#Sample.PlotResults()
		ReloadDeltaPhi(Sample)
	print("End of Simulation Process")
	Sample.PlotResults()

	# Save Results
	savetxt("Resistance10V(0)", (T, R, C), newline="\n")
	WriteToFile(Sample, "V=10(5)T=66000s")

	plt.show()
	return 0
	
if __name__ == '__main__':
	Run(VoltageOnMetal, TimeStep, GridPointsNumber)
