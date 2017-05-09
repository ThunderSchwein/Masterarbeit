# Konstantin Murasov
"""
Solver.py contains a set of tools for the class "Schottky.py"
for solving the Poisson equation on the date, saved in 

"""

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u

def E_Field(Sample):
	E = Sample.E
	Phi = Sample.Phi
	dx = Sample.dx
	
	for i in range(Sample.nx-2) :
		E[i+1] = (Phi[i] - Phi[i+2])/(2*dx)
	E[0] = (Phi[0] - Phi[1])/dx
	E[-1] = (Phi[-2] - Phi[-1])/dx
	
	return E
	
def solverLU(Sample):
	dx = Sample.dx
	nx = Sample.nx
	Phi = Sample.Phi
	Rho = Sample.ChargeDensity
	
	#print("Solver")
	
	# Solution based only on present unbalanced Charge

	Phi[-1] = Phi[-2] = 0	
	for j in range(nx-2):
		Phi[nx-3-j] = 2*Phi[nx-j-2] - Phi[nx-j-1] - Rho[nx-j-2]*dx**2/Sample.Dieelek[nx-j-2]
		
	if(Sample.Bias > abs(u.BarrierHeight)):
		"""
		Phi2 = zeros(nx)
		Phi2[0] = Sample.Bias + u.BarrierHeight - Phi[0]
		R = Resistivity(Sample)
		Rs = sum(R)
		for i in range(nx -1):
			Phi2[i+1] = Phi2[i] - R[i+1]/Rs*Phi2[0] 
		#print(Phi2[-1])
		Sample.Phi = (Phi + Phi2)
		
		"""
		Sample.Phi[0] = Sample.Bias + u.BarrierHeight
		# Inversion
		R = Resistivity(Sample)
		Rs = sum(R)
		#print(Rs, Sample.dx, u.ds)
		for i in range(nx -1):
			Sample.Phi[i+1] = Sample.Phi[i] - R[i+1]/Rs*Phi[0]
		
		
	if(Sample.W >= u.DeviceLength):
		# Depletion region reaches the size of the sample.
		Phi2 = zeros(nx)
		Phi2[0] = Sample.Bias + u.BarrierHeight - Phi[0]
		if(Phi2[0] > 0 ) : return
		R = Resistivity(Sample)
		Rs = sum(R)
		for i in range(nx -1):
			Phi2[i+1] = Phi2[i] - R[i+1]/Rs*Phi2[0] 
		#print(Phi2[-1])
		Sample.Phi = (Phi + Phi2)
	
	#print(Sample.W)
	return
	
def Current(Sample):
	n = 4.4
	A = 1.08e7
	i_0 = 1e-7
	i   = 1e-7
	v = Sample.Bias + u.BarrierHeight
	KbT = u.Kb*u.T
	R = sum(Resistivity(Sample))*Sample.dx
	
	#g = 10
	#Phi2 = zeros(g) 
	#for j in range(g):
	#	Phi2[g-3-j] = 2*Phi2[g-j-2] - Phi2[g-j-1] - 2*Sample.DopingDensity[g-j-2]*Sample.dx**2/Sample.Dieelek[g-j-2]
		
	Barrier = u.BarrierHeight# - abs(Phi2[0])
	#print(Barrier, sum(Sample.DopingDensity[0:5]), Phi2[0])
	
	if(Sample.Bias <= abs(Barrier)): return 0
	
	#if(Sample.Bias <= abs(u.BarrierHeight)): return 0
	
	i_0 = A * (exp((v/(2*n))/KbT)- 1 )*exp(Barrier/(KbT))
	i = A * (exp((v/n - i_0*R)/KbT) -1 )*exp(Barrier/(KbT))
	while(abs(i-i_0) > (0.001*(min(i, i_0)))):
		if(i>i_0):  i_0 = i_0 * 1.003
		if(i<i_0):  i_0 = i_0 / 1.005
		i = A*(exp((v/n - i_0*R)/KbT) -1 )*exp(Barrier/(KbT))
		
	return i	
	
def Resistivity(Sample):
	DP = Sample.DopingDensity
	r = zeros(Sample.nx)
	x = 1e-7
	for i in range(Sample.nx):
		#if(DP[i] == 0) :	r[i] = 20000/Sample.dx
		#else :			
		r[i] = 20000*(x/(x+DP[i]*Sample.dx))/Sample.dx
	return r
#-----------------------------------------------------------
