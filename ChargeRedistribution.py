# Konstantin Murasov
# ChargeRedistribution.py
""" 
This Module calculates the Distribution of charges, which
results from doping migration in a TiO2-Sample, with Data
provided by "Schottky.py", until the solution of the 
Poisson-equation is conform with the boundary conditions
of the Sample:
Potential(x = 0) = Poisson (x = 0) + BarrierHeight

Please take note that "ReloadDeltaPhi"-Method utilizes a
try-and-error approach. to find the right width of the 
depletion region.
"""

from numpy import *
import matplotlib.pyplot as plt
from copy import deepcopy
import Units as u
from SolverLU import *


def Reload(Sample):
	# In case, there are multiple methods for calcualtion of
	# the charge density, they can be switched via this function
	# by commenting the unwanted one out.

	Reload0(Sample)
	return
	
def Reload0(Sample) :
	# This method is redistributing the electrons, until the
	# solution of the poisson eqaution is conform with the
	# Sample boundary conditions.
	
	# The measure for the quality of the solution is the difference
	# between the boundary value at x = 0 calculated by "Solver.py"
	# and the boundary valie dictated by the applied bias:
	# Sample.Phi[0]>(u.BarrierHeight+Sample.Bias
	# The usage of factors (u.BarrierHeight*0.99875+Sample.Bias)
	# allows to control the precision of the calculation, deeming 
	# the values outside this range as not sufficently precise.
	
	# print() - comments can be enabled to monitor the changes in the
	# depletion region width value and possibly spot an error.
	
	# shortcut to the length of the X-Axis
	nx = Sample.nx
	
	# Resetting the Charge distribution
	Sample.ChargeDensity = zeros(nx)
	
	
	if(Sample.Bias < -u.BarrierHeight):
		# Case 1 : the depletion region is present.

		
		if((Sample.Bias > -u.BarrierHeight*0.99) & (Sample.Bias < -u.BarrierHeight*1.01)) :
					# Case 1.1 :  Sample.Bias close to -u.BarrierHeight :
					# In this case no depletion region is formed and thus
					# Sample.ChargeDensity = 0 for all x.
					
					Sample.W = 1 	# (Never set to zero !)
					# Flat potential, as the bias and the built-in potential (BarrierHeight) cancel each other.
					Sample.Phi = zeros(nx)
					# No unbalanced charges in the Sample.
					Sample.ChargeDensity = zeros(nx)
				
		if(Sample.W >= u.DeviceLength):
			# Case 1.2 : The width of the depletion egon excedds the size of the Sample (DeviceLength)
			
			Sample.ChargeDensity = Sample.DopingDensity*u.Titanium_DopingCharge
			Sample.W = u.DeviceLength
			# Solver takes care of the arising potential.
			solverLU(Sample)
		else:
			# Case 1.3 : Stepwise improvement of the value of the depletion region width.
			# The factors rule the enlargement/minimizing of this width.
			
			factor1 = 1.05
			factor2 = 1.01
			factor3 = 1.002
			factor4 = 1.0004
			
			while ((Sample.Phi[0]>(u.BarrierHeight+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				# If the initially guessed depletion region width is too low, it must be expanded.
				Sample.W = Sample.W*factor1
				#Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				# Solver calculates the resulting potential.
				solverLU(Sample)
				#print(1, Sample.W, Sample.Phi[0])
				
			while (Sample.Phi[0]<(u.BarrierHeight+Sample.Bias)):
				# In the opposite case it must shrink. Do not forget to reset the Sample.ChargeDensity to zero before the calculation.
				Sample.W = Sample.W / factor2
				Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				# Solver calculates the resulting potential.
				solverLU(Sample)
				#print(2, Sample.W, Sample.Phi[0])
				
			while ((Sample.Phi[0]>(u.BarrierHeight*0.99875+Sample.Bias)) & (Sample.W < u.DeviceLength)):
				# Once we are close to finding the value of the depletion region, we must use smaller steps,
				# as well as intoduce the tolerance u.BarrierHeight*0.00125 ( = 0.125 %).
				Sample.W = Sample.W*factor3
				#Sample.ChargeDensity[:] = 0.0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				# Solver calculates the resulting potential.
				solverLU(Sample)
				#print(3, Sample.W, Sample.Phi[0])
				
			while (Sample.Phi[0]<(u.BarrierHeight*1.00125+Sample.Bias)):
				# Slow approach to the final value of the depletion region width. Do not for get to reset Sample.ChargeDensity to zero.
				Sample.W = Sample.W /factor4
				Sample.ChargeDensity[:] = 0
				Sample.ChargeDensity[:int(Sample.W/u.DeviceLength*nx)] = Sample.DopingDensity[:int(Sample.W/u.DeviceLength*nx)]*u.Titanium_DopingCharge
				# Solver calculates the resulting potential.
				solverLU(Sample)
				#print(4, Sample.W, Sample.Phi[0])
			#del DepletionRegionCharge
	
	else:
		# Case 2 : The applied Bias is greater then the built-in potential (BarrierHeight)
		# All of the doping ions are unionized, as the depletion region is annihilated.
		B = 5
		Sample.W = B
		Sample.Phi[0] = Sample.Bias + u.BarrierHeight
		Sample.ChargeDensity[:B] = Sample.DopingDensity[:B]*u.Titanium_DopingCharge
		
	# Solver calculates the resulting potential.
	solverLU(Sample)

	#print(0, Sample.W, Sample.Phi[0])
	# Calculation of the resulting electric field
	Sample.E = E_Field(Sample)

	return
