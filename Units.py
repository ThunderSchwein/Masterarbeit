# KonstantinMurasov
"""
This Module contains the constants needed to calculate the Potential 
and the electric field in a given Sample.

Please take note that the default length unit (ds) is 1 nm.
"""
from numpy import *

# Length Unit [m]
ds = 1e-9	# 1 nm
# Tmperature
T  = 300	# K
# Physical Constants
Kb = 8.617e-5
# Electron Charge [Coulomb]
e =-1.6021766208e-19
# Vacuum permittivity 8.854e-12 As/Vm
e0 = 55.263e-3 # (e-)/(V*nm) 

# Titanum Dioxide Properties

# Metal Properties 
# [Pt]
Platinum_WorkFunktion = -5.2 #eV [ @ 5.0 - 5.5 eV ]

# Semiconductor Properties
# [TiO2]
# Bandgap
Titanium_Bandgap = 3 	#[ev]
# Intrinsic carrier concentration
Titanium_IntrinsicCarrierConcentration = 0.0374*exp(-Titanium_Bandgap/(Kb*T))
# Electron mobility [0.2 cm^2/Vs]
Titanium_ElectronMobility = 0.2e14
Titanium_HoleMobility = 0
# Electron affinity
Titanim_ElectronAffinity = -4.0 	#eV
# Dielectricity constant factor
Titanium_DielectricityFactor = 70	#[ @ 86 - 173 ]
# Doping density [1/nm**3]
Titanium_DopingDensity = 2.2e-4	#[1e23 /m**3]
Titanium_DopingCharge = 2
# Oxygen Vacancy Mobility
Titanium_DopingMobility = 1e4	#[nm**2 / V/s]

# Pt/TiO2 barrierheight
BarrierHeight = -1.1 #eV

#------------------------------------------
# Strungaru Paper

# Filament Length [nm]
DeviceLength = 5000	# 5µm

# Oxygen Vacancy Drift Mobility [nm**2 / Vs]
#Mobility = 1000	# 1e-16 m**2 / Vs
#Mobility = 1e6	# 1e-16 m**2 / Vs

# Electroforming break-even Voltage [V/nm]
#E_min = 1.6e-5
E_min = 0.0
#E_min = 4e-4 	# 4e5 V/m

# Conductivity [Simens = 1/R = Ampere/nm]
Sigma1 = 3.5e-3
Sigma2 = 6e-3
Sigma3 = 1e-4
