#Length Unit [m]
ds = 1e-9	# 1 nm
T  = 300	# K

# Physical Constants
Kb = 8.617e-6

# Electron Charge [Coulomb]
e =-1.6021766208e-19

# Vacuum permittivity 8.854e-12 As/Vm
e0 = 55.263e-3 # (e-)/(V*nm) 

# Titanum Dioxide Properties

# Metal Properties 
# [Pt]
Platinum_WorkFunktion = 5.2 #V [ @ 5.0 - 5.5 eV ]

# Pt/TiO2 barrierheight
BarrierHeight = 1.9 #eV

# Semiconductor Properties
# [TiO2]
# Electron affinity
Titanim_ElectronAffinity = 4.0 	#V
# Dielectricity constant factor
Titanium_DielectricityFactor = 70	#[ @ 86 - 173 ]
# Doping density [1/nm**3]
Titanium_DopingDensity = 1e-4	#[1e23 /m**3]
Titanium_DopingCharge = 2
# Oxygen Vacancy Mobility
Titanium_DopingMobility = 2.5 	#[nm**2 / V/s]

#------------------------------------------
# Strungaru Paper

# Filament Length [nm]
DeviceLength = 5000	# 5µm

# Oxygen Vacancy Drift Mobility [nm**2 / Vs]
Mobility = 1000		# 1e-15 m**2 / VS
# Electroforming break-even Voltage [V/nm]
E_min = 4e-4 	# 4e5 V/m

# Conductivity [Simens = 1/R = Ampere/nm]
Sigma1 = 3.5e-3
Sigma2 = 6e-3
Sigma3 = 1e-3
