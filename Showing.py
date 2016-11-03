from numpy import *
import matplotlib.pyplot as plt

X, Rho, Phi = loadtxt('Example.txt', unpack = True)

plt.figure("Potential")
plt.plot(X,Rho)
plt.plot(X,Phi, color = "G")
#plt.plot(X, X**3, color = "R")
#plt.plot(X, X**4, color = "R")
#plt.plot(X,-sin(2*pi*X)/(4*pi**2), color = "R")
plt.plot(X, exp(X), color = "R")
plt.show()
