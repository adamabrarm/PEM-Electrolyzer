import numpy as np
from numpy import log as ln
import math

#Input Variable
T0 = 273        #Zero Degree Celsius to Kelvin (K)
Tcell = T0+60   #Cell Operating Temperature
NH2O = 1.5        #Flow water input to Electrolyzer (mol/s)

#Operating Condition
Pops = 101325   #Operating Pressure (Pa)
Tamb = 298      #Ambient Temperature (K)

#Defining Thermodynamics Parameters
F = 96485   #Faraday's constant (C/mol)
n = 2       #Electron number from complete reaction
R = 8.314   #Gas Constant (J/molK)
CpH2O = 4.19    #Specific Heat Capacity of Water (J/gK)
MwH2O = 18  #Molecular Weight of H2O (g/mol)
MwO2 = 32   #Molecular weight of O2 (g/mol)
MwH2 = 2    #Molecular weight of H2 (g/mol)

#Mass Balance around electrolyzer
NH2 = NH2O      #All water turned into H2 and O2
MH2 = NH2*MwH2     #Mass flow of H2 (g/s)
NO2 = NH2O/2    #All water turned into H2 and O2
MO2 = NO2*MwO2    #Mass flow of O2 (g/s)
Istack = NH2O*2*F 
print("Istack: ",Istack)
Ncell = 1000    #Number of cell per stack
Nstack = 100    #Number of stack needed
Icell = Istack/Ncell/Nstack #Current per cell

#Calculating Voltage reversible & Overpotential
Vrev = 1.5148-1.5421*1e-03*Tcell+9.523*1e-05*Tcell*ln(Tcell)+9.84*10e-08*Tcell**2   #Voltage reversible calculation
print("Vrev: ",Vrev)
Vact = 0.0514*Icell+0.2798  #Calculating activation overpotential
print("Vact: ",Vact)
Vohm = 0.09*Icell           #Calculating ohmic overpotential
print("Vohm: ",Vohm)
Vovr = Vrev+Vact+Vohm       #Overall voltage
print("Vovr: ",Vovr)
dV = Vact+Vohm              #Overpotential Difference
CellEff = Vrev/Vovr*100     #Cell Efficiency Calculation
print("Cell Efficiency: ",CellEff,"%")

#Heat Balance for Outlet Streams
Qheat = Icell*dV              #Waste heat calculation (W)
print(Qheat)
dT = Qheat/(CpH2O)/(NH2O*MwH2O/Ncell)  #Temperature raise by waste heat (K)
Tout = Tcell+dT
print(Tout-T0)





