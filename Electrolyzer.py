from cmath import asinh, sqrt
import numpy as np
from numpy import log as ln
import math

#Dummy Variable
T0 = 273
Trel = T0+90
CD = 2000

#Operating Condition
Pops = 101325   #Operating Pressure (Pa)
Tamb = 298      #Ambient Temperature (K)

#Defining Thermodynamics Parameters
F = 96485   #Faraday's constant (C/mol)
n = 2       #Electron number from complete reaction
R = 8.314   #Gas Constant (J/molK)

#Defining Electrolyzer Parameter
AlphaAn = 0.5   #Charge Transfer Coefficient for Anode (-)
AlphaCat = 1    #Charge Transfer Coefficient for Cathode (-)
#CDAn = 1e-13    #Limiting Current Density on Anode (A/m^2) not used
#CDCat = 1e-16   #Limiting Current Density on Cathode (A/m^2) not used
CDLim= 20000    #Limiting Current Density Global (A/m^2)
Deltam = 5e-05  #Thickness of membrane (m)
WA = 1          #Water Activity, Pw/Psat or 1 if full humidification
Beta2 = 2       #Concentration overvoltage coefficient

#Calculating Vrev
PH2O = Pops
PO2 = Pops
PH2 = Pops
Vrev = 1.229-8.5*1e-04*(Trel-Tamb)+4.3085*1e-05*Trel*ln(PH2*sqrt(PO2)/PH2O)
#print("Vrev: ",Vrev.real)

#Calculating activation overpotential
CD0 = 1.08*1e-17*math.exp(0.086*Trel)
print("CD0: ",CD0)
Vact = ((AlphaAn+AlphaCat)/(AlphaAn*AlphaCat))*R*Trel/(n*F)*ln(CD/CD0)
print("alpha: ",(AlphaAn+AlphaCat)/(AlphaAn*AlphaCat))
print("ln CD/CD0: ",ln(CD/CD0))
print("RT/2F: ",R*Trel/(n*F))
print("Vact: ",Vact)
#Vact2 = R*Trel/(0.42*F)*asinh(CD/(2*1e-02))+R*Trel/(0.42*F)*asinh(CD/(2*1e-02))
#print("Vact2: ",Vact2)

#Calculating ohmic overpotential
lambdam = 0.043+17.81*WA-39.85*pow(WA,2)+36*pow(WA,3)
Sigmam = (0.005139*lambdam-0.00326)*math.exp(1268*(1/303-1/Tamb))
Vohm = CD*Deltam/Sigmam
#print("Vohm: ",Vohm)

#Calculating ohmic overpotential
Pow10 = 2.8206+0.02953*(Trel-T0)-9.1837*1e-05*pow((Trel-T0),2)+1.4454*1e-07*pow((Trel-T0),3)
Psat = 1.01325*(10**Pow10)
#print("Psat: ",Psat)
Px = PO2/(0.1173*Pops)+Psat/Pops
#print("Px: ",Px)
Beta1 = (8.66*1e-5*Trel-0.068)*Px-1.6*1e-04*Trel+0.54
Vcon = CD*(Beta1*CD/CDLim)**Beta2
#print("Vcon: ",Vcon)

#Overall voltage
Vovr = Vrev.real+Vact/10+Vohm/10+Vcon/10
print("Vovr: ",Vovr)
