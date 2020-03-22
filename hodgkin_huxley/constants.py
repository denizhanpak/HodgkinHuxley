import numpy as np

dt = 0.001

#MODEL CONSTANTS
Cm = 1.0 #uF/cm^3 : Mean = 0.91, Range = 0.8-1.5
VNa = 55#-115 #mV : Mean = -109, Range = -95 to -119
VK = -77#-72#12 # mV : Mean = 11. Range = 9-14
VL = -65#-49#10.613 # mV : Mean = -11, Range = -4 to -22
gNa = 40#120 #mS/cm^2 : Mean = 80/160, Range = 65-90/120-260
gK = 35#36 #mS/cm^2 : Mean = 34, Range = 26-49
gL = 0.3 #mS/cm^2 : DEfault = 0.3, Range = 0-3, 0-26, 13-50
gC = 0.06 #mS/cm^2 : Default = 0.06, Range = 0-0.33
tau_r = 0.5 #ms : characteristic rise time
tau_d = 8 #ms : characteristic decay time
Vsyn = 20 #mV
V0 = -20 #mV

def alpha_n(Vi):
    #return 0.01*(Vi + 50)/(1 - np.exp(-0.1*(Vi + 50)))
    return 0.02*(Vi - 25)/(1 - np.exp(-(Vi - 25)/9))

def beta_n(Vi):
    #return 0.125*np.exp((Vi+60)/80)
    return -0.002*(Vi - 25)/(1 - np.exp((Vi - 25)/9))

def alpha_m(Vi):
    #return 0.1*(Vi + 35)/(1 - np.exp(-0.1*(Vi + 35)))
    return 0.182*(Vi + 35)/(1 - np.exp(-(Vi + 35)/9))

def beta_m(Vi):
    #return 4*np.exp((Vi + 60)/18)
    return -0.124*(Vi + 35)/(1 - np.exp((Vi + 35)/9))

def alpha_h(Vi):
    #return 0.07*np.exp((Vi + 60)/20)
    return 0.25*np.exp(-(Vi + 90)/12)

def beta_h(Vi):
    #return 1/(1 + np.exp(-0.1*(Vi +30)))
    return 0.25*np.exp((Vi + 62)/6)/np.exp((Vi + 90)/12)
