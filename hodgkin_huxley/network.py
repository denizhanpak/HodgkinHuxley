import random
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.signal import find_peaks
from .constants import *
#from hodgkin_huxley import visualizer

def x0(V, gating): #The equilibrium functions x0 for gating variable x=m,n,h
    if gating=='n':
        return alpha_n(V)/(alpha_n(V) + beta_n(V))
    elif gating=='m':
        return alpha_m(V)/(alpha_m(V) + beta_m(V))
    elif gating=='h':
        return alpha_h(V)/(alpha_h(V) + beta_h(V))
    
def show_x0():
    Vi = np.arange(-100,100,0.1)
    n0 = x0(Vi, 'n')
    m0 = x0(Vi, 'm')
    h0 = x0(Vi, 'h')
    plt.figure()
    plt.plot(Vi, n0, label='n0')
    plt.plot(Vi, m0, label='m0')
    plt.plot(Vi, h0, label='h0')
    plt.axvline(x=-63, color='r', linestyle='--', label='Resting potential')
    plt.title("Equilibrium functions of the gating variable")
    plt.xlabel("Membrane voltage (mV)")
    plt.legend()
    plt.show()

def adjmat(size, p=1, q=1, r=2, dim=2, show=True): #default is p=1 and r=2
    directed_sw = nx.navigable_small_world_graph(size, p, q, r, dim) #A navigable small-world graph is a directed grid with additional long-range connections that are chosen randomly.
    #visualizer.colored_circle(directed_sw)
    aspl = nx.average_shortest_path_length(directed_sw)
    acc = nx.average_clustering(directed_sw)
    print("Average shortest path length: " + str(aspl))
    print("Average clustering coefficient: " + str(acc))
    aij = nx.to_numpy_matrix(directed_sw)
    #ADDING INHIBTORY NEURONS
    inh = random.sample(range(size**2), k=math.ceil(0.15 * size**2)) #15% OF NEURONS ARE INHIBITORY
    for i in inh:
        aij[i, :] *= -1
    if show:
        plt.matshow(aij)
    return aij, directed_sw #Returns the graph adjacency matrix as a NumPy matrix.
    


#This class defines the network. It is initialized with the size of the networks and ?
#The nodes are individual HodgkinHuxley neurons where the value of each neuron is stored as an element of a vector.
class HodgkinHuxley():
    def __init__(self, size, time):
        self.size = size
        self.time = time
        self.I = np.zeros(size) #External input current
        self.V = np.zeros(size) #Membrane voltage
        self.n = np.zeros(size) #K activation var.
        self.m = np.zeros(size) #NA activation var.
        self.h = np.zeros(size) #Na inactivation var.
        self.r = np.zeros(size) #Fraction of bond receptors with presynaptic neurons
        self.a = np.zeros((size,size)) #Synaptic adjacency matrix
        self.voltage = np.zeros((size, int(time/dt)))
        self.rt = np.zeros((size, int(time/dt)))
        #ADD TIME RECORDINGS
        
    def initialize(self, p=1, q=3, r=2.5, dim=2):
        self.V = 100*np.random.rand(self.size) - 80  # -63*np.ones(self.size)
        self.n = x0(self.V, 'n')
        self.m = x0(self.V, 'm')
        self.h = x0(self.V, 'h')
        self.r = np.random.rand(self.size) # (self.V + 80)/100
        self.a, self.network = adjmat(int(np.sqrt(self.size)), p, q, r, dim)
    
    def step(self, i):
        #RECORD CURRENT STATE
        self.voltage[:, i] = self.V
        self.rt[:, i] = self.r
        
        INa = gNa*(self.m**3)*self.h*(self.V - VNa)
        IK = gK*(self.n**4)*(self.V - VK)
        IL = gL*(self.V - VL)
        
        drdt = (1/tau_r - 1/tau_d)*(1 - self.r)/(1 + np.exp(-self.V + V0)) - self.r/tau_d
        self.r = drdt*dt + self.r
#        weight = self.r*self.a
#        weight = np.asarray(weight)
#        weight = weight.reshape((self.size,))
#        Isyn = gC*weight*(Vsyn - self.V)

        Isyn = gC*(Vsyn - self.V)*self.r*self.a
        Isyn = np.asarray(Isyn)
        Isyn = Isyn.reshape((self.size,))

        dVdt = (-INa -IK -IL + Isyn)/Cm
        dndt = -(alpha_n(self.V) + beta_n(self.V))*self.n + alpha_n(self.V)
        dmdt = -(alpha_m(self.V) + beta_m(self.V))*self.m + alpha_m(self.V)
        dhdt = -(alpha_h(self.V) + beta_h(self.V))*self.h + alpha_h(self.V)

        self.V = dVdt*dt + self.V
        self.n = dndt*dt + self.n
        self.m = dmdt*dt + self.m
        self.h = dhdt*dt + self.h

    def simulate(self):
        for t in range(int(self.time/dt)):
            self.step(t)

    def trace(self):
        t = np.arange(0,self.time,dt)
        plt.figure()
        plt.subplot(411)
        plt.plot(t, self.voltage[0,:])
        plt.title("Neuron 1")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.subplot(412)
        plt.plot(t, self.voltage[1,:])
        plt.title("Neuron 2")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.subplot(413)
        plt.plot(t, self.voltage[2,:])
        plt.title("Neuron 3")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.subplot(414)
        plt.plot(t, self.voltage[3,:])
        plt.title("Neuron 4")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.suptitle("Neural traces")
        plt.tight_layout()
        plt.show()

    def show_r(self):
        plt.figure()
        t = np.arange(0,self.time,dt)
        for r in self.rt:
            plt.plot(t, r)
        plt.title("r_ij(t)")
        plt.xlabel("Time (ms)")
        plt.show()

    def showvoltage(self):
        plt.figure()
        t = np.arange(0,self.time,dt)
        for v in self.voltage:
            plt.plot(t, v)
        plt.title("V_i(t)")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.show()

#    def findPeaks(self):
#
#
#    def raster(self):

