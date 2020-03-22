import sys
from network import *
    
if __name__=='__main__':
    #My suspicions (proven wrong) are for getting complex patterns: Need 
    #clustering > 0.1, Inhibitory neurons must have out-degree > 0.15*netsize
    #randomize initial voltages?
    #CPG like structure?
    #Take away r_ij?
    #Add stimulating current?

    if len(sys.argv) >= 4 or "-h" in sys.argv or "--help" in sys.argv:
        print("usage: hodgkin_huxley.py [NETSIZE=25] [DURATION=200]\nNETSIZE will be rounded up to nearest square.")
        exit()
    if len(sys.argv) >= 1:
        netsize = 25
        duration = 200
        if len(sys.argv) >= 2:
            netsize = sys.argv[1]
            if len(sys.argv) == 3:
                duration = sys.argv[2]

    print("Building network of size %s for a duration of %s" % (netsize, duration))
    hh = HodgkinHuxley(netsize, duration)
    hh.initialize(1, 3, 2.5, 2)
    hh.simulate()
    hh.showOutput()
    #hh.show_r()
    #show_x0()
