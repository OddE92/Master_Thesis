import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import numpy as np
import itertools as it

basename = "Data/Compare/Ee"
contname = "_Bpercent"

Earr = [1, 10, 100, 1000]

# Counters: [0] = +-0.5R_L, [1] = +-2.5, [2] = +-5, [3] = +- 10
count   =  np.zeros((4,5,4), dtype=int)


l = 0
for E in range (15, 19):

    k = 0

    for percent in ["1", "10", "25", "50", "100"]:    

        strGCT      = basename + str(E) + contname + percent + "_GCT.dat"
        strExact    = basename + str(E) + contname + percent + ".dat"

        dataGCT     = np.loadtxt(strGCT)
        dataExact   = np.loadtxt(strExact)

        tempdr = np.sqrt( (dataGCT[:,0]-dataExact[:,0])**2 + (dataGCT[:,1]-dataExact[:,1])**2 + (dataGCT[:,2]-dataExact[:,2])**2 )
       
        tempdr = abs(tempdr - dataGCT[:, 3])

        for x in tempdr:
            if x <= 0.5:
                count[l,k,0] += 1
            if x <= 2.5  and x > 0.5:
                count[l,k,1] += 1
            if x <= 5.0  and x > 2.5:
                count[l,k,2] += 1
            if x <= 10.0 and x > 5.0:
                count[l,k,3] += 1

        k += 1
    l += 1


np.savetxt('Data/Particle_count/E15.dat', count[0])
np.savetxt('Data/Particle_count/E16.dat', count[1])
np.savetxt('Data/Particle_count/E17.dat', count[2])
np.savetxt('Data/Particle_count/E18.dat', count[3])


            

