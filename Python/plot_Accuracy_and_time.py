import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

basename = "Data/Compare/Ee17_Bpercent"

dr = np.empty((0,1000), float)             # Difference in position
R_larmor = np.empty((0,1000), float)
dt = np.empty((0,1000), float)
posGCT = np.empty((0,1000), float)
posExact = np.empty((0,1000), float)
drR_larmor = np.empty((0,1000), float)

drR_larmorAverage = []
drAverage = []
R_larmorAverage = []
dtAverage = []

n = np.linspace(1, 1000, num=1000)

print dr
for percent in ["1", "10", "25", "50", "100"]:

    strGCT      = basename + percent + "_GCT.dat"
    strExact    = basename + percent + ".dat"

    dataGCT     = np.loadtxt(strGCT)
    dataExact   = np.loadtxt(strExact)

    temp1   = np.sqrt( dataGCT[:,0]**2 + dataGCT[:,1]**2 + dataGCT[:,2]**2 )
    temp2   = np.sqrt( dataExact[:,0]**2 + dataExact[:,1]**2 + dataExact[:,2]**2 )
    tempdr  = np.sqrt( (dataGCT[:,0]-dataExact[:,0])**2 + (dataGCT[:,1]-dataExact[:,1])**2 + (dataGCT[:,2]-dataExact[:,2])**2 )
    tempdt  = dataExact[:,3] / dataGCT[:,4]

    dr          = np.append(dr, [tempdr], axis=0 )
    R_larmor    = np.append(R_larmor, [dataGCT[:, 3]], axis=0)
    dt          = np.append(dt, [tempdt], axis=0)
    posGCT      = np.append(posGCT, [temp1], axis=0)
    posExact    = np.append(posExact, [temp2], axis=0)

    drR_larmor = np.append(drR_larmor, [tempdr/dataGCT[:,3]], axis=0)
    drR_larmorAverage = np.append(drR_larmorAverage, np.mean(tempdr/dataGCT[:,3]))

    drAverage = np.append(drAverage, np.mean(tempdr))
    R_larmorAverage = np.append(R_larmorAverage, np.mean(dataGCT[:,3]))
    dtAverage = np.append(dtAverage, np.mean(tempdt))



####### Figure 1 = position compared to R_L
fig1 = plt.figure(1)
ax = fig1.gca()

ax.plot(n, dr[0,:], 'x')
ax.plot([1, 1000], [drAverage[0], drAverage[0]], '-', color='green')
ax.plot(n, R_larmor[0,:], '-', color='red')
ax.plot([1, 1000], [R_larmorAverage[0], R_larmorAverage[0]], '-', color='black')


####### Figure 2 = dr/R_L with average
fig2 = plt.figure(2)
ax2 = fig2.gca()
ax2.set_title(r"Amplitude of $\vec{dr}$ divided by $R_{Larmor}$")
ax2.set_ylabel("pc")
ax2.set_xlabel("Magnetic field instance #")

ax2.plot(n, drR_larmor[0,:],'x', color='green')
ax2.plot([1,1000], [drR_larmorAverage[0],drR_larmorAverage[0]] , '-', color='black')

####### Figure 3 = Time with average
fig3 = plt.figure(3)
ax3 = fig3.gca()
ax3.set_title("Time difference in simulations")

ax3.set_ylim(0, 50)

ax3.plot(n, dt[0,:], 'x')
ax3.plot([1,1000], [dtAverage[0],dtAverage[0]], '-', color='black')

plt.show()

    

