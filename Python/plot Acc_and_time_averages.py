import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import numpy as np
import itertools as it

basename = "Data/Compare/Ee"
contname = "_Bpercent"

drR_larmorAverage = []
drAverage = []
R_larmorAverage = []
dtAverage = []

drR_larmorAverage_E = np.empty((0,5))
drAverage_E = np.empty((0,5))
R_larmorAverage_E = np.empty((0,5))
dtAverage_E = np.empty((0,5))


Earr = [1, 10, 100, 1000]

for E in range (15, 19):
    for percent in ["1", "10", "25", "50", "100"]:

        strGCT      = basename + str(E) + contname + percent + "_GCT.dat"
        strExact    = basename + str(E) + contname + percent + ".dat"

        dataGCT     = np.loadtxt(strGCT)
        dataExact   = np.loadtxt(strExact)

        temp1   = np.sqrt( dataGCT[:,0]**2 + dataGCT[:,1]**2 + dataGCT[:,2]**2 )
        temp2   = np.sqrt( dataExact[:,0]**2 + dataExact[:,1]**2 + dataExact[:,2]**2 )
        tempdr  = np.sqrt( (dataGCT[:,0]-dataExact[:,0])**2 + (dataGCT[:,1]-dataExact[:,1])**2 + (dataGCT[:,2]-dataExact[:,2])**2 )
        tempdt  = dataExact[:,3] / dataGCT[:,4]

        drR_larmorAverage = np.append(drR_larmorAverage, np.mean(tempdr/dataGCT[:,3]))

        drAverage = np.append(drAverage, np.mean(tempdr))

        R_larmorAverage = np.append(R_larmorAverage, np.mean(dataGCT[:,3]))
        
        dtAverage = np.append(dtAverage, np.mean(tempdt))
    
    drR_larmorAverage_E = np.append(drR_larmorAverage_E, [drR_larmorAverage], axis=0)
    drAverage_E = np.append(drAverage_E, [drAverage], axis=0)
    dtAverage_E = np.append(dtAverage_E, [dtAverage], axis=0)
    R_larmorAverage_E = np.append(R_larmorAverage_E, [R_larmorAverage], axis=0)

    drR_larmorAverage = []
    drAverage = []
    R_larmorAverage = []
    dtAverage = []

################################################################
# 
# Figure 1 plots dr/R_L averages
#
################################################################
matplotlib.rc('xtick', labelsize=24)
matplotlib.rc('ytick', labelsize=24)


fig1 = plt.figure(1)
ax = fig1.gca()

marker = it.cycle((r'$\dag$', '1', '+', '2', 'x'))
color = it.cycle(("#FFA500", 'g', 'r', 'k', 'b'))

for i in range(4, -1, -1):
    ax.plot(Earr, drR_larmorAverage_E[:,i], linestyle='', marker=marker.next(), markersize=22, color=color.next())

ax.legend(  [r"$B_{\mathrm{turb}}/B_0 = 1.00$", r"$B_{\mathrm{turb}}/B_0 = 0.50$", r"$B_{\mathrm{turb}}/B_0 = 0.25$",
             r"$B_{\mathrm{turb}}/B_0 = 0.10$", r"$B_{\mathrm{turb}}/B_0 = 0.01$"], 
             fontsize='28', frameon=False)


#ax.set_title(r"Average length of $dr$ proportional to average $R_L$ as functions of energy and $B_0/B_{\mathrm{turb}}$", 
#                fontsize='22', y=1.03)
ax.set_xlabel(r"Energy [PeV]", fontsize='28')
ax.set_ylabel(r"dr/$\mathrm{R_L}$", fontsize='28')

################################################################
# 
# Figure 2 plots timeExact/timeGCT-averages
#
################################################################
fig2 = plt.figure(2)
ax2 = fig2.gca()

for i in range(4, -1, -1):
    ax2.plot(Earr, dtAverage_E[:,i], linestyle='', marker=marker.next(), markersize=22, color=color.next())

ax2.legend( [r"$B_{\mathrm{turb}}/B_0 = 1.00$", r"$B_{\mathrm{turb}}/B_0 = 0.50$", r"$B_{\mathrm{turb}}/B_0 = 0.25$",
             r"$B_{\mathrm{turb}}/B_0 = 0.10$", r"$B_{\mathrm{turb}}/B_0 = 0.01$"], 
             fontsize='28', frameon=False)

#ax2.set_title(r"Average runtime of the exact solution proportional to the average runtime of the GC-solution", 
#                fontsize='18', y=1.03)
ax2.set_xlabel(r"Energy [PeV]", fontsize='28')
ax2.set_ylabel(r"t/$\mathrm{\tau}$", fontsize='28')


ax.set_yscale("log")
ax.set_xscale("log")

ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

ax2.set_xscale("log")
ax2.xaxis.set_major_formatter(ticker.ScalarFormatter())

fig1.set_size_inches(30, 22.5)
fig1.savefig("Figures/Final/Average_dr_Rlarmor.png", bbox_inches='tight')

fig2.set_size_inches(30, 22.5)
fig2.savefig("Figures/Final/Average_runtime.png", bbox_inches='tight')

plt.show()
