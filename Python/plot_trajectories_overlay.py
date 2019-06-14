import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import gaussian_kde
from numpy.fft import fft
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

matplotlib.rc('xtick', labelsize=24)
matplotlib.rc('ytick', labelsize=24)

trajectory = np.loadtxt('Data/trajectory.dat')
GC = np.loadtxt('Data/trajectory_GC.dat')

#fig1 = plt.figure(1)
#ax = fig1.gca(projection='3d')

#ax.plot3D(trajectory[:,0], trajectory[:,1], trajectory[:,2], color='blue')
#ax.plot3D(GC[:,0], GC[:,1], GC[:,2], color='red')

#ax.set_title(r'Comparison of GC (Red) and Trajectory, T = $5\cdot10^3$ years, E = $10^{17}$ eV')
#ax.set_xlabel('x [pc]', fontsize=28)
#ax.set_ylabel('y [pc]', fontsize=28)
#ax.set_zlabel('z [pc]', fontsize=28)

#ax.set_xlim(min(trajectory[:,0]), max(trajectory[:,0]))
#ax.set_ylim(min(trajectory[:,0]), max(trajectory[:,0]))
#ax.set_zlim(min(trajectory[:,0]), max(trajectory[:,0]))

fig2, ax2 = plt.subplots(1, 2)

ax2[0].plot(trajectory[:,0], trajectory[:,1], color='blue')
ax2[0].plot(GC[:,0], GC[:,1], color='red')
ax2[0].set_xlabel('x [pc]', fontsize=28)
ax2[0].set_ylabel('y [pc]', fontsize=28)

ax2[1].plot(trajectory[:,0], trajectory[:,2], color='blue')
ax2[1].plot(GC[:,0], GC[:,2], color='red')
ax2[1].set_xlabel('x [pc]', fontsize=28)
ax2[1].set_ylabel('z [pc]', fontsize=28)

fig2.subplots_adjust(wspace=0.3)


plt.show()