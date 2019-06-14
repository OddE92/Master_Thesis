import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import gaussian_kde
from numpy.fft import fft
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

#trajectory = np.loadtxt('Data/trajectory_GC.dat')                               # Reads the data from the RK4 approx
#trajectory = np.loadtxt('data/RK4_approx_stepsizectrl.dat')
#trajectory2 = np.loadtxt('data/B_average.dat')
#mfield = np.loadtxt('data/test_mf.dat')
#test_rng = np.loadtxt('data/test_rng.dat')

# fig1 = plt.figure(1)
# ax = fig1.gca(projection='3d')
# ax.plot3D(trajectory[:,0], trajectory[:,1], trajectory[:,2])

# ax.set_title('Trajectory BS step size control')
# ax.set_xlabel('x [pc]')
# ax.set_ylabel('y [pc]')
# ax.set_zlabel('z [pc]')

#Plot multiple trajectories
tr1 = np.loadtxt('Data/trajectory_GC.dat')
tr2 = np.loadtxt('Data/trajectory.dat')

fig2 = plt.figure(3)
fig2.suptitle('Sample of Trajectories')
ax = fig2.add_subplot(1, 2, 1, projection='3d')
first = ax.plot3D(tr1[:,0], tr1[:,1], tr1[:,2])

ax.set_xlim(min(tr2[:,0]), max(tr2[:,0]))
ax.set_ylim(min(tr2[:,0]), max(tr2[:,0]))
ax.set_zlim(min(tr2[:,0]), max(tr2[:,0]))

ax = fig2.add_subplot(1, 2, 2, projection='3d')
first = ax.plot3D(tr2[:,0], tr2[:,1], tr2[:,2])

ax.set_xlim(min(tr2[:,0]), max(tr2[:,0]))
ax.set_ylim(min(tr2[:,0]), max(tr2[:,0]))
ax.set_zlim(min(tr2[:,0]), max(tr2[:,0]))

# ax = fig2.add_subplot(2, 3, 3, projection='3d')
# first = ax.plot3D(tr3[:,0], tr3[:,1], tr3[:,2])

# ax = fig2.add_subplot(2, 3, 4, projection='3d')
# first = ax.plot3D(tr4[:,0], tr4[:,1], tr4[:,2])

# ax = fig2.add_subplot(2, 3, 5, projection='3d')
# first = ax.plot3D(tr5[:,0], tr5[:,1], tr5[:,2])


#ax[0,0].plot3D(tr1[:,0], tr1[:,1], tr1[:,2], color='red')
#ax[0,1].plot3D(tr2[:,0], tr2[:,1], tr2[:,2], color='black')
#ax[0,2].plot3D(tr3[:,0], tr3[:,1], tr3[:,2], color='red')
#ax[1,0].plot3D(tr4[:,0], tr4[:,1], tr4[:,2], color='red')
#ax[1,2].plot3D(tr5[:,0], tr5[:,1], tr5[:,2], color='red')



plt.show()