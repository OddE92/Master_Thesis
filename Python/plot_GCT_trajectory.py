import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import gaussian_kde
from numpy.fft import fft
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

trajectory = np.loadtxt('Data/trajectory_GC.dat')                               # Reads the data from the RK4 approx
#trajectory = np.loadtxt('Data/GC_trajectory_ODE.dat')
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')
ax.plot3D(trajectory[:,0], trajectory[:,1], trajectory[:,2])

ax.set_title('Trajectory BS step size control')
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')
ax.set_zlabel('z [pc]')

#ax.set_xlim(min(trajectory[:,0]), max(trajectory[:,0]))
#ax.set_ylim(min(trajectory[:,0]), max(trajectory[:,0]))
#ax.set_zlim(min(trajectory[:,0]), max(trajectory[:,0]))

#Plot multiple trajectories
# tr1 = np.loadtxt('data/RK4_approx_stepsizectrl_nonturb.dat')
# tr2 = np.loadtxt('data/RK4_approx_stepsizectrl_b+bturb_LM10.dat')
# tr3 = np.loadtxt('data/RK4_approx_stepsizectrl_b+bturb_LM150.dat')
# tr4 = np.loadtxt('data/RK4_approx_stepsizectrl_turb_LM10.dat')
# tr5 = np.loadtxt('data/RK4_approx_stepsizectrl_turb_LM150.dat')

# fig2 = plt.figure(3)
# fig2.suptitle('Sample of Trajectories')
# ax = fig2.add_subplot(2, 3, 1, projection='3d')
# first = ax.plot3D(tr1[:,0], tr1[:,1], tr1[:,2])

# ax = fig2.add_subplot(2, 3, 2, projection='3d')
# first = ax.plot3D(tr2[:,0], tr2[:,1], tr2[:,2])

# ax = fig2.add_subplot(2, 3, 3, projection='3d')
# first = ax.plot3D(tr3[:,0], tr3[:,1], tr3[:,2])

# ax = fig2.add_subplot(2, 3, 4, projection='3d')
# first = ax.plot3D(tr4[:,0], tr4[:,1], tr4[:,2])

# ax = fig2.add_subplot(2, 3, 5, projection='3d')
# first = ax.plot3D(tr5[:,0], tr5[:,1], tr5[:,2])


#ax[0,0].plot3D(tr1[:,0], tr1[:,1], tr1[:,2], color='red')
#ax[0,1].plot3D(tr2[:,0], tr2[:,1], tr2[:,2], color='red')
#ax[0,2].plot3D(tr3[:,0], tr3[:,1], tr3[:,2], color='red')
#ax[1,0].plot3D(tr4[:,0], tr4[:,1], tr4[:,2], color='red')
#ax[1,2].plot3D(tr5[:,0], tr5[:,1], tr5[:,2], color='red')


#Plot <B^2>
#fig2 = plt.figure(2)
#plt.plot(trajectory2)



#Plot mean field value
#f,(ax1, ax2) = plt.subplots(2)
#f.suptitle('Comparison of Turbulent Field Components \n for $L_\mathrm{max} = 10\mathrm{pc}$ and $L_\mathrm{max} = 150\mathrm{pc}$ ')

#x = np.arange(0, 500, 0.1)
#l = len(mfield)/2

#ax1.plot(x, mfield[0:l, 0], color='red')
#ax1.plot(x, mfield[0:l, 1], color='green')
#ax1.plot(x, mfield[0:l, 2], color='blue')

#ax1.set_ylabel('B [$\mu$G]')
#ax1.set_xlabel('x [pc]')
#ax1.legend(('$B_x$', '$B_y$', '$B_z$'), loc="upper right")
#ax1.set_title('$L_\mathrm{max} = 10$pc')

#ax2.plot(x, mfield[l:len(mfield), 0], color='red')
#ax2.plot(x, mfield[l:len(mfield), 1], color='green')
#ax2.plot(x, mfield[l:len(mfield), 2], color='blue')

#ax2.set_ylabel('B [$\mu$G]')
#ax2.set_xlabel('x [pc]')
#ax2.legend(('$B_x$', '$B_y$', '$B_z$'), loc="upper right")
#ax2.set_title('$L_\mathrm{max} = 150$pc')
#print(np.mean(mfield))

#Plot scatter of RNG-pulls
#fig4 = plt.figure(4)
#plt.scatter(test_rng[:,1], test_rng[:,0])



plt.show()