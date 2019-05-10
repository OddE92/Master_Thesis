import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

bfield = np.loadtxt('Data/spiral_MF.dat')
bfield2 = np.loadtxt('Data/spiral_MF2.dat')

colors =  np.arctan2(bfield[:,2], bfield[:,3])      # Decides color based on the xy-plane angle

norm = Normalize()
norm.autoscale(colors)                              # Normalizes the colors to the given range of "colors"

fig = plt.figure()
ax = fig.gca()

colormap = cm.get_cmap('seismic')
ax.quiver(bfield[:,0], bfield[:,1], bfield[:,2], bfield[:,3], scale=50, color='red')
#ax.quiver(bfield2[:,0], bfield2[:,1], bfield2[:,2], bfield2[:,3], scale=50, color='blue')

plt.show()