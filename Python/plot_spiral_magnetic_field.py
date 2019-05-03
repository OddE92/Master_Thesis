import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

Bx = np.zeros((500, 500))
By = np.zeros((500, 500))
x = np.linspace(-500, 500, 500)
y = np.linspace(-500, 500, 500)


colors =  np.arctan2(y, x)                                      # Decides color based on the xy-plane angle

norm = Normalize()
norm.autoscale(colors)                              # Normalizes the colors to the given range of "colors"

theta = np.deg2rad(11.5)

si = np.sin(theta)
ci = np.cos(theta)

for x1 in range (0, 1000):
    for y1 in range (0, 1000):
        phi = np.arctan2(y1-500, x1-500)

        Bx[x1/2][y1/2] = si*np.cos(phi) - ci*np.sin(phi)
        By[x1/2][y1/2] = si*np.sin(phi) + ci*np.cos(phi)


fig = plt.figure()
ax = fig.gca()

Bx = np.multiply(Bx, 10)
By = np.multiply(By, 10)


colormap = cm.get_cmap('seismic')
ax.quiver(x, y, Bx, By, color=colormap(norm(colors)), units='width')

plt.show()