import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
#from pylab import figure
from pathlib import Path
import sys

if len(sys.argv)>1:
    directory_in_str=sys.argv[1]
else:
    directory_in_str="./"

E_file=directory_in_str + 'exact.dat'
E = np.loadtxt(E_file)

ny, nx = E.shape
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

f0=directory_in_str + 'out_000000.dat'
M = np.loadtxt(f0)

fig = plt.figure(figsize=plt.figaspect(0.33))

# subplot 1 - calculated
ax = fig.add_subplot(1, 3, 1)
surf = ax.pcolormesh(M)

fig.colorbar(surf)
surf.set_clim(-1.0,1.0)
ax.set_title("calculated solution") 

# subplot 2 - exact solution

ax2 = fig.add_subplot(1, 3, 2)
surf2 = ax2.pcolormesh(E)

fig.colorbar(surf2)
surf2.set_clim(-1.0,1.0)
ax2.set_title("exact solution")

# subplot 3 - error

ax3 = fig.add_subplot(1, 3, 3)
surf3 = ax3.pcolormesh(E)

fig.colorbar(surf3)
surf3.set_clim(-1e-8,1e-8)
ax3.set_title("exact solution")
pathlist = Path(directory_in_str).glob('out_*.dat')
for path in pathlist:
    
    M = np.loadtxt(str(path))
    surf.set_array(M.ravel())

    Error = M-E
    surf3.set_array(Error.ravel())

    # redraw to update plot with new data
    plt.draw() 

    # image output and closing of plots
    out_str=str(path)
    out_str+= '.png'
    plt.savefig(out_str)

