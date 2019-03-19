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

if nx > 100:
    XS=round(nx/100)
else:
    XS=1

if ny > 100:
    YS=round(ny/100)
else:
    YS=1

pathlist = Path(directory_in_str).glob('out_*.dat')
for path in pathlist:

    M = np.loadtxt(str(path))
    fig = plt.figure(figsize=plt.figaspect(0.5))
    # subplot 1 - calculated
    ax = fig.add_subplot(1, 2, 1)
    #surf = ax.plot_surface(X,Y,M, rstride=YS, cstride=XS, cmap='viridis', linewidth=0, antialiased=False)
    surf = ax.pcolor(M)
    
    fig.colorbar(surf)
    surf.set_clim(-1.0,1.0)
    ax.set_title("calculated solution") 
    # subplot 2 - exact solution
    
    ax2 = fig.add_subplot(1, 2, 2)
    #surf2 = ax2.plot_surface(X,Y,E, rstride=YS, cstride=XS, cmap='viridis', linewidth=0, antialiased=False)
    surf2 = ax2.pcolor(E)
    
    fig.colorbar(surf2)
    surf2.set_clim(-1.0,1.0)
    ax2.set_title("exact solution")

    # image output and closing of plots
    out_str=str(path)
    out_str+= '.png'
    plt.savefig(out_str)
    plt.close("all")

