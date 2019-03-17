import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import figure
from pathlib import Path

directory_in_str='./'

E = np.loadtxt("exact.dat")

pathlist = Path(directory_in_str).glob('out*.dat')
for path in pathlist:

    M = np.loadtxt(str(path))
    ny, nx = M.shape
    x = np.linspace(-1, 1, nx)
    y = np.linspace(-1, 1, ny)
    X, Y = np.meshgrid(x, y)
    fig = plt.figure(figsize=plt.figaspect(0.5))
    # subplot 1 - calculated
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X,Y,M, rstride=1, cstride=1, cmap='viridis', linewidth=0, antialiased=False)
    ax.set_zlim(0,1.0)
    
    fig.colorbar(surf)
    surf.set_clim(0,1.0)
    ax.set_title("calculated solution") 
    # subplot 2 - exact solution
    
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf2 = ax2.plot_surface(X,Y,E, rstride=1, cstride=1, cmap='viridis', linewidth=0, antialiased=False)
    ax2.set_zlim(0,1.0)
    
    fig.colorbar(surf2)
    surf2.set_clim(0,1.0)
    ax2.set_title("exact solution")

    # image output and closing of plots
    out_str=str(path)
    out_str+= '.png'
    plt.savefig(out_str)
    plt.close("all")

