from __future__ import print_function, absolute_import, division 
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma
from pyKratos import *


def PlotContour(Nodes, variable,name):
    nnodes = len(Nodes)
    x = []
    y = []
    z = []
    for node in Nodes:
        x.append(node.coordinates[0])
        y.append(node.coordinates[1])
        z.append(node.GetSolutionStepValue(variable, 0))

    x = np.array(x)
    y = np.array(y)
    xmin = np.amin(x)
    xmax = np.amax(x)
    ymin = np.amin(y)
    ymax = np.amax(y)
    print("xmin = ",xmin," xmax = ",xmax)
    print("ymin = ",ymin," ymax = ",ymax)
    #print(x,y,z)
    
    # define grid.
    xi = np.linspace(xmin, xmax, 100)
    yi = np.linspace(ymin, ymax, 100)
    # grid the data.
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear') #nearest linear cubic
        
    # contour the gridded data, plotting dots at the randomly spaced data
    # points.
    CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.jet)
    plt.colorbar()  # draw colorbar
    # plot data points.
    plt.scatter(x, y, marker='o', c='b', s=5)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title('griddata test (%d points)' % nnodes)
    plt.axes().set_aspect('equal', 'datalim')

    plt.savefig(name)
    plt.close()
    #plt.ion()
    #plt.show()
    #plt.close()
