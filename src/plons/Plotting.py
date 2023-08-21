import matplotlib.pyplot as plt
import numpy as np

import plons.SmoothingKernelScript    as sk
import plons.PhysicalQuantities       as pq
import plons.ConversionFactors_cgs    as cgs

def plotSlice(ax, X, Y, smooth, observable, logplot=False, cmap = plt.cm.get_cmap('inferno'), clim=(None, None)):
    ax.set_aspect('equal')
    ax.set_facecolor('k')

    if logplot:
        obs = np.log10(smooth[observable]+1e-99)
    else:
        obs = smooth[observable]
    axPlot = ax.pcolormesh(X/cgs.au, Y/cgs.au, obs, cmap=cmap, vmin=clim[0], vmax = clim[1])

    cbar = plt.colorbar(axPlot, ax = ax, location='right', fraction=0.0471, pad=0.01)
    if logplot:
        cbar.set_label("log("+observable+")")
    else:
        cbar.set_label(observable)
    return cbar

def plotSink(ax, dumpData, setup, rotate=False):
    if rotate: circleAGB = plt.Circle((-np.linalg.norm(dumpData['posAGB'])/cgs.au, 0.), setup["wind_inject_radius"], transform=ax.transData._b, color="black", zorder=10)
    else: circleAGB = plt.Circle(dumpData['posAGB']/cgs.au, setup["wind_inject_radius"], color="black", zorder=10)
    ax.add_artist(circleAGB)
    
    if not setup['single_star']:
        if rotate: circleComp = plt.Circle((np.linalg.norm(dumpData['posComp'])/cgs.au, 0.), setup["rAccrComp"], transform=ax.transData._b, color="black", zorder=10)
        else: circleComp = plt.Circle(dumpData['posComp']/cgs.au, setup["rAccrComp"], color="black", zorder=10)
        ax.add_artist(circleComp)

        if setup['triple_star']:
            if rotate: circleComp_in = plt.Circle((np.linalg.norm(dumpData['posComp_in'])/cgs.au, 0.), setup["rAccrComp_in"], transform=ax.transData._b, color="black", zorder=10)
            else: circleComp_in = plt.Circle(dumpData['posComp_in']/cgs.au, setup["rAccrComp_in"], color="black", zorder=10)
            ax.add_artist(circleComp_in)

def SlicePlot2D(ax, dumpData, setup, n = 200, xlim=(-30, 30), ylim=(-30, 30), zlim=None, rotate = False,
                observable="rho", logplot=True, cmap = plt.cm.get_cmap('inferno'), clim = (-17, -14)):
    x = np.linspace(xlim[0], xlim[1], n)*cgs.au
    y = np.linspace(ylim[0], ylim[1], n)*cgs.au
    X, Y = np.meshgrid(x, y)
    Z    = np.zeros_like(X)
    
    if rotate:
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1]) # Calculate the angle around which to rotate
        X_rot, Y_rot, Z_rot = sk.rotateMeshAroundZ(theta, X, Y, Z)                        # Rotate the mesh grid before computing the smoothed data. In this way the image will be constructed in the rotated frame
        smooth = sk.smoothMesh(X_rot, Y_rot, Z_rot, dumpData, [observable])               # Smooth the data with the rotated mesh
    else:
        smooth = sk.smoothMesh(X, Y, Z, dumpData, [observable])

    cbar = plotSlice(ax, X, Y, smooth, observable, logplot=logplot, cmap=cmap, clim = clim)
    plotSink(ax, dumpData, setup, rotate)

    return ax, cbar