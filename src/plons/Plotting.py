from typing import Dict, Tuple, Any, Optional

import numpy        as np
import numpy.typing as npt
import matplotlib
import matplotlib.pyplot   as plt

import plons.SmoothingKernelScript    as sk
import plons.PhysicalQuantities       as pq
import plons.ConversionFactors_cgs    as cgs

def plotSlice(ax: plt.Axes,
              X: npt.NDArray[np.single],
              Y: npt.NDArray[np.single],
              smooth: Dict[str, npt.NDArray[np.single]],
              observable: str,
              logplot: bool = False,
              cmap: matplotlib.colors.Colormap = plt.colormaps['inferno'],
              clim: Tuple[Optional[float], Optional[float]] = (None, None)) -> matplotlib.colorbar.Colorbar:
    """Plot a property given a grid and smoothed data ontop of the grid

    Args:
        ax (plt.Axes): axis of figure on which you want to plot the slice
        X (npt.NDArray[np.single]): X values in meshgrid which you want to plot
        Y (npt.NDArray[np.single]): Y values in meshgrid which you want to plot
        smooth (Dict[str, npt.NDArray[np.single]]): Dictionary pointing at smoothed values in meshgrid which you want to plot
        observable (str): Name of the observable you want to plot, corresponding to the name in the smooth directory
        logplot (bool, optional): plot in log scale?. Defaults to False.
        cmap (matplotlib.colors.Colormap, optional): Colormap to use. Defaults to plt.colormaps['inferno'].
        clim (Tuple[Optional[float], Optional[float]], optional): limits for the colorbar. Defaults to (None, None).

    Returns:
        colorbar.Colorbar: Colorbar
    """

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

def plotSink(ax: plt.Axes,
             dumpData: Dict[str, Any],
             setup: Dict[str, Any],
             rotate: bool = False) -> matplotlib.patches.Circle:
    """Plot the sink particles on the axis

    Args:
        ax (plt.Axes): axis of figure on which you want to plot the slice
        dumpData (Dict[str, Any]): Data of the dump file
        setup (Dict[str, Any]): Setup of the simulation
        rotate (bool, optional): Should the sinks be rotated. Defaults to False.

    Returns:
        matplotlib.patches.Circle: Circle objects of the sink particles
    """
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
            return circleAGB, circleComp, circleComp_in
        return circleAGB, circleComp
    return circleAGB

def SlicePlot2D(ax: plt.Axes,
                dumpData: Dict[str, Any],
                setup: Dict[str, Any],
                n: int = 200,
                xlim: tuple[float, float] = (-30, 30),
                ylim: tuple[float, float] = (-30, 30),
                rotate: bool = False,
                observable: str = "rho",
                logplot: bool = True,
                cmap: matplotlib.colors.Colormap = plt.colormaps['inferno'],
                clim: tuple[float, float] = (-17, -14)) -> matplotlib.colorbar.Colorbar:
    """Plot a property given xlims and ylims

    Args:
        ax (plt.Axes):  axis of figure on which you want to plot the slice
        dumpData (Dict[str, Any]):  Data of the dump file
        setup (Dict[str, Any]):  Setup of the simulation
        n (int, optional): amount of points on each axis to plot. Defaults to 200.
        xlim (tuple[float, float], optional): xlimits for the plot. Defaults to (-30, 30).
        ylim (tuple[float, float], optional): ylimits for the plot. Defaults to (-30, 30).
        rotate (bool, optional): should the binary be rotated to lay on the x-axis. Defaults to False.
        observable (str, optional): property to plot. Defaults to "rho".
        logplot (bool, optional): plot in log scale?. Defaults to True.
        cmap (matplotlib.colors.Colormap, optional): colormap to use. Defaults to plt.colormaps['inferno'].
        clim (tuple[float, float], optional): limits of the colormap. Defaults to (-17, -14).

    Returns:
        matplotlib.colorbar.Colorbar: the colorbar in the plot
    """
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

    return cbar