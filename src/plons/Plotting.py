from typing import Dict, Tuple, Any, Optional

import numpy        as np
import numpy.typing as npt
import matplotlib
import matplotlib.pyplot   as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

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
              clim: Tuple[Optional[float], Optional[float]] = (None, None),
              cbar = True) -> matplotlib.colorbar.Colorbar:
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
        cbar (bool, optional): Should a colorbar be plotted? Defaults to True.

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
    ax.set_xlim(X.min()/cgs.au, X.max()/cgs.au)
    ax.set_ylim(Y.min()/cgs.au, Y.max()/cgs.au)

    if cbar == True:
        cbar = plt.colorbar(axPlot, ax = ax, location='right', fraction=0.0471, pad=0.01)
        if logplot:
            cbar.set_label("log("+observable+")")
        else:
            cbar.set_label(observable)
        return cbar
    else:
        return axPlot

def plotSink(ax: plt.Axes,
             dumpData: Dict[str, Any],
             setup: Dict[str, Any],
             rotate: bool = False,
             plane: str = "xy") -> matplotlib.patches.Circle:
    """Plot the sink particles on the axis

    Args:
        ax (plt.Axes): axis of figure on which you want to plot the slice
        dumpData (Dict[str, Any]): Data of the dump file
        setup (Dict[str, Any]): Setup of the simulation
        rotate (bool, optional): Should the sinks be rotated. Defaults to False.
        plane (str, optional): Plane to plot the sinks on. Defaults to "xy".

    Returns:
        matplotlib.patches.Circle: Circle objects of the sink particles
    """
    if rotate: plotAGB = (-np.linalg.norm(dumpData._params['posAGB'])/cgs.au, 0.)
    else: plotAGB = dumpData._params['posAGB'][:2]/cgs.au if plane == "xy" else dumpData._params['posAGB'][[0,2]]/cgs.au
    circleAGB = plt.Circle(plotAGB, setup["wind_inject_radius"], color="black", zorder=10)
    ax.add_patch(circleAGB)
    
    if not setup['single_star']:
        if rotate: plotComp = (np.linalg.norm(dumpData._params['posComp'])/cgs.au, 0.)
        else: plotComp = dumpData._params['posComp'][:2]/cgs.au if plane == "xy" else dumpData._params['posComp'][[0,2]]/cgs.au
        circleComp = plt.Circle(plotComp, setup["rAccrComp"], color="black", zorder=10)
        ax.add_patch(circleComp)

        if setup['triple_star']:
            if rotate: plotComp_in = (np.linalg.norm(dumpData._params['posComp_in'])/cgs.au, 0.)
            else: plotComp_in = dumpData._params['posComp_in'][:2]/cgs.au if plane == "xy" else dumpData._params['posComp_in'][[0,2]]/cgs.au
            circleComp_in = plt.Circle(plotComp_in, setup["rAccrComp_in"], color="black", zorder=10)
            ax.add_patch(circleComp_in)
            return circleAGB, circleComp, circleComp_in
        return circleAGB, circleComp
    return circleAGB

def SlicePlot2D(ax: plt.Axes,
                dumpData: Dict[str, Any],
                setup: Dict[str, Any],
                n: int = 512,
                xlim: tuple[float, float] = (-30, 30),
                ylim: tuple[float, float] = (-30, 30),
                rotate: bool = False,
                plane: str = "xy",
                observable: str = "rho",
                logplot: bool = True,
                quiver: bool = False,
                n_arrows: int = 25,
                cmap: matplotlib.colors.Colormap = plt.colormaps['inferno'],
                clim: tuple[float, float] = (-17, -14),
                cbar = True) -> matplotlib.colorbar.Colorbar:
    """Plot a property given xlims and ylims

    Args:
        ax (plt.Axes):  axis of figure on which you want to plot the slice
        dumpData (Dict[str, Any]):  Data of the dump file
        setup (Dict[str, Any]):  Setup of the simulation
        n (int, optional): amount of points on each axis to plot. Defaults to 200.
        xlim (tuple[float, float], optional): xlimits for the plot. Defaults to (-30, 30).
        ylim (tuple[float, float], optional): ylimits for the plot. Defaults to (-30, 30).
        rotate (bool, optional): should the binary be rotated to lay on the x-axis. Defaults to False.
        plane (str, optional): plane to plot the slice in. Defaults to "xy".
        observable (str, optional): property to plot. Defaults to "rho".
        logplot (bool, optional): plot in log scale?. Defaults to True.
        quiver (bool, optional): Should a quiver plot be added? Defaults to False.
        cmap (matplotlib.colors.Colormap, optional): colormap to use. Defaults to plt.colormaps['inferno'].
        clim (tuple[float, float], optional): limits of the colormap. Defaults to (-17, -14).
        cbar (bool, optional): Should a colorbar be plotted? Defaults to True. If 'top', the colorbar is plotted on the top of the axis, if False, no colorbar is plotted.

    Returns:
        matplotlib.colorbar.Colorbar: the colorbar in the plot
    """
    
    if rotate:
        theta = pq.getPolarAngleCompanion(dumpData._params['posComp'][0], dumpData._params['posComp'][1]) # Calculate the angle around which to rotate

    if logplot:
        clim = (10**clim[0], 10**clim[1])

    # Change only the x and y data from cm to au. Nessessary for the rendering
    dumpData.x = dumpData.x / cgs.au
    dumpData.y = dumpData.y / cgs.au
    dumpData.z = dumpData.z / cgs.au
    dumpData.h = dumpData.h / cgs.au

    divider = make_axes_locatable(ax)
    if cbar == 'top':
        cax = divider.append_axes("top", size="5%", pad=0.05)
    elif cbar == True:
        cax = divider.append_axes("right", size="5%", pad=0.05)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax = dumpData.render(
            observable,
            ax=ax,
            rotation=(-180/np.pi*theta if rotate else 0, 0, 90 if plane == "xz" else 0),
            rot_origin=[0, 0, 0],
            xlim=xlim,
            x_pixels=n,
            y_pixels=n,
            ylim=ylim,
            log_scale=logplot,
            cmap=cmap,
            cbar_ax=cax if cbar else None,
            cbar_kws={
                "orientation": "horizontal" if cbar == 'top' else 'vertical',
                'label': f"log({observable})" if logplot else observable
                },
            vmin=clim[0], vmax=clim[1],
            xsec=0,
            normalize=True
            )
        
        if quiver:
            ax = dumpData.arrowplot(
                ("vx", "vy", "vz"),
                ax=ax,
                rotation=(-180/np.pi*theta, 0, 0) if rotate else (0, 0, 0),
                rot_origin=[0, 0, 0],
                xlim=xlim,
                ylim=ylim,
                x_arrows=n_arrows, y_arrows=n_arrows,
                qkey=False,
                normalize=True
            )
    
    # revert the changes to the dumpData so that other plots are not affected
    dumpData.x = dumpData.x * cgs.au
    dumpData.y = dumpData.y * cgs.au
    dumpData.z = dumpData.z * cgs.au
    dumpData.h = dumpData.h * cgs.au

    colorbar = ax.get_images()[0].colorbar
    if cbar == 'top':
        colorbar.ax.xaxis.set_label_position('top')
        colorbar.ax.xaxis.set_ticks_position('top')
    if cbar == False:
        colorbar.remove()

    plotSink(ax, dumpData, setup, rotate, plane)

    return colorbar