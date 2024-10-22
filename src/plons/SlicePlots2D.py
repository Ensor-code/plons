#Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import os
import math

# import plons scripts
import plons.SmoothingKernelScript    as sk
import plons.ConversionFactors_cgs    as cgs
import plons.Plotting                 as plot

# import certain things from packages
from mpl_toolkits.axes_grid1    import make_axes_locatable
from matplotlib.ticker          import MultipleLocator

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

roundToInteger  = True
velocity_vec    = True


# colormap per parameter
cm =   {'rho':   plt.colormaps['inferno'],
        'speed': plt.colormaps['CMRmap'],
        #'Tgas': plt.colormaps['RdYlGn'],
        'Tgas':  plt.colormaps['hot'], #nipy_spectral
        'tau':   plt.colormaps['viridis_r'],
        'tauL':  plt.colormaps['viridis'],
        'kappa': plt.colormaps['Spectral_r'],
        'Tdust': plt.colormaps['Spectral_r'],
        'Gamma': plt.colormaps['Spectral_r']
        }

# label name per parameter
name = {'rho':   r'$\log \, \rho$ [g$\,$cm$^{-3}$]',
        'speed': r'$v$ [km/s]',
        'Tgas':  r'$\log \, T$ [K]',
        'tau':   r'$\tau$ [/]',
        'tauL':  r'$\tau_L$ [/]',
        'kappa': r'$\kappa$ [g/cm$^3$]',
        'Tdust': r'$T_{\rm eq}$ [K]',
        'Gamma': r'$\Gamma$ [/]'
        }

# label name per parameter
log =  {'rho':   True,
        'speed': False,
        'Tgas':  True,
        'tau':   False,
        'tauL':  False,
        'kappa': False,
        'Tdust': False,
        'Gamma': False
        }

'''
Load the smoothing kernel data
'''
def smoothData(dumpData, setup, observables, theta = 0., zoom = 1, nneighb = 10, n_grid = 200, n_grid_vec = 25):
    pixCoord = sk.getPixels('z', n_grid, 'comp', dumpData, (setup['bound']) * cgs.au * np.sqrt(2.) / 2. / zoom)
    results_sph_sl_z = sk.getSmoothingKernelledPix(nneighb, dumpData, observables, sk.rotatePixCoordAroundZ(theta, pixCoord))
    X1, Y1, Z1, results_sph_sl_z = sk.convertToMesh(pixCoord, results_sph_sl_z, observables)

    pixCoord = sk.getPixels('z', n_grid_vec, 'comp', dumpData, (setup['bound']) * cgs.au * np.sqrt(2.) / 2. / zoom)
    results_sph_sl_z_vec = sk.getSmoothingKernelledPix(nneighb, dumpData, ['vx', 'vy', 'vz'], sk.rotatePixCoordAroundZ(theta, pixCoord))
    VX1, VY1, VZ1, results_sph_sl_z_vec = sk.convertToMesh(pixCoord, results_sph_sl_z_vec, ['vx', 'vy', 'vz'])
    results_sph_sl_z_vec = sk.rotateVelocityAroundZ(-theta, results_sph_sl_z_vec)
    
    pixCoord = sk.getPixels('y', n_grid, 'comp', dumpData, (setup['bound']) * cgs.au * np.sqrt(2.) / 2. / zoom)
    results_sph_sl_y = sk.getSmoothingKernelledPix(nneighb, dumpData, observables, sk.rotatePixCoordAroundZ(theta, pixCoord))
    X2, Y2, Z2, results_sph_sl_y = sk.convertToMesh(pixCoord, results_sph_sl_y, observables)
    
    pixCoord = sk.getPixels('y', n_grid_vec, 'comp', dumpData, (setup['bound']) * cgs.au * np.sqrt(2.) / 2. / zoom)
    results_sph_sl_y_vec = sk.getSmoothingKernelledPix(nneighb, dumpData, ['vx', 'vy', 'vz'], sk.rotatePixCoordAroundZ(theta, pixCoord)) 
    VX2, VY2, VZ2, results_sph_sl_y_vec = sk.convertToMesh(pixCoord, results_sph_sl_y_vec, ['vx', 'vy', 'vz'])
    results_sph_sl_y_vec = sk.rotateVelocityAroundZ(-theta, results_sph_sl_y_vec)

    smooth = {  'smooth_z'     :  results_sph_sl_z,
                'x_z'          :  X1,
                'y_z'          :  Y1,
                'smooth_y'     :  results_sph_sl_y,
                'x_y'          :  X2,
                'z_y'          :  Z2
                }

    smooth_vec = {  'smooth_z' :  results_sph_sl_z_vec,
                    'x_z'      :  VX1,
                    'y_z'      :  VY1,
                    'smooth_y' :  results_sph_sl_y_vec,
                    'x_y'      :  VX2,
                    'z_y'      :  VZ2
                    }

    return smooth, smooth_vec


'''
Make figure with the xy(orbital plane) slice plot of log density [g/cm^3]. 
    - smooth            is smoothing kernel data
    - zoom              zoom factor of the plot
    - rhoMin/rhoMax     determine the min and max of the colorbar
    - dumpData      data of dump file
    - setup         setup data
    - rAccComp      accretion radius of the companion
'''
def densityPlot(smooth, zoom, limits, dumpData, setup, orbital=True, cmap=plt.colormaps['inferno']):
    
    fig, ax = plt.subplots(1, figsize=(7, 7))

    ax.set_aspect('equal')
    ax.set_facecolor('k')

    axPlot = None
    if orbital:
        dataRho = np.log10(smooth[zoom]['smooth_z']["rho"]+1e-99)
        axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.au, smooth[zoom]['y_z'] / cgs.au,
                            dataRho, cmap=cmap, vmin=limits[0], vmax=limits[1],
                            rasterized=True)
    else:
        dataRho = np.log10(smooth[zoom]['smooth_y']["rho"]+1e-99)
        axPlot = ax.pcolormesh(smooth[zoom]['x_y'] / cgs.au, smooth[zoom]['z_y'] / cgs.au,
                            dataRho, cmap=cmap, vmin=limits[0], vmax=limits[1],
                            rasterized=True)

    plot.plotSink(ax, dumpData, setup, rotate=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="7%", pad=0.15)
    cbar = plt.colorbar(axPlot, cax=cax)
    cbar.set_label(r'$\log \, \rho$ [g$\,$cm$^{-3}$]', fontsize=24)
    cbar.ax.tick_params(labelsize=20)

    lim = (setup['bound'] * np.sqrt(2.)/2.) / zoom
    if roundToInteger:
        if lim >= 2.:
            lim = math.floor(lim / 2.) * 2.
        elif lim >= 1.:
            lim = 1.

    ax.axis([-lim, lim, -lim + 0.001, lim])
    ax.xaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.xaxis.set_minor_locator(MultipleLocator(lim / 8.))
    ax.yaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.yaxis.set_minor_locator(MultipleLocator(lim / 8.))
    ax.set_xlabel(r"$x$ [AU]", fontsize=22)
    if orbital: ax.set_ylabel(r"$y$ [AU]", fontsize=22)
    else:       ax.set_ylabel(r"$z$ [AU]", fontsize=22)

    ax.tick_params(labelsize=20)

    return fig


'''
Makes one slice plot 

INPUT
    - ax        is given subplot
    - par       is the name of the parameter
    - mi/ma     colorbar limits
    - smooth    smoothed data
    - zoom      zoom factor for the plot
'''

def onePlot(fig, ax, par, limits, smooth, smooth_vec, zoom, dumpData, setup, plane, cbar=True, velocity_vec = False):
    axPlot = None
    if plane == 'z':
        if log[par]:
            data = np.log10(smooth[zoom]['smooth_z'][par]+1e-99)
        else:
            data = smooth[zoom]['smooth_z'][par]
            
        axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.au, smooth[zoom]['y_z'] / cgs.au,
                                data, cmap=cm[par], vmin=limits[0], vmax=limits[1],
                                rasterized=True, shading="nearest")

        if par == 'speed' and velocity_vec:
            vx = smooth_vec[zoom]['smooth_z']['vx']
            vy = smooth_vec[zoom]['smooth_z']['vy']
            normaliseVectorLength = np.hypot(vx, vy)
            ax.quiver(smooth_vec[zoom]['x_z'] / cgs.au, smooth_vec[zoom]['y_z'] / cgs.au,
                        vx / normaliseVectorLength, vy / normaliseVectorLength, scale_units="dots", scale=0.05)

    if plane == 'y':
        if log[par]:
            data = np.log10(smooth[zoom]['smooth_y'][par]+1e-99)
        else:
            data = smooth[zoom]['smooth_y'][par]

        axPlot = ax.pcolormesh(smooth[zoom]['x_y'] / cgs.au, smooth[zoom]['z_y'] / cgs.au, data,
                                cmap=cm[par], vmin=limits[0], vmax=limits[1], rasterized=True, shading="nearest")

        if par == 'speed' and velocity_vec:
            vx = smooth_vec[zoom]['smooth_y']['vx']
            vz = smooth_vec[zoom]['smooth_y']['vz']
            normaliseVectorLength = np.hypot(vx, vz)
            ax.quiver(smooth_vec[zoom]['x_y'] / cgs.au, smooth_vec[zoom]['z_y'] / cgs.au,
                        vx / normaliseVectorLength, vz / normaliseVectorLength, scale_units="dots", scale=0.05)

    plot.plotSink(ax, dumpData, setup, rotate=True)

    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="7%", pad=0.15)
        cbar = fig.colorbar(axPlot, cax=cax)
        cbar.set_label(name[par], fontsize=24)
        cbar.ax.tick_params(labelsize=20)

    ax.tick_params(labelsize=20)
    ax.set_aspect('equal')
    ax.set_facecolor('k')

    lim = (setup['bound'] * np.sqrt(2.) / 2.) / zoom
    if roundToInteger:
        if lim >= 2.:
            lim = math.floor(lim / 2.) * 2.
        elif lim >= 1.:
            lim = 1.
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim + 0.001, lim)

    ax.xaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.xaxis.set_minor_locator(MultipleLocator(lim / 8.))
    ax.yaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.yaxis.set_minor_locator(MultipleLocator(lim / 8.))


'''
Make plot with 6 pannels:
    left:    xy-plane (= orbital plane = face-on plane)
    right:   xz-plane (= edge-on plane)
    upper:   log density [g/cm^3]
    middle:  speed [km/s]
    lower:   log temperature [K].

INPUT:
    - smooth        smoothed data
    - zoom          zoom factor for the plot
    - rhoMin/rhoMax limits of colorbar of density plots
    - vmax          upper limit of colorbar velocity plots
    - dumpData      data of dump file
    - setup         setup data
    - run           run number [str]
    - loc           output directory
'''


def allPlots(smooth, smooth_vec, zoom, limits, dumpData, setup, observables):

    fig = plt.figure(figsize=(10 + 0.35*5, 5*len(observables)))
    spec = fig.add_gridspec(ncols=2, nrows=len(observables), width_ratios=[1., 1.], height_ratios=[1.]*len(observables))

    axs = []
    for i in range(len(observables)):
        axs.append(fig.add_subplot(spec[i, 0]))
        axs.append(fig.add_subplot(spec[i, 1]))

    for i in range(len(observables)):
        if observables[i] == 'speed':
            onePlot(fig, axs[2*i],   observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, 'z', cbar=False, velocity_vec = velocity_vec)
            onePlot(fig, axs[2*i+1], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, 'y', cbar=True, velocity_vec = velocity_vec)
        else:
            onePlot(fig, axs[2*i],   observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, 'z', cbar=False)
            onePlot(fig, axs[2*i+1], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, 'y', cbar=True)

    for i in range(len(observables)-1):
        axs[2*i].set_xticklabels([])
        axs[2*i+1].set_xticklabels([])
    for i in range(len(observables)):
        axs[2*i].set_ylabel(r"$y$ [AU]", fontsize=22)
        axs[2*i+1].set_ylabel(r"$z$ [AU]", fontsize=22)
        axs[2*i+1].set_yticklabels([])
    axs[-1].set_xlabel(r"$x$ [AU]", fontsize=22)
    axs[-2].set_xlabel(r"$x$ [AU]", fontsize=22)

    fig.subplots_adjust(wspace=-0.1, hspace=0.15)

    return fig
