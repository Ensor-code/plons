#Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import os
import math

# import plons scripts
import plons.SmoothingKernelScript    as sk
import plons.ConversionFactors_cgs    as cgs
import plons.PhysicalQuantities       as pq
import plons.Plotting                 as plot

# import certain things from packages
from mpl_toolkits.axes_grid1    import make_axes_locatable
from matplotlib.ticker          import MultipleLocator

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

roundToInteger  = True
round_bounds    = False
velocity_vec    = True
minInfLog       = -300
minInf          = 10.**(-minInfLog)
maxInfLog       = 300
maxInf          = 10.**maxInfLog
sigma_bounds_u  = 3.
sigma_bounds_l  = 2.


# colormap per parameter
cm =   {'rho':   plt.cm.get_cmap('inferno'),
        'speed': plt.cm.get_cmap('CMRmap'),
        #'Tgas': plt.cm.get_cmap('RdYlGn'),
        'Tgas':  plt.cm.get_cmap('hot'), #nipy_spectral
        'tau':   plt.cm.get_cmap('viridis_r'),
        'tauL':  plt.cm.get_cmap('viridis'),
        'kappa': plt.cm.get_cmap('Spectral_r'),
        'Tdust': plt.cm.get_cmap('Spectral_r'),
        'Gamma': plt.cm.get_cmap('Spectral_r')
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
def densityPlot(smooth, zoom, limits, dumpData, setup, orbital=True, cmap=plt.cm.get_cmap('inferno')):
    
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


def makeLimits(observables, smooth, zoom):
    limits = {}
    for observable in observables:
        limits[observable] = {}
        limits[observable][zoom] = [0.,0.]

    if "rho" in observables:
        limits["rho"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_y']["rho"]), log=True, round=round_bounds)
    if "speed" in observables:
        limits["speed"][zoom] = findBounds(smooth[zoom]['smooth_y']["speed"], log=False, round=round_bounds)
        limits["speed"][zoom][0] = max(limits["speed"][zoom][0], 0.)
    if "Tgas" in observables:
        limits["Tgas"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_z']["Tgas"]), log=True, round=round_bounds)
    if "kappa" in observables:
        limits["kappa"][zoom] = findBounds(smooth[zoom]['smooth_y']["kappa"], log=False, round=round_bounds)
    if "Gamma" in observables:
        limits["Gamma"][zoom] = findBounds(smooth[zoom]['smooth_y']["Gamma"], log=False, round=round_bounds)
    if "tau" in observables:
        limits["tau"][zoom] = findBounds(smooth[zoom]['smooth_y']["tau"], log=True, round=round_bounds)
        limits["tau"][zoom][0] = max(limits["tau"][zoom][0], 0.)

    return limits

'''
main definition
'''
def SlicePlots(run, loc, dumpData, setup, number = -1, zoomin = [1, 5, 10], 
               observables = ['rho', 'Tgas', 'speed'], limits = False,
               nneighb = 10, n_grid = 200, n_grid_vec = 25, printout=False):
    if limits: customRanges = True
    else:      customRanges = False
    theta=0
    if not setup["single_star"]:
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
    
    if printout: print('     Calculating the smoothing kernels. This may take a while, please wait...')
    smooth = {}
    smooth_vec = {}
    for zoom in zoomin:
        smooth[zoom], smooth_vec[zoom] = smoothData(dumpData, setup, observables, theta, zoom, nneighb, n_grid, n_grid_vec)

        if not customRanges: 
            limits = makeLimits(observables, smooth, zoom)    
        if printout:
            print("          Ranges of Parameters: zoom = "+str(zoom))
            if "rho"   in observables: print("          rhoMin,   rhoMax   = {0:10.5f}, {1:10.5f}".format(limits["rho"][zoom][0], limits["rho"][zoom][1]))
            if "speed" in observables: print("          vMin,     vMax     = {0:10.5f}, {1:10.5f}".format(limits["speed"][zoom][0], limits["speed"][zoom][1]))
            if "Tgas"  in observables: print("          TMin,     TMax     = {0:10.5f}, {1:10.5f}".format(limits["Tgas"][zoom][0], limits["Tgas"][zoom][1]))
            if "kappa" in observables: print("          kappaMin, kappaMax = {0:10.5f}, {1:10.5f}".format(limits["kappa"][zoom][0], limits["kappa"][zoom][1]))
            if "Gamma" in observables: print("          GammaMin, GammaMax = {0:10.5f}, {1:10.5f}".format(limits["Gamma"][zoom][0], limits["Gamma"][zoom][1]))
            if "tau"   in observables: print("          tauMin,   tauMax   = {0:10.5f}, {1:10.5f}".format(limits["tau"][zoom][0], limits["tau"][zoom][1]))

        # Make plots
        fig = densityPlot(smooth, zoom, limits["rho"][zoom], dumpData, setup, orbital=True)
        saveFig(fig, loc, '2Dplot_density_orbital', zoom, number)
        fig = densityPlot(smooth, zoom, limits["rho"][zoom], dumpData, setup, orbital=False)
        saveFig(fig, loc, '2Dplot_density_meridional', zoom, number)
        fig = allPlots(smooth, smooth_vec, zoom, limits, dumpData, setup, observables)
        struct = ""
        for i in observables:
            struct+=i
        saveFig(fig, loc, '2Dplot_'+struct, zoom, number)
        plt.close()
        if printout: print('          Slice plots (zoom factor = ' + str(zoom) + ') model ' + str(run) + ' ready and saved!\n')

def saveFig(fig, loc, name, zoom, number):
    if number == -1:
        fig.savefig(os.path.join(loc, 'png/'+name+'_zoom{0:01d}.png'.format(zoom)), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(loc, 'pdf/'+name+'_zoom{0:01d}.pdf'.format(zoom)), dpi=300, bbox_inches="tight")
    else:
        fig.text(0.5, 0.9, "Dumpfile {0:05d}".format(number), size=28)
        fig.savefig(os.path.join(loc, 'animation/'+name+'_zoom{0:01d}_{1:04d}.png'.format(zoom,number)), dpi=200, bbox_inches="tight")
    plt.close()
    
def findBounds(data, log = False, round = False):
    filtered_data = data
    result = np.zeros(2)

    if (-np.inf in data) or (np.inf in data):
        min = minInf
        max = maxInf
        if log == True:
            min = minInfLog
            max = maxInfLog

        filtered_data = filtered_data[np.logical_and(min <= data, data <= max)]

    if np.nan in data:
        filtered_data = filtered_data[not np.isnan(data)]

    if (0. in data) and (log == False) :
        filtered_data = filtered_data[filtered_data != 0.]

    mean = np.mean(filtered_data)
    std  = np.std(filtered_data)

    result[0] = mean - sigma_bounds_l * std
    result[1] = mean + sigma_bounds_u * std

    if round == True:
        result = np.round(result)

    return result