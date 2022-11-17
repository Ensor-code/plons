#Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import os
import math

# import own scripts
import SmoothingKernelScript    as sk
import ConversionFactors_cgs    as cgs
import PhysicalQuantities       as pq

# import certain things from packages
from matplotlib                 import colors
from astropy                    import constants
from mpl_toolkits.axes_grid1    import AxesGrid

from matplotlib                 import rcParams, rc
from mpl_toolkits.axes_grid1    import make_axes_locatable
from matplotlib.ticker          import MultipleLocator
import matplotlib
matplotlib.use("Agg")

#rc('font', family='serif')
#rc('text', usetex=True)

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

n_grid          = 500
n_grid_vec      = 25
mesh            = True
roundToInteger  = True
round_bounds    = False
customRanges    = True
velocity_vec    = True
printRanges     = True
minInfLog       = -300
minInf          = 10.**(-minInfLog)
maxInfLog       = 300
maxInf          = 10.**maxInfLog
sigma_bounds_u  = 3.
sigma_bounds_l  = 2.
nneighb = 20


'''
Load the smoothing kernel data
'''
def smoothData(dumpData, setup, theta, observables, zoom=1):
    print('          Calculating zoom = '+str(zoom), end='\r')
    results_sph_sl_z, x1, y1, z1  = sk.getSmoothingKernelledPix(n_grid, nneighb, dumpData, observables, 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2. / zoom, theta, mesh)
    results_sph_sl_z_vec, x1_vec, y1_vec, z1_vec  = sk.getSmoothingKernelledPix(n_grid_vec, nneighb, dumpData, ['vx', 'vy', 'vz'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2. / zoom, theta, mesh, vec=True)
    results_sph_sl_y, x2, y2, z2  = sk.getSmoothingKernelledPix(n_grid, nneighb, dumpData, observables, 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2. / zoom, theta, mesh)
    results_sph_sl_y_vec, x2_vec, y2_vec, z2_vec  = sk.getSmoothingKernelledPix(n_grid_vec, nneighb, dumpData, ['vx', 'vy', 'vz'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2. / zoom, theta, mesh, vec=True)

    if (setup['single_star']):
        xcomp, ycomp = 0,0
        smooth = {  'smooth_z'     :  results_sph_sl_z,
                    'x_z'          :  x1,
                    'y_z'          :  y1,
                    'smooth_y'     :  results_sph_sl_y,
                    'x_y'          :  x2,
                    'z_y'          :  z2
                    }

        smooth_vec = {  'smooth_z' :  results_sph_sl_z_vec,
                        'x_z'      :  x1_vec,
                        'y_z'      :  y1_vec,
                        'smooth_y' :  results_sph_sl_y_vec,
                        'x_y'      :  x2_vec,
                        'z_y'      :  z2_vec
                        }
    else:
        xcomp = dumpData['posComp'][0]
        ycomp = dumpData['posComp'][1]
        smooth = {  'smooth_z'     :  results_sph_sl_z,
                    'x_z'          :  x1,
                    'y_z'          :  y1,
                    'smooth_y'     :  results_sph_sl_y,
                    'x_y'          :  planeCoordinates(n_grid, x2, y2, xcomp, ycomp),
                    'z_y'          :  z2
                    }
        smooth_vec = {  'smooth_z' :  results_sph_sl_z_vec,
                        'x_z'      :  x1_vec,
                        'y_z'      :  y1_vec,
                        'smooth_y' :  results_sph_sl_y_vec,
                        'x_y'      :  planeCoordinates(n_grid_vec, x2_vec, y2_vec, xcomp, ycomp),
                        'z_y'      :  z2_vec
                        }

    return smooth, smooth_vec


'''
Make figure with the xy(orbital plane) slice plot of log density [g/cm^3]. 
    - smooth            is smoothing kernel data
    - zoom              zoom factor of the plot
    - rhoMin/rhoMax     determine the min and max of the colorbar
    - dumpData      data of dump file
    - setup         setup data
    - run           run number [str]
    - loc           output directory
    - rAccComp      accretion radius of the companion
'''
def densityPlot(smooth, zoom, limits, dumpData, setup, run, loc, rAccComp, rAccComp_in, number = -1):

    cm_rho  = plt.cm.get_cmap('inferno')
    fig, ax = plt.subplots(1, figsize=(7, 7))

    ax.set_aspect('equal')
    ax.set_facecolor('k')

    axPlot = None
    dataRho = np.log10(smooth[zoom]['smooth_z']["rho"])
    if mesh == False:
        axPlot = ax.scatter(smooth[zoom]['x_z']/cgs.AU_cm(), smooth[zoom]['y_z']/cgs.AU_cm(),
                            s=5, c=dataRho,cmap=cm_rho,vmin=limits[0], vmax = limits[1],
                            rasterized=True)
    else:
        axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(),
                               dataRho, cmap=cm_rho, vmin=limits[0], vmax=limits[1],
                               rasterized=True)

    if setup['single_star']== False:
        xAGB  = dumpData['posAGB' ][0] / cgs.AU_cm()
        yAGB  = dumpData['posAGB' ][1] / cgs.AU_cm()
        xcomp = dumpData['posComp'][0] / cgs.AU_cm()
        ycomp = dumpData['posComp'][1] / cgs.AU_cm()

        circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black", zorder=10)
        ax.add_artist(circleAGB)
        circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
        ax.add_artist(circleComp)
        
        if setup['triple_star']==True:
            xcomp_in = dumpData['posComp_in' ][0] / cgs.AU_cm()
            ycomp_in = dumpData['posComp_in' ][1] / cgs.AU_cm()
            circleComp_in = plt.Circle((xcomp_in, ycomp_in), rAccComp_in, transform=ax.transData._b, color="black", zorder=10)
            ax.add_artist(circleComp_in)
            

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="7%", pad=0.15)
    cbar = plt.colorbar(axPlot, cax=cax)
    cbar.set_label(r'$\log \, \rho$ [g$\,$cm$^{-3}$]', fontsize=24)
    cbar.ax.tick_params(labelsize=20)

    lim = (setup['bound'] * np.sqrt(2.)/2.) / zoom
    if roundToInteger: lim = math.floor(lim / 2.) * 2.
    ax.axis([-lim, lim, -lim + 0.001, lim])
    ax.xaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.xaxis.set_minor_locator(MultipleLocator(lim / 8.))
    ax.yaxis.set_major_locator(MultipleLocator(lim / 2.))
    ax.yaxis.set_minor_locator(MultipleLocator(lim / 8.))
    ax.set_xlabel(r"$x$ [AU]", fontsize=22)
    ax.set_ylabel(r"$y$ [AU]", fontsize=22)

    ax.tick_params(labelsize=20)

    if number == -1:
        fig.savefig(os.path.join(loc, 'png/2Dplot_density_zoom{0:01d}.png'.format(zoom)), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(loc, 'pdf/2Dplot_density_zoom{0:01d}.pdf'.format(zoom)), dpi=300, bbox_inches="tight")
    else:
        fig.text(0.5, 0.9, "Dumpfile {0:05d}".format(number), size=28)
        fig.savefig(os.path.join(loc, 'animation/2Dplot_density_zoom{0:01d}_{1:04d}.png'.format(zoom, int(
            number / setup['nfulldump']))), dpi=200,
                    bbox_inches="tight")
    plt.close()


'''
Makes one slice plot 

INPUT
    - ax        is given subplot
    - par       is the name of the parameter
    - mi/ma     colorbar limits
    - smooth    smoothed data
    - zoom      zoom factor for the plot
'''

def onePlot(fig, ax, par, limits, smooth, smooth_vec, zoom, dumpData, setup, axs, plane, rAccComp, rAccComp_in, velocity_vec = False):
    # colormap per parameter
    cm = {'rho': plt.cm.get_cmap('inferno'),
          'speed': plt.cm.get_cmap('CMRmap'),
          #'Tgas': plt.cm.get_cmap('RdYlGn'),
          'Tgas': plt.cm.get_cmap('hot'), #nipy_spectral
          'tau': plt.cm.get_cmap('viridis_r'),
          'kappa': plt.cm.get_cmap('Spectral_r'),
          'Gamma': plt.cm.get_cmap('Spectral_r')
          }

    # label name per parameter
    name = {'rho': r'$\log \, \rho$ [g$\,$cm$^{-3}$]',
            'speed': '$v$ [km/s]',
            'Tgas': r'$\log \, T$ [K]',
            'tau': r'$\tau [/]$',
            'kappa': r'$\kappa$ [g/cm$^3$]',
            'Gamma': r'$\Gamma$ [/]'
            }
    logtau = False

    xAGB = dumpData['posAGB'][0] / cgs.AU_cm()
    yAGB = dumpData['posAGB'][1] / cgs.AU_cm()
    zAGB = dumpData['posAGB'][2] / cgs.AU_cm()


    if setup['single_star'] == False:
        xcomp = dumpData['posComp'][0] / cgs.AU_cm()
        ycomp = dumpData['posComp'][1] / cgs.AU_cm()
        zcomp = dumpData['posComp'][2] / cgs.AU_cm()
        if setup['triple_star']==True:
            xcomp_in = dumpData['posComp_in'][0] / cgs.AU_cm()
            ycomp_in = dumpData['posComp_in'][1] / cgs.AU_cm()
            zcomp_in = dumpData['posComp_in'][2] / cgs.AU_cm()
            

    gamma = setup["gamma"]
    axPlot = None
    if plane == 'z':
        if (not logtau and par == 'tau') or par in ('kappa', 'Gamma', 'speed'):
            data = smooth[zoom]['smooth_z'][par]
        else:
            data = np.log10(smooth[zoom]['smooth_z'][par])
            
        if mesh == False:
            axPlot = ax.scatter(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(), s=5,
                                c=data, cmap=cm[par], vmin=limits[0], vmax=limits[1],
                                rasterized=True)
        else:
            axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(),
                                    data, cmap=cm[par], vmin=limits[0], vmax=limits[1],
                                    rasterized=True, shading="nearest")

            if par == 'speed' and velocity_vec:
                vx = smooth_vec[zoom]['smooth_z']['vx']
                vy = smooth_vec[zoom]['smooth_z']['vy']
                normaliseVectorLength = np.hypot(vx, vy)
                ax.quiver(smooth_vec[zoom]['x_z'] / cgs.AU_cm(), smooth_vec[zoom]['y_z'] / cgs.AU_cm(),
                            vx / normaliseVectorLength, vy / normaliseVectorLength, scale_units="dots", scale=0.05)

    if plane == 'y':
        if (not logtau and par == 'tau') or par in ('kappa', 'Gamma', 'speed'):
            data = smooth[zoom]['smooth_y'][par]
        else:
            data = np.log10(smooth[zoom]['smooth_y'][par])
        if mesh == False:
            axPlot = ax.scatter(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(), s=5,
                                c=data, cmap=cm[par], vmin=limits[0], vmax=limits[1],
                                rasterized=True)

        else:
            axPlot = ax.pcolormesh(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(), data,
                                    cmap=cm[par], vmin=limits[0], vmax=limits[1], rasterized=True, shading="nearest")

            if par == 'speed' and velocity_vec:
                vx = smooth_vec[zoom]['smooth_y']['vx']
                vz = smooth_vec[zoom]['smooth_y']['vz']
                normaliseVectorLength = np.hypot(vx, vz)
                ax.quiver(smooth_vec[zoom]['x_y'] / cgs.AU_cm(), smooth_vec[zoom]['z_y'] / cgs.AU_cm(),
                            vx / normaliseVectorLength, vz / normaliseVectorLength, scale_units="dots", scale=0.05)

    # plot the position of the AGB star and comp in the edge-on plane & make colorbar
    if gamma <= 1.:
        if ax == axs[2] or ax == axs[4]:
            ax.set_ylabel(r"$z$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black",
                                   zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b,
                                        color="black", zorder=10)
                ax.add_artist(circleComp)
                if setup['triple_star']==True:
                    circleComp_in = plt.Circle((xcomp_in, zcomp_in), rAccComp_in, transform=ax.transData._b,
                                        color="black", zorder=10)
                    ax.add_artist(circleComp_in)                    

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="7%", pad=0.15)
            cbar = fig.colorbar(axPlot, cax=cax)
            cbar.set_label(name[par], fontsize=24)
            cbar.ax.tick_params(labelsize=20)

        # plot the position of the AGB star and comp in the face-on plane
        if ax == axs[1] or ax == axs[3]:
            ax.set_ylabel(r"$y$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black",
                                   zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b,
                                        color="black", zorder=10)
                ax.add_artist(circleComp)
                if setup['triple_star']==True:
                    pos = np.array([xcomp_in,ycomp_in,zcomp_in])
                    circleComp_in = plt.Circle((xcomp_in, ycomp_in), rAccComp_in, transform=ax.transData._b,
                                        color="black", zorder=10)
                    ax.add_artist(circleComp_in) 
                
                

    else:
        if ax == axs[2] or ax == axs[4] or ax == axs[6]:
            ax.set_ylabel(r"$z$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black", zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
                ax.add_artist(circleComp)
                if setup['triple_star']==True:
                    circleComp_in = plt.Circle((xcomp_in, zcomp_in), rAccComp_in, transform=ax.transData._b,
                                        color="black", zorder=10)
                    ax.add_artist(circleComp_in)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="7%", pad=0.15)
            cbar = fig.colorbar(axPlot, cax=cax)
            cbar.set_label(name[par], fontsize=24)
            cbar.ax.tick_params(labelsize=20)

        # plot the position of the AGB star and comp in the face-on plane
        if ax == axs[1] or ax == axs[3] or ax == axs[5]:
            ax.set_ylabel(r"$y$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black", zorder=10)
            #circleAGB = plt.Circle((xAGB, yAGB), max(setup["wind_inject_radius"],setup["primary_Reff"]), transform=ax.transData._b, color="black", zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
                #circleComp = plt.Circle((xcomp, ycomp), rAccComp, transform=ax.transData._b, color="black", zorder=10)
                ax.add_artist(circleComp)
                if setup['triple_star']==True:
                    circleComp_in = plt.Circle((xcomp_in, ycomp_in), rAccComp_in, transform=ax.transData._b,
                                        color="black", zorder=10)
                    ax.add_artist(circleComp_in)

    ax.tick_params(labelsize=20)
    ax.set_aspect('equal')
    ax.set_facecolor('k')

    lim = (setup['bound'] * np.sqrt(2.) / 2.) / zoom
    if roundToInteger: lim = math.floor(lim / 2.) * 2.
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


def allPlots(smooth, smooth_vec, zoom, limits, dumpData, setup, run, loc, rAccComp, rAccComp_in, observables, number = - 1):
    gamma = setup["gamma"]

    fig = None
    if gamma <= 1.:
        fig = plt.figure(figsize=(10 + 0.35*5, 10))
        spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[1., 1.], height_ratios=[1, 1])

        ax1 = fig.add_subplot(spec[0, 0])
        ax2 = fig.add_subplot(spec[0, 1])
        ax3 = fig.add_subplot(spec[1, 0])
        ax4 = fig.add_subplot(spec[1, 1])

        axs = {1: ax1,
               2: ax2,
               3: ax3,
               4: ax4
               }

        # the temperature colorbar limits may have to be changed...
        onePlot(fig, ax1, 'rho', limits["rho"][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, rAccComp_in)
        onePlot(fig, ax2, 'rho', limits["rho"][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, rAccComp_in)
        onePlot(fig, ax3, 'speed', limits["speed"][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, rAccComp_in, velocity_vec = velocity_vec)
        onePlot(fig, ax4, 'speed', limits["speed"][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, rAccComp_in, velocity_vec = velocity_vec)

        ax2.set_yticklabels([])
        ax4.set_yticklabels([])

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xlabel(r"$x$ [AU]", fontsize=22)
        ax4.set_xlabel(r"$x$ [AU]", fontsize=22)

    else:
        fig = plt.figure(figsize=(10 + 0.35*5, 15))
        spec = fig.add_gridspec(ncols=2, nrows=len(observables), width_ratios=[1., 1.], height_ratios=[1.]*len(observables))

        ax1 = fig.add_subplot(spec[0, 0])
        ax2 = fig.add_subplot(spec[0, 1])
        ax3 = fig.add_subplot(spec[1, 0])
        ax4 = fig.add_subplot(spec[1, 1])
        ax5 = fig.add_subplot(spec[2, 0])
        ax6 = fig.add_subplot(spec[2, 1])

        axs = {1: ax1,
               2: ax2,
               3: ax3,
               4: ax4,
               5: ax5,
               6: ax6
               }

        for i in range(len(observables)):
            if observables[i] == 'speed':
                onePlot(fig, axs[2*i+1], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, rAccComp_in, velocity_vec = velocity_vec)
                onePlot(fig, axs[2*i+2], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, rAccComp_in, velocity_vec = velocity_vec)
                axs[2*i+2].set_yticklabels([])
            else:
                onePlot(fig, axs[2*i+1], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, rAccComp_in)
                onePlot(fig, axs[2*i+2], observables[i], limits[observables[i]][zoom], smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, rAccComp_in)
                axs[2*i+2].set_yticklabels([])

        for i in range(len(observables)-1):
            axs[2*i+1].set_xticklabels([])
            axs[2*i+2].set_xticklabels([])
        axs[2*len(observables)-1].set_xlabel(r"$x$ [AU]", fontsize=22)
        axs[2*len(observables)  ].set_xlabel(r"$x$ [AU]", fontsize=22)

    fig.subplots_adjust(wspace=-0.1, hspace=0.15)
    if number == -1:
        plt.savefig(os.path.join(loc, 'png/2Dplot_DensTempTau_zoom{0:01d}.png'.format(zoom)), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(loc, 'pdf/2Dplot_DensTempTau_zoom{0:01d}.pdf'.format(zoom)), bbox_inches="tight")
        print('          Slice plots (zoom factor = ' + str(zoom) + ') model ' + str(run) + ' ready and saved!\n')

    else:
        fig.text(0.5, 0.9, "Dumpfile {0:05d}".format(number), size=28)
        plt.savefig(os.path.join(loc, 'animation/2Dplot_DensSpeedTau_zoom{0:01d}_{1:04d}.png'.format(zoom, int(number/setup['nfulldump']))), dpi=200,
                    bbox_inches="tight")
    plt.close()

'''
main definition
'''


def SlicePlots(run, loc, dumpData, setup, number = -1, zoomin = [1,2,5,10], observables = ['rho', 'Tgas', 'speed']):
    print('')
    print('(1)  Start calculations for slice plots...')

    rAccComp_in = 0
    theta=0
    rAccComp = 0
    if not setup["single_star"]:
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
        rAccComp = setup['rAccrComp']
        if rAccComp <= 0.05 * max(setup["wind_inject_radius"],setup["primary_Reff"]):
            rAccComp = 0.05 * max(setup["wind_inject_radius"],setup["primary_Reff"])
        if setup['triple_star']:
            rAccComp_in = setup['rAccrComp_in']

    limits = {}
    for observable in observables:
        limits[observable] = {}
        for zoom in zoomin:
            limits[observable][zoom] = [0.,0.]
            
    # Ranges are the same as for '/gamma/non_isowind/gamma1.4/global'
    if "global" in run:
        if "rho" in observables:
            limits["rho"][1]  = [-24.39280, -16.86249]
            limits["rho"][2]  = [-20.80625, -17.63237]
            limits["rho"][5]  = [-19.71264, -16.90193]
            limits["rho"][10] = [-19.71264, -16.90193]

        if "speed" in observables:
            limits["speed"][1]  = [5.96568, 38.44218]
            limits["speed"][2]  = [5.99680, 33.89827]
            limits["speed"][5]  = [3.17029, 25.16178]
            limits["speed"][10] = [3.17029, 25.16178]

        if "Tgas" in observables:
            limits["Tgas"][1]  = [1.26680, 4.65751]
            limits["Tgas"][2]  = [2.58142, 4.07818]
            limits["Tgas"][5]  = [2.52443, 4.53982]
            limits["Tgas"][10] = [2.52443, 4.53982]

        if "tau" in observables:
            limits["tau"][1]  = [0, 0.17]
            limits["tau"][2]  = [0, 0.17]
            limits["tau"][5]  = [0, 0.17]
            limits["tau"][10] = [0, 0.17]

        if "kappa" in observables:
            limits["kappa"][1]  = [0., 3.]
            limits["kappa"][2]  = [0., 3.]
            limits["kappa"][5]  = [0., 3.]
            limits["kappa"][10] = [0., 3.]

        if "Gamma" in observables:
            limits["Gamma"][1]  = [0., 1.]
            limits["Gamma"][2]  = [0., 1.]
            limits["Gamma"][5]  = [0., 1.]
            limits["Gamma"][10] = [0., 1.]

    # Ranges are the same as for '/gamma/non_isowind/gamma1.4/local'
    elif "local" in run:
        if "rho" in observables:
            limits["rho"][1]  = [-19.96524, -16.86532]
            limits["rho"][2]  = [-19.57196, -15.85257]
            limits["rho"][5]  = [-18.47570, -14.90456]
            limits["rho"][10] = [-18.47570, -14.90456]

        if "speed" in observables:
            limits["speed"][1]  = [2.54833, 21.86302]
            limits["speed"][2]  = [0.03749, 18.39068]
            limits["speed"][5]  = [0.00000, 22.53676]
            limits["speed"][10] = [0.00000, 22.53676]

        if "Tgas" in observables:
            limits["Tgas"][1]  = [2.26195, 4.49702]
            limits["Tgas"][2]  = [2.14475, 4.92498]
            limits["Tgas"][5]  = [1.97842, 5.32584]
            limits["Tgas"][10] = [1.97842, 5.32584]

        if "tau" in observables:
            limits["tau"][1]  = [0, 0.17]
            limits["tau"][2]  = [0, 0.17]
            limits["tau"][5]  = [0, 0.17]
            limits["tau"][10] = [0, 0.17]

        if "kappa" in observables:
            limits["kappa"][1]  = [0., 3.]
            limits["kappa"][2]  = [0., 3.]
            limits["kappa"][5]  = [0., 3.]
            limits["kappa"][10] = [0., 3.]

        if "Gamma" in observables:
            limits["Gamma"][1]  = [0., 1.]
            limits["Gamma"][2]  = [0., 1.]
            limits["Gamma"][5]  = [0., 1.]
            limits["Gamma"][10] = [0., 1.]
    
    elif customRanges:
        if "rho" in observables:
            limits["rho"][1]  = [-21, -14]
            limits["rho"][2]  = [-21, -14]
            limits["rho"][5]  = [-21, -14]
            limits["rho"][10] = [-21, -14]

        if "speed" in observables:
            limits["speed"][1]  = [0., 20.]
            limits["speed"][2]  = [0., 20.]
            limits["speed"][5]  = [0., 20.]
            limits["speed"][10] = [0., 20.]

        if "Tgas" in observables:
            limits["Tgas"][1]  = [1., 4.]
            limits["Tgas"][2]  = [1., 4.]
            limits["Tgas"][5]  = [1., 4.]
            limits["Tgas"][10] = [1., 4.]

        if "tau" in observables:
            limits["tau"][1]  = [0, 1]
            limits["tau"][2]  = [0, 1]
            limits["tau"][5]  = [0, 1]
            limits["tau"][10] = [0, 1]

        if "kappa" in observables:
            limits["kappa"][1]  = [0., 3.]
            limits["kappa"][2]  = [0., 3.]
            limits["kappa"][5]  = [0., 3.]
            limits["kappa"][10] = [0., 3.]

        if "Gamma" in observables:
            limits["Gamma"][1]  = [0., 1.]
            limits["Gamma"][2]  = [0., 1.]
            limits["Gamma"][5]  = [0., 1.]
            limits["Gamma"][10] = [0., 1.]

    print('     Calculating the smoothing kernels. This may take a while, please wait...')
    smooth = {}
    smooth_vec = {}
    for zoom in zoomin:
        smooth[zoom], smooth_vec[zoom] = smoothData(dumpData, setup, theta, observables, zoom)

        if not "global" in run and not "local" in run and not customRanges:
            if "rho" in observables:
                limits["rho"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_y']["rho"]), log=True, round=round_bounds)
            if "speed" in observables:
                limits["speed"][zoom] = findBounds(smooth[zoom]['smooth_y']["speed"], log=False, round=round_bounds)
                limits["speed"][zoom][0] = max(limits["speed"][zoom][0], 0.)
            if "Tgas" in observables:
                limits["Tgas"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_z']["Tgas"]), log=True, round=round_bounds)
            if "tau" in observables:
                limits["tau"][zoom] = findBounds(smooth[zoom]['smooth_y']["tau"], log=True, round=round_bounds)
                limits["tau"][zoom][0] = max(limits["tau"][zoom][0], 0.)
            if "tau" in observables:
                limits["kappa"][zoom] = findBounds(smooth[zoom]['smooth_y']["kappa"], log=False, round=round_bounds)
            if "tau" in observables:
                limits["Gamma"][zoom] = findBounds(smooth[zoom]['smooth_y']["Gamma"], log=False, round=round_bounds)

        if printRanges:
            print("          Ranges of Parameters: zoom = "+str(zoom))
            if "rho"   in observables: print("          rhoMin,   rhoMax   = {0:10.5f}, {1:10.5f}".format(limits["rho"][zoom][0], limits["rho"][zoom][1]))
            if "speed" in observables: print("          vMin,     vMax     = {0:10.5f}, {1:10.5f}".format(limits["speed"][zoom][0], limits["speed"][zoom][1]))
            if "Tgas"  in observables: print("          TMin,     TMax     = {0:10.5f}, {1:10.5f}".format(limits["Tgas"][zoom][0], limits["Tgas"][zoom][1]))
            if "kappa" in observables: print("          kappaMin, kappaMax = {0:10.5f}, {1:10.5f}".format(limits["kappa"][zoom][0], limits["kappa"][zoom][1]))
            if "Gamma" in observables: print("          GammaMin, GammaMax = {0:10.5f}, {1:10.5f}".format(limits["Gamma"][zoom][0], limits["Gamma"][zoom][1]))
            if "tau"   in observables: print("          tauMin,   tauMax   = {0:10.5f}, {1:10.5f}".format(limits["tau"][zoom][0], limits["tau"][zoom][1]))

        # Make plots
        densityPlot(smooth, zoom, limits["rho"][zoom], dumpData, setup, run, loc, rAccComp, rAccComp_in, number = number)
        allPlots(smooth, smooth_vec, zoom, limits, dumpData, setup, run, loc, rAccComp, rAccComp_in, observables, number = number)

    
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

def planeCoordinates(n, x, y, x_comp, y_comp):
    r = np.zeros(shape=(n, n))
    dot_product = x * x_comp + y * y_comp
    r[dot_product < 0] = -np.hypot(x[dot_product < 0], y[dot_product < 0])
    r[dot_product >= 0] = np.hypot(x[dot_product >= 0], y[dot_product >= 0])

    return r
