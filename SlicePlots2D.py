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

radius_AGB      = 1.15 #AU
n_grid          = 500
n_grid_vec      = 50
mesh            = True
roundToInteger  = True
round_bounds    = False
customRanges    = False
velocity_vec    = True
printRanges     = False
minInfLog       = -300
minInf          = 10.**(-minInfLog)
maxInfLog       = 300
maxInf          = 10.**maxInfLog
sigma_bounds_u  = 3.
sigma_bounds_l  = 2.


'''
Load the smoothing kernel data
'''
def smoothData(dumpData, setup, theta):
    #comment the ones you don't need, it takes long to run, so don't run if not nescessary
    print('     Calculating the smoothing kernels. This may take a while, please wait...')
    #zoom = 1
    results_sph_sl_z, x1, y1, z1  = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2., theta, mesh)
    results_sph_sl_z_vec, x1_vec, y1_vec, z1_vec  = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2., theta, mesh, vec=True)
    results_sph_sl_y, x2, y2, z2  = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2., theta, mesh)
    results_sph_sl_y_vec, x2_vec, y2_vec, z2_vec  = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 2., theta, mesh, vec=True)

    #zoom = 2
    results_sph_sl_zZ2,x1Z2,y1Z2,z1Z2 = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 4., theta, mesh)
    results_sph_sl_zZ2_vec, x1Z2_vec, y1Z2_vec, z1Z2_vec = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 4., theta, mesh, vec=True)
    results_sph_sl_yZ2,x2Z2,y2Z2,z2Z2 = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 4., theta, mesh)
    results_sph_sl_yZ2_vec, x2Z2_vec, y2Z2_vec, z2Z2_vec = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 4., theta, mesh, vec=True)

    #zoom = 5
    results_sph_sl_zZ5,x1Z5,y1Z5,z1Z5 = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 10., theta, mesh)
    results_sph_sl_zZ5_vec, x1Z5_vec, y1Z5_vec, z1Z5_vec = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'z', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 10., theta, mesh, vec=True)
    results_sph_sl_yZ5,x2Z5,y2Z5,z2Z5 = sk.getSmoothingKernelledPix(n_grid, 20, dumpData, ['rho', 'temp', 'speed'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 10., theta, mesh)
    results_sph_sl_yZ5_vec, x2Z5_vec, y2Z5_vec, z2Z5_vec = sk.getSmoothingKernelledPix(n_grid_vec, 20, dumpData, ['vx', 'vy', 'vz'], 'comp', 'y', (setup['bound']) * cgs.AU_cm() * np.sqrt(2.) / 10., theta, mesh, vec=True)

    xcomp = dumpData['posComp'][0]
    ycomp = dumpData['posComp'][1]
    smooth_zoom1 = {'smooth_z'     :  results_sph_sl_z,
                    'x_z'          :  x1,
                    'y_z'          :  y1,
                    'smooth_y'     :  results_sph_sl_y,
                    'x_y'          :  planeCoordinates(n_grid, x2, y2, xcomp, ycomp),
                    'z_y'          :  z2
                    }

    smooth_zoom2 = {'smooth_z'     :  results_sph_sl_zZ2,
                    'x_z'          :  x1Z2,
                    'y_z'          :  y1Z2,
                    'smooth_y'     :  results_sph_sl_yZ2,
                    'x_y'          :  planeCoordinates(n_grid, x2Z2, y2Z2, xcomp, ycomp),
                    'z_y'          :  z2Z2
                    }

    smooth_zoom5 = {'smooth_z'     :  results_sph_sl_zZ5,
                    'x_z'          :  x1Z5,
                    'y_z'          :  y1Z5,
                    'smooth_y'     :  results_sph_sl_yZ5,
                    'x_y'          :  planeCoordinates(n_grid, x2Z5, y2Z5, xcomp, ycomp),
                    'z_y'          :  z2Z5
                    }

    smooth_zoom1_vec = {'smooth_z': results_sph_sl_z_vec,
                        'x_z': x1_vec,
                        'y_z': y1_vec,
                        'smooth_y': results_sph_sl_y_vec,
                        'x_y': planeCoordinates(n_grid_vec, x2_vec, y2_vec, xcomp, ycomp),
                        'z_y': z2_vec
                        }

    smooth_zoom2_vec = {'smooth_z': results_sph_sl_zZ2_vec,
                        'x_z': x1Z2_vec,
                        'y_z': y1Z2_vec,
                        'smooth_y': results_sph_sl_yZ2_vec,
                        'x_y': planeCoordinates(n_grid_vec, x2Z2_vec, y2Z2_vec, xcomp, ycomp),
                        'z_y': z2Z2_vec
                        }

    smooth_zoom5_vec = {'smooth_z': results_sph_sl_zZ5_vec,
                        'x_z': x1Z5_vec,
                        'y_z': y1Z5_vec,
                        'smooth_y': results_sph_sl_yZ5_vec,
                        'x_y': planeCoordinates(n_grid_vec, x2Z5_vec, y2Z5_vec, xcomp, ycomp),
                        'z_y': z2Z5_vec
                        }

    smooth = { 1: smooth_zoom1,
               2: smooth_zoom2,
               5: smooth_zoom5
               }
    smooth_vec = { 1: smooth_zoom1_vec,
                   2: smooth_zoom2_vec,
                   5: smooth_zoom5_vec
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
def densityPlot(smooth, zoom, rhoMin, rhoMax, dumpData, setup, run, loc, rAccComp, number = -1):

    cm_rho  = plt.cm.get_cmap('inferno')
    fig, ax = plt.subplots(1, figsize=(7, 7))

    ax.set_aspect('equal')
    ax.set_facecolor('k')

    axPlot = None
    dataRho = np.log10(smooth[zoom]['smooth_z']["rho"])
    if mesh == False:
        axPlot = ax.scatter(smooth[zoom]['x_z']/cgs.AU_cm(), smooth[zoom]['y_z']/cgs.AU_cm(),
                            s=5, c=dataRho,cmap=cm_rho,vmin=rhoMin, vmax = rhoMax,
                            rasterized=True)
    else:
        axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(),
                               dataRho, cmap=cm_rho, vmin=rhoMin, vmax=rhoMax,
                               rasterized=True)

    if setup['single_star']== False:
        xAGB  = dumpData['posAGB' ][0] / cgs.AU_cm()
        yAGB  = dumpData['posAGB' ][1] / cgs.AU_cm()
        xcomp = dumpData['posComp'][0] / cgs.AU_cm()
        ycomp = dumpData['posComp'][1] / cgs.AU_cm()

        circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), radius_AGB, transform=ax.transData._b, color="black", zorder=10)
        circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
        ax.add_artist(circleAGB)
        ax.add_artist(circleComp)

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

def onePlot(fig, ax, par, mi, ma, smooth, smooth_vec, zoom, dumpData, setup, axs, plane, rAccComp, velocity_vec = False):
    # colormap per parameter
    cm = {'rho': plt.cm.get_cmap('inferno'),
          'speed': plt.cm.get_cmap('Spectral'),
          'temp': plt.cm.get_cmap('RdYlGn')
          }

    # label name per parameter
    name = {'rho': r'$\log \, \rho$ [g$\,$cm$^{-3}$]',
            'temp': r'$\log \, T$ [K]',
            'speed': '$v$ [km/s]'
            }

    xAGB = dumpData['posAGB'][0] / cgs.AU_cm()
    yAGB = dumpData['posAGB'][1] / cgs.AU_cm()
    zAGB = dumpData['posAGB'][2] / cgs.AU_cm()

    if setup['single_star'] == False:
        xcomp = dumpData['posComp'][0] / cgs.AU_cm()
        ycomp = dumpData['posComp'][1] / cgs.AU_cm()
        zcomp = dumpData['posComp'][2] / cgs.AU_cm()

    gamma = setup["gamma"]
    axPlot = None
    if not par == 'speed':
        if plane == 'z':
            data = np.log10(smooth[zoom]['smooth_z'][par])
            if mesh == False:
                axPlot = ax.scatter(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(), s=5,
                                    c=data, cmap=cm[par], vmin=mi, vmax=ma,
                                    rasterized=True)
            else:
                axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(),
                                       data, cmap=cm[par], vmin=mi, vmax=ma,
                                       rasterized=True, shading="nearest")

        if plane == 'y':
            data = np.log10(smooth[zoom]['smooth_y'][par])
            if mesh == False:
                axPlot = ax.scatter(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(), s=5,
                                    c=data, cmap=cm[par], vmin=mi, vmax=ma,
                                    rasterized=True)

            else:
                axPlot = ax.pcolormesh(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(), data,
                                       cmap=cm[par], vmin=mi, vmax=ma, rasterized=True, shading="nearest")


    if par == 'speed':
        if plane == 'z':
            data = smooth[zoom]['smooth_z'][par] * cgs.cms_kms()
            if mesh == False:
                axPlot = ax.scatter(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(), s=5,
                                    c=data, cmap=cm[par], vmin=mi, vmax=ma,
                                    rasterized=True)
            else:
                axPlot = ax.pcolormesh(smooth[zoom]['x_z'] / cgs.AU_cm(), smooth[zoom]['y_z'] / cgs.AU_cm(),
                                       data, cmap=cm[par], vmin=mi, vmax=ma,
                                       rasterized=True, shading="nearest")

                if velocity_vec:
                    vx = smooth_vec[zoom]['smooth_z']['vx'] * cgs.cms_kms()
                    vy = smooth_vec[zoom]['smooth_z']['vy'] * cgs.cms_kms()
                    normaliseVectorLength = np.hypot(vx, vy)
                    ax.quiver(smooth_vec[zoom]['x_z'] / cgs.AU_cm(), smooth_vec[zoom]['y_z'] / cgs.AU_cm(),
                              vx / normaliseVectorLength, vy / normaliseVectorLength, scale_units="dots", scale=0.2)

        if plane == 'y':
            data = smooth[zoom]['smooth_y'][par] * cgs.cms_kms()
            if mesh == False:
                axPlot = ax.scatter(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(), s=5,
                                    c=data, cmap=cm[par], vmin=mi, vmax=ma,
                                    rasterized=True)
            else:
                axPlot = ax.pcolormesh(smooth[zoom]['x_y'] / cgs.AU_cm(), smooth[zoom]['z_y'] / cgs.AU_cm(),
                                       data, cmap=cm[par], vmin=mi, vmax=ma,
                                       rasterized=True, shading="nearest")

                if velocity_vec != None:
                    vx = smooth_vec[zoom]['smooth_y']['vx'] * cgs.cms_kms()
                    vz = smooth_vec[zoom]['smooth_y']['vz'] * cgs.cms_kms()
                    normaliseVectorLength = np.hypot(vx, vz)
                    ax.quiver(smooth_vec[zoom]['x_y'] / cgs.AU_cm(), smooth_vec[zoom]['z_y'] / cgs.AU_cm(),
                              vx / normaliseVectorLength, vz / normaliseVectorLength, scale_units="dots", scale=0.2)


    # plot the position of the AGB star and comp in the edge-on plane & make colorbar
    if gamma <= 1.:
        if ax == axs[2] or ax == axs[4]:
            ax.set_ylabel(r"$z$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), radius_AGB, transform=ax.transData._b, color="black",
                                   zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b,
                                        color="black", zorder=10)
                ax.add_artist(circleComp)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="7%", pad=0.15)
            cbar = fig.colorbar(axPlot, cax=cax)
            cbar.set_label(name[par], fontsize=24)
            cbar.ax.tick_params(labelsize=20)

        # plot the position of the AGB star and comp in the face-on plane
        if ax == axs[1] or ax == axs[3]:
            ax.set_ylabel(r"$y$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), radius_AGB, transform=ax.transData._b, color="black",
                                   zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b,
                                        color="black", zorder=10)
                ax.add_artist(circleComp)

    else:
        if ax == axs[2] or ax == axs[4] or ax == axs[6]:
            ax.set_ylabel(r"$z$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), radius_AGB, transform=ax.transData._b, color="black", zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
                ax.add_artist(circleComp)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="7%", pad=0.15)
            cbar = fig.colorbar(axPlot, cax=cax)
            cbar.set_label(name[par], fontsize=24)
            cbar.ax.tick_params(labelsize=20)

        # plot the position of the AGB star and comp in the face-on plane
        if ax == axs[1] or ax == axs[3] or ax == axs[5]:
            ax.set_ylabel(r"$y$ [AU]", fontsize=22)
            circleAGB = plt.Circle((-np.hypot(xAGB, yAGB), 0.), radius_AGB, transform=ax.transData._b, color="black", zorder=10)
            ax.add_artist(circleAGB)
            if setup['single_star'] == False:
                circleComp = plt.Circle((np.hypot(xcomp, ycomp), 0.), rAccComp, transform=ax.transData._b, color="black", zorder=10)
                ax.add_artist(circleComp)

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


def allPlots(smooth, smooth_vec, zoom, rhoMin, rhoMax, vmin, vmax, Tmin, Tmax, dumpData, setup, run, loc, rAccComp, number = - 1):
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
        onePlot(fig, ax1, 'rho', rhoMin, rhoMax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp)
        onePlot(fig, ax2, 'rho', rhoMin, rhoMax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp)
        onePlot(fig, ax3, 'speed', vmin, vmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, velocity_vec = velocity_vec)
        onePlot(fig, ax4, 'speed', vmin, vmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, velocity_vec = velocity_vec)

        ax2.set_yticklabels([])
        ax4.set_yticklabels([])

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xlabel(r"$x$ [AU]", fontsize=22)
        ax4.set_xlabel(r"$x$ [AU]", fontsize=22)

    else:
        fig = plt.figure(figsize=(10 + 0.35*5, 15))

        spec = fig.add_gridspec(ncols=2, nrows=3, width_ratios=[1., 1.], height_ratios=[1, 1, 1])

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

        # the temperature colorbar limits may have to be changed...
        onePlot(fig, ax1, 'rho', rhoMin, rhoMax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp)
        onePlot(fig, ax2, 'rho', rhoMin, rhoMax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp)
        onePlot(fig, ax3, 'speed', vmin, vmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp, velocity_vec = velocity_vec)
        onePlot(fig, ax4, 'speed', vmin, vmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp, velocity_vec = velocity_vec)
        onePlot(fig, ax5, 'temp', Tmin, Tmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'z', rAccComp)
        onePlot(fig, ax6, 'temp', Tmin, Tmax, smooth, smooth_vec, zoom, dumpData, setup, axs, 'y', rAccComp)

        ax2.set_yticklabels([])
        ax4.set_yticklabels([])
        ax6.set_yticklabels([])

        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xticklabels([])
        ax4.set_xticklabels([])
        ax5.set_xlabel(r"$x$ [AU]", fontsize=22)
        ax6.set_xlabel(r"$x$ [AU]", fontsize=22)


    fig.subplots_adjust(wspace=-0.1, hspace=0.15)
    if number == -1:
        plt.savefig(os.path.join(loc, 'png/2Dplot_DensSpeedTemp_zoom{0:01d}.png'.format(zoom)), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(loc, 'pdf/2Dplot_DensSpeedTemp_zoom{0:01d}.pdf'.format(zoom)), bbox_inches="tight")
        print('         Slice plots (zoom factor = ' + str(zoom) + ') model ' + str(run) + ' ready and saved!')

    else:
        fig.text(0.5, 0.9, "Dumpfile {0:05d}".format(number), size=28)
        plt.savefig(os.path.join(loc, 'animation/2Dplot_DensSpeedTemp_zoom{0:01d}_{1:04d}.png'.format(zoom, int(number/setup['nfulldump']))), dpi=200,
                    bbox_inches="tight")
    plt.close()

'''
main definition
'''


def SlicePlots(run, loc, dumpData, setup, number = -1):
    print('')
    print('(1)  Start calculations for slice plots...')

    # Make sliceplots
    Mdot = setup['Mdot']
    bound = setup['bound']

    theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
    smooth, smooth_vec = smoothData(dumpData, setup, theta)
    rAccComp = setup['rAccrComp']
    if rAccComp <= 0.05 * radius_AGB: rAccComp = 0.05 * radius_AGB

    # Bounds
    rhoMin1 = 0.
    rhoMax1 = 0.
    rhoMin2 = 0.
    rhoMax2 = 0.
    rhoMin5 = 0.
    rhoMax5 = 0.

    vMin1 = 0.
    vMax1 = 0.
    vMin2 = 0.
    vMax2 = 0.
    vMin5 = 0.
    vMax5 = 0.

    TMin1 = 0.
    TMax1 = 0.
    TMin2 = 0.
    TMax2 = 0.
    TMin5 = 0.
    TMax5 = 0.

    # Ranges are the same as for '/gamma/non_isowind/gamma1.4/global'
    if "global" in run:
        rhoMin1, rhoMax1 = -24.39280, -16.86249
        rhoMin2, rhoMax2 = -20.80625, -17.63237
        rhoMin5, rhoMax5 = -19.71264, -16.90193

        vMin1, vMax1 = 5.96568, 38.44218
        vMin2, vMax2 = 5.99680, 33.89827
        vMin5, vMax5 = 3.17029, 25.16178

        TMin1, TMax1 = 1.26680, 4.65751
        TMin2, TMax2 = 2.58142, 4.07818
        TMin5, TMax5 = 2.52443, 4.53982

    # Ranges are the same as for '/gamma/non_isowind/gamma1.4/local'
    elif ("local" in run):
        rhoMin1, rhoMax1 = -19.96524, -16.86532
        rhoMin2, rhoMax2 = -19.57196, -15.85257
        rhoMin5, rhoMax5 = -18.47570, -14.90456

        vMin1, vMax1 = 2.54833, 21.86302
        vMin2, vMax2 = 0.03749, 18.39068
        vMin5, vMax5 = 0.00000, 22.53676

        TMin1, TMax1 = 2.26195, 4.49702
        TMin2, TMax2 = 2.14475, 4.92498
        TMin5, TMax5 = 1.97842, 5.32584


    elif customRanges:
        rhoMin1, rhoMax1 = -21, -16
        rhoMin2, rhoMax2 = -21, -16
        rhoMin5, rhoMax5 = -21, -16

        vMin1, vMax1 = 0., 30.
        vMin2, vMax2 = 0., 30.
        vMin5, vMax5 = 0., 30.

        TMin1, TMax1 = 1., 4.
        TMin2, TMax2 = 1., 4.
        TMin5, TMax5 = 1., 4.

    else:
        rhoMin1, rhoMax1 = findBounds(np.log10(smooth[1]['smooth_y']["rho"]), log=True, round=round_bounds)
        rhoMin2, rhoMax2 = findBounds(np.log10(smooth[2]['smooth_y']["rho"]), log=True, round=round_bounds)
        rhoMin5, rhoMax5 = findBounds(np.log10(smooth[5]['smooth_y']["rho"]), log=True, round=round_bounds)

        vMin1, vMax1 = findBounds(smooth[1]['smooth_y']["speed"] * cgs.cms_kms(), log=False, round=round_bounds)
        vMin2, vMax2 = findBounds(smooth[2]['smooth_y']["speed"] * cgs.cms_kms(), log=False, round=round_bounds)
        vMin5, vMax5 = findBounds(smooth[5]['smooth_y']["speed"] * cgs.cms_kms(), log=False, round=round_bounds)
        vMin1 = max(vMin1, 0.)
        vMin2 = max(vMin2, 0.)
        vMin5 = max(vMin5, 0.)

        TMin1, TMax1 = findBounds(np.log10(smooth[1]['smooth_z']["temp"]), log=True, round=round_bounds)
        TMin2, TMax2 = findBounds(np.log10(smooth[2]['smooth_z']["temp"]), log=True, round=round_bounds)
        TMin5, TMax5 = findBounds(np.log10(smooth[5]['smooth_z']["temp"]), log=True, round=round_bounds)

    if printRanges:
        print("Ranges of Parameters:\n")
        print("rhoMin1, rhoMax1 = {0:10.5f}, {1:10.5f}".format(rhoMin1, rhoMax1))
        print("rhoMin2, rhoMax2 = {0:10.5f}, {1:10.5f}".format(rhoMin2, rhoMax2))
        print("rhoMin5, rhoMax5 = {0:10.5f}, {1:10.5f}".format(rhoMin5, rhoMax5))
        print("")
        print("vMin1, vMax1 = {0:10.5f}, {1:10.5f}".format(vMin1, vMax1))
        print("vMin2, vMax2 = {0:10.5f}, {1:10.5f}".format(vMin2, vMax2))
        print("vMin5, vMax5 = {0:10.5f}, {1:10.5f}".format(vMin5, vMax5))
        print("")
        print("TMin1, TMax1 = {0:10.5f}, {1:10.5f}".format(TMin1, TMax1))
        print("TMin2, TMax2 = {0:10.5f}, {1:10.5f}".format(TMin2, TMax2))
        print("TMin5, TMax5 = {0:10.5f}, {1:10.5f}".format(TMin5, TMax5))
        print("")
        print('     Start making the slice plot figures, please wait..')
        print('')

    # Make plots
    densityPlot(smooth, 1, rhoMin1, rhoMax1, dumpData, setup, run, loc, rAccComp, number = number)
    densityPlot(smooth, 2, rhoMin2, rhoMax2, dumpData, setup, run, loc, rAccComp, number = number)
    densityPlot(smooth, 5, rhoMin5, rhoMax5, dumpData, setup, run, loc, rAccComp, number = number)

    allPlots(smooth, smooth_vec, 1, rhoMin1, rhoMax1, vMin1, vMax1, TMin1, TMax1, dumpData, setup, run, loc, rAccComp, number = number)
    allPlots(smooth, smooth_vec, 2, rhoMin2, rhoMax2, vMin2, vMax2, TMin2, TMax2, dumpData, setup, run, loc, rAccComp, number = number)
    allPlots(smooth, smooth_vec, 5, rhoMin5, rhoMax5, vMin5, vMax5, TMin5, TMax5, dumpData, setup, run, loc, rAccComp, number = number)

    
def findBounds(data, log = False, round = False):
    filtered_data = data
    result = np.zeros(2)

    if (-np.inf in data) or (np.inf in data):
        min = minInf
        max = maxInf
        if log == True:
            min = minInfLog
            max = maxInfLog

        filtered_data = filtered_data[np.logical_and(minInfLog <= data, data <= maxInfLog)]

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