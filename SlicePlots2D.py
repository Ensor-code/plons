#Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import os

# import own scripts
import SmoothingKernelScript    as sk
import ConversionFactors_cgs    as cgs

# import certain things from packages
from matplotlib                 import colors
from astropy                    import constants
from mpl_toolkits.axes_grid1    import AxesGrid

from matplotlib                 import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size':   11})
rcParams.update({'figure.dpi': 200})
#rc('font', family='serif')
#rc('text', usetex=True)

# ignore warnings
import warnings
warnings.filterwarnings("ignore")


'''
Load the smoothing kernel data
'''
def smoothData(dumpData,setup):
    #comment the ones you don't need, it takes long to run, so don't run if not nescessary
    print('     Calculating the smoothing kernels. This may take a while, please wait...')
    #zoom = 1
    results_sph_sl_z  ,x1  ,y1   ,z1  = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','z',(setup['bound'])*cgs.AU_cm() )
    results_sph_sl_y  ,x2  ,y2   ,z2  = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','y',(setup['bound'])*cgs.AU_cm() )
    #zoom = 2
    results_sph_sl_zZ2,x1Z2,y1Z2,z1Z2 = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','z',(setup['bound'])*cgs.AU_cm()/2)
    results_sph_sl_yZ2,x2Z2,y2Z2,z2Z2 = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','y',(setup['bound'])*cgs.AU_cm()/2)
    #zoom = 5
    results_sph_sl_zZ5,x1Z5,y1Z5,z1Z5 = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','z',(setup['bound'])*cgs.AU_cm()/5)
    results_sph_sl_yZ5,x2Z5,y2Z5,z2Z5 = sk.getSmoothingKernelledPix(300,20,dumpData,['rho','temp','speed'], 'comp','y',(setup['bound'])*cgs.AU_cm()/5)
   
    smooth_zoom1 = {'smooth_z'     :  results_sph_sl_z,
                    'x_z'          :  x1,
                    'y_z'          :  y1,
                    'smooth_y'     :  results_sph_sl_y,
                    'x_y'          :  x2,
                    'z_y'          :  z2
        }
    
    smooth_zoom2 = {'smooth_z'     :  results_sph_sl_zZ2,
                    'x_z'          :  x1Z2,
                    'y_z'          :  y1Z2, 
                    'smooth_y'     :  results_sph_sl_yZ2,
                    'x_y'          :  x2Z2,
                    'z_y'          :  z2Z2
        }
    
    smooth_zoom5 = {'smooth_z'     :  results_sph_sl_zZ5,
                    'x_z'          :  x1Z5,
                    'y_z'          :  y1Z5, 
                    'smooth_y'     :  results_sph_sl_yZ5,
                    'x_y'          :  x2Z5,
                    'z_y'          :  z2Z5
        }
    
    smooth = { 1: smooth_zoom1,
               2: smooth_zoom2,
               5: smooth_zoom5
        }

    return smooth


'''
Make figure with the xy(orbital plane) slice plot of log density [g/cm^3]. 
    - smooth            is smoothing kernel data
    - zoom              zoom factor of the plot
    - rhoMin/rhoMax     determine the min and max of the colorbar
    - dumpData      data of dump file
    - setup         setup data
    - run           run number [str]
    - loc           output directory
'''
def densityPlot(smooth, zoom, rhoMin, rhoMax, dumpData, setup, run, loc):

    cm_rho  = plt.cm.get_cmap('viridis')
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(9.8, 8)
        
    ax.axis('equal')
    ax.set_facecolor('k')
    axPlot = ax.scatter(smooth[zoom]['x_z']/cgs.AU_cm(),smooth[zoom]['y_z']/cgs.AU_cm(),s=5,c=np.log10(smooth[zoom]['smooth_z']['rho']),cmap=cm_rho,vmin=rhoMin, vmax = rhoMax)

    if setup['single_star']== False:
        xAGB  = dumpData['posAGB' ][0] / cgs.AU_cm()
        yAGB  = dumpData['posAGB' ][1] / cgs.AU_cm()
        xcomp = dumpData['posComp'][0] / cgs.AU_cm()
        ycomp = dumpData['posComp'][1] / cgs.AU_cm()
        ax.plot(xAGB ,yAGB , 'ko', markersize =   4*zoom ,label = 'AGB')
        ax.plot(xcomp,ycomp, 'ro', markersize = 1.5*zoom ,label = 'comp')   

    cbar = plt.colorbar(axPlot, ax = ax)
    cbar.set_label('log density [g/cm$^3$]', fontsize = 25)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xlabel('x [AU]', fontsize = 18)
    ax.set_ylabel('y [AU]', fontsize = 18)
      
    lim = (setup['bound']*np.sqrt(2)/2)/zoom
    ax.axis([-lim,lim,-lim,lim])
    ax.tick_params(labelsize=14)

    fig.tight_layout()
    fig.savefig(os.path.join(loc, 'png/2Dplot_density_zoom'+str(zoom)) + ".png",dpi = 300)
    fig.savefig(os.path.join(loc, 'pdf/2Dplot_density_zoom'+str(zoom)) + ".pdf")

    print('         Density slice plot (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')


'''
Makes one slice plot 

INPUT
    - ax        is given subplot
    - par       is the name of the parameter
    - mi/ma     colorbar limits
    - smooth    smoothed data
    - zoom      zoom factor for the plot
'''
def onePlot(ax, par, mi, ma, smooth, zoom, dumpData, setup, axs, plane):
    
    # colormap per parameter
    cm = {'rho'  : plt.cm.get_cmap('viridis' ),
          'temp' : plt.cm.get_cmap('afmhot'  ),
          'speed': plt.cm.get_cmap('PuBuGn_r')
          }
    # label name per parameter
    name = {'rho'  : 'log density [g/cm$^3$]',
            'temp' : 'log temperature [K]' ,
            'speed': 'speed [km/s]'
          }
            
    xAGB = dumpData['posAGB'][0]/cgs.AU_cm()
    yAGB = dumpData['posAGB'][1]/cgs.AU_cm()
    zAGB = dumpData['posAGB'][2]/cgs.AU_cm()

    if setup['single_star']== False:
        xcomp = dumpData['posComp'][0]/cgs.AU_cm()
        ycomp = dumpData['posComp'][1]/cgs.AU_cm()
        zcomp = dumpData['posComp'][2]/cgs.AU_cm()

 
    ax.axis('equal')
    ax.set_facecolor('k')
    if not par == 'speed':
        if plane == 'z':
            axPlot = ax.scatter(smooth[zoom]['x_z']/cgs.AU_cm(),smooth[zoom]['y_z']/cgs.AU_cm(),s=5,c=np.log10(smooth[zoom]['smooth_z'][par]),cmap=cm[par],vmin=mi, vmax = ma)
            ax.set_ylabel('y [AU]')
        if plane == 'y':
            axPlot = ax.scatter(smooth[zoom]['x_y']/cgs.AU_cm(),smooth[zoom]['z_y']/cgs.AU_cm(),s=5,c=np.log10(smooth[zoom]['smooth_y'][par]),cmap=cm[par],vmin=mi, vmax = ma)
            ax.set_ylabel('z [AU]')
    
    if par == 'speed':
        if plane == 'z':
            axPlot = ax.scatter(smooth[zoom]['x_z']/cgs.AU_cm(),smooth[zoom]['y_z']/cgs.AU_cm(),s=5,c=(smooth[zoom]['smooth_z'][par]*cgs.cms_kms()),cmap=cm[par],vmin=mi, vmax = ma)
            ax.set_ylabel('y [AU]')
        if plane == 'y':
            axPlot = ax.scatter(smooth[zoom]['x_y']/cgs.AU_cm(),smooth[zoom]['z_y']/cgs.AU_cm(),s=5,c=(smooth[zoom]['smooth_y'][par]*cgs.cms_kms()),cmap=cm[par],vmin=mi, vmax = ma)
            ax.set_ylabel('z [AU]')
    
    
    
    # plot the position of the AGB star and comp in the edge-on plane & make colorbar 
    if ax == axs[2] or ax == axs[4] or ax == axs[6]:
        ax.plot(xAGB, zAGB, 'ko', markersize = 4*zoom,label = 'AGB')
        if setup['single_star'] == False:
            ax.plot(xcomp, zcomp, 'ro', markersize = 1.5*zoom,label = 'comp')  
        cbar = plt.colorbar(axPlot, ax = ax)
        cbar.set_label(name[par], fontsize = 18)
        
    if ax == axs[5] or ax == axs[6]:
        ax.set_xlabel('x [AU]')
    
    # plot the position of the AGB star and comp in the face-on plane
    if ax == axs[1] or ax == axs[3] or ax == axs[5]:
        ax.plot(xAGB,yAGB, 'ko', markersize = 4*zoom,label = 'AGB')
        if setup['single_star'] == False:
            ax.plot(xcomp,ycomp, 'ro', markersize = 1.5*zoom,label = 'comp')  
            
    lim = (setup['bound']-30)/zoom
    ax.axis([-lim,lim,-lim,lim])


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
def allPlots(smooth, zoom, rhoMin, rhoMax, vmax, dumpData, setup, run, loc):

    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6))= plt.subplots(3, 2,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [0.81,1]})
    fig.set_size_inches(12, 14.3)
    
    axs = {1: ax1,
           2: ax2,
           3: ax3,
           4: ax4,
           5: ax5,
           6: ax6
        }
          
    # the temperature colorbar limits may have to be changed...
    onePlot(ax1,'rho'  , rhoMin, rhoMax, smooth, zoom, dumpData, setup, axs, 'z')
    onePlot(ax2,'rho'  , rhoMin, rhoMax, smooth, zoom, dumpData, setup, axs, 'y')
    onePlot(ax5,'temp' , 1.7   , 5.2   , smooth, zoom, dumpData, setup, axs, 'z')
    onePlot(ax6,'temp' , 1.7   , 5.2   , smooth, zoom, dumpData, setup, axs, 'y')
    onePlot(ax3,'speed', 0     , vmax  , smooth, zoom, dumpData, setup, axs, 'z')
    onePlot(ax4,'speed', 0     , vmax  , smooth, zoom, dumpData, setup, axs, 'y')
    
    ax1.axes.get_xaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax4.axes.get_xaxis().set_visible(False)
    
    fig.tight_layout()
    fig.subplots_adjust(wspace = 0.2,hspace = 0.005)

    fig.savefig(os.path.join(loc, 'png/2Dplot_DensSpeedTemp_zoom'+str(zoom)) + '.png',dpi = 300)
    fig.savefig(os.path.join(loc, 'pdf/2Dplot_DensSpeedTemp_zoom'+str(zoom)) + '.pdf')
    print('         Slice plots (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')


'''
main definition
'''
def SlicePlots(run,loc, dumpData, setup):
    print('')
    print('(1)  Start calculations for slice plots...')

    # Make sliceplots
    Mdot  = setup['Mdot']
    bound = setup['bound']
    rhoMin = 0
    rhoMax = 0

    #Set the limits of the density plot colorbars, scales exactly with Mdot if no cooling is included
    if Mdot > 1e-5:
        rhoMin = -18.5
        rhoMax = -13.5
    elif Mdot < 5e-7:
        rhoMin = -21
        rhoMax = -16
    elif 5e-7 <= Mdot <= 1e-5:
        rhoMin = -20
        rhoMax = -15
    vmax = 26

    smooth = smoothData(dumpData,setup)

    print('     Start making the slice plot figures, please wait..')
    print('')
    # Make plots
    densityPlot( smooth, 1, rhoMin, rhoMax, dumpData, setup, run, loc)
    densityPlot( smooth, 2, rhoMin, rhoMax, dumpData, setup, run, loc)
    densityPlot( smooth, 5, rhoMin, rhoMax, dumpData, setup, run, loc)
    allPlots(    smooth, 1, rhoMin, rhoMax, vmax, dumpData, setup, run, loc)
    allPlots(    smooth, 2, rhoMin, rhoMax, vmax, dumpData, setup, run, loc)
    allPlots(    smooth, 5, rhoMin, rhoMax, vmax, dumpData, setup, run, loc)
