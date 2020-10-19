#Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import os

# import own scripts
import smoothingKernelScript    as sk
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
    #Uncomment the ones you need, takes long to run, so don't run if not nescessary
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
Make figure with the x-y(orbital plane) slice plot of log(rho[g/cm3]). Takes 'zoom'factor as input
'''
def densityPlot(smooth,zoom,rhoMin,rhoMax,vmax,dumpData, setup,run, loc):

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
        ax.plot(xAGB ,yAGB , 'ko', markersize = 4  ,label = 'AGB')
        ax.plot(xcomp,ycomp, 'ro', markersize = 2,label = 'comp')   

    cbar = plt.colorbar(axPlot, ax = ax)
    cbar.set_label('log density [cm/g$^3$]', fontsize = 25)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xlabel('x [AU]', fontsize = 18)
    ax.set_ylabel('y [AU]', fontsize = 18)
      
    lim = (setup['bound']-30)/zoom
    ax.axis([-lim,lim,-lim,lim])
    ax.tick_params(labelsize=14)

    fig.tight_layout()
    fig.savefig(loc+'1Plot_'+str(run)+'_zoom'+str(zoom)+'.png',dpi = 300)
    
    print('Density slice plot (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')

'''
Makes one of the plots that will be combined in one figure
'''
def onePlot(ax, par, mi, ma, smooth, zoom, dumpData, setup, axs, plane):
    
    # colormap per parameter
    cm = {'rho'  : plt.cm.get_cmap('viridis' ),
          'temp' : plt.cm.get_cmap('afmhot'  ),
          'speed': plt.cm.get_cmap('PuBuGn_r')
          }
    # label name per parameter
    name = {'rho'  : 'log density [cm/g$^3$]',
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
Make figure with the x-y(left) and x-z(right) slice plots of log(rho[g/cm3]), log(T[K]) and |v|[km/s]. Takes 'zoom'factor as input
'''
def allPlots(smooth, zoom, rhoMin, rhoMax, vmax,  bound, dumpData, setup, run, loc):

    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6))= plt.subplots(3, 2,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [0.81,1]})
    fig.set_size_inches(12, 14.4)
    
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
    
    fig.savefig(loc +str(run)+'_zoom'+str(zoom)+'.png',dpi = 300)
    
    print('         Slice plots (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')



def SlicePlots(run,loc, dumpData, setup):
    print('')
    print('(1)  Start calculations for slice plots...')

    # Make sliceplots
    Mdot  = setup['Mdot']
    bound = setup['bound']   

    #Set the limits of the colorbars
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
    densityPlot( smooth, 1, rhoMin, rhoMax, bound, dumpData, setup, run, loc)
    densityPlot( smooth, 2, rhoMin, rhoMax, bound, dumpData, setup, run, loc)
    densityPlot( smooth, 5, rhoMin, rhoMax, bound, dumpData, setup, run, loc)
    allPlots(    smooth, 1, rhoMin, rhoMax, vmax, bound, dumpData, setup, run, loc)
    allPlots(    smooth, 2, rhoMin, rhoMax, vmax, bound, dumpData, setup, run, loc)
    allPlots(    smooth, 5, rhoMin, rhoMax, vmax, bound, dumpData, setup, run, loc)

'''   
NOTE: uncomment and change if we add Hill sphere plots
'''

## Make figure with the x-y(left) and x-z(right) Hill sphere slice plots of log(rho[g/cm3]) and |v|[km/s].

# def HillsphPlot(modelname,results_sph_sl_y,x2,z2,zoom,rhoMin,rhoMax):
#     cm_rho  = plt.cm.get_cmap('inferno')
#     cm_T    = plt.cm.get_cmap('afmhot')
#     cm_v    = plt.cm.get_cmap('viridis')
#     cm_vtvv = plt.cm.get_cmap('seismic')
# #     fig = plt.figure(figsize=(11, 14))
#     fig, ((ax1,ax2),(ax3,ax4))= plt.subplots(2, 2,  gridspec_kw={'height_ratios':[1,1],'width_ratios': [1,1]})
#     fig.set_size_inches(10, 8)
    
#     #CHANGE 200 TO OUTER BOUNDARY
#     lim = setup['bound']/zoom        

#     def onePlotRH(axis,par,name,colormap,mi,ma,x,y,xlabel,ylabel):
#         axi = axis
#         axi.axis('equal')
#         axi.set_facecolor('k')
#         axPlot = axi.scatter(x/cgs.AU_cm(),y/cgs.AU_cm(),s=11,c=par,cmap=colormap,vmin=mi, vmax = ma)
#         cbar = plt.colorbar(axPlot, ax = axi)
#         cbar.set_label(name, fontsize = 18)
#         if axi == ax3 or axi == ax4:
#             axi.set_xlabel(xlabel)
#         if axi == ax1 or axi == ax3:
#             axi.set_ylabel(ylabel)
#         axi.axis([-lim,lim,-lim,lim])
#         axi.plot(0,0, 'ro', markersize = 5,label = 'comp')
 
#     onePlotRH(ax1,np.log10(results_sph_sl_y['rho']),'log($\\rho$[g/cm$^3$])', cm_rho, rhoMin, rhoMax, x2,z2,'x[AU]', 'z[AU]')
#     onePlotRH(ax2,np.log10(results_sph_sl_y['temp']),'log($T$[K])', cm_T, 3, 6,x2,z2, 'x[AU]', 'z[AU]')
#     onePlotRH(ax3,results_sph_sl_y['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, 26, x2,z2, 'x[AU]', 'z[AU]')
#     onePlotRH(ax4,results_sph_sl_y['vtvv'],'$(v_t/v)^2$', cm_vtvv, 0, 1, x2, z2, 'x[AU]', 'z[AU]' )
#     zoom = '_RH'
#     fig.savefig(loc+'2DSliceplots/HillSphere_'+str(run)+'Z'+str(zoom)+'.png',dpi = 300)
#     print('Hill sphere plot model '+str(run)+' ready and saved!')


    # #Load data for Hill sphere
    # xcomp = dumpData['posComp'][0]/cgs.AU_cm()
    # ycomp = dumpData['posComp'][1]/cgs.AU_cm()
    # zcomp = dumpData['posComp'][2]/cgs.AU_cm()
    # hillsph = dumpData['rHill']/cgs.AU_cm()
    # # bound = dumpData['boundary']/cgs.AU_cm()
    # # zoom = int(bound/hillsph)
    # # dataHS = ld.LoadData_cgs_inner(str(run),hillsph,xcomp*AU,ycomp*AU,zcomp*AU)
    # #Make Hill sphere plots 
    # dumpDataHS  = dataHS #rcomp is now (0,0,0)
    # dumpData    = data
    # xAGB = dumpData['posAGB'][0]/cgs.AU_cm()
    # yAGB = dumpData['posAGB'][1]/cgs.AU_cm()
    # zAGB = dumpData['posAGB'][2]/cgs.AU_cm()
    # xcomp = dumpData['posComp'][0]/cgs.AU_cm()
    # ycomp = dumpData['posComp'][1]/cgs.AU_cm()
    # zcomp = dumpData['posComp'][2]/cgs.AU_cm()
    # Mdot  = setup['Mdot']
    # v_ini = setup['v_ini'] *cgs.cms_kms()
    # # hillsph = dumpData['rHill']/cgs.AU_cm()
    # # bound = setup['boundary']/cgs.AU_cm()
    # # zoom = int(bound/hillsph)

    # results_sph_sl_yHS,x2HS,y2HS,z2HS = sk.getSmoothingKernelledPix(400,20,dumpDataHS,['rho','temp','speed','vtvv'], 'comp','y',10*AU)

    # #Change these limits if the plots don't look nice
    # if Mdot > 1e-5:
    #     if v_ini == 5:
    #         rhoMin = -14.5
    #         rhoMax = -11.5
    #     elif v_ini == 20:
    #         rhoMin = -15.5
    #         rhoMax = -12.5

    # elif Mdot < 5e-7:
    #     if v_ini == 5:
    #         rhoMin = -17
    #         rhoMax = -14
    #     elif v_ini == 20:
    #         rhoMin = -18
    #         rhoMax = -15

    # elif 5e-7 <= Mdot <= 1e-5:
    #     rhoMin = -16
    #     rhoMax = -13

    # HillsphPlot('M'+str(run),results_sph_sl_yHS,x2HS,z2HS,zoom,rhoMin,rhoMax)
