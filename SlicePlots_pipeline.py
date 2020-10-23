#Import packages
import numpy as np
import matplotlib.pyplot as plt
import os

# import own scripts
import smoothingKernelScript as sk
import ConversionFactors_cgs as cgs

# import certain things from packages
from matplotlib                 import colors
from astropy                    import constants
from mpl_toolkits.axes_grid1    import AxesGrid

from matplotlib                 import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size':   11})
rcParams.update({'figure.dpi': 200})
rc('font', family='serif')
rc('text', usetex=True)

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

#Definitions needed to make sliceplots:

#Load the smoothing kernel data
def smoothData(dumpData,setup):
    #Uncomment the ones you need, takes long to run, so don't run if not nescessary
    print('     Calculating the smoothing kernels. This may take a while, please wait...')
    #zoom = 1
    results_sph_sl_z,x1,y1,z1 = sk.getSmoothingKernelledPix(400,20,dumpData,['rho','temp','speed'], 'comp','z',setup['bound']*cgs.AU_cm())
    results_sph_sl_y,x2,y2,z2 = sk.getSmoothingKernelledPix(400,20,dumpData,['rho','temp','speed'], 'comp','y',setup['bound']*cgs.AU_cm())
    #zoom =2
    results_sph_sl_zZ2,x1Z2,y1Z2,z1Z2 = sk.getSmoothingKernelledPix(400,20,dumpData,['rho','temp','speed'], 'comp','z',setup['bound']*cgs.AU_cm()/2)
    results_sph_sl_yZ2,x2Z2,y2Z2,z2Z2 = sk.getSmoothingKernelledPix(400,20,dumpData,['rho','temp','speed'], 'comp','y',setup['bound']*cgs.AU_cm()/2)
    #zoom = 5
    results_sph_sl_zZ5,x1Z5,y1Z5,z1Z5 = sk.getSmoothingKernelledPix(500,20,dumpData,['rho','temp','speed'], 'comp','z',setup['bound']*cgs.AU_cm()/5)
    results_sph_sl_yZ5,x2Z5,y2Z5,z2Z5 = sk.getSmoothingKernelledPix(500,20,dumpData,['rho','temp','speed'], 'comp','y',setup['bound']*cgs.AU_cm()/5)

    smooth = {'sph_sl_z'     : results_sph_sl_z,
             'sph_sl_zZ2'    : results_sph_sl_zZ2,
             'sph_sl_zZ5'    : results_sph_sl_zZ5,
             'sph_sl_y'      : results_sph_sl_y,
             'sph_sl_yZ2'    : results_sph_sl_yZ2,
             'sph_sl_yZ5'    : results_sph_sl_yZ5,
             'xz'            : x1 ,
             'yz'            : y1 ,
             'xy'            : x2 ,
             'zy'            : z2 ,
        }

    return smooth

# Makes one of the plots that will be combined in one figure
def onePlot(ax,par,name,colormap,mi,ma,x,y,xlabel,ylabel,lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6):

    xAGB = dumpData['posAGB'][0]/cgs.AU_cm()
    yAGB = dumpData['posAGB'][1]/cgs.AU_cm()
    zAGB = dumpData['posAGB'][2]/cgs.AU_cm()

    if setup['single_star']== False:
        xcomp = dumpData['posComp'][0]/cgs.AU_cm()
        ycomp = dumpData['posComp'][1]/cgs.AU_cm()
        zcomp = dumpData['posComp'][2]/cgs.AU_cm()

 
    ax.axis('equal')
    ax.set_facecolor('k')
    axPlot = ax.scatter(x/cgs.AU_cm(), y/cgs.AU_cm(), s=11, c=par, cmap=colormap, vmin=mi, vmax = ma)
    if ax == ax2 or ax == ax4 or ax == ax6:
        ax.plot(xAGB,zAGB, 'bo', markersize = 4,label = 'AGB')
        if setup['single_star']== False:
            ax.plot(xcomp,zcomp, 'ro', markersize = 3,label = 'comp')  
        cbar = plt.colorbar(axPlot, ax = ax)
        cbar.set_label(name, fontsize = 18)
    if ax == ax5 or ax == ax6:
        ax.set_xlabel(xlabel)
    if ax == ax1 or ax == ax3 or ax == ax5:
        ax.plot(xAGB,yAGB, 'bo', markersize = 4,label = 'AGB')
        if setup['single_star']== False:
            ax.plot(xcomp,ycomp, 'ro', markersize = 3,label = 'comp')  
    ax.axis([-lim,lim,-lim,lim])
    ax.set_ylabel(ylabel)


# Make figure with the x-y(left) and x-z(right) slice plots of log(rho[g/cm3]), log(T[K]) and |v|[km/s]. Takes 'zoom'factor as input
def allPlots(modelname, smooth, zoom, rhoMin, rhoMax, vmax, bound, dumpData, setup, run, loc):

    cm_rho  = plt.cm.get_cmap('inferno')
    cm_T    = plt.cm.get_cmap('afmhot')
    cm_v    = plt.cm.get_cmap('viridis')
#     fig = plt.figure(figsize=(11, 14))
    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6))= plt.subplots(3, 2,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [0.8,1]})
    fig.set_size_inches(10, 14)
    lim = bound/zoom  
          
    # the temperature colorbar limits may have to be changed...
    onePlot(ax1,np.log10(smooth['sph_sl_z']['rho']),'log($\\rho[$g/cm$^3]$)', cm_rho, rhoMin, rhoMax, smooth['xz'],smooth['yz'],'x[AU]', 'y[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    onePlot(ax2,np.log10(smooth['sph_sl_y']['rho']),'log($\\rho[$g/cm$^3]$)', cm_rho, rhoMin, rhoMax, smooth['xy'],smooth['zy'],'x[AU]', 'z[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    onePlot(ax3,np.log10(smooth['sph_sl_z']['temp']),'log$(T$[K])', cm_T, 1.7, 5.2,smooth['xz'],smooth['yz'],'x[AU]', 'y[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    onePlot(ax4,np.log10(smooth['sph_sl_y']['temp']),'log$(T$[K])', cm_T, 1.7, 5.2,smooth['xy'],smooth['zy'], 'x[AU]', 'z[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    onePlot(ax5,smooth['sph_sl_z']['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, vmax, smooth['xz'],smooth['yz'], 'x[AU]', 'y[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    onePlot(ax6,smooth['sph_sl_y']['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, vmax, smooth['xy'],smooth['zy'], 'x[AU]', 'z[AU]',lim,dumpData,setup,ax1,ax2,ax3,ax4,ax5,ax6)
    fig.savefig(loc+'2DSliceplots/'+str(run)+'Z'+str(zoom)+'.png',dpi = 300)
    print('         Slice plots (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')

# Make figure with the x-y(orbital plane) slice plot of log(rho[g/cm3]). Takes 'zoom'factor as input
def densityPlot(modelname,smooth, zoom,rhoMin,rhoMax,vmax,bound,dumpData, setup,run, loc):

    cm_rho  = plt.cm.get_cmap('inferno')
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(9.8, 8)
        
    ax.axis('equal')
    ax.set_facecolor('k')
    axPlot = ax.scatter(smooth['xz']/cgs.AU_cm(),smooth['yz']/cgs.AU_cm(),s=11,c=np.log10(smooth['sph_sl_z']['rho']),cmap=cm_rho,vmin=rhoMin, vmax = rhoMax)
    if setup['single_star']== False:
        xAGB = dumpData['posAGB'][0]/cgs.AU_cm()
        yAGB = dumpData['posAGB'][1]/cgs.AU_cm()
        xcomp = dumpData['posComp'][0]/cgs.AU_cm()
        ycomp = dumpData['posComp'][1]/cgs.AU_cm()
        ax.plot(xAGB,yAGB, 'bo', markersize = 5,label = 'AGB')
        ax.plot(xcomp,ycomp, 'ro', markersize = 5,label = 'comp')   

    cbar = plt.colorbar(axPlot, ax = ax)
    cbar.set_label('log($\\rho$[g/cm$^3]$)', fontsize = 25)
    cbar.ax.tick_params(labelsize=14)
    ax.set_xlabel('x[AU]', fontsize = 18)
    ax.set_ylabel('y[AU]', fontsize = 18)
    lim = bound/zoom  
    ax.axis([-lim,lim,-lim,lim])
    ax.tick_params(labelsize=14)

    fig.savefig(loc+'2DSliceplots/1Plot_'+str(run)+'Z'+str(zoom)+'.png',dpi = 300)
    print('         Density slice plot (zoom factor = '+str(zoom)+') model '+str(run)+' ready and saved!')


def SlicePlots(run,loc, dumpData, setup):
    print('')
    print('(1)  Start calculations for slice plots...')
     # Make new directory to save figures in, unless directory already exists
    try:
        os.mkdir(loc+'2DSliceplots/')
    except OSError:
        print('')

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
    densityPlot('M'+str(run),smooth,1,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)
    densityPlot('M'+str(run),smooth,2,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)
    densityPlot('M'+str(run),smooth,5,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)
    allPlots('M'+str(run),smooth,1,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)
    allPlots('M'+str(run),smooth,2,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)
    allPlots('M'+str(run),smooth,5,rhoMin,rhoMax,vmax,bound,dumpData,setup,run, loc)

   
 # # Uncomment and change if we add Hill sphere plots

 # Make figure with the x-y(left) and x-z(right) Hill sphere slice plots of log(rho[g/cm3]) and |v|[km/s].

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
