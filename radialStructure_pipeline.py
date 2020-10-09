#!/usr/bin/env python
# coding: utf-8

#Import nescessary packages and the python code to load the data (as ld)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os


# import own scripts
import smoothingKernelScript as sk
import ConversionFactors_cgs as cgs

# import certain things from packages
from matplotlib    import colors
from astropy       import constants
from mpl_toolkits.axes_grid1 import AxesGrid

from matplotlib    import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size':   11})
rcParams.update({'figure.dpi': 200})
rc('font', family='serif')
rc('text', usetex=True)

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

#Definitions used:
      
def getParamsLine(results_line,rR, minT):
    indicesToKeep = np.where(np.log10(results_line['temp']) > minT)
    R = rR[indicesToKeep]
    temp  = np.log10(results_line['temp'][indicesToKeep])
    speed = (results_line['speed'][indicesToKeep]) *1e-5
    rho   = np.log10(results_line['rho'][indicesToKeep])
    return(rho,speed,temp,R)

#definition used in plotParR, plots radial structure of given parameter (log(rho)/|v|/T) on the x- and y-axis
def oneRadialStructurePlot(parX,parY,parZ, X, Y, Z, parName, axis, parMin, parMax, xcomp, xAGB, bound):

    axis.plot((X/cgs.AU_cm()),parX, '.', color = 'royalblue', label = 'x-axis', markersize = 1)
    axis.plot((Y/cgs.AU_cm()),parY, '.', color = 'firebrick', label = 'y-axis', markersize = 1)
    axis.plot((Z/cgs.AU_cm()),parZ, '.', color = 'goldenrod', label = 'z-axis', markersize = 1)

    #CHANGE TO OUTER BOUNDARY
    axis.set_xlim(-bound,bound)
    axis.vlines(xcomp, parMin, parMax ,'k', linestyle = 'dashed' ,linewidth = 1)
    axis.vlines(xAGB, parMin, parMax,'k', linestyle = 'solid', linewidth = 1)

    axis.set_ylabel(parName, fontsize = 27)
    axis.tick_params(labelsize=20)


def radialStructPlots(run,directory, dumpData, setup):
    print('')
    print('(2)  Start calculations for the radial structure plots.')
    try:
        os.mkdir(directory+'RadialStructurePlots/')
    except OSError:
        print('')


    # Define legend
    xAxis  = mlines.Line2D([],[], color = 'royalblue', label='x-axis')
    yAxis  = mlines.Line2D([],[], color = 'firebrick', label='y-axis')
    zAxis  = mlines.Line2D([],[], color = 'goldenrod', label='z-axis')
    comp   = mlines.Line2D([],[], color = 'k', linestyle = 'dashed',linewidth = 1, label = 'comp')
    AGB    = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 1, label='AGB')
    AGBs   = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 2, label='AGB')
    handles1 = [xAxis,yAxis,zAxis,comp,AGB]
    handles2 = [xAxis,yAxis,zAxis,AGBs]


    #plots radial structure of log(rho), |v| and T on the x- and y-axis
        
    if setup['single_star'] == True:
        xcomp = 0
        xAGB  = 0
        handl = handles2
        
    else:       
        xcomp = dumpData['posComp'][0]/cgs.AU_cm()
        xAGB  = dumpData['posAGB'][0]/cgs.AU_cm()
        handl = handles1
    
    fig, (ax1,ax2,ax3)= plt.subplots(3, 1,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [1]})
    fig.set_size_inches(8, 18)
    
    #Remove data of AGB and companion
    dataToUse = {}
    dataToUse['rho']      = dumpData['rho'][:-2]       
    dataToUse['temp']     = dumpData['temp'][:-2]      
    dataToUse['speed']    = dumpData['speed'][:-2]     
    dataToUse['mass']     = dumpData['mass'][:-2]
    dataToUse['position'] = dumpData['position'][:-2]
    dataToUse['h']        = dumpData['h'][:-2]
    
    #calculate smoothed data around one axis
    results_line_X,xX,yX,zX = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_x',setup['bound']*cgs.AU_cm())
    results_line_Y,xY,yY,zY = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_y',setup['bound']*cgs.AU_cm())
    results_line_Z,xZ,yZ,zZ = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_z',setup['bound']*cgs.AU_cm())

    
    parX = getParamsLine(results_line_X, xX,1.8)
    parY = getParamsLine(results_line_Y, yY,1.8)
    parZ = getParamsLine(results_line_Z, zZ,1.8)
    X = parX[3]
    Y = parY[3]
    Z = parZ[3]
    
    
    Mdot  = setup['Mdot']
    vini  = setup['v_ini'] 
    
    #Select limits 
    if Mdot > 1e-5:
        rhoMin = -21.5
        rhoMax = -10.5
    elif Mdot < 5e-7:
        rhoMin = -24
        rhoMax = -13
    elif 5e-7 <= Mdot <= 1e-5:
        rhoMin = -23
        rhoMax = -12

    vmin = -3
    vmax = 60
    Tmin = 1.7
    Tmax = 5.7

    
    oneRadialStructurePlot(parX[0],parY[0], parZ[0], X, Y, Z, 'log($\\rho$[g/cm$^3$])', ax1, rhoMin, rhoMax, xcomp, xAGB, setup['bound'])
    ax1.set_title('v = '+ str(vini)+ 'km/s', fontsize = 33)#, Mdot ='+ str(Mdot)+ '$M_\odot$/yr, ecc = ' +str(ecc))
    oneRadialStructurePlot(parX[1],parY[1],parZ[1], X, Y, Z, '$|v|$[km/s]', ax2, vmin, vmax,  xcomp, xAGB,setup['bound'])
    # ax2.set_ylim(0,vmax+5)
    ax2.legend(handles = handl, loc = 'upper left', fontsize = 20) #only for M18 and M16
    oneRadialStructurePlot(parX[2],parY[2], parZ[2], X, Y, Z, 'log($T$[K])', ax3, Tmin, Tmax, xcomp, xAGB, setup['bound'])    
    ax3.set_xlabel('r[AU]', fontsize = 27)
                   
    fig.savefig(directory+'RadialStructurePlots/1D_'+str(run)+'.png')
    print('     Radial structure plot model '+str(run)+' ready and saved!')



