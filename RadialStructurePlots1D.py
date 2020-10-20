import numpy                    as np
import matplotlib.pyplot        as plt
import matplotlib.lines         as mlines
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

    
def getParamsLine(results_line,rR, minT):
    indicesToKeep = np.where(np.log10(results_line['temp']) > minT)
    R     = rR[indicesToKeep]
    temp  = (results_line['temp' ][indicesToKeep])
    speed = (results_line['speed'][indicesToKeep]) * cgs.cms_kms()
    rho   = (results_line['rho'  ][indicesToKeep])
    
    return(rho,speed,temp,R)

'''
definition used in plotParR, plots radial structure of given parameter (log(rho)/|v|/T) on the x- and y-axis
'''
def oneRadialStructurePlot(parX,parY,parZ, X, Y, Z, parName, color, axis, parMin, parMax, xcomp, xAGB, bound):
    
    
    axis.plot((X/cgs.AU_cm()),parX, ls = '-' , color = color, label = 'x-axis', markersize = 0.3, lw = 0.8)
    axis.plot((Y/cgs.AU_cm()),parY, ls = ':' , color = color, label = 'y-axis', markersize = 0.3, lw = 0.8)
    axis.plot((Z/cgs.AU_cm()),parZ, ls = '-.', color = color, label = 'z-axis', markersize = 0.3, lw = 0.8)

    axis.set_xlim(-bound,bound)
    #axis.vlines(xcomp, parMin, parMax ,'k', linestyle = 'dashed' ,linewidth = 0.5)
    #axis.vlines(xAGB , parMin, parMax ,'k', linestyle = 'solid', linewidth = 0.5)
    axis.set_ylim(parMin,parMax)


    axis.set_ylabel(parName, fontsize = 13)
    axis.tick_params(labelsize=10)


def radialStructPlots(run,loc, dumpData, setup):
    print('')
    print('(2)  Start calculations for the radial structure plots.')

    # Define legend
    xAxis  = mlines.Line2D([],[], ls = '-' , color = 'k', label='x-axis', markersize = 0.3, lw = 0.8)
    yAxis  = mlines.Line2D([],[], ls = ':' , color = 'k', label='y-axis', markersize = 0.3, lw = 0.8)
    zAxis  = mlines.Line2D([],[], ls = '-.', color = 'k', label='z-axis', markersize = 0.3, lw = 0.8)
    #comp   = mlines.Line2D([],[], color = 'royalblue', linestyle = 'dashed',linewidth = 1, label = 'comp')
    #AGB    = mlines.Line2D([],[], color = 'royalblue', linestyle = 'solid', linewidth = 1, label = 'AGB' )
    handles1 = [xAxis,yAxis,zAxis]#,AGB, comp]
    handles2 = [xAxis,yAxis,zAxis]#,AGB]


    #plots radial structure of log(rho), |v| and T on the x- and y-axis
        
    if setup['single_star'] == True:
        xcomp = 0
        xAGB  = 0
        handl = handles2
        
    else:       
        xcomp = dumpData['posComp'][0]/cgs.AU_cm()
        xAGB  = dumpData['posAGB' ][0]/cgs.AU_cm()
        handl = handles1
    
    fig = plt.figure(figsize=(5, 9))
    
    #fig = plt.figure(figsize=(4.5, 10))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    
    
    #Remove data of AGB and companion
    dataToUse = {}
    dataToUse['rho'     ] = dumpData['rho'     ][:-2]       
    dataToUse['temp'    ] = dumpData['temp'    ][:-2]      
    dataToUse['speed'   ] = dumpData['speed'   ][:-2]     
    dataToUse['mass'    ] = dumpData['mass'    ][:-2]
    dataToUse['position'] = dumpData['position'][:-2]
    dataToUse['h'       ] = dumpData['h'       ][:-2]
    
    #calculate smoothed data around one axis
    results_line_X,xX,yX,zX = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_x',setup['bound']*cgs.AU_cm())
    results_line_Y,xY,yY,zY = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_y',setup['bound']*cgs.AU_cm())
    results_line_Z,xZ,yZ,zZ = sk.getSmoothingKernelledPix(10000,20,dataToUse,['rho','temp','speed'], 'comp','line_z',setup['bound']*cgs.AU_cm())

    
    parX = getParamsLine(results_line_X, xX,2.3)
    parY = getParamsLine(results_line_Y, yY,1.8)
    parZ = getParamsLine(results_line_Z, zZ,1)
    X = parX[3]
    Y = parY[3]
    Z = parZ[3]
    
    
    Mdot  = setup['Mdot' ]
    vini  = setup['v_ini'] 
    
    #Select limits 
    if   Mdot <= 1e-5:
        rhoMin = 10**(-22.5)
        rhoMax = 10**(-13)
    elif Mdot >= 5e-7:
        rhoMin = 10**(-18.5)
        rhoMax = 10**(-11.5)
    elif 5e-7 < Mdot < 1e-5:
        rhoMin = 10**(-23)
        rhoMax = 10**(-12)

    vmin = 0
    vmax = 50
    Tmin = 1e1
    Tmax = 1e6


    #ax1.set_title('v = '+ str(vini)+ 'km/s', fontsize = 33)#, Mdot ='+ str(Mdot)+ '$M_\odot$/yr, ecc = ' +str(ecc))
    
    # Plots
    oneRadialStructurePlot(parX[0],parY[0], parZ[0], X, Y, Z, 'density [cm/g$^3$]','royalblue', ax1, rhoMin, rhoMax, xcomp, xAGB, setup['bound'])
    oneRadialStructurePlot(parX[1],parY[1], parZ[1], X, Y, Z, 'speed [km/s]'      ,'firebrick', ax2, vmin  , vmax  , xcomp, xAGB, setup['bound'])
    oneRadialStructurePlot(parX[2],parY[2], parZ[2], X, Y, Z, 'temperature [K]'   ,'goldenrod', ax3, Tmin  , Tmax  , xcomp, xAGB, setup['bound'])
    
    # Plot make up
    ax1.legend(handles = handl, fontsize = 12, loc = 'upper right')
    ax1.set_yscale('log')
    ax3.set_yscale('log')
    ax1.axes.get_xaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax3.set_xlabel('$r$ [AU]', fontsize = 16)

    
    fig.tight_layout()
    fig.subplots_adjust(wspace = 0.005,hspace = 0.01)
    fig.savefig(loc+str(run)+'_1Dplot_radialStructure.png')
    
    
    #Make text file with info to make plots
    title = loc+str(run)+'_data_1D_radialStructure.txt'
    with open (title,'w') as f:
        f.write('Model '+str(run)+'\n')
        f.write('Data to make radial structure plots yourself:')
        f.write('\n')
        names = ['X [cm]', 'Y[cm]', 'Z[cm]', 'Rho(x) [cm/g$^3$]', 'Rho(y) [cm/g$^3$]', 'Rho(z) [cm/g$^3$]', '|v|(x) [km/s]', '|v|(y) [km/s]', '|v|(z) [km/s]', 'T(x) [K]', 'T(y) [K]', 'T(z) [K]']
        f.write("{: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34}".format(*names))
        col_format = "{:<35}" * 12 + "\n"   # 7 left-justfied columns with 15 character width
        f.write('\n')
        for i in zip(X,Y,Z,parX[0],parY[0],parZ[0],parX[1],parY[1],parZ[1],parX[2],parY[2],parZ[2]):
            f.write(col_format.format(*i))
    
    print('     Radial structure plot model '+str(run)+' ready and saved!')



