import numpy                    as np
import matplotlib.pyplot        as plt
import os

# import plons scripts
import plons.SmoothingKernelScript    as sk
import plons.ConversionFactors_cgs    as cgs
import plons.PhysicalQuantities       as pq

# import certain things from packages
from matplotlib.ticker          import MultipleLocator
from scipy.signal import argrelextrema

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

n_grid = 10000
round_bounds = False
    
def getParamsLine(results_line, vec1, minT, gamma, vec2 = [], dumpData = None):
    indicesToKeep = np.where(np.log10(results_line['Tgas']) > minT)
    if gamma <= 1.:
        indicesToKeep = np.where(results_line['rho'] > 0.)

    R     = vec1[indicesToKeep]
    if len(vec2) != 0:
        x_comp = 0.
        y_comp = 0.
        if 'posComp' in dumpData:
            x_comp = dumpData['posComp'][0]
            y_comp = dumpData['posComp'][1]
        R = lineCoordinates(len(indicesToKeep[0]), vec1[indicesToKeep], vec2[indicesToKeep], x_comp, y_comp)

    temp  = (results_line['Tgas' ][indicesToKeep])
    speed = (results_line['speed'][indicesToKeep])
    rho   = (results_line['rho'  ][indicesToKeep])

    return(rho,speed,temp,R)

'''
definition used in plotParR, plots radial structure of given parameter (log(rho)/|v|/T) on the x- and y-axis
'''
def oneRadialStructurePlot(parX,parZ, X, Z, parName, axis, parMin, parMax, rcomp, bound, limX):
    axis.set_xlim(limX, bound)
    # Uncomment if you want both sides
    # axis.set_xlim(-bound,bound)
    axis.vlines(rcomp, parMin, parMax, 'black', linestyle='--', linewidth=1., label=r"$x_\mathrm{comp}$", zorder = 10)
    
    axis.plot(X/cgs.au, parX, color = 'C0', label = '$x$-axis', lw = 1)
    axis.plot(Z/cgs.au, parZ, color = 'C1', label = '$z$-axis', lw = 1)

    # ind_maxima = argrelextrema(parX,np.greater)
    # Xmax = X[ind_maxima]
    # RhoMax = parX[ind_maxima]
    # ind_maxima = argrelextrema(RhoMax,np.greater)
    # Xmax = Xmax[ind_maxima]
    # RhoMax = RhoMax[ind_maxima]
    # axis.plot(Xmax/cgs.au,RhoMax,'*','r')

    axis.set_ylim(parMin,parMax)

    axis.set_ylabel(parName, fontsize=22)
    axis.tick_params(labelsize=18)
    axis.tick_params(axis='x', pad=9)

def radialStructPlots(run,loc, dumpData, setup):
    #plots radial structure of log(rho), |v| and T on the x- and y-axis
    fig = None
    gamma = setup["gamma"]

    if gamma <= 1.:
        fig = plt.figure(figsize=(10, 10))

        # fig = plt.figure(figsize=(4.5, 10))
        ax1 = plt.subplot(311)
        ax2 = plt.subplot(312)

        # Remove data of AGB and companion
        dataToUse = {}
        dataToUse['rho'] = dumpData['rho'][:-2]
        dataToUse['Tgas'] = dumpData['Tgas'][:-2]
        dataToUse['speed'] = dumpData['speed'][:-2]
        dataToUse['mass'] = dumpData['mass'][:-2]
        dataToUse['position'] = dumpData['position'][:-2]
        dataToUse['h'] = dumpData['h'][:-2]

        # calculate smoothed data around one axis
        theta = 0.
        if setup['single_star'] == False:
            theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])

        pixCoordX      = sk.getPixels('line_x', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
        results_line_X = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], sk.rotatePixCoordAroundZ(theta, pixCoordX))
        pixCoordZ      = sk.getPixels('line_z', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
        results_line_Z = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], pixCoordZ)

        parX = getParamsLine(results_line_X, pixCoordX.transpose()[0], 1., gamma, pixCoordX.transpose()[1], dumpData)
        parZ = getParamsLine(results_line_Z, pixCoordZ.transpose()[2], 1., gamma)

        X = parX[3]
        Z = parZ[3]

        # Bounds
        
        rhoMinX, rhoMaxX = np.min(parX[0]), np.max(parX[0])
        rhoMinZ, rhoMaxZ = np.min(parZ[0]), np.max(parZ[0])
        rhoMin           = 0.1 * min(rhoMinX[rhoMinX>0], rhoMinZ[rhoMinZ>0])
        rhoMax           = 10 * max(rhoMaxX, rhoMaxZ)

        vMinX, vMaxX = np.min(parX[1]), np.max(parX[1])
        vMinZ, vMaxZ = np.min(parZ[1]), np.max(parZ[1])
        vMin = min(vMinX, vMinZ)
        vMin = 0.9 * vMin
        vMax = 1.1 * max(vMaxX, vMaxZ)

        TMinX, TMaxX = np.min(parX[2]), np.max(parX[2])
        TMinZ, TMaxZ = np.min(parZ[2]), np.max(parZ[2])
        TMin = 0.9 * min(TMinX, TMinZ)
        TMax = 1.1 * max(TMaxX, TMaxZ)

        # ax1.set_title('v = '+ str(vini)+ 'km/s', fontsize = 33)#, Mdot ='+ str(Mdot)+ '$M_\odot$/yr, ecc = ' +str(ecc))
        # Plots
        limX = setup['wind_inject_radius']
        posAGB = np.hypot(dumpData['posAGB'][0], dumpData['posAGB'][1])
        X += posAGB
        posComp = posAGB
        if setup['single_star'] == False:
            posComp = (np.hypot(dumpData['posComp'][0], dumpData['posComp'][1])) / cgs.au
            if setup['triple_star']==True:
                posComp = [posComp,(np.hypot(dumpData['posComp_in'][0],dumpData['posComp_in'][1])) / cgs.au]
            print('radial position companion(s): ', posComp)


        bound = setup['bound']
        # bound = 500

        oneRadialStructurePlot(parX[0],parZ[0], X, Z, r'$\rho$ [g$\,$cm$^{-3}$]', ax1, rhoMin, rhoMax,
                               posComp, bound, limX)
        oneRadialStructurePlot(parX[1], parZ[1], X, Z, r'$v$ [km/s]', ax2, vMin, vMax, posComp, bound,
                               limX)

        # Plot make up
        ax1.xaxis.set_major_locator(MultipleLocator(bound / 5.))
        ax1.xaxis.set_minor_locator(MultipleLocator(bound / 30.))
        ax2.xaxis.set_major_locator(MultipleLocator(bound / 5.))
        ax2.xaxis.set_minor_locator(MultipleLocator(bound / 30.))

        ax1.legend(loc='upper center', ncol=5, bbox_to_anchor=[0., 0.3, 1., 1.], prop={'size': 20}, labelspacing=2)
        ax1.set_yscale('log')
        #ax2.set_yscale('log')
        ax2.set_xlabel('$r$ [AU]', fontsize=22)
        ax1.set_xticklabels([])

        '''
        # Make text file with info to make plots
        title = os.path.join(loc, 'txt/data_1D_radialStructure.txt')
        with open(title, 'w') as f:
            f.write('Model ' + str(run) + '\n')
            f.write('Data to make radial structure plots yourself:')
            f.write('\n')
            names = ['X [cm]', 'Z[cm]', 'Rho(x) [g/cm$^3$]', 'Rho(z) [g/cm$^3$]', '|v|(x) [km/s]', '|v|(z) [km/s]',
                     'T(x) [K]', 'T(z) [K]']
            f.write("{: <34} {: <34} {: <34}  {: <34} {: <34} {: <34}".format(*names))
            col_format = "{:<35}" * 6 + "\n"  # 7 left-justfied columns with 15 character width
            f.write('\n')
            for i in zip(X, Z, parX[0], parZ[0], parX[1], parZ[1]):
                f.write(col_format.format(*i))
        '''

    else:
        fig = plt.figure(figsize=(15,10))

        #fig = plt.figure(figsize=(4.5, 10))
        ax1 = plt.subplot(311)
        ax2 = plt.subplot(312)
        ax3 = plt.subplot(313)

        #Remove data of AGB and companion
        dataToUse = {}
        dataToUse['rho'     ] = dumpData['rho'     ][:-2]
        dataToUse['Tgas'    ] = dumpData['Tgas'    ][:-2]
        dataToUse['speed'   ] = dumpData['speed'   ][:-2]
        dataToUse['mass'    ] = dumpData['mass'    ][:-2]
        dataToUse['position'] = dumpData['position'][:-2]
        dataToUse['h'       ] = dumpData['h'       ][:-2]

        #calculate smoothed data around one axis
        theta = 0.
        if setup['single_star'] == False:
            theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])

        pixCoordX      = sk.getPixels('line_x', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
        results_line_X = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], sk.rotatePixCoordAroundZ(theta, pixCoordX))
        pixCoordZ      = sk.getPixels('line_z', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
        results_line_Z = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], pixCoordZ)

        gamma = setup["gamma"]
        parX = getParamsLine(results_line_X, pixCoordX.transpose()[0], 1., gamma, pixCoordX.transpose()[1], dumpData)
        parZ = getParamsLine(results_line_Z, pixCoordZ.transpose()[2], 1., gamma)

        X = parX[3]
        Z = parZ[3]

        # Bounds
        rhoMinX, rhoMaxX = np.min(parX[0]), np.max(parX[0])
        rhoMinZ, rhoMaxZ = np.min(parZ[0]), np.max(parZ[0])
        rhoMin           = 0.1 * min(rhoMinX[rhoMinX>0], rhoMinZ[rhoMinZ>0])
        rhoMax           = 10 * max(rhoMaxX, rhoMaxZ)
        #print('rhomin en max zijn: ',rhoMin,rhoMax)

        vMinX, vMaxX = np.min(parX[1]), np.max(parX[1])
        vMinZ, vMaxZ = np.min(parZ[1]), np.max(parZ[1])
        vMin         = min(vMinX,vMinZ)
        vMin         = 0.9 * vMin
        vMax         = 1.1*max(vMaxX, vMaxZ)
        #print(vMin,vMax)

        TMinX, TMaxX = np.min(parX[2]), np.max(parX[2])
        TMinZ, TMaxZ = np.min(parZ[2]), np.max(parZ[2])
        TMin         = 0.9 * min(TMinX, TMinZ)
        TMax         = 1.1 * max(TMaxX, TMaxZ)

        #ax1.set_title('v = '+ str(vini)+ 'km/s', fontsize = 33)#, Mdot ='+ str(Mdot)+ '$M_\odot$/yr, ecc = ' +str(ecc))
        # Plots
        limX = setup['wind_inject_radius']
        posAGB = np.hypot(dumpData['posAGB'][0], dumpData['posAGB'][1])
        X += posAGB
        posComp = posAGB
        if setup['single_star'] == False:
            posComp = (np.hypot(dumpData['posComp'][0], dumpData['posComp'][1])+posAGB) / cgs.au
            if setup['triple_star']==True:
                posComp = [posComp,(np.hypot(dumpData['posComp_in'][0],dumpData['posComp_in'][1])+posAGB) / cgs.au]
            print('radial position 2 companions: ', posComp)
        bound = setup['bound']
        # bound = 500

        oneRadialStructurePlot(parX[0], parZ[0], X, Z, r'$\rho$ [g$\,$cm$^{-3}$]', ax1, rhoMin, rhoMax, posComp,  bound, limX)
        oneRadialStructurePlot(parX[1], parZ[1], X, Z, r'$v$ [km/s]', ax2, vMin  , vMax  , posComp, bound, limX)
        oneRadialStructurePlot(parX[2], parZ[2], X, Z, r'$T$ [K]', ax3, TMin  , TMax  , posComp, bound, limX)

        # Plot make up
        ax1.xaxis.set_major_locator(MultipleLocator(bound / 5.))
        ax1.xaxis.set_minor_locator(MultipleLocator(bound / 30.))
        ax2.xaxis.set_major_locator(MultipleLocator(bound / 5.))
        ax2.xaxis.set_minor_locator(MultipleLocator(bound / 30.))
        ax3.xaxis.set_major_locator(MultipleLocator(bound / 5.))
        ax3.xaxis.set_minor_locator(MultipleLocator(bound / 30.))

        ax1.legend(loc = 'upper center', ncol=5, bbox_to_anchor=[0., 0.3, 1., 1.], prop={'size': 20}, labelspacing=2)
        ax1.set_yscale('log')
        #ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax3.set_xlabel('$r$ [AU]', fontsize=22)
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])

        '''
        # Make text file with info to make plots
        title = os.path.join(loc, 'txt/data_1D_radialStructure.txt')
        with open (title,'w') as f:
            f.write('Model '+str(run)+'\n')
            f.write('Data to make radial structure plots yourself:')
            f.write('\n')
            names = ['X [cm]', 'Z[cm]', 'Rho(x) [g/cm$^3$]', 'Rho(z) [g/cm$^3$]', '|v|(x) [km/s]', '|v|(z) [km/s]', 'T(x) [K]', 'T(z) [K]']
            f.write("{: <34} {: <34} {: <34}  {: <34} {: <34} {: <34} {: <34} {: <34}".format(*names))
            col_format = "{:<35}" * 8 + "\n"   # 7 left-justfied columns with 15 character width
            f.write('\n')
            for i in zip(X,Z,parX[0],parZ[0],parX[1],parZ[1],parX[2],parZ[2]):
                f.write(col_format.format(*i))
        '''
    print('     Radial structure plot model '+str(run)+' ready and saved!')

    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.savefig(os.path.join(loc, 'png/1Dplot_radialStructure.png'), bbox_inches="tight")
    fig.savefig(os.path.join(loc, 'pdf/1Dplot_radialStructure.pdf'), bbox_inches="tight")

def lineCoordinates(n, x, y, x_comp, y_comp):
    r = np.zeros(shape=(n))
    dot_product = x * x_comp + y * y_comp
    r[dot_product < 0] = -np.hypot(x[dot_product < 0], y[dot_product < 0])
    r[dot_product >= 0] = np.hypot(x[dot_product >= 0], y[dot_product >= 0])

    return r


def plotRho_r2(run,loc, dumpData, setup):

    fig = plt.figure(figsize=(10, 4.5))
    ax1 = plt.subplot(111)

    #Remove data of AGB and companion
    dataToUse = {}
    dataToUse['rho'     ] = dumpData['rho'     ][:-2]
    dataToUse['position'] = dumpData['position'][:-2]

    #calculate smoothed data around one axis
    theta = 0.
    if setup['single_star'] == False:
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])

    pixCoordX      = sk.getPixels('line_x', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
    results_line_X = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], sk.rotatePixCoordAroundZ(theta, pixCoordX))
    pixCoordZ      = sk.getPixels('line_z', n_grid, 'comp', dumpData,  setup['bound'] * cgs.au)
    results_line_Z = sk.getSmoothingKernelledPix(20, dumpData, ['rho', 'Tgas', 'speed'], pixCoordZ)

    gamma = setup["gamma"]
    parX = getParamsLine(results_line_X, pixCoordX.transpose()[0], 1., gamma, pixCoordX.transpose()[1], dumpData)
    parZ = getParamsLine(results_line_Z, pixCoordZ.transpose()[2], 1., gamma)

    X = parX[3]
    Z = parZ[3]


    #ax1.set_title('v = '+ str(vini)+ 'km/s', fontsize = 33)#, Mdot ='+ str(Mdot)+ '$M_\odot$/yr, ecc = ' +str(ecc))
    # Plots
    limX = setup['wind_inject_radius']
    posAGB = np.hypot(dumpData['posAGB'][0], dumpData['posAGB'][1])
    X += posAGB
    posComp = posAGB
    if setup['single_star'] == False:
        posComp = (np.hypot(dumpData['posComp'][0], dumpData['posComp'][1])+posAGB) / cgs.au
        if setup['triple_star']==True:
            posComp = [posComp,(np.hypot(dumpData['posComp_in'][0],dumpData['posComp_in'][1])+posAGB) / cgs.au]
        print('radial position 2 companions: ', posComp)
    bound = setup['bound']
    # bound = 500

    parX = parX[0]*X**2
    parZ = parZ[0]*Z**2
    # Bounds
    MinX, MaxX = np.min(parX), np.max(parX)
    MinZ, MaxZ = np.min(parZ), np.max(parZ)
    parMin     = min(MinX[MinX>0], MinZ[MinZ>0])
    parMax     = 10 * max(MaxX, MaxZ)

    oneRadialStructurePlot(parX, parZ, X, Z, r'$\rho * r^2$ [g$\,$cm$^{-1}$]', ax1, parMin, parMax, posComp,  bound, limX)

    # Plot make up
    ax1.xaxis.set_major_locator(MultipleLocator(bound / 5.))
    ax1.xaxis.set_minor_locator(MultipleLocator(bound / 30.))

    ax1.legend(loc = 'upper center', ncol=5, bbox_to_anchor=[0., 0.3, 1., 1.], prop={'size': 20}, labelspacing=2)
    ax1.set_yscale('log')
    ax1.set_xlabel('$r$ [AU]', fontsize=22)


    print('     Radial structure plot model '+str(run)+' ready and saved!')

    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.savefig(os.path.join(loc, 'png/1Dplot_rho_r2.png'), bbox_inches="tight")
    fig.savefig(os.path.join(loc, 'pdf/1Dplot_rho_r2.pdf'), bbox_inches="tight")