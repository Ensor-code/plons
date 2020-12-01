import os
import math                     as math
import numpy                    as np
import matplotlib.pyplot        as plt
# own scripts
import ConversionFactors_cgs    as cgs
import GeometricalFunctions     as gf
import Tools                    as tl
import PhysicalQuantities       as pq
# ignore warnings
import warnings
warnings.filterwarnings("ignore")


'''
Double bin in a dataset 'rho' according to two sets (so in 2D):
    First bin in function of 'x',
    then bin every x-bin in function of 'xr'.
Calculate per bin the average 'rho'.
Finally to get the average 'rho' in a tube, select the data within
'x' < sma (= semi-major axis) and take the average per 'xr' 

RETURN:
    'tube' is a 2D list containing
        [0] the 'x'   values to plot
        [1] the 'rho' values to plot
'''
def getBinnedTube(rho, x, xr, sma, bound):    
    
    # bin density & xr ifo x
    binned_rho, nb_bins_rho = tl.bin_data(rho,x)
    
    binned_xr, nb_bins_xr   = tl.bin_data(xr, x)
    
    # bin density_binned ifo xr --> double binned 
    double_binned = {}
    
    for key in binned_rho:
        if len( binned_xr[key] ) != 0 and key < (bound - bound/20):
            innerbin, nb_bins   = tl.bin_data(binned_rho[key], binned_xr[key])
            double_binned[key]  = innerbin
            
    # for the binned density ifo x, not every xr-bin(x) will contain a particle --> empty lists
    # remove these empty bins from the binned_density
    empties = []
    # find empty bins
    for key1 in double_binned:
        for key2 in double_binned[key1]:
            if len(double_binned[key1][key2]) == 0:
                empties.append([key1,key2])
    # delete the empty bins from the dict
    for i in range(len(empties)):
        del double_binned[empties[i][0]][empties[i][1]]
    
    
    # calculate the average density per (x,xr)-bin
    rho_means = {}
    for key1 in double_binned:
        rho_means[key1] = []
        for key2 in double_binned[key1]:
            mean = np.mean(double_binned[key1][key2])
            rho_means[key1].append(mean)
            
    # calculate the average density in a tube of 2*sma around chosen axis & make two datasets for plotting
    r_tube   = []
    rho_tube = []
    for key in rho_means:
        r_tube.append(key)
        rho_tube.append(np.mean(rho_means[key][:int(round(sma*2))]))
        
    rho_tube_sm = tl.smoothen(rho_tube,2)
        
    tube = [r_tube, rho_tube_sm]
        
    return tube  
  
    
'''
Get the data for calculating the tubes.
Calculate the tube data:
    - z_tube:    a tube around the z-axis
    - orb_tube:  a tube containing the averaged data over the orbital plane
    
RETURN:
    Two 2D lists with the data in the tubes.
'''
def getTubeData_orbPlane(setup,data):
    
    bound = setup['bound'  ]        # [au]
    
    if setup['single_star'] == False:
        sma = setup['sma_ini']      # [au]
        
    if setup['single_star'] == True:
        sma = bound/100
    
    pos = data['position'].transpose()
    
    x   = pos[0]  /cgs.AU_cm()      # [au]    
    y   = pos[1]  /cgs.AU_cm()      # [au]
    z   = pos[2]  /cgs.AU_cm()      # [au]
    
    xr  = gf.calc_r_2D(y,z)         # [au]
    yr  = gf.calc_r_2D(x,z)         # [au]
    zr  = gf.calc_r_2D(x,y)         # [au]
    
    rho = data['rho']
    
    z_tube   = getBinnedTube(rho, z, zr, sma, bound)
    orb_tube = getBinnedTube(rho, zr, z, sma, bound)
    
    return z_tube, orb_tube
    
    
'''
Get the data for calculating the tubes.
Calculate the tube data:
    - x_tube:    a tube around the x-axis
    - y_tube:    a tube around the y-axis
    
        => take average of these two tubes in function of the radial coordinate
            to come to a representation for the orbital plane
    
RETURN:
    2D list with the data in the tube.
'''
def getTubeData_xy(setup,data):
    
    bound = setup['bound'  ]        # [au]
    
    if setup['single_star'] == False:
        sma = setup['sma_ini']      # [au]
        
    if setup['single_star'] == True:
        sma = bound/100
    
    pos = data['position'].transpose()
    
    x   = pos[0]  /cgs.AU_cm()      # [au]    
    y   = pos[1]  /cgs.AU_cm()      # [au]
    z   = pos[2]  /cgs.AU_cm()      # [au]
    
    xr  = gf.calc_r_2D(y,z)         # [au]
    yr  = gf.calc_r_2D(x,z)         # [au]
    zr  = gf.calc_r_2D(x,y)         # [au]
    
    rho = data['rho']
    
    x_tube   = getBinnedTube(rho, x, xr, sma, bound)
    y_tube   = getBinnedTube(rho, y, yr, sma, bound)

    
    xy_tube = []
    for i in range(len(x_tube[1])):
        xy_tube.append( (x_tube[1][i] + y_tube[1][i] ) /2 )
        
    orb_tube = [x_tube[0], xy_tube]
    
    return orb_tube

    
    
    
'''
MAIN def of this script.

Get the tube data and plots:
    - FIG1: for the average in the whole orbital plane      --> better
    - FIG2: for the average over the x- and y-tube
'''
def main_tube(run, outloc, setup, data):
    
    print('')
    print('(7)  Start calculations for the radial structure plots, using tubes.')
    
        
    ### --- fig using the orb plane 'tube'
    z_tube, orb_tube = getTubeData_orbPlane(setup,data)
    
    fig1 = plt.figure(figsize=(4,3.5))
    ax1  = plt.subplot(111)

    ax1.plot(orb_tube[0], orb_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'orbital plane')
    ax1.plot(z_tube[0]  , z_tube[1]   ,c='c',              lw = 0.8, label = 'polar axis')

    ax1.set_yscale('log')
    ax1.set_xlabel('$r$ [AU]', fontsize = 9)
    ax1.set_ylabel('density [cm/g$^3$]', fontsize = 9)
    ax1.tick_params(labelsize=7)

    ax1.legend(fontsize = 7, loc = 'upper right')
    fig1.tight_layout()
    
    plt.savefig(outloc+str(run)+'_1D_tube_orb.png', dpi = 200)
    
    
    ### --- fig using mean of the x and y tube
    
    #xy_tube = getTubeData_xy(setup,data)
        
    #fig2 = plt.figure(figsize=(4,3.5))
    #ax2  = plt.subplot(111)

    #ax2.plot(xy_tube[0] , xy_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'mean $x$ & $y$ axis')
    #ax2.plot(z_tube[0]  , z_tube[1]  ,c='c',              lw = 0.8, label = 'polar axis')

    #ax2.set_yscale('log')
    #ax2.set_xlabel('$r$ [AU]', fontsize = 9)
    #ax2.set_ylabel('density [cm/g$^3$]', fontsize = 9)
    #ax2.tick_params(labelsize=7)

    #ax2.legend(fontsize = 7, loc = 'upper right')
    #fig2.tight_layout()
    
    #plt.savefig(outloc+str(run)+'_1D_tube_xy.png', dpi = 200)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
