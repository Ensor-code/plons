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


   
    
def getTubeData_orbPlane(setup,data):
    
    if setup['single_star'] == False:
        sma = setup['sma_ini']
        
    if setup['single_star'] == True:
        sma = setup['bound']/100
    
    pos = data['position'].transpose()
    
    x   = pos[0]  /cgs.AU_cm()      # [au]    
    y   = pos[1]  /cgs.AU_cm()      # [au]
    z   = pos[2]  /cgs.AU_cm()      # [au]
    
    xr  = gf.calc_r_2D(y,z)         # [au]
    yr  = gf.calc_r_2D(x,z)         # [au]
    zr  = gf.calc_r_2D(x,y)         # [au]
    
    rho = data['rho']
    
    z_tube   = getBinnedTube(rho, z, zr, sma)
    orb_tube = getBinnedTube(rho, zr, z, sma)
    
    return orb_tube, z_tube
    
    
    
def getTubeData_xy(setup,data):
    
    if setup['single_star'] == False:
        sma = setup['sma_ini']
        
    if setup['single_star'] == True:
        sma = setup['bound']/100
    
    pos = data['position'].transpose()
    
    x   = pos[0]  /cgs.AU_cm()      # [au]    
    y   = pos[1]  /cgs.AU_cm()      # [au]
    z   = pos[2]  /cgs.AU_cm()      # [au]
    
    xr  = gf.calc_r_2D(y,z)         # [au]
    yr  = gf.calc_r_2D(x,z)         # [au]
    zr  = gf.calc_r_2D(x,y)         # [au]
    
    rho = data['rho']
    
    x_tube = getBinnedTube(rho, x, xr, sma)
    y_tube = getBinnedTube(rho, y, yr, sma)

    
    xy_tube = []
    for i in range(len(x_tube[1])):
        xy_tube.append( (x_tube[1][i] + y_tube[1][i] ) /2 )
        
    orb_tube = [x_tube[0], xy_tube]
    
    return orb_tube

    
    
def getBinnedTube(rho, x, xr, sma):    
    
    # bin density & xr ifo x
    binned_rho, nb_bins_rho = tl.bin_data(rho,x)
    
    binned_xr, nb_bins_xr   = tl.bin_data(xr, x)
    
    # bin density_binned ifo xr --> double binned 
    double_binned = {}
    
    for key in binned_rho:
        if len( binned_xr[key] ) != 0:
            innerbin, nb_bins   = tl.bin_data(binned_rho[key], binned_xr[key])
            double_binned[key]  = innerbin
            
            
    # for the binned density ifo x, not every xr-bin(x) will contain a particle --> empty list
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
    
    
    
def main_tube(run, outloc, setup, data):
    
    print('')
    print('(7)  Start calculations for the radial structure plots, using tubes.')
    
    orb_tube, z_tube = getTubeData_orbPlane(setup,data)
    xy_tube          = getTubeData_xy(setup,data)
    
    
    # fig using the orb plane 'tube'
    fig1 = plt.figure(figsize=(4,3.5))
    ax1  = plt.subplot(111)

    ax1.plot(orb_tube[0], orb_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'orbital plane')
    ax1.plot(z_tube[0]  , z_tube[1]   ,c='c',              lw = 0.8, label = 'polar axis')

    ax1.set_yscale('log')
    ax1.set_xlabel('$r$ [AU]', fontsize = 9)
    ax1.set_ylabel('density [cm/g$^3$]', fontsize = 9)
    ax1.tick_params(labelsize=7)

    ax1.legend(fontsize = 7, loc = 'upper right')
    
    plt.savefig(outloc+str(run)+'1D_tube_orb.png', dpi = 200)
    
    
    
    # fig using mean of the x and y tube
    fig2 = plt.figure(figsize=(4,3.5))
    ax2  = plt.subplot(111)

    ax2.plot(xy_tube[0] , xy_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'mean $x$ & $y$ axis')
    ax2.plot(z_tube[0]  , z_tube[1]  ,c='c',              lw = 0.8, label = 'polar axis')

    ax2.set_yscale('log')
    ax2.set_xlabel('$r$ [AU]', fontsize = 9)
    ax2.set_ylabel('density [cm/g$^3$]', fontsize = 9)
    ax2.tick_params(labelsize=7)

    ax2.legend(fontsize = 7, loc = 'upper right')
    
    plt.savefig(outloc+str(run)+'1D_tube_xy.png', dpi = 200)
    
    
    
    
    
    
    
    
    
    
    
    
    
