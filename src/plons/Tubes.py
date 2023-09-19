import os
import math                     as math
import numpy                    as np
import matplotlib.pyplot        as plt
# own scripts
import plons.ConversionFactors_cgs    as cgs
import plons.GeometricalFunctions     as gf
import plons.Tools                    as tl
# ignore warnings
import warnings
warnings.filterwarnings("ignore")


'''
Double bin in a dataset 'rho' according to two sets (so in 2D):
    First bin as a function of 'x',
    then bin every x-bin as a function of 'xr'.
Calculate per bin the average 'rho'.
Finally to get the average 'rho' in a tube, select the data within
'x' < sma (= semi-major axis) and take the average per 'xr'

RETURN:
    'tube' is a 2D list containing
        [0] the 'x'   values to plot
        [1] the 'rho' values to plot        # density in SI-units
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

    tube = [np.array(r_tube), np.array(rho_tube_sm)] #*cgs.gcm3_kgm3()]

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

    x   = np.abs( pos[0]  /cgs.au )     # [au]
    y   = np.abs( pos[1]  /cgs.au )     # [au]
    z   = np.abs( pos[2]  /cgs.au )     # [au]

    #xr  = np.abs( gf.calc_r_2D(y,z)    )     # [au]
    #yr  = np.abs( gf.calc_r_2D(x,z)    )     # [au]
    zr  = np.abs( gf.calc_r_2D(x,y)    )     # [au]       Distance to polar axis

    rho = data['rho']

    z_tube   = getBinnedTube(rho, z, zr, sma, bound)      #  Tube  for zr <= 2 sma, binned according to all z (:= r) values
    orb_tube = getBinnedTube(rho, zr, z, sma, bound)      # 'Tube' for z  <= 2 sma, binned according to all zr(:= r) values

    return z_tube, orb_tube


'''
Get the data for calculating the tubes.
Calculate the tube data:
    - x_tube:    a tube around the x-axis
    - y_tube:    a tube around the y-axis

        => take average of these two tubes as a function of the radial coordinate
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

    x   = np.abs( pos[0]  /cgs.au )     # [au]
    y   = np.abs( pos[1]  /cgs.au )     # [au]
    z   = np.abs( pos[2]  /cgs.au )     # [au]

    xr  = np.abs( gf.calc_r_2D(y,z)    )     # [au]
    yr  = np.abs( gf.calc_r_2D(x,z)    )     # [au]
    #zr  = np.abs( gf.calc_r_2D(x,y)    )     # [au]

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

    Mdot  = setup['Mdot' ]
    # Select ylimits as in 1Dplots
    #if   Mdot <= 5e-7:
        #rhoMin = 10**(-21)
        #rhoMax = 10**(-14)

    #elif Mdot >= 1e-5:
        #rhoMin = 10**(-18.5)
        #rhoMax = 10**(-11.5)

    #elif 5e-7 < Mdot < 1e-5:
        #rhoMin = 10**(-20)
        #rhoMax = 10**(-13)

    print('')
    print('(7)  Start calculations for the radial structure plots, using tubes.')


    ### --- fig using the orb plane 'tube'
    z_tube, orb_tube = getTubeData_orbPlane(setup,data)

    fig1 = plt.figure(figsize=(4,3))
    ax1  = plt.subplot(111)

    ax1.plot(orb_tube[0], orb_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'orbital plane')
    ax1.plot(z_tube[0]  , z_tube[1]   ,c='c',              lw = 0.8, label = 'polar axis')

    ax1.set_yscale('log')
    ax1.set_xlabel('$r$ [AU]', fontsize = 8)
    ax1.set_xlim([orb_tube[0][0],orb_tube[0][-1]])
    ax1.set_ylabel('mean density [g/cm$^3$]', fontsize = 8)
    ax1.tick_params(labelsize=7)

    ax1.legend(fontsize = 7, loc = 'upper right')
    fig1.tight_layout()

    plt.savefig(os.path.join(outloc, 'png/1D_tube_orb.png'), dpi = 200)
    plt.savefig(os.path.join(outloc, 'pdf/1D_tube_orb.pdf'))

    ### --- fig using the orb plane 'tube' & plotting the difference
    diff = []
    r = []

    for i in range(len(orb_tube[0])):
        for j in range(len(z_tube[0])):
            if orb_tube[0][i] == z_tube[0][j]:
                r.append(orb_tube[0][i])
                diff.append(np.log10(orb_tube[1][i]) - np.log10(z_tube[1][j]))
    diff= tl.smoothen(diff,2)


    ax_extra = ax1.twinx()
    ax_extra.plot(r,  diff            ,c='grey'        , lw = 0.3, label = 'difference' )
    ax_extra.plot([r[0],r[-1]], [0,0] ,c='grey', ls=':', lw = 0.4)
    color = 'tab:grey'
    ax_extra.set_ylabel('difference', fontsize = 9, color=color)
    ax_extra.set_ylim([-0.6,2.1])
    ax_extra.tick_params(axis='y',labelsize=7)

    ax1.legend(fontsize = 6, loc = 'upper center')
    # ax_extra.legend(fontsize = 5, loc = 'upper right')
    fig1.tight_layout()

    plt.savefig(os.path.join(outloc, 'png/1D_tube_orb_diff.png'), dpi = 200)
    plt.savefig(os.path.join(outloc, 'pdf/1D_tube_orb_diff.pdf'))


    ### --- fig using the orb plane 'tube' ONLY plotting the difference
    fig3 = plt.figure(figsize=(4,3.5))
    ax3  = plt.subplot(111)
    ax3.plot(r,  diff ,c='r', ls = '-'  , lw = 0.8, label = 'niet log')

    ax3.set_xlabel('$r$ [AU]', fontsize = 9)
    ax3.set_title('Mean density orb plane - mean density polar axis [g/cm$^3$]', fontsize = 9)
    ax3.tick_params(labelsize=7)

    fig3.tight_layout()

    plt.savefig(os.path.join(outloc, 'png/1D_tube_difference.png'), dpi = 200)
    plt.savefig(os.path.join(outloc, 'pdf/1D_tube_difference.pdf'))


    ### --- put data in a text file
    title = os.path.join(outloc, 'txt/data_1D_tube.txt')
    with open (title,'w') as f:
        f.write('Model '+str(run)+'\n')
        f.write('Data to make tube plots yourself:')
        f.write('\n')
        names = ['r [au]', 'MeanRhoOrbpl [g/cm$^3$]', 'MeanPolAx [g/cm$^3$]' , 'Difference logrho (OrbPl - PolAx)']
        f.write("{: <34} {: <34} {: <34} {: <34} ".format(*names))
        col_format = "{:<35}" * 4 + "\n"   # 4 left-justfied columns with 35 character width
        f.write('\n')
        for i in zip(z_tube[0], orb_tube[1], z_tube[1], (np.log(orb_tube[1]) - np.log(z_tube[1])) ):
            f.write(col_format.format(*i))



    ### --- fig using mean of the x and y tube

    #xy_tube = getTubeData_xy(setup,data)

    #fig2 = plt.figure(figsize=(4,3.5))
    #ax2  = plt.subplot(111)

    #ax2.plot(xy_tube[0] , xy_tube[1] ,c='k', ls = '--'  , lw = 0.8, label = 'mean $x$ & $y$ axis')
    #ax2.plot(z_tube[0]  , z_tube[1]  ,c='c',              lw = 0.8, label = 'polar axis')

    #ax2.set_yscale('log')
    #ax2.set_xlabel('$r$ [AU]', fontsize = 9)
    #ax2.set_ylabel('mean density [g/cm$^3$]', fontsize = 9)
    #ax2.tick_params(labelsize=7)

    #ax2.legend(fontsize = 7, loc = 'upper right')
    #fig2.tight_layout()

    #plt.savefig(os.path.join(outloc, 'png/1D_tube_xy.png'), dpi = 200)
    # plt.savefig(os.path.join(outloc, 'pdf/1D_tube_xy.pdf'))
    print('     Plot for radial structure, using tubes, of model '+str(run)+' ready and saved!')

















