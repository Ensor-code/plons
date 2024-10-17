import matplotlib.pyplot as plt
import numpy as np
import os
from typing import Dict, List

import plons.SlicePlots2D       as sl
import plons.PhysicalQuantities as pq


# ***************************************
# ************* Slice Plots *************
# ***************************************

round_bounds    = False
minInfLog       = -300
minInf          = 10.**(-minInfLog)
maxInfLog       = 300
maxInf          = 10.**maxInfLog
sigma_bounds_u  = 3.
sigma_bounds_l  = 2.
customRanges = True

'''
main definition
'''
def SlicePlots(run, loc, dumpData, setup, number = -1, zoomin = [1, 5, 10], 
               observables = ['rho', 'Tgas', 'speed'],
               nneighb = 10, n_grid = 200, n_grid_vec = 25, printout=False):
    theta=0
    if not setup["single_star"]:
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
    
    if printout: print('     Calculating the smoothing kernels. This may take a while, please wait...')
    smooth = {}
    smooth_vec = {}
    for zoom in zoomin:
        smooth[zoom], smooth_vec[zoom] = sl.smoothData(dumpData, setup, observables, theta, zoom, nneighb, n_grid, n_grid_vec)

        if not customRanges:
            limits = makeLimits(observables, smooth, zoom)   
        else:
            limits = customRanges(observables, zoomin) 
        if printout:
            print("          Ranges of Parameters: zoom = "+str(zoom))
            if "rho"     in observables: print("          rhoMin,   rhoMax   = {0:10.5f}, {1:10.5f}".format(limits["rho"][zoom][0], limits["rho"][zoom][1]))
            if "speed"   in observables: print("          vMin,     vMax     = {0:10.5f}, {1:10.5f}".format(limits["speed"][zoom][0], limits["speed"][zoom][1]))
            if "Tgas"    in observables: print("          TMin,     TMax     = {0:10.5f}, {1:10.5f}".format(limits["Tgas"][zoom][0], limits["Tgas"][zoom][1]))
            if "kappa"   in observables: print("          kappaMin, kappaMax = {0:10.5f}, {1:10.5f}".format(limits["kappa"][zoom][0], limits["kappa"][zoom][1]))
            if "Gamma_E" in observables: print("          GammaMin, GammaMax = {0:10.5f}, {1:10.5f}".format(limits["Gamma_E"][zoom][0], limits["Gamma_E"][zoom][1]))
            if "tau"     in observables: print("          tauMin,   tauMax   = {0:10.5f}, {1:10.5f}".format(limits["tau"][zoom][0], limits["tau"][zoom][1]))

        # Make plots
        fig = sl.densityPlot(smooth, zoom, limits["rho"][zoom], dumpData, setup, orbital=True)
        saveFig(fig, loc, '2Dplot_density_orbital', zoom, number)
        fig = sl.densityPlot(smooth, zoom, limits["rho"][zoom], dumpData, setup, orbital=False)
        saveFig(fig, loc, '2Dplot_density_meridional', zoom, number)
        fig = sl.allPlots(smooth, smooth_vec, zoom, limits, dumpData, setup, observables)
        struct = ""
        for i in observables:
            struct+=i
        saveFig(fig, loc, '2Dplot_'+struct, zoom, number)
        plt.close()
        if printout: print('          Slice plots (zoom factor = ' + str(zoom) + ') model ' + str(run) + ' ready and saved!\n')

def saveFig(fig, loc, name, zoom, number):
    if number == -1:
        fig.savefig(os.path.join(loc, 'png/'+name+'_zoom{0:01d}.png'.format(zoom)), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(loc, 'pdf/'+name+'_zoom{0:01d}.pdf'.format(zoom)), dpi=300, bbox_inches="tight")
    else:
        fig.text(0.5, 0.9, "Dumpfile {0:05d}".format(number), size=28)
        fig.savefig(os.path.join(loc, 'animation/'+name+'_zoom{0:01d}_{1:04d}.png'.format(zoom,number)), dpi=200, bbox_inches="tight")
    plt.close()

def makeLimits(observables, smooth, zoom):
    limits = {}
    for observable in observables:
        limits[observable] = {}
        limits[observable][zoom] = [0.,0.]

    if "rho" in observables:
        limits["rho"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_y']["rho"]), log=True, round=round_bounds)
    if "speed" in observables:
        limits["speed"][zoom] = findBounds(smooth[zoom]['smooth_y']["speed"], log=False, round=round_bounds)
        limits["speed"][zoom][0] = max(limits["speed"][zoom][0], 0.)
    if "Tgas" in observables:
        limits["Tgas"][zoom] = findBounds(np.log10(smooth[zoom]['smooth_z']["Tgas"]), log=True, round=round_bounds)
    if "kappa" in observables:
        limits["kappa"][zoom] = findBounds(smooth[zoom]['smooth_y']["kappa"], log=False, round=round_bounds)
    if "Gamma_E" in observables:
        limits["Gamma_E"][zoom] = findBounds(smooth[zoom]['smooth_y']["Gamma_E"], log=False, round=round_bounds)
    if "tau" in observables:
        limits["tau"][zoom] = findBounds(smooth[zoom]['smooth_y']["tau"], log=True, round=round_bounds)
        limits["tau"][zoom][0] = max(limits["tau"][zoom][0], 0.)

    return limits
    
def findBounds(data, log = False, round = False):
    filtered_data = data
    result = np.zeros(2)

    if (-np.inf in data) or (np.inf in data):
        min = minInf
        max = maxInf
        if log == True:
            min = minInfLog
            max = maxInfLog

        filtered_data = filtered_data[np.logical_and(min <= data, data <= max)]

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

def customRanges(observables, zoomin):
    limits: Dict[str, Dict[int, List]] = {}
    for observable in observables:
        limits[observable] = {}
        for zoom in zoomin:
           limits[observable][zoom] = [0.,0.]
    if customRanges:
        if "rho" in observables:
            limits["rho"][1]  = [-19, -14]
            limits["rho"][2]  = [-19, -14]
            limits["rho"][4]  = [-17, -14]
            limits["rho"][5]  = [-17, -14]
            limits["rho"][10] = [-19, -14]
            limits["rho"][20] = [-18, -13]

        if "speed" in observables:
            limits["speed"][1]  = [0., 20.]
            limits["speed"][2]  = [0., 20.]
            limits["speed"][5]  = [0., 20.]
            limits["speed"][10] = [0., 20.]
            limits["speed"][20] = [0., 20.]

        if "Tgas" in observables:
            limits["Tgas"][1]  = [1., 4.]
            limits["Tgas"][2]  = [1., 4.]
            limits["Tgas"][5]  = [1., 4.]
            limits["Tgas"][10] = [1., 4.]
            limits["Tgas"][20] = [1., 4.]

        if "tau" in observables:
            limits["tau"][1]  = [0, 1]
            limits["tau"][2]  = [0, 1]
            limits["tau"][5]  = [0, 1]
            limits["tau"][10] = [0, 1]
            limits["tau"][20] = [0, 1]

        if "Tdust" in observables:
            limits["Tdust"][1]  = [0, 2000]
            limits["Tdust"][2]  = [0, 2000]
            limits["Tdust"][5]  = [0, 2000]
            limits["Tdust"][10] = [0, 2000]
            limits["Tdust"][20] = [0, 2000]

        if "tauL" in observables:
            limits["tauL"][1]  = [0, 0.05]
            limits["tauL"][2]  = [0, 0.05]
            limits["tauL"][5]  = [0, 0.05]
            limits["tauL"][10] = [0, 0.05]
            limits["tauL"][20] = [0, 0.05]

        if "kappa" in observables:
            limits["kappa"][1]  = [0., 3.]
            limits["kappa"][2]  = [0., 3.]
            limits["kappa"][5]  = [0., 3.]
            limits["kappa"][10] = [0., 3.]
            limits["kappa"][20] = [0., 3.]

        if "Gamma_E" in observables:
            limits["Gamma_E"][1]  = [0., 1.]
            limits["Gamma_E"][2]  = [0., 1.]
            limits["Gamma_E"][5]  = [0., 1.]
            limits["Gamma_E"][10] = [0., 1.]
            limits["Gamma_E"][20] = [0., 1.]
    return limits


