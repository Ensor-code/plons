#import datetime                     as dt
import numpy                        as np
import os
import warnings
warnings.filterwarnings("ignore")
import argparse                     as ap
import matplotlib.pyplot            as plt
from typing import Dict, List

# import plons scripts
import plons.RadialStructurePlots1D       as rs
import plons.SlicePlots2D                 as sl
import plons.CMF_meanRho                  as cmf
import plons.OrbitalEvolution             as ov
import plons.PhysicalQuantities           as pq
import plons.LoadData                     as ld
import plons.TerminalVelocity             as tmv
import plons.Profiles1D                   as dp
import userSettings                       as us
import plons.ArchimedianSpiral            as ars


round_bounds    = False
minInfLog       = -300
minInf          = 10.**(-minInfLog)
maxInfLog       = 300
maxInf          = 10.**maxInfLog
sigma_bounds_u  = 3.
sigma_bounds_l  = 2.

#print('------------------START:', dt.datetime.now(),'---------------------')
print('')

options = { '1': '(1) 2D slice plots',
            '2': '(2) 1D line plots & tube plots',
            '3': '(3) velocity related quantities and Morphological parameters',
            '4': '(4) Cummulative mass fraction',
            '5': '(5) Orbital evolution',
            '6': '(6) 1D spherical model',
            '7': '(7) ArchimedianSpiral'
            }


# ***************************************
# ************ USER INPUT ***************
# ***************************************

p = ap.ArgumentParser(description="This pipeline reduces PHANTOM output to usable plots and datasets.", formatter_class=ap.RawTextHelpFormatter)
p.add_argument("-m", "--models", action="store", default='', type = str, nargs='*', dest = "models", help = "The integers of models you want to run, split with spaces. To see the corresponding integers, run without this option")
helpoptions = "The integers of options you want to run:"
for option in options:
    helpoptions+='\n'+options[option]
p.add_argument("-o", "--options", action="store", default='', type = str, nargs='*', dest = "options", help = helpoptions)
p.add_argument("-z", "--zoom", action="store", default='1 5 10', type = str, nargs='*', dest = "zoom")
args = p.parse_args()
runModels = args.models
runParts = args.options
zoomin = np.array(args.zoom.split(), dtype=int)

def run_main(outputloc,runParts,numbers, models):
    for number in numbers:
        run = models[int(number)][1]
        print('------ MODEL '+str(number)+': '+run+' ------')
        saveloc = os.path.join(outputloc, run)
        try:
            os.makedirs(os.path.join(saveloc, 'png'))
            os.makedirs(os.path.join(saveloc, 'txt'))
            os.makedirs(os.path.join(saveloc, 'pdf'))
            os.makedirs(os.path.join(saveloc, 'animation'))
        except OSError:
            pass

        print('')
        print('Data is loading...')
        [setup, dumpData, sinkData, outerData] = LoadData_cgs(run, loc, userSettingsDictionary, bound=bound, factor=factor, runPart=runParts)
        print('All data is loaded and ready to use.')
        print('')
        for part in runParts:
            if part == '0':
                for i in options:
                    runPart(i, run, saveloc, dumpData, setup, sinkData, outerData)
            else:
                runPart(part, run, saveloc, dumpData, setup, sinkData, outerData)
        print('')

def searchModels(loc, prefix):
    result = []
    for path, directories, files in os.walk(loc):
        dumpFiles = list(filter(lambda x: (prefix+'_' in x and not '.tmp' in x), files))
        if len(dumpFiles) != 0:
            slicedString = path.replace(loc, "")
            result.append([path, slicedString])
    return sorted(result)

def runPart(part, run, saveloc, dumpData, setup, sinkData, outerData):
    if part == '1':
        print('')
        print('(1)  Start calculations for slice plots...')
        if customRanges:
            SlicePlots(run, saveloc, dumpData, setup, zoomin=zoomin, observables=observables, limits=limits, printout=True)
        else:
            SlicePlots(run, saveloc, dumpData, setup, zoomin=zoomin, observables=observables, printout=True)

    if part == '2':
        print('')
        print('(2)  Start calculations for the radial structure plots.')
        rs.radialStructPlots(run, saveloc, dumpData, setup)

    if part == '3':
        print('')
        print('(3) Start calculations for terminal velocity...')
        tmv.main_terminalVelocity(setup, dumpData, sinkData, saveloc, run)

    if part == '4':
        print('')
        print('(4)  Start calculations for the cummulative mass fraction and mean density plots...')
        if setup['single_star'] == True:
            cmf.CMF_meanRho(run, saveloc, dumpData, setup, factor)
        else:
            cmf.CMF_meanRho(run, saveloc, outerData, setup, factor)

    if part == '5':
        print('')
        if setup['single_star']:
            print('(5)  A single model has no orbit, and thereby no orbital evolution.')
            print('     The orbital evolution part is therefore skipped.')
        else:
            print('(5)  Start calculations for orbital evolution...')
            ov.orbEv_main(run, saveloc, sinkData, setup)

    if part == '6':
        print('')
        print('(6)  Start calculating for the 1D spherical plots')
        dp.profiles_main(run, loc, saveloc, dumpData, setup)

    if part == '7':
        print('')
        print('(7)  Archimedian spiral')
        ars.ArchimedianSpiral(run, saveloc, setup)

def LoadData_cgs(run, loc, userSettings, bound = None, factor = -1, number = -1, runPart = [0]):
    dir       = os.path.join(loc, run)
    setup     = ld.LoadSetup(dir, userSettings["prefix"])

    # Pick either last dump file or user chosen file
    if number == -1: index = ld.lastFullDumpIndex(dir, userSettings["prefix"], setup)
    else: index = number
    fileName       = dir+'/{0:s}_{1:05d}'.format(userSettings["prefix"], index)
    dumpData  = ld.LoadFullDump(fileName, setup)

    if len(set(runPart)-set([1, 2, 4, 6, 7]))>0:
        sinkData  = ld.LoadSink(dir, userSettings['prefix'], setup["icompanion_star"])
        if bound == None:
            bound = setup['bound']
        if factor > 0:
            outerData = ld.LoadDump_outer(factor, bound, setup, dumpData)
        else: outerData = None

    # save the final specifics of the AGB star to dumpData
    dumpData['maccrAGB' ] = sinkData['maccrAGB'   ][-1]

    # save the final specifics of the companion to dumpData
    if not setup['single_star']:
        dumpData['maccrComp'] = sinkData['maccrComp'  ][-1]
        dumpData['v_orbAGB' ] = sinkData['v_orbAGB_t' ][-1]
        dumpData['v_orbComp'] = sinkData['v_orbComp_t'][-1]
        dumpData['sma_fi'   ] = dumpData['rAGB'] + dumpData['rComp']        # [cm]

    if setup['triple_star']:
        dumpData['maccrComp_in'] = sinkData['maccrComp_in'  ][-1]
        dumpData['v_orbComp_in'] = sinkData['v_orbComp_in_t'][-1]

    return setup, dumpData, sinkData, outerData



# ***************************************
# ************* Slice Plots *************
# ***************************************

'''
main definition
'''
def SlicePlots(run, loc, dumpData, setup, number = -1, zoomin = [1, 5, 10], 
               observables = ['rho', 'Tgas', 'speed'], limits = False,
               nneighb = 10, n_grid = 200, n_grid_vec = 25, printout=False):
    if limits: customRanges = True
    else:      customRanges = False
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
        if printout:
            print("          Ranges of Parameters: zoom = "+str(zoom))
            if "rho"   in observables: print("          rhoMin,   rhoMax   = {0:10.5f}, {1:10.5f}".format(limits["rho"][zoom][0], limits["rho"][zoom][1]))
            if "speed" in observables: print("          vMin,     vMax     = {0:10.5f}, {1:10.5f}".format(limits["speed"][zoom][0], limits["speed"][zoom][1]))
            if "Tgas"  in observables: print("          TMin,     TMax     = {0:10.5f}, {1:10.5f}".format(limits["Tgas"][zoom][0], limits["Tgas"][zoom][1]))
            if "kappa" in observables: print("          kappaMin, kappaMax = {0:10.5f}, {1:10.5f}".format(limits["kappa"][zoom][0], limits["kappa"][zoom][1]))
            if "Gamma" in observables: print("          GammaMin, GammaMax = {0:10.5f}, {1:10.5f}".format(limits["Gamma"][zoom][0], limits["Gamma"][zoom][1]))
            if "tau"   in observables: print("          tauMin,   tauMax   = {0:10.5f}, {1:10.5f}".format(limits["tau"][zoom][0], limits["tau"][zoom][1]))

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
    if "Gamma" in observables:
        limits["Gamma"][zoom] = findBounds(smooth[zoom]['smooth_y']["Gamma"], log=False, round=round_bounds)
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

def customRanges(observables):
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

        if "Gamma" in observables:
            limits["Gamma"][1]  = [0., 1.]
            limits["Gamma"][2]  = [0., 1.]
            limits["Gamma"][5]  = [0., 1.]
            limits["Gamma"][10] = [0., 1.]
            limits["Gamma"][20] = [0., 1.]
    return limits




print('')
print('-------------------------------------------------------')
print('     Welcome to PLONS!'        )
print('-------------------------------------------------------')
print('')
print('The PLOtting tool for Nice Simulations.')
print('')


# Initialise user settings or load them
userSettingsFilePath = os.path.join( os.getcwd(), "userSettings.txt")
if not os.path.isfile(userSettingsFilePath) or os.stat(userSettingsFilePath).st_size == 0: us.create(userSettingsFilePath)

userSettingsDictionary = us.load(userSettingsFilePath)
prefix = userSettingsDictionary["prefix"]
loc = userSettingsDictionary["data_location"]
outputloc = userSettingsDictionary["pipeline_output_location"]
if "observables" in userSettingsDictionary:
    observables = userSettingsDictionary["observables"]
else: observables = ['rho', 'Tgas', 'speed']
limits = customRanges(observables)

# Which parts do you want to run?
factor      = 3   # the without inner, is without r< factor * sma
bound       = None
foundModels = searchModels(loc, prefix)
if runModels == '':
    print("The following models within %s "
        "have been found based on the prefix '%s' "%(loc, prefix))
    print('Enter the numbers of the models that you would like to analyse, split multiple models by a space (q or exit to quit, a or all to run all models):')
    for i in range(len(foundModels)):
        if foundModels[i][1] == "": print("\t(%d) %s"%(i, foundModels[i][0]))
        else: print("\t(%d) /...%s"%(i, foundModels[i][1]))

    print()
    runModels = str(input( '  >>>   '))
    runModels = runModels.split()

if any(model in ('q', 'exit') for model in runModels):
    print('')
    print('Program is stopped by the user!!')
    print('')

    print('')
    #print('------------------END:', dt.datetime.now(),'---------------------')
else:
    if any(model in ('a', 'all') for model in runModels):
        models = range(len(foundModels))
    else:
        models = runModels
    if runParts == '':
        print('')
        print('Which components of the pipeline do you want to run?')
        print()
        print('     (1) 2D slice plots of the global structure of the last dump full of the model.')
        print('     (2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.')
        print('     (3) Information about the velocity related quantities of the model + Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.')
        print('     (4) Cummulative mass fraction as a function of the polar coordinate theta.')
        print('     (5) Information of the orbital evolution.')
        print('     (6) 1D spherical profiles for single star models.')
        print()
        print('Choose from 0 to 6, where 0 means \'all\', split multiple components by a space (q or exit to quit):')
        runParts = input('  >>>   ')
        runParts = runParts.split()
    if any(part in ('q', 'exit') for part in runParts):
        print('')
        print('Program is stopped by the user!!')
        print('')

        print('')
        #print('------------------END:', dt.datetime.now(),'---------------------')
    else:
        for i in range(len(runParts)):
            if runParts[i] in ('0', 0):
                for option in options: print(options[option])
            else:
                print(options[runParts[i]])

        print('')
        print('')
        print('It takes some time, so sit back and relax!')
        print('')
        print('')


        try:
            os.makedirs(outputloc)
        except OSError:
            pass

        run_main(outputloc, runParts, models, foundModels)

        print('')
        print('')
        print('--------')
        print('')
        print('')
        print('The pipeline is finished! You can find all plots and reduced data in the following directory:')
        print(str(outputloc))
        print('')

        print('')
        #print('------------------END:', dt.datetime.now(),'---------------------')
