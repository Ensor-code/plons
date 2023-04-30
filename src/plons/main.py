#import datetime                     as dt
from distutils.command.config import dump_file
import numpy                        as np
import sys
import os
import warnings
warnings.filterwarnings("ignore")
import argparse                     as ap
import time

# import plons scripts
import plons.RadialStructurePlots1D       as rs
import plons.SlicePlots2D                 as sl
import plons.CMF_meanRho                  as cmf
import plons.OrbitalEvolution             as ov
import plons.LoadDataPHANTOM              as ld
import plons.TerminalVelocity             as tmv
import plons.Profiles1D                   as dp
import plons.userSettings                 as us
import plons.ArchimedianSpiral            as ars

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
        run = models[int(number)][1].replace('/','',1)
        print('------ MODEL '+str(number)+': '+run+' ------')
        saveloc = os.path.join(outputloc, run)
        print(saveloc, run,loc)
        try:
            os.makedirs(os.path.join(saveloc, 'png'))
            os.makedirs(os.path.join(saveloc, 'txt'))
            os.makedirs(os.path.join(saveloc, 'pdf'))
            os.makedirs(os.path.join(saveloc, 'animation'))
        except OSError:
            pass
        start_time = time.time()
        print('')
        print('Data is loading...')
        [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound, userSettingsDictionary)
        print('All data is loaded and ready to use.')
        print("Loading took %s seconds" % np.round(time.time() - start_time,2))

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
        if zoomin != '':
            for i in range(len(zoomin)): zoomin[i] = int(zoomin[i])
            start_time2 = time.time()
            sl.SlicePlots(run, saveloc, dumpData, setup, zoomin=zoomin, observables=observables)
            print("Sliceplot took %s seconds" % np.round(time.time() - start_time2,2))
        else:
            start_time3 = time.time()
            sl.SlicePlots(run, saveloc, dumpData, setup, observables=observables)
            print("Sliceplot took %s seconds" % np.round(time.time() - start_time3,2))

    if part == '2':
        print('')
        print('(2)  Start calculations for the radial structure plots.')
        start_time4 = time.time()
        rs.radialStructPlots(run, saveloc, dumpData, setup)
        print("Radialstructure took %s seconds" % np.round(time.time() - start_time4,2))
    if part == '3':
        print('')
        print('(3) Start calculations for terminal velocity...')
        start_time5 = time.time()
        tmv.main_terminalVelocity(setup, dumpData, sinkData, saveloc, run)
        print("Terminal velocity took %s seconds" % np.round(time.time() - start_time5,2))
    if part == '4':
        print('')
        print('(4)  Start calculations for the cummulative mass fraction and mean density plots...')
        if setup['single_star'] == True:
            start_timecmf = time.time()
            cmf.CMF_meanRho(run, saveloc, dumpData, setup, factor)
            print("CMF took %s seconds" % np.round(time.time() - start_timecmf,2))
        else:
            start_time6 = time.time()
            cmf.CMF_meanRho(run, saveloc, outerData, setup, factor)
            print("CMF took %s seconds" % np.round(time.time() - start_time6,2))
    if part == '5':
        print('')
        if setup['single_star']:
            print('(5)  A single model has no orbit, and thereby no orbital evolution.')
            print('     The orbital evolution part is therefore skipped.')
        else:
            print('(5)  Start calculations for orbital evolution...')
            start_time7 = time.time()
            ov.orbEv_main(run, saveloc, sinkData, setup)
            print("Orbits took %s seconds" % np.round(time.time() - start_time7,2))

    if part == '6':
        print('')
        print('(6)  Start calculating for the 1D spherical plots')
        start_time8 = time.time()
        dp.profiles_main(run, loc, saveloc, dumpData, setup)
        print("1D profiles took %s seconds" % np.round(time.time() - start_time8,2))

    if part == '7':
        print('')
        print('(7)  Archimedian spiral')
        start_time9 = time.time()
        ars.ArchimedianSpiral(run, saveloc, setup)
        print("Archimedes spiral took %s seconds" % np.round(time.time() - start_time9,2))

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
phantom_dir = userSettingsDictionary["hard_path_to_phantom"]
if "observables" in userSettingsDictionary:
    observables = userSettingsDictionary["observables"]
else: observables = ['rho', 'Tgas', 'speed']

customRanges = True
if customRanges:
    limits = {}
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
        print('     (4) Cummulative mass fraction in function of the polar coordinate theta.')
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
