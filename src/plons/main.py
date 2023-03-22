#import datetime                     as dt
from distutils.command.config import dump_file
import numpy                        as np
import sys
import os
import warnings
warnings.filterwarnings("ignore")
import argparse                     as ap

# own scripts
import plons.RadialStructurePlots1D       as rs
import plons.SlicePlots2D                 as sl
import plons.CMF_meanRho                  as cmf
import plons.OrbitalEvolution             as ov
import plons.LoadDataPHANTOM              as ld
import plons.TerminalVelocity             as tmv
import plons.Tubes                        as tb
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
p.add_argument("-z", "--zoom", action="store", default='', type = str, nargs='*', dest = "zoom", help = helpoptions)
args = p.parse_args()
runModels = args.models
runParts = args.options
zoomin = args.zoom

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
        [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound, userSettingsDictionary)
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
        if zoomin != '':
            for i in range(len(zoomin)): zoomin[i] = int(zoomin[i])
            sl.SlicePlots(run, saveloc, dumpData, setup, zoomin=zoomin, observables=observables)
        else:
            sl.SlicePlots(run, saveloc, dumpData, setup, observables=observables)
        
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

print('')
print('-------------------------------------------------------')
print('     Welcome to the very first PHANTOM pipeline!'        )
print('-------------------------------------------------------')
print('')
print('This pipeline reduces PHANTOM output to usable plots and datasets.')
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
