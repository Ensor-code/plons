import datetime                     as dt
from distutils.command.config import dump_file
import numpy                        as np
import sys
import os
import warnings
warnings.filterwarnings("ignore")

# own scripts
import RadialStructurePlots1D       as rs
import SlicePlots2D                 as sl
import CMF_meanRho                  as cmf
import OrbitalEvolution             as ov
import LoadDataPHANTOM              as ld
import TerminalVelocity             as tmv
import Tubes                        as tb
import userSettings                 as us
import LoadDump                     as lodu

print('------------------START:', dt.datetime.now(),'---------------------')
print('')

options = { '0': '(1) 2D slice plots \n(2) 1D line/tube plots \n(3) Terminal velocity \n(4) Morphological parameters\n(5) Cummulative mass fraction\n(6) Orbital evolution\n(7) Tube plots ', 
            '1': '(1) 2D slice plots', 
            '2': '(2) 1D line plots & tube plots',
            '3': '(3) velocity related quantities',
            '4': '(4) Morphological parameters',
            '5': '(5) Cummulative mass fraction',
            '6': '(6) Orbital evolution',
            '7': '(7) Tube plots',
            '8': '(8) Animation 2D slice plots'
            }


def run_main(outputloc,runParts,numbers, models):
    for number in numbers:
        run = models[int(number)][1][1:]
        print('---- MODEL '+run+' ----')
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
                # (1) 2D slice plots
                sl.SlicePlots(run, saveloc, dumpData, setup)
                # (2) 1D line plots
                rs.radialStructPlots(run, saveloc, dumpData, setup)
                # (3) and (4) terminal velocity, eta, Qp
                tmv.main_terminalVelocity(setup, dumpData, sinkData, saveloc, run)
                # (5) cummulative mass fraction
                if setup['single_star'] == True:
                    cmf.CMF_meanRho(run, saveloc, dumpData, setup, factor)
                else:
                    cmf.CMF_meanRho(run, saveloc, outerData, setup, factor)
                # (6) orbital evolution
                ov.orbEv_main(run, saveloc, sinkData, setup)
                # (7) tube plots
                tb.main_tube(run, saveloc, setup, dumpData)
                
            if part == '1':
                # (1) 2D slice plots
                sl.SlicePlots(run, saveloc, dumpData, setup)
                
            if part == '2':
                # (2) 1D line plots
                rs.radialStructPlots(run, saveloc, dumpData, setup)

            if part == '3' or part == '4':  
                # (3) and (4) terminal velocity, eta, Qp
                tmv.main_terminalVelocity(setup, dumpData, sinkData, saveloc, run)
                
            if part == '5':
                # (5) cummulative mass fraction
                if setup['single_star'] == True:
                    cmf.CMF_meanRho(run, saveloc, dumpData, setup, factor)
                else:
                    cmf.CMF_meanRho(run, saveloc, outerData, setup, factor)
                    
            if part == '6':
                # (6) orbital evolution
                ov.orbEv_main(run, saveloc, sinkData, setup)
                
            if part == '7':
                # (7) tube plots
                tb.main_tube(run, saveloc, setup, dumpData)
        print('')

def searchModels(loc, prefix):
    result = []
    for path, directories, files in os.walk(loc):
        dumpFiles = list(filter(lambda x: prefix+'_' in x, files))
        if len(dumpFiles) != 0:
            slicedString = path.replace(loc, "")
            result.append([path, slicedString])

    return result

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


# Which parts do you want to run?
factor      = 3   # the without inner, is without r< factor * sma
bound       = None
foundModels = searchModels(loc, prefix)
print("The following models within %s "
      "have been found based on the prefix '%s' "%(loc, prefix))
print('Enter the numbers of the models that you would like to analyse, split multiple models by a space (q to quit):')
for i in range(len(foundModels)):
    if foundModels[i][1] == "": print("\t(%d) %s"%(i, foundModels[i][0]))
    else: print("\t(%d) /...%s"%(i, foundModels[i][1]))

print()
runModels = str(input( '  >>>   '))
if runModels == 'q':
    print('')
    print('Program is stopped by the user!!')
    print('')

    print('')
    print('------------------END:', dt.datetime.now(),'---------------------')
else:
    models    = runModels.split()
    print('')
    print('Which components of the pipeline do you want to run?')
    print()
    print('     (1) 2D slice plots of the global structure of the last dump full of the model.')
    print('     (2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.')
    print('     (3) Information about the velocity related quantities of the model.')
    print('     (4) Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.')
    print('     (5) Cummulative mass fraction in function of the polar coordinate theta.')
    print('     (6) Information of the orbital evolution.')
    print('     (7) Tube plots for classifying EDEs/flattenings.')
    print()
    print('Choose from 0 to 7, where 0 means \'all\', split multiple components by a space (q to quit):')
    part = input('  >>>   ')
    if part == 'q':
        print('')
        print('Program is stopped by the user!!')
        print('')

        print('')
        print('------------------END:', dt.datetime.now(),'---------------------')
    else:
        runParts = part.split()
        for i in range(len(runParts)):
            print(options[runParts[i]])

        if '3' in runParts and '4' in runParts:
            runParts.remove('4')

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
        print('------------------END:', dt.datetime.now(),'---------------------')
