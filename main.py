import datetime                     as dt
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


print('------------------START:', dt.datetime.now(),'---------------------')
print('')

options = { '0': '(1) 2D slice plots \n(2) 1D line/tube plots \n(3) Terminal velocity \n(4) Morphological parameters\n(5) Cummulative mass fraction\n(6) Orbital evolution\n(7) Tube plots ', 
            '1': '(1) 2D slice plots', 
            '2': '(2) 1D line plots & tube plots',
            '3': '(3) velocity related quantities',
            '4': '(4) Morphological parameters',
            '5': '(5) Cummulative mass fraction',
            '6': '(6) Orbital evolution',
            '7': '(7) Tube plots'
            }


def run_main(outputloc,part,Numbers):
    for run in Numbers:
        print('---- MODEL '+run+' ----')
        saveloc = outputloc+run+'/' 
        
        try:
            os.mkdir(outputloc+run+'/')
        except OSError:
            print('')
            
        print('')
        print('Data is loading...')
        [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound)
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
                    cmf.CMF_meanRho(run, saveloc, dumpData, setup)
                else:
                    cmf.CMF_meanRho(run, saveloc, outerData, setup)
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
                    cmf.CMF_meanRho(run, saveloc, dumpData, setup)
                else:
                    cmf.CMF_meanRho(run, saveloc, outerData, setup)
                    
            if part == '6':
                # (6) orbital evolution
                ov.orbEv_main(run, saveloc, sinkData, setup)
                
            if part == '7':
                # (7) tube plots
                tb.main_tube(run, saveloc, setup, dumpData)
        print('')
        
        
    

    

print('')
print('-------------------------------------------------------')
print('     Welcome to the very first PHANTOM pipeline!'        )
print('-------------------------------------------------------')
print('')
print('This pipeline reduces PHANTOM output to usable plots and datasets.')
print('It returns:')
print('     (1) 2D slice plots of the global structure of the last dump of the model.')
print('     (2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.')
print('     (3) Information about the velocity related quantities of the model.')
print('     (4) Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.')
print('     (5) Cummulative mass fraction in function of the polar coordinate theta.')
print('     (6) Information of the orbital evolution.')
print('     (7) Tube plots for classifying EDEs/flattenings.')
print('')
print('')


loc         = None
outputloc   = None
factor      = 3   # the without inner, is without r< factor * sma
bound       = None


if loc == None or outputloc == None:
    print(' !! Before you start, give the directory where: ')
    print('     - the PHANTOM output of the models is saved   (LoadDataPHANTOM.py, line 109)')
    print('     - the output of the pipeline should be saved  (LoadDataPHANTOM.py, line 110)')
    print('')
    print('------------------END:', dt.datetime.now(),'---------------------')
    sys.exit()


# Which parts do you want to run?

print('Which models do you want to run this for? ')
print('Give the number, split multiple models by a space (q to quit):')
runNumbers = str(input( '  >>>   '))
if runNumbers == 'q':
    print('')
    print('Program is stopped by the user!!')
    print('')
    
    print('')
    print('------------------END:', dt.datetime.now(),'---------------------')
else:
    numbers    = runNumbers.split()
    print('')
    print('Which components of the pipeline do you want to run?')
    print('Choose from 0 to 6, where 0 means \'all\', split multiple components by a space (q to quit):')
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
            os.mkdir(outputloc)
        except OSError:
            print('')


        run_main(outputloc, runParts, numbers)

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
