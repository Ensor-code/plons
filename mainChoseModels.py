import datetime                     as dt
import numpy                        as np
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


print('------------------START:', dt.datetime.now(),'---------------------')
print('')

options = { '0': '(1) 2D slice plots \n(2) 1D line plots \n(3) Terminal velocity \n(4) Morphological parameters\n(5) Cummulative mass fraction\n(6) Orbital evolution ', 
            '1': '(1) 2D slice plots', 
            '2': '(2) 1D line plots',
            '3': '(3) Terminal velocity',
            '4': '(4) Morphological parameters',
            '5': '(5) Cummulative mass fraction',
            '6': '(6) Orbital evolution'
            }


def run_main(part,Numbers):
    for run in Numbers:
        print('---- MODEL '+run+' ----')

        print('')
        print('Data is loading...')
        [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound)
        print('All data is loaded and ready to use.')
        print('')

        if part == '0':
            # (1) 2D slice plots
            sl.SlicePlots(run, outputloc, dumpData, setup)
            # (2) 1D line plots
            rs.radialStructPlots(run, outputloc, dumpData, setup)
            # (3) and (4) terminal velocity, eta, Qp
            tmv.main_terminalVelocity(setup, dumpData, outputloc, run)
            # (5) cummulative mass fraction
            if setup['single_star'] == True:
                cmf.CMF_meanRho(run, outputloc, dumpData, setup)
            else:
                cmf.CMF_meanRho(run, outputloc, outerData, setup)
            # (6) orbital evolution
            ov.orbEv_main(run, outputloc, sinkData, setup)
            
        if part == '1':
            # (1) 2D slice plots
            sl.SlicePlots(run, outputloc, dumpData, setup)
            
        if part == '2':
            # (2) 1D line plots
            rs.radialStructPlots(run, outputloc, dumpData, setup)

        if part == '3' or part == '4':  
            # (3) and (4) terminal velocity, eta, Qp
            tmv.main_terminalVelocity(setup, dumpData, sinkData, outputloc, run)
            
        if part == '5':
            # (5) cummulative mass fraction
            if setup['single_star'] == True:
                cmf.CMF_meanRho(run, outputloc, dumpData, setup)
            else:
                cmf.CMF_meanRho(run, outputloc, outerData, setup)
                
        if part == '6':
            # (6) orbital evolution
            ov.orbEv_main(run, outputloc, sinkData, setup)
        
        

print('')
print('-------------------------------------------------------')
print('     Welcome to the very first PHANTOM pipeline!'        )
print('-------------------------------------------------------')
print('')
print('This pipeline reduces PHANTOM output to usable plots and datasets.')
print('It returns:')
print('     (1) 2D slice plots of the global structure of the last dump of the model.')
print('     (2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.')
print('     (3) Information about the terminal velocity of the model.')
print('     (4) Quantitative measurement of the degree of aspherical morphology: morphological parameters eta and Qp.')
print('     (5) Cummulative mass fraction in function of the polar coordinate theta.')
print('     (6) Information of the orbital evolution.')

# which parts do you want to run?
print('')
print('')
print('Which components of the pipeline do you want to run?')
print('Choose from 0 to 6, where 0 means \'all\', split multiple components by a space:')
part = input('  >>>   ')
runParts = part.split() 
print('For which models do you want to run this? ')
print('Give the modelnames, split multiple models by a space:')
runNumbers = str(input( '  >>>   '))
Numbers    = runNumbers.split()

for i in range(len(runParts)):
    print(options[runParts[i]])
    
if '3' in runParts and '4' in runParts:
    runParts.remove('4')


#print(options[part])
# # Dit komt nog als inputs
# run       = 'M16A'
loc       = '/home/user/Documents/thesis/modellen_pc/'
outputloc = '/home/user/Documents/phantom pipeline/'
factor    = 10   # the without inner, is without r< factor * sma
bound     = None
    
print('')
print('')
print('It takes some time, so sit back and relax!')




try:
    os.mkdir(outputloc)
except OSError:
    print('')


print('')


for i in range(len(runParts)):
    run_main(runParts[i],Numbers)

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
