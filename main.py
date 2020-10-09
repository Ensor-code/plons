import numpy                            as np
import os
import warnings
warnings.filterwarnings("ignore")

# own scripts
import radialStructure_pipeline         as rs
import SlicePlots_pipeline              as sl
import CMF_meanRho_pipeline             as cmf
import orbitalEvolution_pipeline        as ov
import LoadDataPHANTOM                  as d
import TerminalVelocity                 as tmv


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
print('     (4) Quantitative measurement of the degree of aspherical morphology (eta and Qp).')
print('     (5) Cummulative mass fraction in function of the polar coordinate theta.')
print('     (6) Information of the orbital evolution.')
print('')
print('')
print('It takes some time, so sit back and relax!')

# Dit komt nog als inputs
loc         = '/home/silke/Documents/Univ/Master/THESIS/Models/phantom_Masterthesis/desktop_run'
run         = '59'
outputloc   = '/home/silke/Documents/Univ/PhD/Pipeline/testOutput/'
factor      = 3   # the without inner, is without r< factor * sma
bound       = None


try:
    os.mkdir(outputloc)
except OSError:
    print('')


print('')
print(' ---- MODEL '+run+' ----')

print('')
print('Data is loading...')
[setup, dumpData, sinkData, outerData] = d.LoadData_cgs(run, loc, factor, bound)
print('All data is loaded and ready to use.')


# ov.orbEv(run,loc)

sl.SlicePlots(run, outputloc, dumpData, setup)
rs.radialStructPlots(run, outputloc, dumpData, setup)
tmv.main_terminalVelocity(setup, dumpData, outputloc, run)

if setup['single_star'] == True:
	cmf.CMF_meanRho(run, outputloc, dumpData, setup)
else:
	cmf.CMF_meanRho(run, outputloc, outerData, setup)


#print('setup keys:')
## print('')
#for key in setup:
	#print(key)
#print('')

#print('dumpData keys:')
## print('')
#for key in dumpData:
	#print(key)
#print('')


#print('sinkData keys:')
#for key in sinkData:
	#print(key)
#print('')


#print('outerData keys:')
#for key in outerData:
	#print(key)	
#print('')
