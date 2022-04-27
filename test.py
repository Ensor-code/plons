import numpy                    as np
import math                     as math
import os
import ConversionFactors_cgs    as cgs
import PhysicalQuantities       as pq
import LoadDataPHANTOM          as ld
#import LoadSink                 as snk
import TerminalVelocity         as tmv
import sys
sys.path.append('/home/mats/codes/phantom/scripts')
sys.path.append('/home/matse/codes/phantom/scripts')
from readPhantomDump import *

dump = read_dump('../rayTracer/wind_00010')

# print(dump["quantities"]["hfact"])
for i in dump["blocks"][1]["data"]["m"]: print(i)


#print(pq.getTemp(10e10,10e3,1.43))
#print(pq.getPeriod(1.5*cgs.Msun_gram(),1*cgs.Msun_gram(),3)*cgs.sec_year())
#print(pq.getOrbitalVelocity(3/cgs.sec_year(),3)*cgs.cms_kms())
#print('')


# loc = '/home/silke/Documents/Univ/Master/THESIS/Models/phantom_Masterthesis/desktop_run'
# run = '39'
# outputloc = '/home/silke/Documents/Univ/PhD/Pipeline/testOutput/'

'''
Test load data
'''
#print('Binary star')
#print('')
#print('data from model '+run+' is loading.')
# print('Load data')
# print('')
# setup, dump, sink, outer = ld.LoadData_cgs(run,loc, 5, 200)
#print('loading model complete.')
#print('')
#print(setup['sma_ini'])
#print(setup['period_ini']*cgs.sec_year())
#print(dump['period_fi' ]*cgs.sec_year())
#print(dump['rho'][1000])
#print(dump['posAGB'])
#print(sink['posAGB'][-1])
#print(outer['rho'][1000])
##sinkData = snk.LoadSink_cgs(run, loc)
##print(sinkData['posAGB'])
##print(sinkData['posAGB'][1])

#print('')
#print('')
#print('Single star')
#print('')

#loc = '/home/silke/Documents/Univ/Master/THESIS/Models/phantom_Masterthesis/desktop_run'
#run = 'O'
#print('data from model '+run+' is loading.')
#setup, dump, sink, outer = ld.LoadData_cgs(run,loc, 5, 200)
#print('loading model complete.')
#print('')
#print(setup['v_ini'])
##print(dump['period_fi' ]*cgs.sec_year())
#print(dump['rho'][-1])
#print(dump['posAGB'])
#print(sink['posAGB'][-1])

# '''
# Test terminal velocities
# '''

# tmv.main_terminalVelocity(setup, dump, outputloc, run)

