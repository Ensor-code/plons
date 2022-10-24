import numpy                    as np
import math                     as math
import os
# own scripts
import PhysicalQuantities       as pq
import ConversionFactors_cgs    as cgs
import LoadDump                 as dmp
import LoadSink                 as snk
import LoadSetup                as stp


'''
Load the data for a general model
    * INPUT
        - 'run'     [str]   number of the model 
        - 'loc'     [str]   directory of the output data from PHANTOM 
        - 'factor'  [float] factor of the semi-major axis that you can to discard in order to consider the selfsimilar part of the model 
        - 'bound'   [float] outer boundary you want to consider [AU] if you want the setup outer boundary, use None
    * RETURN: four different dictionaries for a binary/single model containing:           
        - 'setup'           setup information of the model
        - 'dumpData'        data from the last full dump (last time step of the model)
        - 'sinkData'        data of the two sink particles in function of time
        - 'outerData'       data from the last dump in a chosen range (None for a single model)
'''
def LoadData_cgs(run, loc, factor, bound, userSettingsDictionary, number = -1):
    
    setup       = stp.LoadSetup(run, loc, userSettingsDictionary)
    single_star = setup['single_star']
    
    if single_star == False:
        dumpData, sinkData, outerData = LoadData_binary_cgs(run, loc, factor, bound, setup, userSettingsDictionary, number)
        
        return setup, dumpData, sinkData, outerData
        
    if single_star == True:
        dumpData, sinkData = LoadData_single_cgs(run, loc, setup, userSettingsDictionary)
        
        return setup, dumpData, sinkData, None
    
    

'''
Load the data for a binary star model
'''
def LoadData_binary_cgs(run, loc, factor, bound, setup, userSettingsDictionary, number):
    
    dumpData  = dmp.LoadDump_cgs(run, loc, setup, userSettingsDictionary, number)
    sinkData  = snk.LoadSink_cgs(run, loc, setup, userSettingsDictionary)
    if bound == None:
        bound = setup['bound']
    outerData = dmp.LoadDump_outer_cgs(run, loc, factor, bound, setup, dumpData)
    
    # save the final specifics of the AGB star and companion to dumpData
    #dumpData['posAGB'   ] = sinkData['posAGB'     ][-1]
    #dumpData['posComp'  ] = sinkData['posComp'    ][-1]
                                                  
    dumpData['velAGB'   ] = sinkData['velAGB'     ][-1]
    dumpData['velComp'  ] = sinkData['velComp'    ][-1]
    
    #dumpData['rAGB'     ] = sinkData['rAGB'       ][-1]
    #dumpData['rComp'    ] = sinkData['rComp'      ][-1]
                                                  
    #dumpData['massAGB'  ] = sinkData['massAGB'    ][-1]
    #dumpData['massComp' ] = sinkData['massComp'   ][-1]
                                                  
    dumpData['maccrAGB' ] = sinkData['maccrAGB'   ][-1]
    dumpData['maccrComp'] = sinkData['maccrComp'  ][-1]
                                                  
    dumpData['period_fi'] = sinkData['period_t'   ][-1]
                                                  
    #dumpData['rHill'    ] = sinkData['rHill_t'    ][-1]
                                                  
    dumpData['v_orbAGB' ] = sinkData['v_orbAGB_t' ][-1]
    dumpData['v_orbComp'] = sinkData['v_orbComp_t'][-1]
    
    dumpData['sma_fi'   ] = dumpData['rAGB'] + dumpData['rComp']        # [cm]
    dumpData['v_orb_fi' ] = pq.getOrbitalVelocity(dumpData['period_fi'], dumpData['sma_fi'] /cgs.AU_cm() )

    if setup['triple_star']==True:
        dumpData['velComp_in'  ] = sinkData['velComp_in'    ][-1]
        dumpData['maccrComp_in'] = sinkData['maccrComp_in'  ][-1]
        dumpData['period_fi_in'] = sinkData['period_t_in'   ][-1]
        dumpData['v_orbComp_in'] = sinkData['v_orbComp_in_t'][-1]


    return dumpData, sinkData, outerData





'''
Load the data for a single star model
'''
def LoadData_single_cgs(run, loc, setup, userSettingsDictionary):
    
    dumpData = dmp.LoadDump_single_cgs(run, loc, setup, userSettingsDictionary)
    sinkData = snk.LoadSink_single_cgs(run, loc, setup, userSettingsDictionary)
    
    # save the final specifics of the AGB star and companion to dumpData
    dumpData['posAGB'   ] = sinkData['posAGB'     ][-1]
    dumpData['velAGB'   ] = sinkData['velAGB'     ][-1]
    dumpData['rAGB'     ] = sinkData['rAGB'       ][-1]
    dumpData['massAGB'  ] = sinkData['massAGB'    ][-1]
    dumpData['maccrAGB' ] = sinkData['maccrAGB'   ][-1]
   
    return dumpData, sinkData



