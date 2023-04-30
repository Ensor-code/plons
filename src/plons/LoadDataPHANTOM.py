import math                     as math
import time
import numpy as np
#import plons scripts
import plons.LoadDump                 as dmp
import plons.LoadSink                 as snk
import plons.LoadSetup                as stp


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
    start_time = time.time()
    setup       = stp.LoadSetup(run, loc, userSettingsDictionary)
    print("Loading setup took %s seconds" % np.round(time.time() - start_time,2))
    start_time = time.time()
    dumpData  = dmp.LoadDump_cgs(run, loc, setup, userSettingsDictionary, number)
    print("Loading dumpData took %s seconds" % np.round(time.time() - start_time,2))
    if setup['single_star']:
        start_timea = time.time()
        sinkData  = snk.LoadSink_single_cgs(run, loc, setup, userSettingsDictionary)
        outerData = None
        print("Loading sinkData took %s seconds" % np.round(time.time() - start_timea,2))
    else:
        start_timeb = time.time()
        sinkData  = snk.LoadSink_cgs(run, loc, setup, userSettingsDictionary)
        print("Loading sinkData took %s seconds" % np.round(time.time() - start_timeb,2))

        if bound == None:
            bound = setup['bound']
        start_time = time.time()
        outerData = dmp.LoadDump_outer_cgs(factor, bound, setup, dumpData)
        print("Loading outerData took %s seconds" % np.round(time.time() - start_time,1))

    # save the final specifics of the AGB star to dumpData
    dumpData['posAGB'   ] = sinkData['posAGB'     ][-1]
    dumpData['velAGB'   ] = sinkData['velAGB'     ][-1]
    dumpData['rAGB'     ] = sinkData['rAGB'       ][-1]
    dumpData['massAGB'  ] = sinkData['massAGB'    ][-1]
    dumpData['maccrAGB' ] = sinkData['maccrAGB'   ][-1]

    # save the final specifics of the companion to dumpData
    if not setup['single_star']:
        dumpData['velComp'  ] = sinkData['velComp'    ][-1]
        dumpData['maccrComp'] = sinkData['maccrComp'  ][-1]
        # dumpData['period_fi'] = sinkData['period_t'   ][-1]
        dumpData['v_orbAGB' ] = sinkData['v_orbAGB_t' ][-1]
        dumpData['v_orbComp'] = sinkData['v_orbComp_t'][-1]
        dumpData['sma_fi'   ] = dumpData['rAGB'] + dumpData['rComp']        # [cm]
        #dumpData['v_orb_fi' ] = pq.getOrbitalVelocity(dumpData['period_fi'], dumpData['sma_fi'] /cgs.au )

        if setup['triple_star']:
            dumpData['velComp_in'  ] = sinkData['velComp_in'    ][-1]
            dumpData['maccrComp_in'] = sinkData['maccrComp_in'  ][-1]
            # dumpData['period_fi_in'] = sinkData['period_t_in'   ][-1]
            dumpData['v_orbComp_in'] = sinkData['v_orbComp_in_t'][-1]
        #hier nog quad doen
        if setup['quadruple_star']:
            dumpData['velComp_in'  ] = sinkData['velComp_in'    ][-1]
            dumpData['maccrComp_in'] = sinkData['maccrComp_in'  ][-1]
            dumpData['v_orbComp_in'] = sinkData['v_orbComp_in_t'][-1]
            #dumpData['posComp_in2'  ] = sinkData['posComp_in2'    ][-1]
            #dumpData['posComp_in'  ] = sinkData['posComp_in'    ][-1]
    return setup, dumpData, sinkData, outerData
