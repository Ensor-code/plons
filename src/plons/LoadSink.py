import numpy                    as np
import math                     as math
import os
import sys
import pandas as pd
import numba as nb

import time
# Import plons scripts
import plons.PhysicalQuantities       as pq
import plons.GeometricalFunctions     as gf
import plons.ConversionFactors_cgs    as cgs
import dask.dataframe as dd

unit_factors = {'t': cgs.cu_time(), # evolution time             [yrs]
                'x': cgs.au,  # position coordinates       [cm]
                'y': cgs.au,
                'z': cgs.au,
                'mass': cgs.Msun, # mass of sph particles      [g]
                'vx': cgs.cu_vel(),
                'vy': cgs.cu_vel(),  # velocity components        [cm/s]
                'vz': cgs.cu_vel(),
                'maccr': cgs.Msun}  # accreted mass              [g]

# Pick last wind evolution file
def findLastWindSinkIndex(runName,userPrefix):
    listevFiles = sortedWindSinkList(userPrefix+'Sink0001N0', runName)
    lastFile = listevFiles[-1]
    t1 = lastFile.lstrip(userPrefix+'Sink0001')
    t2 = t1.lstrip("N")
    t3 = t2.rstrip(".ev")
    lastIndex = int(t3)
    return lastIndex
def sortedWindSinkList(userPrefix, runName):
    return np.sort([x for x in os.listdir(runName) if userPrefix in x and ".ev" in x])
def read_data_from_sink_file(filename):
    # Load in the Dask DataFrame
    data = dd.read_table(filename, skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 11), delim_whitespace=True, header=None, names=['t', 'x', 'y', 'z', 'mass', 'vx', 'vy', 'vz', 'maccr'])
    #Apply the unit factors to the relevant columns using `map_partitions()`
    data = data.map_partitions(lambda df: df.assign(**{name: df[name] * factor for name, factor in unit_factors.items()}))
    return data
def append_data(file_lists):
    data_frames = []
    for file_list in file_lists:
        data_frames.append(read_data_from_sink_file(file_list))
    concatenated_data = dd.concat(data_frames)
    numpy_array = concatenated_data.values.compute()
    t, x, y, z, mass, vx, vy, vz, maccr = numpy_array.T
    return t, x, y, z, mass, vx, vy, vz, maccr


'''
Load the .ev-files from a phantom model
      - This file gives the specifics of the sink particles present in the model (AGB star and companion, not of the sph particles)
        in function of the evolution time of the model. The last entry corresponds to the data from the last dump.
      - Only suited for a binary model
      - Units in cgs

INPUT:
    - 'run'   is the number of the run specifically           [str]
    - 'loc'   is the directory where the model is located     [str]
    - 'setup' is the setup data                               [dict]

RETURN:
    a dictionary containing the data from the sink files (all units in cgs)

'''
def LoadSink_cgs(run, loc, setup, userSettingsDictionary):
    runName = os.path.join(loc, run)
    userPrefix = userSettingsDictionary["prefix"]
    # load the dump file wind_00xxx
    t1, x1, y1, z1, mass1, vx1, vy1, vz1, maccr1 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    t2, x2, y2, z2, mass2, vx2, vy2, vz2, maccr2 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    if setup['triple_star']==True:
        t3, x3, y3, z3, mass3, vx3, vy3, vz3, maccr3 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    else:
        setup['triple_star']==False
    #Check quadruple star
    if os.path.exists(os.path.join(runName, str('%sSink0004N0'%userPrefix+str(1)+'.ev'))):
        print('Quadruple setup')
        setup['quadruple_star']=True
        t3, x3, y3, z3, mass3, vx3, vy3, vz3, maccr3 = 0, 0, 0, 0, 0, 0, 0, 0, 0
        t4, x4, y4, z4, mass4, vx4, vy4, vz4, maccr4 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    else:
        setup['quadruple_star']=False
    numberOfevFiles = findLastWindSinkIndex(runName,userPrefix)
    list1, list2, list3, list4 = [], [], [], []

    for n in range(1, numberOfevFiles+1):
        if setup['triple_star']:
            stars = 3
        elif setup['quadruple_star']:
            stars = 4
        else:
            stars = 2

        for i in range(1, stars+1):
            fileName_sink = os.path.join(runName, f'{userPrefix}Sink00{i:02d}N0{n}.ev')

            if i == 1:
                list1.append(fileName_sink)
            elif i == 2:
                list2.append(fileName_sink)
            elif i == 3:
                list3.append(fileName_sink)
            elif i == 4:
                list4.append(fileName_sink)
    t1, x1, y1, z1, mass1, vx1, vy1, vz1, maccr1 = append_data(list1)
    if stars >=2:
        t2, x2, y2, z2, mass2, vx2, vy2, vz2, maccr2 = append_data(list2)
        if stars >=3:
            t3, x3, y3, z3, mass3, vx3, vy3, vz3, maccr3= append_data(list3)
            if stars ==4:
                t4, x4, y4, z4, mass4, vx4, vy4, vz4, maccr4= append_data(list4)
    # AGB star
    r1 = gf.calc_r(x1, y1, z1)                         # [cm]
    position1 = np.column_stack((x1, y1, z1))
    velocity1 = np.column_stack((vx1, vy1, vz1))
    # companion star
    r2 = gf.calc_r(x2, y2, z2)                         # [cm]
    position2 = np.column_stack((x2, y2, z2))
    velocity2 = np.column_stack((vx2, vy2, vz2))
    if setup['triple_star']==True:
        r3 = gf.calc_r(x3, y3, z3)                         # [cm]
        position3 = np.column_stack((x3, y3, z3))
        velocity3 = np.column_stack((vx3, vy3, vz3))
        period_in       = pq.getPeriod(mass1, mass3, (r1 + r3) /cgs.au )
        rHill_in        = pq.getRHill( abs(r1 + r3), mass3, mass1      )              # [cm]
    if setup['quadruple_star']==True:
        r3 = gf.calc_r(x3, y3, z3)                         # [cm]
        position3 = np.column_stack((x3, y3, z3))
        velocity3 = np.column_stack((vx3, vy3, vz3))
        period_in       = pq.getPeriod(mass1, mass3, (r1 + r3) /cgs.au )
        rHill_in        = pq.getRHill( abs(r1 + r3), mass3, mass1      )              # [cm]
        r4 = gf.calc_r(x4, y4, z4)                         # [cm]
        position4 = np.column_stack((x4, y4, z4))
        velocity4 = np.column_stack((vx4, vy4, vz4))
        #deze nog aanpassen?
        period_in2       = pq.getPeriod(mass1, mass4, (r1 + r4) /cgs.au )
        rHill_in2        = pq.getRHill( abs(r1 + r4), mass4, mass1      )              # [cm]
    orbitalVel_AGB = np.sqrt(np.sum(velocity1[:, :2] ** 2, axis=1))
    orbitalVel_comp = np.sqrt(np.sum(velocity2[:, :2] ** 2, axis=1))
    if setup['triple_star']==True:
        orbitalVel_comp_in = np.sqrt(np.transpose(velocity3)[0]**2+np.transpose(velocity3)[1]**2)
    if setup['quadruple_star']==True:
        orbitalVel_comp_in = np.sqrt(np.transpose(velocity3)[0]**2+np.transpose(velocity3)[1]**2)
        orbitalVel_comp_in2 = np.sqrt(np.transpose(velocity4)[0]**2+np.transpose(velocity4)[1]**2)
    rHill           = pq.getRHill( abs(r1 + r2), mass2, mass1           )         # [cm]
    # output
    #    "_t" stands for the fact that these values are in function of the evolution time, not from the last dump in function of location
    if setup['triple_star']==True:
        data = {'posAGB'      : position1,              # [cm]
                'velAGB'      : velocity1,              # [cm/s]
                'massAGB'     : mass1,                  # [gram]
                'maccrAGB'    : maccr1,                 # [gram]
                'rAGB'        : r1,                     # [cm]
                'posComp'     : position2,
                'velComp'     : velocity2,
                'massComp'    : mass2,
                'maccrComp'   : maccr2,
                'rComp'       : r2,
                'time'        : t1,                     # [yrs]
                #'period_t'    : period,                 # [s]
                'v_orbAGB_t'  : orbitalVel_AGB,         # [cm/s]
                'v_orbComp_t' : orbitalVel_comp,        # [cm/s]
                'rHill_t'     : rHill,                  # [cm]
                'posComp_in'  : position3,
                'velComp_in'  : velocity3,
                'massComp_in' : mass3,
                'maccrComp_in': maccr3,
                'rComp_in'    : r3,
                'rHill_in_t'  : rHill_in,
                'period_t_in' : period_in,
                'v_orbComp_in_t': orbitalVel_comp_in
                }
    elif setup['quadruple_star']==True:
        data = {'posAGB'      : position1,              # [cm]
                'velAGB'      : velocity1,              # [cm/s]
                'massAGB'     : mass1,                  # [gram]
                'maccrAGB'    : maccr1,                 # [gram]
                'rAGB'        : r1,                     # [cm]
                'posComp'     : position2,
                'velComp'     : velocity2,
                'massComp'    : mass2,
                'maccrComp'   : maccr2,
                'rComp'       : r2,
                'time'        : t1,                     # [yrs]
                'v_orbAGB_t'  : orbitalVel_AGB,         # [cm/s]
                'v_orbComp_t' : orbitalVel_comp,        # [cm/s]
                'rHill_t'     : rHill,                  # [cm]
                'posComp_in'  : position3,
                'velComp_in'  : velocity3,
                'massComp_in' : mass3,
                'maccrComp_in': maccr3,
                'rComp_in'    : r3,
                'rHill_in_t'  : rHill_in,
                'period_t_in' : period_in,
                'v_orbComp_in_t': orbitalVel_comp_in,
                #Quad Additional
                'posComp_in2'  : position4,
                'velComp_in2'  : velocity4,
                'massComp_in2' : mass4,
                'maccrComp_in2': maccr4,
                'rComp_in2'    : r4,
                'rHill_in_t2'  : rHill_in2,
                'period_t_in2' : period_in2,
                'v_orbComp_in_t2': orbitalVel_comp_in2
                }

    else: #binary
        data = {'posAGB'      : position1,              # [cm]
                'velAGB'      : velocity1,              # [cm/s]
                'massAGB'     : mass1,                  # [gram]
                'maccrAGB'    : maccr1,                 # [gram]
                'rAGB'        : r1,                     # [cm]
                'posComp'     : position2,
                'velComp'     : velocity2,
                'massComp'    : mass2,
                'maccrComp'   : maccr2,
                'rComp'       : r2,
                'time'        : t1,                     # [yrs]
                #'period_t'    : period,                 # [s]
                'v_orbAGB_t'  : orbitalVel_AGB,         # [cm/s]
                'v_orbComp_t' : orbitalVel_comp,        # [cm/s]
                'rHill_t'     : rHill                  # [cm]
                }


    return data






'''
Load the .ev-files from a phantom model of a SINGLE star
      - This file gives the specifics of the sink particle present in the model (AGB star, not of the sph particles)
        in function of the evolution time of the model. The last entry corresponds to the data from the last dump.
      - Only suited for a single model
      - Units in cgs
'''
def LoadSink_single_cgs(run, loc, setup, userSettingsDictionary):

    runName = os.path.join(loc,run)
    userPrefix = userSettingsDictionary["prefix"]
    fileName = os.path.join(runName, '%sSink0001N01.ev'%userPrefix)

    # load the dump file prefix_00xxx
    t1, x1, y1, z1, mass1, vx1, vy1, vz1, maccr1 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    try:
       (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(fileName, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)


    except OSError:

        print(' ERROR: No sink file found for this model in the current directory!')
        sys.exit()

    # AGB star

    t1     = t1     *  cgs.cu_time()                   # evolution time             [yrs]
    x1     = x1     *  cgs.au                          # position coordinates       [cm]
    y1     = y1     *  cgs.au
    z1     = z1     *  cgs.au
    mass1  = mass1  *  cgs.Msun                        # mass of sph particles      [g]
    vx1    = vx1    *  cgs.cu_vel()                    # velocity components        [cm/s]
    vy1    = vy1    *  cgs.cu_vel()
    vz1    = vz1    *  cgs.cu_vel()
    maccr1 = maccr1 *  cgs.Msun                        # accreted mass              [g]

    r1 = gf.calc_r(x1, y1, z1)                         # [cm]
    position1 = np.array((x1, y1, z1 )).transpose()
    velocity1 = np.array((vx1,vy1,vz1)).transpose()



    data = {'posAGB'      : position1,       # [cm]
            'velAGB'      : velocity1,       # [cm/s]
            'massAGB'     : mass1,           # [gram]
            'maccrAGB'    : maccr1,          # [gram]
            'rAGB'        : r1,              # [cm]
            'time'        : t1,              # [yrs]
            }

    return data
