import numpy                    as np
import math                     as math
import os
# own scripts
import PhysicalQuantities       as pq
import GeometricalFunctions     as gf
import ConversionFactors_cgs    as cgs
import LoadSetup                as stp
import sys




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

    numberOfevFiles = findLastWindSinkIndex(runName)

    for n in range(1,numberOfevFiles+1):
        #print(n)
        fileName_sink1 = os.path.join(runName, str('%sSink0001N0'%userPrefix+str(n)+'.ev'))
        fileName_sink2 = os.path.join(runName, str('%sSink0002N0'%userPrefix+str(n)+'.ev'))

        if setup['triple_star']==True:
            fileName_sink3 = os.path.join(runName, str('%sSink0003N0'%userPrefix+str(n)+'.ev'))

        try:
        # to calculate period, we need masses and sma, so coordinates
            (t1e, x1e,y1e,z1e, mass1e, vx1e,vy1e,vz1e, maccr1e) = np.loadtxt(fileName_sink1, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
            n_file = len(t1e)
            (t2e, x2e,y2e,z2e, mass2e, vx2e,vy2e,vz2e, maccr2e) = np.loadtxt(fileName_sink2, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)[:, :n_file]
            if setup['triple_star']==True:
                (t3e, x3e,y3e,z3e, mass3e, vx3e,vy3e,vz3e, maccr3e) = np.loadtxt(fileName_sink3, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)[:, :n_file]
        except OSError:
            print(' ERROR: No sink files found for this model in the current directory!')
            
        t1     = np.append(t1,t1e)
        x1     = np.append(x1,x1e)
        y1     = np.append(y1,y1e)
        z1     = np.append(z1,z1e)
        mass1  = np.append(mass1, mass1e)
        vx1    = np.append(vx1, vx1e)
        vy1    = np.append(vy1, vy1e)
        vz1    = np.append(vz1, vz1e)
        maccr1 = np.append(maccr1, maccr1e)

        t2     = np.append(t2,t2e)
        x2     = np.append(x2,x2e)
        y2     = np.append(y2,y2e)
        z2     = np.append(z2,z2e)
        mass2  = np.append(mass2, mass2e)
        vx2    = np.append(vx2, vx2e)
        vy2    = np.append(vy2, vy2e)
        vz2    = np.append(vz2, vz2e)
        maccr2 = np.append(maccr2, maccr2e)
        
        if setup['triple_star']==True:
            t3     = np.append(t3,t3e)
            x3     = np.append(x3,x3e)
            y3     = np.append(y3,y3e)
            z3     = np.append(z3,z3e)
            mass3  = np.append(mass3, mass3e)
            vx3    = np.append(vx3, vx3e)
            vy3    = np.append(vy3, vy3e)
            vz3    = np.append(vz3, vz3e)
            maccr3 = np.append(maccr3, maccr3e)

        


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
    
    # companion star

    t2     = t2     *  cgs.cu_time()                   # evolution time             [yrs]                                 
    x2     = x2     *  cgs.au                          # position coordinates       [cm]
    y2     = y2     *  cgs.au
    z2     = z2     *  cgs.au
    mass2  = mass2  *  cgs.Msun                        # mass of sph particles      [g]
    vx2    = vx2    *  cgs.cu_vel()                    # velocity components        [cm/s]
    vy2    = vy2    *  cgs.cu_vel()
    vz2    = vz2    *  cgs.cu_vel()
    maccr2 = maccr2 *  cgs.Msun                        # accreted mass              [g]

    r2 = gf.calc_r(x2, y2, z2)                         # [cm]

    position2 = np.array((x2, y2, z2 )).transpose()
    velocity2 = np.array((vx2,vy2,vz2)).transpose()

    if setup['triple_star']==True:
        #close companion star
        t3     = t3     *  cgs.cu_time()                   # evolution time             [yrs]                                 
        x3     = x3     *  cgs.au                          # position coordinates       [cm]
        y3     = y3     *  cgs.au
        z3     = z3     *  cgs.au
        mass3  = mass3  *  cgs.Msun                        # mass of sph particles      [g]
        vx3    = vx3    *  cgs.cu_vel()                    # velocity components        [cm/s]
        vy3    = vy3    *  cgs.cu_vel()
        vz3    = vz3    *  cgs.cu_vel()
        maccr3 = maccr3 *  cgs.Msun                        # accreted mass              [g]

        r3 = gf.calc_r(x3, y3, z3)                         # [cm]

        position3 = np.array((x3, y3, z3 )).transpose()
        velocity3 = np.array((vx3,vy3,vz3)).transpose()
        
        period_in       = pq.getPeriod(mass1, mass3, (r1 + r3) /cgs.au )
        rHill_in        = pq.getRHill( abs(r1 + r3), mass3, mass1      )              # [cm]

    
    
    
    # orbital information 
    # NOT CORRECT!!!
    #period          = pq.getPeriod(mass1, mass2, setup['sma_ini'] )
    
    #print('period 1',period)
    #print('period 2',setup['period'])
    #periodFixed = setup['period']
    #ONLY ORBITAL VEL OF OUTER COMPANION WILL BE CORRECT IN CASE OF TRIPLE, INNER HAS COMPLICATED ORBITAL VELOCITY
    #orbitalVel_AGB  = pq.getOrbitalVelocity(period, r1     /cgs.au_cm() )
    #orbitalVel_comp = pq.getOrbitalVelocity(period, r2     /cgs.au_cm() )
    #orbotalVel_comp_in = pq.getOrbitalVelocity(period_in, r3     /cgs.au_cm() )
    
    orbitalVel_AGB  = np.sqrt(np.transpose(velocity1)[0]**2+np.transpose(velocity1)[1]**2) 
    orbitalVel_comp = np.sqrt(np.transpose(velocity2)[0]**2+np.transpose(velocity2)[1]**2)
    if setup['triple_star']==True:
        orbitalVel_comp_in = np.sqrt(np.transpose(velocity3)[0]**2+np.transpose(velocity3)[1]**2)
        
        
        
    rHill           = pq.getRHill( abs(r1 + r2), mass2, mass1           )         # [cm]
    
    #print('test orbital velocity:')
    #print(np.shape(np.linalg.norm(np.transpose(velocity1))),np.shape(orbitalVel_AGB))
    #print(np.nanmean(np.linalg.norm(np.transpose(velocity1))),np.nanmean(orbitalVel_AGB))
    #print(np.array(np.linalg.norm(velocity1))[0:20],orbitalVel_AGB[0:20])
    
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
    else:
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
        
    #print(np.shape(data['velAGB']))
    #print(np.shape(data['velAGB'][0]))
    #print(np.shape(np.transpose(data['velAGB'])[0]))

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



# Pick last wind evolution file
def findLastWindSinkIndex(runName):
    listevFiles = sortedWindSinkList('windSink0001N0', runName)
    lastFile = listevFiles[-1]
    t1 = lastFile.lstrip("windSink0001")
    t2 = t1.lstrip("N")
    t3 = t2.rstrip(".ev")
    lastIndex = int(t3)
    return lastIndex

def sortedWindSinkList(userPrefix, runName):
    return np.sort(list(filter(lambda x: ("%s"%userPrefix in x) and (".ev" in x), os.listdir(runName))))

