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

    fileName_sink11 = os.path.join(runName, '%sSink0001N01.ev'%userPrefix)
    fileName_sink21 = os.path.join(runName, '%sSink0002N01.ev'%userPrefix)

    # load the dump file wind_00xxx
    t1, x1, y1, z1, mass1, vx1, vy1, vz1, maccr1 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    t2, x2, y2, z2, mass2, vx2, vy2, vz2, maccr2 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    try:
    # to calculate period, we need masses and sma, so coordinates, we call the parameters xI, I stands for input
        (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(fileName_sink11, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        (t2, x2,y2,z2, mass2, vx2,vy2,vz2, maccr2) = np.loadtxt(fileName_sink21, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        
    except OSError:

        print(' ERROR: No sink files found for this model in the current directory!')
        
    
    '''
    Uncommand this part of the code when your simulation has been paused and restarted.
    Give the model's name.
    '''
    fileName_sink12 = os.path.join(runName, '%sSink0001N02.ev'%userPrefix)
    fileName_sink22 = os.path.join(runName, '%sSink0002N02.ev'%userPrefix)
    t1e, x1e, y1e, z1e, mass1e, vx1e, vy1e, vz1e, maccr1e = 0, 0, 0, 0, 0, 0, 0, 0, 0
    t2e, x2e, y2e, z2e, mass2e, vx2e, vy2e, vz2e, maccr2e = 0, 0, 0, 0, 0, 0, 0, 0, 0
    if run == '59' or run == '52' or run == '48' or run == '41':

        try:
            (t1e, x1e,y1e,z1e, mass1e, vx1e,vy1e,vz1e, maccr1e) = np.loadtxt(fileName_sink12, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
            (t2e, x2e,y2e,z2e, mass2e, vx2e,vy2e,vz2e, maccr2e) = np.loadtxt(fileName_sink22, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        except OSError:

            print(' ERROR: No extra sink files found for this model.')


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


    # AGB star

    t1     = t1     *  cgs.cu_time()                   # evolution time             [yrs]                                
    x1     = x1     *  cgs.AU_cm()                     # position coordinates       [cm]
    y1     = y1     *  cgs.AU_cm()       
    z1     = z1     *  cgs.AU_cm()      
    mass1  = mass1  *  cgs.Msun_gram()                 # mass of sph particles      [g]
    vx1    = vx1    *  cgs.cu_vel()                    # velocity components        [cm/s]
    vy1    = vy1    *  cgs.cu_vel()
    vz1    = vz1    *  cgs.cu_vel()
    maccr1 = maccr1 *  cgs.Msun_gram()                 # accreted mass              [g]

    r1 = gf.calc_r(x1, y1, z1)                         # [cm]

    position1 = np.array((x1, y1, z1 )).transpose()
    velocity1 = np.array((vx1,vy1,vz1)).transpose()
    
    # companion star

    t2     = t2     *  cgs.cu_time()                   # evolution time             [yrs]                                 
    x2     = x2     *  cgs.AU_cm()                     # position coordinates       [cm]
    y2     = y2     *  cgs.AU_cm()       
    z2     = z2     *  cgs.AU_cm()      
    mass2  = mass2  *  cgs.Msun_gram()                 # mass of sph particles      [g]
    vx2    = vx2    *  cgs.cu_vel()                    # velocity components        [cm/s]
    vy2    = vy2    *  cgs.cu_vel()
    vz2    = vz2    *  cgs.cu_vel()
    maccr2 = maccr2 *  cgs.Msun_gram()                 # accreted mass              [g]

    r2 = gf.calc_r(x2, y2, z2)                            # [cm]

    position2 = np.array((x2, y2, z2 )).transpose()
    velocity2 = np.array((vx2,vy2,vz2)).transpose()
    
    
    # orbital information 
    
    period          = pq.getPeriod(mass1, mass2, (r1 + r2) /cgs.AU_cm()  )
    orbitalVel_AGB  = pq.getOrbitalVelocity(period, r1     /cgs.AU_cm()  )
    orbitalVel_comp = pq.getOrbitalVelocity(period, r2     /cgs.AU_cm()  )
    rHill           = pq.getRHill( abs(r1 + r2), mass2, mass1            )         # [cm]
    
    # output
    #    "_t" stands for the fact that these values are in function of the evolution time, not from the last dump in function of location
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
            'period_t'    : period,                 # [s]
            'v_orbAGB_t'  : orbitalVel_AGB,         # [cm/s]
            'v_orbComp_t' : orbitalVel_comp,        # [cm/s]
            'rHill_t'     : rHill                   # [cm]
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
    x1     = x1     *  cgs.AU_cm()                     # position coordinates       [cm]
    y1     = y1     *  cgs.AU_cm()       
    z1     = z1     *  cgs.AU_cm()      
    mass1  = mass1  *  cgs.Msun_gram()                 # mass of sph particles      [g]
    vx1    = vx1    *  cgs.cu_vel()                    # velocity components        [cm/s]
    vy1    = vy1    *  cgs.cu_vel()
    vz1    = vz1    *  cgs.cu_vel()
    maccr1 = maccr1 *  cgs.Msun_gram()                 # accreted mass              [g]

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