import numpy                    as np
import math                     as math
import os
# own scripts
import PhysicalQuantities       as pq
import GeometricalFunctions     as gf
import ConversionFactors_cgs    as cgs
import LoadSetup                as stp





'''
Load the .ev-files from a phantom model
      - This file gives the specifics of the objects only present in the model (not of the sph particles)
        in function of the evolution time of the model. The last entry corresponds to the data from the last dump.
      - Only suited for a binary model
      - Units in cgs
'''
def LoadSink_cgs(run, loc, setup):

    
    runName = loc + run

    fileNInt   = int(setup['tmax'])
    fileNumber = str(fileNInt)

    if fileNInt  < 100:
        fileNumber = str(0) + fileNumber

    # make ascii file of this filenumber    
    fileName  = runName+'/wind_00'+fileNumber +'.ascii'
    
    # load the dump file wind_00xxx
    try:
    # to calculate period, we need masses and sma, so coordinates, we call the parameters xI, I stands for input
        (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(runName+'/windSink0001N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        (t2, x2,y2,z2, mass2, vx2,vy2,vz2, maccr2) = np.loadtxt(runName+'/windSink0002N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        
        
    except OSError:

        print('Copying file from STER...')
        
        os.system("scp -r silkem@copernicus.ster.kuleuven.be:/STER/silkem/THESIS/phantom_Masterthesis/desktop_run"+run+"/windSink000* /home/silke/Documents/Univ/THESIS/Models/phantom_Masterthesis/desktop_run"+runNumber+"/")

        print('Converting file to ascii...')
        
        os.system('cd')
        os.system('cd '+runName)
        #os.system("ssplash to ascii "+runName+"/wind_00"+fileNumber)
        (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(runName+'/windSink0001N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        (t2, x2,y2,z2, mass2, vx2,vy2,vz2, maccr2) = np.loadtxt(runName+'/windSink0002N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)



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
    
    period          = pq.getPeriod(mass1, mass2, (r1 + r2) /cgs.AU_cm()  )
    orbitalVel_AGB  = pq.getOrbitalVelocity(period, r1     /cgs.AU_cm()  )
    orbitalVel_comp = pq.getOrbitalVelocity(period, r2     /cgs.AU_cm()  )
    rHill           = pq.getRHill(abs(r1 + r2), mass2, mass1             )         # [cm]
    
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






# Load the .ev-files for a single phantom model
def LoadSink_single_cgs(run, loc, setup):

    runName = loc + run

    fileNInt   = int(setup['tmax'])    
    fileNumber = str(fileNInt)

    if fileNInt  < 100:
        fileNumber = str(0) + fileNumber

    # make ascii file of this filenumber    
    fileName  = runName+'/wind_00'+fileNumber +'.ascii'
    
    # load the dump file wind_00xxx
    try:
       (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(runName+'/windSink0001N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        
        
    except OSError:

        print('Copying file from STER...')
        
        os.system("scp -r silkem@copernicus.ster.kuleuven.be:/STER/silkem/THESIS/phantom_Masterthesis/desktop_run"+runNumber+"/windSink000* /home/silke/Documents/Univ/THESIS/Models/phantom_Masterthesis/desktop_run"+runNumber+"/")

        print('Converting file to ascii...')
        
        os.system('cd')
        os.system('cd '+runName)
        #os.system("ssplash to ascii "+runName+"/wind_00"+fileNumber)
        (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(runName+'/windSink0001N01.ev', skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        

   

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

    r1 = gf.calc_r(x1, y1, z1)                            # cm

    position1 = np.array((x1, y1, z1 )).transpose()
    velocity1 = np.array((vx1,vy1,vz1)).transpose()
    
    
    
    data = {'posAGB'      : position1,
            'velAGB'      : velocity1,
            'massAGB'     : mass1,
            'maccrAGB'    : maccr1,
            'rAGB'        : r1,
            'time'        : t1,
            }


    
    return data


