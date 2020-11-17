import numpy                    as np
import math                     as math
import os
# own scripts
import PhysicalQuantities       as pq
import GeometricalFunctions     as gf
import ConversionFactors_cgs    as cgs
import LoadSetup                as stp


'''
Loads the final full dump of a phantom model, given the number, in cgs-units 
(converts it to ascii or copies it from the server, if needed)
    Suited for binary model only
    
INPUT:
    - 'run'   is the number of the run specifically           [str]
    - 'loc'   is the directory where the model is located     [str]
    - 'setup' is the setup data                               [dict]
    
RETURNS
    a dictionary containing the data from the dump (all units in cgs)
'''
def LoadDump_cgs(run, loc, setup):
    
    runName = loc + run
    

    fileNInt   = int(setup['tmax'])
    fileNumber = str(fileNInt)

    if fileNInt  < 100:
        fileNumber = str(0) + fileNumber

    # make ascii file of this filenumber    
    fileName  = runName+'/wind_00'+fileNumber +'.ascii'
    
    # load the dump file wind_00xxx
    try:
        '''
        We don't need the information about the AGB star and companion star (two last rows of the dumps).
        To efficiently skips these rows, we load only one parameter, to get the length. All other parameters
        are loaded after without the two last rows.
        '''
        (x,y,z,mass, h, rho, vx,vy,vz, u) = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
     
    except OSError:
        try: 
            print('Converting dump file to ascii...')
            
            os.system("splash to ascii "+runName+"/wind_00"+fileNumber)
            (x,y,z,mass, h, rho, vx,vy,vz, u) = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0,1,2,3,4,5,6,7,8,9),  unpack=True)
     
        
        except OSError:
            print(' ERROR: No dump file found for this model in the current directory!')
            
            
    #rows = len(x)       
    #(x,y,z,mass, h, rho, vx,vy,vz, u) = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0,1,2,3,4,5,6,7,8,9), max_rows = rows-2, unpack=True)
        
    
    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x     [h > 0.0] * cgs.AU_cm()                       # position coordinates          [cm]
    y     = y     [h > 0.0] * cgs.AU_cm()       
    z     = z     [h > 0.0] * cgs.AU_cm()      
    mass  = mass  [h > 0.0] * cgs.Msun_gram()                   # mass of sph particles         [g]
    vx    = vx    [h > 0.0] * cgs.cu_vel()                      # velocity components           [cm/s]
    vy    = vy    [h > 0.0] * cgs.cu_vel()
    vz    = vz    [h > 0.0] * cgs.cu_vel()
    u     = u     [h > 0.0] * cgs.cu_e()                        # specific internal density     [erg/g]
    rho   = rho   [h > 0.0] * cgs.cu_dens()                     # density                       [g/cm^3]
    h     = h     [h > 0.0] * cgs.AU_cm()                       # smoothing length              [cm]
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
    temp  = pq.getTemp(p, rho, setup['mu'])                     # temperature                   [K]
    cs    = pq.getSoundSpeed(p, rho, setup['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates

    
    # Data of the two stars
    # AGB star
    xAGB     = x[-2]
    yAGB     = y[-2]
    zAGB     = z[-2]
    posAGB   = [xAGB, yAGB, zAGB]
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)   
    massAGB  = mass[-2] 
    
    # companion
    xComp    = x[-1]
    yComp    = y[-1]
    zComp    = z[-1]
    posComp  = [xComp, yComp, zComp]
    rComp    = gf.calc_r(xComp, yComp, zComp)
    massComp = mass[-1]
    
    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()
    
    speed = np.linalg.norm(velocity, axis=1)
    mach  = speed/cs
    vtvv  = (vtan/speed)**2      # fraction of the velocity that is tangential: if vtvv > 0.5 -> tangential
    rHill = pq.getRHill(abs(rComp+rAGB),massComp,massAGB)
    
  
    # get normal of the edge-on plane
    #normal_edgeon_plane = gf.getNormalPerpTo(AGBcoord, compCoord, [1,1,0])
    

    
    # output
    data = {'position'      : position[:-2],         # [cm]
            'velocity'      : velocity[:-2],         # [cm/s]
            'h'             : h[:-2],                # [cm]
            'mass'          : mass[:-2],             # [g]
            'rho'           : rho[:-2],              # [g/cm^3]
            'u'             : u[:-2],                # [erg/g]
            'temp'          : temp[:-2],             # [K]    
            'speed'         : speed[:-2],            # [cm/s]
            'mach'          : mach[:-2],             
            'vtvv'          : vtvv[:-2],
            'fileNumber'    : fileNumber,
            'r'             : r[:-2],                # [cm]
            'phi'           : phi[:-2],
            'theta'         : theta[:-2],
            'cs'            : cs[:-2],               # [cm]
            'posAGB'        : posAGB,           # [cm]
            'rAGB'          : rAGB,             # [cm]
            'massAGB'       : massAGB,          # [g]
            'posComp'       : posComp,          # [cm]
            'rComp'         : rComp,            # [cm]
            'massComp'      : massComp,         # [g]
            'rHill'         : rHill             # [cm]
            }
    

    
    return data


'''
Loads the final full dump of a phantom model, given the number, in cgs-units 
(converts it to ascii or copies it from the server, if needed)
    Suited for single model only
    
INPUT:
    - 'run'   is the number of the run specifically           [str]
    - 'loc'   is the directory where the model is located     [str]
    - 'setup' is the setup data                               [dict]
    
RETURNS
    a dictionary containing the data from the dump (all units in cgs)
'''
def LoadDump_single_cgs(run, loc, setup):
    
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
        x = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0), unpack=True)

    except OSError:
        try:
            print('Converting dump file to ascii...')
            
            os.system("splash to ascii "+runName+"/wind_00"+fileNumber)
            x = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0), unpack=True)
        
        except OSError:
            print(' ERROR: No dump file found for this model in the current directory!')
        
    rows = len(x)       
    (x,y,z,mass, h, rho, vx,vy,vz, u) = np.loadtxt(runName+'/wind_00'+str(fileNumber)+'.ascii', skiprows=14, usecols=(0,1,2,3,4,5,6,7,8,9), max_rows = rows-1, unpack=True)
    

    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x     [h > 0.0] * cgs.AU_cm()                       # position coordinates          [cm]
    y     = y     [h > 0.0] * cgs.AU_cm()       
    z     = z     [h > 0.0] * cgs.AU_cm()      
    mass  = mass  [h > 0.0] * cgs.Msun_gram()                   # mass of sph particles         [g]
    vx    = vx    [h > 0.0] * cgs.cu_vel()                      # velocity components           [cm/s]
    vy    = vy    [h > 0.0] * cgs.cu_vel()
    vz    = vz    [h > 0.0] * cgs.cu_vel()
    u     = u     [h > 0.0] * cgs.cu_e()                        # specific internal density     [erg/g]
    rho   = rho   [h > 0.0] * cgs.cu_dens()                     # density                       [g/cm^3]
    h     = h     [h > 0.0] * cgs.AU_cm()                       # smoothing length              [cm]
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
    temp  = pq.getTemp(p, rho, setup['mu'])                     # temperature                   [K]
    cs    = pq.getSoundSpeed(p, rho, setup['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates

    
    
    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()
    
    speed = np.linalg.norm(velocity, axis=1)
    mach  = speed/cs
    vtvv  = (vtan/speed)**2
    
        
    
    # output
    data = {'position'      : position,       # [cm]
            'velocity'      : velocity,       # [cm/s]
            'h'             : h,              # [cm]
            'mass'          : mass,           # [g]
            'rho'           : rho,            # [g/cm^3]
            'u'             : u,              # [erg/g]
            'temp'          : temp,           # [K]    
            'speed'         : speed,          # [cm/s]
            'mach'          : mach,           
            'vtvv'          : vtvv,
            'fileNumber'    : fileNumber,
            'r'             : r,              # [cm]
            'phi'           : phi,            
            'theta'         : theta,
            'cs'            : cs              # [cm]
            }                                 

    
    return data



'''
Cut out the inner part of the wind, since here the wind is not yet self-similar.
      Suited only for binary model
    
INPUT:
    - 'run'    is the number of the run specifically                    [str]
    - 'loc'    is the directory where the model is located              [str]
    - 'factor' is the factor of the semi-major axis you want to cut out [float]
    - 'bound'  is a chosen boundary of the model                        [float]
    - 'setup'  is the setup data                                        [dict]
    - 'dump'   is the original dump data                                [dict]

RETURNS
    a dictionary containing the data from the dump of a chosen outer range (all units in cgs)
'''
def LoadDump_outer_cgs(run, loc, factor, bound, setup, dump):
    

    
    position = dump['position'].transpose()
    velocity = dump['velocity'].transpose()
    
    x     = position[0]
    y     = position[1]
    z     = position[2]
    mass  = dump['mass']
    vx    = velocity[0]
    vy    = velocity[1]
    vz    = velocity[2]
    u     = dump['u']
    rho   = dump['rho']
    h     = dump['h']
    r     = dump['r']

    
    x     = x     [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # cm
    y     = y     [r > factor * setup['sma_ini'] * cgs.AU_cm() ]               
    z     = z     [r > factor * setup['sma_ini'] * cgs.AU_cm() ]               
    mass  = mass  [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # g
    vx    = vx    [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # cm/s
    vy    = vy    [r > factor * setup['sma_ini'] * cgs.AU_cm() ]               
    vz    = vz    [r > factor * setup['sma_ini'] * cgs.AU_cm() ]               
    u     = u     [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # erg/g
    rho   = rho   [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # g/cm^3
    h     = h     [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # cm
    r     = r     [r > factor * setup['sma_ini'] * cgs.AU_cm() ] 
    
    
    x     = x     [r < bound * cgs.AU_cm()]                                 # cm
    y     = y     [r < bound * cgs.AU_cm()]               
    z     = z     [r < bound * cgs.AU_cm()]               
    mass  = mass  [r < bound * cgs.AU_cm()]                                 # g
    vx    = vx    [r < bound * cgs.AU_cm()]                                 # cm/s
    vy    = vy    [r < bound * cgs.AU_cm()]               
    vz    = vz    [r < bound * cgs.AU_cm()]               
    u     = u     [r < bound * cgs.AU_cm()]                                 # erg/g
    rho   = rho   [r < bound * cgs.AU_cm()]                                 # g/cm^3
    h     = h     [r < bound * cgs.AU_cm()]                                 # cm
    r     = r     [r < bound * cgs.AU_cm()] 
    
    
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
    temp  = pq.getTemp(p, rho, setup['mu'])                     # temperature                   [K]
    cs    = pq.getSoundSpeed(p, rho, setup['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates

    
   
    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()
    
    speed = np.linalg.norm(velocity, axis=1)
    mach  = speed/cs
    vtvv  = (vtan/speed)**2      # fraction of the velocity that is tangential: if vtvv > 0.5 -> tangential
  
 
    
    # output
    data = {'position'      : position,      # [cm]
            'velocity'      : velocity,      # [cm/s]
            'h'             : h,             # [cm]
            'mass'          : mass,          # [g]
            'rho'           : rho,           # [g/cm^3]
            'u'             : u,             # [erg/g]
            'temp'          : temp,          # [K]    
            'speed'         : speed,         # [cm/s]
            'mach'          : mach,          
            'vtvv'          : vtvv,
            'r'             : r,             # [cm]
            'phi'           : phi,           
            'theta'         : theta,         
            'cs'            : cs             # [cm]
            }                                
                                             
    

    
    return data


