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
def LoadDump_cgs(run, loc, setup, userSettingsDictionary, number):
    runName = os.path.join(loc,run)
    userPrefix = userSettingsDictionary["prefix"]
    phantom_dir = userSettingsDictionary["hard_path_to_phantom"]
    sys.path.append(phantom_dir+"/scripts")
    from readPhantomDump import read_dump

    # Pick either last dump file or user chosen file
    index = 0
    if number == -1: index = findLastFullDumpIndex(userPrefix, setup, runName)
    else: index = number


    # make filename of this filenumber
    fileName       = runName+'/{0:s}_{1:05d}'.format(userPrefix, index)

    # reading the dumpfile
    dump = read_dump(fileName)
    
    if "tau" in dump["blocks"][0]["data"]: containsTau = True
    else: containsTau=False
    if "Tdust" in dump["blocks"][0]["data"]: containsTemp = True
    else: containsTemp=False

    x = dump["blocks"][0]["data"]["x"]
    y = dump["blocks"][0]["data"]["y"]
    z = dump["blocks"][0]["data"]["z"]
    h = dump["blocks"][0]["data"]["h"]
    mass = np.ones(len(h))*dump["quantities"]["massoftype"][0]
    vx = dump["blocks"][0]["data"]["vx"]
    vy = dump["blocks"][0]["data"]["vy"]
    vz = dump["blocks"][0]["data"]["vz"]
    u = dump["blocks"][0]["data"]["u"]
    if containsTau:  tau  = dump["blocks"][0]["data"]["tau"]
    if containsTemp: temp = dump["blocks"][0]["data"]["Tdust"]
    
    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x                     [h > 0.0] * cgs.AU_cm()       # position coordinates          [cm]
    y     = y                     [h > 0.0] * cgs.AU_cm()       
    z     = z                     [h > 0.0] * cgs.AU_cm()      
    mass  = mass                  [h > 0.0] * cgs.Msun_gram()   # mass of sph particles         [g]
    vx    = vx                    [h > 0.0]                     # velocity components           [cm/s]
    vy    = vy                    [h > 0.0]
    vz    = vz                    [h > 0.0]
    u     = u                     [h > 0.0]                     # specific internal density     [erg/g]
    if containsTemp: temp = temp  [h > 0.0]                     # temperature                   [K]
    else:            temp = pq.getTemp(p, rho, setup['mu'], u)
    if containsTau: tau   = tau   [h > 0.0]                     # optical depth                 
    h     = h                     [h > 0.0] * cgs.AU_cm()       # smoothing length              [cm]
    rho   = pq.getRho(h, dump["quantities"]["hfact"], mass)     # density                       [g/cm^3]
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
    cs    = pq.getSoundSpeed(p, rho, setup['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates

    
    # Data of the two stars
    # AGB star
    xAGB     = dump["blocks"][1]["data"]["x"][0]
    yAGB     = dump["blocks"][1]["data"]["y"][0]
    zAGB     = dump["blocks"][1]["data"]["z"][0]
    posAGB   = [xAGB, yAGB, zAGB]
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)   
    massAGB  = dump["blocks"][1]["data"]["m"][0]
    
    # companion
    xComp    = dump["blocks"][1]["data"]["x"][1]
    yComp    = dump["blocks"][1]["data"]["y"][1]
    zComp    = dump["blocks"][1]["data"]["z"][1]
    posComp  = [xComp, yComp, zComp]
    rComp    = gf.calc_r(xComp, yComp, zComp)
    massComp = dump["blocks"][1]["data"]["m"][1]
    
    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()
    
    speed = np.linalg.norm(velocity, axis=1)
    mach  = speed/cs
    vtvv  = (vtan/speed)**2      # fraction of the velocity that is tangential: if vtvv > 0.5 -> tangential
    rHill = pq.getRHill(abs(rComp+rAGB),massComp,massAGB)
    
  
    # get normal of the edge-on plane
    #normal_edgeon_plane = gf.getNormalPerpTo(AGBcoord, compCoord, [1,1,0])
    

    
    # output
    data = {'position'      : position,              # [cm]
            'velocity'      : velocity,              # [cm/s]
            'h'             : h,                     # [cm]
            'mass'          : mass,                  # [g]
            'rho'           : rho,                   # [g/cm^3]
            'u'             : u,                     # [erg/g]
            'temp'          : temp,                  # [K] 
            'speed'         : speed,                 # [cm/s]
            'mach'          : mach,                  
            'vtvv'          : vtvv,     
            'fileNumber'    : index,
            'r'             : r,                     # [cm]
            'phi'           : phi,     
            'theta'         : theta,     
            'cs'            : cs,                    # [cm]
            'posAGB'        : posAGB,                # [cm]
            'rAGB'          : rAGB,                  # [cm]
            'massAGB'       : massAGB,               # [g]
            'posComp'       : posComp,               # [cm]
            'rComp'         : rComp,                 # [cm]
            'massComp'      : massComp,              # [g]
            'rHill'         : rHill,                 # [cm]
            'vx'            : vx,                    # [cm/s]
            'vy'            : vy,                    # [cm/s]
            'vz'            : vz                     # [cm/s]
            }
    if containsTau: data["tau"] = tau

    
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
def LoadDump_single_cgs(run, loc, setup, userSettingsDictionary):

    runName = os.path.join(loc, run)
    userPrefix = userSettingsDictionary["prefix"]
    phantom_dir = userSettingsDictionary["hard_path_to_phantom"]
    sys.path.append(phantom_dir+"/scripts")
    from readPhantomDump import read_dump

    # Pick last file from model
    lastFullDumpIndexInt = findLastFullDumpIndex(userPrefix, setup, runName)

    # make filename of this filenumber
    fileName  = runName+'/{0:s}_{1:05d}'.format(userPrefix, lastFullDumpIndexInt)

    # reading the dumpfile
    dump = read_dump(fileName)

    if "tau" in dump["blocks"][0]["data"]: containsTau = True
    else: containsTau=False
    if "Tdust" in dump["blocks"][0]["data"]: containsTemp = True
    else: containsTemp=False

    x = dump["blocks"][0]["data"]["x"]
    y = dump["blocks"][0]["data"]["y"]
    z = dump["blocks"][0]["data"]["z"]
    h = dump["blocks"][0]["data"]["h"]
    mass = np.ones(len(h))*dump["quantities"]["massoftype"][0]
    vx = dump["blocks"][0]["data"]["vx"]
    vy = dump["blocks"][0]["data"]["vy"]
    vz = dump["blocks"][0]["data"]["vz"]
    u = dump["blocks"][0]["data"]["u"]
    if containsTau:  tau  = dump["blocks"][0]["data"]["tau"]
    if containsTemp: temp = dump["blocks"][0]["data"]["Tdust"]
        
    
    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x                     [h > 0.0] * cgs.AU_cm()       # position coordinates          [cm]
    y     = y                     [h > 0.0] * cgs.AU_cm()       
    z     = z                     [h > 0.0] * cgs.AU_cm()      
    mass  = mass                  [h > 0.0] * cgs.Msun_gram()   # mass of sph particles         [g]
    vx    = vx                    [h > 0.0] * cgs.cu_vel()      # velocity components           [cm/s]
    vy    = vy                    [h > 0.0] * cgs.cu_vel()
    vz    = vz                    [h > 0.0] * cgs.cu_vel()
    u     = u                     [h > 0.0] * cgs.cu_e()        # specific internal density     [erg/g]
    if containsTemp: temp = temp  [h > 0.0]                     # temperature                   [K]
    else:            temp = pq.getTemp(p, rho, setup['mu'], u)
    if containsTau: tau   = tau   [h > 0.0]                     # optical depth                 
    h     = h                     [h > 0.0] * cgs.AU_cm()       # smoothing length              [cm]
    rho   = pq.getRho(h, dump["quantities"]["hfact"], mass)     # density                       [g/cm^3]
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
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
            'fileNumber'    : lastFullDumpIndexInt,
            'r'             : r,              # [cm]
            'phi'           : phi,            
            'theta'         : theta,
            'cs'            : cs,             # [cm]
            'vx'            : vx,             # [cm/s]
            'vy'            : vy,             # [cm/s]
            'vz'            : vz              # [cm/s]
            }                    
    if containsTau: data["tau"] = tau             

    
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
    temp  = dump['temp']
    if "tau" in dump: tau   = dump['tau']
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
    temp  = temp  [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 # K
    if "tau" in dump: tau = tau   [r > factor * setup['sma_ini'] * cgs.AU_cm() ]                                 
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
    temp  = temp  [r < bound * cgs.AU_cm()]                                 # K
    if "tau" in dump: tau = tau   [r < bound * cgs.AU_cm()]                                 # K
    h     = h     [r < bound * cgs.AU_cm()]                                 # cm
    r     = r     [r < bound * cgs.AU_cm()] 
    
    
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
    # temp  = pq.getTemp(p, rho, setup['mu'], u)                  # temperature                   [K]
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
    if "tau" in dump: data["tau"] = tau
    
    return data

# Pick last full dump file from model
def findLastFullDumpIndex(userPrefix, setup, runName):
    listDumps = sortedDumpList(userPrefix, runName)
    lastFile = listDumps[-1]
    lastDumpIndex = int(lastFile.lstrip("%s_0" % userPrefix))
   
    # Nearest full dump to max
    lastFullDumpIndex = int(int(math.floor(lastDumpIndex / setup['nfulldump'])) * setup['nfulldump']) 
    return lastFullDumpIndex

def findAllFullDumpIndices(userSettingsDictionary, setup, runName):
    userPrefix = userSettingsDictionary['prefix']
    listDumps = sortedDumpList(userPrefix, runName)
    lastFullDumpIndex = findLastFullDumpIndex(userPrefix, setup, runName)
    fullDumpLists = np.arange(0, lastFullDumpIndex + 1, setup['nfulldump'], dtype=int)
    return fullDumpLists

def sortedDumpList(userPrefix, runName):
    return np.sort(list(filter(lambda x: ("%s_"%userPrefix in x) and (not (".asc" in x or ".dat" in x)), os.listdir(runName))))