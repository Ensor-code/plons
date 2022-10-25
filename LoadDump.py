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
    
INPUT:
    - 'run'   is the number of the run specifically           [str]
    - 'loc'   is the directory where the model is located     [str]
    - 'setup' is the setup data                               [dict]
    
RETURNS
    a dictionary containing the data from the dump (all units in cgs)
'''
def LoadDump_cgs(run, loc, setup, userSettingsDictionary, number = -1):

    # Load read_dump script to read PHANTOM dump files
    runName = os.path.join(loc,run)
    userPrefix = userSettingsDictionary["prefix"]
    phantom_dir = userSettingsDictionary["hard_path_to_phantom"]
    sys.path.append(phantom_dir+"/scripts")
    from readPhantomDump import read_dump

    # Pick either last dump file or user chosen file
    if number == -1: index = findLastFullDumpIndex(userPrefix, setup, runName)
    else: index = number

    # make filename of this filenumber
    fileName       = runName+'/{0:s}_{1:05d}'.format(userPrefix, index)

    # reading the dumpfile
    dump = read_dump(fileName)

    unit_dist = dump['units']['udist']
    unit_mass = dump['units']['umass']
    unit_time = dump['units']['utime']
    
    
    unit_velocity = unit_dist/unit_time
    unit_density  = unit_mass/unit_dist**3
    unit_pressure = unit_mass/(unit_dist*unit_time**2)
    unit_ergg     = unit_velocity**2
    unit_energ    = unit_mass*unit_ergg
    unit_opacity  = unit_dist**2/unit_mass

    # Data of the two stars
    # AGB star
    xAGB     = dump["blocks"][1]["data"]["x"][0] * cgs.AU_cm()
    yAGB     = dump["blocks"][1]["data"]["y"][0] * cgs.AU_cm()
    zAGB     = dump["blocks"][1]["data"]["z"][0] * cgs.AU_cm()
    posAGB   = [xAGB, yAGB, zAGB]
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)   
    massAGB  = dump["blocks"][1]["data"]["m"][0] * unit_mass
    lumAGB   = dump["blocks"][1]["data"]["lum"][0] * unit_energ/unit_time

    # companion
    if not setup["single_star"]:
        xComp    = dump["blocks"][1]["data"]["x"][1] * cgs.AU_cm()
        yComp    = dump["blocks"][1]["data"]["y"][1] * cgs.AU_cm()
        zComp    = dump["blocks"][1]["data"]["z"][1] * cgs.AU_cm()
        posComp  = [xComp, yComp, zComp]
        rComp    = gf.calc_r(xComp, yComp, zComp)
        massComp = dump["blocks"][1]["data"]["m"][1] * unit_mass
        rHill = pq.getRHill(abs(rComp+rAGB),massComp,massAGB)

    containsTau = "tau" in dump["blocks"][0]["data"]
    containsTemp = "Tdust" in dump["blocks"][0]["data"]
    bowenDust = "bowen_kmax" in setup
    
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
    if containsTemp:
        temp = dump["blocks"][0]["data"]["Tdust"]
        if (len(temp) == 2*len(h)): temp = temp[1::2] # A problem with loading the files, info: Mats

    filter = h > 0.0
    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x                     [filter] * unit_dist         # position coordinates          [cm]
    y     = y                     [filter] * unit_dist     
    z     = z                     [filter] * unit_dist
    mass  = mass                  [filter] * unit_mass         # mass of sph particles         [g]
    vx    = vx                    [filter] * unit_velocity * cgs.cms_kms()                   # velocity components           [cm/s]
    vy    = vy                    [filter] * unit_velocity * cgs.cms_kms()
    vz    = vz                    [filter] * unit_velocity * cgs.cms_kms()
    u     = u                     [filter] * unit_energ        # specific internal density     [erg/g]
    if containsTemp: temp = temp  [filter]                     # temperature                   [K]
    if containsTau and len(tau) == len(h): tau = tau [filter]  # optical depth                 
    h     = h                     [filter] * unit_dist         # smoothing length              [cm]
    rho   = pq.getRho(h, dump["quantities"]["hfact"], mass)     # density                       [g/cm^3]
    p     = pq.getPressure(rho, u, dump['quantities']['gamma']) # pressureure                   [Ba = 1e-1 Pa]
    if not containsTemp: temp = pq.getTemp(p, rho, dump['quantities']['gamma'], setup['mu'], u) # temperature                [K]
    if bowenDust: kappa = pq.getKappa(temp, setup['kappa_gas'], setup['bowen_delta'], setup['bowen_Tcond'], setup['bowen_kmax'])
    if bowenDust: Gamma = pq.getGamma(kappa, lumAGB, massAGB)
    cs    = pq.getSoundSpeed(p, rho, dump['quantities']['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates

    
    # Data of the two stars
    # AGB star
    xAGB     = dump["blocks"][1]["data"]["x"][0] * cgs.AU_cm()
    yAGB     = dump["blocks"][1]["data"]["y"][0] * cgs.AU_cm()
    zAGB     = dump["blocks"][1]["data"]["z"][0] * cgs.AU_cm()
    posAGB   = [xAGB, yAGB, zAGB]
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)   
    massAGB  = dump["blocks"][1]["data"]["m"][0] * cgs.Msun_gram()
    
    if not setup["single_star"]:
        # companion
        xComp    = dump["blocks"][1]["data"]["x"][1] * cgs.AU_cm()
        yComp    = dump["blocks"][1]["data"]["y"][1] * cgs.AU_cm()
        zComp    = dump["blocks"][1]["data"]["z"][1] * cgs.AU_cm()
        posComp  = [xComp, yComp, zComp]
        rComp    = gf.calc_r(xComp, yComp, zComp)
        massComp = dump["blocks"][1]["data"]["m"][1] * cgs.Msun_gram()
    
    if setup['triple_star']:
        # inner companion
        xComp_in    = dump["blocks"][1]["data"]["x"][2] * cgs.AU_cm()
        yComp_in    = dump["blocks"][1]["data"]["y"][2] * cgs.AU_cm()
        zComp_in    = dump["blocks"][1]["data"]["z"][2] * cgs.AU_cm()
        posComp_in  = [xComp_in, yComp_in, zComp_in]
        rComp_in    = gf.calc_r(xComp_in, yComp_in, zComp_in)
        massComp_in = dump["blocks"][1]["data"]["m"][2] * cgs.Msun_gram()
        
    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()
        
    speed = np.linalg.norm(velocity, axis=1)
    mach  = speed/cs
    vtvv  = (vtan/speed)**2      # fraction of the velocity that is tangential: if vtvv > 0.5 -> tangential
    
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
            'vx'            : vx,                    # [cm/s]
            'vy'            : vy,                    # [cm/s]
            'vz'            : vz                     # [cm/s]
            }

    if not setup["single_star"]:
        data['posComp'    ] = posComp                # [cm]
        data['rComp'      ] = rComp                  # [cm]
        data['massComp'   ] = massComp               # [g]
        data['rHill'      ] = rHill                  # [cm]

    if setup['triple_star']:
        data["posComp_in"]  = posComp_in               # [cm]
        data['rComp_in']    = rComp_in                 # [cm]
        data['massComp_in'] = massComp_in              # [g]
        

    if containsTau:
        data["tau"        ] = tau

    if bowenDust:
        data["kappa"      ] = kappa                  # [cm^2/g]
        data["Gamma"      ] = Gamma
    
    return data

'''
Cut out the inner part of the wind, since here the wind is not yet self-similar.
      Suited only for binary model
    
INPUT:
    - 'factor' is the factor of the semi-major axis you want to cut out [float]
    - 'bound'  is a chosen boundary of the model                        [float]
    - 'setup'  is the setup data                                        [dict]
    - 'dump'   is the original dump data                                [dict]

RETURNS
    a dictionary containing the data from the dump of a chosen outer range (all units in cgs)
'''
def LoadDump_outer_cgs(factor, bound, setup, dump):

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
    if "tau"   in dump: tau   = dump['tau']
    if "kappa" in dump: kappa = dump['kappa']
    if "Gamma" in dump: Gamma = dump['Gamma']
    h     = dump['h']
    r     = dump['r']

    filter = r > factor * setup['sma_ini'] * cgs.AU_cm() and r < bound * cgs.AU_cm()
    x     = x                         [filter]             # cm
    y     = y                         [filter]
    z     = z                         [filter]
    mass  = mass                      [filter]             # g
    vx    = vx                        [filter]             # cm/s
    vy    = vy                        [filter]
    vz    = vz                        [filter]
    u     = u                         [filter]             # erg/g
    rho   = rho                       [filter]             # g/cm^3
    temp  = temp                      [filter]             # K
    if "tau"   in dump: tau   = tau   [filter]
    if "kappa" in dump: kappa = kappa [filter]             # cm^2/g
    if "Gamma" in dump: Gamma = Gamma [filter]
    h     = h                         [filter]             # cm
    r     = r                         [filter]
    
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressureure                   [Ba = 1e-1 Pa]
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
    if "tau"   in dump: data["tau"]   = tau
    if "kappa" in dump: data["kappa"] = kappa
    if "Gamma" in dump: data["Gamma"] = Gamma
    
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
