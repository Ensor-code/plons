import numpy                    as np
import math                     as math
import os
import sys
import time

# Import plons scripts
import plons.PhysicalQuantities       as pq
import plons.GeometricalFunctions     as gf
import plons.ConversionFactors_cgs    as cgs

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

    print('Dump file read')
    unit_dist = dump['units']['udist']
    unit_mass = dump['units']['umass']
    unit_time = dump['units']['utime']

    unit_velocity = unit_dist/unit_time
    unit_density  = unit_mass/unit_dist**3
    unit_pressure = unit_mass/(unit_dist*unit_time**2)
    unit_ergg     = unit_velocity**2
    unit_energ    = unit_mass*unit_ergg
    unit_opacity  = unit_dist**2/unit_mass
    xarray=dump["blocks"][1]["data"]["x"] * unit_dist
    yarray=dump["blocks"][1]["data"]["y"] * unit_dist
    zarray=dump["blocks"][1]["data"]["z"] * unit_dist
    marray=dump["blocks"][1]["data"]["m"] * unit_mass
    # Data of the two stars
    # AGB star
    xAGB     = xarray[0]
    yAGB     = yarray[0]
    zAGB     = zarray[0]
    posAGB   = [xAGB, yAGB, zAGB]
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)
    massAGB  = marray[0]
    lumAGB   = dump["blocks"][1]["data"]["lum"][0] * unit_energ/unit_time
    # companion
    if not setup["single_star"]:
        xComp    = xarray[1]
        yComp    = yarray[1]
        zComp    = zarray[1]
        posComp  = [xComp, yComp, zComp]
        rComp    = gf.calc_r(xComp, yComp, zComp)
        massComp = marray[1]
        rHill = pq.getRHill(abs(rComp+rAGB),massComp,massAGB)

        if setup['triple_star']:
            # inner companion
            xComp_in    = xarray[2]
            yComp_in    = yarray[2]
            zComp_in    = zarray[2]
            posComp_in  = [xComp_in, yComp_in, zComp_in]
            rComp_in    = gf.calc_r(xComp_in, yComp_in, zComp_in)
            massComp_in = marray[2]
        if setup['quadruple_star']:
            #0 is agb, 1 is outer, 2 is close companion, 3 is outer

            xComp_in    = xarray[2]
            yComp_in    = yarray[2]
            zComp_in    = zarray[2]
            posComp_in  = [xComp_in, yComp_in, zComp_in]
            rComp_in    = gf.calc_r(xComp_in, yComp_in, zComp_in)
            massComp_in = marray[2]

            xComp_out    = xarray[3]
            yComp_out    = yarray[3]
            zComp_out    = zarray[3]
            posComp_in2  = [xComp_out, yComp_out, zComp_out]
            rComp_out    = gf.calc_r(xComp_out, yComp_out, zComp_out)
            massComp_out = marray[3]


    containsTau  = "tau" in dump["blocks"][0]["data"]
    containsTauL = "tau_lucy" in dump["blocks"][0]["data"]
    bowenDust    = "bowen_kmax" in setup

    x = dump["blocks"][0]["data"]["x"]
    y = dump["blocks"][0]["data"]["y"]
    z = dump["blocks"][0]["data"]["z"]
    h = dump["blocks"][0]["data"]["h"]
    mass = np.ones(len(h))*dump["quantities"]["massoftype"][0]
    vx = dump["blocks"][0]["data"]["vx"]
    vy = dump["blocks"][0]["data"]["vy"]
    vz = dump["blocks"][0]["data"]["vz"]
    u = dump["blocks"][0]["data"]["u"]

    filter = h > 0.0
    # dir   = '/STER/matse/Magritte/Lucy/'
    # filter[len(np.load(dir+'xtemp.npy')):] = False
    # filter[-1] = False
    # Format the data (select only data with positive smoothing length (h) and convert it to cgs-units
    x     = x                     [filter] * unit_dist          # position coordinates          [cm]
    y     = y                     [filter] * unit_dist
    z     = z                     [filter] * unit_dist
    mass  = mass                  [filter] * unit_mass          # mass of sph particles         [g]
    vx    = vx                    [filter] * unit_velocity / cgs.kms                   # velocity components           [cm/s]
    vy    = vy                    [filter] * unit_velocity / cgs.kms
    vz    = vz                    [filter] * unit_velocity / cgs.kms
    u     = u                     [filter] * unit_ergg          # specific internal density     [erg/g]
    h     = h                     [filter] * unit_dist          # smoothing length              [cm]
    rho   = pq.getRho(h, dump["quantities"]["hfact"], mass)     # density                       [g/cm^3]
    p     = pq.getPressure(rho, u, dump['quantities']['gamma']) # pressureure                   [Ba = 1e-1 Pa]


    if containsTau:
        tau  = dump["blocks"][0]["data"]["tau"][filter]         # optical depth
    if containsTauL:
        tauL  = dump["blocks"][0]["data"]["tau_lucy"][filter]   # Lucy optical depth
    if setup['isink_radiation'] > 1 and setup['iget_tdust'] == 0: temp = dump["blocks"][0]["data"]["Tdust"][filter]
    elif "temperature" in dump["blocks"][0]["data"]: temp = dump["blocks"][0]["data"]["temperature"][filter]
    else: temp = pq.getTemp(dump['quantities']['gamma'], setup['mu'], u) # temperature                [K]
    if bowenDust:
        Tdust = dump["blocks"][0]["data"]["Tdust"][filter]     # temperature                   [K]
    if setup['isink_radiation'] == 0: Gamma = 0.
    if setup['isink_radiation'] == 1: Gamma = np.ones_like(x)*setup['alpha_rad']
    if bowenDust:
        kappa = pq.getKappa(Tdust, setup['kappa_gas'], setup['bowen_delta'], setup['bowen_Tcond'], setup['bowen_kmax'])
        if containsTau:
            Gamma = pq.getGamma(kappa, lumAGB, massAGB, tau)
        else:
            Gamma = pq.getGamma(kappa, lumAGB, massAGB)
        if setup['isink_radiation'] == 3: Gamma = Gamma + setup['alpha_rad']

    cs    = pq.getSoundSpeed(p, rho, dump['quantities']['gamma'])            # speed of sound                [cm/s]
    vtan  = pq.getRadTanVelocity(x,y,vx,vy)                     # tangential velocity           [cm/s]
    r, phi, theta = gf.TransformToSpherical(x,y,z)              # sperical coordinates
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
            'Tgas'          : temp,                  # [K]
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
            'vz'            : vz,                    # [cm/s]
            'Gamma'         : Gamma
            }

    if not setup["single_star"]:
        data['posComp'    ] = posComp                # [cm]
        data['rComp'      ] = rComp                  # [cm]
        data['massComp'   ] = massComp               # [g]
        data['rHill'      ] = rHill                  # [cm]

    if setup['triple_star']:
        data["posComp_in" ] = posComp_in               # [cm]
        data['rComp_in'   ] = rComp_in                 # [cm]
        data['massComp_in'] = massComp_in              # [g]
    if setup['quadruple_star']:#hier nog verderdoen
        data["posComp_in" ] = posComp_in               # [cm]
        data["posComp_in2" ]= posComp_in2
        data['rComp_in'   ] = rComp_in                 # [cm]
        data['massComp_in'] = massComp_in              # [g]

    if containsTau:
        data["tau"        ] = tau

    if containsTauL:
        data["tauL"       ] = tauL

    if bowenDust:
        data["kappa"      ] = kappa                  # [cm^2/g]
        data["Tdust"      ] = Tdust


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
    temp  = dump['Tgas']
    if "tau"   in dump: tau   = dump['tau']
    if "tauL"  in dump: tauL   = dump['tauL']
    if "kappa" in dump: kappa = dump['kappa']
    if "Gamma" in dump: Gamma = dump['Gamma']
    if "Tdust" in dump: Tdust = dump['Tdust']
    h     = dump['h']
    r     = dump['r']

    filter = (r > factor * setup['sma_ini'] * cgs.au) & (r < bound * cgs.au)
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
    if "tauL"  in dump: tauL  = tauL  [filter]
    if "kappa" in dump: kappa = kappa [filter]             # cm^2/g
    if "Gamma" in dump: Gamma = Gamma [filter]
    if "Tdust" in dump: Tdust = Tdust [filter]
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
            'Tgas'          : temp,          # [K]
            'speed'         : speed,         # [cm/s]
            'mach'          : mach,
            'vtvv'          : vtvv,
            'r'             : r,             # [cm]
            'phi'           : phi,
            'theta'         : theta,
            'cs'            : cs             # [cm]
            }
    if "tau"   in dump: data["tau"]   = tau
    if "tauL"  in dump: data["tauL"]  = tauL
    if "kappa" in dump: data["kappa"] = kappa
    if "Gamma" in dump: data["Gamma"] = Gamma
    if "Tdust" in dump: data["Tdust"] = Tdust

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
    with os.scandir(runName) as entries:
        file_list = [entry.name for entry in entries if entry.is_file() and ("%s_" % userPrefix in entry.name) and (not (".asc" in entry.name or ".dat" in entry.name))]
    return np.sort(file_list)
