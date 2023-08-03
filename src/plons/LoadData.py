import math                     as math
import numpy                    as np
import sys

#import plons scripts
import plons.LoadDump                 as dmp
import plons.LoadSink                 as snk
import plons.PhysicalQuantities       as pq
import plons.ConversionFactors_cgs    as cgs
import plons.GeometricalFunctions     as gf

import os


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
def LoadData_cgs(run, loc, userSettingsDictionary, bound = None, factor = -1, number = -1):
    
    dir       = os.path.join(loc, run)
    setup     = LoadSetup(dir, userSettingsDictionary["prefix"])

    # Pick either last dump file or user chosen file
    if number == -1: index = dmp.findLastFullDumpIndex(userSettingsDictionary["prefix"], setup, dir)
    else: index = number
    fileName       = dir+'/{0:s}_{1:05d}'.format(userSettingsDictionary["prefix"], index)
    dumpData  = dmp.LoadDump_cgs(fileName, setup, userSettingsDictionary["hard_path_to_phantom"])

    sinkData  = snk.LoadSink_cgs(dir, userSettingsDictionary['prefix'], setup["icompanion_star"])
    if bound == None:
        bound = setup['bound']
    if factor > 0:
        outerData = dmp.LoadDump_outer_cgs(factor, bound, setup, dumpData)
    else: outerData = None

    # save the final specifics of the AGB star to dumpData
    dumpData['maccrAGB' ] = sinkData['maccrAGB'   ][-1]

    # save the final specifics of the companion to dumpData
    if not setup['single_star']:
        dumpData['maccrComp'] = sinkData['maccrComp'  ][-1]
        dumpData['v_orbAGB' ] = sinkData['v_orbAGB_t' ][-1]
        dumpData['v_orbComp'] = sinkData['v_orbComp_t'][-1]
        dumpData['sma_fi'   ] = dumpData['rAGB'] + dumpData['rComp']        # [cm]

    if setup['triple_star']:
        dumpData['maccrComp_in'] = sinkData['maccrComp_in'  ][-1]
        dumpData['v_orbComp_in'] = sinkData['v_orbComp_in_t'][-1]

    return setup, dumpData, sinkData, outerData


'''
Only load the prefix.in and prefix.setup files to get general information about the phantom model
      Suited for binary and single model
      
INPUT:
    - 'location' is the directory where the model is located     [str]
    
RETURNS
    a dictionary containing the info from the setup files 
        (!! check units, they are not all in SI or cgs)

'''
def LoadSetup(dir, prefix):
    # load the prefix.in & prefix.setup file
    setup = {}
    try:  
        with open(os.path.join(dir,'%s.setup'%prefix), 'r') as f:
            lines = f.readlines()
            for string in lines:
                line = string.split()
                if len(line) != 0:
                    if line[0] != '#':
                        stringName = line[0]

                        # Floats
                        if stringName == 'primary_mass': stringName = 'massAGB_ini'
                        elif stringName == 'secondary_mass': stringName = 'massComp_ini'
                        elif stringName == 'semi_major_axis': stringName = 'sma_ini'
                        elif stringName == 'binary2_a' : stringName = 'sma_in_ini'
                        elif stringName == 'wind_gamma' or stringName == "temp_exponent": stringName = 'gamma'
                        elif stringName == 'eccentricity': stringName = 'ecc'
                        elif stringName == 'binary2_e' : stringName = 'ecc_in'                            
                        elif stringName == 'secondary_racc': stringName = 'rAccrComp'
                        elif stringName == 'accr2b' : stringName = 'rAccrComp_in'
                        elif stringName == 'racc2b' : stringName = 'rAccrComp_in'

                        setup[stringName] = float(line[2])
                        
                        if stringName == 'q2': 
                            stringName = 'massComp_in_ini'    
                            setup[stringName] = float(line[2])*setup['massAGB_ini']

                        # Boolean
                        if stringName == 'icompanion_star':
                            stringName = 'single_star'
                            if int(line[2]) == 0: setup[stringName] = True
                            else: setup[stringName] = False
                            stringName = 'triple_star'
                            if int(line[2]) == 2: setup[stringName] = True
                            else: setup[stringName] = False

    except FileNotFoundError:
        print('')
        print(" ERROR: No %s.setup file found!"%prefix)
        print('')
        exit()

    try:
        with open(os.path.join(dir,'%s.in'%prefix), 'r') as f:
            lines = f.readlines()
            for string in lines:
                line = string.split()
                if len(line) != 0:
                    if line[0] != '#':
                        stringName = line[0]

                        # Strings
                        if stringName == 'logfile': setup[stringName] = str(line[2])
                        elif stringName == 'dumpfile': setup[stringName] = str(line[2])
                        elif stringName == 'twallmax': setup[stringName] = str(line[2])
                        elif stringName == 'dtwallmax': setup[stringName] = str(line[2])

                        # Floats
                        else:
                            if stringName == 'wind_velocity': stringName = 'v_ini'
                            elif stringName == 'outer_boundary': stringName = 'bound'
                            elif stringName == 'wind_mass_rate': stringName = 'Mdot'
                            if line[2]=='F': line[2] = 0
                            setup[stringName] = float(line[2])

    except FileNotFoundError:
        print('')
        print(" ERROR: No %s.setup file found!"%prefix)
        print('')
        exit()
    
    
    # Additional Parameters
    massAGB_ini = setup["massAGB_ini"]
    v_ini = setup["v_ini"]
    if setup["single_star"] == False:
        massComp_ini = setup["massComp_ini"]
        sma = setup["sma_ini"]
        period = pq.getPeriod(massAGB_ini * cgs.Msun, massComp_ini * cgs.Msun, sma)           # [s]
        #v_orb = pq.getOrbitalVelocity(period, sma) * cgs.cms_kms()                           # [km/s]
        Rcap = pq.getCaptureRadius(massComp_ini * cgs.Msun, v_ini * cgs.kms) / cgs.au         # [au]
        if setup['triple_star']==True:
            massComp_in_ini = setup["massComp_in_ini"]        
            sma_in = setup["sma_in_ini"]
            period = pq.getPeriod((massAGB_ini+massComp_in_ini) * cgs.Msun, massComp_ini * cgs.Msun, sma)  # [s]
            period_in = pq.getPeriod(massAGB_ini * cgs.Msun, massComp_in_ini * cgs.Msun, sma_in)           # [s]
            Rcap_in = pq.getCaptureRadius(massComp_in_ini * cgs.Msun, v_ini * cgs.kms) / cgs.au            # [au]
            setup["period_in"] = period_in
            setup["Rcap_in"] = Rcap_in

        setup["period"] = period
        setup["Rcap"] = Rcap
    
    return setup

'''
Loads the final full dump of a phantom model, given the number, in cgs-units 
    
INPUT:
    - 'dir'   is the directory where the model is located     [str]
    - 'setup' is the setup data                               [dict]
    
RETURNS
    a dictionary containing the data from the dump (all units in cgs)
'''
def LoadDump_cgs(fileName, setup, phantom_dir):

    # Load read_dump script to read PHANTOM dump files
    sys.path.append(phantom_dir+"/scripts")
    from readPhantomDump import read_dump

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
    xAGB     = dump["blocks"][1]["data"]["x"][0] * unit_dist
    yAGB     = dump["blocks"][1]["data"]["y"][0] * unit_dist
    zAGB     = dump["blocks"][1]["data"]["z"][0] * unit_dist
    posAGB   = np.array([xAGB, yAGB, zAGB])
    rAGB     = gf.calc_r(xAGB, yAGB, zAGB)   
    massAGB  = dump["blocks"][1]["data"]["m"][0] * unit_mass
    lumAGB   = dump["blocks"][1]["data"]["lum"][0] * unit_energ/unit_time
    vxAGB    = dump["blocks"][1]["data"]["vx"][0] * unit_velocity
    vyAGB    = dump["blocks"][1]["data"]["vy"][0] * unit_velocity
    vzAGB    = dump["blocks"][1]["data"]["vz"][0] * unit_velocity
    velAGB   = np.array([vxAGB, vyAGB, vzAGB])


    # companion
    if not setup["single_star"]:
        xComp    = dump["blocks"][1]["data"]["x"][1] * unit_dist
        yComp    = dump["blocks"][1]["data"]["y"][1] * unit_dist
        zComp    = dump["blocks"][1]["data"]["z"][1] * unit_dist
        posComp  = np.array([xComp, yComp, zComp])
        rComp    = gf.calc_r(xComp, yComp, zComp)
        massComp = dump["blocks"][1]["data"]["m"][1] * unit_mass
        rHill = pq.getRHill(abs(rComp+rAGB),massComp,massAGB)
        vxComp   = dump["blocks"][1]["data"]["vx"][1] * unit_velocity
        vyComp   = dump["blocks"][1]["data"]["vy"][1] * unit_velocity
        vzComp   = dump["blocks"][1]["data"]["vz"][1] * unit_velocity
        velComp  = np.array([vxComp,vyComp,vzComp])


        if setup['triple_star']:
            # inner companion
            xComp_in    = dump["blocks"][1]["data"]["x"][2] * unit_dist
            yComp_in    = dump["blocks"][1]["data"]["y"][2] * unit_dist
            zComp_in    = dump["blocks"][1]["data"]["z"][2] * unit_dist
            posComp_in  = np.array([xComp_in, yComp_in, zComp_in])
            rComp_in    = gf.calc_r(xComp_in, yComp_in, zComp_in)
            massComp_in = dump["blocks"][1]["data"]["m"][2] * unit_mass
            vxComp_in   = dump["blocks"][1]["data"]["vx"][2] * unit_velocity
            vyComp_in   = dump["blocks"][1]["data"]["vy"][2] * unit_velocity
            vzComp_in   = dump["blocks"][1]["data"]["vz"][2] * unit_velocity
            velComp_in  = np.array([vxComp_in,vyComp_in,vzComp_in])


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
    vx    = vx                    [filter] * unit_velocity / cgs.kms                   # velocity components           [cm/s]  !!---> no, in km/s??
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
            'r'             : r,                     # [cm]
            'phi'           : phi,     
            'theta'         : theta,     
            'cs'            : cs,                    # [cm]
            'posAGB'        : posAGB,                # [cm]
            'rAGB'          : rAGB,                  # [cm]
            'massAGB'       : massAGB,               # [g]
            'velAGB'        : velAGB,
            'vx'            : vx,                    # [cm/s]
            'vy'            : vy,                    # [cm/s]
            'vz'            : vz,                    # [cm/s]
            'Gamma'         : Gamma
            }

    if not setup["single_star"]:
        data['posComp'    ] = posComp                # [cm]
        data['rComp'      ] = rComp                  # [cm]
        data['massComp'   ] = massComp               # [g]
        data['velComp'    ] = velComp
        data['rHill'      ] = rHill                  # [cm]

    if setup['triple_star']:
        data["posComp_in" ] = posComp_in               # [cm]
        data['rComp_in'   ] = rComp_in                 # [cm]
        data['massComp_in'] = massComp_in              # [g]
        data['velComp_in' ] = velComp_in
        
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
    
    p     = pq.getPressure(rho, u, setup['gamma'])              # pressure                      [Ba = 1e-1 Pa]
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
    lastFullDumpIndex = findLastFullDumpIndex(userPrefix, setup, runName)
    fullDumpLists = np.arange(0, lastFullDumpIndex + 1, setup['nfulldump'], dtype=int)
    return fullDumpLists

def sortedDumpList(userPrefix, runName):
    return np.sort(list(filter(lambda x: ("%s_"%userPrefix in x) and (not (".asc" in x or ".dat" in x)), os.listdir(runName))))
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
def LoadSink_cgs(dir, userPrefix, icompanion_star = 0):
    
    fileName_sink1 = os.path.join(dir, str('%sSink0001N0'%userPrefix+'1.ev'))
    # companion
    if icompanion_star > 0:
        fileName_sink2 = os.path.join(dir, str('%sSink0002N0'%userPrefix+'1.ev'))
        # inner companion
        if icompanion_star == 2:
            fileName_sink3 = os.path.join(dir, str('%sSink0003N0'%userPrefix+'1.ev'))

    try:
    # to calculate period, we need masses and sma, so coordinates
        (t1, x1,y1,z1, mass1, vx1,vy1,vz1, maccr1) = np.loadtxt(fileName_sink1, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
        if icompanion_star > 0:
            n_file = len(t1)
            (t2, x2,y2,z2, mass2, vx2,vy2,vz2, maccr2) = np.loadtxt(fileName_sink2, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)[:, :n_file]
            if icompanion_star == 2:
                (t3, x3, y3, z3, mass3, vx3, vy3, vz3, maccr3) = np.loadtxt(fileName_sink3, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)[:, :n_file]
    except OSError:
        print(' ERROR: No sink files found for this model in the current directory!')


    numberOfevFiles = findLastWindSinkIndex(dir,userPrefix)
    for n in range(2,numberOfevFiles+1):
        fileName_sink1 = os.path.join(dir, str('%sSink0001N0'%userPrefix+str(n)+'.ev'))
        if icompanion_star > 0:
            fileName_sink2 = os.path.join(dir, str('%sSink0002N0'%userPrefix+str(n)+'.ev'))
            if icompanion_star == 2:
                fileName_sink3 = os.path.join(dir, str('%sSink0003N0'%userPrefix+str(n)+'.ev'))

        try:
        # to calculate period, we need masses and sma, so coordinates
            (t1e, x1e,y1e,z1e, mass1e, vx1e,vy1e,vz1e, maccr1e) = np.loadtxt(fileName_sink1, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)
            if icompanion_star > 0:
                n_file = len(t1e)
                (t2e, x2e,y2e,z2e, mass2e, vx2e,vy2e,vz2e, maccr2e) = np.loadtxt(fileName_sink2, skiprows=1, usecols=(0,1,2,3,4,5,6,7,11), unpack=True)[:, :n_file]
                if icompanion_star == 2:
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

        if icompanion_star > 0:
            t2     = np.append(t2,t2e)
            x2     = np.append(x2,x2e)
            y2     = np.append(y2,y2e)
            z2     = np.append(z2,z2e)
            mass2  = np.append(mass2, mass2e)
            vx2    = np.append(vx2, vx2e)
            vy2    = np.append(vy2, vy2e)
            vz2    = np.append(vz2, vz2e)
            maccr2 = np.append(maccr2, maccr2e)

            if icompanion_star == 2:
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

    # companion
    if icompanion_star > 0:
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
        
        rHill           = pq.getRHill( abs(r1 + r2), mass2, mass1           )         # [cm]

        if icompanion_star == 2:
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
    if icompanion_star > 0:
        orbitalVel_comp = np.sqrt(np.transpose(velocity2)[0]**2+np.transpose(velocity2)[1]**2)
        if icompanion_star == 2:
            orbitalVel_comp_in = np.sqrt(np.transpose(velocity3)[0]**2+np.transpose(velocity3)[1]**2)




    #print('test orbital velocity:')
    #print(np.shape(np.linalg.norm(np.transpose(velocity1))),np.shape(orbitalVel_AGB))
    #print(np.nanmean(np.linalg.norm(np.transpose(velocity1))),np.nanmean(orbitalVel_AGB))
    #print(np.array(np.linalg.norm(velocity1))[0:20],orbitalVel_AGB[0:20])

    # output
    #    "_t" stands for the fact that these values are in function of the evolution time, not from the last dump in function of location
    data = {'posAGB'      : position1,              # [cm]
            'velAGB'      : velocity1,              # [cm/s]
            'massAGB'     : mass1,                  # [gram]
            'maccrAGB'    : maccr1,                 # [gram]
            'rAGB'        : r1,                     # [cm]
            }
    if icompanion_star > 0:
        data['posComp'     ] = position2
        data['velComp'     ] = velocity2
        data['massComp'    ] = mass2
        data['maccrComp'   ] = maccr2
        data['rComp'       ] = r2
        data['time'        ] = t1                     # [yrs]
        data['v_orbAGB_t'  ] = orbitalVel_AGB         # [cm/s]
        data['v_orbComp_t' ] = orbitalVel_comp        # [cm/s]
        data['rHill_t'     ] = rHill                  # [cm]
        if icompanion_star == 2:
            data['posComp_in'    ] = position3,
            data['velComp_in'    ] = velocity3,
            data['massComp_in'   ] = mass3,
            data['maccrComp_in'  ] = maccr3,
            data['rComp_in'      ] = r3,
            data['rHill_in_t'    ] = rHill_in,
            data['period_t_in'   ] = period_in,
            data['v_orbComp_in_t'] = orbitalVel_comp_in

    return data


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
    return np.sort(list(filter(lambda x: ("%s"%userPrefix in x) and (".ev" in x), os.listdir(runName))))