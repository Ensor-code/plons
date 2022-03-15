import numpy                    as np
import math                     as math
import os
# own scripts
import PhysicalQuantities       as pq
import GeometricalFunctions     as gf
import ConversionFactors_cgs    as cgs



'''
Only load the prefix.in and prefix.setup files to get general information about the phantom model
      Suited for binary and single model
      
INPUT:
    - 'run' is the number of the run specifically           [str]
    - 'loc' is the directory where the model is located     [str]
    
RETURNS
    a dictionary containing the info from the setup files 
        (!! check units, they are not all in SI or cgs)

'''
def LoadSetup(run, loc, userSettingsDictionary):
    runName = os.path.join(loc, run)
    userPrefix = userSettingsDictionary["prefix"]

    # load the prefix.in & prefix.setup file
    setup = {}
    try:  
        with open(os.path.join(runName,'%s.setup'%userPrefix), 'r') as f:
            lines = f.readlines()
            for string in lines:
                line = string.split()
                if len(line) != 0:
                    if line[0] != '#':
                        stringName = line[0]

                        # Boolean
                        if stringName == 'icompanion_star':
                            stringName = 'single_star'
                            if int(line[2]) == 0: setup[stringName] = True
                            else: setup[stringName] = False

                        # Floats
                        else:
                            if stringName == 'primary_mass': stringName = 'massAGB_ini'
                            elif stringName == 'secondary_mass': stringName = 'massComp_ini'
                            elif stringName == 'semi_major_axis': stringName = 'sma_ini'
                            elif stringName == 'wind_gamma' or stringName == "temp_exponent": stringName = 'gamma'
                            elif stringName == 'eccentricity': stringName = 'ecc'
                            elif stringName == 'secondary_racc': stringName = 'rAccrComp'

                            setup[stringName] = float(line[2])

    except FileNotFoundError:
        print('')
        print(" ERROR: No %s.setup file found!"%userPrefix)
        print('')
        exit()

    try:
        with open(os.path.join(runName,'%s.in'%userPrefix), 'r') as f:
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
        print(" ERROR: No %s.setup file found!"%userPrefix)
        print('')
        exit()
    
    
    # Additional Parameters
    massAGB_ini = setup["massAGB_ini"]
    massComp_ini = setup["massComp_ini"]
    sma = setup["sma_ini"]
    v_ini = setup["v_ini"]
    period = pq.getPeriod(massAGB_ini * cgs.Msun_gram(), massComp_ini * cgs.Msun_gram(), sma)           # [s]
    v_orb = pq.getOrbitalVelocity(period, sma) * cgs.cms_kms()                                          # [km/s]
    Rcap = pq.getCaptureRadius(massComp_ini * cgs.Msun_gram(), v_ini / cgs.cms_kms()) / cgs.AU_cm()     # [au]

    setup["period"] = period
    setup["Rcap"] = Rcap
    setup["v_orb"] = v_orb

    
    return setup
