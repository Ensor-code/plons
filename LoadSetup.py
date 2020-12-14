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
    try:  
        with open(os.path.join(runName,'%s.setup'%userPrefix), 'r') as f:
            lines = f.readlines()  
            setup = []
            for i in range(len(lines)):
                setup.append(lines[i].split())
    except FileNotFoundError:
        print('')
        print(" ERROR: No %s.setup file found!"%userPrefix)
        print('')
        exit()

    try:
        with open(os.path.join(runName,'%s.in'%userPrefix), 'r') as f:
            lines = f.readlines()  
            infile = []
            for i in range(len(lines)):
                infile.append(lines[i].split())
    except FileNotFoundError:
        print('')
        print(" ERROR: No %s.setup file found!"%userPrefix)
        print('')
        exit()
    
    
    # Get specifics about the model
    is_single   = int(setup[6][2])
    if is_single == 0:
        single_star = True
    if is_single == 1:
        single_star = False


    tmax          = float(infile[9][2])         # last dump in code units
    nfulldump     = float(infile[16][2])        # Every time full dump is produced
    v_ini         = float(infile[58][2])        # initial wind velocity                         [km/s]
    mu            = float(infile[38][2])        # mean molecular weight of the gas
    Mdot          = float(infile[62][2])        # mass loss rate                                [M_sun/yr]
    # !! phantom error: the mass of the sph particles is wrong with a factor 2, to 'make up' for this, the mass loss rate Mdot can be multiplied by 2
    bound         = float(infile[69][2])        # outer boundary of the model                   [au] 
    massAGB_ini   = float(setup[1][2])          # initial mass of the AGB star                  [Msun]
    
    
    if single_star == False:
        massComp_ini  = float(setup[7][2])      # initial mass of the companion                 [Msun]
        rAccrComp     = float(setup[8][2])      # accretion radius of companion                 [au]
        sma           = float(setup[12][2])     # semi-major axis (initial orbital separation)  [au]
        gamma         = float(setup[15][2])     # adiabatic constant        
        ecc           = float(setup[13][2])     # eccentricity
        period        = pq.getPeriod(massAGB_ini*cgs.Msun_gram(),massComp_ini*cgs.Msun_gram(),sma)              # [s]
        v_orb         = pq.getOrbitalVelocity(period, sma)*cgs.cms_kms()                                        # [km/s]
        Rcap          = pq.getCaptureRadius(massComp_ini*cgs.Msun_gram(), v_ini/cgs.cms_kms())/cgs.AU_cm()      # [au]                                   
        
    if single_star == True:
        gamma         = float(setup[8][2])     # adiabatic constant   
    

    
    # output
    if single_star == False:
        setup = {'v_ini'        : v_ini,                # [km/s]
                'massComp_ini'  : massComp_ini,         # [Msun]
                'sma_ini'       : sma,                  # [au]
                'single_star'   : single_star,          # boolean
                'bound'         : bound,                # [au]
                'v_orb'         : v_orb,                # [km/s]
                'tmax'          : tmax,
                'nfulldump'     : nfulldump,
                'mu'            : mu,
                'Mdot'          : Mdot,                 # [Msun/yr]
                'massAGB_ini'   : massAGB_ini,          # [Msun]
                'ecc'           : ecc,
                'gamma'         : gamma,
                'rAccrComp'     : rAccrComp,            # [au]
                'period_ini'    : period,               # [s]
                'Rcap'          : Rcap                  # [AU]
                                
                }
        
    if single_star == True:
        setup = {'v_ini'        : v_ini,               # [km/s]
                'single_star'   : single_star,
                'bound'         : bound,                # [au]
                'mu'            : mu,
                'Mdot'          : Mdot,                 # [Msun/yr]
                'massAGB_ini'   : massAGB_ini,          # [Msun]
                'tmax'          : tmax,
                'gamma'         : gamma
                }

    
    return setup
