import os
import math                     as math
import numpy                    as np
# own scripts
import ConversionFactors_cgs    as cgs
import GeometricalFunctions     as gf
import Tools                    as tl
# ignore warnings
import warnings
warnings.filterwarnings("ignore")

'''
# TODO:
   - index error lijn 90 --> met de twist...
'''


'''
Calculate the terminal velocity of a given model by binning the speed 
of every sph particle according to the radius from the AGB star.
Since there is a spreading in speed, 3 different values for the 
terminal velocity are obtained (min, mean, max).
'''
def getTerminalVelocity(setup, dump):
    
    single_star = setup['single_star']
    if single_star == False:
        sma         = setup['sma_ini']
    outerBound  = int(round( setup['bound']  ))  
    r           = gf.getRadiusCoordinate(dump['position'],dump['posAGB'])/cgs.AU_cm()  # radius [AU] from AGB, not barycentre!
    

    #Prepare for binning.
    rmin = min( r )     # minimum radius in the data set [AU]
    rmax = max( r )     # maximum radius in the data set [AU]
    
    
    binned_speed, number_of_bins = tl.bin_data(dump['speed'],r)
    
    binned_speed_means = {}
    for i in range(number_of_bins):
        if len(binned_speed[i]) != 0:
            binned_speed_means[i] = [min(binned_speed[i]),np.mean(binned_speed[i]),max(binned_speed[i])]
    
    
    # collect per radius bin the minimum, mean and max
    binned_term_speed = {}
    binned_term_speed['min' ] = []
    binned_term_speed['mean'] = []
    binned_term_speed['max' ] = []
    for key in binned_speed_means:
        binned_term_speed['min' ].append(binned_speed_means[key][0])
        binned_term_speed['mean'].append(binned_speed_means[key][1])
        binned_term_speed['max' ].append(binned_speed_means[key][2])
    
    
    
    # cut off unphysical data
    '''
    If the speed is plot in function of r (3D to 1D), 
    we see that for some models an unphysical twist is 
    present in the speed. This mostly happens for models
    with a low mass companion, which does not affect the
    wind enough to influence the speed in the whole 
    regions of the model. Therefore, we look for the 
    position r at which this twist is present, and only
    take into account the data for r < r_twist.
    '''
    binned_term_speed['diff'] = []
    for i in range(len(binned_term_speed['mean'])-1):
        binned_term_speed['diff'].append(np.abs(binned_term_speed['mean'][i]-binned_term_speed['mean'][i+1]))
    differences = binned_term_speed['diff']
    
    twist = []
    #for i in range(5,len(differences)-5):
        #if (differences[i]+differences[i-1] <= differences[i+1]) and (differences[i]+differences[i-1] <= differences[i+2]) and (differences[i]+differences[i-1] <= differences[i+3]) and (differences[i]+differences[i-1] <= differences[i+4]) and (differences[i]+differences[i-1] >= differences[i-2]) and (differences[i]+differences[i-1] >= differences[i-3]) and (differences[i]+differences[i-1] >= differences[i-4]) and (differences[i]+differences[i-1] >= differences[i-5]): 
            #twist.append(i)
    
    if single_star == True:
        orbSep = 100
    
    if len(twist) != 0 and twist[0] > sma:
        #print('twist is there')
        index = tl.find_nearest(r,twist[0]) 
        r = r[index:]
        rmin = min(r)
        rmax = max(r)
        outerBound = int((round(r[index]))-1)
    else:
        index = -1
        

    ### TERMINAL VELOCITY ###
   
    if outerBound > number_of_bins:
        outerBound = number_of_bins
   
    # three options for the terminal velocity in [km/s]
    terminal_speed_min  = np.mean(binned_term_speed['min' ][round(outerBound-0.2*outerBound):int(outerBound)])
    terminal_speed_mean = np.mean(binned_term_speed['mean'][round(outerBound-0.2*outerBound):int(outerBound)])
    terminal_speed_max  = np.mean(binned_term_speed['max' ][round(outerBound-0.2*outerBound):int(outerBound)])
    
    terminal_speed = {'min': terminal_speed_min, 'mean': terminal_speed_mean, 'max': terminal_speed_max }
    
    if single_star == True:
        
        return terminal_speed, binned_term_speed
    
    if single_star == False:
        ### VELOCITY AT COMPANION ###
        
        wind_comp_min  = binned_term_speed['min' ][int(round(sma))]
        wind_comp_mean = binned_term_speed['mean'][int(round(sma))]
        wind_comp_max  = binned_term_speed['max' ][int(round(sma))]
        
        wind_comp = {'min': wind_comp_min, 'mean': wind_comp_mean, 'max': wind_comp_max}

    
        return terminal_speed, binned_term_speed, wind_comp, index
    
    


'''
Calculate the value of the parameter eta = speed / orbital velocity
for different speeds.
    1) terminal velocity
    2) velocity at the semi-major axis
    3) velocity at the semi-major axis of the corresponding single model
    4) vector sum of the velocity at the semi-major axis of the corresponding single model and the orbital velocity
'''
# NOTE Not yet possible to get (3) and (4), see later
def getEta_binary(setup, dump, sinkData, terminal_speed, wind_comp):
    sma = setup['sma_ini']
    
    if setup['ecc'] == 0:
        vOrb_AGB  = dump['v_orbAGB' ]
        vOrb_comp = dump['v_orbComp']

    elif setup['ecc'] >0:
        vOrb_AGB  = [min(sinkData['v_orbAGB_t' ]), np.mean(sinkData['v_orbAGB_t' ]), max(sinkData['v_orbAGB_t' ])]
        vOrb_comp = [min(sinkData['v_orbComp_t']), np.mean(sinkData['v_orbComp_t']), max(sinkData['v_orbComp_t'])]
  
    
    '''
    Uses the terminal velocity to compute eta
    '''
    eta1_min  = terminal_speed['min' ] / vOrb_comp
    eta1_mean = terminal_speed['mean'] / vOrb_comp
    eta1_max  = terminal_speed['max' ] / vOrb_comp
    
    
    eta1 = {'min': eta1_min, 'mean': eta1_mean, 'max': eta1_max}
    #print('')
    #print('eta1 = v_term/v_orb')
    #print(round(eta1_max , 2))
    #print(round(eta1_mean, 2))
    #print(round(eta1_min , 2))

    '''
    Uses the speed at the semi-major axis 
    '''
    
    #print('')
    #print('Wind speed at companion [km/s]:')
    #print(round(wind_speed_max , 2))
    #print(round(wind_speed_mean, 2))
    #print(round(wind_speed_min , 2))
    #print('')
    
    eta2_min  = wind_comp['min' ] / vOrb_comp
    eta2_mean = wind_comp['mean'] / vOrb_comp
    eta2_max  = wind_comp['max' ] / vOrb_comp
    
    eta2 = {'min': eta2_min, 'mean': eta2_mean, 'max': eta2_max}
    
    #print('eta2 = v_w(rcomp)/v_orb')
    #print(round(eta2_max , 2))
    #print(round(eta2_mean, 2))
    #print(round(eta2_min , 2))
    
    return eta1, eta2

    
'''
Give the speed at the location where the companion may be located for a single model.
Needed to calculate (3) and (4) of eta.
'''
# NOTE --> later
def getVelocity_single(binned_term_speed):
     if single_star == True:
        
        orbSep = [2.5,4,9,5]
        
        wind_speed_min  = {}
        wind_speed_mean = {}
        wind_speed_max  = {}
        wind_mass_at_comp = {}
        
        for i in range(len(orbSep)):
            # if the given semi-major axis is not an integer
            if not isinstance(orbSep[i], int):
                # three options for the terminal velocity in [km/s] 
                wind_speed_min[ orbSep[i] ] = ( np.mean([np.mean( binned_term_speed['min' ][int(orbSep[i])]),np.mean(binned_term_speed['min' ][int(orbSep[i]+1)]) ]) )
                wind_speed_mean[orbSep[i] ] = ( np.mean([np.mean( binned_term_speed['mean'][int(orbSep[i])]),np.mean(binned_term_speed['mean'][int(orbSep[i]+1)]) ]) )
                wind_speed_max[ orbSep[i] ] = ( np.mean([np.mean( binned_term_speed['max' ][int(orbSep[i])]),np.mean(binned_term_speed['max' ][int(orbSep[i]+1)]) ]) )
                #wind_mass_at_comp[orbSep[i]] = np.sum(binned_mass_hill[int(round(orbSep[i]))/Msun])
            
            # if the semi-major axis is an integer 
            if isinstance(orbSep[i], int):
                # three options for the terminal velocity in [km/s]
                wind_speed_min[ orbSep[i] ] = ( np.mean( binned_term_speed['min' ][int(orbSep[i])] ) )
                wind_speed_mean[orbSep[i] ] = ( np.mean( binned_term_speed['mean'][int(orbSep[i])] ) )
                wind_speed_max[ orbSep[i] ] = ( np.mean( binned_term_speed['max' ][int(orbSep[i])] ) )
                #wind_mass_at_comp[orbSep[i]] = np.sum(binned_mass_hill[int(round(orbSep[i]))]/Msun)
        
        #print('')
        #print('Wind speed at companion [km/s]:')
        #print(wind_speed_max )
        #print(wind_speed_mean)
        #print(wind_speed_min )
        #print('')
        #print('Wind mass at companion [Msun]:')
        #print(wind_mass_at_comp)
        
        
        
def getQp(setup, dump, wind_comp):
    
    z       = dump['position'].transpose()[2]
    mass    = dump['mass']
    rHill   = dump['rHill']
    sma     = setup['sma_ini'] * cgs.AU_cm()
    
    #print(len(mass))
    
    mass    = mass      [ z <  rHill ]
    r       = dump['r'] [ z <  rHill ]
    z       = z         [ z <  rHill ]
    #print(len(mass))
    mass    = mass      [ z > -rHill ]
    r       = r         [ z > -rHill ]
    z       = z         [ z > -rHill ]
    #print(len(mass))

    mass    = mass    [ r > sma-rHill ]  
    r       = r       [ r > sma-rHill ]
    #print(len(mass))
    mass    = mass    [ r < sma+rHill ]
    r       = r       [ r < sma+rHill ]
    #print(len(mass))
    
    massHill  = sum(mass)
    
    #print(mass)
    #print(massHill)
    
    
    wind_comp_mean   = np.mean([wind_comp['min'],wind_comp['max']])     # mean wind speed at the location of the compnion
    
    v_orb    = np.sqrt( cgs.G() * (dump['massComp']+dump['massAGB']) / (sma) ) 
    Qp       = ( (dump['massComp'] * v_orb)/(massHill * wind_comp_mean) ) 
        
    return Qp, massHill
    
    
    
def main_terminalVelocity(setup, dump, sinkData, outputloc, run):
    
    single_star = setup['single_star']
    
    if single_star == False:
        print('')
        print('(3)  Start calculations for terminal velocity...')
        print('')
        terminal_speed, binned_term_speed, wind_comp, index = getTerminalVelocity(setup, dump)
        print('')
        print('(4)  Start calculations for morphological parameters eta and Qp...')
        print('')
        eta1, eta2 = getEta_binary(setup, dump, sinkData, terminal_speed, wind_comp)
        Qp, massHill = getQp(setup, dump, wind_comp)
        
        
    if single_star == True:
        terminal_speed, binned_term_speed = getTerminalVelocity(setup, dump)
        

    title = outputloc+'vel_'+run+'.txt'
    with open (title,'w') as f:
        f.write('Model '+run+'\n')
        f.write('---------\n')
        f.write('\n')
        #f.write('Wind type:          '+str(data['windtype'])+'\n'  )
        f.write('-  Info about the model:\n')
        f.write('\n')
        f.write('Mass AGB star:             '+str(round(setup['massAGB_ini'],2))+' Msun\n'  )
        f.write('Initial wind velocity:     '+str(round(setup['v_ini'],2))+' km/s\n')
        f.write('Initial mass loss rate:    '+str(round(setup['Mdot'],2))+' Msun/yr\n')
        f.write('\n')
        #f.write('Companion:\n')
        if single_star == False:    
            f.write('Companion mass:            '+str(round(setup['massComp_ini'],2))+' Msun\n'  )
            f.write('Orbital separation:        '+str(round(setup['sma_ini'],2))+' au \n' ) 
        if single_star == True:
            f.write('Single star model, so no companion information.\n')
            f.write('\n')
        f.write('\n')
        f.write('\n')
        #f.write('\n')
        f.write('-  For the terminal velocity and eta, three different values are given.\n')
        f.write('   This is due to the spread on the speed of the model, because of the morphology.\n')
        f.write('   By the method of binning, a minimum, mean and maximum value is computed\n')
        f.write('\n')
        f.write('Terminal velocities [km/s]:'+'\n'      )
        f.write(str(round(terminal_speed['max' ]*cgs.cms_kms(), 3))+'\n' )
        f.write(str(round(terminal_speed['mean']*cgs.cms_kms(), 3))+'\n' )
        f.write(str(round(terminal_speed['min' ]*cgs.cms_kms(), 3))+'\n' )
        f.write('\n')
        if single_star == False:
            if setup['ecc'] == 0:
                f.write('eta1 = v_term/v_orb'+'\n')
                f.write(str(round(eta1['min' ], 2))+'\n')
                f.write(str(round(eta1['mean'], 2))+'\n')
                f.write(str(round(eta1['max' ], 2))+'\n')
                f.write('\n')
                f.write('Wind speed at companion [km/s]:\n')
                f.write(str(round(wind_comp['min']*cgs.cms_kms() , 2))+'\n')
                f.write(str(round(wind_comp['mean']*cgs.cms_kms(), 2))+'\n')
                f.write(str(round(wind_comp['max']*cgs.cms_kms() , 2))+'\n')
                f.write('\n')
                f.write('eta2 = v_w(rcomp)/v_orb\n')
                f.write(str(round(eta2['min' ], 2))+'\n')
                f.write(str(round(eta2['mean'], 2))+'\n')
                f.write(str(round(eta2['max' ], 2))+'\n')
            else:
                f.write('eta1 = v_term/v_orb'+'\n')
                f.write('Apastron:' +'\n')
                f.write(str(round(eta1['min' ][0], 2))+'\n')
                f.write(str(round(eta1['mean'][0], 2))+'\n')
                f.write(str(round(eta1['max' ][0], 2))+'\n')
                f.write('mean:' +'\n')
                f.write(str(round(eta1['min' ][1], 2))+'\n')
                f.write(str(round(eta1['mean'][1], 2))+'\n')
                f.write(str(round(eta1['max' ][1], 2))+'\n')
                f.write('periastron:' +'\n')
                f.write(str(round(eta1['min' ][2], 2))+'\n')
                f.write(str(round(eta1['mean'][2], 2))+'\n')
                f.write(str(round(eta1['max' ][2], 2))+'\n')
                f.write('\n')
                f.write('Wind speed at companion [km/s]:\n')
                f.write(str(round(wind_comp['min']*cgs.cms_kms() , 2))+'\n')
                f.write(str(round(wind_comp['mean']*cgs.cms_kms(), 2))+'\n')
                f.write(str(round(wind_comp['max']*cgs.cms_kms() , 2))+'\n')
                f.write('\n')
                f.write('eta2 = v_w(rcomp)/v_orb\n')
                f.write('Apastron:'+'\n')
                f.write(str(round(eta2['min' ][0], 2))+'\n')
                f.write(str(round(eta2['mean'][0], 2))+'\n')
                f.write(str(round(eta2['max' ][0], 2))+'\n')
                f.write('mean:'+'\n')
                f.write(str(round(eta2['min' ][1], 2))+'\n')
                f.write(str(round(eta2['mean'][1], 2))+'\n')
                f.write(str(round(eta2['max' ][1], 2))+'\n')
                f.write('Periastron:'+'\n')
                f.write(str(round(eta2['min' ][2], 2))+'\n')
                f.write(str(round(eta2['mean'][2], 2))+'\n')
                f.write(str(round(eta2['max' ][2], 2))+'\n')
            f.write('\n')
            f.write('Total mass the companion encounters at r = rcomp, within the Hill radius = Mwind [Msun]:\n')
            f.write(str(massHill/cgs.Msun_gram())+'\n')
            f.write('\n')
            f.write('Qp = (Mcomp*v_orb)/(Mwind*v_wind)\n')
            f.write('   (most of the time given as Q1 = 1e-6*Qp)\n')
            f.write(str(round(Qp)))
            #f.write(str(m_wind_hill/Msun)+' in bin\n')
        #if single_star == True:
            #f.write('Velocities [cm/s] at the different orbital separations:\n')
            #for key in wind_speed_min:
                #f.write('At '+str(key)+' au:\n')
                #f.write(str(round(wind_speed_max[ key], 2))+'\n')
                #f.write(str(round(wind_speed_mean[key], 2))+'\n')
                #f.write(str(round(wind_speed_min[ key], 2))+'\n')
                #f.write('\n')
                
        print('         The output of (3) terminal velocity and (4) morphological parameters is saved in a text file!')
        print('')
    
    
    
    
    
    
    
    
    
    
