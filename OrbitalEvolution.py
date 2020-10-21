# Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import matplotlib.lines         as mlines
import math                     as math
import os
from matplotlib                 import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size' :   12})
rcParams.update({'figure.dpi': 200})
#rc('font', family='serif')
#rc('text', usetex=True)

# import own scripts
import ConversionFactors_cgs    as cgs
import PhysicalQuantities       as pq

'''
Calculates:
    - times at which orbit is in apastron and periastron          (timeApa, timePer)
    - orbital separation at all apastron and periastron passages  (apastronOS, periastronOS)

returns first maximum, last maximum, first minimum, last minumimum
uses sinkData = data
'''
def orbSep_apa_per(data, setup):

    period = setup['period_ini'] * cgs.sec_year()
    time   = data['time'       ]
    ecc    = setup['ecc'       ]
    orbSep = data['rComp'      ] + data['rAGB' ]
    # print('time goes from ', min(time),' to ', max(time))
    # print('period is ', period)
    # print('the amount of periods is max(time)/period: ', math.ceil(max(time)/period))

    apastronOS   = []
    periastronOS = []  
    timeApa      = []
    timePer      = []
    i=0

    while i<math.ceil(max(time)/period):
        # Model starts in apastron passage, so apastron passages are at i*period for i = 0, 1, ...
        tApa     = i*period
        # The exact apastron is defined as where the orbital separation is maximal
        # Find the index of that maximum
        # time will be between tminApa and tmaxApa
        tminApa  = tApa - 0.5
        tmaxApa  = tApa + 0.5
        # index will be between Iapa1 and Iapa2
        Iapa1    = np.abs(time - tminApa).argmin()
        Iapa2    = np.abs(time - tmaxApa).argmin()
        indexApa = (orbSep[Iapa1:Iapa2]).argmax() + Iapa1

        # Periastron passages are at period*(i+0.5) for i = 0, 1, ...
        tPer     = (i+0.5)*period
        # The exact periastron is defined as where the orbital separation is minimal
        # Find the index of that miniimum
        # time will be between tminPer and tmaxPer
        tminPer  = tPer - 0.5
        tmaxPer  = tPer + 0.5
        # index will be between Iper1 and Iper2
        Iper1    = np.abs(time - tminPer).argmin()
        Iper2    = np.abs(time - tmaxPer).argmin()
        indexPer = (orbSep[Iper1:Iper2].argmin()) + Iper1
        
        #Make arrays with all apastron and periastron orbital separations
        apastronOS.append(orbSep[indexApa])
        periastronOS.append(orbSep[indexPer])
        #Make arrays with all apastron and periastron times
        timeApa.append(time[indexApa])#*timeUnit)
        timePer.append(time[indexPer])#*timeUnit)
        
        i = i+1
        
    return apastronOS, periastronOS, timeApa, timePer



'''
Visualises the orbit of the system
Uses sinkData = data
'''
def plot_orbit(data, setup, ax): 
    t  = data['time' ]
    rC = data['rComp']
    rA = data['rAGB' ]

    xc = data['posComp'].transpose()[0] 
    xa = data['posAGB' ].transpose()[0]
    yc = data['posComp'].transpose()[1]
    ya = data['posAGB' ].transpose()[1]
    zc = data['posComp'].transpose()[2]
    za = data['posAGB' ].transpose()[2]
    
    M1   = data['massComp']
    M2   = data['massAGB' ]
    xCOM = (M1*xc + M2 * xa)/(M1+M2)
    yCOM = (M1*yc + M2 * ya)/(M1+M2)
    rCOM = np.sqrt(xCOM**2 + yCOM**2)
    
    ax.axis('equal')
    ax.set_ylim(-4,4)
    
    ax.plot(xc,yc, c = 'r', label = 'companion')
    ax.plot(xa,ya, c = 'b', label = 'AGB'      )
    ax.set_ylabel('$y$ [cm]',fontsize = 15)
    ax.set_xlabel('$x$ [cm]',fontsize = 15)
    ax.plot(0,0,'*', color = 'k', label = 'COM')
    ax.set_title('$e = $'+str(setup['ecc'])+', $a = $'+ str(setup['sma_ini'])+' AU, $q = $'+ str(setup['massAGB_ini']/setup['massComp_ini']), fontsize = 15)

    ax.tick_params(labelsize=12)
    

'''
Calculates:
    - semi-major axis of each orbit
    - time at middle of each orbit
'''
def calcSMA(OrbSep_a_p):
    i = OrbSep_a_p
    sma  = []
    tOrbit = []
    for j in range(len(i[0])):
        # sma is (apastronOS+periastronOS)/2
        sma.append((i[0][j]+i[1][j])/2)
        # tOrbit is mean time of timeApa and timePer
        tOrbit.append((i[2][j]+i[3][j])/2)
    return (sma,tOrbit)



'''
Makes plots of the evolution of the semi-major axis and of the orbital separation at the apastron and periastron passages
Calculates the change in eccentricity and semi-major axis 
Returns text file with data to make the evolution plot and with the calculated delta(e) and delta(a)
'''
def plotChangeOrbSep(info, sinkData, setup, run, loc):#, ylabel, unit, name, title):
    
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    
    # legend
    apastron      = mlines.Line2D([],[], color = 'k', linestyle = 'None'  , marker = '$a$' ,label = 'Apastron', markersize = 8)
    periastron    = mlines.Line2D([],[], color = 'k', linestyle = 'None'  , marker = '$p$', label = 'Periastron', markersize =8)
    sma           = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'sma')
    handles1      = [apastron,sma, periastron]

    
    fig.set_size_inches(7, 5)
    c = 'k'
    
    # the data to plot are the orbital separations at apastron and periastron passages
    toPlot        = info['OrbSep_a_p']
    #Calculate the sma and mean time at every orbit
    [sma,orbit_t] = calcSMA(info['OrbSep_a_p'])
    
    #t_total is the final timestep
    t_total     = max(sinkData['time']) #in years

    #The change in sma per year is: delta a = (delta r_apa + delta r_per) / (2*time)
    delta_a              = (((toPlot[0]-toPlot[0][0] + (toPlot[1]-toPlot[1][0]) )/2)[-1]) / t_total
    ratio_delta_a_per_yr = (delta_a/(setup['sma_ini']*cgs.AU_cm()))
    
    #The change in ecc per year is: delta e = (delta r_apa - delta r_per) / (2*(a+delta a) *time)
    delta_e              = (((toPlot[0]-toPlot[0][0] - (toPlot[1]-toPlot[1][0]) )/(2*setup['sma_ini']*cgs.AU_cm() + 2*delta_a) )[-1])/t_total
    ratio_delta_e_per_yr = (delta_e/setup['ecc'] ) 

    # Plot the change in orbital separation at apastron and periastron passages for each orbit
    ax.plot(orbit_t,toPlot[0]-toPlot[0][0],color = c, marker = '$a$', linestyle = 'dashed', markersize = 10)
    ax.plot(orbit_t,toPlot[1]-toPlot[1][0],color = c, marker = '$p$', linestyle = 'dashed', markersize = 10)
    # Change markers to 'a' and 'p'
    for i in range(len(orbit_t)):
        ax.plot(orbit_t[i],toPlot[0][i]-toPlot[0][0],color = 'white', marker = 'o'  , markersize = 12)
        ax.plot(orbit_t[i],toPlot[0][i]-toPlot[0][0],color = c      , marker = '$a$', markersize = 11)
        ax.plot(orbit_t[i],toPlot[1][i]-toPlot[1][0],color = 'white', marker = 'o', markersize = 12)
        ax.plot(orbit_t[i],toPlot[1][i]-toPlot[1][0],color = c      , marker = '$p$', markersize = 11)
    #Plot the change in sma at each orbit
    ax.plot(orbit_t, sma-sma[0], color = c, linestyle = 'dotted', markersize = 10)
    
    #change xticks to orbit numbers
    ax.tick_params(labelsize=12)     
    period = setup['period_ini'] * cgs.sec_year()
    tickorbit = []
    for p in range(1, int(setup['tmax']/period)+1,1):
        tickorbit.append(str(p))
        
    plt.setp(ax, xticks= orbit_t, xticklabels=tickorbit)
    
    ax.set_xlabel('Orbit', fontsize = 18)
    ax.set_ylabel('$\Delta$Orb sep [cm]', fontsize = 16)
    ax.set_title('Orbital evolution')
    ax.legend(handles = handles1, fontsize = 12)#, loc = 'lower left')
    fig.tight_layout()

    plt.savefig(loc+str(run)+'_evolution_OrbitalSeparation')
    
    #Write text file with usefull info
    title = loc+str(run)+'_data_OrbitalEvolution.txt'
    with open (title,'a') as f:
        f.write('\n')
        if setup['ecc']>0:
            f.write('The change in eccentricity: delta_e/e =         '+ str(round(ratio_delta_e_per_yr,9))+ '/yr'+'\n')
        else:
            f.write('The change in eccentricity: delta_e   =         '+ str(round(delta_e,9))+ '/yr'+'\n')
        f.write('The change in sma:          delta_a/a =         '+ str(round(ratio_delta_a_per_yr,9))+ '/yr'+'\n')
        f.write('\n')
        f.write('To plot the change in eccentricity and change in sma, use the following: '+ '\n')
        f.write('\n')
        f.write('Times of middle of each orbit, to plot on x axis [yrs]: '+ '\n')
        f.write(str(orbit_t)+'\n' )
        f.write('Orbital separation at apastron passages [cm] : '+ '\n')
        f.write(str(toPlot[0])+ '\n')
        f.write('Orbital separation at periastron passages [cm]: '+ '\n')
        f.write(str(toPlot[1])+ '\n')
        f.write('Semi major axis at times of periastron passage [cm]: '+ '\n')
        f.write(str(sma)+ ' \n')



'''
Makes plot of the evolution of mass accretion by the companion
Returns text file with data to make the evolution plot 
'''

def plotMassAccr(setup, sinkData, run, loc):
  # make plot of the mass accretion evolution, very interesting to plot!
    fig = plt.figure(figsize=(8, 5))
    #legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    handles_ap    = [apaLine, perLine]

    #plot the accreted mass in function of time
    plt.plot(sinkData['time'],  sinkData['maccrComp'], color = 'royalblue', linestyle = 'solid')
    
    #Plot vertical lines indicating where there are apastron and periastron passages
    period = setup['period_ini'] * cgs.sec_year()
    j = period/2  # First periastron
    i = 0         # Start at apastron

    maxi = max(sinkData['maccrComp'])
    mini = min(sinkData['maccrComp'])
    for orbit in range(0, int(sinkData['time'][-1]/period)+1):
        plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
        i = i+period
        j = j+period

    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
   
    ax = plt.subplot(111)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Accreted mass [g]', fontsize = 16)

    plt.title('Total accreted mass by the companion', fontsize = 18)
    plt.legend(handles = handles_ap)
    fig.tight_layout()
    plt.savefig(loc+str(run)+'_evolution_Maccr_companion')



'''
Makes plot of the evolution of the orbital velocity of the AGB star and companion 
'''

def plotOrbVel(sinkData, run, loc):
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.suptitle('Orbital velocities (model'+str(run)+')', fontsize = 15)

    ax.plot(sinkData['time'], sinkData['v_orbComp_t'], label ='companion')
    ax.plot(sinkData['time'], sinkData['v_orbAGB_t'], label = 'AGB')
    ax.set_ylabel('$v_{orb}$ [cm/s]', fontsize = 12)

    ax.set_xlabel('time[yrs]', fontsize = 10)

    ax.tick_params(labelsize=10)
          
    plt.legend()
    fig.tight_layout()




'''
Makes plot of the evolution of the orbital radii of the AGB star and companion 
'''
def plotOrbRad(sinkData, run, loc):
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})

    ax.plot(sinkData['time'], sinkData['rComp'], label= 'r comp')
    ax.set_ylabel('$r$ [cm]', fontsize = 12)
    ax.set_title('r comp, r AGB, orb sep (model'+str(run)+')', fontsize = 15)
       
    ax.plot(sinkData['time'], sinkData['rAGB'], label ='r AGB')
    ax.plot(sinkData['time'], sinkData['rAGB']+sinkData['rComp'], label = 'Orb sep')
    ax.set_xlabel('time [yrs]', fontsize = 14)
    ax.tick_params(labelsize=10)

    plt.legend()
    fig.tight_layout()
    plt.savefig(loc+str(run)+'_evolution_rComp_rAGB_orbSep')



'''
Main function executing calculations about orbital evolution
'''

def orbEv_main(run,loc, sinkData, setup):
    print('')
    if setup['single_star'] == True:
        print('(6)  A single model has no orbit, and thereby no orbital evolution.')
        print('     The orbital evolution part is therefore skipped.')
        
    else:
        print('(6)  Start calculations for orbital evolution...')

        # calculate total mass accreted by companion and AGB star, total mass lost by the AGB star and ratio mass accreted to mass lost
        # save it in dictionary 'info'
        info = {}
        info['TotMaC']         = sinkData['maccrComp'][-1] # the total mass accreted by the companion
        info['TotMaA']         = sinkData['maccrAGB' ][-1]  # the total mass accreted by the AGB
        info['MassLostAGB']    = sinkData['massAGB'  ][0] - (sinkData['massAGB'][-1] - info['TotMaA'])
        info['RatioMaC_MLAGB'] = info['TotMaC']/ info['MassLostAGB']

        #Write text file with this info
        title = loc+str(run)+'_data_OrbitalEvolution.txt'
        with open (title,'w') as f:
            f.write('\n')
            f.write('Model '+str(run)+'\n')
            f.write('\n')
            f.write('The total mass accreted by the companion is                                    : '+ str(info['TotMaC'      ]) +' g \n')
            f.write('The total mass accreted by the AGB is                                          : '+ str(info['TotMaA'      ]) +' g \n')
            f.write('The total mass lost by the AGB is                                              : '+ str(info['MassLostAGB' ]) +' g \n')
            f.write('The ratio of the mass accreted by the companion to the mass lost by the AGB is : '+ str(round(info['RatioMaC_MLAGB'],5)))
            f.write('\n')

        #Visualise the orbit of the system
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        fig.set_size_inches(7, 7)
        plot_orbit(sinkData,setup, ax)
        ax.axis('equal')
        ax.legend(fontsize = 15, loc = 'center right')
        fig.tight_layout()
        plt.savefig(loc+str(run)+'_Orbit')


        # Calculate orbital separation and time in apastron and periastron passages
        info['OrbSep_a_p'] = orbSep_apa_per(sinkData, setup) 
        #Make plots of the change in sma and OrbSep 
        plotChangeOrbSep(info, sinkData, setup, run,loc)
        #Plot evolution of the mass accretion by the companion
        plotMassAccr(setup,sinkData, run, loc)
        #Plot evolution of orbital velocities 
        plotOrbVel(sinkData, run, loc)
        #Plot evolution of orbital radii and orbital separation
        plotOrbRad(sinkData, run, loc)


        #Write text file with usefull info
        title = loc+str(run)+'_data_OrbitalEvolution.txt'
        with open (title,'a') as f:
            f.write('\n')
            f.write('To plot mass accretion, orbital velocities and orbital radii: '+ '\n')
            f.write('\n')

            names = ['Time on x-axis [yrs]', 'Total accr mass comp [g]', 'Orbital Radius comp [cm]','Orbital Radius AGB [cm]', 'Orbital separation [cm]', 'Orbital vel comp [cm/s]', 'Orbital vel AGB [cm/s]']
            f.write("{: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34}".format(*names))
            f.write('\n')

            col_format = "{:<35}" * 7 + "\n"   # 7 left-justfied columns with 15 character width
            for i in zip(sinkData['time'], sinkData['maccrComp'],sinkData['rComp'],sinkData['rAGB'], sinkData['rComp']+sinkData['rAGB'], sinkData['v_orbComp_t'],sinkData['v_orbAGB_t'] ):
                    f.write(col_format.format(*i))

        print('     Orbital evolution plots of model '+ run +' ready and saved!')
