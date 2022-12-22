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
import Tools                    as tl

'''
Calculates:
    - times at which orbit is in apastron and periastron          (timeApa, timePer)
    - orbital separation at all apastron and periastron passages  (apastronOS, periastronOS)

returns first maximum, last maximum, first minimum, last minumimum
uses sinkData = data
'''
def orbSep_apa_per(data, setup):

    period = setup['period'] / cgs.year
    time   = data['time'   ]
    ecc    = setup['ecc'   ]
    rc     = data['rComp'  ]
    rAGB   = data['rAGB' ]
    
    
    if setup['triple_star']==True:
        M1   = data['massComp']
        M2   = data['massAGB' ]
        rCOM = (M1*rc + M2 * rAGB)/(M1+M2)
        orbSep = data['rComp_in'] + rCOM
    else:
        orbSep = rc + rAGB

    apastronOS   = []
    periastronOS = []  
    timeApa      = []
    timePer      = []
    i=0

    while i<math.floor(max(time)/period):
        # Model starts in apastron passage, so apastron passages are at i*period for i = 0, 1, ...
        # Not sure if this is true for triples.. :(
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
    #t  = data['time' ]
    #rC = data['rComp']
    #rA = data['rAGB' ]

    xc = data['posComp'].transpose()[0] 
    xa = data['posAGB' ].transpose()[0]
    yc = data['posComp'].transpose()[1]
    ya = data['posAGB' ].transpose()[1]
    #zc = data['posComp'].transpose()[2]
    #za = data['posAGB' ].transpose()[2]
    
    M1   = data['massComp']
    M2   = data['massAGB' ]
    
    ax.axis('equal')
    ax.set_ylim(-4,4)
    
    if setup['triple_star']==True:
        xc_in = data['posComp_in'].transpose()[0] 
        yc_in = data['posComp_in'].transpose()[1] 
        #zc_in = data['posComp_in'].transpose()[2] 
        M_in  = data['massComp_in']
        #CoM
        xCOM_in = (M_in*xc_in + M2 * xa)/(M_in+M2)
        yCOM_in = (M_in*yc_in + M2 * ya)/(M_in+M2)
        MCOM_in = M_in + M2
        xCOM_out = ( MCOM_in*xCOM_in + M1 *xc)/(MCOM_in+M1)
        yCOM_out = ( MCOM_in*yCOM_in + M1 *yc)/(MCOM_in+M1)
        #ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'k', label = 'CoM_inner')
        #ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'k', label = 'CoM_outer', markersize = 5)
        #ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'g', label = 'inner comp')
        #ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'r', label = 'outer comp')
        #FOR SCIENTIST@SCHOOL
        ax.set_facecolor('k')
        ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'white')#, label = 'MM')
        ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'white')#, label = 'MM', markersize = 5)
        ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'yellow')#, label = 'S2')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'yellow')#, label = 'S3')



    else:
        xCOM = (M1*xc + M2 * xa)/(M1+M2)
        yCOM = (M1*yc + M2 * ya)/(M1+M2)
        #rCOM = np.sqrt(xCOM**2 + yCOM**2)
        ax.plot(xCOM/cgs.au,yCOM/cgs.au,'*', color = 'k', label = 'CoM')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'r', label = 'companion')
    
    #ax.plot(xa[1:]/cgs.au,ya[1:]/cgs.au, c = 'b', label = 'AGB'      )
    #FOR SCIENTIST@SCHOOL
    ax.plot(xa[1:]/cgs.au,ya[1:]/cgs.au, c = 'r')#, label = 'AGB'      )

    #ax.set_ylabel('$y$ [cm]',fontsize = 15)
    #ax.set_xlabel('$x$ [cm]',fontsize = 15)
    ax.set_ylabel('$y$ [au]',fontsize = 15)
    ax.set_xlabel('$x$ [au]',fontsize = 15)
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
- Makes plots of the evolution of the semi-major axis and of 
    the orbital separation at the apastron and periastron passages.
- Calculates the change in eccentricity and semi-major axis.
- Returns text file with data to make the evolution plot and 
    with the calculated delta(e) and delta(a).
'''
def plotChangeOrbSep(info, sinkData, setup, run, loc):#, ylabel, unit, name, title):
    
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    
    # Legend
    apastron      = mlines.Line2D([],[], color = 'k', linestyle = 'None'  , marker = '$a$' ,label = 'Apastron', markersize = 8)
    periastron    = mlines.Line2D([],[], color = 'k', linestyle = 'None'  , marker = '$p$', label = 'Periastron', markersize =8)
    sma           = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'sma')
    handles1      = [apastron,sma, periastron]

    
    fig.set_size_inches(7, 5)
    c = 'k'
    
    # The data to plot are the orbital separations at apastron and periastron passages
    toPlot        = info['OrbSep_a_p']
    # Calculate the sma and mean time at every orbit
    [sma,orbit_t] = calcSMA(info['OrbSep_a_p'])
    
    # t_total is the final timestep
    t_total     = max(sinkData['time']) #in years

    # The change in sma per year is: delta a = (delta r_apa + delta r_per) / (2*time)
    delta_a              = (((toPlot[0]-toPlot[0][0] + (toPlot[1]-toPlot[1][0]) )/2)[-1]) / t_total
    ratio_delta_a_per_yr = (delta_a/(setup['sma_ini']*cgs.au))
    
    # The change in ecc per year is: delta e = (delta r_apa - delta r_per) / (2*(a+delta a) *time)
    delta_e              = (((toPlot[0]-toPlot[0][0] - (toPlot[1]-toPlot[1][0]) )/(2*setup['sma_ini']*cgs.au + 2*delta_a) )[-1])/t_total
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
    # Plot the change in sma at each orbit
    ax.plot(orbit_t, sma-sma[0], color = c, linestyle = 'dotted', markersize = 10)
    
    # Change xticks to orbit numbers
    ax.tick_params(labelsize=12)     
    period = setup['period'] / cgs.year
    tickorbit = []
    for p in range(1, int(t_total/period)+1,1):
        tickorbit.append(str(p))
        
    plt.setp(ax, xticks= orbit_t, xticklabels=tickorbit)
    
    ax.set_xlabel('Orbit', fontsize = 18)
    ax.set_ylabel('$\Delta$Orb sep [cm]', fontsize = 16)
    ax.set_title('Orbital evolution')
    ax.legend(handles = handles1, fontsize = 12)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'pdf/evolution_OrbitalSeparation.pdf'))
    plt.savefig(os.path.join(loc, 'png/evolution_OrbitalSeparation.png'))
    
    # Write text file with usefull info
    title = os.path.join(loc, 'txt/data_OrbitalEvolution.txt')
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
'''
def plotMassAccr(setup, sinkData, run, loc):
    # Make plot of the mass accretion evolution, very interesting to plot!
    fig = plt.figure(figsize=(8, 5))
    # Legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')

    # Plot the accreted mass in function of time
    plt.plot(sinkData['time'],  sinkData['maccrComp']/cgs.Msun, color = 'royalblue', linestyle = 'solid')
    if setup['triple_star']==True:
        plt.plot(sinkData['time'],  sinkData['maccrComp_in']/cgs.Msun, color = 'red', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'royalblue', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'r', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
    else:
        handles_ap    = [apaLine, perLine]

    # Plot vertical lines indicating where there are apastron and periastron passages
    period = setup['period'] / cgs.year
    #print('period in years: ',period)
    i = 0         # Start at apastron
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = max(np.max(sinkData['maccrComp']/cgs.Msun),np.max(sinkData['maccrComp_in']/cgs.Msun))
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
            i = i+period
        plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)   
    else:
        maxi = np.max(sinkData['maccrComp']/cgs.Msun)
        j = period/2  # First periastron
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
            plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
            i = i+period
            j = j+period
        plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
   
    ax = plt.subplot(111)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Accreted mass [Msun]', fontsize = 16)

    plt.title('Total accreted mass by the companion', fontsize = 18)
    plt.legend(handles = handles_ap)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_Maccr_companion.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_Maccr_companion.pdf'))



'''
Makes plot of the evolution of mass accretion rate by the companion
returns t_yrs and mass accretion rate to make plot yourself
'''
def plotMassAccrRate(setup, sinkData, run, loc):
    # Make plot of the mass accretion evolution, very interesting to plot!
    fig = plt.figure(figsize=(8, 5))
    # Legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    #print(min(sinkData['time']),max(sinkData['time']), len(sinkData['time']))
    
    # Make empty array for accrRates per year
    accrRates    = []
    if setup['triple_star']==True:
        accrRates_in = []
        # Make array with the years
        factor   = 0.5
    else:
        if setup["sma_ini"] > 20: #wide binary
            factor = 0.5
        else: #close binary
            factor   = 5.0
    yrs          = np.arange(1,int(max(sinkData['time']))+1,1./factor)   # 1 datapunt per 1/factor jaar
    # Make array with the indices of years in sinkData['time']
    indices_yrs  = [0]
    t_yrs        = [sinkData['time'][0]]
    
    #For each year, calculate the mass accreted in that 1 year and add it to the accRates array
    for year in yrs:
        i           = tl.find_nearest(sinkData['time'], year)
        indices_yrs.append(i)
        t_yrs.append(sinkData['time'][i])
        
        #calculate difference in accreted mass between this and the previous year, to find accretion rate per year in this year
        accrRate = factor*(sinkData['maccrComp'][indices_yrs[-1]]-sinkData['maccrComp'][indices_yrs[-2]])/cgs.Msun
        accrRates.append(accrRate)
        if setup['triple_star']==True:
            accrRate_in = factor*(sinkData['maccrComp_in'][indices_yrs[-1]]-sinkData['maccrComp_in'][indices_yrs[-2]])/cgs.Msun
            accrRates_in.append(accrRate_in)
    
    plt.plot(t_yrs[:-1], accrRates,color = 'royalblue', linestyle = 'solid')    
    if setup['triple_star']==True:
        plt.plot(t_yrs[:-1], accrRates_in,color = 'r', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'royalblue', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'r', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
    else:
        handles_ap    = [apaLine, perLine]



    # Plot vertical lines indicating where there are apastron and periastron passages
    period = setup['period'] * cgs.year
    #print('period in years: ',period)
    i = 0         # Start at apastron
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = max(np.max(accrRates),np.max(accrRates_in))
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
            i = i+period
        plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)   
    else:
        maxi = np.max(accrRates)
        j = period/2  # First periastron
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
            plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
            i = i+period
            j = j+period
        plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
   
    #ax = plt.subplot(111)
    plt.xlabel('Time[yrs]', fontsize = 14)
    plt.ylabel('Mass accretion rate [Msun/yr]', fontsize = 14)

    plt.title('Mass accretion rate by the companion', fontsize = 18)
    plt.legend(handles = handles_ap)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_MaccrRate_companion.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_MaccrRate_companion.pdf'))
    
    return(t_yrs[:-1], accrRates)


'''
Makes plot of the evolution of the orbital velocity of the AGB star and companion 
'''
def plotOrbVel(setup,sinkData, run, loc):
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    t    = sinkData['time'][1:]

    ax.plot(t, sinkData['v_orbAGB_t'][1:]/cgs.kms, label = 'AGB')
    ax.plot(t, sinkData['v_orbComp_t'][1:]/cgs.kms, label ='companion')
    

    if setup['triple_star']==True:
        ax.plot(t, sinkData['v_orbComp_in_t'][1:]/cgs.kms, label ='inner companion')

    period = setup['period'] / cgs.year
    i = 0         
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = max(np.max(sinkData['v_orbComp_t'][1:]/cgs.kms),np.max(sinkData['v_orbComp_in_t'][1:]/cgs.kms))
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
            i = i+period
        plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)   
        #handles_ap    = [periodLine,comp_in,comp_out]

    else:
        maxi = np.max(sinkData['v_orbComp_t'][1:]/cgs.kms)
        j = period/2  # First periastron
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
            plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
            i = i+period
            j = j+period
        plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
        #handles_ap    = [apaLine, perLine]
    
   
    ax.set_ylabel('$v_{orb}$ [km/s]', fontsize = 12)

    ax.set_xlabel('time[yrs]', fontsize = 10)

    ax.tick_params(labelsize=10)
    plt.title('Orbital velocities', fontsize = 15)
    plt.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_OrbitalVelocity.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_OrbitalVelocity.pdf'))




'''
Makes 1 plot of the evolution of the orbital radii of the AGB star and companion 
'''
def plotOrbRad(setup,sinkData, run, loc):
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})

    #ax.set_ylabel('$r$ [cm]', fontsize = 12)
    ax.set_ylabel('$r$ [au]', fontsize = 12)
    ax.set_title('orbital radii', fontsize = 15)
    ra   = sinkData['rAGB'][1:] /cgs.au
    t    = sinkData['time'][1:]
    
    
    ax.plot(t, ra, label ='r AGB')
    if setup['triple_star']==True:
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax.plot(t, rc_out, label= 'r outer comp')
        ax.plot(t, rc_in, label= 'r inner comp')
        
        Mc    = sinkData['massComp'][1:]
        Ma    = sinkData['massAGB' ][1:]
        Mc_in = sinkData['massComp_in'][1:]
        #CoM
        rCOM_in = (Mc_in*rc_in + Ma *ra )/(Mc_in+Ma)
        ax.plot(t, rCOM_in+rc_out, label = 'Orb sep outer')
 
        #ax.plot(sinkData['time'], sinkData['rAGB']+sinkData['rComp'], label = 'Orb sep')
    else:
        ax.plot(t, sinkData['rComp'][1:]/cgs.au, label= 'r comp')
        ax.plot(t, ra +sinkData['rComp'][1:]/cgs.au, label = 'Orb sep')
        
    # Plot vertical lines indicating where there are apastron and periastron passages
    # Legend
    #apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    #perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    #periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    period = setup['period'] / cgs.year
    #print('period in years: ',period)
    i = 0         # Start at apastron
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = np.max(rCOM_in+rc_out)
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
            i = i+period
        plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)   
        #handles_ap    = [periodLine,comp_in,comp_out]

    else:
        maxi = np.max(ra +sinkData['rComp'][1:]/cgs.au)
        j = period/2  # First periastron
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
            plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
            i = i+period
            j = j+period
        plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
        #handles_ap    = [apaLine, perLine]

        
    ax.set_xlabel('time [yrs]', fontsize = 14)
    ax.tick_params(labelsize=10)

    plt.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_rComp_rAGB_orbSep.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_rComp_rAGB_orbSep.pdf'))


'''
Makes seperate plots of the evolution of the orbital radii of the AGB star and companion(s)
'''
def plotOrbRadSeperate(setup,sinkData, run, loc):
    if setup['triple_star']==True:
        fig, axs= plt.subplots(2, 2,  gridspec_kw={'height_ratios':[1,1],'width_ratios': [1,1]},figsize = (20,10))
    else:
        fig, axs= plt.subplots(3, 1,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [1]})

    ra   = sinkData['rAGB'][1:] /cgs.au
    t    = sinkData['time'][1:]
    
    
    if setup['triple_star']==True:
        ax1 = axs[0][0]
        ax2 = axs[0][1]
        ax3 = axs[1][0]
        ax4 = axs[1][1]
        ax1.plot(t, ra, label ='r AGB')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax2.plot(t, rc_out, label= 'r outer comp')
        ax2.set_ylabel('r out comp [au]', fontsize = 12)
        ax3.plot(t, rc_in, label= 'r inner comp')
        ax3.set_ylabel('r inn comp [au]', fontsize = 12)        
        
        Mc    = sinkData['massComp'][1:]
        Ma    = sinkData['massAGB' ][1:]
        Mc_in = sinkData['massComp_in'][1:]
        #CoM
        rCOM_in = (Mc_in*rc_in + Ma *ra )/(Mc_in+Ma)

        ax4.plot(t, rCOM_in+rc_out, label = 'Orb sep outer')
        ax4.set_ylabel('orb sep outer [au]', fontsize = 12)        
 
        #ax.plot(sinkData['time'], sinkData['rAGB']+sinkData['rComp'], label = 'Orb sep')
    else:
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax1.plot(t, ra, label ='r AGB')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        ax2.set_ylabel('r comp [au]', fontsize = 12)
        ax2.plot(t, sinkData['rComp'][1:]/cgs.au, label= 'r comp')
        ax3.set_ylabel('Orb Sep [au]', fontsize = 12)
        ax3.plot(t, ra +sinkData['rComp'][1:]/cgs.au, label = 'Orb sep')
        
    # Plot vertical lines indicating where there are apastron and periastron passages
    # Legend
    #apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    #perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    #periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    period = setup['period'] / cgs.year
    #print('period in years: ',period)
    '''
    i = 0         # Start at apastron
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = np.max(rCOM_in+rc_out)
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            ax1.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
            i = i+period
        ax1.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)   
        #handles_ap    = [periodLine,comp_in,comp_out]

    else:
        maxi = np.max(ra +sinkData['rComp'][1:]/cgs.au)
        j = period/2  # First periastron
        for orbit in range(0, int(sinkData['time'][-1]/period)+1):
            ax1.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
            ax1.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
            i = i+period
            j = j+period
        ax1.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
        ax1.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
        #handles_ap    = [apaLine, perLine]

    '''    
    ax3.set_xlabel('time [yrs]', fontsize = 14)
    #ax.tick_params(labelsize=10)

    #plt.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/new_evolution_rComp_rAGB_orbSep.png'))
    plt.savefig(os.path.join(loc, 'pdf/new_evolution_rComp_rAGB_orbSep.pdf'))




'''
Main function executing calculations about orbital evolution
'''
def orbEv_main(run,loc, sinkData, setup):
    print('')
    if setup['single_star']:
        return
        
    else:
        # Calculate total mass accreted by companion and AGB star, total mass lost by the AGB star and ratio mass accreted to mass lost
        # Save it in dictionary 'info'
        info = {}
        info['TotMaC']         = sinkData['maccrComp'][-1] # the total mass accreted by the companion
        info['TotMaA']         = sinkData['maccrAGB' ][-1]  # the total mass accreted by the AGB
        info['MassLostAGB']    = sinkData['massAGB'  ][0] - (sinkData['massAGB'][-1] - info['TotMaA'])
        info['RatioMaC_MLAGB'] = info['TotMaC']/ info['MassLostAGB']
        if setup['triple_star']==True:
            info['TotMaC_in']         = sinkData['maccrComp_in'][-1] # the total mass accreted by the inner companion
            info['RatioMaC_in_MLAGB'] = info['TotMaC_in']/ info['MassLostAGB']

            
        
        # Write text file with this info
        title = os.path.join(loc, 'txt/data_OrbitalEvolution.txt')
        with open (title,'w') as f:
            f.write('\n')
            f.write('Model '+str(run)+'\n')
            f.write('\n')
            f.write('The total mass accreted by the companion is                                    : '+ str(info['TotMaC'      ]) +' g \n')
            f.write('The total mass accreted by the AGB is                                          : '+ str(info['TotMaA'      ]) +' g \n')
            f.write('The total mass lost by the AGB is                                              : '+ str(info['MassLostAGB' ]) +' g \n')
            f.write('The ratio of the mass accreted by the companion to the mass lost by the AGB is : '+ str(round(info['RatioMaC_MLAGB'],5)))
            f.write('\n')
            if setup['triple_star']==True:
                f.write('The total mass accreted by the inner companion is                                    : '+ str(info['TotMaC_in'      ]) +' g \n')
                f.write('The ratio of the mass accreted by the inner companion to the mass lost by the AGB is : '+ str(round(info['RatioMaC_in_MLAGB'],5)))


        # Visualise the orbit of the system
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        fig.set_size_inches(7, 7)
        plot_orbit(sinkData,setup, ax)
        ax.axis('equal')
        #ax.legend(fontsize = 15, loc = 'lower right')
        fig.tight_layout()

        plt.savefig(os.path.join(loc, 'png/orbit_sc@sc.png'))
        plt.savefig(os.path.join(loc, 'pdf/orbit_sc@sc.pdf'))

# werkt nog niet voor triples in elk geval
        '''
        if setup['triple_star']==False and setup['single_star']==False:
            # Calculate orbital separation and time in apastron and periastron passages
            info['OrbSep_a_p'] = orbSep_apa_per(sinkData, setup) 
            # Make plots of the change in sma and OrbSep 
            plotChangeOrbSep(info, sinkData, setup, run,loc)
        '''
        
        # Plot evolution of the mass accretion by the companion
        plotMassAccr(setup,sinkData, run, loc)
        
        # Plot evolution of the mass accretion rate by the companion
        (t_yrs, accrRates) = plotMassAccrRate(setup,sinkData, run, loc)
        #if setup['triple_star']==False:
            # Plot evolution of orbital velocities 
        plotOrbVel(setup,sinkData, run, loc)
        
        # Plot evolution of orbital radii and orbital separation
        #if setup['ecc'] == 0:
        plotOrbRadSeperate(setup,sinkData, run, loc)            
        #else:
        plotOrbRad(setup,sinkData, run, loc)


        # Write text file with usefull info
        title = os.path.join(loc, 'txt/data_OrbitalEvolution.txt')
        with open (title,'a') as f:
            f.write('\n')
            f.write('To plot mass accretion rate per year: '+ '\n')
            f.write('\n')
            names = ['Time in years on x-axis [yrs]', 'Mass accretion rates per year [Msun/yr]']
            f.write("{: <34} {: <34} ".format(*names))
            f.write('\n')
            col_format = "{:<35}" * 2 + "\n"   # 2 left-justfied columns with 35 character width
            for i in zip(t_yrs, accrRates):
                    f.write(col_format.format(*i))
                    
            f.write('\n')
            f.write('To plot mass accretion, orbital velocities and orbital radii: '+ '\n')
            f.write('\n')

            names = ['Time on x-axis [yrs]', 'Total accr mass comp [g]', 'Orbital Radius comp [cm]','Orbital Radius AGB [cm]', 'Orbital separation [cm]', 'Orbital vel comp [cm/s]', 'Orbital vel AGB [cm/s]' ]
            f.write("{: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34} ".format(*names))
            f.write('\n')

            col_format = "{:<35}" * 7 + "\n"   # 7 left-justfied columns with 35 character width
            for i in zip(sinkData['time'], sinkData['maccrComp'],sinkData['rComp'],sinkData['rAGB'], sinkData['rComp']+sinkData['rAGB'], sinkData['v_orbComp_t'][1:],sinkData['v_orbAGB_t'][1:]):
                    f.write(col_format.format(*i))
                    

        print('     Orbital evolution plots of model '+ run +' ready and saved!')
