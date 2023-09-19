# Import packages
import numpy                    as np
import matplotlib.pyplot        as plt
import matplotlib.lines         as mlines
import math                     as math
from scipy.signal import argrelextrema

import os
from matplotlib                 import rcParams
# Change the matplotlib default parameters
rcParams.update({'font.size' :   12})
rcParams.update({'figure.dpi': 200})
#rc('font', family='serif')
#rc('text', usetex=True)

# import plons scripts
import plons.ConversionFactors_cgs    as cgs
import plons.PhysicalQuantities       as pq
import plons.Tools                    as tl

'''
Visualises the orbit of the system
Uses sinkData = data
'''
def plot_orbit(data, setup,loc):

    # Visualise the orbit of the system
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(7, 4)

    xc = data['posComp'].transpose()[0]
    xa = data['posAGB' ].transpose()[0]
    yc = data['posComp'].transpose()[1]
    ya = data['posAGB' ].transpose()[1]
    #zc = data['posComp'].transpose()[2]
    #za = data['posAGB' ].transpose()[2]

    M1   = data['massComp']
    M2   = data['massAGB' ]

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
        ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'k', label = 'CoM_inner')
        ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'k', label = 'CoM_outer', markersize = 5)
        ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'gold', label = 'inner comp')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'crimson', label = 'outer comp')
        #  #FOR SCIENTIST@SCHOOL
        #ax.set_facecolor('k')
        #ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'white')#, label = 'MM')
        #ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'white')#, label = 'MM', markersize = 5)
        #ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'yellow')#, label = 'S2')
        #ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'yellow')#, label = 'S3')

    else:
        xCOM = (M1*xc + M2 * xa)/(M1+M2)
        yCOM = (M1*yc + M2 * ya)/(M1+M2)
        #rCOM = np.sqrt(xCOM**2 + yCOM**2)
        ax.plot(xCOM/cgs.au,yCOM/cgs.au, color = 'navy', label = 'CoM')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'crimson', label = 'companion')

    ax.plot(xa[1:]/cgs.au,ya[1:]/cgs.au, c = 'greenyellow', label = 'AGB'      )
    ax.set_ylabel('$y$ [au]',fontsize = 15)
    ax.set_xlabel('$x$ [au]',fontsize = 15)
    #ax.set_title('$e = $'+str(setup['ecc'])+', $a = $'+ str(setup['sma_ini'])+' AU, $q = $'+ str(setup['massAGB_ini']/setup['massComp_ini']), fontsize = 15)
    ax.set_title('(a)', y=0.9,x=0.95, fontweight='bold', ha='right', fontsize=16)
    ax.tick_params(labelsize=12)
    ax.axis('equal')
    ax.set_ylim(-70,100)
    ax.legend(fontsize = 15, loc = 'upper left')
    fig.tight_layout()
    plt.savefig(os.path.join(loc, 'png/orbit.png'))
    plt.savefig(os.path.join(loc, 'pdf/orbit.pdf'))



'''
Calculates:
    - times at which orbit is in apastron and periastron          (timeApa, timePer)
    - orbital separation at all apastron and periastron passages  (apastronOS, periastronOS)
returns times and values of relative maxima and minima of orbital separation
uses sinkData = data
'''
def orbSep_apa_per(data, setup):

    #data from sink files (so about evolution)
    time   = data['time'   ]
    rc     = data['rComp'  ]
    rAGB   = data['rAGB' ]

    #input value
    ecc    = setup['ecc'   ]

    if setup['triple_star']==True:
        M1   = data['massComp']
        M2   = data['massAGB' ]
        rCOM = (M1*rc + M2 * rAGB)/(M1+M2)
        #orbSep of outer orbit
        orbSep = data['rComp_in'] + rCOM
    else:
        orbSep = rc + rAGB

    #find relative maxima and minima
    #indicesApa, _ = find_peaks(orbSep,200*cgs.au)
    indicesApa   = argrelextrema(orbSep, np.greater, order=1500)  #order indicates how many neighbouring datapoints are considered when checking maxima/minima
    indicesPer   = argrelextrema(orbSep, np.less, order=1500)

    apastronOS   = np.array(orbSep[indicesApa])
    periastronOS = np.array(orbSep[indicesPer])
    timeApa      = np.array(time[indicesApa])
    timePer      = np.array(time[indicesPer])

    #Simulations stars in apastron, so its possible that there is one apastron separation value mores
    if len(apastronOS)>len(periastronOS):
        apastronOS   = apastronOS[:-1]
        timeApa      = timeApa[:-1]

    return apastronOS, periastronOS, timeApa, timePer


'''
Makes plot of the evolution of the orbital separation at apastron and periastron passage
'''
def plotApaPerChange(sinkData, setup, loc):
    fig, (ax1)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(7, 5)
    #Calculate orbital separation and time at apastron and periastron
    [apastronOS, periastronOS, timeApa, timePer] = orbSep_apa_per(sinkData, setup)

    ax1.set_xlabel('time [yrs]', fontsize = 16)
    # Plot the apastron orbital separations
    colorA = 'crimson'
    ax1.plot(timeApa,apastronOS/cgs.au, c=colorA, marker= '$*$', linestyle = 'dashed', markersize = 10)
    ax1.vlines(timeApa, 0.99*np.min(apastronOS)/cgs.au,1.015*np.max(apastronOS)/cgs.au, linestyle = 'dotted',color=colorA , linewidth = 0.5)
    ax1.tick_params(axis='y',labelsize=12, labelcolor = colorA)
    ax1.set_ylabel('Apastron orb sep [au]', fontsize = 16, color=colorA)
    #add second yaxis
    ax2 = ax1.twinx()
    # Plot the apastron orbital separations
    colorP = 'navy'
    ax2.plot(timePer,periastronOS/cgs.au, c=colorP, marker = '$*$', linestyle = 'dashed', markersize = 10)
    ax2.vlines(timePer, 0.99*np.min(periastronOS)/cgs.au,1.015*np.max(periastronOS)/cgs.au, linestyle = 'dotted',color=colorP , linewidth = 0.5)
    ax2.set_ylabel('Periastron orb sep [au]', fontsize = 16, color=colorP)
    ax2.tick_params(axis='y',labelsize=12, labelcolor = colorP)
    #ax1.set_title('Evolution apastron and periastron separation')
    ax1.set_title('(c)', y=0.9,x=0.95, fontweight='bold', ha='right', fontsize=16)
    fig.tight_layout()
    plt.savefig(os.path.join(loc, 'pdf/evolution_ApaPer.pdf'))
    plt.savefig(os.path.join(loc, 'png/evolution_ApaPer.png'))

'''
Makes plot of the evolution of the orbital period
And
Makes plot of the estimates of the evolution of the semi-major axis and eccentricity (are related, no exact values)
Uses
    # ra = a (1+e) = a + ae
    # rp = a (1-e) = a - ae
    # rp + ra = 2a
    # ra - rp = 2ae
'''
def plotEstimate_a_Per(sinkData,setup,loc):
    #Calculate orbital separation and time at apastron and periastron
    [apastronOS, periastronOS, timeApa, timePer] = orbSep_apa_per(sinkData, setup)
    # Estimates a and e, using the orbital separation at consequent apastron and periastron
    est_a = (apastronOS + periastronOS)/2
    est_e = (apastronOS - periastronOS)/(2*est_a)
    # Times corresponding best to estimates are times inbetween the apastron and periastron used
    t = (timeApa+timePer)/2

    # Find indices corresponding to these times
    ind  = []
    for time in t:
        i  = tl.find_nearest(sinkData['time'], time)
        ind.append(i)

    # Estimate period using the estimate of a [au], and the masses of the AGB and companion [gram] at these times t
    est_p     = pq.getPeriod(sinkData['massAGB'][ind], sinkData['massComp'][ind], est_a/cgs.au )/cgs.year

    # Make plot of orbital period, and add linear line for visualisation
    fig_p, ((ax_p))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    ax_p.plot(t, est_p, label = 'orbital period',marker= '$*$', linestyle = 'solid', color = 'navy',markersize = 10)
    ax_p.plot([t[0],t[-1]],[est_p[0],est_p[-1]],'--',label='linear',color ='navy',linestyle='dashed',linewidth = 0.8)
    ax_p.set_ylabel('estimate orbital period [yrs]', fontsize = 12)
    ax_p.set_xlabel('time[yrs]', fontsize = 10)
    ax_p.tick_params(labelsize=10)
    ax_p.grid(which='both')
    ax_p.set_ylim(845,1240)
    ax_p.legend(fontsize = 15, loc = 'upper left')
    #plt.title('Orbital Period', fontsize = 15)
    plt.title('(d)', y=0.9,x=0.95, fontweight='bold', ha='right', fontsize=16)
    fig_p.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_OrbitalPeriod.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_OrbitalPeriod.pdf'))
    plt.close

    #Make plot of estimates of a and e
    fig, (ax1)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(7, 5)
    ax1.set_xlabel('time [yrs]', fontsize = 16)
    # Plot the estimates of a
    color_a = 'crimson'
    ax1.plot(t,est_a/cgs.au, c=color_a, marker= '$*$', linestyle = 'dashed', markersize = 10)
    ax1.vlines(t, 0.99*np.min(est_a)/cgs.au,1.01*np.max(est_a)/cgs.au, linestyle = 'dotted',color=color_a , linewidth = 0.5)
    ax1.tick_params(axis='y',labelsize=12, labelcolor = color_a)
    ax1.set_ylabel('estimate semi-major axis [au]', fontsize = 16, color=color_a)
    # Plot the estimates of e, with seperate y-axis
    ax2 = ax1.twinx()
    color_e = 'navy'
    ax2.plot(t,est_e, c=color_e, marker = '$*$', linestyle = 'dashed', markersize = 10)
    ax2.vlines(t, 0.99*np.min(est_e),1.01*np.max(est_e), linestyle = 'dotted',color=color_e , linewidth = 0.5)
    ax2.set_ylabel('estimate eccentricity', fontsize = 16, color=color_e)
    ax2.tick_params(axis='y',labelsize=12, labelcolor = color_e)
    ax1.set_title('Evolution semi-major axis and eccentricity')
    fig.tight_layout()
    plt.savefig(os.path.join(loc, 'pdf/evolution_a_e.pdf'))
    plt.savefig(os.path.join(loc, 'png/evolution_a_e.png'))


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

    # Plot the accreted mass as a function of time
    plt.plot(sinkData['time'],  sinkData['maccrComp']/cgs.Msun, color = 'crimson', linestyle = 'solid')
    if setup['triple_star']==True:
        plt.plot(sinkData['time'],  sinkData['maccrComp_in']/cgs.Msun, color = 'gold', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'gold', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
    #else:
        handles_ap    = [apaLine, perLine]
    '''
    # Plot vertical lines indicating where there should be apastron and periastron passages
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
    '''
    ax = plt.subplot(111)
    ax.grid(which='both')
    ax.set_ylim(-0.000003,0.000195)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Accreted mass [Msun]', fontsize = 16)
    #plt.title('Total accreted mass by the companion', fontsize = 18)
    plt.title('(a)', y=0.9,x=0.95, fontweight='bold', ha='right', fontsize=16)
    if setup['triple_star']==True:
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

    plt.plot(t_yrs[:-1], accrRates,color = 'crimson', linestyle = 'solid')
    if setup['triple_star']==True:
        plt.plot(t_yrs[:-1], accrRates_in,color = 'gold', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'gold', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
    #else:
        #handles_ap    = [apaLine, perLine]

    '''
    # Plot vertical lines indicating where there are apastron and periastron passages
    period = setup['period'] / cgs.year
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
    '''
    ax = plt.subplot(111)
    ax.grid()
    plt.xlabel('Time[yrs]', fontsize = 14)
    plt.ylabel('Mass accretion rate [Msun/yr]', fontsize = 14)

    plt.title('Mass accretion rate by the companion', fontsize = 18)
    if setup['triple_star']==True:
        plt.legend(handles = handles_ap)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_MaccrRate_companion.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_MaccrRate_companion.pdf'))

    #return(t_yrs[:-1], accrRates)


'''
Makes plot of the evolution of the orbital velocity of the AGB star and companion
'''
def plotOrbVel(setup,sinkData, run, loc):
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    t    = sinkData['time'][1:]

    ax.plot(t, sinkData['v_orbAGB_t'][1:]/cgs.kms, label = 'AGB', c='greenyellow')
    ax.plot(t, sinkData['v_orbComp_t'][1:]/cgs.kms, label ='companion', c='crimson')

    if setup['triple_star']==True:
        ax.plot(t, sinkData['v_orbComp_in_t'][1:]/cgs.kms, label ='inner companion', c='gold')

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

    ax.set_ylabel('$r$ [au]', fontsize = 12)
    #ax.set_title('orbital radii', fontsize = 15)
    ax.set_title('(b)', y=0.9,x=0.95, fontweight='bold', ha='right', fontsize=16)
    ra   = sinkData['rAGB'][1:] /cgs.au
    t    = sinkData['time'][1:]

    ax.plot(t, ra, label ='r AGB', c='greenyellow')
    if setup['triple_star']==True:
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax.plot(t, rc_out, label= 'r outer comp',c='crimson')
        ax.plot(t, rc_in, label= 'r inner comp',c='gold')

        Mc    = sinkData['massComp'][1:]
        Ma    = sinkData['massAGB' ][1:]
        Mc_in = sinkData['massComp_in'][1:]
        #CoM
        rCOM_in = (Mc_in*rc_in + Ma *ra )/(Mc_in+Ma)
        ax.plot(t, rCOM_in+rc_out, label = 'Orb sep outer',c='navy')

        #ax.plot(sinkData['time'], sinkData['rAGB']+sinkData['rComp'], label = 'Orb sep')
    else:
        ax.plot(t, sinkData['rComp'][1:]/cgs.au, label= 'r comp', c='crimson')
        ax.plot(t, ra +sinkData['rComp'][1:]/cgs.au, label = 'Orb sep', c='navy')
    '''
    # Plot vertical lines indicating where there should be apastron and periastron passages
    # Legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
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
    '''
    ax.grid()
    ax.set_xlabel('time [yrs]', fontsize = 14)
    ax.tick_params(labelsize=10)
    ax.legend(fontsize = 15, loc = 'upper left')
    ax.set_ylim(-5,330)
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
        ax1.plot(t, ra, label ='r AGB',c='greenyellow')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax2.plot(t, rc_out, label= 'r outer comp',c='crimson')
        ax2.set_ylabel('r out comp [au]', fontsize = 12)
        ax3.plot(t, rc_in, label= 'r inner comp',c='gold')
        ax3.set_ylabel('r inn comp [au]', fontsize = 12)

        Mc    = sinkData['massComp'][1:]
        Ma    = sinkData['massAGB' ][1:]
        Mc_in = sinkData['massComp_in'][1:]
        #CoM
        rCOM_in = (Mc_in*rc_in + Ma *ra )/(Mc_in+Ma)

        ax4.plot(t, rCOM_in+rc_out, label = 'Orb sep outer',c='navy')
        ax4.set_ylabel('orb sep outer [au]', fontsize = 12)

    else:
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
        ax1.plot(t, ra, label ='r AGB',c='greenyellow')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        ax2.set_ylabel('r comp [au]', fontsize = 12)
        ax2.plot(t, sinkData['rComp'][1:]/cgs.au, label= 'r comp',c='crimson')
        ax3.set_ylabel('Orb Sep [au]', fontsize = 12)
        ax3.plot(t, ra +sinkData['rComp'][1:]/cgs.au, label = 'Orb sep',c='navy')

    '''
    # Plot vertical lines indicating where there should be apastron and periastron passages
    # Legend
    #apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    #perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    #periodLine    = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    period = setup['period'] / cgs.year
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
    ax1.grid()
    ax2.grid()
    ax3.grid()
    #plt.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_rComp_rAGB_orbSep_seperate.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_rComp_rAGB_orbSep_seperate.pdf'))


'''
Main function executing calculations about orbital evolution
'''
def orbEv_main(run,loc, sinkData, setup):
    print('')
    if setup['single_star']:
        return

    else:
        #print('(5)  Start calculations for orbital evolution...')  # already printed in main

        # Visualise the orbit of the system
        plot_orbit(sinkData,setup,loc)

        #For eccentric binaries
        if setup['ecc']>0 and setup['triple_star']==False and setup['single_star']==False:
            # Make plot of estimates of orbital period, a and e, using apastron and periastron passages
            plotEstimate_a_Per(sinkData,setup,loc)
            # Make plot of change in orbital separation at apastron and periastron passages
            plotApaPerChange(sinkData,setup,loc)

        # Plot evolution of the mass accretion by the companion
        plotMassAccr(setup,sinkData, run, loc)

        # Plot evolution of the mass accretion rate by the companion
        plotMassAccrRate(setup,sinkData, run, loc)
        # Plot evolution of orbital velocities
        plotOrbVel(setup,sinkData, run, loc)

        # Plot evolution of orbital radii and orbital separation
        #if setup['ecc'] == 0:
        plotOrbRadSeperate(setup,sinkData, run, loc)
        #else:
        plotOrbRad(setup,sinkData, run, loc)

        '''
        #   UNCOMMENT (PARTS) IF YOU WANT TO WRITE OUT INFO IN TEXT FILE
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

        # Write text file with usefull info; if you want to do this, first save t_yrs and accrRates from plotMassAccrRate function
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
        '''

        print('     Orbital evolution plots of model '+ run +' ready and saved!')
