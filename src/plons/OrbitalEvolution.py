# Import packages
from cProfile import label, run
from re import T
import numpy                    as np
import matplotlib.pyplot        as plt
import matplotlib.lines         as mlines
import math                     as math
from scipy.signal               import argrelextrema

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
    fig.set_size_inches(7, 7)

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
        ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'navy', label = 'CoM_inner')
        ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'navy', label = 'CoM_outer', markersize = 5)
        ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'lime', label = 'inner comp')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'crimson', label = 'outer comp')
        # for figure intro PhD thesis
        # ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'green', label = 'CoM_*', markersize = 5)
        # ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'darkred',linestyle=(0,(5,10)), label = 'B')
        # ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'green', label = 'CoM_A')
        # ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'darkred',linestyle='dotted', label = 'Ab')
        #  #FOR SCIENTIST@SCHOOL
        #ax.set_facecolor('navy')
        #ax.plot(xCOM_in/cgs.au,yCOM_in/cgs.au, color = 'white')#, label = 'MM')
        #ax.plot(xCOM_out/cgs.au,yCOM_out/cgs.au,'*', color = 'white')#, label = 'MM', markersize = 5)
        #ax.plot(xc_in[1:]/cgs.au,yc_in[1:]/cgs.au, c = 'yellow')#, label = 'S2')
        #ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'yellow')#, label = 'S3')

    else:
        xCOM = (M1*xc + M2 * xa)/(M1+M2)
        yCOM = (M1*yc + M2 * ya)/(M1+M2)
        #rCOM = np.sqrt(xCOM**2 + yCOM**2)
        ax.plot(xCOM/cgs.au,yCOM/cgs.au,'*', color = 'navy', label = 'CoM')
        ax.plot(xc[1:]/cgs.au,yc[1:]/cgs.au, c = 'crimson', label = 'companion')

    ax.plot(xa[1:]/cgs.au,ya[1:]/cgs.au, c = 'gold', label = 'AGB'      )
    # ax.plot(xa[1:]/cgs.au,ya[1:]/cgs.au, c = 'darkred',linestyle='-', label = 'Aa'      )
    ax.set_ylabel('$y$ [au]',fontsize = 15)
    ax.set_xlabel('$x$ [au]',fontsize = 15)
    ax.set_title('$e = $'+str(setup['ecc'])+', $a = $'+ str(setup['sma_ini'])+' AU, $q = $'+ str(setup['massAGB_ini']/setup['massComp_ini']), fontsize = 15)
    ax.tick_params(labelsize=12)
    ax.axis('equal')
    ax.legend(fontsize = 15, loc = 'lower right')
    fig.tight_layout()
    plt.savefig(os.path.join(loc, 'png/orbit.png'))
    plt.savefig(os.path.join(loc, 'pdf/orbit.pdf'))


def plotVerticalApaPerLines(setup,sinkData,mini,maxi):
    # Plot vertical lines indicating where there should be apastron and periastron passages
        period = setup['period'] / cgs.year
        i = 0         # Simulation starts at apastron
        mini = 0
        if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
            for orbit in range(0, int(sinkData['time'][-1]/period)+1):
                plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
                i = i+period
            plt.vlines(i,mini, maxi,  linestyle = 'dotted' , linewidth = 0.5)
        else:
            j = period/2  # First periastron
            for orbit in range(0, int(sinkData['time'][-1]/period)):
                plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
                plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
                i = i+period
                j = j+period
            plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)  #end with last apastron
            # plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)



'''
Calculates the instantaneous eccentricity using the Laplace-Runge-Lenz vector:
returns an array containing the evolution of the eccentricity
uses sinkData = data
'''
def plotEccentricity(data, setup, loc):

    #data from sink files (so about evolution)
    r1     = data['posAGB' ]
    v1     = data['velAGB' ]
    M1     = data['massAGB' ]
    r2     = data['posComp' ]
    v2     = data['velComp' ]
    M2     = data['massComp']
    t      = data['time'][1:]
    
    #centre of mass position and velocity
    Mtot  = (M1+M2)*cgs.G

    #position and velocity coordinates relative to each other
    r  = r2 - r1
    v  = v2 - v1
    rnorm = np.sqrt(np.sum(r**2,axis=1))
    v_sqr = np.sum(v**2,axis=1)
    rv    = np.sum(r*v,axis=1)

    a = v_sqr/Mtot-1./rnorm
    ecc = r*a[:, None] - v*rv[:, None]/Mtot[:, None]

    eccentricity = np.sqrt(np.sum(ecc**2,axis=1))

    # trueAnomaly = np.arccos(np.sum(ecc*r,axis=1)/(rnorm*eccentricity))*180/(math.pi)

    print('The eccentricity is about ',np.mean(eccentricity[1:]))
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    maxi = np.max(eccentricity[1:])
    mini = np.min(eccentricity[1:])
    plotVerticalApaPerLines(setup,data,mini,maxi)
    if setup['ecc']>0:
        plt.hlines(setup['ecc'],t[0],t[-1])
    ax.set_xlabel('time [yrs]', fontsize = 16)
    #ax.plot(trueAnomaly[1:], eccentricity[1:], c='gold')
    ax.plot(t, eccentricity[1:], c='gold')
    #ax.plot(t, Mtot[1:]/cgs.Msun, c='gold')
    ax.set_title('Evolution of eccentricity')
    fig.tight_layout()
    plt.savefig(os.path.join(loc, 'pdf/evolution_ecc.pdf'))
    plt.savefig(os.path.join(loc, 'png/evolution_ecc.png'))


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
    # indicesApa   = argrelextrema(orbSep, np.greater, order=1500)  #order indicates how many neighbouring datapoints are considered when checking maxima/minima

    indicesApa   = argrelextrema(orbSep, np.greater, order=5000)  #order indicates how many neighbouring datapoints are considered when checking maxima/minima
    indicesPer   = argrelextrema(orbSep, np.less, order=5000)

    apastronOS   = np.array(orbSep[indicesApa])
    periastronOS = np.array(orbSep[indicesPer])
    timeApa      = np.array(time[indicesApa])
    timePer      = np.array(time[indicesPer])


    #Simulations starts in apastron, so its possible that there is one apastron separation value more
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
    ax1.vlines(timeApa, 0.99*np.min(apastronOS)/cgs.au,1.01*np.max(apastronOS)/cgs.au, linestyle = 'dotted',color=colorA , linewidth = 0.5)
    ax1.tick_params(axis='y',labelsize=12, labelcolor = colorA)
    ax1.set_ylabel('Apastron orb sep [au]', fontsize = 16, color=colorA)
    #add second yaxis
    ax2 = ax1.twinx()
    # Plot the apastron orbital separations
    colorP = 'navy'
    ax2.plot(timePer,periastronOS/cgs.au, c=colorP, marker = '$*$', linestyle = 'dashed', markersize = 10)
    ax2.vlines(timePer, 0.99*np.min(periastronOS)/cgs.au,1.01*np.max(periastronOS)/cgs.au, linestyle = 'dotted',color=colorP , linewidth = 0.5)
    ax2.set_ylabel('Periastron orb sep [au]', fontsize = 16, color=colorP)
    ax2.tick_params(axis='y',labelsize=12, labelcolor = colorP)
    ax1.set_title('Evolution apastron and periastron separation')
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
    plt.title('Orbital Period', fontsize = 15)
    plt.legend()
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
    ax1.vlines(t, np.min(est_a)/cgs.au,np.max(est_a)/cgs.au, linestyle = 'dotted',color=color_a , linewidth = 0.5)
    ax1.tick_params(axis='y',labelsize=12, labelcolor = color_a)
    ax1.set_ylabel('estimate semi-major axis [au]', fontsize = 16, color=color_a)
    # Plot the estimates of e, with seperate y-axis
    ax2 = ax1.twinx()
    color_e = 'navy'
    ax2.plot(t,est_e, c=color_e, marker = '$*$', linestyle = 'dashed', markersize = 10)
    # ax2.vlines(t, 0.99*np.min(est_e),1.01*np.max(est_e), linestyle = 'dotted',color=color_e , linewidth = 0.5)
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
    apaLine       = mlines.Line2D([],[], color = 'navy', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'navy', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'navy', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')

    # Plot the accreted mass as a function of time
    plt.plot(sinkData['time'],  sinkData['maccrComp']/cgs.Msun, color = 'crimson', linestyle = 'solid')
    if setup['triple_star']==True:
        plt.plot(sinkData['time'],  sinkData['maccrComp_in']/cgs.Msun, color = 'lime', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'lime', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
    #else:
        handles_ap    = [apaLine, perLine]

    ax = plt.subplot(111)
    ax.grid(which='both')
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Accreted mass [Msun]', fontsize = 16)

    plt.title('Total accreted mass by the companion', fontsize = 18)
    if setup['triple_star']==True:
        plt.legend(handles = handles_ap,loc='upper right')
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_Maccr_companion.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_Maccr_companion.pdf'))



def reduceToNPeriods(setup,time):
    finalTime = time[-1]
    period = setup['period'] / cgs.year
    fullOrbits = int(finalTime/period)
    tmax = fullOrbits * setup['period'] / cgs.year
    ind_tmax = np.where(time > tmax)[0]
    index_tmax = ind_tmax[0]+1
    time_new = time[:index_tmax]
    # print('finalTimestep',time_new[-1])
    return time_new,index_tmax


'''
Calculate theoretical BHL mass accretion rate
'''
def BHLMassAccrRate(setup,data):
    # If you want to reduce calculations to N full periods
    # if setup['ecc']>0:
    time, index_tmax = reduceToNPeriods(setup,data['time'])

    # else
    # time = data['time'][:-1]
    # index_tmax = -1

    alphaBHL = 0.75
    Mdot     = setup['Mdot'] * cgs.Msun / cgs.year                               # [g/s]
    orbSep   = (data['rComp'][:index_tmax] + data['rAGB'][:index_tmax])                                    # cm   ?
    G        = cgs.G                                                             # [cm^3 g^-1 s^-1]
    Mc       = data['massComp'][:index_tmax]                                                  # [g]
    Md       = data['massAGB'][:index_tmax]
    vorbp    = data['v_orbAGB_t'][:index_tmax]                                                # [cm/s]  ?
    # vorbc    = data['v_orbComp_t']#/cgs.kms                                     # [cm/s]

    if setup['v_ini'] == 5:
        vw = 12.7*cgs.kms
    elif setup['v_ini'] == 10:
        vw = 15.0*cgs.kms
    elif setup['v_ini'] == 20:
        vw = 22.8*cgs.kms
    else:  #Incorrect!
        vw       = np.sqrt(np.mean(vorbp)**2 + (setup['v_ini']*cgs.kms)**2)          # [cm/s]  local wind velocity, wordt berekend in ander script, nu afschatten door np.sqrt(vw^2+vorb^2) En gemiddelde, want klopt niet anders

    print('vw', vw/cgs.kms)
    
    # vtest2       = vorbc + vorbp                                              # Instantaneous relative velocity of stars, exactly the same as v
    #INCORRECT FORMULA...? v        = np.sqrt(G * (Mc + Md)/orbSep)                                     # [cm/s] Instantaneous relative velocity of stars
    # vtest        = np.sqrt(G * (Mc + Md)/orbSep)
    sma      = setup['sma_ini']        *cgs.au                               # [cm]  

    v       = np.sqrt(G * (Mc + Md)*(2/orbSep - 1/sma))                         # !Correct relative orbital velocity in eccentric systems
    print('mean rel orb vel: ',np.mean(v)/cgs.kms)
    
    MaccrBHL = (alphaBHL * (Mdot / orbSep**2)  * ( G * Mc / ( vw**2))**2 * (1 +  (v/vw)**2)**(-1.5) )    #g/s
    MaccrEffBHL = MaccrBHL/Mdot #*100  #%
    ecc      = setup['ecc']
    sma      = setup['sma_ini']        *cgs.au                               # [cm]  
    Mc       = setup['massComp_ini']   *cgs.Msun                             # [g]
    Md       = setup['massAGB_ini']    *cgs.Msun                             # [g]
    v_orb    = np.sqrt(G * (Mc + Md)/sma)                                    # [cm/s]
    print('mean vorb from input parameters: ',v_orb/cgs.kms)    
    print('-------------------')

   # average mass accretion rate formula:
    # MaccrBHL_av    = (alphaBHL * (Mdot/  (sma**2 * np.sqrt(1 - ecc**2))) * ( G * Mc / ( vw**2))**2 * (1 + (v_orb/vw)**2)**(-1.5) )  #g/s
    # average mass accretion efficiency:
    # MaccrEffBHL_av  = MaccrBHL_av / Mdot #*100 # %
    # print('Average BHL accr efficiency (formula)',MaccrEffBHL_av)
    # print('time-averaged mean BHL MaccrEff: ',np.mean(MaccrEffBHL))
    # print('---------------------')

    return time,MaccrBHL #, MaccrEffBHL, MaccrBHL_av,MaccrEffBHL_av


def calculateAverageAccrRate_lastNorbits(timeArray,accrRates,period,n):
    firstTimestep = timeArray[-1] - n * period
    index_firsttimestep = np.where(timeArray > firstTimestep)[0][0]
    meanAccrRate = np.mean(accrRates[index_firsttimestep:])
    return meanAccrRate

def calculateAccretedMassRate_lastNorbits(timeArray,maccrComp,period,n):
    firstTimestep = timeArray[-1] - n * period
    index_firsttimestep = np.where(timeArray > firstTimestep)[0][0]
    totalMassAccr =(maccrComp[-1]-maccrComp[index_firsttimestep])/cgs.Msun
    accretedMassRate = totalMassAccr/(n*period)
    return accretedMassRate

def calc_accrRates(setup,sinkData):
    # If you want to reduce to N periods:
    time, index_tmax = reduceToNPeriods(setup,sinkData['time'])
    # Else
    # time = sinkData['time']
    # index_tmax = -1

    # Make empty array to store accretion rates
    accrRates    = np.array([]) 
    if setup['triple_star']==True:
        accrRates_in = np.array([]) 
        factor   = 0.5
    else:
        if setup["sma_ini"] > 20: #wide binary
            factor = 0.5
        else: #close binary
            factor   = 5.0   #1 datapoint every 1/5 years --> 5 datapoints per year

    # Make array of timesteps on which you calculate mass accretion rate
    timeArray    = np.arange(0,int(time[-1]),1./factor)   # 1 datapoint per year/factor
    # Make array with the indices of selected times in original dataset sinkData['time']
    indices_time  = np.array([0])
    accretedMass_t_previous = 0
    #For each time in timeArray, calculate the mass accreted wrt previous time, and add it to the accRates array
    for t in timeArray:
        i            = tl.find_nearest(sinkData['time'], t)  #find the time datapoint in the original data array
        indices_time = np.append(indices_time,i) 
        accretedMass_t = sinkData['maccrComp'][indices_time[-1]]
        #calculate difference in accreted mass between this and the previous timestep, to find accretion rate
        accrRate  = factor*(accretedMass_t-accretedMass_t_previous)/cgs.Msun   #accreted mass * factor to get accretion rate per year
        accrRates = np.append(accrRates,accrRate)
        accretedMass_t_previous = accretedMass_t

        if setup['triple_star']==True:
            accrRate_in = factor*(sinkData['maccrComp_in'][indices_time[-1]]-sinkData['maccrComp_in'][indices_time[-2]])/cgs.Msun
            accrRates_in = np.append(accrRates_in,accrRate_in)

    period = setup['period'] / cgs.year
    massAccretedEff = calculateAccretedMassRate_lastNorbits(time,sinkData['maccrComp'][:index_tmax],period,2)

    # print('Mass accretion rate over last 2 orbits is ', massAccretedEff)
    Mdot     = setup['Mdot']                                # [Msun/yr]
    print('Mass accreted / last 2 orbital periods /Mdot is ', massAccretedEff/Mdot)
    print('-------------------')


    if setup['triple_star']==True:
        return timeArray,accrRates,accrRates_in
    else:
        return timeArray,accrRates



'''
Makes plot of the evolution of mass accretion rate by the companion
returns timeArray and mass accretion rate to make plot yourself
'''
def plotMassAccrRate(setup, sinkData, loc):
    # Make plot of the mass accretion evolution, very interesting to plot!
    fig = plt.figure(figsize=(8, 5))
    # Legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'navy', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    BHL           = mlines.Line2D([],[], color = 'royalblue', linestyle = 'dotted', linewidth = 0.8, label = 'BHL')
    simulation    = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 1, label = 'simulation')

    if setup['triple_star']==True:
        timeArray,accrRates,accrRates_in = calc_accrRates(setup,sinkData)
        plt.plot(timeArray, accrRates,color = 'crimson', linestyle = 'solid')    
        plt.plot(timeArray[:-1], accrRates_in[:-1],color = 'lime', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'lime', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
        maxi = max(np.max(accrRates),np.max(accrRates_in))

    else:
        timeArray,accrRates = calc_accrRates(setup,sinkData)
        plt.plot(timeArray, accrRates,color = 'crimson', linestyle = 'solid')    
        handles_ap    = [BHL,simulation, apaLine, perLine]
        #Plot BHL mass accretion rate
        time,MaccrBHL = BHLMassAccrRate(setup,sinkData)
        plt.plot(time, MaccrBHL/cgs.Msun * cgs.year,color = 'royalblue', linestyle = 'dotted',linewidth=0.8)    
        maxi = 1.1* max(np.max(accrRates),np.max(MaccrBHL/cgs.Msun * cgs.year))
        
    mini = 0
    plotVerticalApaPerLines(setup,sinkData,mini,maxi)
    
    # ax = plt.subplot(111)
    # ax.grid()
    plt.xlabel('Time [yrs]', fontsize = 14)
    plt.ylabel('Mass accretion rate [Msun/yr]', fontsize = 14)

    # plt.title('Mass accretion rate by the companion', fontsize = 18)
    plt.legend(handles = handles_ap)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_MaccrRate_companion_BHLcomp.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_MaccrRate_companion_BHLcomp.pdf'))
    plt.close()


def plotMassAccrEff(setup, sinkData, loc):
    ### Mass accretion efficiency plot
    Mdot     = setup['Mdot']                                # [Msun/yr]
    
    # Make plot of the mass accretion evolution, very interesting to plot!
    fig = plt.figure(figsize=(8, 5))
    # Legend
    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')
    periodLine    = mlines.Line2D([],[], color = 'navy', linestyle = 'dotted', linewidth = 0.5, label = 'orbital period')
    BHL           = mlines.Line2D([],[], color = 'royalblue', linestyle = 'dotted', linewidth = 0.8, label = 'BHL')
    simulation    = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 1, label = 'simulation')

    if setup['triple_star']==True:
        timeArray,accrRates,accrRates_in = calc_accrRates(setup,sinkData)
        plt.plot(timeArray, accrRates/Mdot,color = 'crimson', linestyle = 'solid')    
        plt.plot(timeArray[:-1], accrRates_in[:-1]/Mdot,color = 'lime', linestyle = 'solid')
        comp_out       = mlines.Line2D([],[], color = 'crimson', linestyle = 'solid', linewidth = 0.5, label = 'outer comp')
        comp_in        = mlines.Line2D([],[], color = 'lime', linestyle = 'solid', linewidth = 0.5, label = 'inner comp')
        handles_ap    = [periodLine,comp_in,comp_out]
        maxi = max(np.max(accrRates),np.max(accrRates_in))/Mdot

        period = setup['period'] / cgs.year
        meanAccrEff = calculateAverageAccrRate_lastNorbits(timeArray,accrRates,period,2)/Mdot
        print('mean accretion eff of outer companion over last 2 orbits is ', meanAccrEff)
        meanAccrEff = calculateAverageAccrRate_lastNorbits(timeArray,accrRates_in,period,2)/Mdot
        print('mean accretion eff of inner companion over last 2 orbits is ', meanAccrEff)
    else:
        timeArray,accrRates = calc_accrRates(setup,sinkData)
        plt.plot(timeArray, accrRates/Mdot,color = 'crimson', linestyle = 'solid')    
        handles_ap    = [BHL,simulation, apaLine, perLine]
        #Plot BHL mass accretion rate
        time,MaccrBHL = BHLMassAccrRate(setup,sinkData)
        # plt.plot(sinkData['time'], MaccrBHL/cgs.Msun * cgs.year,color = 'royalblue', linestyle = 'dotted',linewidth=0.8)    
        plt.plot(time, MaccrBHL/cgs.Msun * cgs.year /Mdot,color = 'royalblue', linestyle = 'dotted',linewidth=0.8)    
        maxi = 1.1* max(np.max(accrRates),np.max(MaccrBHL/cgs.Msun * cgs.year))/Mdot
    
        period = setup['period'] / cgs.year
        meanAccrEff = calculateAverageAccrRate_lastNorbits(timeArray,accrRates,period,2)/Mdot
        meanAccrEffBHL = calculateAverageAccrRate_lastNorbits(time,MaccrBHL/cgs.Msun * cgs.year,period,2)/Mdot
        print('mean accretion eff over last 2 orbits is ', meanAccrEff)
        print('mean BHL accretion eff over last 2 orbits is ', meanAccrEffBHL)


    mini = 0
    plotVerticalApaPerLines(setup,sinkData,mini,maxi)
    
    # ax = plt.subplot(111)
    # ax.grid()
    plt.xlabel('Time [yrs]', fontsize = 14)
    plt.ylabel(r'Mass accretion efficiency $\beta$', fontsize = 14)

    # plt.title('Mass accretion efficiency', fontsize = 18)
    plt.legend(handles = handles_ap)
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_MaccrEff_companion_BHLcomp.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_MaccrEff_companion_BHLcomp.pdf'))
    plt.close()


def plot_vw_vorb(setup,sinkData,loc):
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    t    = sinkData['time'][1:]
    orbSep   = (sinkData['rComp'] + sinkData['rAGB'])                                    # cm   ?
    G        = cgs.G                                                             # [cm^3 g^-1 s^-1]
    Mc       = sinkData['massComp']                                                  # [g]
    Md       = sinkData['massAGB']
    sma      = setup['sma_ini']        *cgs.au                               # [cm]  
    v_rel       = np.sqrt(G * (Mc + Md)*(2/orbSep - 1/sma))
    vorbp    = sinkData['v_orbAGB_t']                                                # [cm/s]  ?
    vorbc    = sinkData['v_orbComp_t']#/cgs.kms                                     # [cm/s]
    v_orb_sum = vorbp+vorbc
    v_orb_oldFormula = np.sqrt(G * (Mc + Md)*(1/orbSep))

    ax.plot(t, (v_rel[1:]/cgs.kms),label ='vorb 2/r - 1/a',linestyle = 'solid')
    ax.plot(t, (v_orb_sum[1:]/cgs.kms),label ='vorb sum',linestyle = 'dotted')
    ax.plot(t, (v_orb_oldFormula[1:]/cgs.kms), label ='vorb 1/r',linestyle = 'dashdot')
    ax.set_ylabel(r'$v_{\rm orb}$', fontsize = 12)

    '''
    print('mean relative orb vel correct formula',np.mean(v_rel/cgs.kms))
    # ax.plot(t, v_rel[1:]/cgs.kms, label ='relative vorb', c='k')
    if setup['v_ini'] == 5:
        vw = 12.7
        # plt.hlines(vw,t[0],t[-1],label='vwind')
        ax.plot(t, (v_rel[1:]/cgs.kms)/vw, label ='vorb/vw', c='k')
    elif setup['v_ini'] == 10:
        vw = 15.0
        ax.plot(t, (v_rel[1:]/cgs.kms)/vw, label ='vorb/vw', c='k')
        # plt.hlines(vw,t[0],t[-1],label='vwind')
    elif setup['v_ini'] == 20:
        vw = 22.8
        # plt.hlines(vw,t[0],t[-1],label='vwind')
        ax.plot(t, (v_rel[1:]/cgs.kms)/vw, label ='vorb/vw', c='k')
    
    print('max rel vorb is ',np.max(v_rel/cgs.kms))
    print('min rel vorb is ',np.min(v_rel/cgs.kms))
    print('vw is ',vw)

    ax.set_ylabel(r'$v_{\rm w}/v_{\rm orb}$', fontsize = 12)
    '''
    ax.set_xlabel('time[yrs]', fontsize = 10)

    ax.tick_params(labelsize=10)
    plt.title('vw'+str(setup['v_ini'])+'e'+str(setup['ecc']), fontsize = 15)
    plt.legend()
    fig.tight_layout()
    # plt.savefig(os.path.join(loc, 'png/evolution_vw_vorb.png'))
    # plt.savefig(os.path.join(loc, 'pdf/evolution_vw_vorb.pdf'))    
    plt.savefig(os.path.join(loc, 'png/diff_vorb.png'))
    plt.savefig(os.path.join(loc, 'pdf/diff_vorb.pdf'))

'''
Makes plot of the evolution of the orbital velocity of the AGB star and companion
'''
def plotOrbVel(setup,sinkData, run, loc):
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    t    = sinkData['time'][1:]
    print('mean orb vel AGB',np.mean(sinkData['v_orbAGB_t'][1:]/cgs.kms))
    print('mean orb vel comp',np.mean(sinkData['v_orbComp_t'][1:]/cgs.kms))

    ax.plot(t, sinkData['v_orbAGB_t'][1:]/cgs.kms, label = 'AGB', c='gold')
    ax.plot(t, sinkData['v_orbComp_t'][1:]/cgs.kms, label ='companion', c='crimson')
    
    if setup['triple_star']==True:
        ax.plot(t, sinkData['v_orbComp_in_t'][1:]/cgs.kms, label ='inner companion', c='lime')
        print('mean orb vel comp_in',np.mean(sinkData['v_orbComp_in_t'][1:]/cgs.kms))

    period = setup['period'] / cgs.year
    i = 0
    mini = 0
    if setup['triple_star']==True: #(not at apastron at t=0, so difficult)
        maxi = max(np.max(sinkData['v_orbComp_t'][1:]),np.max(sinkData['v_orbComp_in_t'][1:]))/cgs.kms
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
    ax.set_title('orbital radii', fontsize = 15)
    ra   = sinkData['rAGB'][1:] /cgs.au
    t    = sinkData['time'][1:]

    ax.plot(t, ra, label ='r AGB', c='gold')
    if setup['triple_star']==True:
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax.plot(t, rc_out, label= 'r outer comp',c='crimson')
        ax.plot(t, rc_in, label= 'r inner comp',c='lime')

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

    ax.grid()
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
        ax1.plot(t, ra, label ='r AGB',c='gold')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        rc_in  = sinkData['rComp_in'][1:] /cgs.au
        rc_out = sinkData['rComp'][1:] /cgs.au
        ax2.plot(t, rc_out, label= 'r outer comp',c='crimson')
        ax2.set_ylabel('r out comp [au]', fontsize = 12)
        ax3.plot(t, rc_in, label= 'r inner comp',c='lime')
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
        ax1.plot(t, ra, label ='r AGB',c='gold')
        ax1.set_ylabel('rAGB [au]', fontsize = 12)
        ax2.set_ylabel('r comp [au]', fontsize = 12)
        ax2.plot(t, sinkData['rComp'][1:]/cgs.au, label= 'r comp',c='crimson')
        ax3.set_ylabel('Orb Sep [au]', fontsize = 12)
        ax3.plot(t, ra +sinkData['rComp'][1:]/cgs.au, label = 'Orb sep',c='navy')


    ax3.set_xlabel('time [yrs]', fontsize = 14)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    #plt.legend()
    fig.tight_layout()

    plt.savefig(os.path.join(loc, 'png/evolution_rComp_rAGB_orbSep_seperate.png'))
    plt.savefig(os.path.join(loc, 'pdf/evolution_rComp_rAGB_orbSep_seperate.pdf'))



'''
Plot angular momentum, time derivative of angular momentum, and specific angular momentum
'''
def calcJ(setup,sinkData):
    Jcomp = np.sqrt((sinkData['J_comp'].transpose()[0])**2 + sinkData['J_comp'].transpose()[1]**2 + sinkData['J_comp'].transpose()[2]**2)

    #To get derivative:
        # If you want to reduce to N periods:
    time, index_tmax = reduceToNPeriods(setup,sinkData['time'])
    # Else
    # time = sinkData['time']
    # index_tmax = -1

    # Make empty array to store accretion rates
    derivJ    = np.array([]) 
    # if setup['triple_star']==True:
    #     accrRates_in = np.array([]) 
    #     factor   = 0.5
    # else:
    if setup["sma_ini"] > 20: #wide binary
        factor = 0.5
    else: #close binary
        factor   = 5.0   #1 datapoint every 1/5 years --> 5 datapoints per year

    # Make array of timesteps on which you calculate mass accretion rate
    timeArray    = np.arange(0,int(time[-1]),1./factor)   # 1 datapoint per year/factor
    # Make array with the indices of selected times in original dataset sinkData['time']
    indices_time  = np.array([0])
    J_t_previous = 0
    #For each time in timeArray, calculate the added J wrt previous time, and add it to the derivJ array
    for t in timeArray:
        i            = tl.find_nearest(sinkData['time'], t)  #find the time datapoint in the original data array
        indices_time = np.append(indices_time,i) 
        J_t = Jcomp[indices_time[-1]]
        # print((J_t-J_t_previous)*factor, ' J accreted per year')
        #calculate difference in accreted mass between this and the previous timestep, to find accretion rate
        Jder  = factor*(J_t-J_t_previous)/cgs.year  #added J * factor to get J rate per year # /cgs.yr to get it per s
        # print(Jder,'J accreted per second' )
        derivJ = np.append(derivJ,Jder)
        J_t_previous = J_t

        # if setup['triple_star']==True:
        #     accrRate_in = factor*(sinkData['maccrComp_in'][indices_time[-1]]-sinkData['maccrComp_in'][indices_time[-2]])/cgs.Msun
        #     accrRates_in = np.append(accrRates_in,accrRate_in)

    timeArray2,accrRates = calc_accrRates(setup,sinkData)  
    accrRates_gs = accrRates*cgs.Msun / (cgs.year)                                     # Msun/yr --> g/s 
    # print(np.mean(accrRates), '?')
    j = derivJ/accrRates_gs                                                            # g cm**2 /s**2    /(g/s) --> cm**2/s
    jKepl = np.sqrt(cgs.G * setup['rAccrComp']*cgs.au * sinkData['massComp'])       # cm^3 g^-1 s^-2 *  cm  * g   --> sqrt(cm^4/s^2) --> cm^2 / s
    jKepl2 = np.sqrt(cgs.G * 0.8*setup['rAccrComp']*cgs.au * sinkData['massComp'])
    # est_Rstar = np.mean(j[-100:-1])**2/(cgs.G * np.mean(sinkData['massComp'][-100:-1])) /cgs.au  # in au
    # print('R_*:',est_Rstar)
    # print('R_*/Racc = ',est_Rstar/setup['rAccrComp'])
    # print('test')
    # jKepl_Rstar = np.sqrt(cgs.G * est_Rstar*cgs.au * sinkData['massComp']) 
    # print(np.mean(jKepl_Rstar[-100:-1]),np.mean(j[-100:-1]))
    jK1mean = np.mean(jKepl[-100:-1])
    jK2mean = np.mean(jKepl2[-100:-1])
    jSimMean = np.mean(j[-100:-1])
    print('jKepl: ',jK1mean,' - ',jK2mean)
    print('jsim/jK: ',jSimMean/jK1mean,' - ',jSimMean/jK2mean)


    return Jcomp,derivJ,j,jKepl,timeArray

def plot_angMom(setup,sinkData,loc):
    Jcomp,derivJ,j,jKepl,timeArray = calcJ(setup,sinkData)

    Jcomp_z = sinkData['J_comp'].transpose()[2] # *1e-3 *(1e-2)**2         #g cm**2 /s --> kg m**2 /s
    fig1 = plt.figure(figsize=(8, 5))
    plt.plot(sinkData['time'],Jcomp_z)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel(r'J_z [g cm$^2$ / s]', fontsize = 16)
    fig1.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_Jz.png'))

    fig4 = plt.figure(figsize=(8,5))
    plt.plot(sinkData['time'],Jcomp)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel(r'J [g cm$^2$ / s]', fontsize = 16)
    fig4.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_J.png'))

    fig5 = plt.figure(figsize=(8,5))
    plt.plot(sinkData['massComp']/cgs.Msun,Jcomp)
    plt.xlabel('M_c [Msun]', fontsize = 16)
    plt.ylabel(r'J [g cm$^2$ / s]', fontsize = 16)
    fig5.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_J_Mc.png'))

    fig2 = plt.figure(figsize=(8,5))
    plt.plot(timeArray,derivJ)
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel(r'dJ/dt [g cm$^2$ / s$^2$]', fontsize = 16)
    mini = 0
    maxi = np.max(derivJ)
    plotVerticalApaPerLines(setup,sinkData,mini,maxi)
    fig2.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_dJdt.png'))

    fig3 = plt.figure(figsize=(8,5))
    plt.plot(timeArray,j,label='j')
    plt.plot(sinkData['time'],jKepl,c='k',label=r'j_K')
    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel(r'j [cm$^2$/s]', fontsize = 16)
    plt.legend()
    mini = 0
    maxi = np.max(j)
    plotVerticalApaPerLines(setup,sinkData,mini,maxi)
    fig3.tight_layout()
    plt.savefig(os.path.join(loc, 'png/evolution_j.png'))




'''
Main function executing calculations about orbital evolution
'''
def orbEv_main(run,loc, sinkData, setup,dumpData):
    print('')
    if setup['single_star']:
        return

    else:
        # Visualise the orbit of the system
        plot_orbit(sinkData,setup,loc)
        plotEccentricity(sinkData,setup,loc)


        #For eccentric binaries
        if setup['ecc']>0.09 and setup['triple_star']==False and setup['single_star']==False:
            # Make plot of estimates of orbital period, a and e, using apastron and periastron passages
            plotEstimate_a_Per(sinkData,setup,loc)
            # Make plot of change in orbital separation at apastron and periastron passages
            plotApaPerChange(sinkData,setup,loc)
        
        if setup['triple_star']==False and setup['single_star']==False:
            plot_angMom(setup,sinkData,loc)
        # Plot evolution of the mass accreted by the companion
        plotMassAccr(setup,sinkData, run, loc)
        # Plot evolution of the mass accretion rate by the companion
        # plotMassAccrRate(setup,sinkData, loc)
        plotMassAccrEff(setup,sinkData,loc)
        # Plot evolution of orbital velocities
        plotOrbVel(setup,sinkData, run, loc)
        #Plot vw vs vorb
        plot_vw_vorb(setup,sinkData,loc)
        '''
        # Plot evolution of orbital radii and orbital separation
        #if setup['ecc'] == 0:
        plotOrbRadSeperate(setup,sinkData, run, loc)
        #else:
        plotOrbRad(setup,sinkData, run, loc)
        

        #   UNCOMMENT (PARTS) IF YOU WANT TO WRITE OUT INFO IN TEXT FILE
        # Calculate total mass accreted by companion and AGB star, total mass lost by the AGB star and ratio mass accreted to mass lost
        # Save it in dictionary 'info'
        info = {}
        info['TotMaC']         = sinkData['maccrComp'][-1] # the total mass accreted by the companion
        info['TotMaA']         = sinkData['maccrAGB' ][-1]  # the total mass accreted by the AGB
        info['MassLostAGB']    = sinkData['massAGB'  ][0] - (sinkData['massAGB'][-1] - info['TotMaA'])


        # info['RatioMaC_MLAGB'] = info['TotMaC']/ info['MassLostAGB']
        if setup['triple_star']==True:
            info['TotMaC_in']         = sinkData['maccrComp_in'][-1] # the total mass accreted by the inner companion
            # info['RatioMaC_in_MLAGB'] = info['TotMaC_in']/ info['MassLostAGB']


        # Write text file with this info
        title = os.path.join(loc, 'txt/data_OrbitalEvolution.txt')
        with open (title,'w') as f:
            f.write('\n')
            f.write('Model '+str(run)+'\n')
            f.write('\n')
            f.write('The total mass accreted by the companion is                                    : '+ str(info['TotMaC'      ]/cgs.Msun) +' Msun \n')
            f.write('The total mass accreted by the AGB is                                          : '+ str(info['TotMaA'      ]/cgs.Msun) +' Msun \n')
            f.write('The total mass lost by the AGB sink particle is                                : '+ str(info['MassLostAGB' ]/cgs.Msun) +' Msun \n')
            f.write('The expected mass lost by the AGB is                                           : '+ str(setup['Mdot']*sinkData['time'][-1]) + 'Msun \n')
            f.write('The total mass in the system withouth the mass loss at the boundary is         : '+ str((info['TotMaC'] + np.sum(dumpData['mass']))/cgs.Msun )+ 'Msun \n')
            f.write('The ratio of total mass accreted by companion to total mass lost by AGB is     : '+ str(round((info['TotMaC']/cgs.Msun)/ (setup['Mdot']*sinkData['time'][-1]),5)))
            f.write('\n')
            if setup['triple_star']==True:
                f.write('The total mass accreted by the inner companion is                                    : '+ str(info['TotMaC_in'      ]/cgs.Msun) +' Msun \n')
                f.write('The ratio of the mass accreted by the inner companion to the mass lost by the AGB is : '+ str(round(info['TotMaC_in']/ (setup['Mdot']*sinkData['time'][-1]),5)))
        '''
        '''
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
