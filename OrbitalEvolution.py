import numpy                    as np
import matplotlib.pyplot        as plt
import matplotlib.lines         as mlines
import math                     as math
import os
from matplotlib                 import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size' :   12})
rcParams.update({'figure.dpi': 200})
rc('font', family='serif')
rc('text', usetex=True)

# own scripts
import ConversionFactors_cgs    as cgs
    

'''
Calculates orbital separation and r of the COM, used sinkData (data)
'''
def calcRadii(data):
    t  = data['time']
    x1 = data['posComp'].transpose()[0] 
    x2 = data['posAGB' ].transpose()[0]
    y1 = data['posComp'].transpose()[1]
    y2 = data['posAGB' ].transpose()[1]
    z1 = data['posComp'].transpose()[2]
    z2 = data['posAGB' ].transpose()[2]
    M1 = data['massComp']
    M2 = data['massAGB' ]
    xCOM = (M1*x1 + M2 * x2)/(M1+M2)
    yCOM = (M1*y1 + M2 * y2)/(M1+M2)
    zCOM = (M1*z1 + M2 * z2)/(M1+M2)
    rCOM = np.sqrt(xCOM**2 + yCOM**2)

    rC = data['rComp']
    rA = data['rAGB']
    OrbSep = rC+rA

    output = {'orbSep': OrbSep,
              'rCOM'  : rCOM }
    
    return output

'''
Calculate first peak vs last peak heights of eccentric models: 
    - first apastron is in time range 7 - 12 (9.3), first periastron in time range 12-17 (13.95)---> 7-17
    - last apastron is in time range 45-50 (46.5), last periastron in time range 48-53   (51.15)---> 44-54

returns first maximum, last maximum, first minimum, last minumimum
uses sinkData = data
'''
def calcPeaksOrbSep(data, setup, radii):

    period = setup['period_ini'] * cgs.sec_year()
    time   = data['time']
    ecc    = setup['ecc']
    orbSep = radii['orbSep']    
    # print('time goes from ', min(time),' to ', max(time))
    # print('period is ', period)
    # print('the amount of periods is max(time)/period: ', math.ceil(max(time)/period))

    apastron   = []
    periastron = []  
    timeApa    = []
    timePer    = []
    i=0

    while i<math.ceil(max(time)/period):
        tApa     = i*period
        tminApa  = tApa - 0.5
        tmaxApa  = tApa + 0.5
        Iapa1    = np.abs(time - tminApa).argmin()
        Iapa2    = np.abs(time - tmaxApa).argmin()
        indexApa = (orbSep[Iapa1:Iapa2]).argmax() + Iapa1

        tPer     = (i+0.5)*period
        tminPer  = tPer - 0.5
        tmaxPer  = tPer + 0.5
        Iper1    = np.abs(time - tminPer).argmin()
        Iper2    = np.abs(time - tmaxPer).argmin()
        indexPer = (orbSep[Iper1:Iper2].argmin()) + Iper1
        
        apastron.append(orbSep[indexApa])
        periastron.append(orbSep[indexPer])
        timeApa.append(time[indexApa])#*timeUnit)
        timePer.append(time[indexPer])#*timeUnit)
        
        i = i+1
        
    return apastron, periastron, timeApa, timePer

'''
Visualises the orbit of the system
'''
def plot_orbit(data, setup, radii, ax): 
    r  = radii 
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
    ax.plot(xa,ya, c = 'b', label = 'AGB')
    ax.set_ylabel('y[cm]',fontsize = 15)
    ax.set_xlabel('x[cm]',fontsize = 15)
    ax.plot(0,0,'*', color = 'k', label = 'COM')
    ax.set_title('$e = $'+str(setup['ecc'])+', $a =$'+ str(setup['sma_ini'])+' AU, $q=$'+ str(setup['massAGB_ini']/setup['massComp_ini']), fontsize = 15)

    ax.tick_params(labelsize=12)

'''
Calculates semi-major axis at timesteps of the input peaks 'infoPeaksOrbSep'
'''
def calcSMA(infoPeaksOrbSep):
    i = infoPeaksOrbSep
#     print(i)
    sma  = []
    tsma = []
    for j in range(len(i[0])):
#         print(j)
        sma.append((i[0][j]+i[1][j])/2)
        tsma.append((i[2][j]+i[3][j])/2)
    return (sma,tsma)

'''
Makes plots of the change of the apastron and periastron values of a selected parameter
'''
def plotChangeOrbSep(info, sinkData, setup, peaksPar, run, loc):#, ylabel, unit, name, title):
    
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    
    # LEGEND
    apastron      = mlines.Line2D([],[], color = 'k', linestyle = 'None',marker = '$a$' ,label = 'Apastron', markersize = 8)
    periastron    = mlines.Line2D([],[], color = 'k', linestyle = 'None',marker = '$p$', label = 'Periastron', markersize =8)
    sma           = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'sma')
    handles1      = [apastron,sma, periastron]


    fig.set_size_inches(7, 5)
    c = 'k'
    toPlot      = info[peaksPar]
    [sma,sma_t] = calcSMA(toPlot)
        
    t_total     = max(sinkData['time']) #in years

    #detla a = detla r_apa + delta r_per / 2 
    delta_a              = ((toPlot[0]-toPlot[0][0] + (toPlot[1]-toPlot[1][0]) )/2)[-1]
    ratio_delta_a_per_yr = (delta_a/(setup['sma_ini']*cgs.AU_cm()))/ t_total

    delta_e              = ((toPlot[0]-toPlot[0][0] - (toPlot[1]-toPlot[1][0]) )/(2*setup['sma_ini']*cgs.AU_cm() + 2*delta_a) )[-1]
    ratio_delta_e_per_yr = (delta_e/setup['ecc'] ) /t_total



    ax.plot(sma_t,toPlot[0]-toPlot[0][0],color = c, marker = '$a$', linestyle = 'dashed', markersize = 10)
    ax.plot(sma_t,toPlot[1]-toPlot[1][0],color = c, marker = '$p$', linestyle = 'dashed', markersize = 10)
    for i in range(len(sma_t)):
        ax.plot(sma_t[i],toPlot[0][i]-toPlot[0][0],color = 'white', marker = 'o', markersize = 12)
        ax.plot(sma_t[i],toPlot[0][i]-toPlot[0][0],color = c, marker = '$a$', markersize = 11)
        ax.plot(sma_t[i],toPlot[1][i]-toPlot[1][0],color = 'white', marker = 'o', markersize = 12)
        ax.plot(sma_t[i],toPlot[1][i]-toPlot[1][0],color = c, marker = '$p$', markersize = 11)

    ax.plot(sma_t, sma-sma[0], color = c, linestyle = 'dotted', markersize = 10)
    
    ax.tick_params(labelsize=12)     
    
    ax.set_xlabel('Orbit', fontsize = 18)
    ax.set_ylabel('$\Delta$Orb sep [AU]', fontsize = 16)
    ax.set_title('Orbital evolution model'+str(run))
    ax.legend(handles = handles1, fontsize = 12)#, loc = 'lower left')
    plt.savefig(loc+'orbEvolution/model'+str(run)+'/ChangeOrbSep_'+str(run))
    
    #Write text file with usefull info
    title = loc+'orbEvolution/model'+str(run)+'/info_OrbEvol_'+str(run)+'.txt'
    with open (title,'a') as f:
        f.write('\n')
        f.write('The change in eccentricity, delta_e/e=         '+ str(round(ratio_delta_e_per_yr,9))+ '/yr'+'\n')
        f.write('The change in sma,          delta_a/a=         '+ str(round(ratio_delta_a_per_yr,9))+ '/yr'+'\n')
        f.write('\n')

        f.write('To plot the change in eccentricity and change in sma, use the following: '+ '\n')
        f.write('\n')
        names = ['Periastron passages [yrs]','Orb Sep Apa [cm]','Orb Sep Per [cm]','sma Per [cm]']
        f.write("{: <34} {: <34} {: <34} {: <34}".format(*names))
        f.write('\n')
        col_format = "{:<35}" * 4 + "\n"   # 7 left-justfied columns with 15 character width

        for i in zip(sma_t, toPlot[0], toPlot[1],sma):
            f.write(col_format.format(*i))
        f.write('\n')


def plotMassAccr(sinkData, run, loc):
  # make plot of the mass accretion evolution, very interesting to plot!
    plt.figure(figsize=(8, 5))
    #to scale:
    maxi = max(sinkData['maccrComp'])
    mini = min(sinkData['maccrComp'])

    plt.plot(sinkData['time'],  sinkData['maccrComp'], color = 'royalblue', linestyle = 'solid')

    j = 9.3/2
    i = 0
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid' , linewidth = 0.5)


    ax = plt.subplot(111)

    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Total accreted mass [g]', fontsize = 16)

    plt.title('Total accreted mass evolution', fontsize = 18)
    plt.savefig(loc+'orbEvolution/model'+str(run)+'/MacrEvolution_'+str(run))



def plotOrbVelEcc(sinkData, run, loc):
    # Plot of orbital velocity of non eccentric models
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    # fig.set_size_inches(16, 10)
    fig.suptitle('Orbital velocities (model'+str(run)+')', fontsize = 15)

    ax.plot(sinkData['time'], sinkData['v_orbComp_t'], label ='companion')
    ax.plot(sinkData['time'], sinkData['v_orbAGB_t'], label = 'AGB')
    ax.set_ylabel('$v_{orb}$ [cm/s]', fontsize = 12)

    ax.set_xlabel('time[yrs]', fontsize = 10)

    ax.tick_params(labelsize=10)
          
    plt.legend()
    plt.savefig(loc+'orbEvolution/model'+str(run)+'/orbVel_'+str(run))



def plotOrbEvEcc(sinkData, run, loc):
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})

    ax.plot(sinkData['time'], sinkData['rComp'], label= 'r comp')
    ax.set_ylabel('$r$[cm]', fontsize = 12)
    ax.set_title('r comp, r AGB, orb sep (model'+str(run)+')', fontsize = 15)
       
    ax.plot(sinkData['time'], sinkData['rAGB'], label ='r AGB')
    ax.plot(sinkData['time'], sinkData['rAGB']+sinkData['rComp'], label = 'Orb sep')
    ax.set_xlabel('time[yrs]', fontsize = 14)
    ax.tick_params(labelsize=10)

    plt.legend()
    plt.savefig(loc+'orbEvolution/model'+str(run)+'/rc,rA,orbSep_'+str(run))


def orbEv_main(run,loc, sinkData, setup):
    print('')
    print('(6)  Start calculations for orbital evolution...')
    if setup['single_star'] == True:
        print('A single model has no orbit, and thereby no orbital evolution.')
        print('The orbital evolution part is skipped.')
    else:
        try:
            os.mkdir(loc+'orbEvolution/')
        except OSError:
            print('')

        try:
            os.mkdir(loc+'orbEvolution/model'+str(run)+'/')
        except OSError:
            print('')

        # calculate orbSep, COM
        # calculate tot mass accreted, total mass lost, ratio mass accreted to mass lost
        # calculate peaks of eccentric models

        info = {}
        radii= calcRadii(sinkData)

        info['TotMaC'] = sinkData['maccrComp'][-1] # the total mass accreted by the companion
        info['TotMaA'] = sinkData['maccrAGB'][-1]  # the total mass accreted by the AGB
        info['MassLostAGB'] = sinkData['massAGB'][0] - (sinkData['massAGB'][-1] - info['TotMaA'])
        info['RatioMaC_MLAGB'] = info['TotMaC']/ info['MassLostAGB']

        #Write text file with usefull info
        title = loc+'orbEvolution/model'+str(run)+'/info_OrbEvol_'+str(run)+'.txt'
        with open (title,'w') as f:
            f.write('\n')
            f.write('Model '+str(run)+'\n')
            f.write('\n')
            f.write('The total mass accreted by the companion is                                    : '+ str(info['TotMaC']) +' g \n')
            f.write('The total mass accreted by the AGB is                                          : '+ str(info['TotMaA']) +' g \n')
            f.write('The total mass lost by the AGB is                                              : '+ str(info['MassLostAGB']) +' g \n')
            f.write('The ratio of the mass accreted by the companion to the mass lost by the AGB is : '+ str(round(info['RatioMaC_MLAGB'],5)))
            f.write('\n')

        #Plots orbit
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        fig.set_size_inches(7, 7)

        plot_orbit(sinkData,setup,radii, ax)
        ax.axis('equal')
        ax.legend(fontsize = 15, loc = 'center right')
        plt.savefig(loc+'orbEvolution/model'+str(run)+'/Orbit_model'+str(run))


        if setup['ecc'] >0:
            # Calculate values of parameters in apastron and periastron 
            info['peaksOrbSep'] = calcPeaksOrbSep(sinkData, setup, radii) 
            #Make the plots of the change in OrbSep and eccentrictiy
            plotChangeOrbSep(info, sinkData, setup, 'peaksOrbSep',run,loc)
            #Plot evolution of the mass accretion
            plotMassAccr(sinkData, run, loc)
            #Plot orbital velocities 
            plotOrbVelEcc(sinkData, run, loc)
            #Plot orbital radii and orbital separation
            plotOrbEvEcc(sinkData, run, loc)


        #Write text file with usefull info
        title = loc+'orbEvolution/model'+str(run)+'/info_OrbEvol_'+str(run)+'.txt'
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

