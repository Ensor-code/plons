#!/usr/bin/env python
# coding: utf-8

# In[1]:

def orbEv(runNumber,directory):
    #import packages
    import numpy                 as np
    import matplotlib.pyplot     as plt
    import matplotlib.lines as mlines
    import loadPhantomData       as ld
    import math as math
    import os
    import ConversionFactors_cgs as cgs

    from matplotlib    import rcParams, rc
    # Change the matplotlib default parameters
    rcParams.update({'font.size':   12})
    rcParams.update({'figure.dpi': 200})
    rc('font', family='serif')
    rc('text', usetex=True)
        
    AU = cgs.AU_cm()                    # cm
    vToKms = 1e-3
    timeUnit = 1.5916423E-01            # yrs per timestep
    
    # r is in AU, 
    # M accr  in solar masses     

    try:
        os.mkdir(directory+'orbEvolution/')
    except OSError:
        print('')

    try:
        os.mkdir(directory+'orbEvolution/model'+str(runNumber)+'/')
    except OSError:
        print('')

    

    # In[2]:


    #load in info from windSink files, get t, x, y, z, accreted mass, radii and orbital velocities of the AGB star and companion
     
    def loadFile(fileNumber,runName):
        try:
            (time, x,y,z, mass, maccr,r ,v ) = np.loadtxt(runName, skiprows = 12, usecols = (0,1,2,3,4, 11, 18, 19), unpack = True)
        except OSError:
            print('Converting file to ascii...')
            os.system('cd')
            os.system('cd /home/user/Documents/thesis/ev_files/')
            os.system('splash to ascii M'+fileNumber+'windSink0002N01.ev')
            (time, x,y,z, mass, maccr,r ,v ) = np.loadtxt(runName, skiprows = 12, usecols = (0,1,2,3,4, 11, 18, 19), unpack = True)

        return(time, x,y,z, mass, maccr,r,v)

    def loadIn(runNumber):
        fileNumber = str(runNumber)
        runNameComp = '/home/user/Documents/thesis/ev_files/M'+fileNumber+'windSink0002N01.ev.ascii'
        runNameAGB  = '/home/user/Documents/thesis/ev_files/M'+fileNumber+'windSink0001N01.ev.ascii'
        [timeC,xC,yC,zC, massC, maccrC, rc, vc] = loadFile(fileNumber,runNameComp)
        [timeA,xA,yA,zA, massA, maccrA, rA, vA] = loadFile(fileNumber,runNameAGB)
        
        #a model that was paused during its running, has 2 .ev files per star
        if fileNumber == '19':
            [timeC2, xC2 ,yC2, zC2, massC2, maccrC2, rc2, vc2] = loadFile(fileNumber,  '/home/user/Documents/thesis/ev_files/M'+fileNumber+'windSink0002N02.ev.ascii')
            [timeA2, xA2, yA2, zA2, massA2, maccrA2, rA2, vA2] = loadFile(fileNumber,  '/home/user/Documents/thesis/ev_files/M'+fileNumber+'windSink0001N02.ev.ascii')
            timeC = np.append(timeC,timeC2)
            xC    = np.append(xC,xC2)
            yC    = np.append(yC,yC2)
            zC    = np.append(zC,zC2)
            massC = np.append(massC, massC2)
            maccrC= np.append(maccrC,maccrC2)
            vc    = np.append(vc, vc2)
            
            xA    = np.append(xA,xA2)
            yA    = np.append(yA,yA2)
            zA    = np.append(zA,zA2)
            massA = np.append(massA, massA2)
            maccrA= np.append(maccrA,maccrA2)
            vA    = np.append(vA,vA2)
        
        

        output = {'time'   : timeC * timeUnit,
                  'xC'     : xC,
                  'yC'     : yC,
                  'zC'     : zC,
                  'xA'     : xA,
                  'yA'     : yA,
                  'zA'     : zA,
                  'Mc'     : massC,
                  'MA'     : massA,
                  'MaC'    :  maccrC,
                  'MaA'    :  maccrA,
                  'vc'     :  vc * vToKms,
                  'vA'     :  vA * vToKms
                 }
        
        return output


    # In[3]:


    # calculates radii, orbital separation and r of the COM
    def calcRadii(fileNumber,data):
        d = data 

        t  = d['time']
        x1 = d['xC'] 
        x2 = d['xA'] 
        y1 = d['yC'] 
        y2 = d['yA'] 
        z1 = d['zC'] 
        z2 = d['zA'] 
        M1 = d['Mc']
        M2 = d['MA']
        xCOM = (M1*x1 + M2 * x2)/(M1+M2)
        yCOM = (M1*y1 + M2 * y2)/(M1+M2)
        zCOM = (M1*z1 + M2 * z2)/(M1+M2)
        rCOM = np.sqrt(xCOM**2 + yCOM**2)

        rC = np.sqrt((x1-xCOM)**2 + (y1-yCOM)**2)
        rA = np.sqrt((x2-xCOM)**2 + (y2-yCOM)**2)
        OrbSep = rC+rA

        Output = {'rC'    :rC,
                  'rA'    :rA,
                  'orbSep': OrbSep,
                  'rCOM'  : rCOM }
        
        return Output


    # In[5]:


    # calculate first peak vs last peak heights of eccentric models: 
    # first apastron is in time range 7 - 12 (9.3), first periastron in time range 12-17 (13.95)---> 7-17
    # last apastron is in time range 45-50 (46.5), last periastron in time range 48-53   (51.15)---> 44-54

    # returns first maximum, last maximum, first minimum, last minumimum
    def calcPeaks(fileNumber, parameter):
        #CHANGE PERIOD
        period = 9.30
        
        time   = data['time']
        ecc    = info['ecc']
        orbSep = radii['orbSep']
        param  = radii[parameter]
        
        apastron   = []
        periastron = []  
        timeApa    = []
        timePer    = []
        i=0
        # CHANGE WHILE LOOP!
        while i<6:
            tApa = i*period
            tminApa = tApa - 0.5
            tmaxApa = tApa + 0.5
            Iapa1 = np.abs(time - tminApa).argmin()
            Iapa2 = np.abs(time - tmaxApa).argmin()
            indexApa = (orbSep[Iapa1:Iapa2]).argmax() + Iapa1

            tPer = (i+0.5)*period
            tminPer = tPer - 0.5
            tmaxPer = tPer + 0.5
            Iper1 = np.abs(time - tminPer).argmin()
            Iper2 = np.abs(time - tmaxPer).argmin()
            indexPer = (orbSep[Iper1:Iper2].argmin()) + Iper1
            
            apastron.append(param[indexApa])
            periastron.append(param[indexPer])
            timeApa.append(time[indexApa]*timeUnit)
            timePer.append(time[indexPer]*timeUnit)
            
            i = i+1
            
        return apastron, periastron, timeApa, timePer


    # In[10]:


    # Load in data and info 
    # calculate radii, orbSep, COM
    # calculate tot mass accreted, total mass lost, ratio mass accreted to mass lost
    # calculate peaks of eccentric models

    data = loadIn(runNumber)
    info = ld.LoadInfo(str(runNumber))
    radii= calcRadii(runNumber,data)

    info['TotMaC'] = data['MaC'][-1]
    info['TotMaA'] = data['MaA'][-1]
    info['MassLostAGB']= data['MA'][0] - (data['MA'][-1] - info['TotMaA'])
    info['RatioMaC_MLAGB'] = info['TotMaC']/ info['MassLostAGB']

    #     values of parameters in apastron and periastron 
    info['peaksOrbSep'] = calcPeaks(runNumber, 'orbSep') 
    info['peaksrC']     = calcPeaks(runNumber, 'rC')


    # In[11]:


    # LEGEND
    apastron      = mlines.Line2D([],[], color = 'k', linestyle = 'None',marker = '$a$' ,label = 'Apastron', markersize = 8)
    periastron    = mlines.Line2D([],[], color = 'k', linestyle = 'None',marker = '$p$', label = 'Periastron', markersize =8)

    sma           = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'sma')

    apaLine       = mlines.Line2D([],[], color = 'k', linestyle = 'solid', linewidth = 0.5, label = 'Apastron')
    perLine       = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', linewidth = 0.5, label = 'Periastron')


    maxgr  = mlines.Line2D([],[], color = 'k', linestyle = 'None', marker = '*', label = 'Growth of peak max')
    mingr  = mlines.Line2D([],[], color = 'k', linestyle = 'None', marker = 'o', label = 'Growth of peak min')

    handles1       =[apastron,sma, periastron]
    handles2       =[apaLine,perLine]


    # In[12]:


    # Visualises the orbit of the system
    def plot_orbit(data,radii,ax):
        d = data 
        r = radii 
        t  = d['time']
        rC = r['rC']
        rA = r['rA']
        
        xc = d['xC']
        yc = d['yC']
        zc = d['zC']
        xa = d['xA']
        ya = d['yA']
        za = d['zA']
        
        M1 = d['Mc']
        M2 = d['MA']
        xCOM = (M1*xc + M2 * xa)/(M1+M2)
        yCOM = (M1*yc + M2 * ya)/(M1+M2)
        rCOM = np.sqrt(xCOM**2 + yCOM**2)
        
        ax.axis('equal')
        ax.set_ylim(-4,4)
        
        ax.plot(xc,yc, c = 'r', label = 'companion')
        ax.plot(xa,ya, c = 'b', label = 'AGB')
        ax.set_ylabel('y[AU]',fontsize = 15)
        ax.set_xlabel('x[AU]',fontsize = 15)
        ax.plot(0,0,'*', color = 'k', label = 'COM')
        ecc = str(info['ecc'])
        ax.set_title('e = '+ecc+', $a = 6$ AU, $q=1.5$', fontsize = 15)

        ax.tick_params(labelsize=12)


    # In[13]:


    #Plots orbits
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    fig.set_size_inches(7, 7)

    plot_orbit(data,radii, ax)
    ax.axis('equal')
    ax.legend(fontsize = 15, loc = 'center right')
    plt.savefig(directory+'orbEvolution/model'+str(runNumber)+'/Orbit_model'+str(runNumber))


    # In[14]:


    #calculates sma at timesteps of the input peaks 'infoPeaksOrbSep'
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


    # In[53]:


    #Makes plots of the change of the apastron and periastron values of a selected parameter
    def plotChangeParam(peaksPar):#, ylabel, unit, name, title):
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        
        fig.set_size_inches(7, 5)
        c='k'
        toPlot = info[peaksPar]
        [sma,sma_t] = calcSMA(toPlot)
            
        t_total = max(data['time'])

        #!CHANGE 6 TO SMA INI
        
        delta_a = ((toPlot[0]-toPlot[0][0] + (toPlot[1]-toPlot[1][0]) )/6 )[-1]
        ratio_delta_a_per_yr = (delta_a/6)/ t_total

        delta_e = ((toPlot[0]-toPlot[0][0] - (toPlot[1]-toPlot[1][0]) )/(12 + 2*delta_a) )[-1]
        ratio_delta_e_per_yr = (delta_e/info['ecc'] ) /t_total

        print('delta_e/e= ', str(round(ratio_delta_e_per_yr,9)), '/yr')
        print('delta_a/a= ', str(round(ratio_delta_a_per_yr,9)), '/yr')


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
        ax.set_title('Orbital evolution model'+str(runNumber))
        ax.legend(handles = handles1, fontsize = 12)#, loc = 'lower left')
        plt.savefig(directory+'orbEvolution/model'+str(runNumber)+'/ChangeOrbSep_'+str(runNumber))
        
        #Write text file with usefull info
        title = directory+'orbEvolution/model'+str(runNumber)+'/info_OrbEvol_'+str(runNumber)+'.txt'
        with open (title,'w') as f:
            f.write('Model '+str(runNumber)+'\n')
            f.write('\n')
            f.write('The total mass accreted by the companion is                                    : '+ str(info['TotMaC']) +'\n')
            f.write('The total mass accreted by the AGB is                                          : '+ str(info['TotMaA']) +'\n')
            f.write('The total mass lost by the AGB is                                              : '+ str(info['MassLostAGB']) +'\n')
            f.write('The ratio of the mass accreted by the companion to the mass lost by the AGB is : '+ str(info['RatioMaC_MLAGB']))
            f.write('\n')
            f.write('\n')
            f.write('The change in eccentricity, delta_e/e=         '+ str(round(ratio_delta_e_per_yr,9))+ '/yr'+'\n')
            f.write('The change in sma,          delta_a/a=         '+ str(round(ratio_delta_a_per_yr,9))+ '/yr'+'\n')
            f.write('\n')
            f.write('To plot the change in eccentricity and change in sma, use the following: \n')
            f.write('Times of periastron passage, to plot on x axis: \n'+ str(sma_t)+'\n' )
            f.write('Orbital separation at apastron passages : \n'+ str(toPlot[0])+ '\n')
            f.write('Orbital separation at periastron passages : \n'+ str(toPlot[1])+ '\n')
            f.write('Semi major axis at times of periastron passage : \n'+ str(sma)+ '\n')


    # In[52]:


    #Make the plots of the change in OrbSep

    #ONLY USEFULL IF ECC>0
    plotChangeParam('peaksOrbSep')


    # In[143]:


    # make plot of the mass accretion evolution, very interesting plot!
    plt.figure(figsize=(8, 5))
    #to scale:
    maxi = max(data['MaC'])
    mini = min(data['MaC'])

    plt.plot(data['time'],  data['MaC'], color = 'royalblue', linestyle = 'solid')

    j = 9.3/2
    i = 0
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)
    i = i+9.3
    j = j+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)
    plt.vlines(j,mini, maxi,  linestyle = 'dotted', linewidth = 0.5)

    i = i+9.3
    plt.vlines(i,mini, maxi,  linestyle = 'solid', linewidth = 0.5)


    ax = plt.subplot(111)

    plt.xlabel('Time[yrs]', fontsize = 16)
    plt.ylabel('Total accreted mass', fontsize = 16)
    #GIVE CORRECT UNIT

    plt.title('Total accreted mass evolution', fontsize = 18)
    plt.savefig(directory+'orbEvolution/model'+str(runNumber)+'/MacrEvolution_'+str(runNumber))


    # In[166]:


    # Plot of orbital evolution of non eccentric models
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})

    ax.plot(data['time'], radii['rC'], label= 'r comp')
    ax.set_ylabel('$r$[AU]', fontsize = 12)
    ax.set_title('r comp, r AGB, orb sep (model'+str(runNumber)+')', fontsize = 15)
       
    ax.plot(data['time'], radii['rA'], label ='r AGB')
    ax.plot(data['time'], radii['orbSep'], label = 'Orb sep')
    ax.set_xlabel('time[yrs]', fontsize = 14)
    ax.tick_params(labelsize=10)

    plt.legend()
    plt.savefig(directory+'orbEvolution/model'+str(runNumber)+'/rc,rA,orbSep_'+str(runNumber))


    # In[163]:


    # Plot of orbital velocity of non eccentric models
    fig, ((ax))= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    # fig.set_size_inches(16, 10)
    fig.suptitle('Orbital velocities (model'+str(runNumber)+')', fontsize = 15)


    ax.plot(data['time'], data['vc'], label ='companion')
    ax.plot(data['time'], data['vA'], label = 'AGB')
    ax.set_ylabel('$v_{orb}$ [km/s]', fontsize = 12)

    ax.set_xlabel('time[yrs]', fontsize = 10)

    ax.tick_params(labelsize=10)
          
    plt.legend()
    plt.savefig(directory+'orbEvolution/model'+str(runNumber)+'/orbVel_'+str(runNumber))
          

