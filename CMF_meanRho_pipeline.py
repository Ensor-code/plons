#!/usr/bin/env python
# coding: utf-8

# In[1]:


# !!

#This script only gives the data for 1 model, does not make plots comparing the models, therefore look at phiScriptNewTry.ipynb
def CMF_meanRho(runNumber,directory):

    #Import packages
    import numpy as np
    import os
    import math  as math
    import healpy as hp
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    # import own scripts
    import loadPhantomData as ld

    import scipy.ndimage.filters

    # import certain things from packages
    from collections   import OrderedDict
    from astropy       import constants


    from matplotlib    import rcParams, rc
    # Change the matplotlib default parameters
    rcParams.update({'font.size':   11})
    rcParams.update({'figure.dpi': 200})
    rc('font', family='serif')
    rc('text', usetex=True)


    # In[2]:


    try:
        os.mkdir(directory+'CMF_meanRho/')
    except OSError:
        print('')


    # In[3]:


    # !!
    #returns index where arrayvalue is closest to input value
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx


    # In[4]:


    # !!
    #smoothens values to plot
    def smoothen(values,sigma):
    	#Sigma chosen random, since trying to set sigma=standard_deviations did not work
    	fit=scipy.ndimage.filters.gaussian_filter(values, sigma)
    	return fit


    # In[5]:


    # !!
    # BINNING definitions
    def create_bins(lower_bound, width, quantity):
    #""" create_bins returns an equal-width (distance) partitioning. 
        #It returns an ascending list of tuples, representing the intervals.
        #A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
        #and i < quantity, satisfies the following conditions:
            #(1) bins[i][0] + width == bins[i][1]
            #(2) bins[i-1][0] + width == bins[i][0] and
                #bins[i-1][1] + width == bins[i][1]
    #"""
        bins = []
        for low in range(lower_bound, lower_bound + quantity*width + 1, width):
            bins.append((low, low+width))
        return bins


    def find_bin(value, bins):
    #""" bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        #binning returns the smallest index i of bins so that
        #bin[i][0] <= value < bin[i][1]
    #"""
    # searches the bin containing the value
        if value < bins[0][0]:
            return 0
        elif value > bins[-1][1]:
            return len(bins) - 1
        
        else:
            for i in range(0, len(bins)):
                if bins[i][0] <= value < bins[i][1]:
                    return i
        
            return -1


    # In[6]:


    # !!
    # Calculates cumulative mass and mean density profiles
    # Makes sure the arrays are in the correct order to plot them wrt to the angle
    def getEverything(data): 
        #if single model, remove AGB star data
        # !! CHANGE TO IF SINGLE and uncomment
    #     if key == '10A' or key == '19A' or key == '16A':
    #         theta = data['theta'][:-1]
    #         mass  = data['mass'][:-1]
    #         rho   = data['rho'][:-1]
    #     else: 
        #for the other models, the data without the inner data was loaded in, so stars already removed
        theta = data['theta']
        mass  = data['mass']
        rho   = data['rho']
        
        TotalMass        = np.sum(mass)
        
        #Make array of theta from pi/2 till pi (Orb plane till edge on where z<0)
        indThetaOrdered  = np.argsort(theta)
    #     print(indThetaOrdered)
        thetaOrdered     = theta[indThetaOrdered]
        indPi2Ordered    = find_nearest(thetaOrdered, np.pi/2)
        thetaPi2_Pi      = thetaOrdered[indPi2Ordered:]
        
        #Make array of theta from pi/2 till 0 (Orb plane till edge on where z>0)
        indThetaReversed = np.flip(indThetaOrdered)
        thetaReversed    = theta[indThetaReversed]
        indPi2Reversed   = find_nearest(thetaReversed,np.pi/2)
        thetaPi2_0       = thetaReversed[indPi2Reversed:]
        
        
        # make arrays with mass for those theta values (same indices)
        massPi2_Pi       = mass[indThetaOrdered][indPi2Ordered:]
        massPi2_0        = mass[indThetaReversed][indPi2Reversed:]    
        
        # make arrays with rho for those theta values (same indices)
        rhoPi2_Pi        = rho[indThetaOrdered][indPi2Ordered:]
        rhoPi2_0         = rho[indThetaReversed][indPi2Reversed:]
            
    #     print('okay till here, start making bins')
        # Make bins of theta values , so for theta going from pi/2 till pi-0 by pi/2+-i
        # to have integers, go from pi/2 * 100 = 157 to pi*100 = 314.15
        # the length of such array is pi/2 *100 
        number_of_bins = int(100 * np.pi/2)
        # create the bins from 157 till 314, with width =1, so we bin thetas with width 1/100
        # the created bin is [(157,158),(158,159),...,(313,314)]
        ThetaBins = create_bins(lower_bound=int(100 * np.pi/2), width=1, quantity=number_of_bins )
        # Now we have bins corresponding to the theta values Pi2_Pi
        # Make array containing the corresponding rho and mass values for the bins:
        
    #     print('start making rho and mass matrices')
        rhoBinned  = {}
        massBinned = {}
        for i in range(len(ThetaBins)):
            rhoBinned[i]  = [0]
            massBinned[i] = [0]

        # link the theta values to a bin, by the function find_bin
        for index in range(len(thetaPi2_Pi)):
            #for each theta from pi2_pi, find index of correct bin
            bin_index = find_bin(100*thetaPi2_Pi[index],ThetaBins)
            #and add mass and rho value for that theta to that bin
            rhoBinned[bin_index].append(rhoPi2_Pi[index])
            massBinned[bin_index].append(massPi2_Pi[index])

        
        # the same for pi2_0, here theta = pi2-x corresponds to thetaBin = pi2+x bin
        for index in range(len(thetaPi2_0)):
            theta = thetaPi2_0[index]
            # this theta = pi2 - x , x = pi2-theta, 
            x = np.pi/2 - theta
            # so theta corresponds to the theta in the bins 
            thetaBin = np.pi/2 + x
            bin_index = find_bin(thetaBin*100, ThetaBins)
            #add mass and rho value for that theta to that bin
            rhoBinned[bin_index].append(rhoPi2_0[index])
            massBinned[bin_index].append(massPi2_0[index])
        
        # calculate the mean rho for each thetaBin,
        # but remove the first 0 that we had to add to make initial array
        # make array such that we dont get errors
    #     rhoMean = [0] * len(rhoBinned)
        rhoMean = []
        for index in range(len(rhoBinned)-1):
            rhoMean.append(np.mean(rhoBinned[index][1:]))
        rhoMeanSmoothened = smoothen(rhoMean,4)
        
        # calculate the mass fraction for theta from pi/2 till pi/2+-pi/2
    #     massAccumulated = [0] * len(massBinned)
    #     totalMassPerTheta = [0] * len(massBinned)
    #     massFraction = [0] * len(massBinned)
        massAccumulated   = []
        totalMassPerTheta = []
        massFraction      = []
        # to find theta value where mass fraction is 50% and 75%:
        index50 = 0
        index75 = 0
        
        #the mass for theta = pi/2
        totalMassPerTheta.append(np.sum(massBinned[0]))
        massAccumulated.append(totalMassPerTheta[0])
        massFraction.append(massAccumulated[0]/TotalMass)
        # for each thetaBin calculate total mass f5or that thetaBin, and add it to the mass fraction 
        for index in range(1,len(massBinned)):
            totalMassPerTheta.append(np.sum(massBinned[index]))
            massAccumulated.append(massAccumulated[index-1]+totalMassPerTheta[index])
            massFraction.append(massAccumulated[index]/TotalMass)
            
            if massFraction[index] < 0.5:
                index50 = index50 +1
            if massFraction[index] < 0.75:
                index75 = index75 +1
                
            
        # theta where 50/75 procent of mass is accumulated is upper boundary of the bin corresponding to index50/75
        phi50 = ThetaBins[index50][1]/100
        phi75 = ThetaBins[index75][1]/100
        
        # rhomean where 50 and 75 procent of mass is accumulated
        rho50   = rhoMeanSmoothened[index50]
        rho75   = rhoMeanSmoothened[index75]
        
        x=[]
        for i in range(len(rhoMean)):
    #         x.append(np.pi/2 + i/(np.pi*100)
            x.append(i/(np.pi*100))
            
    #         x.append((np.pi/2 - i/100)/np.pi)

           
        Results = {'phi50'     :  phi50,
                   'phi75'     : phi75,
                   'rho50'       : rho50,
                   'rho75'       : rho75,
                   'massFraction': massFraction,
                   'meanRho'     : rhoMean,
                   'meanRhoSm'   : rhoMeanSmoothened,
                   'x'           : x,
                  }
            
        return Results


    # In[7]:


    # !!

    # CHANGE TO if single:
    if runNumber == '16A' or runNumber == '19A' or runNumber == '10A':
        data  = ld.LoadData_cgs_single_without_inner(str(runNumber))
        dataL = dataModels[str(runNumber)]
        dataR = dataModels[str(runNumber)]
        print('Model ', str(runNumber), 'is loaded.')

    else:
        data  = ld.LoadData_without_inner(str(runNumber), 'all', 'all')
        dataL = ld.LoadData_without_inner(str(runNumber), 'L', 'all')
        dataR = ld.LoadData_without_inner(str(runNumber), 'R', 'all')
        print('Model ', str(runNumber), 'is loaded.')


    # In[8]:


    # !!
    # Calculation of phi, mean rho, mass fractions, ...
    # get info needed to make plots
    # Perc contains all the mass fractions for which we calculate what theta is

    infoForPlot  = getEverything(data)
    infoForPlotL = getEverything(dataL)
    infoForPlotR = getEverything(dataR)


    # In[20]:


    # !!
    # legend:
    left          = mlines.Line2D([],[], color = 'k', linestyle = 'dashed', label = 'Periastron')
    right         = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'Apastron')
    full          = mlines.Line2D([],[], color = 'k', linestyle = 'solid', label = 'Full')

    phi_50 = mlines.Line2D([],[], color = 'k', linestyle = 'None', fillstyle = 'none', marker = 'o', label = '$\\Phi_{50}$')
    phi_75 = mlines.Line2D([],[], color = 'k', linestyle = 'None', marker = 'o', label = '$\\Phi_{75}$')

    handles      = [left, full, right,  phi_50, phi_75]
    handlesphi   = [phi_50, phi_75]


    # In[11]:


    print('The ratio of the mean density in the orbital plane to the mean density on the polar axis of model', runNumber, 'is', round(infoForPlot['meanRhoSm'][0]/infoForPlot['meanRhoSm'][-1],2))


    # In[82]:


    for key in infoForPlot:
        print(key,' is ',infoForPlot[key])


    # In[12]:


    # !!
    # Makes a plot of the mean density profile of both the apastron and periastron side seperately
    def plotLvsR(ax):
        color = 'k'
    #     normalising_factor = infoForPlot['meanRhoSm'][0]
        
        # plot left
        marker = 'dashed'
        log_ysmooth = np.log10(infoForPlotL['meanRhoSm'] ) 

        phi50x =  0.5-((np.pi-infoForPlotL['phi50'])/np.pi)
        phi50y =  np.log10(infoForPlotL['rho50'] )  
        phi75x =  0.5-((np.pi-infoForPlotL['phi75'])/np.pi)
        phi75y =  np.log10(infoForPlotL['rho75'] )   

        ax.plot(infoForPlotL['x'],(log_ysmooth), color, ls = marker)
        ax.plot(phi50x, phi50y, marker = 'o', color = color, mfc = 'None')
        ax.plot(phi75x, phi75y, marker = 'o', color = color) 

        # plot right
        marker = 'dotted'
        log_ysmooth = np.log10(infoForPlotR['meanRhoSm'] )  

        phi50x =  0.5-((np.pi-infoForPlotR['phi50'])/np.pi)
        phi50y =  np.log10(infoForPlotR['rho50'] )  
        phi75x =  0.5-((np.pi-infoForPlotR['phi75'])/np.pi)
        phi75y =  np.log10(infoForPlotR['rho75'] )   

        ax.plot(infoForPlotR['x'],(log_ysmooth), color, ls = marker)
        ax.plot(phi50x, phi50y, marker = 'o', color = color, mfc = 'None')
        ax.plot(phi75x, phi75y, marker = 'o', color = color) 

        # plot the mean
        marker = 'solid'
        log_ysmooth = np.log10(infoForPlot['meanRhoSm'] )  

        phi50x =  0.5-((np.pi-infoForPlot['phi50'])/np.pi)
        phi50y =  np.log10(infoForPlot['rho50'] )   
        phi75x =  0.5-((np.pi-infoForPlot['phi75'])/np.pi)
        phi75y =  np.log10(infoForPlot['rho75'] )   

        ax.plot(infoForPlot['x'], log_ysmooth, color, ls = marker)
        ax.plot(phi50x, phi50y, marker = 'o', color = color, mfc = 'None')
        ax.plot(phi75x, phi75y, marker = 'o', color = color) 

        ax.set_xlim(0,1./2.)

        ax.set_xlabel('$\\Phi$', fontsize = 13)
        ax.set_ylabel('log($<\\rho>$)', fontsize = 13)    


    # In[14]:


    # !! 
    # mean density profile plots, left and right side
    fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
    plotLvsR(ax)

    plt.setp((ax), xticks= [0,1/4,1/2], xticklabels=['$\pi/2$', '$\pi/4$ $\&$ $3\pi/4$', '0 $\&$ $\pi$'])
    ax.tick_params(labelsize=10)

    ax.set_title('Mean density ($\Phi$) (model '+ str(runNumber)+')', fontsize = 15)

    ax.legend(handles = handles, loc = 'lower right',  fontsize = 8)
    plt.savefig(directory+'CMF_meanRho/meanRhoPlot'+str(runNumber))
    print('mean rho plot of model ',str(runNumber), 'ready and saved')


    # In[40]:


    # !! plot of the cumulative mass fraction for all models
    plt.figure(figsize=(5, 4))

    color = 'k'
    marker = 'solid'
    plt.plot(infoForPlot['x'],infoForPlot['massFraction'][:-1],color = color, linestyle = marker)
    phi50x =  0.5-((np.pi-infoForPlot['phi50'])/np.pi)
    phi75x =  0.5-((np.pi-infoForPlot['phi75'])/np.pi)
    plt.plot(phi50x,0.5, color = color, marker='o', mfc = 'None')
    plt.plot(phi75x,0.75, color = color, marker='o')
    print('phi(50) = ', np.round(infoForPlot['phi50']/np.pi ,3), 'pi, phi(75) = ', np.round(infoForPlot['phi75']/np.pi,3),'pi')



    locs, labels = plt.xticks()  # Get the current locations and labels.
    plt.xticks(np.arange(1/2, 0, step=0.1))  # Set label locations.
    plt.xticks([0,1/4,1/2], ['$\pi/2$', '$\pi/4$ $\&$ $3\pi/4$', '0 $\&$  $\pi$'])  

    plt.xlabel('$\Phi$',fontsize = 13)
    plt.ylabel('$M[\Phi]/M_{tot}$', fontsize = 13)
    plt.legend(handles = handlesphi, loc = 'lower right', fontsize = 12)
    plt.title('Cumulative mass fraction (model '+ str(runNumber) + ')', fontsize = 15)
    plt.savefig(directory+'CMF_meanRho/CMFplot'+str(runNumber))
    print('CMF plot of model ',str(runNumber), 'ready and saved')


    # In[61]:


    #Makes text file with all usefull data
    title = directory+'CMF_meanRho/CMF_meanRho_'+str(runNumber)+'.txt'
    with open (title,'w') as f:
        f.write('Model '+str(runNumber)+'\n')
        f.write('\n')
        f.write('Description parameters:'+'\n')
        f.write('\n')
        f.write('phi50 / phi75: Angle wrt orbital plane at which 50 / 75 percent of total mass is acumulated'+'\n')
        f.write('rho50 / rho75: Mean density corresponding to angle phi50 / phi 75'+'\n')
        f.write('Mass fraction: Array with cumulative mass fraction for increasing angle with the orbital plane.'+'\n')
        f.write('               Array corresponds to phi: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('meanRhoSm:     Smoothened array of mean rho. Array corresponds to phi: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('x:             Array going from 0 to 1/2. Corresponds to phi: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('Apastron side / Periastron side / Full :  gives which part of the data is used: x<0 / x>0 / all.'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('Values:'+'\n')
        f.write('\n')
        f.write('The ratio of the mean density in the orbital plane to the mean density on the polar axis is :'+ str(round(infoForPlot['meanRhoSm'][0]/infoForPlot['meanRhoSm'][-1],3))+'\n')
        f.write('This is a usefull ratio to measure the EDE!'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('phi(50) =                         '+ str(np.round(infoForPlot['phi50']/np.pi ,3))+'pi'+'\n')
        f.write('phi(50) to plot on x-axis:        '+str(np.round(phi50x,3))+'\n')
        f.write('rho(50) to plot on y-axis, Full:  '+str(infoForPlot['rho50'])+'\n')
        f.write('                       Apastron:  '+str(infoForPlotL['rho50'])+'\n')
        f.write('                     Periastron:  '+str(infoForPlotR['rho50'])+'\n')    
        f.write('\n')
        f.write('phi(75) =                         '+str(np.round(infoForPlot['phi75']/np.pi,3))+'pi'+'\n')
        f.write('phi(75) to plot on x-axis:        '+str(np.round(phi75x,3))+'\n')
        f.write('rho(75) to plot on y-axis, Full:  '+str(infoForPlot['rho75'])+'\n')
        f.write('                       Apastron:  '+str(infoForPlotL['rho75'])+'\n')
        f.write('                     Periastron:  '+str(infoForPlotR['rho75'])+'\n')    
        
        f.write('\n')
        f.write('x array, to plot: '+'\n')
        f.write(str(infoForPlot['x'])+'\n')
        
        f.write('\n')
        f.write('Full - mass fraction array: '+'\n')
        f.write(str(infoForPlot['massFraction'][:-1])+'\n')
        f.write('\n')
        f.write('Full - meanRhoSm array: '+'\n')
        f.write(str(infoForPlot['meanRhoSm'])+'\n')
        
        f.write('\n')
        f.write('Apastron - mass fraction array: '+'\n')
        f.write(str(infoForPlotL['massFraction'][:-1])+'\n')
        f.write('\n')
        f.write('Apastron - meanRhoSm array: '+'\n')
        f.write(str(infoForPlotL['meanRhoSm'])+'\n')   
        
        f.write('\n')
        f.write('Periastron - mass fraction array: '+'\n')
        f.write(str(infoForPlotR['massFraction'][:-1])+'\n')
        f.write('\n')
        f.write('Perastron - meanRhoSm array: '+'\n')
        f.write(str(infoForPlotR['meanRhoSm'])+'\n')


    # In[ ]:




