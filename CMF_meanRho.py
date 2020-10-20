import numpy                as np
import os
import matplotlib.pyplot    as plt
import matplotlib.lines     as mlines
# import certain things from packages
from collections            import OrderedDict
from astropy                import constants
from matplotlib             import rcParams, rc
# Change the matplotlib default parameters
rcParams.update({'font.size':   11})
rcParams.update({'figure.dpi': 200})
#rc('font', family='serif')
#rc('text', usetex=True)

# own scripts
import Tools                 as tl

# ignore warnings
import warnings
warnings.filterwarnings("ignore")


'''    
Calculates cumulative mass and mean density profiles
Makes sure the arrays are in the correct order to plot them wrt to the angle: pi/2 to 0 & pi
'''
def getEverything(mass, theta, rho): 
    
    TotalMass        = np.sum(mass)
    
    #Make array of theta from pi/2 till pi (Orb plane till edge on where z<0)
    indThetaOrdered  = np.argsort(theta)
    thetaOrdered     = theta[indThetaOrdered]
    indPi2Ordered    = tl.find_nearest(thetaOrdered, np.pi/2)
    thetaPi2_Pi      = thetaOrdered[indPi2Ordered:]
    
    #Make array of theta from pi/2 till 0 (Orb plane till edge on where z>0)
    indThetaReversed = np.flip(indThetaOrdered)
    thetaReversed    = theta[indThetaReversed]
    indPi2Reversed   = tl.find_nearest(thetaReversed,np.pi/2)
    thetaPi2_0       = thetaReversed[indPi2Reversed:]
    
    
    # make arrays with mass for those theta values (same indices)
    massPi2_Pi       = mass[indThetaOrdered][indPi2Ordered:]
    massPi2_0        = mass[indThetaReversed][indPi2Reversed:]    
    
    # make arrays with rho for those theta values (same indices)
    rhoPi2_Pi        = rho[indThetaOrdered][indPi2Ordered:]
    rhoPi2_0         = rho[indThetaReversed][indPi2Reversed:]
        
    # Make bins of theta values , so for theta going from pi/2 till pi-0 by pi/2+-i
    # to have integers, go from pi/2 * 100 = 157 to pi*100 = 314.15
    # the length of such array is pi/2 *100 
    number_of_bins = int(100 * np.pi/2)
    # create the bins from 157 till 314, with width =1, so we bin thetas with width 1/100
    # the created bin is [(157,158),(158,159),...,(313,314)]
    ThetaBins = tl.create_bins(lower_bound=int(100 * np.pi/2), width=1, quantity=number_of_bins )
    # Now we have bins corresponding to the theta values Pi2_Pi
    # Make array containing the corresponding rho and mass values for the bins:
    
    rhoBinned  = {}
    massBinned = {}
    for i in range(len(ThetaBins)):
        rhoBinned[i]  = [0]
        massBinned[i] = [0]

    # link the theta values to a bin, by the function find_bin
    for index in range(len(thetaPi2_Pi)):
        #for each theta from pi2_pi, find index of correct bin
        bin_index = tl.find_bin(100*thetaPi2_Pi[index],ThetaBins)
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
        bin_index = tl.find_bin(thetaBin*100, ThetaBins)
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
    rhoMeanSmoothened = tl.smoothen(rhoMean,4)
    
    # calculate the mass fraction for theta from pi/2 till pi/2+-pi/2
#     massAccumulated = [0] * len(massBinned)
#     totalMassPerTheta = [0] * len(massBinned)
#     massFraction = [0] * len(massBinned)
    massAccumulated   = []
    totalMassPerTheta = []
    massFraction      = []
    # to find theta value where mass fraction is 25%, 50% and 75%:
    index25 = 0
    index50 = 0
    index75 = 0
    
    # the mass for theta = pi/2
    totalMassPerTheta.append(np.sum(massBinned[0]))
    massAccumulated.append(totalMassPerTheta[0])
    massFraction.append(massAccumulated[0]/TotalMass)
    # for each thetaBin calculate total mass f5or that thetaBin, and add it to the mass fraction 
    for index in range(1,len(massBinned)):
        totalMassPerTheta.append(np.sum(massBinned[index]))
        massAccumulated.append(massAccumulated[index-1]+totalMassPerTheta[index])
        massFraction.append(massAccumulated[index]/TotalMass)
        
        if massFraction[index] < 0.25:
            index25 = index25 +1
        if massFraction[index] < 0.5:
            index50 = index50 +1
        if massFraction[index] < 0.75:
            index75 = index75 +1
            
        
    # theta where 25/50/75 procent of mass is accumulated is upper boundary of the bin corresponding to index25/50/75
    theta25 = ThetaBins[index25][1]/100
    theta50 = ThetaBins[index50][1]/100
    theta75 = ThetaBins[index75][1]/100
    
    # rhomean where 50 and 75 procent of mass is accumulated
    rho25   = rhoMeanSmoothened[index25]
    rho50   = rhoMeanSmoothened[index50]
    rho75   = rhoMeanSmoothened[index75]
    
    x=[]
    for i in range(len(rhoMean)):
#         x.append(np.pi/2 + i/(np.pi*100)
        x.append(i/(np.pi*100))
        
#         x.append((np.pi/2 - i/100)/np.pi)

       
    Results = {'theta25'       : theta25,
               'theta50'       : theta50,
               'theta75'       : theta75,
               'rho25'         : rho25,
               'rho50'         : rho50,
               'rho75'         : rho75,
               'massFraction'  : massFraction,
               'meanRho'       : rhoMean,
               'meanRhoSm'     : rhoMeanSmoothened,
               'x'             : x,
              }
        
    return Results



def calcThetaValues(infoForPlot): 

    theta25x =  0.5-((np.pi-infoForPlot['theta25'])/np.pi)
    theta25y =  np.log10(infoForPlot['rho25'] )
    theta50x =  0.5-((np.pi-infoForPlot['theta50'])/np.pi)
    theta50y =  np.log10(infoForPlot['rho50'] )   
    theta75x =  0.5-((np.pi-infoForPlot['theta75'])/np.pi)
    theta75y =  np.log10(infoForPlot['rho75'] )   

    theta    = { '25x'   : theta25x  ,
                 '25y'   : theta25y  ,
                 '50x'   : theta50x  ,
                 '50y'   : theta50y  ,
                 '75x'   : theta75x  ,
                 '75y'   : theta75y  ,
        }

    return theta

   


def calculateDeltaValues(infoForPlot):
    delta25 = infoForPlot['theta25'] - infoForPlot['theta50']
    delta50 = infoForPlot['theta50'] - infoForPlot['theta75']
    ratio = delta25/delta50
    
    delta = { 'delta25'    : delta25,
              'delta50'    : delta50,
              'ratioDelta' : ratio
        }
    
    return delta


'''
Makes a plot of the mean density profile of both the apastron and periastron side seperately
'''
def plotLvsR(ax, theta, infoForPlot, marker):#, infoForPlotL, infoForPlotR):
    color = 'k'
#     normalising_factor = infoForPlot['meanRhoSm'][0]
    
    log_ysmooth = np.log10(infoForPlot['meanRhoSm'] )  

    ax.plot(infoForPlot['x'], log_ysmooth, color, ls = marker)
    ax.plot(theta['25x'], theta['25y'], marker = 'o', color = 'royalblue', fillstyle = 'right')
    ax.plot(theta['50x'], theta['50y'], marker = 'o', color = 'firebrick', fillstyle = 'none' )
    ax.plot(theta['75x'], theta['75y'], marker = 'o', color = 'goldenrod') 

    ax.set_xlim(0,1./2.)

    ax.set_xlabel('$\\theta$', fontsize = 13)
    ax.set_ylabel('log density [cm/g$^3$]', fontsize = 13)   



def CMF_meanRho(run,outloc, data, setup):
    print('')
    print('(5)  Start calculations for the cummulative mass fraction and mean density plots...')



    # legend for plots
    left          = mlines.Line2D([],[], color = 'k', linestyle = 'dashed', label = 'Periastron')
    right         = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'Apastron')
    full          = mlines.Line2D([],[], color = 'k', linestyle = 'solid' , label = 'Full')
    
    theta_25      = mlines.Line2D([],[], color = 'royalblue', linestyle = 'None', fillstyle = 'right', marker = 'o', label = '$\\theta_{25}$')
    theta_50      = mlines.Line2D([],[], color = 'firebrick', linestyle = 'None', fillstyle = 'none' , marker = 'o', label = '$\\theta_{50}$')
    theta_75      = mlines.Line2D([],[], color = 'goldenrod', linestyle = 'None', marker = 'o', label = '$\\theta_{75}$')

    handles       = [left, full, right, theta_25, theta_50, theta_75]
    handlestheta  = [theta_25, theta_50, theta_75]


    # Calculation of theta, mean rho, mass fractions, ...
    # get info needed to make plots
    # Perc contains all the mass fractions for which we calculate what theta is
    infoForPlot  = getEverything(data['mass'], data['theta'], data['rho'])
    theta  = calcThetaValues(infoForPlot)
    delta  = calculateDeltaValues(infoForPlot)

    # Select apastron and periastron side data in case of binary
    if setup['single_star'] == False:
        x = data['position'].transpose()[0]
        thetaR = data['theta'] [x>0]
        massR  = data['mass']  [x>0]
        rhoR   = data['rho']   [x>0]

        thetaL = data['theta'] [x<0]
        massL  = data['mass']  [x<0]
        rhoL   = data['rho']   [x<0]

        infoForPlotL = getEverything(massL, thetaL, rhoL)
        infoForPlotR = getEverything(massR, thetaR, rhoR)
        
        thetaL = calcThetaValues(infoForPlotL)
        thetaR = calcThetaValues(infoForPlotR)
        
        deltaR = calculateDeltaValues(infoForPlotR)
        deltaL = calculateDeltaValues(infoForPlotL)


        # mean density profile plots, left and right side
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        
        plotLvsR(ax, theta , infoForPlot , 'solid' )
        plotLvsR(ax, thetaL, infoForPlotL, 'dashed')
        plotLvsR(ax, thetaR, infoForPlotR, 'dotted')
        
        plt.setp((ax), xticks= [0,1/4,1/2], xticklabels=['$\pi/2$', '$\pi/4$ $&$ $3\pi/4$', '0 $&$ $\pi$'])
        ax.tick_params(labelsize=10)
        ax.set_title('Model '+ str(run)+': Mean density ', fontsize = 15)
        ax.legend(handles = handles, loc = 'lower right',  fontsize = 8)
        fig.tight_layout()
        plt.savefig(outloc+'meanRhoPlot'+str(run))
        print('     Mean density plot of model',str(run),'ready and saved!')    


    # Plot of the cumulative mass fraction 
    fig = plt.figure(figsize=(5, 6))

    color = 'k'
    marker = 'solid'
    plt.plot(infoForPlot['x'],infoForPlot['massFraction'][:-1],color = color, linestyle = marker)
    
    plt.plot(theta['25x'],0.25, color = 'royalblue', marker='o', fillstyle = 'right')
    plt.plot(theta['50x'],0.5 , color = 'firebrick', marker='o', fillstyle = 'none')
    plt.plot(theta['75x'],0.75, color = 'goldenrod', marker='o')

    plt.legend(handles = handlestheta, loc = 'lower right', fontsize = 12)


    if setup['single_star'] == False:
       
        plt.plot(infoForPlotL['x'],infoForPlotL['massFraction'][:-1],color = color, linestyle = 'dashed')
        plt.plot(thetaL['25x'],0.25, color = 'royalblue', marker='o', fillstyle = 'right')
        plt.plot(thetaL['50x'],0.5 , color = 'firebrick', marker='o', fillstyle = 'none')
        plt.plot(thetaL['75x'],0.75, color = 'goldenrod', marker='o')

        plt.plot(infoForPlotR['x'],infoForPlotR['massFraction'][:-1],color = color, linestyle = 'dotted')
        plt.plot(thetaR['25x'],0.25, color = 'royalblue', marker='o', fillstyle = 'right')
        plt.plot(thetaR['50x'],0.5 , color = 'firebrick', marker='o', fillstyle = 'none')
        plt.plot(thetaR['75x'],0.75, color = 'goldenrod', marker='o')

        plt.legend(handles = handles, loc = 'lower right', fontsize = 12)    


    locs, labels = plt.xticks()  # Get the current locations and labels.
    plt.xticks(np.arange(1/2, 0, step=0.1))  # Set label locations.
    plt.xticks([0,1/4,1/2], ['$\pi/2$', '$\pi/4$ $&$ $3\pi/4$', '0 $&$  $\pi$'])  

    plt.xlabel('$\\theta$',fontsize = 13)
    plt.ylabel('$M[\\theta]/M_{tot}$', fontsize = 13)
    plt.title('Model '+ str(run) + ': Cumulative mass fraction', fontsize = 15)
    fig.tight_layout()
    plt.savefig(outloc+'CMFplot'+str(run))
    print('     Cummulative mass fraction plot of model',str(run), 'ready and saved!')


    #Makes text file with all usefull data
    title = outloc+'CMF_meanRho_'+str(run)+'.txt'
    with open (title,'w') as f:
        f.write('Model '+str(run)+'\n')
        f.write('\n')
        f.write('Description parameters:'+'\n')
        f.write('\n')
        f.write('theta25/50/75:     Angle wrt orbital plane at which 25 / 50 / 75 percent of total mass is accumulated'+'\n')
        f.write('rho25/5075   :     Mean density corresponding to angle theta25/50/75, in unit [g/cm^3]\n')
        f.write('delta        :     = (theta25-theta50)/(theta50-theta75) ratio of difference between thetas. \n')
        f.write('                       If normalised to the corresponding single model, delta gives an indication \n')
        f.write('                       of the strength of the equatorial density enhancement (EDE)\n')
        f.write('Mass fraction:     Array with cumulative mass fraction for increasing angle with the orbital plane.'+'\n')
        f.write('                   Array corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('meanRhoSm:         Smoothened array of mean rho. Array corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('x:                 Array going from 0 to 1/2. Corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('Apastron side / Periastron side / Full :  gives which part of the data is used: x<0 / x>0 / all.'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('Values:'+'\n')
        f.write('\n')
        f.write('The ratio of the mean density in the orbital plane to the mean density on the polar axis is: '+ str(round(infoForPlot['meanRhoSm'][0]/infoForPlot['meanRhoSm'][-1],3))+'\n')
        f.write('This is a usefull ratio to measure the EDE!'+'\n')
        f.write('\n')
        f.write('\n')
        # theta 25%
        f.write('theta(25):\n')
        f.write('   Full = '+ str(np.round(infoForPlot['theta25']/np.pi ,3))+'pi'+'\n')
        if setup['single_star'] == False:
            f.write('   Apa  = '+ str(np.round(infoForPlotL['theta25']/np.pi ,3))+'pi'+'\n')
            f.write('   Per  = '+ str(np.round(infoForPlotR['theta25']/np.pi ,3))+'pi'+'\n')
        f.write('theta(25) to plot on x-axis:\n')
        f.write('   Full:  '+str(np.round(theta['25x'],3))+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(np.round(thetaL['25x'],3))+'\n')
            f.write('    Per:  '+str(np.round(thetaR['25x'],3))+'\n')
        f.write('rho(25) to plot on y-axis:\n')
        f.write('   Full:  '+str(infoForPlot['rho25'])+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(infoForPlotL['rho25'])+'\n')
            f.write('    Per:  '+str(infoForPlotR['rho25'])+'\n')    
        f.write('\n')
        # theta 50%
        f.write('theta(50):\n')
        f.write('   Full = '+ str(np.round(infoForPlot['theta50']/np.pi ,3))+'pi'+'\n')
        if setup['single_star'] == False:
            f.write('   Apa  = '+ str(np.round(infoForPlotL['theta50']/np.pi ,3))+'pi'+'\n')
            f.write('   Per  = '+ str(np.round(infoForPlotR['theta50']/np.pi ,3))+'pi'+'\n')
        f.write('theta(50) to plot on x-axis:\n')
        f.write('   Full:  '+str(np.round(theta['50x'],3))+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(np.round(thetaL['50x'],3))+'\n')
            f.write('    Per:  '+str(np.round(thetaR['50x'],3))+'\n')
        f.write('rho(50) to plot on y-axis:\n')
        f.write('   Full:  '+str(infoForPlot['rho50'])+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(infoForPlotL['rho50'])+'\n')
            f.write('    Per:  '+str(infoForPlotR['rho50'])+'\n')    
        f.write('\n')
        # theta 75%
        f.write('theta(75):\n')
        f.write('   Full = '+str(np.round(infoForPlot['theta75']/np.pi,3))+'pi'+'\n')
        if setup['single_star'] == False:
            f.write('    Apa = '+ str(np.round(infoForPlotL['theta75']/np.pi ,3))+'pi'+'\n')
            f.write('    Per = '+ str(np.round(infoForPlotR['theta75']/np.pi ,3))+'pi'+'\n')
        f.write('theta(75) to plot on x-axis:\n')
        f.write('   Full:  '+str(np.round(theta['75x'],3))+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(np.round(thetaL['75x'],3))+'\n')
            f.write('    Per:  '+str(np.round(thetaR['75x'],3))+'\n')
        f.write('rho(75) to plot on y-axis:s\n')
        f.write('   Full:  '+str(infoForPlot['rho75'])+'\n')
        if setup['single_star'] == False:
            f.write('    Apa:  '+str(infoForPlotL['rho75'])+'\n')
            f.write('    Per:  '+str(infoForPlotR['rho75'])+'\n')    
        f.write('\n')
        # delta
        f.write('delta ratio:\n')
        f.write('   Full = '+str(np.round(delta['ratioDelta'],5))+'\n')
        if setup['single_star'] == False:
            f.write('    Apa = '+ str(np.round(deltaL['ratioDelta'] ,5))+'\n')
            f.write('    Per = '+ str(np.round(deltaR['ratioDelta'] ,5))+'\n')
        f.write('\n')
        f.write('\n')
        if setup['single_star'] == False:
            names = ['x array, to plot', 'Full - mass fract', 'Full - meanRhoSm', 'Apastron - mass fract', 'Apastron - meanRhoSm', 'Periastron - mass fract', 'Periastron - meanRhoSm']
            f.write("{: <34} {: <34} {: <34} {: <34} {: <34} {: <34} {: <34}".format(*names))

            col_format = "{:<35}" * 7 + "\n"   # 7 left-justfied columns with 15 character width
            f.write('\n')
            for i in zip(infoForPlot['x'], infoForPlot['massFraction'][:-1], infoForPlot['meanRhoSm'],infoForPlotL['massFraction'][:-1],infoForPlotL['meanRhoSm'], infoForPlotR['massFraction'][:-1], infoForPlotR['meanRhoSm']):
                f.write(col_format.format(*i))
        else:
            names = ['x array, to plot', 'Mass fraction', 'MeanRhoSm']
            col_format = "{:<35}" * 3 + "\n"   # 7 left-justfied columns with 15 character width
            f.write("{: <34} {: <34} {: <34} ".format(*names))
            f.write('\n')
            for i in zip(infoForPlot['x'], infoForPlot['massFraction'][:-1], infoForPlot['meanRhoSm']):
                f.write(col_format.format(*i))
        




