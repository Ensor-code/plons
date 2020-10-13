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
rc('font', family='serif')
rc('text', usetex=True)

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
    # to find theta value where mass fraction is 50% and 75%:
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
        
        if massFraction[index] < 0.5:
            index50 = index50 +1
        if massFraction[index] < 0.75:
            index75 = index75 +1
            
        
    # theta where 50/75 procent of mass is accumulated is upper boundary of the bin corresponding to index50/75
    theta50 = ThetaBins[index50][1]/100
    theta75 = ThetaBins[index75][1]/100
    
    # rhomean where 50 and 75 procent of mass is accumulated
    rho50   = rhoMeanSmoothened[index50]
    rho75   = rhoMeanSmoothened[index75]
    
    x=[]
    for i in range(len(rhoMean)):
#         x.append(np.pi/2 + i/(np.pi*100)
        x.append(i/(np.pi*100))
        
#         x.append((np.pi/2 - i/100)/np.pi)

       
    Results = {'theta50'       : theta50,
               'theta75'       : theta75,
               'rho50'         : rho50,
               'rho75'         : rho75,
               'massFraction'  : massFraction,
               'meanRho'       : rhoMean,
               'meanRhoSm'     : rhoMeanSmoothened,
               'x'             : x,
              }
        
    return Results



def calcThetaValues(infoForPlot, infoForPlotL, infoForPlotR):
    theta50xL =  0.5-((np.pi-infoForPlotL['theta50'])/np.pi)
    theta50yL =  np.log10(infoForPlotL['rho50'] )  
    theta75xL =  0.5-((np.pi-infoForPlotL['theta75'])/np.pi)
    theta75yL =  np.log10(infoForPlotL['rho75'] ) 

    theta50xR =  0.5-((np.pi-infoForPlotR['theta50'])/np.pi)
    theta50yR =  np.log10(infoForPlotR['rho50'] )  
    theta75xR =  0.5-((np.pi-infoForPlotR['theta75'])/np.pi)
    theta75yR =  np.log10(infoForPlotR['rho75'] )   

    theta50x =  0.5-((np.pi-infoForPlot['theta50'])/np.pi)
    theta50y =  np.log10(infoForPlot['rho50'] )   
    theta75x =  0.5-((np.pi-infoForPlot['theta75'])/np.pi)
    theta75y =  np.log10(infoForPlot['rho75'] )   

    theta    = { '50x'   : theta50x  ,
                 '50y'   : theta50y  ,
                 '75x'   : theta75x  ,
                 '75y'   : theta75y  ,
                 '50xR'  : theta50xR ,
                 '50yR'  : theta50yR ,
                 '75xR'  : theta75xR ,
                 '75yR'  : theta75yR ,
                 '50xL'  : theta50xL ,
                 '50yL'  : theta50yL ,
                 '75xL'  : theta75xL ,
                 '75yL'  : theta75yL
        }

    return theta



def calcThetaValuesSingle(infoForPlot):

    theta50x =  0.5-((np.pi-infoForPlot['theta50'])/np.pi)
    theta50y =  np.log10(infoForPlot['rho50'] )   
    theta75x =  0.5-((np.pi-infoForPlot['theta75'])/np.pi)
    theta75y =  np.log10(infoForPlot['rho75'] )   

    theta    = {'50x'  : theta50x,
                '50y'  : theta50y,
                '75x'  : theta75x,
                '75y'  : theta75y,
        }

    return theta    


'''
Makes a plot of the mean density profile of both the apastron and periastron side seperately
'''
def plotLvsR(ax, theta, infoForPlot, infoForPlotL, infoForPlotR):
    color = 'k'
#     normalising_factor = infoForPlot['meanRhoSm'][0]
    
    # plot left
    marker = 'dashed'
    log_ysmooth = np.log10(infoForPlotL['meanRhoSm'] ) 

    ax.plot(infoForPlotL['x'],(log_ysmooth), color, ls = marker)
    ax.plot(theta['50xL'], theta['50yL'], marker = 'o', color = color, mfc = 'None')
    ax.plot(theta['75xL'], theta['75yL'], marker = 'o', color = color) 

    # plot right
    marker = 'dotted'
    log_ysmooth = np.log10(infoForPlotR['meanRhoSm'] )  


    ax.plot(infoForPlotR['x'],(log_ysmooth), color, ls = marker)
    ax.plot(theta['50xR'], theta['50yR'], marker = 'o', color = color, mfc = 'None')
    ax.plot(theta['75xR'], theta['75yR'], marker = 'o', color = color) 

    # plot the mean
    marker = 'solid'
    log_ysmooth = np.log10(infoForPlot['meanRhoSm'] )  

    ax.plot(infoForPlot['x'], log_ysmooth, color, ls = marker)
    ax.plot(theta['50x'], theta['50y'], marker = 'o', color = color, mfc = 'None')
    ax.plot(theta['75x'], theta['75y'], marker = 'o', color = color) 

    ax.set_xlim(0,1./2.)

    ax.set_xlabel('$\\theta$', fontsize = 13)
    ax.set_ylabel('log($<\\rho>$)', fontsize = 13)   



def CMF_meanRho(run,loc, data, setup):
    print('')
    print('(5)  Start calculations for the cummulative mass fraction and mean density plots...')

    #Make new folder to save plots and text file, unless folder exists
    try:
        os.mkdir(loc+'CMF_meanRho/')
    except OSError:
        print('')


    # legend:
    left          = mlines.Line2D([],[], color = 'k', linestyle = 'dashed', label = 'Periastron')
    right         = mlines.Line2D([],[], color = 'k', linestyle = 'dotted', label = 'Apastron')
    full          = mlines.Line2D([],[], color = 'k', linestyle = 'solid', label = 'Full')

    theta_50       = mlines.Line2D([],[], color = 'k', linestyle = 'None', fillstyle = 'none', marker = 'o', label = '$\\theta_{50}$')
    theta_75       = mlines.Line2D([],[], color = 'k', linestyle = 'None', marker = 'o', label = '$\\theta_{75}$')

    handles      = [left, full, right,  theta_50, theta_75]
    handlestheta   = [theta_50, theta_75]


    # Calculation of theta, mean rho, mass fractions, ...
    # get info needed to make plots
    # Perc contains all the mass fractions for which we calculate what theta is
    infoForPlot  = getEverything(data['mass'], data['theta'], data['rho'])

    #select apastron and periastron side data in case of binary
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
        theta = calcThetaValues(infoForPlot, infoForPlotL, infoForPlotR)


        # mean density profile plots, left and right side
        fig, (ax)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        plotLvsR(ax, theta, infoForPlot, infoForPlotL, infoForPlotR)
        plt.setp((ax), xticks= [0,1/4,1/2], xticklabels=['$\pi/2$', '$\pi/4$ $\&$ $3\pi/4$', '0 $\&$ $\pi$'])
        ax.tick_params(labelsize=10)
        ax.set_title('Mean density ($\\theta$) (model '+ str(run)+')', fontsize = 15)
        ax.legend(handles = handles, loc = 'lower right',  fontsize = 8)
        plt.savefig(loc+'CMF_meanRho/meanRhoPlot'+str(run))
        print('     Mean density plot of model',str(run),'ready and saved!')    
 
    else:
        theta = calcThetaValuesSingle(infoForPlot)


    # plot of the cumulative mass fraction 
    plt.figure(figsize=(5, 4))

    color = 'k'
    marker = 'solid'
    plt.plot(infoForPlot['x'],infoForPlot['massFraction'][:-1],color = color, linestyle = marker)


    if setup['single_star'] == False:
        plt.plot(theta['50x'],0.5, color = color, marker='o', mfc = 'None')
        plt.plot(theta['75x'],0.75, color = color, marker='o')

        plt.plot(infoForPlotL['x'],infoForPlotL['massFraction'][:-1],color = color, linestyle = 'dashed')
        plt.plot(theta['50xL'],0.5, color = color, marker='o', mfc = 'None')
        plt.plot(theta['75xL'],0.75, color = color, marker='o')

        plt.plot(infoForPlotR['x'],infoForPlotR['massFraction'][:-1],color = color, linestyle = 'dotted')

        plt.plot(theta['50xR'],0.5, color = color, marker='o', mfc = 'None')
        plt.plot(theta['75xR'],0.75, color = color, marker='o')

        plt.legend(handles = handles, loc = 'lower right', fontsize = 12)


    else:
        plt.plot(theta['50x'],0.5, color = color, marker='o', mfc = 'None')
        plt.plot(theta['75x'],0.75, color = color, marker='o')
        plt.legend(handles = handlestheta, loc = 'lower right', fontsize = 12)


    locs, labels = plt.xticks()  # Get the current locations and labels.
    plt.xticks(np.arange(1/2, 0, step=0.1))  # Set label locations.
    plt.xticks([0,1/4,1/2], ['$\pi/2$', '$\pi/4$ $\&$ $3\pi/4$', '0 $\&$  $\pi$'])  

    plt.xlabel('$\\theta$',fontsize = 13)
    plt.ylabel('$M[\\theta]/M_{tot}$', fontsize = 13)
    plt.title('Cumulative mass fraction (model '+ str(run) + ')', fontsize = 15)
    plt.savefig(loc+'CMF_meanRho/CMFplot'+str(run))
    print('     Cummulative mass fraction plot of model',str(run), 'ready and saved!')


    #Makes text file with all usefull data
    title = loc+'CMF_meanRho/CMF_meanRho_'+str(run)+'.txt'
    with open (title,'w') as f:
        f.write('Model '+str(run)+'\n')
        f.write('\n')
        f.write('Description parameters:'+'\n')
        f.write('\n')
        f.write('theta50 / theta75: Angle wrt orbital plane at which 50 / 75 percent of total mass is acumulated'+'\n')
        f.write('rho50 / rho75:     Mean density corresponding to angle theta50 / theta 75'+'\n')
        f.write('Mass fraction:     Array with cumulative mass fraction for increasing angle with the orbital plane.'+'\n')
        f.write('                   Array corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('meanRhoSm:         Smoothened array of mean rho. Array corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('x:                 Array going from 0 to 1/2. Corresponds to theta: pi/2(orbital plane) to pi&0(polar axes))'+'\n')
        f.write('Apastron side / Periastron side / Full :  gives which part of the data is used: x<0 / x>0 / all.'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('Values:'+'\n')
        f.write('\n')
        f.write('The ratio of the mean density in the orbital plane to the mean density on the polar axis is :'+ str(round(infoForPlot['meanRhoSm'][0]/infoForPlot['meanRhoSm'][-1],3))+'\n')
        f.write('This is a usefull ratio to measure the EDE!'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('theta(50),                   full = '+ str(np.round(infoForPlot['theta50']/np.pi ,3))+'pi'+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa = '+ str(np.round(infoForPlotL['theta50']/np.pi ,3))+'pi'+'\n')
            f.write('                              Per = '+ str(np.round(infoForPlotR['theta50']/np.pi ,3))+'pi'+'\n')

        f.write('theta(50) to plot on x-axis, full:  '+str(np.round(theta['50x'],3))+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa:  '+str(np.round(theta['50xL'],3))+'\n')
            f.write('                              Per:  '+str(np.round(theta['50xR'],3))+'\n')


        f.write('rho(50) to plot on y-axis, Full:    '+str(infoForPlot['rho50'])+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa:  '+str(infoForPlotL['rho50'])+'\n')
            f.write('                              Per:  '+str(infoForPlotR['rho50'])+'\n')    
        f.write('\n')
        f.write('theta(75),                    full = '+str(np.round(infoForPlot['theta75']/np.pi,3))+'pi'+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa = '+ str(np.round(infoForPlotL['theta75']/np.pi ,3))+'pi'+'\n')
            f.write('                              Per = '+ str(np.round(infoForPlotR['theta75']/np.pi ,3))+'pi'+'\n')
        f.write('theta(75) to plot on x-axis, full:  '+str(np.round(theta['75x'],3))+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa:  '+str(np.round(theta['75xL'],3))+'\n')
            f.write('                              Per:  '+str(np.round(theta['75xR'],3))+'\n')
        f.write('rho(75) to plot on y-axis, Full:    '+str(infoForPlot['rho75'])+'\n')
        if setup['single_star'] == False:
            f.write('                              Apa:  '+str(infoForPlotL['rho75'])+'\n')
            f.write('                              Per:  '+str(infoForPlotR['rho75'])+'\n')    
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
            f.write("{: <34} {: <34} {: <34} ".format(*names))
            f.write('\n')
            for i in zip(infoForPlot['x'], infoForPlot['massFraction'][:-1], infoForPlot['meanRhoSm']):
                f.write(col_format.format(*i))
        




