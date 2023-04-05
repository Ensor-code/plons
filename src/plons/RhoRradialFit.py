import numpy                    as np
import matplotlib.pyplot        as plt
import os

# import plons scripts
import plons.ConversionFactors_cgs    as cgs
import plons.PhysicalQuantities       as pq

# import certain things from packages
from scipy.optimize             import curve_fit
from matplotlib                 import rcParams
# Change the matplotlib default parameters
rcParams.update({'font.size':   11})
rcParams.update({'figure.dpi': 200})
#rc('font', family='serif')
#rc('text', usetex=True)

# ignore warnings
import warnings
warnings.filterwarnings("ignore")



def func1(x,scaling1):
    return np.log10(scaling1 *(x**0))

def fitRhoR(run,loc,dumpData,sinkData,setup):
    # all models names
    models = {'53':['S90slow','s', 'O',230],'59':['S40slow','s', 'O',125],'51':['S25slow','s', 'O',100],'39':['S90fast','s','C',150],'57':['S40fast','s','C',230],'49':['S25fast','s','C',125],
              '52':['P90slow','p', 'O',125],'56':['P40slow','p', 'O',50],'48':['P25slow','p', 'O',35],'41':['P90fast','p','C',170],'58':['P40fast','p','C',120],'50':['P25fast','p','C',100],
            }
    #rComp = dumpData['r'][-2]
    rComp = np.mean(sinkData['rComp']) 
    #print('rcomp [au] = ', rComp/cgs.au)
    
    dataToUse = {}
    dataToUse['rho1'] = dumpData['rho']    [dumpData['r'] < rComp]
    dataToUse['r1'  ] = dumpData['r'  ]    [dumpData['r'] < rComp]
    
    dens = dumpData['rho']    [(dumpData['r'] >= rComp)]
    rad  = dumpData['r'  ]    [(dumpData['r'] >= rComp)]
    
    dataToUse['rho2'] = dens   [(rad < models[run][-1]*cgs.au)]
    dataToUse['r2'  ] = rad    [(rad < models[run][-1]*cgs.au)]
    
    r1 = dataToUse['r1']/cgs.au
    r2 = dataToUse['r2']/cgs.au
    rC = rComp/cgs.au
    rho1 = np.log10(dataToUse['rho1'])
    rho2 = np.log10(dataToUse['rho2'])

        
    fig = plt.figure(figsize=(7, 6))  
    ax1 = plt.subplot(111)
    plt.title(models[run][0])
    
    
    popt1, pcov1 = curve_fit(func1, r1, rho1,p0=[1e-13])
    [scaling1] = popt1
    #print('scaling1 = ',*popt1) 
    
        
    def func2(x,p):
        return  np.log10(scaling1 * (rC/x)**(p))
    def func3(x,p,sc):
        return  np.log10(sc * (rC/x)**(p))

    
    if models[run][1] == 's':
        popt2, pcov2 = curve_fit(func2, r2, rho2,p0=[2])
        [power2] = popt2
    elif models[run][1] == 'p':
        popt2, pcov2 = curve_fit(func3, r2, rho2,p0=[scaling1,2])
        [power2,scaling2] = popt2
    #print('power2 = ',power2) #,' scaling2 =',scaling2)                                       

                                        
    plt.axvline(x = rC,c='k',lw =0.8)
    ax1.scatter(r1,10**rho1, s= 0.1, c='lightgrey')
    ax1.scatter(r2,10**rho2, s= 0.1, c='darkgrey')
    
    x1 = np.linspace(min(r1),max(r1),100)
    x2 = np.linspace(min(r2),max(r2),100)
    
    
    #fit 
    ax1.plot(x1,10**(func1(x1,*popt1)),lw = 0.8 ,c='royalblue', label = 'y$_1$ = '+str(np.round(scaling1,18)) )
    if models[run][1] == 's':
        ax1.plot(x2,10**(func2(x2,*popt2)),lw = 0.8 ,c='firebrick', label = 'y$_2$ = '+str(np.round(scaling1,18))+' * (r$_{\\rm comp}$/x) ^('+str(np.round(power2,3))+')')# + '+str(np.round(scaling2,17))) 
    elif models[run][1] == 'p':
        ax1.plot(x2,10**(func3(x2,*popt2)),lw = 0.8 ,c='firebrick', label = 'y$_2$ = '+str(np.round(scaling2,18))+' * (r$_{\\rm comp}$/x) ^('+str(np.round(power2,3))+')')# + '+str(np.round(scaling2,17))) 
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [au]')
    ax1.set_ylabel('density [g/cm$^3$]')
    plt.legend()
    fig.tight_layout()
    #plt.show()
    #fig.savefig(loc+str(run)+'_fitRadialRho.png')
    fig.savefig('/lhome/silkem/Documents/Pipeline/testOutput/fitRadialRho/'+str(run)+'_fitRadialRho_logFit.png')

    
    
    input_dens = pq.getInputDensity(setup['v_ini'], 2*setup['Mdot'])
    # Make text file with info 
    title = '/lhome/silkem/Documents/Pipeline/testOutput/fitRadialRho/results_fit_density.txt'
    with open (title,'a') as f:
        if models[run][1] == 's':
            f.write('\n'+str(models[run][0])+'      '+str(np.round(power2,3))+'      '+str(np.round(scaling1,19))+'     '+'////////////'+'     '+str(setup['v_ini'])+'     '+str(2*setup['Mdot'])+'      '+str(round(input_dens,17)))
        elif models[run][1] == 'p':
            f.write('\n'+str(models[run][0])+'      '+str(np.round(power2,3))+'      '+str(np.round(scaling1,19))+'     '+str(np.round(scaling2,19))+'     '+str(setup['v_ini'])+'     '+str(2*setup['Mdot'])+'      '+str(round(input_dens,17)))
        
    print('     Fit radial density '+str(run)+' ready and saved!')
'''

def func(x,sc1,p):
    if x<3.6:
        return sc1 * x**0
    elif x>3.6:
        return sc1 * 3.6**p * x**(-p)



def fitRhoR(run,loc,dumpData,sinkData):
    #rComp = dumpData['r'][-2]
    rComp = np.mean(sinkData['rComp'])
    print('rcomp= ', rComp/cgs.au)
    
    dataToUse = {}
    dataToUse['rho'     ] = dumpData['rho'     ][:-2]  
    dataToUse['r'       ] = dumpData['r'       ][:-2]      
 
    r = dataToUse['r']/cgs.au
        
    fig = plt.figure(figsize=(9, 6))  
    ax1 = plt.subplot(111)
    
    popt, pcov = curve_fit(func, r, dataToUse['rho'])
    [scaling1, power] = popt
    print('scaling1 = ',scaling1,' power = ',power) 

                                      
    plt.axvline(x = rComp/cgs.au,c='k')
    ax1.scatter(r,np.log10(dataToUse['rho']), s= 0.3, c='firebrick')
    
    x1 = np.linspace(min(r),max(r),10000)
    
    #fit 
    ax1.scatter(x1,np.log10(func1(x1,*popt)),s=0.5,c='k') #, label = 'y1 = '+str(np.round(scaling1,17)))

    ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_xlabel('r [au]')
    ax1.set_ylabel('log density [cm/g$^3$]')
    plt.legend()
    #fig.savefig(loc+str(run)+'_fitRadialRho.png')
    fig.savefig('/lhome/jolienm/Documents/phantomPipeline/fitRadialRho/'+str(run)+'_fitRadialRho_oneF.png')

'''
    
