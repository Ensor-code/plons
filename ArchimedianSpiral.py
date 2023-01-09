import scipy.integrate       as integrate
import matplotlib.pyplot     as plt
import numpy as np
import os
from ConversionFactors_cgs   import au, Msun, G


def ArchimedianSpiral(run, saveLoc, setup):
    # paramters of the binary
    M_AGB  = setup['massAGB_ini']*Msun
    M_comp = setup['massComp_ini']*Msun
    a      = setup['sma_ini']*au
    e      = setup['ecc']

    Msum = M_AGB+M_comp
    a_AGB  = M_comp/Msum*a
    a_comp = M_AGB/Msum*a
    P = 2.0*np.pi*np.sqrt(a**3/(Msum*G))

    thetaIni = np.pi
    name = 'apa'

    if   run == 'Lucy/High/binary6':
        v_i = 45e5   # km/s
        v_o = 50e5   # km/s 
    elif run == 'Lucy/High/binary9':
        v_i = 45e5   # km/s
        v_o = 50e5   # km/s 
    elif run == 'Lucy/High/binary9Atten':
        v_i = 13e5   # km/s
        v_o = 18e5   # km/s 
    elif run == 'Lucy/High/binary9Lucy':
        v_i = 16e5   # km/s
        v_o = 22e5   # km/s 
    else:
        print(run)
        v_i = 2e5    # km/s
        v_o = 15e5   # km/s 
    
    fig = plt.figure(figsize=(10,10))
    ##
    #
    # wake spiral of the companion
    theta    = np.linspace(thetaIni, 14*np.pi, 20000)
    thetaBSE = np.linspace(thetaIni, 14*np.pi, 20000)
    
    rcomp0 = a_comp*(1.0-e**2)/(1.0+e*np.cos(theta[0]))
    
    rspiralOuter = integrate.odeint(rprime, rcomp0, theta,    args=(P, a_comp, e, v_o,))
    rspiralInner = integrate.odeint(rprime, rcomp0, thetaBSE, args=(P, a_comp, e, v_i,))
    #function needs: (derivative function(y,t), y0, t points at which we want solution, other arguments of function) 
    
    ro = rspiralOuter[:,0]     
    xo = -ro*np.cos(theta)/au  
    yo = ro*np.sin(theta)/au   

    ri = rspiralInner[:,0]     
    xi = -ri*np.cos(thetaBSE)/au  
    yi = ri*np.sin(thetaBSE)/au   
    
    xcomp = xo[0]		
    ycomp = yo[0]		
    
    # Plot location of AGB star
    rAGB = a_AGB*(1.0-e**2)/(1.0+e*np.cos(theta[0]))
    xAGB = rAGB *np.cos(theta[0])/au
    yAGB = -rAGB *np.sin(theta[0])/au
    
    plt.plot(xi, yi, 'k', linestyle = 'dotted',label = 'BSE',linewidth = 1.4)
    plt.plot(xo, yo, 'k-',label = 'FSE',linewidth = 1.4)
    plt.plot(xAGB, yAGB, 'bo', label = 'AGB')
    plt.plot(xcomp, ycomp, 'ro', label = 'companion')
    
    plt.axis('square')
    plt.xlabel('x[au]',fontsize = 24)
    lim = (setup['bound'] * np.sqrt(2.) / 2.)
    lim = round(lim)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.ylabel('y[au]',fontsize = 24)
    plt.legend(fontsize = 20)
    plt.tick_params(labelsize=24)
    #plt.title('e = '+str(e)+', '+str(vel), fontsize = 18)
    #plt.title('v20e00',fontsize = 19)
    fig.savefig(os.path.join(saveLoc, 'png/2Dplot_ArchimedianSpiral.png'), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(saveLoc, 'pdf/2Dplot_ArchimedianSpiral.pdf'), dpi=200, bbox_inches="tight")
        
        
def rprime(r,theta,P,a,e,v):
    omega = 2.0*np.pi*(1.0+e*np.cos(theta))**2/(P*(1.0-e**2)**(3.0/2.0))
    rprime = v/omega 
    return rprime