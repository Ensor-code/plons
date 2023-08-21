import scipy.integrate       as integrate
import matplotlib.pyplot     as plt
import numpy as np
import os
from plons.ConversionFactors_cgs   import au, Msun, G

def velocities(run):
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
    # elif run == 'Lucy2/High/binary6':
    #     v_i = 22e5   # km/s
    #     v_o = 26e5   # km/s 
    else:
        return False
    return v_i, v_o

def ArchimedianSpiral(run, saveLoc, setup):
    thetaIni = np.pi
    
    if velocities(run):
        xi, yi, theta, xo, yo = ArchSpiral(run, setup, thetaIni)
    else:
        # print(run)
        return

    a_AGB  = setup['massAGB_ini']/(setup['massAGB_ini']+setup['massComp_ini'])*setup['sma_ini']*au
    e      = setup['ecc']

    # Plot location of AGB star
    rAGB = a_AGB*(1.0-e**2)/(1.0+e*np.cos(theta[0]))
    xAGB = rAGB *np.cos(theta[0])/au
    yAGB = -rAGB *np.sin(theta[0])/au
    
    fig = plt.figure(figsize=(10,10))
    
    plt.plot(xi, yi, 'k', linestyle = 'dotted',label = 'BSE',linewidth = 1.4)
    plt.plot(xo, yo, 'k-',label = 'FSE',linewidth = 1.4)
    plt.plot(xAGB, yAGB, 'bo', label = 'AGB')
    plt.plot(xo[0], yo[0], 'ro', label = 'companion')
    
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

def ArchSpiral(run, setup, thetaIni = np.pi):
    M_AGB  = setup['massAGB_ini']
    M_comp = setup['massComp_ini']
    a      = setup['sma_ini']*au
    e      = setup['ecc']

    Msum = M_AGB+M_comp
    a_comp = setup['massAGB_ini']/(setup['massAGB_ini']+setup['massComp_ini'])*a
    P = 2.0*np.pi*np.sqrt(a**3/(Msum*Msun*G))

    theta    = np.linspace(thetaIni, 14*np.pi, 20000)
    thetaBSE = np.linspace(thetaIni, 14*np.pi, 20000)
    
    rcomp0 = a_comp*(1.0-e**2)/(1.0+e*np.cos(theta[0]))

    if velocities(run):
        v_i, v_o = velocities(run)
    
    rspiralOuter = integrate.odeint(rprime, rcomp0, theta,    args=(P, a_comp, e, v_o,))
    rspiralInner = integrate.odeint(rprime, rcomp0, thetaBSE, args=(P, a_comp, e, v_i,))
    #function needs: (derivative function(y,t), y0, t points at which we want solution, other arguments of function) 
    
    ro = rspiralOuter[:,0]     
    xo = -ro*np.cos(theta)/au
    yo = ro*np.sin(theta)/au

    ri = rspiralInner[:,0]     
    xi = -ri*np.cos(thetaBSE)/au
    yi = ri*np.sin(thetaBSE)/au

    return xi, yi, theta, xo, yo