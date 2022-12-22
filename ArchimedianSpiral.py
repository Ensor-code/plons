import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
# physical constants
au = 1.495978707e13
Msun = 1.9891e33
G = 6.67408e-8
#
# paramters of the binary
M_AGB  = 1.5*Msun
M_comp = 1.0*Msun
a = 6*au
e = 0.0

Msum = M_AGB+M_comp
a_AGB  = M_comp/Msum*a
a_comp = M_AGB/Msum*a
P = 2.0*np.pi*np.sqrt(a**3/(Msum*G))

def rprime(r,theta,a,e,v):
    omega = 2.0*np.pi*(1.0+e*np.cos(theta))**2/(P*(1.0-e**2)**(3.0/2.0))
    rprime = v/omega 
    return rprime

def plotTheorSpiral(ecc,orbPh,vel):
    if orbPh == '0':
        thetaIni = 0
        name = 'per'
    elif orbPh == '1':
        thetaIni = -2.4 #np.pi/2  
        #thetaIni = -0.7
        name = 'PerToApa'
    elif orbPh == '2':
        thetaIni = np.pi
        name = 'apa'
    elif orbPh == '3':
        thetaIni = -3.9 #3*np.pi/2
        #thetaIni = 0.7
        name = 'ApaToPer'
        
    if vel == 'Lowv':
        #v0   = 5e5    #7.7e5  # km/s
        #vinf = 8.7e5  # km/s
        v_i = 2e5    # km/s
        v_o = 15e5   # km/s 
    elif vel == 'Medv':
        #v0   = 10e5   # 11.5e5 # km/s
        #vinf = 12.5e5 # km/s  
        v_i = 5e5    # km/s
        v_o = 18e5   # km/s  
    elif vel == 'Highv':
        #v0   = 20e5   # 20.8e5 # km/s
        #vinf = 21.3e5 # km/s
        v_i = 11e5   # km/s
        v_o = 24e5   # km/s      
    elif vel =='Highv_cool':
        v_i = 11e5   # km/s
        v_o = 19.5e5 # km/s
        
    e = ecc      
       
    fig = plt.figure(figsize=(10,10))
    ##
    #
    # wake spiral of the companion
    theta = np.linspace(thetaIni, 14*np.pi, 20000)
    thetaBSE = np.linspace(thetaIni, 5.6*np.pi, 20000)
    
    rcomp0 = a_comp*(1.0-e**2)/(1.0+e*np.cos(theta[0]))
    
    rspiralOuter = integrate.odeint(rprime, rcomp0, theta, args=(a_comp, e, v_o,))
    rspiralInner = integrate.odeint(rprime, rcomp0, thetaBSE, args=(a_comp, e, v_i,))
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
    plt.xlim(-120,120)
    #plt.xlim(-50,50)
    plt.ylabel('y[au]',fontsize = 24)
    plt.ylim(-120,120)
    #plt.ylim(-50,50)
    plt.legend(fontsize = 20)
    plt.tick_params(labelsize=24)
    #plt.title('e = '+str(e)+', '+str(vel), fontsize = 18)
    #plt.title('v20e00',fontsize = 19)
    fig.tight_layout()
    #fig.savefig('/lhome/jolienm/Documents/thesisModellen/theorSpirals/'+vel+'_e_'+str(e)+'_orbPh_'+name+'.png')
    # fig.savefig('/lhome/jolienm/Documents/TierModels/R_Aql/cooling/binariesInPaper/theorSpiral/test_z120_'+vel+'_e_'+str(e)+'_orbPh_'+name+'.png')
    plt.show()
        

vel = ['Highv_cool']
ecc = [0.00]
orbPh = ['2']#'0','1','2','3']
for e in ecc:
    for i in orbPh:
        for v in vel:
            plotTheorSpiral(e,i,v)
    
    
    
    