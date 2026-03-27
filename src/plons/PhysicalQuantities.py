import math                     as math
import numpy                    as np
from scipy.optimize import brentq

#import plons scripts
import plons.ConversionFactors_cgs    as cgs
import plons.GeometricalFunctions     as gf


# --- Physical Quantities in cgs ---

'''
Returns the pressure of a particle, given by the specific internal energy u [erg/g], gamma and density [g/cm^3] in 10^-1 Pa = barye = Ba = g cm^-1 s^-2
'''
def getPressure(density, u, gamma):
    return (gamma-1)*density*u

'''
Returns the temperature, given gamma, mu and the internal energy u
'''
def getTemp(gamma, mu, u):
    temp = (gamma-1.) * mu * u * cgs.mH / cgs.kB
    return temp

'''
Returns the density, given smoothing length, hfact (proportionality factor specifying the smoothing length) and the particle mass [g] in g/cm^3
'''
def getRho(hi, hfact, pmassi):
    rhoh = pmassi*hfact**3/hi/hi/hi
    return rhoh

'''
Returns the opacity, given the equilibrium temperature in [K] in cm^2/g
'''
def getKappa(Teq, kappa_gas = 2e-4, bowen_delta = 60., bowen_Tcond = 1500., bowen_max   = 2.7991, lumAGB = 0., massAGB = 0.):
    if bowen_max < 0:
        bowen_max = calcKappaMax(lumAGB, massAGB)
    kappa = bowen_max/(1+np.exp((Teq-bowen_Tcond)/bowen_delta))+kappa_gas
    return kappa

def calcKappaMax(lumAGB, massAGB):
    kappa_max = 0.95*(4*np.pi*cgs.c*cgs.G*massAGB)/(lumAGB)
    return kappa_max

'''
Returns the Eddington factor, the opacity in [cm^2/g], the luminocity [erg/s] and mass [g] of the AGB star, optional optical depth
'''
def getGamma(kappa, lumAGB, massAGB, tau = 0):
    Gamma = kappa*lumAGB*np.exp(-tau)/(4*np.pi*cgs.c*cgs.G*massAGB)
    return Gamma

'''
Returns the speed of sound [km/s] of the local medium: c_s = sqrt(gamma*P/rho), with gamma = cst, P in [Ba] and rho in [g/cm^3]
'''
def getSoundSpeed(pressure, density, gamma):
    soundSpeed = np.sqrt((gamma*pressure)/(density))
    return soundSpeed                                   # cm/s

'''
Returns the orbital period (s) of the binary system via Kepler's second law.
      mass1 and mass2 in gram, orbSep in AU
'''
def getPeriod(mass1,mass2,orbSep):
    period = 2*np.pi*np.sqrt((orbSep*cgs.au)**3/(cgs.G*(mass1+mass2)))
    return period                           # in s

'''
Returns the orbital velocity [cm/s] of the binary system.
      period in seconds, orbSep in AU
'''
def getOrbitalVelocity(period, orbSep):
    v_orb = (2*np.pi*orbSep*cgs.au)/(period)      
    return v_orb                                    # cm/s

'''
Calculate hill sphere for a certain object.
'''
def getRHill(orbSep, mcomp, mAGB):
    rH = orbSep * ((mcomp/(3*mAGB))**(1./3.))
    return rH

'''
Calculate the accretion radius of the companion.
'''
def getCaptureRadius(mcomp, vwind):
    Rcap = (2*cgs.G*mcomp)/(vwind**2)
    return Rcap                         # [cm]

'''
The epsilon parameter gives an indication of morphological classification:
    >> 1 : spherically symmetric
    ~  1 : regular spiral
    << 1 : complex
It is defined as the ratio of enery densities: 
    epsilon = kin_energy / grav_energy 
            = (v_wind**2 * sma) / ((24 * G**3 * M_comp**2 * M_AGB)^(1/3))
'''
def getEpsilon(vwind, sma, mComp, mAGB):
    epsilon = (vwind**2 * sma)/(cgs.G * (24 * (mComp)**2 * mAGB)**(1/3))
    return epsilon


'''
Calculate the polar angle of the companion w.r.t. the x-axis
! Same as phi from TransformToSpherical in GeometricalFunctions (which is for arrays)
'''
def getPolarAngleCompanion(x, y):
    theta = np.abs(np.arctan(y/x))

    # Second Quadrant
    if x < 0. and y > 0.:
        theta = np.pi - theta

    # Third Quadrant
    elif x < 0. and y < 0.:
        theta = np.pi + theta

    # Fourth Quadrant
    elif x > 0. and y < 0.:
        theta = 2.*np.pi - theta

    return theta

'''
Calculates the potential at a given point
'''
def potential(X, Y, M1, M2, pos_AGB, pos_comp, omega):
        G = cgs.G
        r1 = np.sqrt((X - pos_AGB[0])**2 + (Y - pos_AGB[1])**2)
        r2 = np.sqrt((X - pos_comp[0])**2 + (Y - pos_comp[1])**2)
        return -G * M1 / r1 - G * M2 / r2 - 0.5 * omega**2 * (X**2 + Y**2)
        
'''
Calculates the Roche lobes of the primary and secondary by finding phi(L_1)
'''
def getRocheLobes(dumpData):
    
    
    M_AGB  = dumpData._params['massAGB'][0]
    M_comp = dumpData._params['massComp']
    pos_AGB  = np.array(dumpData._params['posAGB'])
    pos_comp = np.array(dumpData._params['posComp'])
    r_AGB  = np.linalg.norm(pos_AGB)
    r_comp = np.linalg.norm(pos_comp)
    a = np.linalg.norm(pos_AGB - pos_comp)
    omega = np.sqrt(cgs.G * (M_AGB + M_comp) / a**3)
    
    def dphi_dx(r):
        '''
        Calculates the derivative of the potential to find the saddle point
        '''
        
        dr1 = r - r_AGB
        dr2 = r - (-r_comp)
        
        return cgs.G * M_AGB/dr1**2 * np.sign(dr1) + cgs.G*M_comp/dr2**2 * np.sign(dr2) - omega**2 * r

    def calculate_roche_lobes():
        '''
        Solve for the saddle point 
        '''
        x_lo = -r_comp * 0.99
        x_hi =  r_AGB  * 0.99

        L1_r      = brentq(dphi_dx, x_lo, x_hi)
        unit_vec  = (pos_comp - pos_AGB) / a
        frac      = (r_AGB - L1_r) / a
        L1_pos = pos_AGB + frac * unit_vec * a
        phi_L1    = potential(L1_pos[0], L1_pos[1], M_AGB, M_comp, pos_AGB, pos_comp, omega)
        
        return L1_pos, phi_L1     
    
    return calculate_roche_lobes()
