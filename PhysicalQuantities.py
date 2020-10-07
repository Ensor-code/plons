import numpy                    as np
import math                     as math
import os
import ConversionFactors_cgs    as cgs
import GeometricalFunctions     as gf


# --- Physical Quantities in cgs ---

# Returns the pressure of a particle, given by the specific internal energy u [erg/g], gamma and density [g/cm^3] in 10^-1 Pa = barye = Ba = g cm^-1 s^-2
def getPressure(density, u, gamma):
    return (gamma-1)*density*u

# Returns the temperature, given the pressure [Ba] and density [g/cm^3] in K (via Ideal Gas Law)
def getTemp(pressure, density, mu):
    temp = ((pressure*mu*cgs.mH())/(density*cgs.kB()))
    return temp

# Returns the speed of sound [km/s] of the local medium: c_s = sqrt(gamma*P/rho), with gamma = cst, P in [Ba] and rho in [g/cm^3]
def getSoundSpeed(pressure, density, gamma):
    soundSpeed = np.sqrt((gamma*pressure)/(density))      
    return soundSpeed                                   # cm/s

# Returns the orbital period (s) of the binary system via Kepler's second law.
#       mass1 and mass2 in gram, orbSep in AU
def getPeriod(mass1,mass2,orbSep):
    period = 2*np.pi*np.sqrt((orbSep*cgs.AU_cm())**3/(cgs.G()*(mass1+mass2)))
    return period                           # in s

# Returns the orbital velocity [cm/s] of the binary system.
#       period in seconds, orbSep in AU
def getOrbitalVelocity(period, orbSep):
    v_orb = (2*np.pi*orbSep*cgs.AU_cm())/(period)      
    return v_orb                                    # cm/s

# Get the radial and tangential velocity 
#   x, y, and phi are lists
def getRadTanVelocity(x,y,v_x,v_y):
    phi = []
    for i in range(len(y)):
        phi.append(math.atan2(y[i],x[i]))
    #v_rad = []
    v_tan = []
    for i in range(len(x)):
        #v_rad.append(np.abs(coordTransf(v_x[i],v_y[i],phi[i])[0]))
        v_tan.append(np.abs(gf.coordTransf(v_x[i],v_y[i],phi[i])[1]))
    return np.array(v_tan)

# Calculate hill sphere for a certain object:
def getRHill(orbSep, mcomp, mAGB):
    rH = orbSep * ((mcomp/(3*mAGB))**(1./3.))
    return rH


