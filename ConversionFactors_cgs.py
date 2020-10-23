
# --- Physical constants in cgs ---


# gravitational constant 
def G():
    G = 6.674e-8                    # [cm^3 g^-1 s^-1]
    return G

# mass of a hydrogen
def mH():
    mH = 1.6735337254999998e-24     # [g]
    return mH

# Boltzmann constant
def kB():
    kB = 1.380649e-16               # [erg K^-1] 
    return kB


# --- Conversion factors ---

# length: astronomical unit [AU] to cm
def AU_cm():
    au = 1.496e+13                  # [cm]
    return au

# mass: solar mass to gram
def Msun_gram():
    Msun = 1.98847e+33              # [g]
    return Msun

# time: seconds to years
def sec_year():
    year = 1/(60*60*24*365)         # [year]
    return year

# speed: cm/s to km/s
def cms_kms():
    kms = 1e-5                      # [cm s^-1]
    return kms


# --- Conversion factors specific to PHANTOM ---
'''
The SPH code PHANTOM calculates the different quantities in code units,
that need to be converted to cgs units. The different conversions are
given below.
( "cu" in the functions stands for "code unit" )
'''

# Converts the time unit to years
def cu_time():
    timeUnit        = 1.5916423E-01     # [yrs]
    return timeUnit

# Converts the velocity in code units to cm s^-1
def cu_vel():
    factor_velocity = 2.9784608E+06     # [cm s^-1]
    return factor_velocity
    
# Converts the specific energy in code units to erg/g
def cu_e():
    factor_energy   = 8.8712277E+12     # [erg/g]
    return factor_energy

# Converts density in code units to g cm^-3
def cu_dens():
    factor_density  = 5.9410314E-07     # [g cm^-3]
    return factor_density



