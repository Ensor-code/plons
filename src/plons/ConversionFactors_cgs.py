
# --- Physical constants in cgs ---

G = 6.674e-8                    # gravitational constant     [cm^3 g^-1 s^-2]
mH = 1.6735337254999998e-24     # mass of a hydrogen         [g]
kB = 1.380649e-16               # Boltzmann constant         [erg K^-1] 
Rg = 8.31446261815324e7         # Gas constant               [erg/K/g]
steboltz   = 5.67051e-5         # Stefan-Boltzmann constant  [erg cm^-2K^-4 s^-1]
c = 2.99792458e10               # speed of light             [cm s^-1]


# --- Conversion factors ---

au = 1.496e+13                  # length: astronomical unit  [cm]
Msun = 1.98847e+33              # mass: solar mass           [g]
Lsun = 3.826e33                 # lumin: solar luminocity    [ergs s^-1]
year = 60*60*24*365             # time: year                 [s]
kms = 1e5                       # speed: km/s                [cm s^-1]


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
    #return timeUnit
    return timeUnit

# Converts the velocity in code units to cm s^-1
def cu_vel():
    factor_velocity = 2.9784608E+06     # [cm s^-1]
    return factor_velocity
    #return 1.

# Converts the specific energy in code units to erg/g
def cu_e():
    factor_energy   = 8.8712277E+12     # [erg/g]
    return factor_energy
    #return 1.

# Converts density in code units to g cm^-3
def cu_dens():
    factor_density  = 5.9410314E-07     # [g cm^-3]
    #return factor_density
    return 1.

def cu_J():
    factorJ = Msun* au**2 / (5.023E6)  #g cm**2 /s
    return factorJ


