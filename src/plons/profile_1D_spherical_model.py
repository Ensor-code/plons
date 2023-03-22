#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 09:26:15 2021

@author: Ward Homan, Lionel Siess

Make sure that the path in the sys.append points to the scripts directory of phantom
"""


import matplotlib.pyplot as plt
from numpy import *
from scipy import constants, interpolate
import userSettings as us
import os
import sys

userSettingsFilePath = os.path.join( os.getcwd(), "userSettings.txt")
if not os.path.isfile(userSettingsFilePath) or os.stat(userSettingsFilePath).st_size == 0: us.create(userSettingsFilePath)
userSettingsDictionary = us.load(userSettingsFilePath,onlyPathToPhantom=True)
sys.path.append(userSettingsDictionary["hard_path_to_phantom"]+'/scripts')

from readPhantomDump import *
from PhysicalConstantsCGS import Rg, steboltz, kboltz, au, gg, c, mass_proton_cgs, solarm


# specifies the references of the x-axis to which the data is plotted
def reference(axis):
    if axis == 'r':
        x_ref    = r_ref
        x_label  = r'$Distance$ [$R_*]$'
        x_in     = 0.
        x_out    = 12.*x_ref
        x_limits = [x_in/x_ref,x_out/x_ref]
    if axis == 'Tgas':
        x_ref    = T_ref
        x_label  = r'$Temperature$ [$T_{inj}]$'
        x_in     = 1.1*Tinj
        x_out    = 0.1*Tinj
        x_limits = [1.1*Tinj/T_ref,0.1*Tinj/T_ref]
    return x_ref, x_label, x_in, x_out, x_limits

# reads the 1D data and computes derived quantities
def read1D(fname):
    # Read the 1D data
    print('Reading 1D data ',fname)
    f       = open(fname,'r')
    headers = f.readline()
    headers = headers.split()[0:]
    if (headers[0] == '#'):
        headers.pop(0)
    data_1D = dict()
    #print(headers)
    m = loadtxt(f)
    for i,column in enumerate(headers):
      data_1D[column] = m[:,i]
    #add new variables
    try:
      data_1D['Tgas']  = data_1D['T']
      data_1D['cs']    = sqrt(data_1D['gamma']*kboltz*data_1D['T']/(data_1D['mu']*mass_proton_cgs))
    except:
      pass
    data_1D['v_esc'] = sqrt(2*gg*Mstar/(data_1D['r']))
    if idust_opacity == 1:
        data_1D['kappa'] = kappa_max*(1 + exp((data_1D['Tgas']-Tcond)/delta))**(-1)
    if iget_tdust == 1:
        data_1D['Tdust'] = Tstar*(Rstar/data_1D['r'])**tdust_exp
    del headers, m, i, column
    return data_1D

# reads the 3D SPH data and computes derived quantities
def readSPH(fdump):
    data3D = read_dump(fdump)
    print('Reading 3D data ',fdump)
    udist, umass, utime = [data3D['units'][u] for u in ['udist', 'umass', 'utime']]
    gam   = data3D['quantities']['gamma']
    hfact = data3D['quantities']['hfact']
    mpart = data3D['quantities']['massoftype'][0]*umass
    x, y, z, vx, vy, vz, h = [array(data3D['blocks'][0]['data'][c],dtype=float)*conv for c, conv in [('x', udist), ('y', udist), ('z', udist), ('vx', udist/utime), ('vy', udist/utime), ('vz', udist/utime), ('h', udist)]]
    r   = sqrt(x**2+y**2+z**2)
    v   = sqrt(vx**2+vy**2+vz**2)
    rho = mpart * (hfact/h)**3
    try:
      u   = data3D['blocks'][0]['data']['u']*(udist/utime)**2
      if idust_opacity == 1:
        T = mu*mass_proton_cgs/kboltz*(gamma-1.)*u
        data3D['blocks'][0]['data']['Tgas']  = T
        data3D['blocks'][0]['data']['kappa'] = kappa_max*(1 + exp((T-Tcond)/delta))**(-1)
      elif idust_opacity == 2:
        T = data3D['blocks'][0]['data']['mu']*mass_proton_cgs/kboltz*(data3D['blocks'][0]['data']['gamma']-1.)*u
        data3D['blocks'][0]['data']['Tgas']  = T
        data3D['blocks'][0]['data']['cs']    = sqrt(gam*kboltz*T/(data3D['blocks'][0]['data']['mu']*mass_proton_cgs))
      else:
        T = mu*mass_proton_cgs/kboltz*(gam-1.)*u
        data3D['blocks'][0]['data']['Tgas'] = T
        data3D['blocks'][0]['data']['cs']   = sqrt(gam*kboltz*T/(mu*mass_proton_cgs))
        data3D['blocks'][0]['data']['u']    = u
      del T, u
    except:
      pass
    data3D['blocks'][0]['data']['r']   = r
    data3D['blocks'][0]['data']['v']   = v
    data3D['blocks'][0]['data']['rho'] = rho
    del x, y, z, vx, vy, vz, h, hfact, mpart, r, v, rho
    return data3D

# Determines the position of the last boundary shell
def r_last_shell ():
    d_tang       = Rinj*2*(constants.golden*sqrt(5))**(-1/2)/(2*res-1)
    d_rad        = wss*d_tang
    r_last_shell = Rinj+N_bound_sphere*d_rad
    return r_last_shell

# Plots the 1D data
def plot1D (data, i, second=False):
    if whichPlot == 'vel':
        if second == False:
            ax1[i][0].plot(data[xAxis]/x_ref, data['v'] /v_ref,    color='black', linestyle='-',  label='v     analytic')
            ax1[i][0].plot(data[xAxis]/x_ref, data['cs']/v_ref,    color='green', linestyle='--', label='sound speed')
            if isink_radiation>0 and alpha_rad>=1:
              pass
            else:
              ax1[i][0].plot(data[xAxis]/x_ref, data['v_esc']/v_ref, color='blue',  linestyle='-',  label='escape vel')
        else:
            ax1[i][0].plot(data[xAxis]/x_ref, data['v']/v_ref,    color='black', linestyle='--')
        ax1[i][0].set_ylabel(r'$v$ [km s$^{-1}$]')

    if whichPlot == 'temp':
        ax1[i][0].plot(data[xAxis]/x_ref, data['T']/T_ref, color='black',  linestyle='-',  label=r'T$_{gas}$  analytic')
        ax1[i][0].set_ylabel('Temperature [K]')
        if iget_tdust > 0:
          ax1[i][0].plot(data[xAxis]/x_ref, data['Tdust']/T_ref, color='magenta',  linestyle='-',  label=r'T$_{dust}$  analytic')
          
    if whichPlot == 'dustcool':
        lns11 = ax1[i][0].plot(data[xAxis]/x_ref, data['T']/T_ref, color='black',  linestyle='-',  label=r'T$_{gas}$  analytic')
        ax1[i][0].set_ylabel('Temperature [K]')
        lns12,lns2 = Dustcooling(data)
        lns1 = lns11 + lns12        
        return lns1, lns2
        
    if whichPlot == 'v&T':
       lns21 = ax2[i][0].plot(data[xAxis]/x_ref, data['T']/T_ref, color='grey',  linestyle='-',  label=r'T$_{gas}$  analytic')
       if iget_tdust > 0:
          lns23 = ax2[i][0].plot(data[xAxis]/x_ref, data['Tdust']/T_ref, color='magenta',  linestyle='-',  label=r'T$_{dust}$  analytic')
          lns21 = lns21 + lns23
       lns22 = ax1[i][0].plot(data[xAxis]/x_ref, data['cs']/v_ref,    color='green', linestyle='--', label='sound speed')
       lns1  = ax1[i][0].plot(data[xAxis]/x_ref, data['v'] /v_ref,    color='black', linestyle='-',  label='v     analytic')
       lns2  = lns21 + lns22
       ax1[i][0].set_ylabel(r'$v$ [km s$^{-1}$]')
       ax2[i][0].set_ylabel('Temperature [K]')
       return lns1, lns2
     
    if whichPlot == 'dust':
        if iget_tdust > 0:
            lns1 = ax1[i][0].plot(data[xAxis]/x_ref, data['Tdust']/T_ref, color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1[i][0].set_ylabel('Temperature [K]')
        else:
            lns1 = ax1[i][0].plot(data[xAxis]/x_ref, data['T']/T_ref,     color='red',  linestyle='-', label=r'T$_{dust}$=T$_{gas}$ analytic')
            ax1[i][0].set_ylabel('Temperature [K]')
        if idust_opacity > 0:
            lns21 = ax2[i][0].plot(data[xAxis]/x_ref, data['kappa'], color='green',  linestyle='-', label=r'$\kappa_\mathrm{d}$ analytic')
            K3    = multiply(data['K3'],data['rho'])/mass_proton_cgs
            K3    = array([1e-95 if k == 0. else k for k in K3])
            lns22 = ax2[i][0].plot(data[xAxis]/x_ref, log10(K3), color='blue',  linestyle='-', label=r'log($\mathcal{K}_3$) analytic')
            lns2  = lns21 + lns22
            ax2[i][0].set_ylabel(r' $\kappa_\mathrm{d}$ [cm$^2$/g]  OR  log($\mathcal{K}_3$) [cm$^{-3}$]')
            return lns1, lns2
        else:
            return lns1

    if whichPlot == 'chem':
        if iget_tdust > 0:
            lns1 = ax1[i][0].plot(data[xAxis]/x_ref, data['Tdust']/T_ref, color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1[i][0].set_ylabel('Temperature [K]')
        else:
            lns1 = ax1[i][0].plot(data[xAxis]/x_ref, data['T']/T_ref,     color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1[i][0].set_ylabel('Temperature [K]')
        if idust_opacity == 2:
            lns21 = ax2[i][0].plot(data[xAxis]/x_ref, data['gamma'],   color='green',   linestyle='-', label=r'$\gamma$ analytic')
            lns22 = ax2[i][0].plot(data[xAxis]/x_ref, data['mu'],      color='blue',    linestyle='-', label=r'$\mu$ analytic')
            lns23 = ax2[i][0].plot(data[xAxis]/x_ref, log10(data['S']), color='magenta', linestyle='-', label=r'log($S$) analytic')
            lns2  = lns21 + lns22 + lns23
            ax2[i][0].set_ylabel(r' $\gamma$  OR  $\mu$  OR  log($S$)')
            return lns1, lns2

    if whichPlot == 'tau':
        ax1[i][0].plot(data[xAxis]/x_ref, data['tau'], color='black',  linestyle='-',  label=r'$\tau$  analytic')
        ax1[i][0].set_ylabel('optical depth')

    if whichPlot == 'tau_lucy':
        ax1[i][0].plot(data[xAxis]/x_ref, data['tau_lucy'], color='black',  linestyle='-',  label=r'$\tau_{Lucy}$  analytic')
        ax1[i][0].set_ylabel('Lucy optical depth')
  
# Plots the 3D SPH data
def plot3D (data, i, second=False):
    if whichPlot == 'vel':
        if second == False:
           ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['v']/v_ref, 'r.', label='v     SPH')
        else:
          ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['v']/v_ref, 'r.')
          
    if whichPlot == 'temp':
        ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tgas']/T_ref, 'r.', label=r'T$_{gas}$  SPH')
        if iget_tdust > 0:
          ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tdust']/T_ref, 'm.', label=r'T$_{dust}$  SPH')
        
    if whichPlot == 'v&T':
      lns3 = ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['v']/v_ref, 'r.', label='v     SPH')
      ax2[i][0]  = ax1[i][0].twinx()
      lns4 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tgas']/T_ref, 'b.', label=r'T$_{gas}$  SPH')
      if iget_tdust > 0:
        lns41 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tdust']/T_ref, 'm.', label=r'T$_{dust}$  SPH')
        lns4 = lns4 + lns41
      return lns3, lns4
    
    if whichPlot == 'dust':
        if iget_tdust > 0:
            lns3 = ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tdust']/T_ref, 'r.', label=r'T$_{dust}$  SPH')
        else:
            lns3 = ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tgas']/T_ref,  'r.', label=r'T$_{dust}$=T$_{gas}$  SPH')
        if idust_opacity >= 1:
            ax2[i][0]=ax1[i][0].twinx()
            lns41 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['kappa'],  'g.', label=r'$\kappa_\mathrm{d}$  SPH')
            K3    = multiply(data['blocks'][0]['data']['K3'],data['blocks'][0]['data']['rho'])/mass_proton_cgs
            K3    = array([1e-95 if k == 0. else k for k in K3])
            lns42 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, log10(K3), 'b.', label=r'log($\mathcal{K}_3$)  SPH')
            lns4  = lns41 + lns42
            return lns3, lns4
        else:
            return lns3

    if whichPlot == 'chem':
        if iget_tdust > 0:
            lns3 = ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tdust']/T_ref, 'r.', label=r'T$_{dust}$  SPH')
        else:
            lns3 = ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['Tgas']/T_ref,  'r.', label=r'T$_{dust}$  SPH')
        if idust_opacity == 2:            
            ax2[i][0]   = ax1[i][0].twinx()
            lns41 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['gamma'],   'g.', label=r'$\gamma$ SPH')
            lns42 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['mu'],      'b.', label=r'$\mu$ SPH')
            lns43 = ax2[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, log10(data['blocks'][0]['data']['S']), 'm.', label=r'log($S$) SPH')
            lns4  = lns41 + lns42 + lns43
            return lns3, lns4
        else:
            print('\n=== WARNING ===')
            print('Nucleation not activated, no chemistry data available')
            print('Please use different plotting option')
            print('Aborting script')
            print('===============\n')
            sys.exit()

    if whichPlot == 'tau':
        ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['tau'], 'r.', label=r'$\tau$  SPH')
        ax1[i][0].set_ylim([0,1.2*max(data['blocks'][0]['data']['tau'])])

    if whichPlot == 'tau_lucy':
        ax1[i][0].plot(data['blocks'][0]['data'][xAxis]/x_ref, data['blocks'][0]['data']['tau_lucy'], 'r.', label=r'$\tau_{Lucy}$  SPH')
        ax1[i][0].set_ylim([0,1.2*max(data['blocks'][0]['data']['tau_lucy'])])


# Plots auxiliary quantities that are not strictly related to the 1D or 3D SPH data
def plot_Auxiliary (data,i):
    if xAxis == 'r':
        lns5 = ax1[i][0].axvline(x=r_last_shell()/Rstar,ymin=0,ymax=3,color='orange',linestyle='-',label='last boundary shell')
        return lns5
    else:
      pass

# Computes and prints the L2 norm of the SPH data with reference to the 1D profile
def L2norm(data1D,data3D,xaxis):
  
  quantity = 'v'
  x        = data1D[xaxis]
  x_sph    = data3D['blocks'][0]['data'][xaxis]
  y        = data1D[quantity]
  y_sph    = data3D['blocks'][0]['data'][quantity]
  
  # Arrange sph points according to bins of 1D solution
  bins          = x
  binsIndices   = digitize(x_sph,bins)
  arranged_data = [ [] for i in range(len(bins))]
  for i in range(len(x_sph)):
    arranged_data[binsIndices[i]].append(y_sph[i])
  
  # Calculate L2 norm for each bin
  L2      = zeros(len(binsIndices))
  for i in range(len(bins)):
    L2j   = 0
    for j in range(len(arranged_data[i])):
      L2j = L2j+(arranged_data[i][j]-y[i])**2 # Sum square of distance
    L2[i] = sqrt(L2j)/1e5  # Calculate square-root and convert to km/s
    
  L2 = array(L2)   # convert list to array
  L2 = L2[L2 != 0] # Remove empty array elements
  L2 = L2[1:-1]    # Discard first and last element, where binning algorithm fails
  print('===============')
  print('Mean over "'+str(xaxis)+'" of L2 norm of "'+str(quantity)+'" = ' +str(round(mean(L2[1:-1]),2))+' km/s')
  print('===============')  

# Computes and prints the standard deviation of the SPH data with reference to the 1D profile  
def GaussianSigma(data1D,data3D,xaxis):
  
  quantity       = 'v'
  quantity_scale = 1e5
  x        = data1D[xaxis]
  x_sph    = data3D['blocks'][0]['data'][xaxis]
  y        = data1D[quantity]/quantity_scale
  y_sph    = data3D['blocks'][0]['data'][quantity]/quantity_scale
  
  interp       = interpolate.interp1d(x, y, kind='cubic', bounds_error=False, fill_value='extrapolate')
  y_sph_interp = interp(x_sph)
  
  sorted_indices      = argsort(x_sph)
  x_sph_sorted        = x_sph[sorted_indices]
  y_sph_sorted        = y_sph[sorted_indices]
  y_sph_interp_sorted = y_sph_interp[sorted_indices]
  
  outer_boundary = argmax(x_sph_sorted > x_out)  
  
  MSE = sum(square(subtract(y_sph_sorted[:outer_boundary],y_sph_interp_sorted[:outer_boundary])))/len(y_sph_sorted[:outer_boundary])
  STD = sqrt(MSE)
  print('===============')
  print('Standard deviation of '+str(quantity)+'_sph w.r.t. '+str(quantity)+'_1D = ' +str(round(STD,2))+' km/s')
  print('===============')  


# Work in progress
def Dustcooling(data1D):
  
  Tprofile = loadtxt('temperatureprofile_AC_0.1um.txt')
  Tprofile = [*zip(*Tprofile)]
  
  radii        = array(Tprofile[0])
  Tdust        = array(Tprofile[1])  
  interp       = interpolate.interp1d(radii, Tdust, kind='cubic', bounds_error=False, fill_value='extrapolate')
  x_1D_profile = data1D[xAxis]/au
  Tdust_interp = interp(x_1D_profile)
  
  Tgas         = data1D['T']
  
  kappa_gas    = 2e-4
  bowen_Cprime = 3.000e-5
  d2g          = 1/100
  rho_sp       = 2
  dust_size    = 1e-5
  heat_transfer_zone = 1e-7
  therm_conduc = 50
  
  CoolingBowen      = (3*Rg/(2*bowen_Cprime))*multiply( divide(data1D['rho'],data1D['mu']) , Tdust_interp-Tgas )
  CoolingHollenbach = multiply(multiply((d2g*kboltz/(2*rho_sp*dust_size))*sqrt(3*kboltz*divide(Tgas,(data1D['mu']*mass_proton_cgs)**3)),data1D['rho']**2),Tdust_interp-Tgas)
  CoolingFourier    = (3*therm_conduc*d2g/(rho_sp*dust_size*heat_transfer_zone))*multiply(data1D['rho'],Tdust_interp-Tgas)
  CoolingStefan     = 4.*steboltz*(Tdust_interp**4-Tgas**4)*kappa_gas
  
  print(max(CoolingBowen))
  print(max(CoolingHollenbach))
  print(max(CoolingFourier))
  print(max(CoolingStefan))
  
  lns12 = ax1[i][0].plot(x_1D_profile*x_ref/au, Tdust_interp, color='blue',  linestyle='-',  label=r'T$_{dust}$ MCFost')
  ax2[i][0]  = ax1[i][0].twinx()
  lns21 = ax2[i][0].plot(x_1D_profile*x_ref/au,CoolingBowen, color='red',  linestyle='--',  label='Cooling Bowen')
  lns22 = ax2[i][0].plot(x_1D_profile*x_ref/au,CoolingFourier, color='brown',  linestyle='--',  label='Cooling Fourier')
  lns23 = ax2[i][0].plot(x_1D_profile*x_ref/au,CoolingHollenbach, color='green',  linestyle='--',  label='Cooling Hollenbach')
  lns24 = ax2[i][0].plot(x_1D_profile*x_ref/au,CoolingStefan, color='purple',  linestyle='--',  label='Cooling StefBoltz')
  lns2  = lns21 + lns22 + lns23 + lns24
  ax1[i][0].axvline(x=1.8, color='orange')
  ax2[i][0].set_ylabel(r'cooling rate $\Lambda$ [erg s$^{-1}$ cm$^{-3}$]')
  
  return lns12,lns2

# ===== GENERAL INPUT PARAMETERS ========================================================================

# Main directory of the data
# mainPath   = ['/STER/matse/PHANTOM_Models/Lucy/Low/single/', '/STER/matse/PHANTOM_Models/Lucy/Low/single_Lucy/', '/STER/matse/PHANTOM_Models/Lucy/Low/single_Atten/']
# mainPath   = ['/STER/matse/PHANTOM_Models/Lucy/High/single/', '/STER/matse/PHANTOM_Models/Lucy/High/single_Lucy/', '/STER/matse/PHANTOM_Models/Lucy/High/single_Atten/']
mainPath   = ['/STER/matse/PHANTOM_Models/Lucy/High/singleLucy/']
modelLabel = 'wind'
dumpNumber = '00100'

# What do you want to plot?
whichPlot = 'tau_lucy'       # vel, temp, v&T, dust, chem, dustcool, tau, tau_lucy
xAxis     = 'r'          # r, Tgas

#Number of vertical subplots, change second dimension of ax array for horizontal plots
Nsub     = len(mainPath)

# ===== SET UP PLOT LAYOUT =============================================================

size     = (9,6)
fig, ax1 = plt.subplots(1,figsize=size,tight_layout=True,squeeze=False,sharex=True)
ax2      = copy(ax1)
i = 0
for j in range(Nsub):

# ===== READ INPUT ========================================================================

  inputFile = mainPath[j] + modelLabel
  # Extract data from parameter card
  wind_param      = read_infile(inputFile)
  Mstar           = wind_param['primary_mass']*solarm
  Rstar           = wind_param['primary_racc']*au
  Tstar           = wind_param['primary_Teff']
  Rinj            = wind_param['wind_inject_radius']*au
  mu              = wind_param['mu']
  idust_opacity   = wind_param['idust_opacity']
  isink_radiation = wind_param['isink_radiation']
  N_bound_sphere  = wind_param['iboundary_spheres']
  res             = wind_param['iwind_resolution']
  wss             = wind_param['wind_shell_spacing']
  try:
      alpha_rad  = wind_param['alpha_rad']
  except:
      alpha_rad  = 0.
  try:
      iget_tdust = wind_param['iget_tdust']
  except:
      iget_tdust = 0.
  try:
      Tinj       = wind_param['wind_temperature']
  except:
      Tinj       = 0.
  
  if idust_opacity == 1:
      kappa_gas   = wind_param['kappa_gas']
      kappa_max   = wind_param['bowen_kmax']
      Tcond       = wind_param['bowen_Tcond']
      delta       = wind_param['bowen_delta']
      mu          = wind_param['mu']
      gamma       = wind_param['wind_gamma']
  if iget_tdust == 1:
      tdust_exp   = wind_param['tdust_exp']

  # specify reference values
  r_ref  = au
  T_ref  = 1.
  v_ref  = 1.e5
  x_ref, x_label, x_in, x_out, x_limits = reference(xAxis)


  # ===== READ DATA ========================================================================

  # Read input and dump file
  data1D  = read1D(mainPath[j] + modelLabel + '_1D.dat')
  #data1D_2  = read1D(mainPath[i] + modelLabel + '_noR_1D.dat')
  dataSPH = readSPH(mainPath[j] + modelLabel + '_' + dumpNumber)
  #dataSPH_2 = readSPH(mainPath[i] + modelLabel + '_noR_' + dumpNumber)

  # ===== GENERATE PLOT ========================================================================

  if whichPlot == 'dust' or whichPlot == 'chem' or whichPlot == 'v&T' or whichPlot == 'dustcool':
      if idust_opacity > 0:
          #lns3,lns4 = plot3D(dataSPH,i)
          lns1,lns2 = plot1D(data1D,i)
      else:
          #lns3      = plot3D(dataSPH,i)
          lns1      = plot1D(data1D,i)
      #lns5 = plot_Auxiliary(data1D,i)
  else:
      plot3D(dataSPH,i)
      #plot3D(dataSPH_2,i,second=True)
      plot1D(data1D,i)
      #plot1D(data1D_2,i,second=True)
      plot_Auxiliary(data1D,i)
  
  #L2norm(data1D,dataSPH,xAxis)
  #GaussianSigma(data1D,dataSPH,xAxis)
  
  #---------- Legend management ----------
  if whichPlot == 'dust' or whichPlot == 'chem' or whichPlot == 'v&T'or whichPlot == 'dustcool':
    if idust_opacity > 0:
        lns = lns1+lns2#+lns3+lns4#+[lns5]
    else:
        lns = lns1#+lns3+[lns5]
    labs = [l.get_label() for l in lns]
    ax1[i][0].legend(lns, labs, bbox_to_anchor=(1, 0.30), loc='lower right',fontsize=13)
  else:
    ax1[i][0].legend(           bbox_to_anchor=(1, 0.2), loc='lower right', fontsize=13)

  #---------- only set xlabels on bottom panel ----------    
  if i == range(Nsub)[-1]:
    ax1[i][0].set_xlabel(x_label)
    
  #---------- Font size ----------
  for item in (  [ax1[i][0].title, ax1[i][0].xaxis.label, ax1[i][0].yaxis.label]
              + ax1[i][0].get_xticklabels() + ax1[i][0].get_yticklabels()
              + [ax2[i][0].title, ax2[i][0].xaxis.label, ax2[i][0].yaxis.label]
              + ax2[i][0].get_xticklabels() + ax2[i][0].get_yticklabels()
              ):
      item.set_fontsize(16)
  
  #---------- Finalize plot ----------
  #txt = ['high resolution, q=20','mid resolution, q=10','low resolution, q=3']
  #ax1[i][0].text(2.6,19,txt[i],fontsize=16)
  ax1[i][0].set_xlim(x_limits)
  #ax2[i][0].set_ylim([0,40])
  #ax2[i][0].set_yscale('log')
  plt.title('BowenDust model',fontsize=16)   


plt.show()
