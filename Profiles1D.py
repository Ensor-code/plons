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

# reads the 1D data and computes derived quantities
def read1D(runName, setup):
    # Read the 1D data
    f       = open(os.path.join(runName,'wind_1D.dat'),'r')
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
    data_1D['v_esc'] = sqrt(2*gg*setup['massAGB_ini']*solarm/(data_1D['r']))
    if setup['isink_radiation'] > 1:
        if setup['iget_tdust'] == 1:
            data_1D['Tdust'] = setup['primary_Teff']*(setup['primary_Reff']/data_1D['r'])**setup['tdust_exp']
        elif setup['iget_tdust'] == 0:
            data_1D['Tdust'] = data_1D['Tgas']
    del headers, m, i, column
    return data_1D

# Plots the 1D data
def plot1D (data, setup, references, whichPlot, ax1, ax2=None, second=False):
    if whichPlot == 'vel':
        if second == False:
            ax1.plot(data[references['x_axis']]/references['x_ref'], data['v'] /references['v_ref'],    color='black', linestyle='-',  label='v     analytic')
            ax1.plot(data[references['x_axis']]/references['x_ref'], data['cs']/references['v_ref'],    color='green', linestyle='--', label='sound speed')
            if setup['isink_radiation']>0 and 'alpha_rad' in setup:
              pass
            else:
              ax1.plot(data[references['x_axis']]/references['x_ref'], data['v_esc']/references['v_ref'], color='blue',  linestyle='-',  label='escape vel')
        else:
            ax1.plot(data[references['x_axis']]/references['x_ref'], data['v']/references['v_ref'],    color='black', linestyle='--')
        ax1.set_ylabel(r'$v$ [km s$^{-1}$]')

    if whichPlot == 'temp':
        ax1.plot(data[references['x_axis']]/references['x_ref'], data['T']/references['T_ref'], color='black',  linestyle='-',  label=r'T$_{gas}$  analytic')
        ax1.set_ylabel('Temperature [K]')
        if setup['isink_radiation'] > 1 and setup['iget_tdust'] > 0:
          ax1.plot(data[references['x_axis']]/references['x_ref'], data['Tdust']/references['T_ref'], color='magenta',  linestyle='-',  label=r'T$_{dust}$  analytic')
          
    if whichPlot == 'dustcool':
        lns11 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['T']/references['T_ref'], color='black',  linestyle='-',  label=r'T$_{gas}$  analytic')
        ax1.set_ylabel('Temperature [K]')
        lns12,lns2 = Dustcooling(data, references, ax1, ax2)
        lns1 = lns11 + lns12        
        return lns1, lns2
        
    if whichPlot == 'v&T':
       lns21 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['T']/references['T_ref'], color='grey',  linestyle='-',  label=r'T$_{gas}$  analytic')
       if setup['iget_tdust'] > 0:
          lns23 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['Tdust']/references['T_ref'], color='magenta',  linestyle='-',  label=r'T$_{dust}$  analytic')
          lns21 = lns21 + lns23
       lns11 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['cs']/references['v_ref'],    color='green', linestyle='--', label='sound speed')
       lns12  = ax1.plot(data[references['x_axis']]/references['x_ref'], data['v'] /references['v_ref'],    color='black', linestyle='-',  label='v     analytic')
       lns1  = lns11 + lns12
       lns2  = lns21
       ax1.set_ylabel(r'$v$ [km s$^{-1}$]')
       ax2.set_ylabel('Temperature [K]')
       return lns1, lns2
     
    if whichPlot == 'dust':
        if setup['isink_radiation'] > 1 and setup['iget_tdust'] > 0:
            lns1 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['Tdust']/references['T_ref'], color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1.set_ylabel('Temperature [K]')
        else:
            lns1 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['T']/references['T_ref'],     color='red',  linestyle='-', label=r'T$_{dust}$=T$_{gas}$ analytic')
            ax1.set_ylabel('Temperature [K]')
        if setup['idust_opacity'] == 2:
            lns21 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['kappa'], color='green',  linestyle='-', label=r'$\kappa_\mathrm{d}$ analytic')
            K3    = multiply(data['K3'],data['rho'])/mass_proton_cgs
            K3    = array([1e-95 if k == 0. else k for k in K3])
            lns22 = ax2.plot(data[references['x_axis']]/references['x_ref'], log10(K3), color='blue',  linestyle='-', label=r'log($\mathcal{K}_3$) analytic')
            lns2  = lns21 + lns22
            ax2.set_ylabel(r' $\kappa_\mathrm{d}$ [cm$^2$/g]  OR  log($\mathcal{K}_3$) [cm$^{-3}$]')
            return lns1, lns2
        elif setup['idust_opacity'] == 1:
            lns2 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['kappa'], color='green',  linestyle='-', label=r'$\kappa_\mathrm{d}$ analytic')
            ax2.set_ylabel(r' $\kappa_\mathrm{d}$ [cm$^2$/g]')
            return lns1, lns2
        else:
            return lns1

    if whichPlot == 'chem':
        if setup['iget_tdust'] > 0:
            lns1 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['Tdust']/references['T_ref'], color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1.set_ylabel('Temperature [K]')
        else:
            lns1 = ax1.plot(data[references['x_axis']]/references['x_ref'], data['T']/references['T_ref'],     color='red',  linestyle='-', label=r'T$_{dust}$ analytic')
            ax1.set_ylabel('Temperature [K]')
        if setup['idust_opacity'] == 2:
            lns21 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['gamma'],   color='green',   linestyle='-', label=r'$\gamma$ analytic')
            lns22 = ax2.plot(data[references['x_axis']]/references['x_ref'], data['mu'],      color='blue',    linestyle='-', label=r'$\mu$ analytic')
            lns23 = ax2.plot(data[references['x_axis']]/references['x_ref'], log10(data['S']), color='magenta', linestyle='-', label=r'log($S$) analytic')
            lns2  = lns21 + lns22 + lns23
            ax2.set_ylabel(r' $\gamma$  OR  $\mu$  OR  log($S$)')
            return lns1, lns2

    if whichPlot == 'tau':
        ax1.plot(data[references['x_axis']]/references['x_ref'], data['tau'], color='black',  linestyle='-',  label=r'$\tau$  analytic')
        ax1.set_ylabel('optical depth')

    if whichPlot == 'tau_lucy':
        ax1.plot(data[references['x_axis']]/references['x_ref'], data['tau_lucy'], color='black',  linestyle='-',  label=r'$\tau_{Lucy}$  analytic')
        ax1.set_ylabel('Lucy optical depth')

    if whichPlot == 'alpha':
        ax1.plot(data[references['x_axis']]/references['x_ref'], data['alpha'], color='black',  linestyle='-',  label=r'$\alpha$  analytic')
        ax1.set_ylabel('Acceleration parameter')

# specifies the references of the x-axis to which the data is plotted
def reference(axis, r_ref = 0, T_ref = 0, Tinj = 0):
    if axis == 'r':
        references = {
            'x_axis'   : axis,
            'x_ref'    : r_ref,
            'x_label'  : r'$Distance$ [$R_*]$',
            'x_in'     : 0.*r_ref,
            'x_out'    : 12.*r_ref,
            'x_limits' : [0.,12.]
        }
    elif axis == 'Tgas':
        references = {
            'x_axis'   : axis,
            'x_ref'    : T_ref,
            'x_label'  : r'$Temperature$ [$T_{inj}]$',
            'x_in'     : 1.1*Tinj,
            'x_out'    : 0.1*Tinj,
            'x_limits' : [1.1*Tinj/T_ref,0.1*Tinj/T_ref]
        }
    return references

# Plots auxiliary quantities that are not strictly related to the 1D or 3D SPH data
def plot_Auxiliary (data, setup, references, ax1):
    if references['x_axis'] == 'r':
        lns5 = ax1.axvline(x=r_last_shell(setup),ymin=0,ymax=3,color='orange',linestyle='-',label='last boundary shell')
        return lns5
    else:
      pass
    
# Determines the position of the last boundary shell
def r_last_shell (setup):
    d_tang       = setup['wind_inject_radius']*2*(constants.golden*sqrt(5))**(-1/2)/(2*setup['iwind_resolution']-1)
    d_rad        = setup['wind_shell_spacing']*d_tang
    r_last_shell = setup['wind_inject_radius']+setup['iboundary_spheres']*d_rad
    return r_last_shell

# Plots the 3D SPH data
def plot3D (dumpData, setup, references, whichPlot, ax1, ax2=None, second=False):
    if whichPlot == 'vel':
        if second == False:
           ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['speed']*1e5/references['v_ref'], 'r.', label='v     SPH')
        else:
          ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['speed']*1e5/references['v_ref'], 'r.')
          
    if whichPlot == 'temp':
        ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tgas']/references['T_ref'], 'r.', label=r'T$_{gas}$  SPH')
        if setup['isink_radiation'] > 1 and setup['iget_tdust'] > 0:
          ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tdust']/references['T_ref'], 'm.', label=r'T$_{dust}$  SPH')
        
    if whichPlot == 'v&T':
      lns3 = ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['speed']*1e5/references['v_ref'], 'r.', label='v     SPH')
      lns4 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tgas']/references['T_ref'], 'b.', label=r'T$_{gas}$  SPH')
      if setup['iget_tdust'] > 0:
        lns41 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tdust']/references['T_ref'], 'm.', label=r'T$_{dust}$  SPH')
        lns4 = lns4 + lns41
      return lns3, lns4
    
    if whichPlot == 'dust':
        if setup['isink_radiation'] > 1 and setup['iget_tdust'] > 0:
            lns3 = ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tdust']/references['T_ref'], 'r.', label=r'T$_{dust}$  SPH')
        else:
            lns3 = ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tgas']/references['T_ref'],  'r.', label=r'T$_{dust}$=T$_{gas}$  SPH')
        if setup['idust_opacity'] == 2:
            ax2=ax1.twinx()
            lns41 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['kappa'],  'g.', label=r'$\kappa_\mathrm{d}$  SPH')
            K3    = multiply(dumpData['K3'],dumpData['rho'])/mass_proton_cgs
            K3    = array([1e-95 if k == 0. else k for k in K3])
            lns42 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], log10(K3), 'b.', label=r'log($\mathcal{K}_3$)  SPH')
            lns4  = lns41 + lns42
            return lns3, lns4
        elif setup['idust_opacity'] == 1:
            lns4 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['kappa'], 'g.', label=r'$\kappa_\mathrm{d}$ SPH')
            ax2.set_ylabel(r' $\kappa_\mathrm{d}$ [cm$^2$/g]')
            return lns3, lns4
        else:
            return lns3

    if whichPlot == 'chem':
        if ['iget_tdust'] > 0:
            lns3 = ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tdust']/references['T_ref'], 'r.', label=r'T$_{dust}$  SPH')
        else:
            lns3 = ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['Tgas']/references['T_ref'],  'r.', label=r'T$_{dust}$  SPH')
        if ['idust_opacity'] == 2:            
            ax2   = ax1.twinx()
            lns41 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['gamma'],   'g.', label=r'$\gamma$ SPH')
            lns42 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['mu'],      'b.', label=r'$\mu$ SPH')
            lns43 = ax2.plot(dumpData[references['x_axis']]/references['x_ref'], log10(dumpData['S']), 'm.', label=r'log($S$) SPH')
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
        ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['tau'], 'r.', label=r'$\tau$  SPH')
        ax1.set_ylim([0,1.2*max(dumpData['tau'])])

    if whichPlot == 'tau_lucy':
        ax1.plot(dumpData[references['x_axis']]/references['x_ref'], dumpData['tauL'], 'r.', label=r'$\tau_{Lucy}$  SPH')
        ax1.set_ylim([0,0.7])

    if whichPlot == 'alpha':
        alpha = 0
        if setup['isink_radiation'] == 1:
            alpha = [setup['alpha_rad']]*len(dumpData[references['x_axis']])
        elif setup['isink_radiation'] == 2:
            alpha = dumpData['Gamma']
        elif setup['isink_radiation'] == 3:
            alpha = dumpData['Gamma'] + setup['alpha_rad']
        ax1.plot(dumpData[references['x_axis']]/references['x_ref'], alpha, 'r.', label=r'$\tau_{Lucy}$  SPH')
        ax1.set_ylabel('Acceleration parameter')
        ax1.set_ylim([0,1.2*max(alpha)])
    
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
def GaussianSigma(data1D,data3D,xaxis,references):
  
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
  
  outer_boundary = argmax(x_sph_sorted > references['x_out'])  
  
  MSE = sum(square(subtract(y_sph_sorted[:outer_boundary],y_sph_interp_sorted[:outer_boundary])))/len(y_sph_sorted[:outer_boundary])
  STD = sqrt(MSE)
  print('===============')
  print('Standard deviation of '+str(quantity)+'_sph w.r.t. '+str(quantity)+'_1D = ' +str(round(STD,2))+' km/s')
  print('===============')  


# Work in progress
def Dustcooling(data1D, references, ax1, ax2):
  
  Tprofile = loadtxt('temperatureprofile_AC_0.1um.txt')
  Tprofile = [*zip(*Tprofile)]
  
  radii        = array(Tprofile[0])
  Tdust        = array(Tprofile[1])  
  interp       = interpolate.interp1d(radii, Tdust, kind='cubic', bounds_error=False, fill_value='extrapolate')
  x_1D_profile = data1D[references['x_axis']]/au
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
  
  lns12 = ax1.plot(x_1D_profile*references['x_ref']/au, Tdust_interp, color='blue',  linestyle='-',  label=r'T$_{dust}$ MCFost')
  ax2  = ax1.twinx()
  lns21 = ax2.plot(x_1D_profile*references['x_ref']/au,CoolingBowen, color='red',  linestyle='--',  label='Cooling Bowen')
  lns22 = ax2.plot(x_1D_profile*references['x_ref']/au,CoolingFourier, color='brown',  linestyle='--',  label='Cooling Fourier')
  lns23 = ax2.plot(x_1D_profile*references['x_ref']/au,CoolingHollenbach, color='green',  linestyle='--',  label='Cooling Hollenbach')
  lns24 = ax2.plot(x_1D_profile*references['x_ref']/au,CoolingStefan, color='purple',  linestyle='--',  label='Cooling StefBoltz')
  lns2  = lns21 + lns22 + lns23 + lns24
  ax1.axvline(x=1.8, color='orange')
  ax2.set_ylabel(r'cooling rate $\Lambda$ [erg s$^{-1}$ cm$^{-3}$]')
  
  return lns12,lns2


# ===== GENERAL INPUT PARAMETERS ========================================================================
def profiles_main(run, loc, saveloc, dumpData, setup):
    if not setup['single_star']: return
    runName = os.path.join(loc,run)
    data1D = read1D(runName, setup)

    # specify reference values
    r_ref  = setup['primary_Reff']*au
    references = reference('r', r_ref=r_ref)
    references['T_ref']  = 1.
    references['v_ref']  = 1.e5

    whichPlots = ['vel', 'temp', 'v&T']
    if setup['idust_opacity'] > 0: whichPlots.append('dust')
    if setup['isink_radiation'] > 1 and setup['iray_resolution'] >=0: 
        if setup['iget_tdust'] == 3: whichPlots.append('tau_lucy')
        else: whichPlots.append('tau')
    if setup['isink_radiation'] > 0: whichPlots.append('alpha')
    
    # ===== GENERATE PLOT ========================================================================
    for whichPlot in whichPlots:
        # ===== SET UP PLOT LAYOUT =============================================================
        fig, ax1 = plt.subplots(1,figsize=(9,6),tight_layout=True,squeeze=False,sharex=True)
        ax1 = ax1[0][0]

        if whichPlot == 'dust' or whichPlot == 'chem' or whichPlot == 'v&T' or whichPlot == 'dustcool':
            ax2 = ax1.twinx()
            if setup['idust_opacity'] > 0:
                lns1,lns2 = plot1D(data1D, setup, references, whichPlot, ax1, ax2)
                lns3,lns4 = plot3D(dumpData, setup, references, whichPlot, ax1, ax2)
                lns = lns1+lns2+lns3+lns4#+[lns5]
            else:
                lns1      = plot1D(data1D, setup, references, whichPlot, ax1, ax2)
                lns3      = plot3D(dumpData, setup, references, whichPlot, ax1, ax2)
                lns = lns1+lns3#+[lns5]
            labs = [l.get_label() for l in lns]
            ax1.legend(lns, labs, bbox_to_anchor=(1, 0.30), loc='lower right',fontsize=13)
            for item in (  [ax2.title, ax2.xaxis.label, ax2.yaxis.label]
                    + ax2.get_xticklabels() + ax2.get_yticklabels()
                        ):
                item.set_fontsize(16)
        else:
            plot3D(dumpData, setup, references, whichPlot, ax1)
            plot1D(data1D, setup, references, whichPlot, ax1)
            plot_Auxiliary(data1D, setup, references, ax1)
            ax1.legend(           bbox_to_anchor=(1, 0.2), loc='lower right', fontsize=13)
        ax1.set_xlabel(references['x_label'])
            
        #---------- Font size ----------
        for item in (  [ax1.title, ax1.xaxis.label, ax1.yaxis.label]
                    + ax1.get_xticklabels() + ax1.get_yticklabels()
                    ):
            item.set_fontsize(16)
        
        #---------- Finalize plot ----------
        #txt = ['high resolution, q=20','mid resolution, q=10','low resolution, q=3']
        #ax1.text(2.6,19,txt,fontsize=16)
        ax1.set_xlim(references['x_limits'])
        #ax2.set_ylim([0,40])
        #ax2.set_yscale('log')
        plt.title('BowenDust model',fontsize=16)
        plt.savefig(os.path.join(saveloc, 'png/1DProfile_'+whichPlot+'.png'), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(saveloc, 'pdf/1DProfile_'+whichPlot+'.pdf'), bbox_inches="tight")
