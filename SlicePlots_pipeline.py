
def SlicePlots(runNumber,directory):
    #Import packages
    import numpy as np
    import math  as math
    import healpy as hp
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import os

    # import own scripts
    import loadPhantomData as ld
    import smoothingKernelScript as sk
    import ConversionFactors_cgs as cgs

    # import certain things from packages
    from matplotlib    import colors
    from astropy       import constants
    from mpl_toolkits.axes_grid1 import AxesGrid

    from matplotlib    import rcParams, rc
    # Change the matplotlib default parameters
    rcParams.update({'font.size':   11})
    rcParams.update({'figure.dpi': 200})
    rc('font', family='serif')
    rc('text', usetex=True)

    AU = cgs.AU_cm()						# cm


    # In[3]:


    try:
        os.mkdir(directory+'2DSliceplots/')
    except OSError:
        print('')


    # In[4]:


    # Make figure with the x-y(left) and x-z(right) slice plots of log(rho[g/cm3]), log(T[K]) and |v|[km/s]. Takes 'zoom'factor as input
    def finalPlots(modelname,results_sph_sl_y,results_sph_sl_z,x1,y1,x2,z2,zoom,rhoMin,rhoMax,vmax):
        cm_rho  = plt.cm.get_cmap('inferno')
        cm_T    = plt.cm.get_cmap('afmhot')
        cm_v    = plt.cm.get_cmap('viridis')
    #     fig = plt.figure(figsize=(11, 14))
        fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6))= plt.subplots(3, 2,  gridspec_kw={'height_ratios':[1,1,1],'width_ratios': [0.8,1]})
        fig.set_size_inches(10, 14)
        
        lim = 200/zoom        

        def onePlot(axis,par,name,colormap,mi,ma,x,y,xlabel,ylabel):
            axi = axis
            axi.axis('equal')
            axi.set_facecolor('k')
            axPlot = axi.scatter(x/AU,y/AU,s=11,c=par,cmap=colormap,vmin=mi, vmax = ma)
            if axi == ax2 or axi == ax4 or axi == ax6:
                axi.plot(xAGB,zAGB, 'bo', markersize = 4,label = 'AGB')
                axi.plot(xcomp,zcomp, 'ro', markersize = 3,label = 'comp')  
                cbar = plt.colorbar(axPlot, ax = axi)
                cbar.set_label(name, fontsize = 18)
            if axi == ax5 or axi == ax6:

                axi.set_xlabel(xlabel)
            if axi == ax1 or axi == ax3 or axi == ax5:
                axi.plot(xAGB,yAGB, 'bo', markersize = 4,label = 'AGB')
                axi.plot(xcomp,ycomp, 'ro', markersize = 3,label = 'comp')  
            axi.axis([-lim,lim,-lim,lim])
            axi.set_ylabel(ylabel)
        # the temperature colorbar limits may have to be changed...
        onePlot(ax1,np.log10(results_sph_sl_z['rho']),'log($\\rho[$g/cm$^3]$)', cm_rho, rhoMin, rhoMax, x1,y1,'x[AU]', 'y[AU]')
        onePlot(ax2,np.log10(results_sph_sl_y['rho']),'log($\\rho[$g/cm$^3]$)', cm_rho, rhoMin, rhoMax, x2,z2,'x[AU]', 'z[AU]')
        onePlot(ax3,np.log10(results_sph_sl_z['temp']),'log$(T$[K])', cm_T, 1.7, 5.2,x1,y1,'x[AU]', 'y[AU]')
        onePlot(ax4,np.log10(results_sph_sl_y['temp']),'log$(T$[K])', cm_T, 1.7, 5.2,x2,z2, 'x[AU]', 'z[AU]')
        onePlot(ax5,results_sph_sl_z['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, vmax, x1,y1, 'x[AU]', 'y[AU]')
        onePlot(ax6,results_sph_sl_y['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, vmax, x2,z2, 'x[AU]', 'z[AU]')
        # fig.show()
        fig.savefig(directory+'2DSliceplots/'+str(runNumber)+'Z'+str(zoom)+'.png',dpi = 300)
        print('Slice plots model '+str(runNumber)+' ready and saved!')


    # In[5]:


    # Make figure with the x-y(orbital plane) slice plot of log(rho[g/cm3]). Takes 'zoom'factor as input

    def finalPlots1(modelname,results_sph_sl_z,x1,y1,zoom,rhoMin,rhoMax,vmax):
        cm_rho  = plt.cm.get_cmap('inferno')

        fig, (ax1)= plt.subplots(1, 1,  gridspec_kw={'height_ratios':[1],'width_ratios': [1]})
        fig.set_size_inches(9.8, 8)
        
        #CHANGE 200 TO OUTER BOUNDARY
        lim = 200/zoom        

        def onePlot(axis,par,name,colormap,mi,ma,x,y,xlabel,ylabel):
            axi = axis
            axi.axis('equal')
            axi.set_facecolor('k')
            axPlot = axi.scatter(x/AU,y/AU,s=11,c=par,cmap=colormap,vmin=mi, vmax = ma)
            axi.plot(xAGB,yAGB, 'bo', markersize = 5,label = 'AGB')
            axi.plot(xcomp,ycomp, 'ro', markersize = 5,label = 'comp')   
            cbar = plt.colorbar(axPlot, ax = axi)
            cbar.set_label(name, fontsize = 25)
            cbar.ax.tick_params(labelsize=14)
            axi.set_xlabel(xlabel, fontsize = 18)

            axi.set_ylabel(ylabel, fontsize = 18)

        onePlot(ax1,np.log10(results_sph_sl_z['rho']),'log($\\rho$[g/cm$^3]$)', cm_rho, rhoMin, rhoMax, x1,y1,'x[AU]', 'y[AU]')

        ax1.tick_params(labelsize=14)

        
        # fig.show()
        fig.savefig(directory+'2DSliceplots/1Plot_'+str(runNumber)+'Z'+str(zoom)+'.png',dpi = 300)
        print('Density slice plot model '+str(runNumber)+' ready and saved!')

        
    # In[6]:


    # Make figure with the x-y(left) and x-z(right) Hill sphere slice plots of log(rho[g/cm3]) and |v|[km/s].

    def HillsphPlot(modelname,results_sph_sl_y,x2,z2,zoom,rhoMin,rhoMax):
        cm_rho  = plt.cm.get_cmap('inferno')
        cm_T    = plt.cm.get_cmap('afmhot')
        cm_v    = plt.cm.get_cmap('viridis')
        cm_vtvv = plt.cm.get_cmap('seismic')
    #     fig = plt.figure(figsize=(11, 14))
        fig, ((ax1,ax2),(ax3,ax4))= plt.subplots(2, 2,  gridspec_kw={'height_ratios':[1,1],'width_ratios': [1,1]})
        fig.set_size_inches(10, 8)
        
        #CHANGE 200 TO OUTER BOUNDARY
        lim = 200/zoom        

        def onePlotRH(axis,par,name,colormap,mi,ma,x,y,xlabel,ylabel):
            axi = axis
            axi.axis('equal')
            axi.set_facecolor('k')
            axPlot = axi.scatter(x/AU,y/AU,s=11,c=par,cmap=colormap,vmin=mi, vmax = ma)
            cbar = plt.colorbar(axPlot, ax = axi)
            cbar.set_label(name, fontsize = 18)
            if axi == ax3 or axi == ax4:
                axi.set_xlabel(xlabel)
            if axi == ax1 or axi == ax3:
                axi.set_ylabel(ylabel)
    #         axi.plot(0,0, 'ro', markersize = 5,label = 'comp')
            axi.axis([-lim,lim,-lim,lim])
            axi.plot(0,0, 'ro', markersize = 5,label = 'comp')
     
        onePlotRH(ax1,np.log10(results_sph_sl_y['rho']),'log($\\rho$[g/cm$^3$])', cm_rho, rhoMin, rhoMax, x2,z2,'x[AU]', 'z[AU]')
        onePlotRH(ax2,np.log10(results_sph_sl_y['temp']),'log($T$[K])', cm_T, 3, 6,x2,z2, 'x[AU]', 'z[AU]')
        onePlotRH(ax3,results_sph_sl_y['speed']*1e-5,'$|v|$[km/s]', cm_v, 0, 26, x2,z2, 'x[AU]', 'z[AU]')
        onePlotRH(ax4,results_sph_sl_y['vtvv'],'$(v_t/v)^2$', cm_vtvv, 0, 1, x2, z2, 'x[AU]', 'z[AU]' )
        # fig.show()
        zoom = '_RH'
        fig.savefig(directory+'2DSliceplots/HillSphere_'+str(runNumber)+'Z'+str(zoom)+'.png',dpi = 300)
        print('Hill sphere plot model '+str(runNumber)+' ready and saved!')


    # In[7]:


    #Load data 
    data = {}
    data[str(runNumber)] = ld.LoadData_cgs(str(runNumber)) 

    print('Finished loading in model ', runNumber)


    # In[8]:


    # # Uncomment if we add Hill sphere plots
    # #Load data for Hill sphere
    # dataHS = {}
    # xcomp = data[str(runNumber)]['compCoord'][0]/AU
    # ycomp = data[str(runNumber)]['compCoord'][1]/AU
    # zcomp = data[str(runNumber)]['compCoord'][2]/AU
    # hillsph = data[str(runNumber)]['rHill']/AU
    # # bound = data[str(runNumber)]['boundary']/AU
    # # zoom = int(bound/hillsph)
    # # dataHS[str(runNumber)] = ld.LoadData_cgs_inner(str(runNumber),hillsph,xcomp*AU,ycomp*AU,zcomp*AU)
    # #Make Hill sphere plots 
    # dataRunHS  = dataHS[str(runNumber)] #rcomp is now (0,0,0)
    # dataRun    = data[str(runNumber)]
    # xAGB = dataRun['AGBcoord'][0]/AU
    # yAGB = dataRun['AGBcoord'][1]/AU
    # zAGB = dataRun['AGBcoord'][2]/AU
    # xcomp = dataRun['compCoord'][0]/AU
    # ycomp = dataRun['compCoord'][1]/AU
    # zcomp = dataRun['compCoord'][2]/AU
    # Mdot  = dataRun['Mdot']
    # v_ini = dataRun['v_ini'] *1e-5
    # hillsph = dataRun['rHill']/AU
    # bound = dataRun['boundary']/AU
    # zoom = int(bound/hillsph)
    # results_sph_sl_yHS,x2HS,y2HS,z2HS = sk.getSmoothingKernelledPix(400,20,dataRunHS,['rho','temp','speed','vtvv'], 'comp','y',10*AU)

    # #Change these limits if the plots don't look nice
    # if Mdot > 1e-5:
    #     if v_ini == 5:
    #         rhoMin = -14.5
    #         rhoMax = -11.5
    #     elif v_ini == 20:
    #         rhoMin = -15.5
    #         rhoMax = -12.5

    # elif Mdot < 5e-7:
    #     if v_ini == 5:
    #         rhoMin = -17
    #         rhoMax = -14
    #     elif v_ini == 20:
    #         rhoMin = -18
    #         rhoMax = -15

    # elif 5e-7 <= Mdot <= 1e-5:
    #     rhoMin = -16
    #     rhoMax = -13

    # HillsphPlot('M'+str(runNumber),results_sph_sl_yHS,x2HS,z2HS,zoom,rhoMin,rhoMax)


    # In[9]:


    # Make sliceplots
    dataRun    = data[str(runNumber)]
    xAGB = dataRun['AGBcoord'][0]/AU
    yAGB = dataRun['AGBcoord'][1]/AU
    zAGB = dataRun['AGBcoord'][2]/AU
    xcomp = dataRun['compCoord'][0]/AU
    ycomp = dataRun['compCoord'][1]/AU
    zcomp = dataRun['compCoord'][2]/AU
    Mdot  = dataRun['Mdot']

    #Set the limits of the colorbars
    if Mdot > 1e-5:
        rhoMin = -18.5
        rhoMax = -13.5
    elif Mdot < 5e-7:
        rhoMin = -21
        rhoMax = -16
    elif 5e-7 <= Mdot <= 1e-5:
        rhoMin = -20
        rhoMax = -15
    vmax = 26

    #Load the smoothing kernel data
    #Uncomment the ones you need, takes long to run, so don't run if not nescessary
    #zoom = 1
    results_sph_sl_z,x1,y1,z1 = sk.getSmoothingKernelledPix(400,20,dataRun,['rho','temp','speed'], 'comp','z',200*AU)
    results_sph_sl_y,x2,y2,z2 = sk.getSmoothingKernelledPix(400,20,dataRun,['rho','temp','speed'], 'comp','y',200*AU)
    #zoom =2
    results_sph_sl_zZ2,x1Z2,y1Z2,z1Z2 = sk.getSmoothingKernelledPix(400,20,dataRun,['rho','temp','speed'], 'comp','z',100*AU)
    results_sph_sl_yZ2,x2Z2,y2Z2,z2Z2 = sk.getSmoothingKernelledPix(400,20,dataRun,['rho','temp','speed'], 'comp','y',100*AU)
    #zoom = 5
    #     results_sph_sl_zZ5,x1Z5,y1Z5,z1Z5 = sk.getSmoothingKernelledPix(300,20,dataRun,['rho','temp','speed'], 'comp','z',40*AU)
    #     results_sph_sl_yZ5,x2Z5,y2Z5,z2Z5 = sk.getSmoothingKernelledPix(300,20,dataRun,['rho','temp','speed'], 'comp','y',40*AU)

    #Make plots
    finalPlots1('M'+str(runNumber),results_sph_sl_z, x1,y1,1,rhoMin,rhoMax,vmax)
    finalPlots1('M'+str(runNumber),results_sph_sl_zZ2, x1Z2,y1Z2,2,rhoMin,rhoMax,vmax)
    finalPlots('M'+str(runNumber),results_sph_sl_y,results_sph_sl_z,x1,y1,x2,z2,1,rhoMin,rhoMax,vmax)
    finalPlots('M'+str(runNumber),results_sph_sl_yZ2,results_sph_sl_zZ2,x1Z2,y1Z2,x2Z2,z2Z2,2,rhoMin,rhoMax,vmax)

