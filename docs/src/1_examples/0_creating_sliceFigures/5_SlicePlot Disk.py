import numpy                          as np
import matplotlib.pyplot              as plt
import matplotlib.animation           as ani
import os

# import plons scripts
import plons
import plons.SmoothingKernelScript    as sk
import plons.PhysicalQuantities       as pq
import plons.ConversionFactors_cgs    as cgs
import plons.Plotting                 as plot


def makeMovie(prefix, loc, outputloc, model, firstframe,lastframe):

    dump = os.path.join(loc,"wind_00001")

    setup     = plons.LoadSetup(loc, prefix)
    dumpData  = plons.LoadFullDump(dump, setup)

    fig, ax = plt.subplots(1)#, figsize=(10, 6))
    ax.set_aspect('equal')
    ax.set_facecolor('k')

    n = 500
    # Filmpjes by Mats: n=500, dx = 10 au, mijn HI cooling dx = 20
    dx = 10*cgs.au
    x = np.linspace(-dx, dx, n)
    x = np.linspace(-7*cgs.au, 12*cgs.au, n)  #voor eccentrische accrSchijf
    # filmpjes Mats: verhouding 1.5
    y = np.linspace(-dx/1.5, dx/1.5, n)
    # y = np.linspace(-dx, 2*dx, n)

    X, Y = np.meshgrid(x, y)
    Z    = np.zeros_like(X)
    theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
    X_rot, Y_rot, Z_rot = sk.rotateMeshAroundZ(theta, X, Y, Z)
    smooth_rot = sk.smoothMesh(X_rot, Y_rot, Z_rot, dumpData, ['rho'])

    #Filmpjes Mats: inferno; -17 -- -12
    # mesh = ax.pcolormesh(X/cgs.au, Y/cgs.au, np.log10(smooth_rot["rho"]+1e-99), cmap=plt.colormaps['gist_heat'], vmin=-19, vmax = -14)
    mesh = ax.pcolormesh(X/cgs.au, Y/cgs.au, np.log10(smooth_rot["rho"]+1e-99), cmap=plt.colormaps['gist_heat'], vmin=-18, vmax = -11.5)
    ax.set_xlim(x[0]/cgs.au, x[-1]/cgs.au)
    ax.set_ylim(y[0]/cgs.au, y[-1]/cgs.au)
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')

        

    cbar = plt.colorbar(mesh, ax = ax, location='right', fraction=0.0471, pad=0.01)
    cbar.set_label(r'log density $\, \rm{/(g \cdot cm^{-3})}$')

    circleAGB, circleComp = plot.plotSink(ax, dumpData, setup)

    plt.tight_layout()

    def animate(frame):
        dump = os.path.join(loc,prefix+f"_%05i"%frame)
        dumpData  = plons.LoadFullDump(dump, setup)
        theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
        X_rot, Y_rot, Z_rot = sk.rotateMeshAroundZ(theta, X, Y, Z)
        smooth_rot = sk.smoothMesh(X_rot, Y_rot, Z_rot, dumpData, ['rho'])

        mesh.set_array(np.log10(smooth_rot["rho"]+1e-99))
        circleAGB.center = -np.linalg.norm(dumpData['posAGB'])/cgs.au, 0, 0
        circleComp.center = np.linalg.norm(dumpData['posComp'])/cgs.au, 0, 0
        print(frame, end="\r")

    anim = ani.FuncAnimation(fig, animate, frames=range(firstframe, lastframe+1), interval=100)

    anim.save(os.path.join(outputloc,model+'_last4_Orbits_diffColbar.mp4'), writer='ffmpeg')


prefix = "wind"
outputloc = "/lhome/jolienm/Documents/TierModels/finalModelsAccrDisks/plonsMovies/"
# for full last orbit
firstframe = 350
lastframe = 600
model = 'v20e50'
loc = "/STER/hydroModels/jolienm/finalModelsAccrDisks/"+model
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)
model = 'v05e50'
loc = "/STER/hydroModels/jolienm/finalModelsAccrDisks/"+model
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)
model = 'v10e50'
loc = "/STER/hydroModels/jolienm/finalModelsAccrDisks/"+model
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)
'''
## HI cooling binaries
prefix = "wind"
outputloc = "/lhome/jolienm/Documents/TierModels/finalModelsHIcooling/plonsMovies/"
# for full last orbit
firstframe = 1
lastframe = 600

loc = "/STER/hydroModels/jolienm/finalModelsHIcoolingBinaries/v10e50/"
model = 'v10e50'
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)

loc = "/STER/hydroModels/jolienm/finalModelsHIcoolingBinaries/v05e50/"
model = 'v05e50'
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)

loc = "/STER/hydroModels/jolienm/finalModelsHIcoolingBinaries/v20e50/"
model = 'v20e50'
makeMovie(prefix,loc,outputloc,model,firstframe,lastframe)
'''