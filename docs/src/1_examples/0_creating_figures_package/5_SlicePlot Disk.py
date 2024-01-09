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

prefix = "wind"
loc = "/STER/hydroModels/jolienm/modelsLionel/finalAccrDisks/v10e50_T3000_res8_racc01"
dump = os.path.join(loc,"wind_00001")

setup     = plons.LoadSetup(loc, prefix)
dumpData  = plons.LoadFullDump(dump, setup)

fig, ax = plt.subplots(1, figsize=(10, 6))
ax.set_aspect('equal')
ax.set_facecolor('k')

n = 500
dx = 10*cgs.au
x = np.linspace(-dx, dx, n)
y = np.linspace(-dx/1.5, dx/1.5, n)
X, Y = np.meshgrid(x, y)
Z    = np.zeros_like(X)
theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])
X_rot, Y_rot, Z_rot = sk.rotateMeshAroundZ(theta, X, Y, Z)
smooth_rot = sk.smoothMesh(X_rot, Y_rot, Z_rot, dumpData, ['rho'])

mesh = ax.pcolormesh(X/cgs.au, Y/cgs.au, np.log10(smooth_rot["rho"]+1e-99), cmap=plt.cm.get_cmap('inferno'), vmin=-17, vmax = -12)
ax.set_xlim(x[0]/cgs.au, x[-1]/cgs.au)
ax.set_ylim(y[0]/cgs.au, y[-1]/cgs.au)

cbar = plt.colorbar(mesh, ax = ax, location='right', fraction=0.0471, pad=0.01)
cbar.set_label("log(rho)")

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

anim = ani.FuncAnimation(fig, animate, frames=range(1, 301), interval=100)
anim.save('animation.mp4', writer='ffmpeg')
