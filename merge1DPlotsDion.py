import matplotlib.pyplot as plt
import numpy as np
import os
import os.path
import LoadDataPHANTOM as ld
import sys
sys.path.insert(1, '/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/mine')
from fortranRoutines import interpolate3D, coolingRate, findLowerIndicesCoolingTable
import math
from matplotlib.ticker import MultipleLocator


AU                  = 1.495E13
mu                  = 2.381
R                   = 8.314E7
kB                  = 1.38E-16
r_orbit_comp        = 4.
atomic_mass_unit    = 1.67E-24
mass_molecules      = np.array([28.01, 18.01528, 27.0253]) * atomic_mass_unit
abundance_molecules = np.array([1E-4, 5E-5, 3E-7])

colorCool = "C2"
colorAdia = "C1"
colorIso = "C0"
xmin = 1.1
xmax = 30
outputPath  = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output/mergedPlots"
pathCDTable = r"files/table_cd.dat"
pathRadCool = r"files/radcool_all.dat"
pathAscii   = r"files/wind_00078.ascii"

def main():
    pathCoolLocal  = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output/cooling/local/txt/data_1D_radialStructure.txt"
    pathAdiaLocal = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output/gamma/non_isowind/gamma1.4/ward/local/txt/data_1D_radialStructure.txt"
    pathIsoLocal  = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output/gamma/isowind/gamma1/ward/local/txt/data_1D_radialStructure.txt"

    fileCool = np.loadtxt(pathCoolLocal, skiprows=3)
    fileAdia = np.loadtxt(pathAdiaLocal, skiprows=3)
    fileIso  = np.loadtxt(pathIsoLocal , skiprows=3)

    radCoolTable = loadRadCoolTable()
    CDTable      = loadCDTable()

    pressureGradient(fileCool, fileAdia, fileIso)
    dumpData = loadDump()
    lambda_cool_x = findLambda(fileCool, dumpData, "x", radCoolTable, CDTable)
    lambda_cool_z = findLambda(fileCool, dumpData, "z", radCoolTable, CDTable)

    print("Average cooling x: %10.5e"%np.mean(10.**lambda_cool_x[lambda_cool_x > -999.]))
    print("Average cooling z: %10.5e"%np.mean(10.**lambda_cool_z[lambda_cool_z > -999.]))

    plot1DX(fileCool, fileAdia, fileIso, lambda_cool_x)
    plot1DZ(fileCool, fileAdia, fileIso, lambda_cool_z)

def plot1DX(fileCool, fileAdia, fileIso, lambda_cool):
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(18, 10))

    # Ax 1: Density
    ax1.plot(fileAdia[:,0]/AU, fileAdia[:,2], color=colorAdia, linestyle="-", label="Adiabatic")
    ax1.plot(fileIso[:,0]/AU ,  fileIso[:,2], color=colorIso, linestyle="-", label="Isothermal")
    ax1.plot(fileCool[:,0]/AU, fileCool[:,2], color=colorCool, linestyle="-", label="DDWC")

    ax1.axvline(4., color="black", linestyle="--", label=r"$x_\mathrm{comp}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(1E-21, 1E-10)
    ax1.legend(loc="upper center", bbox_to_anchor=[0.04, 0.3, 2.1, 1.], ncol=4, fontsize=22, labelspacing=2)

    # Ax 2: Temperature
    ax2.plot(fileCool[:, 0] / AU, fileCool[:, 6], color=colorCool, linestyle="-")
    ax2.plot(fileAdia[:, 0] / AU, fileAdia[:, 6], color=colorAdia, linestyle="-")
    ax2.axhline(2000., color=colorIso)

    ax2.axvline(4., color="black", linestyle="--")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(1E2, 1E7)

    # Ax 3: Velocity
    ax3.plot(fileCool[:, 0] / AU, fileCool[:, 4], color=colorCool, linestyle="-")
    ax3.plot(fileAdia[:, 0] / AU, fileAdia[:, 4], color=colorAdia, linestyle="-")
    ax3.plot(fileIso[:, 0] / AU, fileIso[:, 4], color=colorIso, linestyle="-")
    ax3.axvline(4., color="black", linestyle="--")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(1E0, 3E2)

    # Ax 4: Cooling Rate
    ax4.plot(fileCool[:, 0] / AU, lambda_cool, color=colorCool, linestyle="-")

    ax4.axvline(4., color="black", linestyle="--")
    ax4.set_xscale("log")
    # ax4.set_yscale("log")
    ax4.set_xlim(xmin, xmax)
    ax4.set_ylim(-24, -16)
    ax4.yaxis.set_major_locator(MultipleLocator(2))
    ax4.yaxis.set_minor_locator(MultipleLocator(0.5))

    ax1.set_ylabel(r"$\rho$ [g$\,$cm$^{-3}$]", fontsize=22)
    ax2.set_ylabel(r"$T$ [K]", fontsize=22)
    ax3.set_ylabel(r"$v$ [km/s]", fontsize=22)
    ax4.set_ylabel("$\log \ \Lambda / \mathrm{erg} \cdot \mathrm{s}^{-1} \cdot \mathrm{cm}^{-3}$", fontsize=22)

    ax3.set_xlabel("$x$ [AU]", fontsize=22)
    ax4.set_xlabel("$x$ [AU]", fontsize=22)

    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

    ax1.tick_params(labelsize=18)
    ax2.tick_params(labelsize=18)
    ax3.tick_params(pad=9, labelsize=18)
    ax4.tick_params(pad=9, labelsize=18)

    plt.subplots_adjust(hspace=0.1)

    #plt.show()
    fig.savefig(os.path.join(outputPath, "pdf/hydro_x.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outputPath, "png/hydro_x.png"), bbox_inches="tight")

def plot1DZ(fileCool, fileAdia, fileIso, lambda_cool):
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(18, 10))

    z_OL = 31

    # Ax 1: Density
    ax1.plot(fileAdia[:, 1] / AU, fileAdia[:, 3], color=colorAdia, linestyle="-", label="Adiabatic")
    ax1.plot(fileIso[:, 1] / AU, fileIso[:, 3], color=colorIso, linestyle="-", label="Isothermal")
    ax1.plot(fileCool[:, 1] / AU, fileCool[:, 3], color=colorCool, linestyle="-", label="DDWC")
    #ax1.axvline(z_OL, color="black", linestyle="--", label=r"$z_\mathrm{OL}$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(1E-21, 1E-10)
    ax1.legend(loc="upper center", bbox_to_anchor=[0.04, 0.3, 2.1, 1.], ncol=4, fontsize=22, labelspacing=2)

    # Ax 2: Temperature
    ax2.axhline(2000., color=colorIso)
    ax2.plot(fileCool[:, 1] / AU, fileCool[:, 7], color=colorCool, linestyle="-")
    ax2.plot(fileAdia[:, 1] / AU, fileAdia[:, 7], color=colorAdia, linestyle="-")
    #ax2.axvline(z_OL, color="black", linestyle="--")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(1E2, 1E7)

    # Ax 3: Velocity
    ax3.plot(fileCool[:, 1] / AU, fileCool[:, 5], color=colorCool, linestyle="-")
    ax3.plot(fileAdia[:, 1] / AU, fileAdia[:, 5], color=colorAdia, linestyle="-")
    ax3.plot(fileIso[:, 1] / AU, fileIso[:, 5], color=colorIso, linestyle="-")
    #ax3.axvline(z_OL, color="black", linestyle="--")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(1E0, 3E2)

    # Ax 4: Cooling Rate
    ax4.plot(fileCool[:, 1] / AU, lambda_cool, color=colorCool, linestyle="-")

    #ax4.axvline(z_OL, color="black", linestyle="--")
    ax4.set_xscale("log")
    # ax4.set_yscale("log")
    ax4.set_xlim(xmin, xmax)
    ax4.set_ylim(-24, -16)
    ax4.yaxis.set_major_locator(MultipleLocator(2))
    ax4.yaxis.set_minor_locator(MultipleLocator(0.5))


    ax1.set_ylabel(r"$\rho$ [g$\,$cm$^{-3}$]", fontsize=22)
    ax2.set_ylabel(r"$T$ [K]", fontsize=22)
    ax3.set_ylabel(r"$v$ [km/s]", fontsize=22)
    ax4.set_ylabel("$\log \ \Lambda / \mathrm{erg} \cdot \mathrm{s}^{-1} \cdot \mathrm{cm}^{-3}$", fontsize=22)

    ax3.set_xlabel("$z$ [AU]", fontsize=22)
    ax4.set_xlabel("$z$ [AU]", fontsize=22)

    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

    ax1.tick_params(labelsize=18)
    ax2.tick_params(labelsize=18)
    ax3.tick_params(pad=9, labelsize=18)
    ax4.tick_params(pad=9, labelsize=18)


    plt.subplots_adjust(hspace=0.1)

    fig.savefig(os.path.join(outputPath, "pdf/hydro_z.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outputPath, "png/hydro_z.png"), bbox_inches="tight")

def pressureGradient(fileCool, fileAdia, fileIso):
    fig, ax = plt.subplots(figsize=(10, 5))
    z_OL = 2.75

    # Gradients in x-direction
    p_cool_x = R / mu * fileCool[:, 2] * fileCool[:, 6]
    p_cool_grad_x = np.gradient(p_cool_x, fileCool[:, 0])

    p_adia_x = R / mu * fileAdia[:, 2] * fileAdia[:, 6]
    p_adia_grad_x = np.gradient(p_adia_x, fileAdia[:, 0])

    p_iso_x = R / mu * fileIso[:, 2] * 2000.
    p_iso_grad_x = np.gradient(p_iso_x, fileIso[:, 0])

    # Gradients in z-direction
    p_cool_z = R / mu * fileCool[:, 3] * fileCool[:, 7]
    p_cool_grad_z = np.gradient(p_cool_z, fileCool[:, 1])

    p_adia_z = R / mu * fileAdia[:, 3] * fileAdia[:, 7]
    p_adia_grad_z = np.gradient(p_adia_z, fileAdia[:, 1])

    p_iso_z = R / mu * fileIso[:, 3] * 2000.
    p_iso_grad_z = np.gradient(p_iso_z, fileIso[:, 1])

    ax.plot(fileAdia[:, 1] / AU, p_adia_grad_z, color=colorAdia, linestyle="-", label="Adiabatic")
    ax.plot(fileIso[:, 1] / AU, p_iso_grad_z, color=colorIso, linestyle="-", label="Isothermal")
    ax.plot(fileCool[:, 1] / AU, p_cool_grad_z, color=colorCool, linestyle="-", label="DDWC")

    ax.axvline(z_OL, color="black", linestyle="--", label=r"$z_\mathrm{OL}$")
    ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.set_xlim(xmin, xmax)
    # ax.set_ylim(1E-12, 1E-4)
    ax.set_ylabel(r"$\nabla p$ [Ba cm$^{-1}$]")
    ax.set_xlabel(r"$z$ [AU]")

    z_iso_ol  = fileIso[:, 1][np.logical_and((fileIso[:, 1] / AU) <= z_OL, fileIso[:, 1] / AU >= xmin)]
    z_adia_ol = fileAdia[:, 1][np.logical_and((fileAdia[:, 1] / AU) <= z_OL, fileAdia[:, 1] / AU >= xmin)]
    z_cool_ol = fileCool[:, 1][np.logical_and((fileCool[:, 1] / AU) <= z_OL, fileCool[:, 1] / AU >= xmin)]

    gradp_iso_ol  = p_iso_grad_z[np.logical_and((fileIso[:, 1] / AU) <= z_OL, (fileIso[:, 1] / AU) >= xmin)]
    gradp_adia_ol = p_adia_grad_z[np.logical_and((fileAdia[:, 1] / AU) <= z_OL, fileAdia[:, 1] / AU >= xmin)]
    gradp_cool_ol = p_cool_grad_z[np.logical_and((fileCool[:, 1] / AU) <= z_OL, fileCool[:, 1] / AU >= xmin)]

    average_pgrad_iso  = np.trapz(gradp_iso_ol, z_iso_ol) / ((xmax - xmin) * AU)
    average_pgrad_adia = np.trapz(gradp_adia_ol, z_adia_ol) / ((xmax - xmin) * AU)
    average_pgrad_cool = np.trapz(gradp_cool_ol, z_cool_ol) / ((xmax - xmin) * AU)

    print("nabla p isothermal: %10.5e"%average_pgrad_iso)
    print("nabla p adiabatic : %10.5e"%average_pgrad_adia)
    print("nabla p DDWC      : %10.5e"%average_pgrad_cool)

    ax.legend(loc="center right", bbox_to_anchor=[0., 0., 1., 1.], labelspacing=2)
    fig.savefig(os.path.join(outputPath, "pdf/pgrad_z.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outputPath, "png/pgrad_z.png"), bbox_inches="tight")

def loadRadCoolTable():
    file = np.loadtxt(pathRadCool)
    nx = int(file[:, 0][-1])
    ny = int(file[:, 1][-1])
    nz = int(file[:, 2][-1])

    data = np.zeros(shape=(nx, ny, nz, 6))
    for line in file:
        i, j, k, T, n_H, N, lambda_CO, lambda_H2O, lambda_HCN = line
        data[int(i) - 1][int(j) - 1][int(k) - 1] = [T, n_H, N, lambda_CO, lambda_H2O, lambda_HCN]

    return data

def loadCDTable():
    file = np.loadtxt(pathCDTable)
    nw = int(file[:, 0][-1])
    nx = int(file[:, 1][-1])
    ny = int(file[:, 2][-1])
    nz = int(file[:, 3][-1])

    data = np.zeros(shape=(nw, nx, ny, nz, 5))
    for line in file:
        i, j, k, l, r_B, vth, m, a, N = line
        data[int(i) - 1][int(j) - 1][int(k) - 1][int(l) - 1] = [r_B, vth, m, a, N]

    return data

def loadDump():
    run = "/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline/files"
    saveloc = "/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output"
    factor = 3
    bound = 30
    userSettingsDictionary = {}
    userSettingsDictionary["prefix"] = "wind"
    userSettingsDictionary[
        "data_location"] = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline"
    userSettingsDictionary[
        "pipeline_output_location"] = r"/Users/diond/Documents/PyCharmProjects/University/2nd Master/Master Thesis/pipeline_output"

    print("load data")
    setup, dumpData, sinkData, outerData = ld.LoadData_cgs(run, run, factor, bound, userSettingsDictionary)

    return dumpData

def findLambda(fileCool, dumpData, line, radCoolTable, CDTable):
    n0, m, v = fitting(dumpData)

    r_B     = None
    T       = None
    vth_CO  = None
    vth_H2O = None
    vth_NCN = None

    if line == "x":
        r_B = np.abs(fileCool[:, 0])
        T   = fileCool[:, 6]

    else:
        r_B = np.abs(fileCool[:,1])
        T = np.abs(fileCool[:, 7])

    vth_CO  = np.sqrt(2. * kB * T / mass_molecules[0]) / v
    vth_CO[vth_CO >= 2.] = 2.
    vth_H2O = np.sqrt(2. * kB * T / mass_molecules[1]) / v
    vth_H2O[vth_H2O >= 2.] = 2.
    vth_HCN = np.sqrt(2. * kB * T / mass_molecules[2]) / v
    vth_HCN[vth_HCN >= 2.] = 2.


    N_CO  = n0 * abundance_molecules[0] * columnDensity(CDTable, r_B,  vth_CO, m, 35, 101, 5, 7)
    N_H2O = n0 * abundance_molecules[1] * columnDensity(CDTable, r_B, vth_H2O, m, 35, 101, 5, 7)
    N_HCN = n0 * abundance_molecules[2] * columnDensity(CDTable, r_B, vth_HCN, m, 35, 101, 5, 7)

    nH = np.zeros(len(r_B))
    nH[r_B < r_orbit_comp * AU] = n0
    nH[r_B >= r_orbit_comp * AU] = n0 * (r_orbit_comp * AU/r_B[r_B >= r_orbit_comp * AU])**m
    nH = np.log10(nH)
    lambda_net  = np.ones(len(N_CO)) * -999.
    for i in range(len(N_CO)):
        lambda_CO = -999.
        lambda_H2O = -999.
        lambda_HCN = -999.

        if N_CO[i] != 0.:
            lambda_CO  = coolingRate(radCoolTable, np.log10(T[i]), np.log10(nH[i]), np.log10(N_CO[i]), 39, 39, 39, "CO")

        if N_H2O[i] != 0.:
            lambda_H2O = coolingRate(radCoolTable, np.log10(T[i]), np.log10(nH[i]), np.log10(N_H2O[i]), 39, 39, 39, "H2O")

        if N_HCN[i] != 0.:
            lambda_HCN = coolingRate(radCoolTable, np.log10(T[i]), np.log10(nH[i]), np.log10(N_HCN[i]), 39, 39, 39, "HCN")

        lambda_net[i] = np.log10( 10.** lambda_CO + 10.**lambda_H2O + 10.**lambda_HCN )

    return lambda_net

def fitting(dumpData):
    # Density Profile
    x = dumpData["position"][:, 0]
    y = dumpData["position"][:, 1]
    z = dumpData["position"][:, 2]
    r = np.sqrt( x**2. + y**2. + z**2.)

    n = dumpData["rho"] / (mu * 1.67E-24)
    r_fit_noCompZone = r[r_orbit_comp * AU < r]
    r_fit_noCompZone = np.log(r_fit_noCompZone)
    rho_fit_noCompZone = n[r_orbit_comp * AU < r]
    rho_fit_noCompZone = np.log(rho_fit_noCompZone)
    n_fit_noCompZone = len(r_fit_noCompZone)

    r_fit_compZone = r[r_orbit_comp * AU >= r]
    r_fit_compZone = np.log(r_fit_compZone)
    rho_fit_compZone = n[r_orbit_comp * AU >= r]
    n_fit_compZone = len(r_fit_compZone)

    # ln y = ln a + b ln r
    b_fit = -(n_fit_noCompZone * np.sum(r_fit_noCompZone * rho_fit_noCompZone) - np.sum(r_fit_noCompZone) * np.sum(
        rho_fit_noCompZone)) / (n_fit_noCompZone * np.sum(r_fit_noCompZone ** 2.) - np.sum(r_fit_noCompZone) ** 2.)
    a_fit = np.mean(rho_fit_compZone)

    # Velocity fit
    vx = dumpData["vx"]
    vy = dumpData["vy"]
    vz = dumpData["vz"]
    v  = np.sqrt(vx**2. + vy**2. + vz**2.)[:-2]

    v_fit = np.mean(v[np.logical_and(1.5 * r_orbit_comp * AU  <= r, r <= 0.9 * xmax * AU)])

    return [a_fit, b_fit, v_fit]

def findLowerIndicesCDTable(CDTable, r_B, vth, m):
    i = 0
    j = 0
    k = 0
    l = 0

    vth0 = CDTable[0][0][0][0][1]
    m0   = CDTable[0][0][0][0][2]
    a0   = CDTable[0][0][0][0][3]

    dv  = 0.02
    dm  = CDTable[0][0][1][0][2] - m0

    i = np.zeros(len(r_B), dtype=int)
    j = np.floor(vth / dv)
    j[np.logical_and(vth >= 0.001, vth < 0.02)] = 0
    j[np.logical_and(vth >= 1.999, vth < 2.)] = 100
    j[vth == 2.] = 101
    k = np.floor(  (m - m0) / dm)
    l = int(r_orbit_comp - a0)

    for h in range(len(r_B)):
        index = np.argmax(np.log10(r_B[h] / AU) < CDTable[:, 0, 0, l, 0]) - 1
        if (np.log10(r_B[h] / AU) > CDTable[0, 0, 0, l, 0]) and (np.log10(r_B[h] / AU) < CDTable[-1, 0, 0, l, 0]):
            i[h] = index

        else:
            i[h] = -1

    return [i, j, k, l]

def columnDensity(CDTable, r_B, vth, m, imax, jmax, kmax, lmax):
    i0_array, j0_array, k0, l0 = findLowerIndicesCDTable(CDTable, r_B, vth, m)

    i1 = 0
    j1 = 0
    k1 = 0

    N_interpol = np.zeros(len(i0_array))

    for p in range(len(N_interpol)):
        i0 = int(i0_array[p])
        j0 = int(j0_array[p])
        k0 = int(k0)
        l0 = int(l0)

        if (i0 >= 0) and (j0 >= 0) and (k0 >= 0):
            if i0 >= imax: i1 = i0
            else: i1 = i0 + 1
            if j0 >= jmax: j1 = j0
            else: j1 = j0 + 1
            if k0 >= kmax: k1 = k0
            else: k1 = k0 + 1

            r_B0  = CDTable[i0][j0][k0][l0][0]
            vth0 = CDTable[i0][j0][k0][l0][1]
            m0   = CDTable[i0][j0][k0][l0][2]

            r_B1  = CDTable[i1][j1][k1][l0][0]
            vth1 = CDTable[i1][j1][k1][l0][1]
            m1   = CDTable[i1][j1][k1][l0][2]

            x_d = 0.
            y_d = 0.
            z_d = 0.

            if i1 - i0 != 0: x_d = (np.log10(r_B[p]/AU) - r_B0) / (r_B1 - r_B0)
            if j1 - j0 != 0: y_d = (vth[p] - vth0) / (vth1 - vth0)
            if k1 - k0 != 0: z_d = (m - m0) / (m1 - m0)

            data = CDTable[i0: i1 + 1, j0: j1 + 1, k0: k1 + 1, l0, 4]
            N_interpol[p] = interpolate3D(x_d, y_d, z_d, data)

    return N_interpol
main()