import os, os.path
import numpy as np
import RadialStructurePlots1D       as rs
import SlicePlots2D                 as sl
import CMF_meanRho                  as cmf
import OrbitalEvolution             as ov
import LoadDataPHANTOM              as ld
import TerminalVelocity             as tmv
import Tubes                        as tb
import userSettings                 as us
import LoadDump                     as lodu

from mpi4py import MPI
comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()
nproc = comm.Get_size()


def main():
    # Initialise user settings or load them
    userSettingsFilePath = os.path.join( os.getcwd(), "userSettings.txt")
    if not os.path.isfile(userSettingsFilePath) or os.stat(userSettingsFilePath).st_size == 0: us.create(userSettingsFilePath)

    run = "cooling/local"
    userSettingsDictionary = us.load(userSettingsFilePath)
    prefix = userSettingsDictionary["prefix"]
    loc = userSettingsDictionary["data_location"]
    outputloc = userSettingsDictionary["pipeline_output_location"]
    saveloc = os.path.join(outputloc, run)
    factor = 3
    bound = None
    try:
        os.makedirs(os.path.join(saveloc, 'animation'))
    except OSError:
        pass

    print("Load data")
    [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound, userSettingsDictionary)
    pathData = os.path.join(loc, run)
    listFiles = os.listdir(pathData)
    listFiles = list(filter(lambda file: file != ".DS_Store", listFiles))
    listFiles = np.sort(listFiles)

    fullDumpIndices = lodu.findAllFullDumpIndices(userSettingsDictionary, setup, pathData)
    for i in range(len(fullDumpIndices[rank::nproc])):
        i_proc = rank + i * nproc
        index = fullDumpIndices[i_proc]
        print("\n")
        print("-------------------------------------------")
        print("Processing dump file {0:03d} out of {1:03d}".format(i, len(fullDumpIndices)))
        print("-------------------------------------------")
        print("\n")
        [setup, dumpData, sinkData, outerData] = ld.LoadData_cgs(run, loc, factor, bound,
                                                                 userSettingsDictionary, number = index)
        sl.SlicePlots(run, saveloc, dumpData, setup, number = index)

main()