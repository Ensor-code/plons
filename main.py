import radialStructure_pipeline as rs
import SlicePlots_pipeline as sl
import CMF_meanRho_pipeline as cmf
import orbitalEvolution_pipeline as ov

runNumber = 21
loc = '/home/user/Documents/phantom pipeline/'

ov.orbEv(runNumber,loc)
rs.radialStructPlots(runNumber,loc)
sl.SlicePlots(runNumber,loc)
cmf.CMF_meanRho(runNumber,loc)
