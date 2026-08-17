import sys as sys
import os as os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import pandas as pd
import multiprocessing as mp
from multiprocessing import freeze_support

import time

from . import SoFiA_Driver as SD
from . import SetFileLocations as SFL
from . import GalaxyFitParameters as GFP
from . import RunWRKP as RW
from . import ReadWRKPFit as RWF
from . import RunInitialFitDCP as RIF
from . import MakeBootstrapSample as MBS
from . import Bootstrap_Error_Analysis as BEA
from . import Bootstrap_Outputs as BO
from . import GeometryCorrection as GC
from . import BootstrapModelPlot as BMP
from . import ExtendedSDProfileCalc as ESDPC
from . import ExtractScalingParams as ESP
from . import BootstrapBoxPlot as BBP



def GalaxyFit():
    print("This is the WRKP Fit Driver Script")
    
    #   Get the various program and file locations
    GeneralDict,KeyRTParams=SFL.GetFileAndFolderLocationAndNames()
    #   Figure out the input parameter file
    ParamFile=GFP.GetRuntimeArguments()
    #   Resolve the config file's own absolute path BEFORE changing directory,
    #   then chdir into it. Every relative path inside a config (CubeName,
    #   TargFolder, etc.) is authored relative to the config file's own
    #   location -- this lets the driver be invoked from any directory instead
    #   of requiring the caller to already be sitting in that folder.
    ParamFile=os.path.abspath(ParamFile)
    os.chdir(os.path.dirname(ParamFile))
    #   Read in the runtime parameters
    RTParams=GFP.ReadParamFile(ParamFile)
    #   Overwrite possible default parameters with the runtime parameters
    GeneralDict,GalaxyDict=GFP.OverwriteDefaults(GeneralDict,RTParams,KeyRTParams)
    #   Make sure all the variables are formatted correctly for use
    GFP.CheckParamTypes(GalaxyDict,KeyRTParams)
    #   Now use the default WRKP file locations that have already been specified and load them in
    GeneralDict=RW.LoadDefaultWRKPFiles(GeneralDict)
    
    #   Assuming that both the general dictionaries and galaxy dictionaries have the correct entries, we can run set up the WRKP input files
    #       Before beginning the analysis we'll copy the galaxy dictionary names to variables that will be used for a specific run.  This allows the running commands to be more general.
    GalaxyDict=GFP.SetCurrentRunVariables(GalaxyDict)
    
    #       before running the fitter we'll make a folder for the galaxy.  Currently the Pipeline will make this as well, but this helps to avoid errors if the fitter fails
    TargFolder=GalaxyDict['TargFolderU']+GalaxyDict['ObjNameU']+"/"
    os.makedirs(TargFolder,exist_ok=True)

    #   The initial/anchor fit: UseDCP==1 runs it entirely through the JS/DCP
    #       pipeline (RunInitialFitDCP.py -- no Fortran involvement at all,
    #       cube/beam geometry derived directly from the raw FITS bytes);
    #       UseDCP==0 (fortran-local) is completely unchanged below.
    start = time.time() #   Time how long the fit takes
    if GalaxyDict.get('UseDCP', 0) == 1:
        #   No SoFiA catalogue text file needed -- RunInitialFitDCP passes
        #       PA_EstimateU/Inc_EstimateU directly in the JSON payload
        #       instead of round-tripping them through a file only Fortran
        #       would read.
        GalaxyDict['BestFitModel']=RIF.RunInitialFit(GeneralDict,GalaxyDict)
    else:
        #       The first step is create a 'catalogue' file with the initial estimates for the inclination and position angle
        GalaxyDict['SoFiAShapeFile']=SD.WriteSoFiACatFileForWRKP(GalaxyDict)
        #   Now run WRKP
        GalaxyDict=RW.RunWRKP(GeneralDict,GalaxyDict,0)
        #   With WRKP now completed, the fit must be ingested
        #       It is possible that the initial fit will have failed.  In that case, we won't be able to read the output file, so we should end the program here
        try:
            GalaxyDict['BestFitModel']=RWF.ReadWRKPOutputFile(GeneralDict,GalaxyDict)
        except:
            GalaxyDict['BestFitModel']={'FITAchieved':False}
    end = time.time()
    InitialFitTime=(end-start)

    if GalaxyDict['BestFitModel']['FITAchieved']==False:
        print("No best fit model made")
        exit()

    #   Write out the total time of the fit
    TimeFile=GalaxyDict['TargFolderU']+GalaxyDict['ObjNameU']+"/"+GalaxyDict['ObjNameU']+"_FitTimeCheck.csv"
    #f=open(TimeFile,'w')
    #Str="Initial Fit Runtime (in seconds) = "+str(InitialFitTime)+"\n"
    #f.write(Str)
    #f.close()
    
    
    #   The next step is to estimate the uncertainties for the model parameters by running a number of bootstrap fits.
    #   Before we loop through all the bootstraps, we need set a few naming conventions
    GalaxyDict=BEA.SetErrorAnalysisVariables(GalaxyDict)
    #       It's also necessary to load in the SoFiA template file
    GeneralDict['SoFiATemplate']=SD.LoadSoFiATemplate(GeneralDict['SoFiATemplateFile'])

    #       Now the loop over all bootstraps can be started
    BootstrapModels=[None]*GalaxyDict['nBootstraps']
    #           Set the number of processors to use for bootstrapping
    nProcessors=GalaxyDict['nProcessors_Bootstraps']
    #           Split the bootstrap fits over some number of parameters
    
    #for i in range(GalaxyDict['nBootstraps']):
    #    BootstrapModels[i]=BootstrapRunStep(i,GeneralDict,GalaxyDict)
    #   Time the bootstrap loop as a whole -- directly comparable between the
    #       DCP and local-pool paths, regardless of which one runs.
    BootstrapLoopStart = time.time()
    if GalaxyDict.get('UseDCP', 0) == 1:
        from . import RunBootstrapsDCP as RBD
        BootstrapModels = RBD.RunBootstrapsDCP(GeneralDict, GalaxyDict)
    else:
        pool = mp.Pool(processes=nProcessors)
        BootstrapModels = pool.starmap(
            BootstrapRunStep,
            [(i, GeneralDict, GalaxyDict) for i in range(GalaxyDict['nBootstraps'])],
            chunksize=1)
        pool.close(); pool.join()
        #   Per-realization timing breakdown (see Bootstrap_Error_Analysis.
        #   GetBootstrapModel) -- pulled off each model dict and written
        #   separately so StoreBootstrappedModels_CSV/EstimateUncertainties
        #   FromBootstraps see the exact same model shape as before this was
        #   added. Same file name/schema RunBootstrapsDCP.py writes for the
        #   DCP path, so run_both.js can compare the two directly.
        BootstrapTimings = [m.pop('Timing', {'realizationIndex': i}) for i, m in enumerate(BootstrapModels)]
        BO.StoreBootstrapTimings_JSON(GalaxyDict, BootstrapTimings)
    BootstrapLoopTime = time.time() - BootstrapLoopStart

    #   Get the uncertaintites on each parameter
    BEA.EstimateUncertaintiesFromBootstraps(GeneralDict,GalaxyDict,BootstrapModels)
    #   Adjust the PA to reach the global PA
    GalaxyDict['BestFitModel']=GC.GetGlobalPositionAngle(GalaxyDict)
    #       With the uncertainties calculated we can calculate an extend Surface density profile
    GalaxyDict=ESDPC.CalcExtendedSDProfile(GalaxyDict)
    #   Extract RHI and VHI for the model
    GalaxyDict=ESP.CalcScalingParams(GalaxyDict,BootstrapModels)
    

    #   Output a text file with the errors and uncertainties
    BO.WriteBootstrappedFitOutputFile_Text(GalaxyDict)
    
    #   Make diagnostic plot from the bootstrapped model
    BMP.MakeBootstrapModelPlot(GalaxyDict,GeneralDict)
    #   Make a box plot
    BBP.BootstrapBoxPlot(GalaxyDict,GeneralDict,BootstrapModels)
    
    #   Save the PV diagrams to Fits files
    BO.SavePVDiagrams(GalaxyDict,GeneralDict)
    
    #   Save the bootstraps to a csv file
    BO.StoreBootstrappedModels_CSV(GalaxyDict,BootstrapModels)

    #   Keep track of the total runtime
    end = time.time()
    TotalFitTime=(end-start)
    f=open(TimeFile,'a')
    #Str="Total Fit Runtime (in seconds) = "+str(TotalFitTime)+"\n"
    #f.write(Str)
    #f.close()
    CheckPoints=["Initial Fit","Bootstrap Loop"]
    RunTimes=[InitialFitTime,BootstrapLoopTime]
    #   When using DCP, also break the bootstrap loop down into three
    #       checkpoints -- resample+SoFiA (now distributed, round 1), local
    #       prep (fixture-only refit, one native call per realization now
    #       instead of three), and the fit dispatch itself (round 2) -- all
    #       set on GalaxyDict by RunBootstrapsDCP.RunBootstrapsDCP.
    if 'DCP_LocalPrepTime' in GalaxyDict:
        if 'DCP_ResampleTime' in GalaxyDict:
            #   Round 1 (resample+SoFiA) and the local fixture-only-fit now
            #   run pipelined/overlapped per realization (RunBootstrapsDCP.
            #   RunBootstrapsDCP), not as two fully serial phases -- so
            #   "DCP Local Prep (native, overlaps above)" is native-call CPU
            #   time, NOT wall-clock-additive to the row above it or to
            #   Total; it's informational only. "DCP Resample+SoFiA+Local
            #   Prep" is the one that's actually additive toward Total.
            CheckPoints+=["DCP Resample+SoFiA+Local Prep"]
            RunTimes+=[GalaxyDict['DCP_ResampleTime']]
        CheckPoints+=["DCP Local Prep (native, overlaps above)","DCP Dispatch"]
        RunTimes+=[GalaxyDict['DCP_LocalPrepTime'],GalaxyDict['DCP_DispatchTime']]
    CheckPoints.append("Total")
    RunTimes.append(TotalFitTime)
    TimeDict={'CheckPoint':CheckPoints,'RunTime':RunTimes}
    DF=pd.DataFrame.from_dict(TimeDict)
    DF.to_csv(TimeFile,index=False)
    
    #   Once the bootstrap run is done, remove the bootstrap and SoFiA folders
    ClnCmd="rm -r "+GalaxyDict['BootstrapFolder']+" "+GalaxyDict['SoFiAFolder']
    os.system(ClnCmd)
   

def BootstrapRunStep(step,GeneralDict,GalaxyDict):
    BootstrapModel=BEA.GetBootstrapModel(GeneralDict,GalaxyDict,step)
    return BootstrapModel
