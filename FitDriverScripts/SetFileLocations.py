import os as os

def GetFileAndFolderLocationAndNames():
    #   Get the current directory for future writing
    CurrDir=os.getcwd()
    #   Get the WRKP directory
    WRKPDir= os.path.dirname(os.path.realpath(__file__))
    #   Adjust the WRKP directory path slightly to point one directory up due to this file being located in WRKP/FitDriverScripts
    WRKPDir=WRKPDir.rsplit('/', 1)[0]
    #   Set the SingleGalaxyFitter and BootstrapSampler paths
    FitterExecPath=WRKPDir+"/Programs/SingleGalaxyFitter"
    BootstrapExecPath=WRKPDir+"/Programs/BootStrapSampler"
    #   Set the SoFiA path to the third_party install
    SoFiAExecPath=WRKPDir+"/third_party/SoFiA-2-master_2_5_1/sofia"
    #   Set the default SoFiA template file
    SoFiATemplateFile=WRKPDir+"/third_party/SoFiA-2-master_2_5_1/template_par_file.par"
    #   Set the default WRKP Main file
    WRKP_GeneralMainIn=WRKPDir+"/Inputs/SingleFitInput_Base.in"
    #   Set the default WRKP Options file
    WRKP_GeneralOptionsIn=WRKPDir+"/Inputs/SingleGalaxyTestFittingOptions_Base.txt"
    #   DCP bootstrap distribution: everything DCP-related lives under a fixed
    #       location, js/app/ (the JS/DCP port's browser+CLI application
    #       layer -- see js/ARCHITECTURE.md for the full layout). These never
    #       vary per-galaxy, so they're computed here rather than requiring
    #       every config to repeat an absolute, machine-specific path.
    #       FixtureFile is where FullSingleGalaxyFit moves the initial fit's
    #       fixture dump when UseDCP=1, so it has to stay derived from the
    #       same DCPJobDir.
    DCPDir=WRKPDir+"/js/app"
    DCPJobDir=DCPDir+"/DCPjobData"
    #   One DCP job, one slice per bootstrap realization -- each slice runs
    #       the entire resample+SoFiA+InitialAnalysis+fit pipeline (see
    #       RunBootstrapsDCP.py's module docstring). Replaces the prior
    #       two-launcher pair (bootstrap-resample-launcher.js +
    #       bootstrap-fit-launcher.js), kept in the tree only as historical
    #       reference for the old two-round-trip design.
    DCPRealizationLauncher=DCPDir+"/bootstrap-realization-launcher.js"
    FixtureFile=DCPJobDir+"/diskfit_fixture.json"
    #   The initial fit's best-fit model cube, dumped by the Fortran fitter
    #       as a side artifact (dense flux array, same fixed-location
    #       reasoning as FixtureFile) -- round 1's static payload needs this
    #       to compute the resampling residual (observed - model).
    ModelCubeFile=DCPJobDir+"/model_cube_bestfit.json"
    #   Opt-in resampling RNG seed for run_both/debug comparisons (see
    #   RunBootstrapsDCP.BuildResamplePayload, MakeBootstrapSample.
    #   WriteBootstrapFile). 0 means "no debug seed" -- every existing config
    #   that doesn't set this gets exactly today's behaviour. Not in
    #   KeyRTParamsList (optional, not required): a config only needs to set
    #   BootstrapSeed=<value> to override this default via OverwriteDefaults'
    #   own name-matching mechanism, same as any other GeneralDict default.
    BootstrapSeed=0
    #   Set the general file dictionary as the set of local variables thus far
    FileDict=locals()
    #   Set the list of runtime parameters that must be set in the RT Params  file
    #       It is necessary to list each variable with it's type i.e. [Variable Name, Type] as the types will be checked in GalaxyFitParamters.CheckParamTypes
    KeyRTParamsList=[['CubeName',str]
    ,['MaskName',str]
    ,['PA_Estimate',float]
    ,['Inc_Estimate',float]
    ,['nBootstraps',int]
    ,['ObjName',str]
    ,['TargFolder',str]
    ,['nProcessors_Bootstraps',int]
    ,['UseDCP',int]]
    return FileDict,KeyRTParamsList
    
