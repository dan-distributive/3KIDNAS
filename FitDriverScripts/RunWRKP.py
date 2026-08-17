import sys as sys
import os as os
import copy as copy

def LoadDefaultWRKPFiles(GeneralDict):
    print("Loading in default WRKP control files")
    
    #   Set the default file names using the general dictionary
    Infile=GeneralDict['WRKP_GeneralMainIn']
    OptionsFile=GeneralDict['WRKP_GeneralOptionsIn']
    
    #   Open up the default main input file and read in all the lines and store them
    MainIn=open(Infile,"r")
    Main_In_Lines=MainIn.readlines()
    MainIn.close()
    #   Do the same for the options file
    OptionsIn=open(OptionsFile,"r")
    WRKP_Options_Lines=OptionsIn.readlines()
    OptionsIn.close()
    #   Store all the default line information into the general dictionary and return it
    GeneralDict['MainWRKPInputLines']=Main_In_Lines
    GeneralDict['MainWRKPOptionsLines']=WRKP_Options_Lines
    return GeneralDict
    
def NameWRKPRunFiles(GalaxyDict,BootstrapSwitch):
    if BootstrapSwitch==0:
        GalaxyDict['MainInputFile']=GalaxyDict['ObjNameU']+"_WRKP_Input.txt"
        GalaxyDict['FittingOptionsFile']=GalaxyDict['ObjNameU']+"_WRKP_Options.txt"
        GalaxyDict['LogFile']=GalaxyDict['ObjNameU']+"_FittingLog.txt"
    
    return GalaxyDict

def WriteWRKPMainFile(WorkingMainLines,GalaxyDict,FixtureOnlySwitch=0,ProbeSwitch=0,TraceSwitch=0,DumpFixtureSwitch=0):

    #   Set the name of the cube to be fit
    WorkingMainLines[1]=GalaxyDict['CubeNameU']+"\n"
    #   Set the name of the fitting options file
    WorkingMainLines[3]=GalaxyDict['FittingOptionsFile']+"\n"
    #   Set the name of the object
    WorkingMainLines[15]=GalaxyDict['ObjNameU']+"\n"
    #   Set the name of the folder that will contain the outputs
    WorkingMainLines[13]=GalaxyDict['TargFolderU']+"\n"

    #   Set up the portion of the file to deal with SoFiA shape estimates
    WorkingMainLines[8]=str(1)+"\n"
    if len(WorkingMainLines) <= 18:
        WorkingMainLines.insert(10,GalaxyDict['SoFiAShapeFile']+"\n")
    else:
        WorkingMainLines[10]=GalaxyDict['SoFiAShapeFile']+"\n"

    #   Append the fixture-only switch (1 = stop right after DumpFittingFixture,
    #       skipping the optimizer). Optional/trailing so older callers and the
    #       Fortran reader (which defaults to 0 if the line is absent) are unaffected.
    WorkingMainLines.append("#\tFixture-only switch (1 = stop after fixture dump)\n")
    WorkingMainLines.append(str(FixtureOnlySwitch)+"\n")

    #   Append the probe switch (1 = evaluate the objective function at
    #       parameter vectors read from probe_params.txt instead of running
    #       the optimizer -- see PipelineGlobals.f). Same optional/trailing
    #       backward-compatible pattern as FixtureOnlySwitch.
    WorkingMainLines.append("#\tProbe switch (1 = evaluate probe_params.txt, skip optimizer)\n")
    WorkingMainLines.append(str(ProbeSwitch)+"\n")

    #   Append the trace switch (1 = print one line per TiltedRingModelComparison
    #       call -- call counter, idum, chi2, PA -- see PipelineGlobals.f).
    #       Same optional/trailing backward-compatible pattern.
    WorkingMainLines.append("#\tTrace switch (1 = print per-call idum/chi2/PA trace)\n")
    WorkingMainLines.append(str(TraceSwitch)+"\n")

    #   Append the dump-fixture switch (1 = dump diskfit_fixture.json /
    #       model_cube_bestfit.json -- the state the DCP/JS legs consume).
    #       Only ever True for the one initial fit on a UseDCP=1 config (see
    #       RunWRKP's own computation of this) -- everything else, including
    #       fortran-local's own per-realization bootstrap fits, leaves this 0
    #       so it stops uselessly (and racily, across nProcessors parallel
    #       realizations) re-dumping both files on every single fit call.
    #       Same optional/trailing backward-compatible pattern.
    WorkingMainLines.append("#\tDump-fixture switch (1 = dump diskfit_fixture.json/model_cube_bestfit.json)\n")
    WorkingMainLines.append(str(DumpFixtureSwitch)+"\n")

    #   Now that we have the main input file contents written, we can write out the file
    mFile=open(GalaxyDict['MainInputFile'],'w')
    for x in WorkingMainLines:
        mFile.write(x)
    mFile.close()
    
def WriteWRKPOptionsFile(WorkingOptionsLines,GalaxyDict,BootstrapSwitch):
    
    #   Set the name of the cube to be used for initial estimates beyond the inclination and position angle
    WorkingOptionsLines[9]=GalaxyDict['CubeNameU']+"\n"
    #   And give it the name of the mask file
    WorkingOptionsLines[10]=GalaxyDict['MaskNameU']+"\n"
    #   If we are doing a bootstrap, we want to set the number of rings to the original number
    if BootstrapSwitch ==1:
        nRTarg=len(GalaxyDict['BestFitModel']['R'])
        #print("Bootstrap WRKP options", BootstrapSwitch)
        WorkingOptionsLines[30]=str(nRTarg)+"\n"
        #print(nRTarg,WorkingOptionsLines[28])
        RArr=GalaxyDict['BestFitModel']['R']
        RStr=[None]*len(RArr)
        for i in range(len(RArr)):
            RStr[i]=str(RArr[i])
        RStr=",".join(RStr)
        RStr+="\n"
        BootStrLabel="#  If using a grid, supply the radius array in arcsecs\n"
        WorkingOptionsLines.insert(31,BootStrLabel)
        WorkingOptionsLines.insert(32,RStr)
    
    #   Now write the contents to the temporary options file
    oFile=open(GalaxyDict['FittingOptionsFile'],'w')
    for x in WorkingOptionsLines:
        oFile.write(x)
    oFile.close()



def RunWRKP(GeneralDict,GalaxyDict,BSSwitch):
    #   Only announced for the one-off initial fit -- BSSwitch=1 fires once
    #   per bootstrap realization (up to nBootstraps times) and the per-
    #   realization summary line in Bootstrap_Error_Analysis.GetBootstrapModel
    #   already covers it.
    if BSSwitch == 0:
        print("Running WRKP")

    #   First copy the general Main lines to a local set
    WorkingMainLines=copy.copy(GeneralDict['MainWRKPInputLines'])
    #   Next name the temporary WRKP Files for the run
    GalaxyDict=NameWRKPRunFiles(GalaxyDict,0)
    #   Now write up the main input file
    #   WRKP_TRACE_DEBUG=1 (diagnostic only, mirrors the JS side's own
    #   TRACE_DEBUG env var convention -- FullModelComparison.js) enables
    #   Fortran's per-eval TRACE print for bisecting where JS and Fortran's
    #   optimizer trajectories first diverge, for the SAME seeded bootstrap
    #   fit. Never on by default.
    TraceSwitch = 1 if os.environ.get('WRKP_TRACE_DEBUG') == '1' else 0
    #   Only the ONE initial fit (BSSwitch=0) needs diskfit_fixture.json/
    #   model_cube_bestfit.json, and only when a DCP/JS leg will actually
    #   consume them (UseDCP=1) -- BSSwitch=1 is exclusively fortran-local's
    #   own per-realization bootstrap fits (UseDCP=1 configs never call
    #   RunWRKP with BSSwitch=1 at all; they go through RunBootstrapsDCP
    #   instead), which used to unconditionally re-dump both files on every
    #   realization for no reason (nothing on that path reads them back).
    DumpFixtureSwitch = 1 if (BSSwitch == 0 and GalaxyDict.get('UseDCP', 0) == 1) else 0
    WriteWRKPMainFile(WorkingMainLines,GalaxyDict,TraceSwitch=TraceSwitch,DumpFixtureSwitch=DumpFixtureSwitch)
    #   Once this is done, we can get rid of the WorkingMainLines
    del WorkingMainLines
    #   Now we'll do a similar process for the fitting options:
    #       First copy general options lines to a local variable
    WorkingOptionsLines=copy.copy(GeneralDict['MainWRKPOptionsLines'])
    #       And next we'll write the specific fitting options for this run
    WriteWRKPOptionsFile(WorkingOptionsLines,GalaxyDict,BSSwitch)
    #   As before, we'll clean by getting rid of the working set of options lines
    del WorkingOptionsLines
    
    #   Now we can run WRKP
    #       Start with making the runtime command. The Fortran binary itself
    #       prints a lot (initial-vector dumps, per-guess chi2, param-limit
    #       checks, etc.) -- fine to see once for the initial fit (BSSwitch=0),
    #       but multiplied by nBootstraps realizations (and further
    #       interleaved across nProcessors parallel workers) it floods the
    #       terminal. For a bootstrap realization, redirect it to GalaxyDict
    #       ['LogFile'] instead (already unique per realization, since
    #       SetBootstrapVariableValues sets ObjNameU to "..._Bootstrap_<Step>"
    #       before this runs) -- nothing is lost, it's just on disk instead of
    #       on screen if a specific realization needs a closer look later.
    WRKPCmd=GeneralDict['FitterExecPath']+" "+GalaxyDict['MainInputFile']
    if BSSwitch == 1:
        WRKPCmd += " > "+GalaxyDict['LogFile']+" 2>&1"
    #       And now run it
    os.system(WRKPCmd)
    #   And clean up the input files
    ClnCmd="rm "+GalaxyDict['MainInputFile']+" " + GalaxyDict['FittingOptionsFile'] + " " +GalaxyDict['SoFiAShapeFile']
    os.system(ClnCmd)
    #   Finally return the galaxy dictionary
    return GalaxyDict
