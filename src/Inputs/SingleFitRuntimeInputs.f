cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting general
c       inputs for fitting a single galaxy
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SingleFitInputMod

      use PipelineGlobals
      use FittingOptionsMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating a mock cube
      subroutine SingleFitIn(MaskSwitch)
      implicit none
      character(500) TopLevelInputFile
      integer MaskSwitch
      integer ios


      print*, "Getting the cube inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Pipeline Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the name of the folder containing everything
      read(10,*)
      read(10,'(A)') DataCubeFile
c           Get the name of the fitting options file
      read(10,*)
      read(10,'(A)') FittingOptionsFile
c           Get the mask switch
      read(10,*)
      read(10,*) MaskSwitch

      read(10,*)
      if(MaskSwitch .eq. 1) then
        read(10,'(A)') MaskName
      endif
c       See whether there is a catalogue file to load in
      read(10,*)
      read(10,*) PFlags%CatFlag
      read(10,*)
      if(PFlags%CatFlag .eq. 1) then
        read(10,'(A)') CatalogueFile
      endif



c       Get a value for the random seed
      read(10,*)
      read(10,*) idum
c       Read in the name of the main output folder to store different fits in
      read(10,*)
      read(10,'(A)') MainOutputFolder

c       Read in the name of the object itself to store different fits in
      read(10,*)
      read(10,'(A)') GalaxyDict%GalaxyName
c       Read in the cube noise (this is temporary until I implement a noise calculation
      read(10,*)
      read(10,*) GalaxyDict%RMS

c       Optionally read a fixture-only switch (stop after DumpFittingFixture,
c           skip the optimizer). Absent for older input files -- default to 0
c           (normal behaviour) if the line isn't there.
      FixtureOnlySwitch=0
      read(10,*,IOSTAT=ios)
      if (ios .eq. 0) then
        read(10,*,IOSTAT=ios) FixtureOnlySwitch
      endif

c       Optionally read a probe switch (evaluate the objective function at
c           specific parameter vectors instead of running the optimizer --
c           see PipelineGlobals.f). Same backward-compatible pattern.
      ProbeSwitch=0
      read(10,*,IOSTAT=ios)
      if (ios .eq. 0) then
        read(10,*,IOSTAT=ios) ProbeSwitch
      endif

c       Optionally read a trace switch (print a per-evaluation debug line --
c           call counter, idum, chi2, PA -- from TiltedRingModelComparison.
c           Investigative-only; see FullModelComparison.f). Same
c           backward-compatible pattern.
      TraceSwitch=0
      read(10,*,IOSTAT=ios)
      if (ios .eq. 0) then
        read(10,*,IOSTAT=ios) TraceSwitch
      endif

c       Optionally read a dump-fixture switch (1 = dump diskfit_fixture.json/
c           model_cube_bestfit.json for the DCP/JS legs to consume; see
c           PipelineGlobals.f). Same backward-compatible pattern.
      DumpFixtureSwitch=0
      read(10,*,IOSTAT=ios)
      if (ios .eq. 0) then
        read(10,*,IOSTAT=ios) DumpFixtureSwitch
      endif

      close(10)


      call FittingOptionsIn()

      GalaxyDict%DataCubeFile=trim(DataCubeFile)

      return
      end subroutine
cccccccc

      end module
