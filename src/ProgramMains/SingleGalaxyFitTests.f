ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This code is meant to test routines for fitting a single galaxy
c       using the routines that are built into the WALLABY Kinematic Pipeline
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program SingleGalaxyTester


      use PipelineGlobals
      use SingleFitInputMod
      use DataCubeInputMod
      use PreGalaxyAnalysisMod
      use PostGalaxyAnalysisMod
      use DataCubeOutputsMod

      use FitOutputMod

      use MakeMomentMapsMod

      use CommonConsts

      use ParameterVectorToTiltedRingMod
      implicit none

      integer MaskSwitch,DataSwitch
      integer i
      character(100) TestMapName

      integer j,k

      real SDConv

      Type(Beam2D) TempBeam

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "Single Galaxy Fit Tests"

      Version=1
c       Get the runtime inputs
      call SingleFitIn(MaskSwitch)
c       Load the cube
      DataSwitch=0
      call ReadFullDataCube(ObservedDC,ObservedBeam
     &          ,DataCubeFile,DataSwitch)
      print*, "Total DC brightness", sum(ObservedDC%Flux)
c       Allocate the beam
      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
c           The actual data cube must be in units of Jy/pixel
      call DCBrightnessConversion(ObservedDC,ObservedBeam)
      print*, "Total DC brightness after unit conversions"
     &              , sum(ObservedDC%Flux)


c       If using a mask, load it in
      if(MaskSwitch .eq. 1) then
c        print*, "Reading Mask File"
        call ReadFullDataCube(DataCubeMask,TempBeam
     &          ,MaskName,MaskSwitch)
c        print*, "Done mask input"
        GalaxyDict%MaskFile=MaskName
c               Set the cube flux to either 1 or zero
        do i=0,DataCubeMask%DH%nPixels(0)-1
            do j=0,DataCubeMask%DH%nPixels(1)-1
                do k=0,DataCubeMask%DH%nChannels-1
                    if(DataCubeMask%Flux(i,j,k) .ge. 1.e-7) then
                        DataCubeMask%Flux(i,j,k)=1.
                    else
                        DataCubeMask%Flux(i,j,k)=0.
                    endif
                enddo
            enddo
        enddo
      elseif(MaskSwitch .eq. 3) then
c           Otherwise set the mask equal to the observed cube
        DataCubeMask=ObservedDC
        DataCubeMask%Flux=1.
        GalaxyDict%MaskFile=DataCubeFile
      endif

      GalaxyFit=>GalaxyFit_Simple
      OutputFit=>OutputBestFit_Simple
c      GalaxyFit=>BBaroloFit

c       Set up a catalogue object to match the general pipeline structure
      call MakeCatalogue()

c       Call the pre-galaxy analysis for the object
      call PreGalaxyAnalysis(SCatLocal%Objects(0))
      if (DumpFixtureSwitch .eq. 1) then
        call DumpFittingFixture('diskfit_fixture.json')
      endif

      if (FixtureOnlySwitch .eq. 1) then
        print*, "FixtureOnlySwitch set -- stopping after fixture dump"
        stop
      endif

c      print*, "initial vector", PVIni%Param(0:PVIni%nParams-1)
c      Call Galaxy fit for the test cube
      call GalaxyFit(SCatLocal%Objects(0))
      if (DumpFixtureSwitch .eq. 1) then
        call DumpBestFitModelCube('model_cube_bestfit.json')
      endif
c   Finally do any extra analysis to clean up the results
      call PostGalaxyAnalysis()
c       Output the single galaxy test fit
      call system("mkdir "//trim(MainOutputFolder))
      call OutputFit(SCatLocal%Objects(0))

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccc
      subroutine MakeCatalogue()
      use PipelineGlobals
      use SoFiACatalogueMod
      implicit none

c      print*, "Make catalogue entry"

      SCatLocal%nObjects=1
      call SoFiACatalogueAllocation(SCatLocal)

      call CatalogueEntryDefinition(SCatLocal%Objects(0))

      SCatLocal%SourceName="MCGSuite"

      if(PFlags%CatFlag .eq. 1) then
        open(10,file=trim(CatalogueFile),status='old')
        read(10,*)
        read(10,*) SCatLocal%Objects(0)%EllipseInc
        read(10,*)
        read(10,*) SCatLocal%Objects(0)%EllipsePA
        close(10)
      endif
    

      return
      end subroutine
cccccc

ccccc
      subroutine CatalogueEntryDefinition(CatObj)
      use SoFiACatalogueMod
      use CommonConsts
      use PipelineGlobals
      implicit none
      Type(CatalogueItem), INTENT(INOUT) :: CatObj

c       Set the catalogue object ID
      CatObj%ObID=0
c       Set the catalogue object name using the input name
      CatObj%ObjName=trim(GalaxyDict%GalaxyName)
c       Set the catalogue object RMS using the input RMS
      CatObj%RMS=GalaxyDict%RMS


      return
      end subroutine
cccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   DumpFittingFixture
c
c   Dumps all state needed to run GalaxyFit_Simple in JS.
c   Call after PreGalaxyAnalysis() and before GalaxyFit().
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine DumpFittingFixture(filename)
      use PipelineGlobals
      use DataCubeMod
      use CommonConsts
      implicit none

      character(*), INTENT(IN) :: filename
      integer i, j, l, ii, jj, kk
      integer, parameter :: unit=42
      character(5) BoolStr

      print*, "Dumping fitting fixture to: ", trim(filename)

      open(unit, file=trim(filename),
     &     status='replace', action='write')

      write(unit,'(A)') '{'

c     PVIni
      write(unit,'(A,I6,A)') '  "nParams": ',
     &     PVIni%nParams, ','
      write(unit,'(A)') '  "param": ['
      do i=0, PVIni%nParams-1
        if(i .lt. PVIni%nParams-1) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         PVIni%Param(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         PVIni%Param(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramLowerLims": ['
      do i=0, PVIni%nParams-1
        if(i .lt. PVIni%nParams-1) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         PVIni%ParamLowerLims(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         PVIni%ParamLowerLims(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramUpperLims": ['
      do i=0, PVIni%nParams-1
        if(i .lt. PVIni%nParams-1) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         PVIni%ParamUpperLims(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         PVIni%ParamUpperLims(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramRange": ['
      do i=0, PVIni%nParams-1
        if(i .lt. PVIni%nParams-1) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         PVIni%ParamRange(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         PVIni%ParamRange(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "cyclicSwitch": ['
      do i=0, PVIni%nParams-1
        if(i .lt. PVIni%nParams-1) then
          write(unit,'(A,I2,A)') '    ',
     &         PVIni%CyclicSwitch(i), ','
        else
          write(unit,'(A,I2)') '    ',
     &         PVIni%CyclicSwitch(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

c     TR_FittingOptions
      write(unit,'(A,I6,A)') '  "nRings": ',
     &     TR_FittingOptions%nRings, ','
      write(unit,'(A,I6,A)') '  "nFittedParamsTotal": ',
     &     TR_FittingOptions%nFittedParamsTotal, ','

      write(unit,'(A)') '  "constParams": ['
      do i=0, 12
        if (TR_FittingOptions%ConstParams(i)) then
          BoolStr='true'
        else
          BoolStr='false'
        endif
        if(i .lt. 12) then
          write(unit,'(A,A,A)') '    ',
     &         trim(BoolStr), ','
        else
          write(unit,'(A,A)') '    ',
     &         trim(BoolStr)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "fixedParams": ['
      do i=0, 12
        if (TR_FittingOptions%FixedParams(i)) then
          BoolStr='true'
        else
          BoolStr='false'
        endif
        if(i .lt. 12) then
          write(unit,'(A,A,A)') '    ',
     &         trim(BoolStr), ','
        else
          write(unit,'(A,A)') '    ',
     &         trim(BoolStr)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramLowerLims13": ['
      do i=0, 12
        if(i .lt. 12) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         TR_FittingOptions%ParamLowerLims(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         TR_FittingOptions%ParamLowerLims(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramUpperLims13": ['
      do i=0, 12
        if(i .lt. 12) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         TR_FittingOptions%ParamUpperLims(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         TR_FittingOptions%ParamUpperLims(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "paramRange13": ['
      do i=0, 12
        if(i .lt. 12) then
          write(unit,'(A,ES20.10E3,A)') '    ',
     &         TR_FittingOptions%ParamRange(i), ','
        else
          write(unit,'(A,ES20.10E3)') '    ',
     &         TR_FittingOptions%ParamRange(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

      write(unit,'(A)') '  "cyclicSwitch13": ['
      do i=0, 12
        if(i .lt. 12) then
          write(unit,'(A,I2,A)') '    ',
     &         TR_FittingOptions%CyclicSwitch(i), ','
        else
          write(unit,'(A,I2)') '    ',
     &         TR_FittingOptions%CyclicSwitch(i)
        endif
      enddo
      write(unit,'(A)') '  ],'

c     RadialProfiles
      write(unit,'(A)') '  "radialProfiles": ['
      do i=0, TR_FittingOptions%nRings-1
        write(unit,'(A)') '    {'
        write(unit,'(A,ES20.10E3,A)') '      "rmid": ',
     &       TR_FittingOptions%RadialProfiles(i)%Rmid, ','
        write(unit,'(A,ES20.10E3,A)') '      "rwidth": ',
     &       TR_FittingOptions%RadialProfiles(i)%Rwidth, ','
        write(unit,'(A,ES20.10E3,A)') '      "centPos0": ',
     &       TR_FittingOptions%RadialProfiles(i)%CentPos(0), ','
        write(unit,'(A,ES20.10E3,A)') '      "centPos1": ',
     &       TR_FittingOptions%RadialProfiles(i)%CentPos(1), ','
        write(unit,'(A,ES20.10E3,A)') '      "inclination": ',
     &       TR_FittingOptions%RadialProfiles(i)%Inclination, ','
        write(unit,'(A,ES20.10E3,A)') '      "positionAngle": ',
     &       TR_FittingOptions%RadialProfiles(i)%PositionAngle,
     &       ','
        write(unit,'(A,ES20.10E3,A)') '      "vSys": ',
     &       TR_FittingOptions%RadialProfiles(i)%VSys, ','
        write(unit,'(A,ES20.10E3,A)') '      "vRot": ',
     &       TR_FittingOptions%RadialProfiles(i)%VRot, ','
        write(unit,'(A,ES20.10E3,A)') '      "vRad": ',
     &       TR_FittingOptions%RadialProfiles(i)%VRad, ','
        write(unit,'(A,ES20.10E3,A)') '      "vDisp": ',
     &       TR_FittingOptions%RadialProfiles(i)%VDisp, ','
        write(unit,'(A,ES20.10E3,A)') '      "vvert": ',
     &       TR_FittingOptions%RadialProfiles(i)%Vvert, ','
        write(unit,'(A,ES20.10E3,A)') '      "dvdz": ',
     &       TR_FittingOptions%RadialProfiles(i)%dvdz, ','
        write(unit,'(A,ES20.10E3,A)') '      "sigUse": ',
     &       TR_FittingOptions%RadialProfiles(i)%SigUse, ','
        write(unit,'(A,ES20.10E3,A)') '      "z0": ',
     &       TR_FittingOptions%RadialProfiles(i)%z0, ','
        write(unit,'(A,ES20.10E3)') '      "zGradiantStart": ',
     &       TR_FittingOptions%RadialProfiles(i)%zGradiantStart
        if(i .lt. TR_FittingOptions%nRings-1) then
          write(unit,'(A)') '    },'
        else
          write(unit,'(A)') '    }'
        endif
      enddo
      write(unit,'(A)') '  ],'

c     ObservedDC header
      write(unit,'(A)') '  "observedDC": {'
      write(unit,'(A,I6,A)') '    "nPixelsX": ',
     &     ObservedDC%DH%nPixels(0), ','
      write(unit,'(A,I6,A)') '    "nPixelsY": ',
     &     ObservedDC%DH%nPixels(1), ','
      write(unit,'(A,I6,A)') '    "nChannels": ',
     &     ObservedDC%DH%nChannels, ','
      write(unit,'(A,ES20.10E3,A)') '    "pixelSizeX": ',
     &     ObservedDC%DH%PixelSize(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "pixelSizeY": ',
     &     ObservedDC%DH%PixelSize(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "channelSize": ',
     &     ObservedDC%DH%ChannelSize, ','
      write(unit,'(A,ES20.10E3,A)') '    "startX": ',
     &     ObservedDC%DH%Start(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "startY": ',
     &     ObservedDC%DH%Start(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "startV": ',
     &     ObservedDC%DH%Start(2), ','
      write(unit,'(A,ES20.10E3,A)') '    "refLocX": ',
     &     ObservedDC%DH%RefLocation(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "refLocY": ',
     &     ObservedDC%DH%RefLocation(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "refLocV": ',
     &     ObservedDC%DH%RefLocation(2), ','
      write(unit,'(A,ES20.10E3,A)') '    "refValX": ',
     &     ObservedDC%DH%RefVal(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "refValY": ',
     &     ObservedDC%DH%RefVal(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "refValV": ',
     &     ObservedDC%DH%RefVal(2), ','
      write(unit,'(A,ES20.10E3,A)') '    "uncertainty": ',
     &     ObservedDC%DH%Uncertainty, ','
      write(unit,'(A,I10,A)') '    "nValid": ',
     &     ObservedDC%DH%nValid, ','

c     Valid indices
      write(unit,'(A)') '    "flattendValidIndices": ['
      do l=0, ObservedDC%DH%nValid-1
        if(l .lt. ObservedDC%DH%nValid-1) then
          write(unit,'(A,I10,A)') '      ',
     &         ObservedDC%FlattendValidIndices(l), ','
        else
          write(unit,'(A,I10)') '      ',
     &         ObservedDC%FlattendValidIndices(l)
        endif
      enddo
      write(unit,'(A)') '    ],'

c     Valid flux values only
      write(unit,'(A)') '    "validFlux": ['
      do l=0, ObservedDC%DH%nValid-1
        call ThreeDIndxCalc(
     &       ObservedDC%FlattendValidIndices(l),
     &       ObservedDC%DH, ii, jj, kk)
        if(l .lt. ObservedDC%DH%nValid-1) then
          write(unit,'(A,ES20.10E3,A)') '      ',
     &         ObservedDC%Flux(ii,jj,kk), ','
        else
          write(unit,'(A,ES20.10E3)') '      ',
     &         ObservedDC%Flux(ii,jj,kk)
        endif
      enddo
      write(unit,'(A)') '    ]'
      write(unit,'(A)') '  },'

c     ObservedBeam
      write(unit,'(A)') '  "observedBeam": {'
      write(unit,'(A,ES20.10E3,A)') '    "beamSigma0": ',
     &     ObservedBeam%BeamSigmaVector(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "beamSigma1": ',
     &     ObservedBeam%BeamSigmaVector(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "beamSigma2": ',
     &     ObservedBeam%BeamSigmaVector(2), ','
      write(unit,'(A,I6,A)') '    "nRadialCells": ',
     &     ObservedBeam%nRadialCells, ','
      write(unit,'(A,ES20.10E3,A)') '    "sigmaLengths": ',
     &     ObservedBeam%SigmaLengths, ','
      write(unit,'(A,ES20.10E3,A)') '    "pixelSizeX": ',
     &     ObservedBeam%PixelSize(0), ','
      write(unit,'(A,ES20.10E3,A)') '    "pixelSizeY": ',
     &     ObservedBeam%PixelSize(1), ','
      write(unit,'(A,ES20.10E3,A)') '    "beamMajorAxis": ',
     &     ObservedBeam%BeamMajorAxis, ','
      write(unit,'(A,ES20.10E3,A)') '    "beamMinorAxis": ',
     &     ObservedBeam%BeamMinorAxis, ','

c     Kernel
      write(unit,'(A)') '    "kernel": ['
      do i=-ObservedBeam%nRadialCells,ObservedBeam%nRadialCells
        do j=-ObservedBeam%nRadialCells,ObservedBeam%nRadialCells
          if(i .eq. ObservedBeam%nRadialCells .and.
     &       j .eq. ObservedBeam%nRadialCells) then
            write(unit,'(A,ES20.10E3)') '      ',
     &           ObservedBeam%Kernel(i,j)
          else
            write(unit,'(A,ES20.10E3,A)') '      ',
     &           ObservedBeam%Kernel(i,j), ','
          endif
        enddo
      enddo
      write(unit,'(A)') '    ]'
      write(unit,'(A)') '  },'

c     Scalars
      write(unit,'(A,ES20.10E3,A)') '  "ftol": ', ftol, ','
      write(unit,'(A,I12,A)') '  "idum": ', idum, ','
      write(unit,'(A,I6,A)') '  "linearLogSDSwitch": ',
     &     PFlags%Linear_Log_SDSwitch, ','
      write(unit,'(A,I6,A)') '  "likelihoodSwitch": ',
     &     PFlags%LikelihoodSwitch, ','
c     Model cloud-generation settings (FittingOptionsInputs.f reads these
c     from the fitting-options file per run -- they are NOT fixed constants,
c     so a JS worker rebuilding the model must read them here rather than
c     hardcoding a guess).
      write(unit,'(A,I6,A)') '  "cmode": ',
     &     ModelTiltedRing%cmode, ','
      write(unit,'(A,ES20.10E3,A)') '  "cloudBaseSurfDens": ',
     &     ModelTiltedRing%CloudBaseSurfDens, ','
      write(unit,'(A,I6)') '  "nRingsPerBeam": ',
     &     TR_FittingOptions%nRingsPerBeam

      write(unit,'(A)') '}'
      close(unit)
      print*, "Fixture dump complete."

      return
      end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   DumpBestFitModelCube
c
c   Dumps the post-fit model cube flux to JSON.
c   Call after GalaxyFit() completes.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine DumpBestFitModelCube(filename)
      use PipelineGlobals
      use DataCubeMod
      implicit none

      character(*), INTENT(IN) :: filename
      integer i, j, k, l
      integer, parameter :: unit=44
      integer nCells

      print*, "Dumping best-fit model cube to: ", trim(filename)

      nCells = ModelDC%DH%nPixels(0)
     &       * ModelDC%DH%nPixels(1)
     &       * ModelDC%DH%nChannels

      open(unit, file=trim(filename),
     &     status='replace', action='write')

      write(unit,'(A)') '{'
      write(unit,'(A,I6,A)') '  "nPixelsX": ',
     &     ModelDC%DH%nPixels(0), ','
      write(unit,'(A,I6,A)') '  "nPixelsY": ',
     &     ModelDC%DH%nPixels(1), ','
      write(unit,'(A,I6,A)') '  "nChannels": ',
     &     ModelDC%DH%nChannels, ','
      write(unit,'(A,I10,A)') '  "nCells": ',
     &     nCells, ','
      write(unit,'(A)') '  "flux": ['
      l = 0
      do i=0, ModelDC%DH%nPixels(0)-1
        do j=0, ModelDC%DH%nPixels(1)-1
          do k=0, ModelDC%DH%nChannels-1
            l = l + 1
            if(l .lt. nCells) then
              write(unit,'(A,ES20.10E3,A)') '    ',
     &             ModelDC%Flux(i,j,k), ','
            else
              write(unit,'(A,ES20.10E3)') '    ',
     &             ModelDC%Flux(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(unit,'(A)') '  ]'
      write(unit,'(A)') '}'

      close(unit)
      print*, "Best-fit model cube dump complete."

      return
      end subroutine
