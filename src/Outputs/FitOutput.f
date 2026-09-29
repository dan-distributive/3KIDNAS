cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for writing out a
c       data cube to a file.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FitOutputMod
      use DataCubeOutputsMod
      use DataCubeMod
      use BeamMod
      use SoFiACatalogueMod
      use PipelineGlobals

      use ParameterVectorToTiltedRingMod
      use TiltedRingGenerationMod
      use FillDataCubeWithTiltedRingMod
      use CubeKernelConvolutionMod

      use UnitConvertMod

      use ParameterVectorToTiltedRingMod

      PROCEDURE(GeneralOutputInterface),POINTER :: OutputFit =>null()

      ABSTRACT INTERFACE
        subroutine GeneralOutputInterface(CatItem)
            import :: CatalogueItem
            Type(CatalogueItem),INTENT(IN) :: CatItem
        END subroutine GeneralOutputInterface
      END INTERFACE


      contains

ccccc
c
      subroutine OutputBestFit_Simple(CatItem)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      character(700) OutputFolder, Script


c           Make the output folder for the file
      call MakeModelOutputFolder(OutputFolder,CatItem)
c           Make a PSF convolved cube for the best fitting model
      call OutputCube(OutputFolder,trim(CatItem%ObjName))
c           Make a residual cube and output it
      call OutputResidualCube(OutputFolder,trim(CatItem%ObjName))
c           Output the tilted ring parameters
c      call OutputTiltedRingParams(OutputFolder)

c           Output the standard model output from  the proto pipeline
      call StandardModelOutput(OutputFolder,CatItem,2)

c           Write out the model moment maps
c      call WriteModelMaps(OutputFolder,CatItem)

c           Write out the model PV maps
c      call WritePVMaps(OutputFolder,CatItem)

c           Write out a file containing the set of fitting flags
      call WriteFlagFile(OutputFolder,CatItem)

c           For testing purposes, write out a file containing the
c               first fit -- for that we need to get the TR parameters
c               from the first fit parameter vector as well as the initial estimate
c      call ParamToTiltedRing(PV_FirstFit,ModelTiltedRing
c     &          ,TR_FittingOptions)
c      call StandardModelOutput(OutputFolder,CatItem,1)

c      print*, "Output initial Param",PVIni%Param(0:PVIni%nParams-1)
      call ParamToTiltedRing(PVIni,ModelTiltedRing
     &          ,TR_FittingOptions)
      call StandardModelOutput(OutputFolder,CatItem,0)
 
c           Make a diagnostic plot
c      call GenerateDiagnosticPlot()

      call CorrectRADEC()

      return
      end subroutine
ccccc


ccccc
      subroutine CorrectRADEC()
      use PipelineGlobals
      implicit none
      integer i, j,len,indx
      character(2000) PlotCmd
      character(8) ValStr
      character(500) CodePath,PythonPath
      character(9) CodeDelim

      write(ValStr,'(I4)') Version

      call GetArg(0,CodePath)

      CodeDelim="P"
      indx = SCAN(trim(CodePath),CodeDelim,.True.)
      PythonPath = CodePath(1:indx-2)
     &          //"/src/PythonScripts/GeometryFix.py "


      PlotCmd="python3 "//trim(PythonPath)//" "
     &          //trim(GalaxyDict%OutputFolder)//" "
     &          //trim(GalaxyDict%GalaxyName)//" "
     &          //trim(ValStr)//" "

      PlotCmd=trim(PlotCmd)//" "
     &          //trim(GalaxyDict%DataCubeFile)//" "


c      print*, trim(PlotCmd)
      call system(trim(PlotCmd))

      return
      end subroutine
ccccc

ccccc
c
      subroutine MakeModelOutputFolder(OutputFolder,CatItem)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(700),INTENT(INOUT):: OutputFolder
      character(700) Script
c           Make the output folder for the file
      OutputFolder=trim(MainOutputFolder)//trim(CatItem%ObjName)
      print*, "Moving outputs to ", trim(MainOutputFolder)
      print*, "Making an output folder ", trim(OutputFolder)
     &          ,trim(CatItem%ObjName)
      Script="mkdir "//trim(OutputFolder)
      call system(Script)

      GalaxyDict%OutputFolder=OutputFolder
      return
      end subroutine
ccccccc


cccccc
c
      subroutine OutputCube(OutputFolder,Name)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(*),INTENT(IN) :: Name
      real pixelarea
      character(500) CubeFile,BaseCubeFile
      character(3) VersStr
      character(10) format_string
      real BeamArea
      real SpecNoise
      character(64) EnvVal
      integer EnvLen

c           The best fitting tilted ring model should have been made in Galaxy fit, so we
c           don't need to do the conversion of a parameter vector to tilted ring parameters

c       Make the tilted ring model
      pixelarea=abs(ModelDC%DH%PixelSize(0)*ModelDC%DH%PixelSize(1))
c      BeamArea=ObservedBeam%BeamMajorAxis*ObservedBeam%BeamMinorAxis
      BeamArea=2.*Pi*abs((ObservedBeam%BeamSigmaVector(0)
     &                      *ObservedBeam%BeamSigmaVector(1)))

      SpecNoise=ObservedDC%DH%Uncertainty
     &             *abs(ObservedDC%DH%ChannelSize)
c           One-off diagnostic (Fortran-vs-JS idum-divergence proof, Dan
c               2026-08-17): if set, forces idum to a fixed, externally
c               supplied value right before this final resynthesis call,
c               so the SAME idum can be injected on both platforms and the
c               resulting model cubes compared with every other variable
c               (parameters, idum) controlled for. No effect unless
c               TRACE_OVERRIDE_IDUM is set.
      call MaybeOverrideIdum(idum)
      call MaybeTraceRingFields(ModelTiltedRing)
      block
        use BasicRanNumGen
        character(64) EnvValC
        integer EnvLenC, IDUnit, IDPI
        integer FIdum2, FIv(32), FIy
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvValC,
     &            EnvLenC)
        if (EnvLenC .gt. 0) then
          call GetRan2State(FIdum2,FIv,FIy)
          open(newunit=IDUnit, file="OutputIdumTraceF.txt",
     &        status="unknown", position="append")
          write(IDUnit,'(A,I0)') "idum ", idum
          write(IDUnit,'(A,I0)') "idum2 ", FIdum2
          write(IDUnit,'(A,I0)') "iy ", FIy
          do IDPI=1,32
            write(IDUnit,'(A,I0,A,I0)') "iv",IDPI," ",FIv(IDPI)
          enddo
          close(IDUnit)
        endif
      end block
      call BuildTiltedRingModel(ModelTiltedRing,idum,SpecNoise
     &          ,ObservedDC,ObservedBeam)
c       Create the point-source data cube
      call FillDataCubeWithTiltedRing(ModelDC,ModelTiltedRing)
c      print*, "Filled DC", sum(ModelDC%Flux),pixelarea
c           One-off diagnostic (Fortran-vs-JS pre/post-convolution
c               isolation, Dan 2026-08-18): checksum the model cube right
c               before beam convolution, so a residual divergence can be
c               attributed to particle generation/binning (present here
c               already) vs convolution (introduced after this point).
c               Gated on TRACE_OVERRIDE_IDUM.
      call MaybeTracePreConvChecksum(ModelDC)
c           One-off diagnostic (Dan probe, 2026-09-27): dump the FULL
c               pre-convolution model cube to its own FITS file (not just
c               a checksum) so it can be pixel-diffed directly against
c               JS's own pre-convolution snapshot, isolating whether the
c               particle-generation/binning stage or the beam-convolution
c               stage is where the two platforms' cubes actually start to
c               differ. Gated on TRACE_DUMP_PRECONV so it's silent by
c               default; writes alongside the real AverageModel cube.
      call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal,
     &          EnvLen)
      if (EnvLen .gt. 0) then
        call WriteDataCubeToFITS(ModelDC,ObservedBeam
     &      ,trim(OutputFolder)//"/PreConvModel.fits","Test")
      endif
c        Convolve the cube with the beam
c           Note that it is assumed that the real beam kernel has already been calculated
      call CubeBeamConvolution(ModelDC,ObservedBeam)
c           One-off diagnostic (Dan probe, 2026-09-27): dump the cube
c               immediately after convolution, BEFORE the BeamArea
c               rescale, to isolate whether CubeBeamConvolution itself
c               (the FFTW-based step) introduces the divergence, or
c               whether it's the scalar BeamArea multiply right after it.
c               Gated on TRACE_DUMP_PRECONV, same as the pre-conv dump.
      if (EnvLen .gt. 0) then
        call WriteDataCubeToFITS(ModelDC,ObservedBeam
     &      ,trim(OutputFolder)//"/PostConvPreScaleModel.fits","Test")
      endif
c           Because the cube is in units of Jy/pixel, convert back to Jy/beam
      ModelDC%Flux=ModelDC%Flux*BeamArea
c       Write the cube to a fits file
      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,format_string) Version

      BaseCubeFile=trim(GalaxyDict%GalaxyName)//"_AverageModel_v"
     &      //trim(VersStr)//".fits"
      CubeFile=trim(OutputFolder)//"/"//trim(BaseCubeFile)

      ModelDC%DH%FType='intensity'
      ModelDC%DH%PixelSize=ModelDC%DH%PixelSize/3600.

      call WriteDataCubeToFITS(ModelDC,ObservedBeam,CubeFile
     &              ,Name)

      GalaxyDict%BestFitCubeFile=trim(BaseCubeFile)



      return
      end subroutine
ccccc


cccccc
c
      subroutine OutputResidualCube(OutputFolder,Name)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(*),INTENT(IN)::Name
      character(500) CubeFile
      Type(DataCube) ResidCube
      real BeamArea

c      The model DC and the actual cube should already be in memory
      ResidCube=ModelDC
c       The Observed Cube is still in Jy/pixel so convert back to Jy/beam
c      BeamArea=abs(ObservedBeam%BeamMajorAxis
c     &              *ObservedBeam%BeamMinorAxis)
      BeamArea=2.*Pi*abs((ObservedBeam%BeamSigmaVector(0)
     &                      *ObservedBeam%BeamSigmaVector(1)))

c      print*, "Output beam Area", BeamArea,sum(ObservedDC%Flux)
c     &              ,sum(ObservedDC%Flux*BeamArea)
      ObservedDC%Flux=ObservedDC%Flux*BeamArea

c       Calculate the difference cube
      ResidCube%Flux=ObservedDC%Flux-ModelDC%Flux
      print*, "Flux Checks", sum(ObservedDC%Flux)
     &              ,sum(ModelDC%Flux),sum(ResidCube%Flux)

c       Write the cube to a fits file
      CubeFile=trim(OutputFolder)//"/"//"DifferenceCube.fits"
      call WriteDataCubeToFITS(ResidCube,ObservedBeam,CubeFile,Name)

      return
      end subroutine
ccccc

ccccc
c
      subroutine OutputTiltedRingParams(OutputFolder)
      implicit none
      character(500),INTENT(IN):: OutputFolder
      character(500) ParamsFile

      integer i
      real SDConv, BeamArea,BeamPixels
      real CentPix(0:1),PAOut
      real SOut,SDTemp

c       The surface density units of the TR fitter are natively Jy/pixel.
c           The standard TR output should be Jy km/s arcsecond^-2.
c           Set the SDConversion factor to
      SDConv=abs(ObservedDC%DH%ChannelSize)
     &          /(abs(ObservedDC%DH%PixelSize(0))
     &          *abs(ObservedDC%DH%PixelSize(1)))



      ParamsFile=trim(OutputFolder)//"/"//"TiltedRingParams.txt"
      open(10,file=trim(ParamsFile),status='replace')
      write(10,*) "#    Best fitting Tilted Ring Params"
      write(10,*) "#    Generated by the WRKP"
      write(10,*) "#    Fitting options"
      write(10,*) "#    Core Code algorithm"
      if(PFlags%CoreCodeSwitch .eq. 1) then
        write(10,*) "WRKP internal algorithm"
      elseif(PFlags%CoreCodeSwitch .eq. 2) then
        write(10,*) "3DBarolo"
      endif
      write(10,*) "#    Fitting options"
      write(10,*) "cmode= ", ModelTiltedRing%cmode
      write(10,*) "cdens= ", ModelTiltedRing%CloudBaseSurfDens
      write(10,*) "rings/beam= ", TR_FittingOptions%nRingsPerBeam

      write(10,*) "#    Best Fit"
      write(10,*) "Likelihood=", PVModel%BestLike

      write(10,*) "#    Rmid (arcsec)    Rwidth(arcsec)"
     &          //"    Xcent (pixels)"
     &          //"    Ycent (pixels)"
     &          //"    Inc (degrees)"
     &          //"    PA (degrees)   VSys (km/s)"
     &          //"    VRot (km/s)    VRad (km/s)"
     &          //"    Vvert (km/s)    VDisp (km/s)"
     &          //"    dvdz (km/s)    Sigma (Jy km/s arcsec^-2"
     &          //"    z0  ('')   zGradStart ('')   "


      do i=0, ModelTiltedRing%nRings-1
        PAOut=(ModelTiltedRing%R(i)%PositionAngle*180./Pi)
        PAOut=PAOut-90.
100     continue
        if(PAOut .lt. 0.) then
            PAOut=PAOut+360.
            goto 100
        elseif(PAOut .gt. 360.) then
            PAOut=PAOut-360.
            goto 100
        endif
        if(PFlags%Linear_Log_SDSwitch .eq. 0) then
            SDTemp=ModelTiltedRing%R(i)%SigUse
        elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
            SDTemp=10.**(ModelTiltedRing%R(i)%SigUse)
        endif


c        print*, "Final SD and VRot"
c     &          ,ModelTiltedRing%R(i)%Rmid
c     &          ,ModelTiltedRing%R(i)%VRot
c     &          ,ModelTiltedRing%R(i)%Sigma,SDConv

        write(10,*) ModelTiltedRing%R(i)%Rmid
     &                  *abs(ObservedDC%DH%PixelSize(0))
     &          , ModelTiltedRing%R(i)%Rwidth
     &                  *abs(ObservedDC%DH%PixelSize(0))
     &          , ModelTiltedRing%R(i)%CentPos
     &          , ModelTiltedRing%R(i)%Inclination*180./Pi
     &          , PAOut
     &          , ModelTiltedRing%R(i)%VSys, ModelTiltedRing%R(i)%VRot
     &          , ModelTiltedRing%R(i)%VRad, ModelTiltedRing%R(i)%Vvert
     &          , ModelTiltedRing%R(i)%VDisp, ModelTiltedRing%R(i)%dvdz
     &          , SDTemp*SDConv
     &          , ModelTiltedRing%R(i)%z0
     &          , ModelTiltedRing%R(i)%zGradiantStart


      enddo
      close(10)

      return
      end subroutine
ccccc


ccccccc
c
      subroutine StandardModelOutput(OutputFolder,CatItem
     &                  ,FitNum)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500),INTENT(IN):: OutputFolder
      integer,INTENT(IN) :: FitNum
      character(500) ParamsFile
      character(10) VersStr, format_string
      character(100) DateStr

      character(300) OutStr

c     BUG FIX (2026, flagged by Dan): these were character(8), matching the
c     original F8.2/F5.2/F8.3 formats below -- since GalaxyDict['BestFitModel']
c     is the ONLY handoff into the bootstrap-resample geometry calculation
c     (MakeBootstrapSample.WriteBootstrapFile / RunBootstrapsDCP.ComputeBsCent,
c     both read Model['XCENTER']/['VSYS']/etc back from THIS text file),
c     truncating to 2-3 decimal digits here was silently discarding real4
c     precision (~7 significant digits) before every bootstrap realization's
c     resample center was computed from it. Confirmed by direct comparison:
c     a UseDCP=1 run whose own initial fit skips this text file (JS/DCP,
c     js/bootstrap-realization-launcher.js's runInitialFit)
c     seeded its bootstrap loop from full precision while a UseDCP=0 run
c     seeded from this truncated text -- individual bootstrap realizations'
c     PA_model diverged by up to 6.6 degrees between the two, an effect size
c     four to five orders of magnitude larger than the ~1e-7-relative
c     perturbation already shown elsewhere in this codebase to be enough to
c     flip a near-tied Nelder-Mead simplex to a different local minimum.
c     Widened to character(16) here, paired with F16.6 (a format that fills
c     the ENTIRE declared width, not just widened in place) on every write
c     into ValStr/ErrStr in this subroutine -- ValStr/ErrStr are reused
c     across many fields per call, so a format narrower than the declared
c     width would leave undefined trailing bytes from a previous write,
c     corrupting output. ReadWRKPFit.py's parser already splits on
c     whitespace (GeoLineAssign/NoiseLineAssign), not fixed columns, so
c     widening this needs no change on the Python side.
c     BUG FIX (2026-09-28, Dan): F16.6 (6 decimal places) is NOT sufficient
c     to round-trip a real(4) value losslessly through Python's DOUBLE-
c     PRECISION consumer (RunBootstrapsDCP.ComputeBsCent's
c     VCenter=DeltaV/dV+RefChan arithmetic operates directly on the parsed
c     text, never re-quantizing to float32) -- confirmed directly: VSys_kin
c     wrote "5745.641602" for the true value 5745.6416015625, a ~4.4e-7
c     truncation that is small enough to still round-trip to the SAME
c     nearest float32 bit pattern (explaining why this looked "good enough"
c     when the F8.2->F16.6 fix was first made), but is NOT reabsorbed by
c     ComputeBsCent's float64 division, producing a genuine, confirmed
c     ~1 ULP-scale CentV divergence between the Fortran-local and JS-local
c     pipelines' bootstrap-resample geometry -- traced via a direct
c     BS_Cent%CentV comparison after independently proving both platforms'
c     amoeba trajectory AND raw best-fit parameter vector bit-identical, so
c     this text round-trip was the only remaining place precision could be
c     lost. X_kin/Y_kin happened to round-trip losslessly through the same
c     F16.6 format by luck of their specific decimal digits (rounding down
c     at the 7th place instead of up), not because F16.6 is actually safe.
c     Switched to list-directed (format-free) internal writes, matching the
c     already-correct, already-verified-exact approach RawGeom_v1.txt uses
c     for PA/Inc. Widened the buffer accordingly (list-directed real4
c     output needs more than 16 characters for some values/signs).
      character(30) ValStr,ErrStr
      character(18) PreambleStr
c     BUG FIX (2026-09-29, Dan): same class of bug as the geometry fields
c     above (F16.6 -> list-directed) but for the radial-profile table --
c     Rad/VRot were F8.2 (2 decimal places) and SD_kin was G9.2 (TWO
c     SIGNIFICANT FIGURES, e.g. "2.9" for a true value of 2.9300459...).
c     ExtractScalingParams.py's RHI extraction interpolates against this
c     exact SD_model/Rad profile (read back via ReadWRKPFit.ProfileLineAssign,
c     shared code path for both the Fortran-local and JS-local legs), so
c     quantizing the profile to ~2 significant digits before RHI's
c     threshold-crossing interpolation was amplifying into RHI_AS diffs of
c     up to 4.45% between the two legs -- confirmed directly: Fortran's own
c     BootstrapFits.csv showed "SD_model" values like "2.9, 2.0, 0.31, 1.2"
c     while JS's own (never quantized) profile showed "2.9300459036646607,
c     2.0146256636770046, 0.31207922029800184, 1.1935189127727546" for the
c     SAME realization. Widened to character(30) (list-directed real4
c     output needs more than 20 characters for some values/signs, same
c     reasoning as ValStr/ErrStr above).
      character(30) RadialProfStr(6)

      character(8) date
      character(10) time
      character(5) zone
      
      integer timeArray(8)
      integer i
      real PAOut
      real BeamPixels,SDConv1
      real SDTemp
      real CentValAS(2),RA,DEC

      call date_and_time(date,time,zone,timeArray)




c       The machine SD units are Jy/pixel and we want them in M_sol/pc^2
c           First get the conversion to Jy/arcsec^2
c           This calculation requires that the pixelsize arrays be in arcseconds
      SDConv1=abs(ObservedDC%DH%ChannelSize)
     &          /(abs(ObservedDC%DH%PixelSize(0))
     &          *abs(ObservedDC%DH%PixelSize(1)))

      print*, "SD Convert factor", SDConv1


      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,format_string) Version

c      print*, "Obj Name", trim(CatItem%ObjName)

      call NameParamFile(OutputFolder,CatItem,FitNum,ParamsFile)
c      print*,trim(ParamsFile)

      open(10,file=trim(ParamsFile),status='replace')

c       Write the Preamble
      write(10,'(a,a)') "Object:  ",trim(CatItem%ObjName)
      write(10,'(a,a)') "Source:  ",trim(SCatLocal%SourceName)
      DateStr="Date:    "//date(7:8)//"-"//date(5:6)//"-"//date(1:4)
      write(10,'(a)') trim(DateStr)
      write(10,'(a,a)') "Version: ",trim(VersStr)

c       Write out the noise measurements
c           First the RMS--which needs to be converted
      write(10,*) " "
      PreambleStr="RMS (mJy/beam)"
      write(ValStr,*)DBLE(ObservedDC%DH%Uncertainty
     &              *ObservedBeam%BeamAreaPixels*1000.)
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       Next the integrated S/N
      PreambleStr="SN_Integrated "
      write(ValStr,*)DBLE(ObservedDC%DH%SN_Int)
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the peak S/N
      PreambleStr="SN_Peak "
      write(ValStr,*)DBLE(ObservedDC%DH%SN_Peak)
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the average S/N
      PreambleStr="SN_Avg "
      write(ValStr,*)DBLE(ObservedDC%DH%SN_Avg)
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)
c       And the median S/N
      PreambleStr="SN_Median "
      write(ValStr,*)DBLE(ObservedDC%DH%SN_Median)
      OutStr=PreambleStr//"    "//ValStr
      write(10,'(a)') trim(OutStr)


c       Write out the geometric model parameters

      write(10,*) " "
      write(10,'(a)') "Geometry Paramters"
      write(10,'(a)') "Param Name            Value    Error"

c      print*, "Reference values",ObservedDC%DH%RefVal
c     &                  ,ObservedDC%DH%RefLocation
c       Get the RA and DEC value for the central point
      do i=0,1
        CentValAS(i+1)=(ModelTiltedRing%R(0)%CentPos(i)
     &              -ObservedDC%DH%RefLocation(i))
     &                  *ObservedDC%DH%PixelSize(i)
     &              +ObservedDC%DH%RefVal(i)
      enddo
c      print*, "Central Position AS", CentValAS
      call ArcSecToDegrees(CentValAS(1),RA)
      call ArcSecToDegrees(CentValAS(2),DEC)
c      print*, "Central Position Deg", RA,DEC

      do i=0, 7
        if( i.eq. 0) then
            PreambleStr="X_kin (pixels)"
            write(ValStr,*)DBLE(ModelTiltedRing%R(0)%CentPos(0))
            write(ErrStr,*)0.0
        elseif(i .eq. 1) then
            PreambleStr="Y_kin (pixels)"
            write(ValStr,*)DBLE(ModelTiltedRing%R(0)%CentPos(1))
            write(ErrStr,*)0.0
        elseif(i .eq. 2) then
            PreambleStr="RA_kin (degrees)"
            write(ValStr,*)DBLE(RA)
            write(ErrStr,*)0.0
        elseif(i .eq. 3) then
            PreambleStr="DEC_kin (degrees)"
            write(ValStr,*)DBLE(DEC)
            write(ErrStr,*)0.0
        elseif(i .eq. 4) then
            PreambleStr="Inc_kin (degrees)"
            write(ValStr,*)DBLE(ModelTiltedRing%R(0)%Inclination
     &                  *180./Pi)
            write(ErrStr,*)0.0
        elseif(i .eq. 5) then
            PreambleStr="PA_kin (degrees)"
c           MIRRORED IN JS: this -90 deg convention offset + [0,360) wrap is
c           duplicated in js/bootstrap-realization-launcher.js's
c           toKinematicPA() (the JS worker has no Fortran binary available at
c           runtime, so it can't call back into this routine). If you change
c           this transform, update that function too.
            PAOut=(ModelTiltedRing%R(0)%PositionAngle*180./Pi)
            PAOut=PAOut-90.
 100        continue
            if(PAOut .lt. 0.) then
                PAOut=PAOut+360.
                goto 100
            elseif(PAOut .gt. 360.) then
                PAOut=PAOut-360.
                goto 100
            endif
            write(ValStr,*)DBLE(PAOut)
            write(ErrStr,*)0.0
        elseif(i .eq. 6) then
            PreambleStr="VSys_kin (km/s)"
            write(ValStr,*)DBLE(ModelTiltedRing%R(0)%VSys)
            write(ErrStr,*)0.0
        elseif(i .eq. 7) then
            PreambleStr="VDisp_kin (km/s)"
            write(ValStr,*)DBLE(ModelTiltedRing%R(0)%VDisp)
            write(ErrStr,*)0.0
        endif

        OutStr=PreambleStr//"    "//ValStr//" "//ErrStr
        write(10,'(a)') trim(OutStr)
      enddo

c       One-off diagnostic turned permanent fix (Dan, 2026): bootstrap
c       resampling's geometry (WriteBootstrapFile/computeBsCent) used to
c       re-derive PA by reading PA_kin back OUT of this file (the F16.6
c       degrees value just written above, itself already round-tripped
c       through a -90/wrap/+90 "kinematic PA" display convention) and
c       reversing that convention in Python -- two independent, lossy
c       re-derivations (this one through F16.6 text, JS's through its own
c       in-memory but still degrees-and-back recomputation) that have no
c       reason to agree bit-for-bit. Traced directly via a COORDTRACE/
c       ROTTRACE hex bisection to a 1-ULP PA difference propagating into
c       every cell's coordinate transform via fd_cos/fd_sin(-PA). Fix:
c       write the RAW, pre-"kinematic" angle (no -90 offset, no [0,360)
c       wrap) to its own small companion file, full list-directed
c       precision, so WriteBootstrapFile/computeBsCent can consume it
c       directly instead of each re-deriving their own approximation.
c       GATED ON FitNum.eq.2 (the converged/"AvgModel" call) ONLY --
c       StandardModelOutput is ALSO called with FitNum=0 further down in
c       OutputBestFit_Simple, AFTER ModelTiltedRing has been overwritten
c       with PVIni (the INITIAL GUESS, not the fit) for the "IniEstimate"
c       output. An earlier, unguarded version of this write used a fixed
c       filename with no FitNum check, so that SECOND call silently
c       clobbered this file with the initial guess's PA instead of the
c       converged fit's -- confirmed directly by tagging a diagnostic
c       print with FitNum+ObjName: FitNum=2 printed 2.8587160110 (matches
c       GalaxyFit.f's own FINALVEC for this run); FitNum=0 printed
c       2.9892427921 (the initial guess) and was the one left on disk.
      if (FitNum .eq. 2) then
        open(11,file=trim(OutputFolder)//"/"//trim(CatItem%ObjName)
     &            //"_RawGeom_v1.txt",status='replace')
        write(11,*) ModelTiltedRing%R(0)%PositionAngle
        write(11,*) ModelTiltedRing%R(0)%Inclination
        close(11)
      endif

c       Write out the radial profiles
      write(10,*)" "

      write(10,'(a)') " Radial Profiles"
      write(10,'(a)')" Rad        VROT_kin    e_VRot_kin"
     &           //"       e_VROT_kin,inc      SD_kin      e_SD_kin"
      write(10,'(a)') "('')       (km/s)         (km/s)     "
     &              //"(km/s)             (Msol/pc^2)      (Msol/pc^2)"

      do i=0,ModelTiltedRing%nRings-1
c       Get the Surface density in units of M_sol/pc^2
        if(PFlags%Linear_Log_SDSwitch .eq. 0) then
            SDTemp=ModelTiltedRing%R(i)%SigUse
        elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
            SDTemp=10.**(ModelTiltedRing%R(i)%SigUse)
        endif

        SDTemp=SDTemp*SDConv1
        call JyAS2_To_MSolPc2(SDTemp,SDTemp)

        print*, "Rmid",ModelTiltedRing%R(i)%Rmid
     &                  ,abs(ObservedDC%DH%PixelSize(0))
        write(RadialProfStr(1),*)DBLE(ModelTiltedRing%R(i)%Rmid
     &                  *abs(ObservedDC%DH%PixelSize(0)))
        write(RadialProfStr(2),*)DBLE(ModelTiltedRing%R(i)%VRot)
        write(RadialProfStr(3),*)0.0D0
        write(RadialProfStr(4),*)0.0D0
        write(RadialProfStr(5),*)DBLE(SDTemp)
        write(RadialProfStr(6),*)0.0D0
       OutStr=RadialProfStr(1)//RadialProfStr(2)//RadialProfStr(3)
     &          //RadialProfStr(4)//RadialProfStr(5)//RadialProfStr(6)
        write(10,'(a)') trim(OutStr)

      enddo

      close(10)

      return
      end subroutine
ccccccc

ccccc
      subroutine NameParamFile(OutputFolder,CatItem
     &              ,FitNum,ParamsFile)
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500),INTENT(IN):: OutputFolder
      integer,INTENT(IN) :: FitNum
      character(500),INTENT(INOUT):: ParamsFile
      
      character(8) VersStr,format_string

      if(Version .le. 10) then
        format_string = "(I1)"
      endif
      write(VersStr,trim(format_string)) Version
    

      if(FitNum .eq. 2) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_AvgModel_v"
     &          //trim(VersStr)//".txt"
      elseif(FitNum .eq. 1) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_FirstModel_v"
     &          //trim(VersStr)//".txt"
      elseif(FitNum .eq. 0) then
        GalaxyDict%BestFitModelFile=trim(CatItem%ObjName)
     &          //"_IniEstimate_v"
     &          //trim(VersStr)//".txt"
      endif
      ParamsFile=trim(OutputFolder)//"/"
     &          //trim(GalaxyDict%BestFitModelFile)
      return
      end subroutine
ccccccc


cccccccc
      subroutine GenerateDiagnosticPlot()
      use PipelineGlobals
      implicit none

      integer i, j,len
      character(2000) PlotCmd

      character(8) ValStr
      character(500) MaskTemp,BaseMaskFile
      character(500) CodePath,PythonPath
      character(9) CodeDelim
      integer indx

      write(ValStr,'(I4)') Version

c           The MaskFile string needs to be adjusted in case there are any spaces
c           First get the length of the trimmed file name
c       Initialize the new string
      if(MomentMapSwitch .eq. 0) then
        BaseMaskFile=trim(GalaxyDict%MaskFile)
        MaskTemp=GalaxyDict%MaskFile(1:1)
      elseif(MomentMapSwitch .eq. 3) then
        BaseMaskFile=trim(SecondaryMaskFile)
        MaskTemp=SecondaryMaskFile(1:1)
      endif
      len=len_trim(BaseMaskFile)
      print*, "hmmm", len
      print*, trim(BaseMaskFile)
c           Loop through the trimmed file name
      do i=2,len
c           Replace spaces with '\ ' for terminal commands
        if(BaseMaskFile(i:i) .eq. " ") then
            MaskTemp=trim(MaskTemp)//"\"    !   Only add a '\' here as the next trim would remove the space
        else
            j=len_trim(MaskTemp)
c           Now double check if the last character is a '\', which needs a space to trail it
            if( MaskTemp(j:j) .eq. "\") then
                MaskTemp=trim(MaskTemp)//" "//BaseMaskFile(i:i)
            else
                MaskTemp=trim(MaskTemp)//BaseMaskFile(i:i)
            endif
        endif
      enddo


c      call GetCWD(CodePath)
      call GetArg(0,CodePath)
      print*, "Current Code Path ", trim(CodePath)

      CodeDelim="P"
      indx = SCAN(trim(CodePath),CodeDelim,.True.)
      PythonPath = CodePath(1:indx-2)
     &          //"/src/PythonScripts/BestFitDiagnosticPlot.py "

      print*, "Trial Python Path ", indx,trim(PythonPath)

      PlotCmd="python3.9 "//trim(PythonPath)//" "
     &          //trim(GalaxyDict%OutputFolder)//" "
     &          //trim(GalaxyDict%GalaxyName)//" "
     &          //trim(ValStr)//" "

      if(MomentMapSwitch .eq. 0) then
        PlotCmd=trim(PlotCmd)//" "
     &          //trim(GalaxyDict%DataCubeFile)//" "
      elseif(MomentMapSwitch .eq. 3) then
        PlotCmd=trim(PlotCmd)//" "
     &          //trim(SecondaryCubeFile)//" "
      endif

      PlotCmd=trim(PlotCmd) //" "//trim(MaskTemp)


      write(ValStr,'(F8.2)') GalaxyDict%Distance
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)

      write(ValStr,'(F8.3)') GalaxyDict%Size
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)

      write(ValStr,'(F5.2)') GalaxyDict%logMass
      PlotCmd=trim(PlotCmd)//" "//trim(ValStr)
      print*, trim(PlotCmd)
      call system(trim(PlotCmd))

    

      return
      end subroutine
cccccccc


cccccc
c
      subroutine WriteModelMaps(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(1) MomNum
      integer i
      character(500) MomMapName
      character(8) ValStr

      write(ValStr,'(I1)') Version

      print*, "Writing Model maps"

      do i=0,2
        write(MomNum,'(I1)') i
        print*, MomNum

        MomMapName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_Mom"//trim(MomNum)
     &              //"_v"//trim(ValStr)
     &              //".fits"
        print*, trim(MomMapName)
        call WriteDCSliceToFITS(ModelMaps,ObservedBeam
     &              ,MomMapName,i,trim(CatItem%ObjName))
      enddo

      return
      end subroutine
ccccccc




cccccc
c
      subroutine WritePVMaps(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500) PVMapName
      character(20) Suffix
      integer i


      do i=0,1
        if(i .eq. 0) then
            Suffix="MajorAxisPV"
        elseif(i .eq. 1) then
            Suffix="MinorAxisPV"
        endif
        PVMapName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_"//trim(Suffix)//".fits"
      call WriteDCSliceToFITS(ModelPVMaps,ObservedBeam
     &              ,PVMapName,i,trim(CatItem%ObjName))
      enddo

      return
      end subroutine
ccccccc





cccccc
c
      subroutine WriteFlagFile(OutputFolder,CatItem)
      implicit none

      character(*), INTENT(IN) :: OutputFolder
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(500) FlagFileName
      character(8) ValStr

      write(ValStr,'(I1)') Version

      FlagFileName=trim(OutputFolder)//"/"
     &              //trim(CatItem%ObjName)
     &              //"_Flags_v"//trim(ValStr)//".txt"

      open(10,file=FlagFileName,status='replace')
c      write(10,'(a)') "# Goodness Of Fit  type (1==chi, 3 == log10(chi)"
c      write(10,*) PFlags%LikelihoodSwitch


      write(10,'(a)') "# Goodness Of Fit "
      write(10,*) GalaxyDict%GoodnessOfFit

      write(10,'(a)') "# RMS noise "
      write(10,*) ObservedDC%DH%Uncertainty

      write(10,'(a)') "# Total cube Cells "
      write(10,*) GalaxyDict%nCells

      write(10,'(a)') "# Normalized Goodness Of Fit "
      write(10,*) GalaxyDict%Norm_GoodnessOfFit

      write(10,'(a)') "# Size Flag (1-> <2beams, 2-> >10 beams)"
      write(10,*) GalaxyDict%Flags%SizeFlag

      write(10,'(a)') "# Pre-analysis center Flag (1-> no convergence)"
      write(10,*) GalaxyDict%Flags%CenterFlag

      close(10)

      return
      end subroutine
ccccccc

ccccccc
c       MaybeTraceRingFields: one-off diagnostic (Fortran-vs-JS ring-
c           geometry divergence proof, Dan 2026-08-17) -- dumps ring 0's
c           full field set (not just the 13 fields TRACE_OVERRIDE_PARAMS
c           actually controls -- constant/fixed fields like VDisp, VRad,
c           Vvert, dvdz, z0, zGradiantStart are pulled from this fit's own
c           TiltedRingFittingOptions defaults, not from the param vector,
c           and were never forced to match between platforms). Gated on
c           TRACE_OVERRIDE_IDUM so it only fires during a deliberate
c           controlled-comparison run.
      subroutine MaybeTraceRingFields(TR)
      use TiltedRingGenerationMod
      implicit none
      Type(TiltedRingModel), INTENT(IN) :: TR
      character(64) EnvVal
      integer EnvLen

      call get_environment_variable("TRACE_OVERRIDE_IDUM",EnvVal,
     &          EnvLen)
      if (EnvLen .gt. 0) then
        print '(A,11F14.6)', 'RINGTRACE',
     &      TR%R(0)%Rmid,TR%R(0)%Rwidth,
     &      TR%R(0)%Inclination,TR%R(0)%PositionAngle,
     &      TR%R(0)%VSys,TR%R(0)%VRot,TR%R(0)%VRad,
     &      TR%R(0)%VDisp,TR%R(0)%Vvert,TR%R(0)%dvdz,
     &      TR%R(0)%SigUse
        print '(A,4F14.6)', 'RINGTRACE2',
     &      TR%R(0)%z0,TR%R(0)%zGradiantStart,
     &      TR%R(0)%CentPos(0),TR%R(0)%CentPos(1)
      endif

      return
      end subroutine
ccccccc

ccccccc
c       MaybeTracePreConvChecksum: one-off diagnostic (Fortran-vs-JS
c           pre/post-convolution isolation, Dan 2026-08-18) -- prints a
c           sum/min/max checksum of the model cube's flux right before
c           beam convolution runs. Gated on TRACE_OVERRIDE_IDUM.
      subroutine MaybeTracePreConvChecksum(DC)
      use DataCubeMod
      implicit none
      Type(DataCube), INTENT(IN) :: DC
      character(64) EnvVal
      integer EnvLen

      call get_environment_variable("TRACE_OVERRIDE_IDUM",EnvVal,
     &          EnvLen)
      if (EnvLen .gt. 0) then
        print '(A,ES25.17,1X,ES25.17,1X,ES25.17)', 'PRECONVTRACE',
     &      sum(DC%Flux), minval(DC%Flux), maxval(DC%Flux)
        print '(A,8F16.8)', 'PRECONVPIX',
     &      DC%Flux(14,23,89),DC%Flux(14,23,90),
     &      DC%Flux(14,24,89),DC%Flux(14,24,90),
     &      DC%Flux(15,23,89),DC%Flux(15,23,90),
     &      DC%Flux(15,24,89),DC%Flux(15,24,90)
      endif

      return
      end subroutine
ccccccc

ccccccc
c       MaybeOverrideIdum: one-off diagnostic (Fortran-vs-JS idum-
c           divergence proof, Dan 2026-08-17) -- reads TRACE_OVERRIDE_IDUM
c           from the environment and, if set to a valid integer,
c           overwrites idum with it. No effect otherwise (leaves idum
c           untouched, including on a blank/unset/unparseable value).
      subroutine MaybeOverrideIdum(idum)
      use BasicRanNumGen
      implicit none
      integer, INTENT(INOUT) :: idum
      character(64) EnvVal
      integer EnvLen, IOStat, NewIdum
      integer FullIdum2, FullIv(32), FullIy

      call get_environment_variable("TRACE_OVERRIDE_IDUM",EnvVal,
     &          EnvLen)
      if (EnvLen .gt. 0) then
        read(EnvVal(1:EnvLen),*,IOSTAT=IOStat) NewIdum
        if (IOStat .eq. 0) then
          print*, "TRACE using overridden idum for resynthesis:",
     &          NewIdum
          idum=NewIdum
c               See random.f's module header for why a mid-run idum
c                   reset alone doesn't reset gasdev's separate cached-
c                   spare-value state -- must clear it explicitly too,
c                   for a true apples-to-apples reset against JS's
c                   always-fresh makeRng().
          call ResetGasdevCache()
        endif
      else
c           One-off diagnostic (Dan probe, 2026-09-27): print the NATURAL
c               (un-overridden) idum right before it's used for output
c               resynthesis, so it can be captured and fed into JS's own
c               TRACE_OVERRIDE_IDUM for a controlled cross-platform test
c               -- without disturbing Fortran's own run at all. Gated on
c               WRKP_TRACE_DEBUG so it's silent by default.
        call get_environment_variable("WRKP_TRACE_DEBUG",EnvVal,
     &            EnvLen)
        if (EnvLen .gt. 0) then
          print*, "TRACE natural idum before resynthesis:", idum
          call GetRan2State(FullIdum2,FullIv,FullIy)
          print*, "TRACE natural idum2 before resynthesis:",FullIdum2
          print*, "TRACE natural iy before resynthesis:",FullIy
          print*, "TRACE natural iv before resynthesis:",FullIv
        endif
      endif

      return
      end subroutine
ccccccc

      end module
