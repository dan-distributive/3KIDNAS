cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the general procedures for modeling
c       a galaxy once the pre-analysis is completed.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GalaxyFitMod
      use ParameterVectorMod
      use DataCubeMod
      use SoFiACatalogueMod
      use BeamMod
      use CalcBeamKernelMod
      use DownhillSimplexMod
      use FullModelComparisonMod
      implicit none

      PROCEDURE(GeneralFitInterface),POINTER :: GalaxyFit =>null()

      ABSTRACT INTERFACE
        subroutine GeneralFitInterface(CatItem)
            import :: CatalogueItem
            Type(CatalogueItem),INTENT(IN) :: CatItem
        END subroutine GeneralFitInterface
      END INTERFACE

      contains
cccccccc
c
      subroutine GalaxyFit_Simple(CatItem)
      use PipelineGlobals
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      real,ALLOCATABLE :: paramGuesses(:,:),chiArray(:)
      real chi2
      integer i,j,iter,StrictEstimate
      character(500) SimplexDumpPath
      integer SimplexDumpLen

      logical FitFlag

      character(100) FullEnsambleOutput

      print*, "Fitting Galaxy",PID

      print*, "Initial PV"
      do i=0, PVIni%nParams-1
        print*, i,PVini%Param(i)
      enddo
c
c       Set up the beam
c
c      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
      call Calculate2DBeamKernel(ObservedBeam,ObservedDC%DH%PixelSize)

c       Allocate a model datacube with the same dimensions as the observed cube
      ModelDC%DH=ObservedDC%DH

      call AllocateDataCube(ModelDC)
c       Copy the initial parameter guess into the model parameter vector
      PVModel%nParams=PVIni%nParams
      call AllocateParamVector(PVModel)
      PVModel%Param=PVIni%Param
      PVModel%ParamLowerLims=PVIni%ParamLowerLims
      PVModel%ParamUpperLims=PVIni%ParamUpperLims
      PVModel%CyclicSwitch=PVIni%CyclicSwitch
      PVModel%ParamRange=PVIni%ParamRange

c           Set up the array of parameter guesses and chi^2 values needed
      ALLOCATE(paramGuesses(PVModel%nParams+1,
     &          PVModel%nParams))
      ALLOCATE(chiArray(PVModel%nParams+1))
c           Create an array of parameter guesses


      if (TraceSwitch .eq. 1) then
        print*, "IDUM BEFORE FIRST CALL",idum
      endif
      call TiltedRingModelComparison(PVModel%Param,chi2)
      print*, "Initial model fit", chi2
      block
        character(64) EnvVal5
        integer EnvLen5, CIUnit
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal5,
     &            EnvLen5)
        if (EnvLen5 .gt. 0) then
          open(newunit=CIUnit, file="Chi2IniTraceF.txt",
     &        status="unknown", position="append")
          write(CIUnit,'(A,ES27.19)') "chi2Ini ", DBLE(chi2)
          close(CIUnit)
        endif
      end block

c       ProbeSwitch: evaluate the objective function at specific parameter
c           vectors (read from probe_params.txt) with a fixed, reset RNG
c           seed each time, instead of running the optimizer. Used to
c           directly compare Fortran's chi2 against the JS port's at
c           specific points -- see PipelineGlobals.f.
      if (ProbeSwitch .eq. 1) then
        call RunObjectiveProbe(PVModel%nParams)
        print*, "ProbeSwitch set -- stopping after probe evaluation"
        stop
      endif

      IniGuessWidth=1.
      ftol=0.005

c           Make the set of initial guesses
      StrictEstimate=1
      call MakeParamGuessArray(PVModel,ParamGuesses
     &          ,PVIni%nParams,idum,IniGuessWidth,StrictEstimate)

      call DownhillSimplexRun(PVModel%nParams
     &                  ,paramGuesses,chiArray)


      PV_FirstFit%nParams=PVIni%nParams
      call AllocateParamVector(PV_FirstFit)
      PV_FirstFit%Param(0:PVModel%nParams-1)=
     &              ParamGuesses(1,1:PVModel%nParams)
      PV_FirstFit%BestLike=chiArray(1)

      IniGuessWidth=0.25
      ftol=ftol/5.

c           Make the set of initial guesses
      StrictEstimate=0
      call MakeParamGuessArray(PVModel,ParamGuesses
     &          ,PVIni%nParams,idum,IniGuessWidth,StrictEstimate)

c           One-off diagnostic (Dan, 2026): dump the FULL pass-2 starting
c           simplex (all nParams+1 vertices, full float32 precision) to a
c           file, so the exact same starting simplex can be force-fed into
c           the JS port's amoeba directly -- a controlled experiment to
c           tell apart "amoeba itself has a bug" from "the two platforms
c           start from very slightly different points and Nelder-Mead's
c           well-known sensitivity to that does the rest" for a
c           realization where JS burns way more iterations than Fortran
c           (or fails to converge) on what should be an equivalent fit.
c           No effect unless FORTRAN_SIMPLEX_DUMP_PATH is set. Suffixed
c           with CatItem%ObjName (e.g. "..._Bootstrap_4") so a whole
c           parallel batch of realizations each gets its own dump file
c           instead of racing to overwrite one shared path.
      call get_environment_variable("FORTRAN_SIMPLEX_DUMP_PATH",
     &          SimplexDumpPath,SimplexDumpLen)
      if (SimplexDumpLen .gt. 0) then
        open(77,file=trim(SimplexDumpPath)//"."//trim(CatItem%ObjName),
     &          status='replace')
        do i=1,PVModel%nParams+1
          do j=1,PVModel%nParams
            write(77,'(Z8.8)') transfer(ParamGuesses(i,j),0)
          enddo
        enddo
        close(77)
      endif

c      print*, "PV Model",PVModel%Param

      call DownhillSimplexRun(PVModel%nParams
     &                  ,paramGuesses,chiArray)

c           One-off diagnostic (Dan, 2026): dump the ACTUAL final returned
c           vector (PVModel%Param right after the SECOND/refined-pass
c           DownhillSimplexRun) -- not inferred from the per-evaluation
c           TRACE/FULLVEC history, which doesn't distinguish "the vertex
c           the simplex settled on" from "whatever the last/lowest-chi2
c           individual evaluation happened to probe" (Nelder-Mead's last
c           call is often a rejected probe, and the historical chi2-
c           minimum trial isn't always the final kept vertex either -- see
c           the very next comment block's own "bisection paradox" note for
c           a real, related, already-fixed bug of this exact class). This
c           is the single source of truth for "what did the optimizer
c           actually decide."
      if (TraceSwitch .eq. 1) then
        print '(A)', "FINALVEC"
        do i=0,PVModel%nParams-1
          print '(A,I3,1X,Z8.8,1X,ES17.9)', "FINALVECPARAM",i+1
     &          ,transfer(PVModel%Param(i),0),PVModel%Param(i)
        enddo
      endif

c      print*, PVModel%Param
c           TEMPORARY WORK
c           convert the best model to a TR model
      call ParamToTiltedRing(PVModel,ModelTiltedRing
     &          ,TR_FittingOptions)

c           BUG FIX (Dan, 2026-08-18): ParamToTiltedRing only refreshes
c               SigUse; the Sigma field that particle-flux generation
c               actually reads (Ring_CalcParticleFlux_Basic) is derived
c               from SigUse via this same Linear_Log_SDSwitch check, but
c               that conversion otherwise only happens inline inside
c               TiltedRingModelComparison (FullModelComparison.f), which
c               never runs again after this point in the fit. Without
c               this, ModelTiltedRing%Sigma is left holding whatever
c               amoeba's last internal trial evaluation (not necessarily
c               the converged vertex -- Nelder-Mead's last function call
c               is often a rejected reflection/contraction/shrink probe)
c               happened to set it to, so the final output model cube
c               (built in FitOutput.f's OutputCube, which assumes
c               ModelTiltedRing is already correct and does no further
c               conversion of its own) silently uses a stale, wrong
c               surface density for every ring. Confirmed via a
c               Fortran-vs-JS "bisection paradox" investigation: particle
c               positions/velocities matched bit-exactly end to end, yet
c               every particle's flux differed by a large, ring-specific
c               (non-constant) factor -- traced to exactly this gap.
      if(PFlags%Linear_Log_SDSwitch .eq. 0) then
        ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma=
     &      ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%SigUse
      elseif(PFlags%Linear_Log_SDSwitch .eq. 1) then
        ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%Sigma=
     &   10.**(ModelTiltedRing%R(0:ModelTiltedRing%nRings-1)%SigUse)
      endif

c       Deallocate the model vector at the end
c      call DeAllocateParamVector(PVModel)
      DEALLOCATE(chiArray)
      DEALLOCATE(paramGuesses)
      return
      end subroutine
ccccccccc

ccccccc
c       RunObjectiveProbe: reads probe_params.txt (format below) and
c       evaluates TiltedRingModelComparison once per probe vector, with
c       idum reset to a fixed seed before each call so results are directly
c       comparable to a paired JS-side evaluation using the same seed.
c
c       probe_params.txt format:
c         line 1: nProbes
c         line 2: seed (integer, negative -- reset before EVERY probe)
c         line 3..3+nProbes-1: nParams space-separated values per line
c
c       Output: one line per probe, "PROBE_RESULT <index> <chi2>"
      subroutine RunObjectiveProbe(nParams)
      use PipelineGlobals
      implicit none
      integer,INTENT(IN) :: nParams
      integer nProbes,seed,p,k,ios
      real,ALLOCATABLE :: probeParams(:)
      real chi2

      ALLOCATE(probeParams(nParams))

      open(unit=77,file='probe_params.txt',status='old'
     &          ,action='read',IOSTAT=ios)
      if (ios .ne. 0) then
        print*, "RunObjectiveProbe: could not open probe_params.txt"
        return
      endif

      read(77,*) nProbes
      read(77,*) seed

      do p=1,nProbes
        read(77,*) (probeParams(k),k=1,nParams)
        idum=seed
        call TiltedRingModelComparison(probeParams,chi2)
        print*, "PROBE_RESULT",p,chi2
      enddo

      close(77)
      DEALLOCATE(probeParams)

      return
      end subroutine
ccccccccc



ccccccc
c
      subroutine DownhillSimplexRun(nParams,paramGuesses,chiArray)
      use PipelineGlobals
      implicit none
      integer,INTENT(IN) :: nParams
      real,INTENT(INOUT):: paramGuesses(nParams+1,nParams)
      real,INTENT(INOUT) :: chiArray(nParams+1)
      integer i,iter
      real chi2
      logical FitFlag


      print*, "ftol ChecK", ftol, ObservedDC%DH%Uncertainty

c           Make the set of initial guesses
c      call MakeParamGuessArray(PVModel,ParamGuesses
c     &          ,PVIni%nParams,idum,IniGuessWidth)

c       Get the goodness of fit for each guess
      chiArray=0.
      do i=1,PVModel%nParams+1
        PVModel%Param(0:PVModel%nParams-1)=
     &                  ParamGuesses(i,1:PVModel%nParams)
        call TiltedRingModelComparison(PVModel%Param,chi2)
        chiArray(i)=chi2
        print*, PID,i,paramGuesses(i,:),chi2
      enddo
      block
        character(64) EnvVal6
        integer EnvLen6, SXUnit
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal6,
     &            EnvLen6)
        if (EnvLen6 .gt. 0) then
          open(newunit=SXUnit, file="SimplexTraceF.txt",
     &        status="unknown", position="append")
          do i=1,PVModel%nParams+1
            write(SXUnit,'(A,I0,A,ES27.19)') "vertex",i-1," ",
     &          DBLE(chiArray(i))
          enddo
          write(SXUnit,'(A)') "---"
          close(SXUnit)
        endif
      end block

c       Run the downhill simplex
      call amoeba(paramGuesses,chiArray
     &                  ,PVModel%nParams+1
     &                  ,PVModel%nParams
     &                  ,PVModel%nParams,ftol
     &                  ,TiltedRingModelComparison,iter,PID
     &                  ,FitFlag)

c           NOTE: DownhillSimplexRun is shared by BOTH passes -- this
c               fires twice per fit (pass1 then pass2); use call ORDER
c               (2nd = pass2) to tell them apart, there's no per-pass
c               label available at this shared subroutine's scope.
      print '(A,I6,1X,L1,1X,ES17.9)', "ITER_F",iter,FitFlag
     &          ,chiArray(1)

c       Store the best model in the PVModel object
      PVModel%BestLike=chiArray(1)
      PVModel%Param(0:PVModel%nParams-1)=
     &              ParamGuesses(1,1:PVModel%nParams)

c       One-off diagnostic (Fortran-vs-JS model-cube divergence
c           investigation): dump the converged parameter vector at full
c           REAL(4) precision, unambiguously labelled so it can be grepped
c           straight out of the trace log and fed into the JS side's own
c           model-cube resynthesis for a direct apples-to-apples compare.
      if (TraceSwitch .eq. 1) then
        print*, "CONVERGED_VECTOR", PVModel%Param(0:PVModel%nParams-1)
     &          ,PVModel%BestLike
      endif
      block
        character(64) EnvValA
        integer EnvLenA, PVUnit, PVj
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvValA,
     &            EnvLenA)
        if (EnvLenA .gt. 0) then
          open(newunit=PVUnit, file="PvBestTraceF.txt",
     &        status="unknown", position="append")
          do PVj=0,PVModel%nParams-1
            write(PVUnit,'(A,I0,A,ES27.19)') "p",PVj," ",
     &          DBLE(PVModel%Param(PVj))
          enddo
          write(PVUnit,'(A)') "---"
          close(PVUnit)
        endif
      end block

      return
      end subroutine
cccccccc

ccccccc
      subroutine MakeParamGuessArray(PredictedPV,ParamGuesses
     &                      ,ndim,idum,lambda,StrictEstimate)

      use BasicRanNumGen
      implicit none

      integer idum
      integer ndim
      Type(ParameterVector) PredictedPV
      real ParamGuesses(ndim+1,ndim)
      real lambda,lambdaPar

      integer i,j,k
      integer counter
      integer StrictEstimate,ParCounter,AcceptBeyondLimits

      print*, "Param Limits Check",idum
      do i=0, ndim-1
        print*, i, PredictedPV%ParamLowerLims(i)
     &          ,PredictedPV%ParamUpperLims(i)
     &          ,PredictedPV%Param(i)
     *          ,PredictedPV%ParamRange(i)
      enddo

      counter=0
      print*, "Param Guess Array Creation"
      do k=1,ndim+1
        if(k .eq. 1) then
            ParamGuesses(k,1:ndim)=PredictedPV%Param(0:ndim-1)
        else
            do i=1,ndim
                ParCounter=0
                AcceptBeyondLimits=0
                j=i-1
c100             lambdaPar=lambda*(PredictedPV%ParamUpperLims(j)
c     &                      -PredictedPV%ParamLowerLims(j))
 100            lambdaPar=lambda*PredictedPV%ParamRange(j)
                lambdaPar=(2*ran2(idum)-1.)*lambdaPar
                ParamGuesses(k,i)=PredictedPV%Param(j)+lambdaPar
                counter=counter+1
                ParCounter=ParCounter+1
c                   If the strict estimate is relaxed, accept the parameter regardless of the limits after 200 tries
                if(StrictEstimate .eq. 0) then
                    if(ParCounter .ge. 200) AcceptBeyondLimits=1
                endif
c       print*, i,ParamGuesses(k,i),PredictedPV%ParamLowerLims(j)
c     &                  ,PredictedPV%ParamUpperLims(j),counter
c     &              ,j,PredictedPV%CyclicSwitch(j)
c               For cyclic parameters, adjust them to be inside the range
200             continue
                if(PredictedPV%CyclicSwitch(j) .eq. 1) then
                    if (ParamGuesses(k,i) .lt.
     &                   PredictedPV%ParamLowerLims(j)) then
                            ParamGuesses(k,i)=ParamGuesses(k,i)
     &                          +PredictedPV%ParamUpperLims(j)
                        goto 200
                    elseif(ParamGuesses(k,i) .gt.
     &                   PredictedPV%ParamUpperLims(j)) then
                            ParamGuesses(k,i)=ParamGuesses(k,i)
     &                          -PredictedPV%ParamUpperLims(j)
                        goto 200
                    endif
                endif

                if(AcceptBeyondLimits .eq. 0) then
                    if(ParamGuesses(k,i) .lt.
     &                   PredictedPV%ParamLowerLims(j)) goto 100
                    if(ParamGuesses(k,i) .gt.
     &                   PredictedPV%ParamUpperLims(j)) goto 100
                endif
            enddo
        endif
c        print*, k,ParamGuesses(k,:)
      enddo
      print*, "Done param guess array creation"

      return
      end subroutine
cccccc

      end module
