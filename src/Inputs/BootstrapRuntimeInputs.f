cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting general
c       inputs for the pipeline
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module BootstrapInputMod

      use BootstrapGlobals
      use DataCubeInputMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating bootstrap samples
      subroutine BootstrapIn()
      implicit none
      character(500) TopLevelInputFile
      integer timearray(3)
 

      print*, "Getting the bootstrap generator inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Bootstrap Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the name of the observed cube
      read(10,*)
      read(10,'(A)') ObservedCubeFile
c           Get the name of the model cube file
      read(10,*)
      read(10,'(A)') ModelCubeFile
c           Get the name of the cube mask file
c      read(10,*)
c      read(10,'(A)') MaskCubeFile
c           Get the base name for the output files
      read(10,*)
      read(10,'(A)') BaseOutName

c       Get the size of the blocks in terms of beams and channels
      read(10,*)
c      read(10,*) SpatialBlackSize,VelBlockSize
      read(10,*) VelBlockSize
c       Get the center of the cube for geometry based resampling
      read(10,*)
      read(10,*) BS_Cent%CentX,BS_Cent%CentY,BS_Cent%CentV
     &              ,BS_Cent%PA
      block
        character(64) EnvVal9
        integer EnvLen9, BCUnit
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal9,
     &            EnvLen9)
        if (EnvLen9 .gt. 0) then
          open(newunit=BCUnit, file="BsCentTraceF.txt",
     &        status="unknown", position="append")
          write(BCUnit,'(A,ES27.19)') "CentX ", DBLE(BS_Cent%CentX)
          write(BCUnit,'(A,ES27.19)') "CentY ", DBLE(BS_Cent%CentY)
          write(BCUnit,'(A,ES27.19)') "CentV ", DBLE(BS_Cent%CentV)
          write(BCUnit,'(A,ES27.19)') "PA ", DBLE(BS_Cent%PA)
          close(BCUnit)
        endif
      end block
c           BS_Cent%PA is already written in radians by MakeBootstrapSample.
c           WriteBootstrapFile -- do not convert again here (was a real,
c           confirmed double degrees->radians conversion bug).

c           Get the base name for the outputs
c      read(10,*)
c      read(10,'(A)') BaseOutName
c           Get the base name of the objects
c      read(10,*)
c      read(10,'(A)') OutputFolder
c       Get the number of bootstrap samples to make
c      read(10,*)
c      read(10,*) nBootstrap
c       Get the noise in mJy
c      read(10,*)
c      read(10,*) RMS
c       Get the random seed
      read(10,*)
      read(10,*) idum


      close(10)

c           If the seed is positive, set it by the time array
      if(idum .ge. 0) then
        print*, "using time to generate random seed"
        call itime(timeArray)
        idum=abs(timeArray(1)*idum)
        idum=-int(idum*timeArray(2))-timeArray(3)
      endif


      return
      end subroutine
cccccccc

ccccccc
      subroutine LoadCubesForBootstrap()
      implicit none
      integer MaskSwitch
      Type(Beam2D) TempBeam

      MaskSwitch=0

      print*, "Loading in data cube"
      call ReadFullDataCube(ObservedCube,ObservedBeam
     &                      ,ObservedCubeFile,MaskSwitch)


c      print*, "Loading in model cube"
      call ReadFullDataCube(ModelCube,ObservedBeam
     &                      ,ModelCubeFile,MaskSwitch)
      block
        character(64) EnvValB
        integer EnvLenB, VXUnit, VXi, VXj, VXk
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvValB,
     &            EnvLenB)
        if (EnvLenB .gt. 0) then
          open(newunit=VXUnit, file="VoxelTraceF.txt",
     &        status="unknown", position="append")
          write(VXUnit,'(A,ES27.19)') "obs_29_11_2 ",
     &        DBLE(ObservedCube%Flux(29,11,2))
          write(VXUnit,'(A,ES27.19)') "model_29_11_2 ",
     &        DBLE(ModelCube%Flux(29,11,2))
          write(VXUnit,'(A,ES27.19)') "diff_29_11_2 ",
     &        DBLE(ObservedCube%Flux(29,11,2)
     &            -ModelCube%Flux(29,11,2))
          close(VXUnit)
          open(newunit=VXUnit, file="ModelFullF.txt",
     &        status="unknown")
          do VXk=0,ModelCube%DH%nChannels-1
            do VXj=0,ModelCube%DH%nPixels(1)-1
              do VXi=0,ModelCube%DH%nPixels(0)-1
                write(VXUnit,'(ES17.9)')
     &              DBLE(ModelCube%Flux(VXi,VXj,VXk))
              enddo
            enddo
          enddo
          close(VXUnit)
        endif
      end block
c      call Allocate_Beam2D(ObservedBeam,ObservedCube%DH%nPixels)
c      call DCBrightnessConversion(ObservedCube,ObservedBeam)
c      call DCBrightnessConversion(ModelCube,ObservedBeam)

c      call ReadFullDataCube(MaskCube,TempBeam
c     &                      ,MaskCubeFile,MaskSwitch)
      return
      end subroutine
cccccccc


      end module BootstrapInputMod
