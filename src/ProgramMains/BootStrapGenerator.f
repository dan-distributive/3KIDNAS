ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This is the main routine for a code to generate multiple bootstrap
c       samples from some input cube.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program CubeBootStrapGenerator
      use BootstrapInputMod
      use BootstrapGlobals
      use DataCubeInputMod
      use GenBootstrapMod
      use DataCubeMod
      use DataCubeOutputsMod
      use CalcBeamKernelMod
      use BeamMod

      use CubeDiffMod
      use DataCubeOutputsMod

      use GenBootstrapMod

      use FlippingBootstrapMod


      implicit none

      character(500) OutName
      integer i

      real IncTemp,PATemp,Center(2),VSysTemp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "Boot Strap Sampler"
c           Get the runtime inputs
      call BootstrapIn()
c       Load in the relevant cubes
      call LoadCubesForBootstrap()

c       Make a bootstrap sample
c      call GenBootstrapSample()
      call GenFlipBootstrapSample()

      block
        character(64) EnvVal8
        integer EnvLen8, RSUnit, RFi, RFj, RFk
        double precision RSSum
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal8,
     &            EnvLen8)
        if (EnvLen8 .gt. 0) then
          RSSum=sum(dble(BootstrapCube%Flux),
     &        MASK=(BootstrapCube%Flux.eq.BootstrapCube%Flux))
          open(newunit=RSUnit, file="ResampleTraceF.txt",
     &        status="unknown", position="append")
          write(RSUnit,'(A,ES27.19)') "sum ", RSSum
          write(RSUnit,'(A,ES27.19)') "px000 ",
     &        DBLE(BootstrapCube%Flux(0,0,0))
          write(RSUnit,'(A,ES27.19)') "pxlast ",
     &        DBLE(BootstrapCube%Flux(BootstrapCube%DH%nPixels(0)-1,
     &            BootstrapCube%DH%nPixels(1)-1,
     &            BootstrapCube%DH%nChannels-1))
          close(RSUnit)
          open(newunit=RSUnit, file="ResampleFullF.txt",
     &        status="unknown")
          do RFk=0,BootstrapCube%DH%nChannels-1
            do RFj=0,BootstrapCube%DH%nPixels(1)-1
              do RFi=0,BootstrapCube%DH%nPixels(0)-1
                write(RSUnit,'(ES17.9)')
     &              DBLE(BootstrapCube%Flux(RFi,RFj,RFk))
              enddo
            enddo
          enddo
          close(RSUnit)
        endif
      end block


c       Output the resampled cube
      BootstrapCube%DH%PixelSize=
     &              BootstrapCube%DH%PixelSize/3600.
c      BootstrapCube%DH%RefLocation(0:2)=
c     &                  BootstrapCube%DH%RefLocation(0:2)+1

      SampleFile=trim(BaseOutName)//".fits"
      call WriteDataCubeToFITS(BootstrapCube
     &          ,ObservedBeam,SampleFile
     &          ,trim(BaseOutName))


      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



