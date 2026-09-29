cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CalcBeamKernelMod
      use BeamMod
      use CommonConsts


      implicit none

      contains

ccccccc
      subroutine Calculate2DBeamKernel(B,PixelSizes)
      implicit none
      Type(Beam2d), INTENT(INOUT):: B
      real,INTENT(IN):: PixelSizes(0:1)

      integer i, j, nTot
      real x,y,R2
      real xp,yp,cpa,spa
      real cellArea
      real fd_cos, fd_sin, fd_exp
      external fd_cos, fd_sin, fd_exp
c      print*, "About to calculate kernel"
c      print*, "Calculating Kernel", B%BeamSigmaVector,B%nRadialCells

      cpa=fd_cos(-B%BeamSigmaVector(2))     !Will need to use a negative rotation for the transform
      spa=fd_sin(-B%BeamSigmaVector(2))     ! to the beam's axes
      do i=-B%nRadialCells,B%nRadialCells
        do j=-B%nRadialCells,B%nRadialCells
c               Get the positions in x and y for the particle pixel
c            x=real(i)*PixelSizes(0)
c            y=real(j)*PixelSizes(1)
            x=real(i)
            y=real(j)
c               Now rotate these to the major axis
            xp=x*cpa-y*spa
            yp=x*spa+y*cpa
c               Get the normalized squared radius
c           BUG FIX (Dan, 2026): X**2. not guaranteed bit-identical to
c               X*X in gfortran -- see PhysCoordTransform.f's matching
c               fix/comment. Verbatim upstream code -- reported upstream,
c               also fixed here.
            R2=(xp/B%BeamSigmaVector(0))*(xp/B%BeamSigmaVector(0))
     &              +(yp/B%BeamSigmaVector(1))*(yp/B%BeamSigmaVector(1))
c               Calculate the kernel value for this cell
            B%Kernel(i,j)=1./sqrt(2.*Pi*B%BeamSigmaVector(0)
     &                  *B%BeamSigmaVector(1))
     &                  *fd_exp(-R2/2)
        enddo
      enddo
c           Diagnostic (Dan probe, 2026-09-28): raw pre-normalization
c               kernel center value + sum, to isolate whether the
c               Gaussian evaluation itself or the sum()/renormalization
c               step is where Fortran and JS first diverge.
      block
        character(64) EnvVal3
        integer EnvLen3
        integer KTUnit2
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal3,
     &            EnvLen3)
        if (EnvLen3 .gt. 0) then
          open(newunit=KTUnit2, file="KernelTraceF.txt",
     &        status="unknown", position="append")
          write(KTUnit2,'(A,ES27.19)') "rawcenter ",
     &        DBLE(B%Kernel(0,0))
          write(KTUnit2,'(A,ES27.19)') "rawsum ",
     &        DBLE(sum(B%Kernel))
          close(KTUnit2)
        endif
      end block
c           Renormalize so that the sum of the kernel is 1.
c      cellArea=PixelSizes(0)*PixelSizes(1)
c      B%Kernel=B%Kernel/(sum(B%Kernel)*cellArea)
c      print*, sum(B%kernel)
      B%Kernel=B%Kernel/(sum(B%Kernel))
c      print*, sum(B%kernel)

c           One-off diagnostic (Dan probe, 2026-09-28): dump the
c               real-space kernel's checksum + corner/center sample
c               values to a dedicated file (not stdout -- bootstrap
c               sub-process stdout redirection makes prints unreliable
c               to capture) to check whether the REAL-SPACE Gaussian
c               kernel (before any FFT) already differs between
c               platforms -- the forward r2c FFT itself is now proven
c               bit-exact in isolation, so if this kernel differs,
c               that's upstream of the FFT, in the fd_cos/fd_sin/fd_exp
c               construction.
      block
        character(64) EnvVal2
        integer EnvLen2
        integer KTUnit
        call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal2,
     &            EnvLen2)
        if (EnvLen2 .gt. 0) then
          open(newunit=KTUnit, file="KernelTraceF.txt",
     &        status="unknown", position="append")
          write(KTUnit,'(A,ES27.19)') "pa ",
     &        DBLE(B%BeamSigmaVector(2))
          write(KTUnit,'(A,ES27.19)') "sigma0 ",
     &        DBLE(B%BeamSigmaVector(0))
          write(KTUnit,'(A,ES27.19)') "sigma1 ",
     &        DBLE(B%BeamSigmaVector(1))
          write(KTUnit,'(A,ES27.19)') "cpa ", DBLE(cpa)
          write(KTUnit,'(A,ES27.19)') "spa ", DBLE(spa)
          write(KTUnit,'(A,ES27.19)') "sum ", DBLE(sum(B%Kernel))
          write(KTUnit,'(A,ES27.19)') "center ",
     &        DBLE(B%Kernel(0,0))
          write(KTUnit,'(A,ES27.19)') "corner ",
     &        DBLE(B%Kernel(-B%nRadialCells,-B%nRadialCells))
          write(KTUnit,'(A,ES27.19)') "r1c1 ",
     &        DBLE(B%Kernel(1,1))
          write(KTUnit,'(A,ES27.19)') "r2c3 ",
     &        DBLE(B%Kernel(2,3))
          write(KTUnit,'(A,I0)') "nRadialCells ",
     &        B%nRadialCells
          close(KTUnit)
        endif
      end block

      return
      end subroutine
cccccccc


cccccccc
      subroutine CalculateComplex2DKernel(B)
      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f'

      Type(Beam2d) B
      double precision,ALLOCATABLE :: PaddedKernel(:,:)
      double precision,ALLOCATABLE::WrappedPaddedKernel(:,:)
      integer i,j,k,l,CentKernel(2)
      integer*8 ArrPlan_r2c


c       Allocate  padded and padded & wrapped arrays
      ALLOCATE(PaddedKernel(B%PaddedSize(1),B%PaddedSize(2)))
      ALLOCATE(WrappedPaddedKernel(B%PaddedSize(1),B%PaddedSize(2)))

c       Set up padded kernel
      PaddedKernel=0.
      do i=1, 2*B%nRadialCells+1
        do j=1, 2*B%nRadialCells+1
            k=i-B%nRadialCells-1
            l=j-B%nRadialCells-1
            PaddedKernel(i,j)=B%Kernel(k,l)
c            print*, i,j,k,l
        enddo
      enddo

c           Wrap the padded kernel
      CentKernel=(2*B%nRadialCells+1)/2
      call MakeWrappedArray(B%PaddedSize,CentKernel
     &              ,PaddedKernel,WrappedPaddedKernel)      !In this file.

c           Make the fftw plan to get the complex kernel
      call dfftw_plan_dft_r2c_2d(ArrPlan_r2c,B%PaddedSize(1)
     &                      ,B%PaddedSize(2)
     &                      ,PaddedKernel
     &                      ,B%ComplexKernel,FFTW_ESTIMATE
     &                      ,FFTW_PRESERVE_INPUT)           !FFTW3 Routine
c           Get the complex kernel
      call dfftw_execute_dft_r2c(ArrPlan_r2c, WrappedPaddedKernel
     &                          , B%ComplexKernel)          !FFTW3 Routine

c       Get rid of the fftw plan
      call dfftw_destroy_plan(ArrPlan_r2c)

c       Deallocate the unneeded arrays
      DEALLOCATE(PaddedKernel)
      DEALLOCATE(WrappedPaddedKernel)

c       Note that the complex kernel has been created
      B%ComplexKernelCreated=.True.

      return
      end subroutine

cccccccccc




ccccccc
c           Wrap a double precision array such that the central value is at (1,1) in the new array
      subroutine MakeWrappedArray(SA,CV,Arr,WrappedArr)
      implicit none

      integer,INTENT(IN) :: SA(2),CV(2)     !The size of the array and the 'central values'
      double precision,INTENT(IN) :: Arr(SA(1),SA(2))
      double precision,INTENT(INOUT) :: WrappedArr(SA(1),SA(2))

      integer i,j,k,l

c           Loop over all cells
      do i=1, SA(1)
        do j=1, SA(2)
c               Get the indices of the wrapped array
            k=i-CV(1)
            l=j-CV(2)
c               Wrap in x if k<=0
            if(k .le. 0) then
                k=SA(1)+k
            endif
c               Wrap in y if l<=0
            if(l .le. 0) then
                l=SA(2)+l
            endif
c            print*, i,j,k,l
            WrappedArr(k,l)=Arr(i,j)
        enddo
      enddo


      return
      end subroutine
ccccccccccc


      end module
