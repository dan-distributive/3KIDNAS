cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       convolve some 2D array with a kernel
c
c           The routines assume that the complex kernel has
c           already been calculated
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TwoDConvolutionMod




      implicit none

      contains

ccccccc
      subroutine Convolve2D(Arr,SizeArray,ComplexKernel,SC,SizePad
     &          ,ConvolvedArray)

      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f'

      integer,INTENT(IN) :: SizeArray(2)
      integer,INTENT(IN) :: SizePad(2)

      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))
      real,INTENT(OUT) :: ConvolvedArray(SizeArray(1),SizeArray(2))

      double precision,ALLOCATABLE :: PaddedArray(:,:),PaddedKernel(:,:)
      double precision,ALLOCATABLE :: RealConvolve(:,:)
      double precision, ALLOCATABLE :: WrappedPaddedKernel(:,:)

      integer i,j,k,l

      integer,INTENT(IN) :: SC(2)       !The size of the complex arrays
c      integer SizeComplex(2)
      double complex,ALLOCATABLE :: ComplexArr(:,:)
      double complex, INTENT(IN) :: ComplexKernel(SC(1),SC(2))
      double complex,ALLOCATABLE :: ComplexConvolve(:,:)

c           One-off diagnostic (Dan probe, 2026-09-27): dump the forward
c               FFT's raw double-precision complex output, first call
c               only, so it can be diffed bin-for-bin against JS's own
c               rdft2R2cSync output for the SAME channel -- isolates
c               whether the forward r2c FFT itself already differs
c               (native vs wasm FFTW codelet selection) before the
c               kernel multiply or inverse FFT ever run.
      integer, SAVE :: ConvolveCallCount = 0
      character(64) EnvVal
      integer EnvLen
      character(8192) WisdomStrBuf
      integer WisdomBufLen

      integer*8 ArrPlan_r2c,ArrPlan_c2r
ccccc
c           Allocate Arrays
      call AllocateConvolutionArrays(SizeArray,Arr
     &              ,SizePad,SC,PaddedArray
     &              ,RealConvolve
     &              ,ComplexArr
     &              ,ComplexConvolve)           !In this file

      call SetupPaddedArrays(SizeArray,SizePad
     &                      ,Arr
     &                      ,PaddedArray)       !In this file


c       Make the real-to-complex fftw plan
      call dfftw_plan_dft_r2c_2d(ArrPlan_r2c,SizePad(1),SizePad(2)
     &                      ,PaddedArray
     &                      ,ComplexArr,FFTW_ESTIMATE
     &                      ,FFTW_PRESERVE_INPUT)               !FFTW3 Routine
c       Make the complex-to-real fftw plan
c           BUG FIX (Dan, 2026): dropped FFTW_PRESERVE_INPUT here to match
c               the WASM build's c2r_2d plan, which cannot honor it for
c               this transform size (fftw_c2r_2d_wasm's own comment: plan
c               creation returns NULL the moment the flag is added, for
c               this exact size, under this build's reduced generic
c               codelet set -- forced to plain FFTW_ESTIMATE there
c               already). The two sides were asking FFTW for different
c               guarantees on the same nominal transform, letting the
c               planner pick different internal algorithms even though
c               both call fftw_plan_dft_c2r_2d. Safe here: the input array
c               (ComplexConvolve, passed at execute time below -- note
c               it's NOT the same array the plan was created with, a
c               separate FFTW "new-array execute") is never read again
c               after dfftw_execute_dft_c2r, only deallocated, so nothing
c               depends on FFTW leaving it undestroyed.
      call dfftw_plan_dft_c2r_2d(ArrPlan_c2r,SizePad(1),SizePad(2)
     &              ,ComplexArr,PaddedArray
     &              ,FFTW_ESTIMATE)         !FFTW3 Routine

c           Do the fft transform of the padded image
      call dfftw_execute_dft_r2c(ArrPlan_r2c, PaddedArray, ComplexArr)  !FFTW3 Routine
      ConvolveCallCount=ConvolveCallCount+1
      call get_environment_variable("TRACE_DUMP_PRECONV",EnvVal,
     &          EnvLen)
      if (EnvLen .gt. 0 .and. ConvolveCallCount .eq. 27) then
        print*, "PLANPRINT SizePad", SizePad(1), SizePad(2)
        call dfftw_print_plan(ArrPlan_r2c)
        print*, ""
        WisdomBufLen = 8192
        call fftw_export_wisdom_cstr(WisdomStrBuf, WisdomBufLen)
        print*, "WISDOMSTART"
        print*, trim(WisdomStrBuf)
        print*, "WISDOMEND"
        print '(A,2ES27.19)', "FFTFORWARDTRACE bin(1,1)",
     &      DBLE(ComplexArr(1,1)),DIMAG(ComplexArr(1,1))
        print '(A,2ES27.19)', "FFTFORWARDTRACE bin(2,1)",
     &      DBLE(ComplexArr(2,1)),DIMAG(ComplexArr(2,1))
        print '(A,2ES27.19)', "FFTFORWARDTRACE bin(1,2)",
     &      DBLE(ComplexArr(1,2)),DIMAG(ComplexArr(1,2))
        print '(A,2ES27.19)', "FFTFORWARDTRACE bin(5,7)",
     &      DBLE(ComplexArr(5,7)),DIMAG(ComplexArr(5,7))
        print '(A,2ES27.19)', "FFTFORWARDTRACE sum",
     &      DBLE(sum(ComplexArr)),DIMAG(sum(ComplexArr))
      endif

c           Convolve the transformed arrays
      do i=1, SizePad(1)/2+1
        do j=1, SizePad(2)
            ComplexConvolve(i,j)=ComplexArr(i,j)*ComplexKernel(i,j)
        enddo
      enddo

c           Do the backwards fft transform fo the convolved array
      call dfftw_execute_dft_c2r(ArrPlan_c2r, ComplexConvolve
     &                  ,RealConvolve)                          !FFTW3 Routine
c           Normalize the real array
      RealConvolve=RealConvolve/SizePad(1)/SizePad(2)

c           Set the convolved array (of the correct size) to the real-convolve array
      do i=1, SizeArray(1)
        do j=1, SizeArray(2)
            ConvolvedArray(i,j)=RealConvolve(i,j)
        enddo
      enddo
c           Destroy the fftw plans
      call dfftw_destroy_plan(ArrPlan_r2c)      !FFTW3 Routine
      call dfftw_destroy_plan(ArrPlan_c2r)      !FFTW3 Routine

c       Free up the padded array space
      DEALLOCATE(PaddedArray)
      DEALLOCATE(ComplexArr)
      DEALLOCATE(ComplexConvolve)
      DEALLOCATE(RealConvolve)

      return
      end subroutine
cccccccc


cccccc
c           This function allocates all the different arrays needed
      subroutine AllocateConvolutionArrays(SizeArray
     &              ,Arr
     &              ,SizePad,SizeComplex,PaddedArray
     &              ,RealConvolve
     &              ,ComplexArr
     &              ,ComplexConvolve)
      implicit none
      integer,INTENT(IN) :: SizeArray(2), SizeComplex(2)
      integer,INTENT(IN) :: SizePad(2)

      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))

      double precision,ALLOCATABLE,INTENT(OUT) :: PaddedArray(:,:)
      double precision,ALLOCATABLE,INTENT(OUT) :: RealConvolve(:,:)

      double complex,ALLOCATABLE,INTENT(OUT) :: ComplexArr(:,:)
      double complex,ALLOCATABLE,INTENT(OUT) :: ComplexConvolve(:,:)


c           Allocate the real/double precision arrays
      ALLOCATE(PaddedArray(SizePad(1),SizePad(2)))
      ALLOCATE(RealConvolve(SizePad(1),SizePad(2)))
c           Allocate the complex arrays
      ALLOCATE(ComplexArr(SizeComplex(1),SizeComplex(2)))
      ALLOCATE(ComplexConvolve(SizeComplex(1),SizeComplex(2)))

      return
      end subroutine
cccccccc

ccccccccc
c           This routine initializes the padded arrays
      subroutine SetupPaddedArrays(SizeArray,SP
     &                      ,Arr
     &                      ,PaddedArray)
      implicit none

      integer, INTENT(IN) :: SizeArray(2), SP(2)
      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))
      double precision,INTENT(INOUT)::PaddedArray(SP(1),SP(2))

      integer i, j
c          Make the padded array
      PaddedArray=0.
      do i=1, SizeArray(1)
        do j=1, SizeArray(2)
            PaddedArray(i,j)=Arr(i,j)
        enddo
      enddo

      return
      end subroutine
cccccccc



      end module
