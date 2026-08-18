cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       build a tilted ring model full of particles.
c          The routines assume the tilted ring has been allocated
c           and the individual ring parameters have been set.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FillDataCubeWithTiltedRingMod
      use TiltedRingMod
      use DataCubeMod
      use ParticleMod


      implicit none
      contains

ccccc
      subroutine FillDataCubeWithTiltedRing(DC,TR)
      implicit none
      Type(DataCube), INTENT(INOUT) :: DC
      Type(TiltedRingModel), INTENT(IN) :: TR

      integer i,j,k,CellIndex(0:2)
      logical InBounds

      Type(Particle) PTest

      logical, save :: HeaderTraced = .false.
      character(64) EnvVal
      integer EnvLen
      character(64) BinEnvVal
      integer BinEnvLen

c       One-off diagnostic (Fortran-vs-JS model-cube divergence
c           investigation, Dan 2026-08-17): dump the model cube's header
c           fields the very first time this runs, so they can be diffed
c           directly against the JS port's equivalent -- per-particle
c           channel binning below depends directly on DC%DH%Start(2) and
c           DC%DH%ChannelSize, so if these don't actually match between
c           platforms, every particle's channel assignment would shift
c           systematically. Gated on TRACE_OVERRIDE_IDUM being set so it
c           only fires during a deliberate controlled-comparison run, not
c           every normal evaluation.
      if (.not. HeaderTraced) then
        call get_environment_variable("TRACE_OVERRIDE_IDUM",EnvVal,
     &          EnvLen)
        if (EnvLen .gt. 0) then
          print '(A,3I5,1X,I5,1X,3F16.6,1X,3F16.6,1X,F16.6)',
     &      'HEADERTRACE',DC%DH%nPixels(0),DC%DH%nPixels(1),
     &      DC%DH%nChannels,0,
     &      DC%DH%Start(0),DC%DH%Start(1),DC%DH%Start(2),
     &      DC%DH%RefVal(0),DC%DH%RefVal(1),DC%DH%RefVal(2),
     &      DC%DH%ChannelSize
          HeaderTraced = .true.
        endif
      endif

c      print*, "Filling DataCube"
c           Make sure the initial cube is empty
c      print*, "Cube center", DC%DH%Start(0:2)
      DC%Flux=0.
c           One-off diagnostic (Fortran-vs-JS "bisection paradox"
c               investigation, Dan 2026-08-18): particle generation is now
c               proven bit-exact end to end (BISECTTRACE), yet the
c               pre-convolution cube still diverges (PRECONVTRACE). That
c               narrows the residual to this binning step -- either the
c               per-particle Flux value (never directly compared before
c               now) or the RoundForBinStability-protected cell-index
c               decision. BINTRACE prints every particle that lands in one
c               of the already-tracked pixels (x in {14,15}, y in
c               {23,24}, ch in {89,90}), with its flux and the running
c               cell total after the add, so the two platforms' logs can
c               be diffed particle-by-particle. Gated on
c               TRACE_OVERRIDE_IDUM; fires every evaluation like
c               RINGTRACE/PRECONVTRACE -- grab the last block.
      call get_environment_variable("TRACE_OVERRIDE_IDUM",BinEnvVal,
     &          BinEnvLen)
c           Loop through all particles
      do i=0,TR%nRings-1
c        print*, "Filling data cube with ring", i
c     &          , sum(TR%R(i)%P(0:TR%R(i)%nParticles-1)%Flux)
        do j=0, TR%R(i)%nParticles-1
c                   Get the index for the cell where the particle lives
            call FindParticleCellLocation(TR%R(i)%P(j),DC,CellIndex)
c                   Check that the particle is inside the cube
            call CheckIfInCube(CellIndex,DC,InBounds)
c                   If it is in the cube, add the flux to the proper cell
            if(InBounds) then
                DC%Flux(CellIndex(0),CellIndex(1),CellIndex(2))=
     &                  DC%Flux(CellIndex(0),CellIndex(1),CellIndex(2))+
     &                  TR%R(i)%P(j)%Flux
                if (BinEnvLen .gt. 0) then
                  if ((CellIndex(0).eq.14 .or. CellIndex(0).eq.15) .and.
     &                (CellIndex(1).eq.23 .or. CellIndex(1).eq.24) .and.
     &                (CellIndex(2).eq.89 .or. CellIndex(2).eq.90)) then
                    print '(A,1X,I2,1X,I6,1X,3I5,1X,F14.8,1X,F14.8)',
     &                  'BINTRACE', i, j, CellIndex(0), CellIndex(1),
     &                  CellIndex(2), TR%R(i)%P(j)%Flux,
     &                  DC%Flux(CellIndex(0),CellIndex(1),CellIndex(2))
                  endif
                endif
            endif
        enddo
      enddo

c      print*, "Filled data cube flux check", sum(DC%Flux)
c      print*, "1st part flux", TR%R(0)%P(0)%Flux

c      print*, "TR Center Position Check"
c      PTest%ProjectedPos(0:1)=TR%R(0)%CentPos(0:1)
c      PTest%ProjectedVel(2)=TR%R(0)%VSys
c      print*, "Center Position",PTest%ProjectedPos(0:1)
c     &          ,PTest%ProjectedVel(2)
c      call FindParticleCellLocation(PTest,DC,CellIndex)
c      print*, "Cell Index Calc", CellIndex
c      print*, "X Position Check", PTest%ProjectedPos(0)
c     &          , DC%DH%Start(0), DC%DH%PixelSize(0)
c     &          ,PTest%ProjectedPos(0)-DC%DH%Start(0)
c     & ,(PTest%ProjectedPos(0)-DC%DH%Start(0))/DC%DH%PixelSize(0)
c     &          ,PTest%ProjectedPos(0)/3600., DC%DH%Start(0)/3600.

c      print*, "Y Position Check", PTest%ProjectedPos(1)
c     &          , DC%DH%Start(1), DC%DH%PixelSize(1)
c     &          ,PTest%ProjectedPos(1)-DC%DH%Start(1)
c     & ,(PTest%ProjectedPos(1)-DC%DH%Start(1))/DC%DH%PixelSize(1)
c     &          ,PTest%ProjectedPos(1)/3600., DC%DH%Start(1)/3600.
      

      return
      end subroutine
cccccc

cccccc
c       RoundForBinStability: same technique as SingleRingGeneration.f's
c           RoundForParticleStability / GenerateBootstrap.f's
c           RoundForInterpStability (masks off the low 12 bits of x's
c           float32 mantissa, keeping the top 11 of 23 bits -- ~0.05%
c           relative precision), applied here to a per-particle
c           nearest-pixel rounding decision instead.
c
c           Root cause: P%ProjectedPos/P%ProjectedVel (the product of a
c           chain of float32 trig/rotation arithmetic during particle
c           generation) can differ from the JS port by a handful of
c           float32 ULPs -- the same already-known, unavoidable
c           cross-platform difference documented at
c           RoundForParticleStability. Normally harmless, but
c           CellIndex=int(Pos+0.5) is a hard round-to-nearest-pixel
c           decision, applied independently to EVERY particle (tens of
c           thousands per ring): whenever one particle's position sits
c           within that noise-distance of a half-integer boundary,
c           Fortran and the JS port round it into different, adjacent
c           cells. Since each particle's full flux is added (not
c           interpolated) to whichever single cell it lands in, one
c           straddling particle yanks a chunk of flux out of one cell and
c           into its neighbor on only one platform -- confirmed directly
c           (2026-08-17): with particle COUNTS already verified identical
c           per ring between two independent Fortran/JS initial fits, and
c           the model cube isolated as the only diverging input (observed
c           cube confirmed bit-identical) to a downstream bootstrap-
c           resampling pixel diff, this per-particle binning round is the
c           remaining discrete decision capable of producing that
c           divergence at this scale.
      real function RoundForBinStability(x)
      implicit none
      real, INTENT(IN) :: x
      integer(4) ix
      integer(4), parameter :: MASK = int(z'FFFFF000',4)
      ix = transfer(x, 0)
      ix = iand(ix, MASK)
      RoundForBinStability = transfer(ix, 0.0)
      return
      end function
cccccccc

cccccc
      subroutine FindParticleCellLocation(P,DC,CellIndex)
      implicit none
      Type(Particle), INTENT(IN) :: P
      Type(DataCube), INTENT(IN) :: DC
      integer, INTENT(OUT) :: CellIndex(0:2)
      integer i,j

      do i=0,1
c        CellIndex(i)=int((P%ProjectedPos(i)-DC%DH%Start(i))
c     &                  /(DC%DH%PixelSize(i)) )
        CellIndex(i)=int(RoundForBinStability(
     &                  (P%ProjectedPos(i))+0.5))
      enddo
      CellIndex(2)=int(RoundForBinStability(
     &                  (P%ProjectedVel(2)-DC%DH%Start(2))
     &                  /(DC%DH%ChannelSize) +0.5))

      i=2
c      if(CellIndex(i) .eq. 100) then
c        print*, P%ProjectedVel(2), DC%DH%Start(2)
c     &          ,DC%DH%ChannelSize, DC%DH%RefLocation(2)
c     &          ,DC%DH%RefVal(2)
c     &          ,DC%DH%Start(2)+100.*DC%DH%ChannelSize
c     &          ,DC%DH%Start(2)+99.*DC%DH%ChannelSize
cc     &          ,DC%Channels(CellIndex(i)-2:CellIndex(i)+1)
c      endif
c      print*, "hmmm", CellIndex,DC%DH%Start
c     &          ,P%ProjectedPos(0:1),P%ProjectedVel(2)


      return
      end subroutine
ccccccc

ccccccc
      subroutine CheckIfInCube(CellIndex,DC,InBounds)
      implicit none
      Type(DataCube), INTENT(IN) :: DC
      integer, INTENT(IN) :: CellIndex(0:2)
      logical, INTENT(OUT) :: InBounds
      integer i

      InBounds=.True.
c           Check the spatial dimensions
      do i=0,1
        if(CellIndex(i) .lt. 0 .or.
     &                  CellIndex(i) .ge. DC%DH%nPixels(i)) then
            InBounds=.False.
        endif
      enddo
c           Check the velocity
      if(CellIndex(2) .lt. 0 .or.
     &                  CellIndex(2) .ge. DC%DH%nChannels) then
        InBounds=.False.
      endif
c      print*, "Inbounds", InBounds, CellIndex
c     &              ,DC%DH%nPixels, DC%DH%nChannels


      return
      end subroutine
cccccccc


      end module
