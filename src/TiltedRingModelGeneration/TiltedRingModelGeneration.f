cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       build a tilted ring model full of particles.
c          The routines assume the tilted ring has been allocated
c           and the individual ring parameters have been set.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingGenerationMod
      use TiltedRingMod
      use SingleRingGenerationMod
      use DataCubeMod
      use BeamMod

      implicit none
      contains

ccccc
      subroutine BuildTiltedRingModel(TR,idum,Noise,DC,BUse)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      integer,INTENT(INOUT) :: idum   !A seed for the random number generation
      Type(DataCube), INTENT(IN) :: DC
      Type(Beam2D), INTENT(IN) :: BUse
      real, INTENT(IN) :: Noise

      integer i
      real AvgChanPerPix

c      print*, "Building Tilted Ring"
c      print*, "Surf Dens Sanity", TR%R(1)%Sigma

      do i=0,TR%nRings-1
c           BUse's beam major axis and the ring's own VRot/VDisp/Rmid set
c               how many spectral channels this ring's flux is smeared
c               across -- feeds Ring_CalcNumParticles below so fast/wide
c               rings get proportionally more particles (upstream commit
c               76ade48, "Changed how the number of particles in each ring
c               is calculated").
        call CalcAvgChanPerPix(i,TR,DC,BUse,AvgChanPerPix)
        call Ring_CalcNumParticles(TR%R(i)        !Get the # of particles in the ring
     &              , TR%cmode,TR%CloudBaseSurfDens,Noise
     &              , AvgChanPerPix)
        call Ring_ParticleAllocation(TR%R(i))               !Allocate the particle array
        call Ring_ParticleGeneration(TR%R(i),idum)          !Generate the particles in the ring
      enddo

      return
      end subroutine
cccccc

cccccc
c       CalcAvgChanPerPix: how many spectral channels a ring's own
c           velocity spread (rotation + dispersion) covers, relative to
c           how many beam-widths wide the ring is -- rings with a steep
c           velocity gradient or wide dispersion smear their flux across
c           more channels per spatial pixel, so need proportionally more
c           particles to sample adequately. Floored at DDisp+1 so a ring
c           with ~zero rotation still gets at least dispersion-width
c           coverage.
      subroutine CalcAvgChanPerPix(ringIndx,TR,DC
     &              ,BUse,AvgChanPerPix)
      implicit none
      integer, INTENT(IN) :: ringIndx
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      Type(DataCube), INTENT(IN) :: DC
      Type(Beam2D), INTENT(IN) :: BUse
      real, INTENT(INOUT) :: AvgChanPerPix

      real DV, DR,DDisp

      DV=2.*TR%R(ringIndx)%VRot/abs(DC%DH%ChannelSize)
      DDisp=2.*sqrt(2.)*TR%R(ringIndx)%VDisp/abs(DC%DH%ChannelSize)
      DR=2.*TR%R(ringIndx)%Rmid

      AvgChanPerPix=(BUse%BeamMajorAxis/DR)*(DV+DDisp)
      if (AvgChanPerPix .lt. DDisp+1) then
        AvgChanPerPix=DDisp+1.
      endif

      return
      end subroutine
c
cccccccc


      end module
