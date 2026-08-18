cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SingleRingGenerationMod
      use ParticleMod
      use TiltedRingMod


      implicit none

      contains

ccccccc
c       RoundForParticleStability: masks off the low 12 bits of x's
c           float32 mantissa (keeping the top 11 of 23 bits -- ~0.05%
c           relative precision) before Ring_CalcNumParticles truncates it
c           to an integer particle count.
c
c           Root cause this exists for: R%Sigma (and everything derived
c           from it, e.g. DensMultiplications below) can differ between
c           Fortran and the JS port by a handful of float32 ULPs (~5e-7
c           relative) -- an already-known, unavoidable cross-platform
c           floating-point difference (compiler codegen, not a bug --
c           see UPSTREAM_SYNC.md) that is normally harmless everywhere
c           else in this pipeline, since it stays a smoothly tiny
c           difference in whatever value carries it. But
c           Ring_CalcNumParticles's int(...)+1 is a hard truncation, and a
c           continuous value that happens to land within noise-distance of
c           an integer boundary truncates to *different whole integers* on
c           the two platforms. That single-particle difference then
c           consumes a different number of ran2()/gasdev() draws,
c           permanently desyncing idum for the rest of the fit --
c           confirmed directly: Fortran and the JS port ran in bit-exact
c           lockstep (identical idum every evaluation) until this exact
c           failure mode first triggered, then diverged permanently and
c           the JS optimizer never recovered.
c
c           Masking gives ~900x safety margin over the observed ~5e-7
c           relative noise, at the cost of ~0.05% precision on a Monte
c           Carlo sampling-density parameter -- scientifically
c           irrelevant, and it must be applied identically on both
c           platforms (see hexF32/roundForParticleStability in
c           SingleRingGeneration.js) for the masked value itself to still
c           agree.
      real function RoundForParticleStability(x)
      implicit none
      real, INTENT(IN) :: x
      integer(4) ix
      integer(4), parameter :: MASK = int(z'FFFFF000',4)
      ix = transfer(x, 0)
      ix = iand(ix, MASK)
      RoundForParticleStability = transfer(ix, 0.0)
      return
      end function
cccccccc

ccccccc
      subroutine Ring_CalcNumParticles(R,cmode,CloudSurfDens,Noise
     &          ,AvgChannelsPerPix)
      use CommonConsts
      use PipelineGlobals
      implicit none
      Type(Ring), INTENT(INOUT):: R
      integer,INTENT(IN) ::cmode
      real, INTENT(IN) :: CloudSurfDens
      real,INTENT(IN) :: Noise
      real,INTENT(IN) :: AvgChannelsPerPix
      real Pixel_Ring
      real DensMultiplications
      real Rl,Rh

c           Get the rough area of the ring in pixels
c      Pixel_Ring=2.*Pi*R%Rmid*R%Rwidth /pixelarea      !Original in 3DBarolo GalMod
      Rl=R%Rmid-R%Rwidth/2.
      Rh=R%Rmid+R%Rwidth/2.
c      Pixel_Ring=Pi*(Rh**2.-Rl**2)/pixelarea
      Pixel_Ring=Pi*(Rh**2.-Rl**2)
c       Calculate a term to get the rough number of clouds per pixel area.  It is normalized by
c           the surface density so that each particle has roughly the
c           same amount of flux.
c
c           DELIBERATELY on the pre-upstream-commit-76ade48 formula for
c           now (not using Noise/AvgChannelsPerPix below) -- Nathan's own
c           attached example config used cdens=10, but his email text
c           says the actual default is cdens=100 (10 looks like a typo in
c           the example, not the intended default); revisit once that's
c           confirmed and re-verify Fortran/JS agreement before switching
c           this back on. Noise/AvgChannelsPerPix stay threaded through
c           Ring_CalcNumParticles's signature (unused here) so re-enabling
c           is a one-line change, not a re-plumb -- see
c           TiltedRingModelGeneration.f's CalcAvgChanPerPix, already
c           wired and confirmed bit-matching between Fortran and the JS
c           port.
      DensMultiplications=CloudSurfDens
     &               *((R%Sigma)**real(cmode))
      R%nParticles=int(RoundForParticleStability(DensMultiplications
     &          *Pixel_Ring))+1
      if (TraceSwitch .eq. 1) then
        print*, "NPTRACE",R%Rmid,R%Sigma,DensMultiplications,
     &          Pixel_Ring,AvgChannelsPerPix,R%nParticles
        print '(A,Z8.8,1X,Z8.8,1X,Z8.8,1X,Z8.8,1X,Z8.8)',
     &      "NPHEX",transfer(R%Sigma,0),transfer(Noise,0),
     &      transfer(DensMultiplications,0),
     &      transfer(Pixel_Ring,0),transfer(AvgChannelsPerPix,0)
      endif
c      print*, "Single Ring PA",R%PositionAngle*180./Pi
c      print*, "Single Ring Center",R%CentPos
c      print*, "Number of particles", R%nParticles,R%Sigma
c      print*, R%Rwidth,R%Rmid,R%Sigma,CloudSurfDens,cmode
c     &          ,Pixel_Ring,DensMultiplications, R%PositionAngle*180./Pi

      return
      end subroutine
cccccccc

ccccccc
      subroutine Ring_ParticleGeneration(R,idum)
      use BasicRanNumGen
      use CommonConsts
      use PipelineGlobals, only: TraceSwitch
      implicit none
      Type(Ring), INTENT(INOUT) :: R     !The ring that we are generating particles for
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      integer i
      real Rmin,Rmax, Area       !The minimum and maximum radii of the ring and the ring area
      character(64) BisectEnvVal
      integer BisectEnvLen

      call get_environment_variable("TRACE_OVERRIDE_IDUM",
     &          BisectEnvVal,BisectEnvLen)

c           To randomly sample the ring area we need Rmin and Rmax
      Rmin=R%Rmid-R%Rwidth/2.       !We need Rmin and Rmax for the ring first
      Rmax=R%Rmid+R%Rwidth/2.
      Area=Pi*(Rmax**2.-Rmin**2)
c      print*, "Single Ring Area Check", Area,Rmin,Rmax
c           One-off diagnostic ("bisection paradox" investigation, Dan
c               2026-08-18): BINTRACE showed every particle landing in the
c               same cell on both platforms but with a different constant
c               Flux value for ring 3 -- since Flux=Sigma*Area/nParticles
c               and RINGTRACE only ever dumps ring 0's fields, this prints
c               every ring's own Rmid/Rwidth/Sigma/nParticles/Area so the
c               actual numeric input responsible can be identified
c               directly, gated on TRACE_OVERRIDE_IDUM like the others.
      if (BisectEnvLen .gt. 0) then
        print '(A,ES25.17,1X,ES25.17,1X,ES25.17,1X,I12,1X,ES25.17)',
     &      'RINGGEOMTRACE', R%Rmid, R%Rwidth, R%Sigma,
     &      R%nParticles, Area
      endif


c      print*, "Ring Params", R%CentPos,R%Inclination, R%PositionAngle
c     &          ,R%VSys, R%VRot, R%VRad, R%VDisp,R%Vvert
c     &          ,R%dvdz,R%Sigma,R%z0,R%zGradiantStart
c      print*, "Consistency Check", R%VSys,R%VRad
c           Now go through all the particles
      do i=0, R%nParticles-1
        call Ring_ParticlePosSelect(R,idum, Rmin,Rmax,i)
        call ParticlePosProject(R%P(i), R%Inclination, R%PositionAngle)
        call ParticlePos_NewCenter(R%P(i), R%CentPos)
        call Ring_CalcParticle_VSys(R,i,idum)
        call Ring_CalcParticleFlux_Basic(R,i, Area)
        if (TraceSwitch.eq.1 .and. i.lt.5) then
          print*, "PARTTRACE",R%Rmid,i,idum,R%P(i)%Pos(0),
     &            R%P(i)%Pos(1),R%P(i)%Pos(2),
     &            R%P(i)%ProjectedPos(0),R%P(i)%ProjectedPos(1),
     &            R%P(i)%ProjectedVel(2)
        endif
c           One-off diagnostic (Fortran-vs-JS gasdev-desync bisection,
c               Dan 2026-08-18): checkpoint every 200th particle (plus the
c               very last one) across the FULL ring, not just the first 5
c               PARTTRACE covers -- to localize exactly where a desync
c               (re)starts, if one does, deeper into a ring's particle
c               loop. Gated on TRACE_OVERRIDE_IDUM.
        if (BisectEnvLen .gt. 0 .and.
     &          (mod(i,200).eq.0 .or. i.eq.R%nParticles-1)) then
          print '(A,F10.6,1X,I6,1X,I12,1X,F14.8)', 'BISECTTRACE',
     &        R%Rmid,i,idum,R%P(i)%ProjectedVel(2)
        endif

c        print*, "Random Check", R%P(i)%AngPos, R%P(i)%Pos
c        print*, "Projected Random Check", R%P(i)%ProjectedPos(0:1)
c        print*, "PRojected Velocity", R%P(i)%ProjectedVel(2)
      enddo

c      print*, "particle Flux ChecK", sum(R%P(0:R%nParticles-1)%Flux)


      return
      end subroutine
cccccccc


ccccccc
      subroutine Ring_ParticlePosSelect(R,idum, Rmin, Rmax, PartID )
      use BasicRanNumGen
      use CommonConsts
      implicit none
      Type(Ring), INTENT(INOUT) :: R     !The ring that we are generating particles for
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      integer, INTENT(IN) :: PartID
      real, INTENT(IN) :: Rmin,Rmax
      real RR, Theta, Z     !Temporary cylindrical coordinates of each particle
      real fd_sin, fd_cos, fd_atanh
      external fd_sin, fd_cos, fd_atanh

c           First get a random radius -- however this is not uniform in R.  We want
c           equal area sampling we'll be using a sqrt distribution.
      RR=ran2(idum)*(Rmax**2.-Rmin**2.)             !Sample a uniform squared radius in the correct range
      RR=sqrt(RR+Rmin**2.)          !Add on the minimum radius squared and take the sqrt
c           Next get a random angle
      Theta=ran2(idum)*2.*Pi
c           Finally get a height from a sech^2 distribution (use atanh on a random distrubtion and scale
c                   by the height
      Z=fd_atanh((2.*ran2(idum)-1.))*R%z0
c      print*, "Z", Z,R%z0

c           Store the cylindrical coordinates
      R%P(PartID)%AngPos(0)=RR
      R%P(PartID)%AngPos(1)=Theta
      R%P(PartID)%AngPos(2)= Z
c           Store the Cartesian Coordinates
      R%P(PartID)%Pos(0)=RR*fd_cos(Theta)
      R%P(PartID)%Pos(1)=RR*fd_sin(Theta)
      R%P(PartID)%Pos(2)= Z

      return
      end subroutine
ccccccccc

cccccccc
      subroutine ParticlePosProject(P,Inclination,PositionAngle)
      implicit none
      Type(Particle), INTENT(INOUT) :: P
      real, INTENT(IN) :: Inclination, PositionAngle
      real cpa,spa,XTemp,YTemp
      real fd_sin, fd_cos
      external fd_sin, fd_cos
c       First get the inclined positions in the major axis frame (aligned with X-axis)
c           in the temporary arrays
      XTemp=P%Pos(0)
      YTemp=P%Pos(1)*fd_cos(Inclination)-P%Pos(2)*fd_sin(Inclination)
c       Now rotate the positions by the position angle
      cpa=fd_cos(PositionAngle)        !Get the cos and sin values
      spa=fd_sin(PositionAngle)
c       Calculate the coordinates in the projected array
      P%ProjectedPos(0)=XTemp*cpa-YTemp*spa
      P%ProjectedPos(1)=XTemp*spa+YTemp*cpa

c      print*, P%Pos(0:1), P%AngPos(0:2)
c     &          , P%ProjectedPos(0:1)
c      print*, P%ProjectedPos, PositionAngle*180./3.14, XTemp,YTemp
      return
      end subroutine
cccccccccccc

cccccccc
      subroutine ParticlePos_NewCenter(P,NewCent)
      implicit none
      Type(Particle), INTENT(INOUT) :: P
      real, INTENT(IN) :: NewCent(0:1)

c       Shift the projected position to lie about the new center
c      print*, "Initial Pos", P%ProjectedPos,NewCent
      P%ProjectedPos(0)=P%ProjectedPos(0)+NewCent(0)
      P%ProjectedPos(1)=P%ProjectedPos(1)+NewCent(1)
c      print*, "Final Pos", P%ProjectedPos
      return
      end subroutine
cccccccccccc

cccccc
      subroutine Ring_CalcParticle_VSys(R,PartID,idum)
      use BasicRanNumGen
      implicit none
      Type(Ring), INTENT(INOUT) :: R
      integer, INTENT(IN) :: PartID
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      real VRotP,VRadP,VVertP,VDispP
      real V_FromRotation, V_FromRadial,V_FromVertical
      real delR,delZ
      real cTheta,sTheta
      integer DerivativeSide
      real fd_sin, fd_cos
      external fd_sin, fd_cos

c           First get the particle's rotational, radial, and vertical motions
c           at the current radial positions

c               To do the linear interpolation, figure out which side you're on
c      delR=(R%P(PartID)%AngPos(0)-R%Rmid)
c      if(delR .le. 0) then
c        DerivativeSide=0
c      else
c        DerivativeSide=1
c      endif

c      delZ=abs(R%P(PartID)%AngPos(2))-R%zGradiantStart
      VRotP=R%VRot !+ R%dVRot_dR(DerivativeSide)*delR
      VRadP=R%VRad !+ R%dVRad_dR(DerivativeSide)*delR
      VDispP=R%VDisp !+ R%dVDisp_dR(DerivativeSide)*delR
c           Check if the particle is far enough above the plane that V_rot is affected
      if(abs(R%P(PartID)%AngPos(2)) .gt. R%zGradiantStart) then
        VRotP=VRotP-R%dvdz*delZ
      endif
c
c           Now get the sin and cos of the angle
c
      cTheta=fd_cos(R%P(PartID)%AngPos(1))
      sTheta=fd_sin(R%P(PartID)%AngPos(1))
c
c           Next get the observed velocity due to the main components
c
      V_FromRotation=VRotP*cTheta*fd_sin(R%Inclination)
      V_FromRadial=VRadP*sTheta*fd_sin(R%Inclination)
      V_FromVertical=R%Vvert*fd_cos(R%Inclination)         !Note the Vvert is the global vertical motion of a ring
c        Sum up the various components
      R%P(PartID)%ProjectedVel(2)=R%VSys+V_FromRotation
     &                      +V_FromRadial+V_FromVertical
c
c           Finally get a random velocity from the dispersion and add it in to the projected velocity
c
      R%P(PartID)%ProjectedVel(2)=R%P(PartID)%ProjectedVel(2)
     &                      +gasdev(idum)*VDispP

c      print*, "hmmm", VDispP, R%VDisp

      return
      end subroutine
ccccccccc

cccccc
      subroutine Ring_CalcParticleFlux_Basic(R,PartID,Area)
      implicit none
      Type(Ring), INTENT(INOUT) :: R
      integer, INTENT(IN) :: PartID
      real,INTENT(IN):: Area
      real SigmaP, del R
      integer DerivativeSide
c           Get the local Sigma
      SigmaP=R%Sigma

      R%P(PartID)%Flux=SigmaP*Area/real(R%nParticles)


      return
      end subroutine
cccccccc

      end module
