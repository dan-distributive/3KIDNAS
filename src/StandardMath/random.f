
      module BasicRanNumGen
c           gasdev's cached-spare-value state, moved from a function-local
c               SAVE variable to module level (2026-08-18, Dan) so
c               ResetGasdevCache below can clear it from outside gasdev --
c               needed for the Fortran-vs-JS idum-divergence proof
c               diagnostic (MaybeOverrideIdum in FitOutput.f): idum only
c               stays negative until the FIRST ran2() call inside it, which
c               happens (via position generation) BEFORE gasdev is ever
c               reached, so gasdev's own `if (idum.lt.0) iset=0` never
c               actually fires except on the true first-ever call of the
c               whole program run -- meaning a mid-run idum reset alone
c               does NOT reset this cache the way a fresh JS rng object
c               does. Purely a location change -- identical behavior for
c               every existing call site.
      integer :: GasdevIset = 0
      real :: GasdevGset = 0.0
c           ran2's own shuffle-table state (idum2, iv(32), iy), moved from
c               function-local SAVE to module level for the SAME reason
c               and by the SAME pattern as GasdevIset/GasdevGset above
c               (Dan probe, 2026-09-27): a bare idum integer is NOT
c               ran2's full state -- transplanting only idum into a fresh
c               JS rng object (which never re-triggers its own `idum<=0`
c               reinit for a positive value) leaves iv/iy/idum2 at their
c               construction-time defaults, disconnected from any real
c               evolved stream. GetRan2State/SetRan2State below let a
c               caller (MaybeOverrideIdum) read/write the FULL state for
c               a genuine apples-to-apples transplant test against JS.
      integer :: Ran2Idum2 = 123456789
      integer :: Ran2Iv(32) = 0
      integer :: Ran2Iy = 0
      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran0
c
c     This generates a random number using the Park and Miller
c     method of Press et al.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran0(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IA,IM,IQ,IR,MASK
c
      real AM
      parameter(IA=16897,IM=2147483647, AM=1./IM)
      parameter(IQ=127773, IR=2836, MASK=123459876)

      integer k
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran1
c
c     This generates a random number using the Park and Miller
c     method with a Bays-Durham shuffle of Press et al.
c
c     Needs a negative idum to initialize
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran1(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IA,IM,IQ,IR,NTAB,NDIV
c
      real AM,EPS,RNMX
      parameter(IA=16807,IM=2147483647, AM=1./IM)
      parameter(IQ=127773, IR=2836, NTAB=32,NDIV=1+(IM-1)/NTAB)
      parameter(EPS=1.2e-7,RNMX=1.-EPS)

      integer j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv/NTAB*0/, iy/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum.le.0 .or. iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if(idum.lt.0) idum=idum+IM
            if(j .le.NTAB) iv(j)=idum
 11      enddo
         iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum .lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      if(ran1 .lt. 0.) print*, 'bug in ran1'
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran2
c
c    This function generates a random number using L'Ecuyer method
c     from Press et al.
c
c     Needs a negative idum to initialize
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran2(idum)
c     
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer IM1, IM2, IMM1, IA1, IA2, IQ1
      integer IQ2,IR1,IR2,NTAB,NDIV
c
      real AM,EPS,RNMX
c
      parameter(IM1=2147483563,IM2=2147483399,AM=1./IM1)
      parameter(IMM1=IM1-1)
      parameter(IA1=40014,IA2=40692, IQ1=53668, IQ2=52774,IR1=12211)
      parameter(IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB)
      parameter(EPS=1.2e-7)
      parameter(RNMX=1.-EPS)
c
      integer idum2,j,k,iv(NTAB),iy
c           Local copies of the module-level state, not SAVE'd here anymore
c               -- see the module header comment on Ran2Idum2/Ran2Iv/Ran2Iy
c               for why. Copy in, run the UNCHANGED algorithm below, copy
c               back out before every return. Purely a storage-location
c               change, identical in spirit to the GasdevIset move above.
      idum2=Ran2Idum2
      iv=Ran2Iv
      iy=Ran2Iy
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum .le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 12 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0) idum=idum+IM1
            if(j .le. NTAB) iv(j)=idum
 12      enddo
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ2)-k*IR1
      if(idum .lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if(idum2 .lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy .lt.1) iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      if(ran2 .le. 0.) print*, 'bug in ran2'
      Ran2Idum2=idum2
      Ran2Iv=iv
      Ran2Iy=iy
      return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ran3
c
c     Another method of generating random numbers from Press et al.
c
c     Needs a negative idum to initialize
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function ran3(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      integer MBIG,MSEED,MZ
      real FAC
c
      parameter(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
c
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(56)
      SAVE iff,inext,inextp,ma
      DATA iff/0/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=abs(MSEED-abs(idum))
         mj=mod(mj,MBIG)
         ma(55)=mj
         mk=1
         do 13 i=1,54
            ii=mod(21*i,56)
            ma(ii)=mk
            mk=mj-mk
            if(mk .lt. MZ) mk=mk+MBIG
            mj=ma(ii)
 13      enddo
         do 14 k=1,4
            do 15 i=1, 56
               ma(i)=ma(i)-ma(1+mod(i+30,56))
               if(ma(i) .lt. MZ) ma(i)=ma(i)+MBIG
 15         enddo
 14      enddo
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext .eq. 56) inext=1
      inextp=inextp+1
      if(inextp .eq. 56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj .lt. MZ) mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      end function
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Gasdev
c
c     This function calculates a gaussian deviate with a zero mean
c     and a sigma of 1.
c
c     It can be found in Press et al.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function gasdev(idum)
c
      implicit none
c
      integer, INTENT(INOUT) :: idum
      real fac, rsq, v1,v2
      real fd_log
      external fd_log
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (idum.lt. 0) GasdevIset=0
      if(GasdevIset .eq. 0) then
c           BUG FIX (Dan, 2026): X**2. (real exponent literal) is not
c               guaranteed bit-identical to X*X in gfortran -- confirmed
c               directly (see PhysCoordTransform.f's matching fix/comment
c               for the reproduction). The JS port already used v1*v1/
c               v2*v2 (no real**real ambiguity to port), so this was a
c               one-sided Fortran divergence risk in the RNG's own hot
c               path. Verbatim upstream code (byte-identical to
c               NateDeg/3KIDNAS) -- reported upstream, also fixed here.
 1      v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1*v1+v2*v2
        if(rsq .ge. 1. .or. rsq .eq. 0.) goto 1
        fac=sqrt(-2.*fd_log(rsq)/rsq)
        GasdevGset=v1*fac
        gasdev=v2*fac
        GasdevIset=1
      else
        gasdev=GasdevGset
        GasdevIset=0
      endif
      return
      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       ResetGasdevCache: one-off diagnostic (Fortran-vs-JS
c           idum-divergence proof, Dan 2026-08-18) -- forces GasdevIset
c           back to 0, so a caller can guarantee gasdev's next call
c           computes a fresh pair instead of returning a stale cached
c           spare left over from earlier in the program's run. See this
c           module's header comment for why a mid-run idum reset alone
c           doesn't achieve this.
      subroutine ResetGasdevCache()
      implicit none
      GasdevIset=0
      return
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       GetRan2State/SetRan2State: one-off diagnostic (Dan probe,
c           2026-09-27) -- read/write ran2's FULL shuffle-table state
c           (idum2, iv(32), iy), not just the bare idum scalar. See this
c           module's header comment on Ran2Idum2/Ran2Iv/Ran2Iy for why a
c           bare idum transplant into JS is insufficient (JS's own ran2
c           only reinitializes iv/iy/idum2 when idum<=0, so a positive
c           override idum leaves them at construction-time defaults,
c           disconnected from any real evolved stream).
      subroutine GetRan2State(idum2Out,ivOut,iyOut)
      implicit none
      integer, INTENT(OUT) :: idum2Out,ivOut(32),iyOut
      idum2Out=Ran2Idum2
      ivOut=Ran2Iv
      iyOut=Ran2Iy
      return
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module

