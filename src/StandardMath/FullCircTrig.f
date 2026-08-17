c       This file contains routines that are useful for calculating inverse
c           trig functions over a full circle

      module FullCircTrig
      use CommonConsts
      contains
c

ccccccc
c       The full circle arc-tangent function
      subroutine FullCircATan(X,Y,Theta)
      implicit none
      real, INTENT(IN) :: X,Y
      real, INTENT(OUT) :: Theta

      real fd_atan
      external fd_atan

c           BUG FIX: native atan() was not guaranteed bit-identical to
c               JS's Math.atan() (confirmed: real disagreements for
c               ordinary inputs) -- unlike sin/cos, V8's Math.atan is NOT
c               already a bit-exact fdlibm derivative. Switched to a new
c               fdlibm atan port (fd_atan/fdAtan) verified bit-exact
c               against each other across 500k random values (Dan, 2026).
      Theta=fd_atan(Y/X)

      if(X .ge. 0. .and. Y .ge. 0.) then
        Theta=Theta
      elseif(X .lt. 0. .and. Y .ge. 0.) then
        Theta=Theta+Pi
      elseif(X .lt. 0. .and. Y .lt. 0.) then
        Theta=Theta+Pi
      elseif(X .ge. 0. .and. Y .lt. 0.) then
        Theta=Theta+2.*Pi
      endif


      return
      end subroutine
cccccccc

      end module

