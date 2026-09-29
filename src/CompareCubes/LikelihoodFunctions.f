cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for
c       calculating the likelihood of some model values
c       being drawn from a set of observations with a given uncertainty
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module LikelihoodMod

      implicit none

      PROCEDURE(LikelihoodFuncInterface),POINTER :: LikePoint =>null()

      ABSTRACT INTERFACE
        subroutine LikelihoodFuncInterface(Likelihood
     &              ,nElements,Model,Obs,Uncertainties)
            IMPLICIT NONE
            real, INTENT(OUT) :: Likelihood
            integer, INTENT(IN) :: nElements
            real, INTENT(IN) :: Model(0:nElements),Obs(0:nElements)
            real, INTENT(IN) :: Uncertainties(0:nElements)

        END subroutine LikelihoodFuncInterface
      END INTERFACE



      contains
cccccccc
c           This routine calculates the tradition chi^2 statistic
      subroutine Chi2Calc(chi2,nElements,Model,Obs,Uncertainties)
      implicit none
      real, INTENT(OUT) :: chi2
      integer, INTENT(IN) :: nElements
      real, INTENT(IN) :: Model(0:nElements),Obs(0:nElements)
      real, INTENT(IN) :: Uncertainties(0:nElements)
      integer i

c      print*, "Calculating chi^2",nElements,Uncertainties(0)

c           BUG FIX (Dan, 2026): X**2. (real exponent literal) is not
c               guaranteed bit-identical to X*X in gfortran -- see
c               PhysCoordTransform.f's matching fix/comment for the
c               reproduction. This is the main per-pixel chi^2 objective
c               function, called on every single optimizer evaluation --
c               the hottest path in the whole fit. JS's port already used
c               direct multiplication. Verbatim upstream code -- reported
c               upstream, also fixed here.
      chi2=0.
      do i=0,nElements-1
        chi2=chi2+(Model(i)-Obs(i))*(Model(i)-Obs(i))
     &          /(Uncertainties(i)*Uncertainties(i))
c        print*, "chi2 sum",i, chi2,Model(i),Obs(i),Uncertainties(i)
      enddo

      return
      end subroutine
cccccc



cccccccc
c
      subroutine LogChi2Calc(chi2,nElements,Model,Obs,Uncertainties)
      implicit none
      real, INTENT(OUT) :: chi2
      integer, INTENT(IN) :: nElements
      real, INTENT(IN) :: Model(0:nElements),Obs(0:nElements)
      real, INTENT(IN) :: Uncertainties(0:nElements)
      integer i

c      print*, "Calculating log chi^2"
      call Chi2Calc(chi2,nElements,Model,Obs,Uncertainties)
      chi2=log10(chi2)
      return
      end subroutine
cccccc


cccccccc
c           This routine calculates the tradition chi^2 statistic
      subroutine Chi2Calc_logElements(chi2,nElements,Model
     &                  ,Obs,Uncertainties)
      implicit none
      real, INTENT(OUT) :: chi2
      integer, INTENT(IN) :: nElements
      real, INTENT(IN) :: Model(0:nElements),Obs(0:nElements)
      real, INTENT(IN) :: Uncertainties(0:nElements)
      real,parameter :: Small=1.e-10
      real M,O

      integer i

c      print*, "Calculating chi^2"

      chi2=0.
      do i=0,nElements-1
        M=Model(i)
        O=Obs(i)
        if(M .le. Small) M=Small
        if(O .le. Small) O=Small

c           BUG FIX (Dan, 2026): X**2. not guaranteed bit-identical to
c               X*X in gfortran -- see PhysCoordTransform.f's matching
c               fix/comment. Verbatim upstream code -- reported upstream,
c               also fixed here.
        chi2=chi2+(log10(M)-log10(O))*(log10(M)-log10(O))
     &                  /(Uncertainties(i)*Uncertainties(i))
c        print*, "chi2Elements sum", i, chi2, M, O
c     &              , log10(M), log10(O), Uncertainties(i)
      enddo
c      print*, "Like Sanity", sum(Model),sum(Obs), Uncertainties(0)

      return
      end subroutine
cccccc

      end module
