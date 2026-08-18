cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for
c       calculating the difference in flux between 2 cubes
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CubeDiffMod
      use DataCubeMod




      contains
cccccccc
c
      subroutine ConstructDiffCube(Cube1,Cube2,DiffCube)
      implicit none
      Type(DataCube),INTENT(IN) :: Cube1,Cube2
      Type(DataCube),INTENT(INOUT):: DiffCube
      integer i, j, k

c       First check that the difference cube hasn't been allocated
      if(ALLOCATED(DiffCube%Flux)) then
c           If it is allocated, deallocate it
        call DeAllocateDataCube(DiffCube)
      endif
c       Now set the cube header to match the initial cube
      DiffCube%DH=Cube1%DH
c           And allocate the difference
      call AllocateDataCube(DiffCube)
c       Finally get the difference between the cubes
      do i=0, Cube1%DH%nPixels(0)-1
        do j=0,Cube1%DH%nPixels(1)-1
            do k=0, Cube1%DH%nChannels-1
                DiffCube%Flux(i,j,k)=Cube1%Flux(i,j,k)-Cube2%Flux(i,j,k)
                if ((i.eq.14 .or. i.eq.15) .and.
     &              (j.eq.23 .or. j.eq.24) .and.
     &              (k.eq.89 .or. k.eq.90)) then
                  print '(A,I3,I3,I3,1X,3F12.6)', 'SRCTRACE',i,j,k,
     &              Cube1%Flux(i,j,k),Cube2%Flux(i,j,k),
     &              DiffCube%Flux(i,j,k)
                endif
            enddo
        enddo
      enddo



      return
      end subroutine

cccccc

      end module
