*******************************************************
* FILE: Geometry.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Geometry
            implicit none

      contains
* BasicGrid


      subroutine Grid2DUni(Lth,Np,Ni)

        use CustomDouble

* Grid2DUni calculates the node positions in a bi-dimensional setting
* Uni means mesh is uniformly spaced
* For this assignment, symmetry was applied. So mesh is defined only to first quadrant (both X and Y coordinates are positive)

*     N: number of nodes in each direction
*        1,2 indices refer to x and y respectively

*     Lth: Lenght along each direction
*        1,2 indices refer to x and y respectively

* Np: array containing node positions in X axis (to be calculated)
* Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*     The first two indices locate the node indexes  - I for x axis and J for y axis
*     The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

* Ni: array containing interface positions in X axis (to be calculated)
* Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*     The first two indices locate the node indexes  - I for x axis and J for y axis
*     The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

* Del: array containing distance between neighboring nodes in all axis(to be calculated)
* Index ranges from 1 to N(1) - 1 in X direction, from 1 to N(2) - 1 in Y direction.
*     The first two indices locate the node indexes  - I for x axis and J for y axis
*     The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

* Delta: array containing distance between parallel interfaces of a volume element in all axis (to be calculated)
* Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*     The first two indices locate the node indexes  - I for x axis and J for y axis
*     The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I

        integer :: Nx, Ny
        real(dp) :: Lth(1:2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2),DX,DY

        real(dp) :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        Nx = N(1)
        Ny = N(2)

*       Mesh for X direction
        ! Symmetry is being considered here. Only a quarter of domain is being modelled.
        if (Nx .GT. 1) then
          DX = (Lth(1) / 2.0_dp) / float(Nx - 1)
        else
          goto 100
        end if

        Np(1,1:999,1) = 0.0_dp
        Ni(1,1:999,1) = 0.0_dp
        Del(1,1:999,1) = DX

        do I = 2, Nx - 1, 1
          Del(I,1:999,1) = DX
          Np(I,1:999,1) = Np(I-1,1:999,1) + DX
          Ni(I,1:999,1) = (Np(I,1:999,1) + Np(I-1,1:999,1)) / 2.0_dp
        end do

        Np(Nx,1:999,1) = Lth(1) / 2.0_dp
        Ni(Nx,1:999,1) = (Np(Nx,1:999,1) + Np(Nx-1,1:999,1)) / 2.0_dp
        Ni(Nx+1,1:999,1) = Lth(1) / 2.0_dp

        do I = 1, Nx, 1
          Delta(I,1:999,1) = Ni(I+1,1:999,1) - Ni(I,1:999,1)
        end do

*       Mesh for Y direction
        ! Symmetry is being considered here. Only a quarter of domain is being modelled.
        if (Ny .GT. 1) then
          DY = (Lth(2) / 2.0_dp) / float(Ny - 1)
        else
          goto 100
        end if

        Np(1:999,1,2) = 0.0_dp
        Ni(1:999,1,2) = 0.0_dp
        Del(1:999,1,2) = DY

        do I = 2, Ny - 1, 1
          Del(1:999,I,2) = DY
          Np(1:999,I,2) = Np(1:999,I-1,2) + DY
          Ni(1:999,I,2) = (Np(1:999,I,2) + Np(1:999,I-1,2)) / 2.0_dp
        end do

        Np(1:999,Ny,2) = Lth(2) / 2.0_dp
        Ni(1:999,Ny,2) = (Np(1:999,Ny,2) + Np(1:999,Ny-1,2)) / 2.0_dp
        Ni(1:999,Ny+1,2) = Lth(2) / 2.0_dp

        do I = 1, Ny, 1
          Delta(1:999,I,2) = Ni(1:999,I+1,2) - Ni(1:999,I,2)
        end do

100   continue
      end subroutine Grid2DUni

      end module Geometry























*      subroutine Grid1DUni(Nx,Lx,NPx,NIx,DelX,DeltaX)
** Grid1DUni calculates the node positions in a one-dimensional setting
** Uni means mesh is uniformly spaced
** Interface is located at half distance from both neighbor nodes, except for
** boundary nodes, where node and interface coincides
*
** Nx: number of nodes in direction X
** Lx: length of direction X
** NPx: array containing node positions in X axis (to be calculated)
** NIx: array containing interface positions in X axis (to be calculated)
** DelX: array containing distance between neighboring nodes (to be calculated)
** DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)
*
** Nodes: index ranges from 0 to Nx-1. There are Nx nodes
** Interfaces: index ranges from 0 to Nx. There is one interface more than the number of nodes
** Distances between nodes (DelX): index ranges from 0 to Nx-2. One less than the number of nodes
** Distance between interfaces (DeltaX): index ranges from 0 to Nx-1. The same number of nodes (but one less than the number of interfaces)
*
*          implicit none
*
*          logical :: debugmode
*          common /dbgMode/ debugmode
*
*          integer :: Nx, I
*          real(dp) :: NPx(1:999),NIx(1:999), Lx, DX
*          real(dp) :: DelX(1:999), DeltaX(1:999)
*
**         Mesh for X direction
*          if (Nx .GT. 1) then
*            DX = Lx / float(Nx - 1)
*          else
*            goto 100
*          end if
*
*          NPx(1) = 0.0_dp
*          NIx(1) = 0.0_dp
*          DelX(1) = DX
*
**         Loop trough all nodes, except the last one at Nx index
*          do I = 2, Nx - 1, 1
*            DelX(I) = DX
*            NPx(I) = NPx(I-1) + DX
*            NIx(I) = (NPx(I) + NPx(I-1)) / 2.0_dp
*          end do
*
*          NPx(Nx) = Lx
*          NIx(Nx) = (NPx(Nx) + NPx(Nx-1)) / 2.0_dp
*          NIx(Nx+1) = Lx
*
*          do I = 1, Nx, 1
*            DeltaX(I) = NIx(I+1) - NIx(I)
*          end do
*
*          if (debugmode) then
*            call ShGeom1D(Nx,Lx,NPx,NIx,DelX,DeltaX)
*          end if
*
*100     continue
*
*      end subroutine Grid1DUni
*


*      subroutine Grid3DUni(N,Lth,NPx,NIx,NPy,NIy,NPz,NIz,Del,DeltaX,DeltaY,DeltaZ)
*
**     Original signature
**     subroutine Grid3DUni(N,Lth,NPx,NIx,NPy,NIy,NPz,NIz,DelX,DeltaX,DelY,DeltaY,DelZ,DeltaZ)
*
** Grid3DUni calculates the node positions in a tri-dimensional setting
** Uni means mesh is uniformly spaced
*
**     N: number of nodes in each direction
**        1,2,3 indices refer to x,y,z respectively
*
**     Lth: Lenght along each direction
**          1,2,3 indices refer to x,y,z respectively
** NPx: array containing node positions in X axis (to be calculated)
** NIx: array containing interface positions in X axis (to be calculated)
*
** NPy: array containing node positions in Y axis (to be calculated)
** NIy: array containing interface positions in Y axis (to be calculated)
*
** NPz: array containing node positions in Z axis (to be calculated)
** NIz: array containing interface positions in Z axis (to be calculated)
*
** DelX: array containing distance between neighboring nodes in X axis(to be calculated)
** DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)
*
** DelY: array containing distance between neighboring nodes in Y axis(to be calculated)
** DeltaY: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)
*
** DelZ: array containing distance between neighboring nodes in Y axis(to be calculated)
** DeltaZ: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)
*
** Nodes: index ranges from 0 to N-1. There are Nx nodes
** Interfaces: index ranges from 0 to N. There is one interface more than the number of nodes
** Distances between nodes (DelX,DelY,DelZ): index ranges from 0 to N-2. One less than the number of nodes
** Distance between interfaces (DeltaY,DeltaY,DeltaZ): index ranges from 0 to N-1. The same number of nodes (but one less than the number of interfaces)
*
*        logical :: debugmode
*        common /dbgMode/ debugmode
*
*        integer :: N(1:3), I, J, K
*        real(dp) :: Lth(1:3)
*        real(dp) :: NPx(0:999), NIx(0:999),DX
*        real(dp) :: DelX(0:999), DeltaX(0:999)
*        real(dp) :: Del(0:999,0:999,0:999,1:3)
*
*        real(dp) :: NPy(0:999), NIy(0:999),DY
*        real(dp) :: DelY(0:999), DeltaY(0:999)
*
*        real(dp) :: NPz(0:999), NIz(0:999),DZ
*        real(dp) :: DelZ(0:999), DeltaZ(0:999)
*
*        real(dp) :: X,Y,Z
*
**       Mesh for X direction
*        if (N(1) .GT. 1) then
*          DX = Lth(1) / float(N(1) - 1)
*        else
*          goto 100
*        end if
*
*        NPx(0) = 0
*        NIx(0) = 0
*
**TODO To be deleted soon
*        DelX(0) = DX
*
*        Del(0,0:999,0:999,1) = DX
*
*        do I = 1, N(1) - 2, 1
*
**TODO To be deleted soon
*          DelX(I) = DX
*
*
*          Del(I,0:999,0:999,1) = DX
*          NPx(I) = NPx(I-1) + DX
*          NIx(I) = (NPx(I) + NPx(I-1)) / 2.0_dp
*        end do
*
*        NPx(N(1)-1) = Lth(1)
*        NIx(N(1)-1) = (NPx(N(1)-1) + NPx(N(1)-2)) / 2.0_dp
*        NIx(N(1)) = Lth(1)
*
*        do I = 0, N(1) - 1, 1
*          DeltaX(I) = NIx(I+1) - NIx(I)
*        end do
*
**       Mesh for Y direction
*        if (N(2) .GT. 1) then
*          DY = Lth(2) / float(N(2) - 1)
*        else
*          goto 100
*        end if
*
*        NPy(0) = 0
*        NIy(0) = 0
*
**TODO To be deleted soon
*        DelY(0) = DY
*
*
*        Del(0:999,0,0:999,2) = DY
*
*        do I = 1, N(2) - 2, 1
*
**TODO To be deleted soon
*          DelY(I) = DY
*
*
*
*          Del(0:999,I,0:999,2) = DY
*          NPy(I) = NPy(I-1) + DY
*          NIy(I) = (NPy(I) + NPy(I-1)) / 2.0_dp
*        end do
*
*        NPy(N(2)-1) = Lth(2)
*        NIy(N(2)-1) = (NPy(N(2)-1) + NPy(N(2)-2)) / 2.0_dp
*        NIy(N(2)) = Lth(2)
*
*        do I = 0, N(2) - 1, 1
*          DeltaY(I) = NIy(I+1) - NIy(I)
*        end do
*
**       Mesh for Z direction
*        if (N(3) .GT. 1) then
*          DZ = Lth(3) / float(N(3) - 1)
*        else
*          goto 100
*        end if
*
*        NPz(0) = 0
*        NIz(0) = 0
*
**TODO To be deleted soon
*        DelZ(0) = DZ
*
*
*
*        Del(0:999,0:999,0,3) = DZ
*
*        do I = 1, N(3) - 2, 1
*
**TODO To be deleted soon
*          DelZ(I) = DZ
*
*
*
*
*          Del(0:999,0:999,I,3) = DZ
*          NPz(I) = NPz(I-1) + DZ
*          NIz(I) = (NPz(I) + NPz(I-1)) / 2.0_dp
*        end do
*
*        NPz(N(2)-1) = Lth(2)
*        NIz(N(2)-1) = (NPz(N(2)-1) + NPz(N(2)-2)) / 2.0_dp
*        NIz(N(2)) = Lth(2)
*
*        do I = 0, N(2) - 1, 1
*          DeltaZ(I) = NIz(I+1) - NIz(I)
*        end do
*
*500     FORMAT(' ', A,I3,A,I3,A,I3,A,F5.3,A,F5.3,A,F5.3,A)
*        if (debugmode) then
*          do I = 0, N(1)-1,1
*            X = NPx(I)
*            do J = 0, N(2)-1,1
*              Y = NPy(J)
*              do K=0,N(3)-1,1
*                Z = NPz(K)
*      write(*,500) "N(",I,",",J,",",K,")=(",X,",",Y,",",Z,")"
*              end do
*              write(*,*)
*            end do
*            write(*,*)
*          end do
*        end if
*
*100   continue
*
*      end subroutine Grid3DUni

* RefinedGrid





* Output
*      subroutine ShGeom1D(Nx,Lx,NPx,NIx,DelX,DeltaX)
*        implicit none
*
*        logical :: debugmode
*        common /dbgMode/ debugmode
*
*        integer :: Nx, I
*        real(dp) :: NPx(1:999),NIx(1:999), Lx
*        real(dp) :: DelX(1:999), DeltaX(1:999)
*
*100     FORMAT(' ', A,I1,A,F6.3)
*200     FORMAT(' ', A,I2,A,F7.3)
*300     FORMAT(' ', A,I3,A,F8.3)
*        if (debugmode) then
*          do I = 1, Nx,1
*            if (I .LT. 10) then
*              write(*,100) "NPx(", I,")=", NPx(I)
*              write(*,100) "NIx(", I,")=", NIx(I)
*              write(*,100) "DelX(",I,")=", DelX(I)
*              write(*,100) "DeltaX(",I,")=", DeltaX(I)
*            end if
*            if ((I .GE. 10) .AND. (I .LT. 100)) then
*              write(*,200) "NPx(", I,")=", NPx(I)
*              write(*,200) "NIx(", I,")=", NIx(I)
*              write(*,200) "DelX(",I,")=", DelX(I)
*              write(*,200) "DeltaX(",I,")=", DeltaX(I)
*            end if
*            if (I .GE. 100) then
*              write(*,300) "NPx(", I,")=", NPx(I)
*              write(*,300) "NIx(", I,")=", NIx(I)
*              write(*,300) "DelX(",I,")=", DelX(I)
*              write(*,300) "DeltaX(",I,")=", DeltaX(I)
*            end if
*
*            write(*,*) ""
*          end do
*          write(*,300) "NIx(", Nx,")=", NIx(Nx)
*
*        end if
*
*      end subroutine ShGeom1D
