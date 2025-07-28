*******************************************************
* FILE: Geometry.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Geometry
        implicit none

      contains

      subroutine Grid2DUni(Lth,Np,Ni)

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
        double precision :: Lth(1:2)
        double precision :: Np(1:999,1:999,1:2), NI(1:999,1:999,1:2), DX, DY

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        Nx = N(1)
        Ny = N(2)

*       Mesh for X direction
        ! Symmetry is being considered here. Only a quarter of domain is being modelled.
        if (Nx .GT. 1) then
          DX = (Lth(1) / 2.0D0) / float(Nx - 1)
        else
          goto 100
        end if

        Np(1,1:999,1) = 0.0D0
        Ni(1,1:999,1) = 0.0D0
        Del(1,1:999,1) = DX

        do I = 2, Nx - 1, 1
          Del(I,1:999,1) = DX
          Np(I,1:999,1) = Np(I-1,1:999,1) + DX
          Ni(I,1:999,1) = (Np(I,1:999,1) + Np(I-1,1:999,1)) / 2.0D0
        end do

        Np(Nx,1:999,1) = Lth(1) / 2.0D0
        Ni(Nx,1:999,1) = (Np(Nx,1:999,1) + Np(Nx-1,1:999,1)) / 2.0D0
        Ni(Nx+1,1:999,1) = Lth(1) / 2.0D0

        do I = 1, Nx, 1
          Delta(I,1:999,1) = Ni(I+1,1:999,1) - Ni(I,1:999,1)
        end do

*       Mesh for Y direction
        ! Symmetry is being considered here. Only a quarter of domain is being modelled.
        if (Ny .GT. 1) then
          DY = (Lth(2) / 2.0D0) / float(Ny - 1)
        else
          goto 100
        end if

        Np(1:999,1,2) = 0.0D0
        Ni(1:999,1,2) = 0.0D0
        Del(1:999,1,2) = DY

        do I = 2, Ny - 1, 1
          Del(1:999,I,2) = DY
          Np(1:999,I,2) = Np(1:999,I-1,2) + DY
          Ni(1:999,I,2) = (Np(1:999,I,2) + Np(1:999,I-1,2)) / 2.0D0
        end do

        Np(1:999,Ny,2) = Lth(2) / 2.0D0
        Ni(1:999,Ny,2) = (Np(1:999,Ny,2) + Np(1:999,Ny-1,2)) / 2.0D0
        Ni(1:999,Ny+1,2) = Lth(2) / 2.0D0

        do I = 1, Ny, 1
          Delta(1:999,I,2) = Ni(1:999,I+1,2) - Ni(1:999,I,2)
        end do

100   continue
      end subroutine Grid2DUni

      end module Geometry
