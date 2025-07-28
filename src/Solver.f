*******************************************************
* FILE: Solver.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************

      module Solver
        use Coefficients
        use Properties

        implicit none

      contains

*     Jacobi algorithm
      subroutine Jacobi(Delta_t,tol,alpha)
        implicit none

        integer :: N(1:2)
        common /gridsize/ N

        double precision :: maxDeltaT

        double precision :: tol, total, alpha, resid

        integer :: I, J, Nx, Ny

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        double precision :: Delta_t

        double precision :: Np(1:999,1:999,1:2),Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        double precision :: TitV(0:999,0:999)

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)

        character(LEN=50) :: name
        logical :: last
        name = "Jacobi"

        Nx = N(1)
        Ny = N(2)

        Told(0,:) = Told(2,:)
        Told(:,0) = Told(:,2)
        Told(Nx+1,:) = 0.0D0
        Told(:,Ny+1) = 0.0D0

        ToldMax = Tmax
        ToldMin = Tmin
        ToldAvg = Tavg

*       Setting temperature to initial temperatures, stored in Told
        TitV(:,:) = Told(:,:)
        Tnew(:,:) = Told(:,:)

*       Implementation of Jacobi Algorithm
        resid = 100.0D0
        iter = 0
        do while (resid .GT. tol)

*         Calculate coefficients for a 2D problem
          call coeff2D(Aw,Ae,An,As,Ap,b,Delta_t)

          TitV(0,:) = TitV(2,:)
          TitV(:,0) = TitV(:,2)

*         Solving linear system
          do I = 1, Nx , 1
            do J = 1, Ny , 1
              total = Aw(I,J) * TitV(I-1,J)
              total = total + Ae(I,J) * TitV(I+1,J)
              total = total + As(I,J) * TitV(I,J-1)
              total = total + An(I,J) * TitV(I,J+1)
              total = total + b(I,J)

              Tnew(I,J) = TitV(I,J) + alpha * (total / Ap(I,J) - TitV(I,J))

              maxDeltaT = abs(TitV(I,J) - Tnew(I,J))

              if((abs(TitV(I,J) - Tnew(I,J))) .GT. maxDeltaT) then
                maxDeltaT = abs(TitV(I,J) - Tnew(I,J))
              end if

            end do
          end do

*         Monitoring convergence
          resid = calcResid(Aw,Ae,An,As,Ap,b)
          TitV(:,:) = Tnew(:,:)
          iter = iter + 1
          totalIter = totalIter + 1
          last = .FALSE.
          call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)

        end do

        last = .TRUE.
        call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)
        Told(:,:) = Tnew(:,:)

      end subroutine Jacobi




      subroutine GaSe2D(Delta_t,tol,alpha)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: maxDeltaT

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        double precision :: TitV(0:999,0:999)
        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        double precision :: Delta_t

        double precision :: Np(1:999,1:999,1:2),Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx, Ny

        double precision :: tol, temp, total,alpha,resid

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)

        logical :: last
        character(LEN=50) :: name
        name = "Gauss-Seidel"

*       Boundary conditions
        Nx = N(1)
        Ny = N(2)

        Told(0,:) = Told(2,:)
        Told(:,0) = Told(:,2)
        Told(Nx+1,:) = 0.0D0
        Told(:,Ny+1) = 0.0D0

        ToldMax = Tmax
        ToldMin = Tmin
        ToldAvg = Tavg

*       Setting temperature to initial temperatures, stored in Told
        TitV(:,:) = Told(:,:)
        Tnew(:,:) = Told(:,:)

*       Implementation of Gauss-Seidel Algorithm (properly)
        resid = 100.0D0
        iter = 0
        do while (resid .GT. tol)

*         Calculate coefficients for a 2D problem
          call coeff2D(Aw,Ae,An,As,Ap,b,Delta_t)

          maxDeltaT = 0.0D0

*         Solving linear system
          do I = 1, Nx , 1
            do J = 1, Ny, 1
              total = Aw(I,J) * Tnew(I-1,J)
              total = total + Ae(I,J) * Tnew(I+1,J)
              total = total + As(I,J) * Tnew(I,J-1)
              total = total + An(I,J) * Tnew(I,J+1)
              total = total + b(I,J)
              temp = Tnew(I,J)

              if (alpha .NE. 1.0D0) then
                Tnew(I,J) = temp + alpha * (total / Ap(I,J) - temp)
              else
                Tnew(I,J) = total / Ap(I,J)
              end if

              if((abs(temp - Tnew(I,J))) .GT. maxDeltaT) then
                maxDeltaT = abs(temp - Tnew(I,J))
              end if

            end do
          end do

*         Monitoring convergence
          resid = calcResid(Aw,Ae,An,As,Ap,b)
          TitV(:,:) = Tnew(:,:)
          iter = iter + 1
          totalIter = totalIter + 1
          last = .FALSE.
          call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)
        end do

        last = .TRUE.
        call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)
        Told(:,:) = Tnew(:,:)

      end subroutine GaSe2D




*     TDMA
      subroutine TDMA2D(Delta_t,tol,alpha)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx, Ny

        double precision :: resid, tol, temp, total,alpha,Delta_t

        double precision :: TitV(0:999,0:999)
        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        double precision :: maxDeltaT

        double precision :: Tw, Te, Ts, Tn

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)

*       Solve each horizontal lines (nodes along X axis are solved at each line)
*       Temperatures above and below the line being calculated are considered to be known
        double precision :: P(0:999),Q(0:999)

        double precision :: Np(1:999,1:999,1:2),Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        logical :: last
        character(LEN=50) :: name
        name = "TDMA 2D (without ADI)"

        maxDeltaT = 0.0D0

*       Boundary conditions
        Nx = N(1)
        Ny = N(2)

        Told(0,:) = 0.0D0
        Told(:,0) = 0.0D0
        Told(Nx+1,:) = 0.0D0
        Told(:,Ny+1) = 0.0D0

        ToldMax = Tmax
        ToldMin = Tmin
        ToldAvg = Tavg

*       Setting temperature to initial temperatures, stored in Told
        TitV(:,:) = Told(:,:)
        Tnew(:,:) = Told(:,:)

        iter = 0
        temp = 1
        resid= 1.0D5
        do while (resid .GT. tol)

*         Calculate coefficients for a 2D problem
          call coeff2D(Aw,Ae,An,As,Ap,b,Delta_t)

          do J = 1, Ny, 1

*           define P1, Q1
            P(1) = Ae(1,J) / Ap(1,J)

            Tw = Tnew(0,J)
            Tn = Tnew(1,J+1)
            Ts = Tnew(1,J-1)

            total = Aw(1,J) * Tw
            total = total + An(1,J) * Tn
            total = total + As(1,J) * Ts
            total = total + b(1,J)

            Q(1) = total / Ap(1,J)

*           define Pi, Qi for i = 2 up to Nx-1
            do I = 2, Nx-1, 1

              P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))

              Tn = Tnew(I,J+1)
              Ts = Tnew(I,J-1)

              total = An(I,J) * Tn
              total = total + As(I,J) * Ts
              total = total + b(I,J)
              total = total + Aw(I,J) * Q(I-1)
              Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))

            end do

*           define Pi, Qi for i = Nx, the last unknown
            P(Nx) = 0

            Tn = Tnew(Nx,J+1)
            Ts = Tnew(Nx,J-1)
            Te = Tnew(Nx+1,J)

            total = An(Nx,J) * Tn
            total = total + As(Nx,J) * Ts
            total = total + Ae(Nx,J) * Te
            total = total + b(Nx,J)
            total = total + Aw(Nx,J) * Q(Nx-1)

            Q(Nx) = total / (Ap(Nx,J) - Aw(Nx,J) * P(Nx-1))

*           Calculate Tnew for the last unknown node
            temp = Tnew(Nx,J)
            Tnew(Nx,J) = temp + alpha *(Q(Nx) - temp)
            maxDeltaT = abs(temp - Tnew(Nx,J))


            do I = Nx - 1, 1, -1
              temp = Tnew(I,J)
              Tnew(I,J) = temp + alpha * (P(I) * Tnew(I+1,J) + Q(I) - temp)

              if((abs(temp - Tnew(I,J))) .GT. maxDeltaT) then
                maxDeltaT = abs(temp - Tnew(I,J))
              end if
            end do
          end do


*         Monitoring convergence
          resid = calcResid(Aw,Ae,An,As,Ap,b)
          TitV(:,:) = Tnew(:,:)
          iter = iter + 1
          totalIter = totalIter + 1
          last = .FALSE.
          call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)

        end do

        last = .TRUE.
        call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)
        Told(:,:) = Tnew(:,:)

      end subroutine TDMA2D





      subroutine TDMA2dADI(Delta_t,tol,alpha)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: maxDeltaT

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        double precision :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx, Ny

        double precision :: resid, tol, temp, total,alpha,Delta_t

        double precision :: TitV(0:999,0:999)
        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        double precision :: Tw, Te, Ts, Tn

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)

*       Solve each horizontal lines (nodes along X axis are solved at each line)
*       Temperatures above and below the line being calculated are considered to be known
        double precision :: P(0:999),Q(0:999)

        double precision :: Np(1:999,1:999,1:2),Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        logical :: last
        character(LEN=50) :: name
        name = "TDMA 2D with ADI"

*       Boundary conditions
        Nx = N(1)
        Ny = N(2)

        Told(0,:) = 0.0D0
        Told(:,0) = 0.0D0
        Told(Nx+1,:) = 0.0D0
        Told(:,Ny+1) = 0.0D0

        ToldMax = Tmax
        ToldMin = Tmin
        ToldAvg = Tavg

*       Setting temperature to initial temperatures, stored in Told
        TitV(:,:) = Told(:,:)
        Tnew(:,:) = Told(:,:)

        iter = 0
        temp = 1
        resid= 100.0D0
        do while (resid .GT. tol)

*         Calculate coefficients for a 2D problem
          call coeff2D(Aw,Ae,An,As,Ap,b,Delta_t)

*         Scan the domain from low Y (lower J) to high Y (higher J)
          if (mod(iter,4) .EQ. 0) then

            do J = 1, Ny, 1

*             define P1, Q1
              P(1) = Ae(1,J) / Ap(1,J)

              Tw = Tnew(0,J)
              Tn = Tnew(1,J+1)
              Ts = Tnew(1,J-1)

              total = Aw(1,J) * Tw
              total = total + An(1,J) * Tn
              total = total + As(1,J) * Ts
              total = total + b(1,J)

              Q(1) = total / Ap(1,J)

*             define Pi, Qi for I = 2 up to Nx - 1
              do I = 2, Nx - 1, 1
                P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))

                Tn = Tnew(I,J+1)
                Ts = Tnew(I,J-1)

                total = An(I,J) * Tn
                total = total + As(I,J) * Ts
                total = total + b(I,J)
                total = total + Aw(I,J) * Q(I-1)

                Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))
              end do

*             define Pi, Qi for I = Nx, the last unknown
              P(Nx) = 0

              Tn = Tnew(Nx,J+1)
              Ts = Tnew(Nx,J-1)
              Te = Tnew(Nx+1,J)

              total = An(Nx,J) * Tn
              total = total + As(Nx,J) * Ts
              total = total + Ae(Nx,J) * Te
              total = total + b(Nx,J)
              total = total + Aw(Nx,J) * Q(Nx-1)

              Q(Nx) = total / (Ap(Nx,J) - Aw(Nx,J) * P(Nx-1))

*             Calculate Tnew for the last unknown node
              temp = Tnew(Nx,J)
              Tnew(Nx,J) = temp + alpha * (Q(Nx) - temp)

              do I = Nx - 1, 1, -1
                temp = Tnew(I,J)
                Tnew(I,J) = temp + alpha* (P(I) * Tnew(I+1,J) + Q(I) - temp)
              end do
            end do

          end if

*         Scan the domain from high Y (higher J) to low Y (lower J)
          if (mod(iter,4) .EQ. 1) then

            do J = Ny - 2, 1, -1

*             define P1, Q1
              P(1) = Ae(1,J) / Ap(1,J)

              Tw = Tnew(0,J)
              Tn = Tnew(1,J+1)
              Ts = Tnew(1,J-1)

              total = Aw(1,J) * Tw
              total = total + An(1,J) * Tn
              total = total + As(1,J) * Ts
              total = total + b(1,J)

              Q(1) = total / Ap(1,J)

*             define Pi, Qi for I = 2 up to Nx - 1
              do I = 2, Nx - 1, 1
                P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))

                Tn = Tnew(I,J+1)
                Ts = Tnew(I,J-1)

                total = An(I,J) * Tn
                total = total + As(I,J) * Ts
                total = total + b(I,J)
                total = total + Aw(I,J) * Q(I-1)

                Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))
              end do

*             define Pi, Qi for I = Nx, the last unknown
              P(Nx) = 0

              Tn = Tnew(Nx,J+1)
              Ts = Tnew(Nx,J-1)
              Te = Tnew(Nx+1,J)

              total = An(Nx,J) * Tn
              total = total + As(Nx,J) * Ts
              total = total + Ae(Nx,J) * Te
              total = total + b(Nx,J)
              total = total + Aw(Nx,J) * Q(Nx-1)

              Q(Nx) = total / (Ap(Nx,J) - Aw(Nx,J) * P(Nx-1))

*             Calculate Tnew for the last unknown node
              temp = Tnew(Nx,J)
              Tnew(Nx,J) = temp + alpha * (Q(Nx) - temp)

              do I = Nx - 1, 1, -1
                temp = Tnew(I,J)
                Tnew(I,J) = temp + alpha* (P(I) * Tnew(I+1,J) + Q(I) - temp)
              end do
            end do

          end if

*         Scan the domain from low X (lower I) to high X (higher I)
          if (mod(iter,4) .EQ. 2) then

            do I = 1, Nx, 1

*             define P1, Q1
              P(1) = An(I,1) / Ap(I,1)

                Ts = Tnew(I,0)
                Te = Tnew(I+1,1)
                Tw = Tnew(I-1,1)

                total = As(I,1) * Ts
                total = total + Ae(I,1) * Te
                total = total + Aw(I,1) * Tw
                total = total + b(I,1)

              Q(1) = total / Ap(I,1)

*             define Pi, Qi for J = 2 up to Ny - 1
              do J = 2, Ny - 1, 1
                P(J) = An(I,J) / (Ap(I,J) - As(I,J) * P(J-1))

                  Te = Tnew(I+1,J)
                  Tw = Tnew(I-1,J)

                  total = Ae(I,J) * Te
                  total = total + Aw(I,J) * Tw
                  total = total + b(I,J)
                  total = total + As(I,J) * Q(J-1)

                Q(J) = total / (Ap(I,J) - As(I,J) * P(J-1))
              end do

*             define Pi, Qi for J = Ny, the last unknown
              P(Ny) = 0

              Tn = Tnew(I,Ny+1)
              Tw = Tnew(I-1,Ny)
              Te = Tnew(I+1,Ny)

              total = An(I,Ny) * Tn
              total = total + Aw(I,Ny) * Tw
              total = total + Ae(I,Ny) * Te
              total = total + b(I,Ny)
              total = total + As(I,Ny) * Q(Ny-1)

              Q(Ny) = total / (Ap(I,Ny) - As(I,Ny) * P(Ny-1))

*             Calculate Tnew for the last unknown node
              temp = Tnew(I,Ny)
              Tnew(I,Ny) = temp + alpha * (Q(Ny) - temp)

              do J = Ny - 1, 1, -1
                temp = Tnew(I,J)
                Tnew(I,J) = temp + alpha*(P(J) * Tnew(I,J+1) + Q(J) - temp)
              end do

            end do

          end if

*         Scan the domain from high X (higher I) to low X (lower I)
          if (mod(iter,4) .EQ. 3) then

            do I = Nx - 2, 1, -1

*             define P1, Q1
              P(1) = An(I,1) / Ap(I,1)

                Ts = Tnew(I,0)
                Te = Tnew(I+1,1)
                Tw = Tnew(I-1,1)

                total = As(I,1) * Ts
                total = total + Ae(I,1) * Te
                total = total + Aw(I,1) * Tw
                total = total + b(I,1)

              Q(1) = total / Ap(I,1)

*             define Pi, Qi for J = 2 up to Ny - 1
              do J = 2, Ny - 1, 1
                P(J) = An(I,J) / (Ap(I,J) - As(I,J) * P(J-1))

                  Te = Tnew(I+1,J)
                  Tw = Tnew(I-1,J)

                  total = Ae(I,J) * Te
                  total = total + Aw(I,J) * Tw
                  total = total + b(I,J)
                  total = total + As(I,J) * Q(J-1)

                Q(J) = total / (Ap(I,J) - As(I,J) * P(J-1))
              end do

*             define Pi, Qi for J = Ny, the last unknown
              P(Ny) = 0

              Tn = Tnew(I,Ny+1)
              Tw = Tnew(I-1,Ny)
              Te = Tnew(I+1,Ny)

              total = An(I,Ny) * Tn
              total = total + Aw(I,Ny) * Tw
              total = total + Ae(I,Ny) * Te
              total = total + b(I,Ny)
              total = total + As(I,Ny) * Q(Ny-1)

              Q(Ny) = total / (Ap(I,Ny) - As(I,Ny) * P(Ny-1))

*             Calculate Tnew for the last unknown node
              temp = Tnew(I,Ny)
              Tnew(I,Ny) = temp + alpha * (Q(Ny) - temp)
              maxDeltaT = abs(temp - Tnew(I,J))


              do J = Ny - 1, 1, -1
                temp = Tnew(I,J)
                Tnew(I,J) = temp + alpha*(P(J) * Tnew(I,J+1) + Q(J) - temp)

                if((abs(temp - Tnew(I,J))) .GT. maxDeltaT) then
                  maxDeltaT = abs(temp - Tnew(I,J))
                end if
              end do

            end do

          end if

*         Monitoring convergence
          resid = calcResid(Aw,Ae,An,As,Ap,b)
          TitV(:,:) = Tnew(:,:)
          iter = iter + 1
          totalIter = totalIter + 1
          last = .FALSE.
          call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)

        end do

        last = .TRUE.
        call convMonitor(name,resid,tol,maxDeltaT,last,Delta_t)
        Told(:,:) = Tnew(:,:)

      end subroutine TDMA2dADI




      function calcResid(Aw,Ae,An,As,Ap,b)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        double precision :: SteelChem(1:6)
        common /steelprop/ SteelChem

        double precision :: Tliq

        double precision :: calcResid
        double precision :: resid

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx, Ny

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        integer :: Idummy, Jdummy


        double precision :: sumT

        double precision :: Aw(1:999,1:999), Ae(1:999,1:999)
        double precision :: An(1:999,1:999), As(1:999,1:999)
        double precision :: Ap(1:999,1:999), b(1:999,1:999)
        double precision :: total

        Nx = N(1)
        Ny = N(2)

        Tliq = Tliquidus()

        Idummy = Imin
        Jdummy = Jmin

        sumT = 0.0D0
        resid = 0.0D0
        Tmax = 0.0D0
        Tmin = 1.0D5
        do I = 1, Nx , 1
          do J = 1, Ny, 1
            total = Aw(I,J) * Tnew(I-1,J)
            total = total + Ae(I,J) * Tnew(I+1,J)
            total = total + As(I,J) * Tnew(I,J-1)
            total = total + An(I,J) * Tnew(I,J+1)
            total = total + b(I,J)

            sumT = sumT + Tnew(I,J)

            if ((Tnew(I,J) - Tmax) .GT. 1.0D-8) then
              Tmax = Tnew(I,J)
            end if

            if ((Tmin - Tnew(I,J)) .GT. 1.0D-8) then
              Tmin = Tnew(I,J)
              Imin = I
              Jmin = J
            end if

            resid = resid + (Ap(I,J) * Tnew(I,J) - total)*(Ap(I,J) * Tnew(I,J) - total)
          end do
        end do

        resid = sqrt(resid) / (Nx * Ny)
        calcResid = resid
        Tavg = sumT / (Nx * Ny)

      end function calcResid





      subroutine convMonitor(name,resid,tol,maxDeltaT,last,tstep)

        implicit none

        character(LEN=50) :: name
        double precision :: resid, tol,tstep
        logical :: last

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 50000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Equivalent Specific heat
*          6: Thermal conductivity
*          7: Liquid fraction

*       ShThick
*       First index: 1 for X axis, 2 for Y axis
*       Second index: curTime
        double precision :: TmRs(1:81,0:50000,1:7), ShThick(1:3,0:50000)
        common /transient/ TmRs, ShThick

        integer :: maxX, maxY, Xi(9), Yi(9), tsSavPos
        common /storeLimits/ maxX, maxY, Xi, Yi, tsSavPos

        double precision :: maxDeltaT

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        double precision :: TmaxCh, TminCh, TavgCh

        double precision :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        integer :: tsShLog
        common /convergence/ tsShLog

130     format(' ', A,I12)
230     format(' ', A,F7.2, A,F7.2, A)
330     format(' ', A,E12.6)
430     format(' ', A,F16.6, A)

        TmaxCh = Tmax - ToldMax
        TminCh = Tmin - ToldMin
        TavgCh = Tavg - ToldAvg

*       Monitoring convergence
        if ((mod(iter,1000) .EQ. 0) .AND. (last .EQV. .FALSE.)) then
          print*
          print*,"=================================================================================================="
          print*,"Intermediary Step"
          print*,"--------------------------------------------------------------------------------------------------"
          print 430,"CURRENT TIME STEP:         ", tstep, " seconds."
          print*
          print 430,"CURRENT TIME:              ", curTime, " seconds."
          print*
          print*,"Solver method: ", trim(name)
          print 130,"CONVERGENCE Monitor: Iterations = ", iter
          print 130,"Total iterations:                 ", totalIter
          print*,""
          print 330,"Root-Mean Square Residue:            ", resid
          print 330,"Tolerance:                           ", tol
          print 330,"Max difference between iterations:   ", maxDeltaT
          print*,""
          print 230,"Maximum Temperature:                           ", Tmax, " K - Change from last time step: ", TmaxCh, " K"
          print 230,"Minimum Temperature:                           ", Tmin, " K - Change from last time step: ", TminCh, " K"
          print 230,"Average(algebraic, non-weigthed) Temperature : ", Tavg, " K - Change from last time step: ", TavgCh, " K"
          print*,""
          print 430,"Liquidus Temperature: ", Tliquidus(), " K"
          print 430,"Solidus Temperature:  ", Tsolidus(), " K"
          print*,"=================================================================================================="
          print*," "
          print*," "
        end if

        if ((mod(iter,1000) .GT. 0) .AND. (last .EQV. .TRUE.)) then

            if ((totalIter - tsShLog) .GT. 1000) then
              print*
              print*,"=================================================================================================="
              print*,"Final Step"
              print*,"--------------------------------------------------------------------------------------------------"
              print 430,"CURRENT TIME STEP:         ", tstep, " seconds."
              print*
              print 430,"CURRENT TIME:              ", curTime, " seconds."
              print*
              print*,"Solver method: ", trim(name)
              print 130,"CONVERGENCE Monitor: Iterations = ", iter
              print 130,"Total iterations:                 ", totalIter
              print*,""
              print 330,"Root-Mean Square Residue:            ", resid
              print 330,"Tolerance:                           ", tol
              print 330,"Max difference between iterations:   ", maxDeltaT
              print*,""
              print 230,"Maximum Temperature:                           ", Tmax, " K - Change from last time step: ", TmaxCh, " K"
              print 230,"Minimum Temperature:                           ", Tmin, " K - Change from last time step: ", TminCh, " K"
              print 230,"Average(algebraic, non-weigthed) Temperature : ", Tavg, " K - Change from last time step: ", TavgCh, " K"
              print*,""
              print 430,"Liquidus Temperature: ", Tliquidus(), " K"
              print 430,"Solidus Temperature:  ", Tsolidus(), " K"
              print*,"=================================================================================================="
              print*," "
              print*," "

              tsShLog = totalIter
            end if
        end if
      end subroutine convMonitor



      subroutine calcOldStats()

        implicit none

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        integer :: Idummy, Jdummy

        double precision :: maxT,minT,avgT,sumT

        maxT = 0.0D0
        minT = 1.0D5
        sumT = 0.0D0

        Idummy = Imin
        Jdummy = Jmin

        do I = 1, N(1), 1

          do J = 1, N(2), 1

            sumT = sumT + Told(I,J)

            if (Told(I,J) .GT. maxT) then
              maxT = Told(I,J)
            end if

            if(Told(I,J) .LT. minT) then
              minT = Told(I,J)
              Imin = I
              Jmin = J
            end if

          end do
        end do

        avgT = sumT / (N(1) * N(2))

        ToldMax = maxT
        ToldMin = minT
        ToldAvg = avgT

      end subroutine calcOldStats


      function ShellThick(Ip,axis)
        implicit none

        double precision :: ShellThick
        ! axis = 1 > X axis
        ! axis = 2 > Y axis
        integer ::  Ip, axis

        integer :: I, J
        integer :: Nx, Ny

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        double precision :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempres/ Tnew, Told

        double precision :: Np(1:999,1:999,1:2), NI(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        ! XposLo, YposLo: marks the index of position with liquid fraction just below 0.01
        ! XposHi, XposHi: marks the index of position with liquid fraction just above 0.01
        integer :: XposLo, YposLo
        integer :: XposHi, YposHi

        double precision :: dummy, m

        double precision :: near_f_Lo, near_f_Hi

        double precision :: liqFrac
        logical :: skin

        Nx = N(1)
        Ny = N(2)

        near_f_Lo = -1.0D0
        near_f_Hi = 1.0D0
        XposHi = 1
        YposHi = 1
        XposLo = 1
        YposLo = 1

        skin = .FALSE.

        if (axis .EQ. 1) then

          ! Scan Y axis for a X constant
          do J = Ny, 1, -1
            liqFrac = f_Liq(Tnew(Ip,J))

            if (liqFrac .GE. near_f_Lo) then
              if ((liqFrac - 0.05D0) .LE. 0.0D0) then
                near_f_Lo = liqFrac

                XposLo = Ip
                YposLo = J
                skin = .TRUE.
              end if
            end if

            if (liqFrac .LE. near_f_Hi) then
              if ((liqFrac - 0.05D0) .GT. 0.0D0) then

                near_f_Hi = liqFrac
                XposHi = Ip
                YposHi = J

              end if
            end if

          end do

        end if

        if (axis .EQ. 2) then

          ! Scan X axis for a Y constant
          do I = Nx, 1, -1
            liqFrac = f_Liq(Tnew(I,Ip))

            if (liqFrac .GE. near_f_Lo) then
              if ((liqFrac - 0.05D0) .LE. 0.0D0) then
                near_f_Lo = liqFrac
                XposLo = I
                YposLo = Ip
                skin = .TRUE.
              end if
            end if

            if (liqFrac .LE. near_f_Hi) then
              if ((liqFrac - 0.05D0) .GT. 0.0D0) then
                near_f_Hi = liqFrac
                XposHi = I
                YposHi = Ip
              end if
            end if

          end do

        end if

        ShellThick = 0.0D0

        if (axis .EQ. 1) then
          if (skin .EQV. .TRUE.) then
            m = (near_f_Lo - near_f_Hi)
            m = m / (Np(XposLo,YposLo,2) - Np(XposHi,YposHi,2))
            dummy = (0.05D0 - near_f_Lo) / m
            dummy = Np(XposLo,YposLo,2) + dummy
            ShellThick = Np(Nx,Ny,2) - dummy
          end if
        end if

        if (axis .EQ. 2) then
          if (skin .EQV. .TRUE.) then
            m = (near_f_Lo - near_f_Hi)
            m = m / (Np(XposLo,YposLo,1) - Np(XposHi,YposHi,1))
            dummy = (0.05D0 - near_f_Lo) / m
            dummy = Np(XposLo,YposLo,1) + dummy
            ShellThick = Np(Nx,Ny,1) - dummy
          end if
        end if

      end function ShellThick

      end module Solver
