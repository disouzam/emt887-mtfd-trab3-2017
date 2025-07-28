*******************************************************
* FILE: Results.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Results
        use Properties

        implicit none

      contains

      subroutine print_res1D(Nx,Temp)
        use CustomDouble

        implicit none

        integer :: Nx, I
        real(dp) :: Temp(0:200,0:0)

        print*, "========================================================"
        print*, "Results for 1D calculation will be printed out on screen"
        print*, "========================================================"
        print*,""
        print*,""

        do I = 1, Nx, 1
          print*,"T(",I,",0) = ", Temp(I,0)
        end do

      end subroutine print_res1D




      subroutine print_res2D(name,Temp,Np,curTime)
        use CustomDouble
        implicit none

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J

        real(dp) :: curTime
        real(dp) :: Temp(0:200,0:200)
        real(dp) :: Np(1:200,1:200,1:2)
        character(LEN=50) :: name

        print*, "========================================================"
        print*, "Results for 2D calculation will be printed out on screen"
        print*, "              ", name, "                                "
        print*, "Time = ", curTime, " seconds"
        print*, "========================================================"
        print*,""
        print*,""

100     FORMAT(' ', A,I3,A,F7.2)
200     FORMAT('',F7.2,A)

        J = 1
        do I = 1, N(1), 1
          write(*,200,advance='no')  Np(I,J,1),"  "
        end do
        print*, ""

        I = 1
        do J = N(2), 1, -1
          write(*,100) "J = ", J, "  ============ Y = ", Np(I,J,2)
          print*, ""

          do I = 1, N(1), 1
            write(*,200,advance='no')  Temp(I,J),"  "
          end do
          print*, ""
          print*, ""
        end do
      end subroutine print_res2D




      subroutine print_k_T2D(name,Temp,Np)
        use CustomDouble

        implicit none

        integer :: N(1:2)
        common /gridsize/ N

        integer ::  I, J

        real(dp) :: Temp(0:200,0:200)
        real(dp) :: Np(1:200,1:200,1:2)
        character(LEN=50) :: name

        print*, "==============================================================="
        print*, "Thermal conductivity for 2D grid will be printed out on screen"
        print*, "              ", name, "                                "
        print*, "==============================================================="
        print*,""
        print*,""

100     FORMAT(' ', A,I3,A,F7.2)
200     FORMAT('',F7.2,A)

        J = 1
        do I = 1, N(1), 1
          write(*,200,advance='no')  Np(I,J,1),"  "
        end do
        print*, ""

        I = 1
        do J = N(2), 1, -1
          write(*,100) "J = ", J, "  ============ Y = ", Np(I,J,2)
          print*, ""

          do I = 1, N(1), 1
            write(*,200,advance='no')  k_Tn(Temp(I,J)),"  "
          end do
          print*, ""
          print*, ""
        end do
      end subroutine print_k_T2D




      subroutine compSol2D(N,TGaSe,TTDMA)
        use CustomDouble

        implicit none

        integer :: N(1:2), I, J
        integer :: Imax,Jmax
        real(dp) :: maxdiff, temp
        real(dp) :: TGaSe(0:200,0:200)
        real(dp) :: TTDMA(0:200,0:200)

        maxdiff = 0
        Imax = 0
        Jmax = 0
        do I = 1, N(1) - 1, 1
          do J = 1, N(2) - 1, 1
            temp = abs(TGaSe(I,J)-TTDMA(I,J))
            if (temp .GT. maxdiff) then
              maxdiff = temp
              Imax = I
              Jmax = J
            end if
          end do
        end do

        print*,"Max difference = ", maxdiff
        print*,"Position =", Imax,",",Jmax
        print*,"Temp Gauss Seidel b: ", TGaSe(Imax,Jmax)
        print*,"Temp TDMA: ", TTDMA(Imax,Jmax)

      end subroutine compSol2D



      subroutine res2Dmod(name,N,Temp,NPx,NPy)
        use CustomDouble

        implicit none

        integer :: N(1:3), I, J, interN, totalN
        real(dp) :: Temp(0:200,0:200)
        real(dp) :: NPx(1:200), NPy(1:200)
        character(LEN=50) :: name

        interN = 0
        totalN = 11
        do while (totalN .NE. N(1))
          interN = interN + 1
          totalN = 11 + interN * 10
        end do

        print*, "=================================================================="
        print*, "  Partial Results for 2D calculation will be printed out on screen"
        print*, "                        ", name, "                                "
        print*, "=================================================================="
        print*,""
        print*,""

100     FORMAT(' ', A,I3,A,F7.2)
200     FORMAT('',F7.2,A)

        do I = 1, N(1), interN + 1
          write(*,200,advance='no')  NPx(I),"  "
        end do
        print*, ""

        do J = N(2), 1, - interN - 1
          write(*,100) "J = ", J, "  ============ Y = ", NPy(J)
          print*, ""

          do I = 1, N(1), interN + 1
            write(*,200,advance='no')  Temp(I,J),"  "
          end do
          print*, ""
          print*, ""
        end do
      end subroutine res2Dmod


      end module Results
