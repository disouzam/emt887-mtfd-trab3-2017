*******************************************************
* FILE: Coefficients.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Coefficients
        use Properties

        implicit none

      contains

      subroutine coeff2D(Aw,Ae,An,As,Ap,b,Delta_t)

        use CustomDouble

* Calculates matrix coefficients for 2D problems
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        real(dp) :: tMould, tS1, tS2, tS3,CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf
        common /boundCondLen/ mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf

*       Stores heat transfer coefficient for secondary cooling region
        real(dp) :: H

*       Stores the heat flow rate across interface in mould faces
        real(dp) :: Qh

        real(dp) :: Delta_t

        real(dp) :: Np(1:999,1:999,1:2), NI(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        real(dp) :: dens
        common /density/ dens

        real(dp) :: DelMinus, DelPlus, Ti, Tip1, Ap0, Ap1, Cp_New, Cp_Old, hcalc

        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Nx,Ny

        real(dp) :: Del(1:999,1:999,1:2), Delta(1:999,1:999,1:2)
        common / mesh/ Del, Delta

        real(dp) :: Aw(1:999,1:999), Ae(1:999,1:999)
        real(dp) :: An(1:999,1:999), As(1:999,1:999)
        real(dp) :: Ap(1:999,1:999), b(1:999,1:999)

        character(LEN=50) :: section

        real(dp) :: dummy, maxDumm, maxCp

        Nx = N(1)
        Ny = N(2)

        ! Implementation according Pereira (2004)
        ! Modelamento matem�tico do escoamento turbulento, da transfer�ncia de calor e da solidifica��o no distribuidor
        ! e na m�quina de lingotamento cont�nuo
        ! Rodrigo Ottoni da Silva Pereira - PPGEM - UFMG - 2004 - Disserta��o de mestrado
        dens = 7000.0_dp

        maxDumm = 1.0E50_dp
        maxCp = 0.0_dp

*        if (debugmode .EQV. .TRUE.) then
*          print*,"Density = ", dens
*        end if

        Qh = Q_mould(curTime)

        section = "Left Face"
*       Equations for points at left face (excluding corners at top and bottom)
        I = 1
        do J = 2, Ny - 1, 1

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)

          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

          Aw(I,J) = Ae(I,J)
          Tnew(I-1,J) = Tnew(I+1,J)

          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)

          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)

          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
        end do

        section = "Right Face"
*       Equations for points at right face (excluding corners at top and bottom)
        I = Nx
        do J = 2, Ny - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)

          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Ae(I,J) = 0.0_dp
          dummy = Ae(I,J)

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)

          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)

          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          if (curTime .GT. tMould) then
            H = h_SecCol(curTime,Tnew(I,J))
            hcalc = h_Total(H,Tf,Tnew(I,J))
          else
            hcalc = 0.0_dp
          end if


          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,2)

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,2) * Tf
          b(I,J) = b(I,J) + Qh * Delta(I,J,2)
        end do


        section = "Bottom Face"
*       Equations for points at bottom face (excluding left and right corners)
        J = 1
        do I = 2, Nx - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)

          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)

          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)

          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

          As(I,J) = An(I,J)
          Tnew(I,J-1) = Tnew(I,J+1)

          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
        end do


        section = "Top Face"
*       Equations for points at top face (excluding left and right corners)
        J = Ny
        do I = 2, Nx - 1, 1

          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)

          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)

          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)

          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          An(I,J) = 0.0_dp

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          if (curTime .GT. tMould) then
            H = h_SecCol(curTime,Tnew(I,J))
            hcalc = h_Total(H,Tf,Tnew(I,J))
          else
            hcalc = 0.0_dp
          end if

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,1)

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,1) * Tf
          b(I,J) = b(I,J) + Qh * Delta(I,J,1)
        end do

        section = "Top Left Corner"
        J = Ny
*       Equation for top left corner
        do I = 1, 1, 1
          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)

          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

          Aw(I,J) = Ae(I,J)
          Tnew(I-1,J) = Tnew(I+1,J)

          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)

          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          An(I,J) = 0.0_dp

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          if (curTime .GT. tMould) then
            H = h_SecCol(curTime,Tnew(I,J))
            hcalc = h_Total(H,Tf,Tnew(I,J))
          else
            hcalc = 0.0_dp
          end if

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,1)

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,1) * Tf
          b(I,J) = b(I,J) + Qh * Delta(I,J,1)
        end do

        section = "Top Right Corner"
        J = Ny
*       Equation for top right corner
        do I = Nx, Nx, 1
          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)

          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J,2) - Np(I,J-1,2)
          DelPlus = Np(I,J,2) - Ni(I,J,2)
          Ti = Tnew(I,J-1)
          Tip1 = Tnew(I,J)

          As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

          dummy = As(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Ae(I,J) = 0.0_dp
          An(I,J) = 0.0_dp

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          if (curTime .GT. tMould) then
            H = h_SecCol(curTime,Tnew(I,J))
            hcalc = h_Total(H,Tf,Tnew(I,J))
          else
            hcalc = 0.0_dp
          end if

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,2)
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,1)

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,2) * Tf
          b(I,J) = b(I,J) + hcalc * Delta(I,J,1) * Tf
          b(I,J) = b(I,J) + Qh * Delta(I,J,2)
          b(I,J) = b(I,J) + Qh * Delta(I,J,1)
        end do

        section = "Bottom Left Corner"
        J = 1
*       Equation for bottom left corner
        do I = 1, 1, 1
          DelMinus = Ni(I+1,J,1) - Np(I,J,1)
          DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I+1,J)

          Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

          Aw(I,J) = Ae(I,J)
          Tnew(I-1,J) = Tnew(I+1,J)

          dummy = Ae(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)

          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

          As(I,J) = An(I,J)
          Tnew(I,J-1) = Tnew(I,J+1)

          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
        end do

        section = "Bottom Right Corner"
        J = 1
*       Equation for bottom right corner
        do I = Nx, Nx, 1
          DelMinus = Ni(I,J,1) - Np(I-1,J,1)
          DelPlus = Np(I,J,1) - Ni(I,J,1)
          Ti = Tnew(I-1,J)
          Tip1 = Tnew(I,J)

          Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

          dummy = Aw(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          DelMinus = Ni(I,J+1,2) - Np(I,J,2)
          DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
          Ti = Tnew(I,J)
          Tip1 = Tnew(I,J+1)

          An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

          As(I,J) = An(I,J)
          Tnew(I,J-1) = Tnew(I,J+1)

          dummy = An(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          Ae(I,J) = 0.0_dp

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
          Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          Cp_New = Cp_Eq(Tnew(I,J))
          !if (Cp_New .GT. maxCp) maxCp = Cp_New
          Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

          if (curTime .GT. tMould) then
            H = h_SecCol(curTime,Tnew(I,J))
            hcalc = h_Total(H,Tf,Tnew(I,J))
          else
            hcalc = 0.0_dp
          end if

          Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1
          Ap(I,J) = Ap(I,J) + hcalc * Delta(I,J,2)

          dummy = Ap(I,J)
          if (dummy .GT. maxDumm) then
            print*,"Error"
            ! read*
          end if

          b(I,J) = Ap0 * Told(I,J)
          b(I,J) = b(I,J) + hcalc * Delta(I,J,2) * Tf
          b(I,J) = b(I,J) + Qh * Delta(I,J,2)
        end do

        section = "Inner Points"
*       Equations for inner points
        do I = 2, Nx - 1, 1
          do J = 2, Ny - 1, 1

            DelMinus = Ni(I,J,1) - Np(I-1,J,1)
            DelPlus = Np(I,J,1) - Ni(I,J,1)
            Ti = Tnew(I-1,J)
            Tip1 = Tnew(I,J)

            Aw(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I-1,J,1) ) * Delta(I,J,2)

            dummy = Aw(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
              ! read*
            end if

            DelMinus = Ni(I+1,J,1) - Np(I,J,1)
            DelPlus = Np(I+1,J,1) - Ni(I+1,J,1)
            Ti = Tnew(I,J)
            Tip1 = Tnew(I+1,J)

            Ae(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,1) ) * Delta(I,J,2)

            dummy = Ae(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
              ! read*
            end if

            DelMinus = Ni(I,J,2) - Np(I,J-1,2)
            DelPlus = Np(I,J,2) - Ni(I,J,2)
            Ti = Tnew(I,J-1)
            Tip1 = Tnew(I,J)

            As(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J-1,2) ) * Delta(I,J,1)

            dummy = As(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
              ! read*
            end if

            DelMinus = Ni(I,J+1,2) - Np(I,J,2)
            DelPlus = Np(I,J+1,2) - Ni(I,J+1,2)
            Ti = Tnew(I,J)
            Tip1 = Tnew(I,J+1)

            An(I,J) = (k_Inter(Ti,Tip1,DelMinus,DelPlus)  / Del(I,J,2) ) * Delta(I,J,1)

            dummy = An(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
              ! read*
            end if

          Cp_Old = Cp_Eq(Told(I,J))
*          Cp_Old = Cp_Eq(Tnew(I,J))
            Ap0 = (dens * Cp_Old * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

            Cp_New = Cp_Eq(Tnew(I,J))
            !if (Cp_New .GT. maxCp) maxCp = Cp_New
            Ap1 = (dens * Cp_New * Delta(I,J,1) * Delta(I,J,2)) / Delta_t

            Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J) + Ap1

            dummy = Ap(I,J)
            if (dummy .GT. maxDumm) then
              print*,"Error"
              ! read*
            end if

            b(I,J) = Ap0 * Told(I,J)

          end do
        end do

*107     format(' ', A, 1F12.6)
*        print 107,"MaxCp =           ", maxCp

      end subroutine coeff2D

      end module Coefficients
