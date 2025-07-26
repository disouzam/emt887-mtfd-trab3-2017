*******************************************************
* FILE: Properties.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Properties

            use CustomDouble

            implicit none

            ! Melting temperature in Kelvin
            real(dp), parameter :: Tmelt = 1809.15_dp

      contains

      function k_Tn(Tn)


        implicit none

        real(dp) :: k_Tn
        real(dp) :: Tn

*       Implementation for AISI 1018 steel according to
*         KIM, Y.;FAROUK, B.; KEVERIAN, J. A mathematical model for Thermal Analysis of Thin Strip Casting of Low Carbon Steel
*         Journal of Engineering for Industry, February 1991, v. 113
*        k_Tn = 65.214_dp
*        k_Tn = k_Tn + 2.715E-2_dp * Tn
*        k_Tn = k_Tn - 1.628E-4_dp * Tn * Tn
*        k_Tn = k_Tn + 1.39E-7_dp * Tn * Tn * Tn
*        k_Tn = k_Tn - 3.041E-11_dp * Tn * Tn * Tn * Tn


        ! Implementation according Pereira (2004)
        ! Modelamento matem�tico do escoamento turbulento, da transfer�ncia de calor e da solidifica��o no distribuidor
        ! e na m�quina de lingotamento cont�nuo
        ! Rodrigo Ottoni da Silva Pereira - PPGEM - UFMG - 2004 - Disserta��o de mestrado
        if (Tn .LT. 1000.15_dp) then

          k_Tn = 64.1354_dp - 0.0427_dp * Tn

        else

          k_Tn = 16.9952_dp + 0.0115_dp * Tn

        end if

      end function k_Tn

      function k_Mushy(Tn)

        implicit none

        real(dp) :: k_Mushy
        real(dp) :: Tn
        real(dp) :: fLiq
        real(dp) :: kS, C

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiquidFraction/ mode

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        k_Mushy = 0.0_dp
        C = 5.0_dp

        if (mode .EQ. 1) then
          fLiq = f_Liq_Lever(Tn)
          kS = k_Tn(Tn)

          k_Mushy = kS * (1.0_dp + (C - 1.0_dp) * fLiq * fLiq)

        elseif (mode .EQ. 2) then
          fLiq = f_Linear(Tn)
          kS = k_Tn(Tn)

          k_Mushy = kS * (1.0_dp + (C - 1.0_dp) * fLiq * fLiq)
        end if

      end function k_Mushy

      function k_Inter(Ti,Tip1,DelMinus,DelPlus)
        implicit none

        real(dp) :: k_Inter
        real(dp) :: Ti, Tip1,DelMinus,DelPlus
        real(dp) :: dummy,DelTotal

        DelTotal = DelMinus + DelPlus

        dummy = (DelMinus / DelTotal) * (1.0_dp / k_Mushy(Ti))
        dummy = dummy + (DelPlus / DelTotal) * (1.0_dp / k_Mushy(Tip1))

        dummy = 1.0_dp / dummy

        k_Inter = dummy

      end function k_Inter

      function Cp_T(Tn)

        implicit none

        real(dp) :: Cp_T
        real(dp) :: Tn
        real(dp) :: Tliq

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        ! Implementation for AISI 1018 steel according to
        ! KIM, Y.;FAROUK, B.; KEVERIAN, J. A mathematical model for Thermal Analysis of Thin Strip Casting of Low Carbon Steel
        ! Journal of Engineering for Industry, February 1991, v. 113

*        if (Tn .LT. 1033.0_dp) then
*
*          Cp_T = 2.368_dp
*          Cp_T = Cp_T - 1.492E-2_dp * Tn
*          Cp_T = Cp_T + 4.107E-5_dp * Tn * Tn
*          Cp_T = Cp_T - 4.696E-8_dp * Tn * Tn * Tn
*          Cp_T = Cp_T + 1.953E-11_dp * Tn * Tn * Tn * Tn
*          Cp_T = 1000.0_dp * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K
*
*        elseif (Tn .LT. 1200.0_dp) then
*
*          Cp_T = 7.802_dp
*          Cp_T = Cp_T - 5.278E-3_dp * Tn
*          Cp_T = Cp_T - 3.676E-6_dp * Tn * Tn
*          Cp_T = Cp_T + 1.388E-9_dp * Tn * Tn * Tn
*          Cp_T = Cp_T + 1.031E-12_dp * Tn * Tn * Tn * Tn
*          Cp_T = 1000.0_dp * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K
*
*        elseif (Tn .LT. 1776.0_dp) then
*
*          Cp_T = 0.703_dp
*          Cp_T = 1000.0_dp * Cp_T ! Conversion from KJ/Kg.K to J/Kg.K
*
*        end if

        ! Implementation according Pereira (2004)
        ! Modelamento matem�tico do escoamento turbulento, da transfer�ncia de calor e da solidifica��o no distribuidor
        ! e na m�quina de lingotamento cont�nuo
        ! Rodrigo Ottoni da Silva Pereira - PPGEM - UFMG - 2004 - Disserta��o de mestrado

        Tliq = Tliquidus()

        if (Tn .LE. 1020.0_dp) then

          Cp_T = 2368.0_dp
          Cp_T = Cp_T - 14.92_dp * Tn
          Cp_T = Cp_T + 4.107E-2_dp * Tn * Tn
          Cp_T = Cp_T - 4.696E-5_dp * Tn * Tn * Tn
          Cp_T = Cp_T + 1.953E-8_dp * Tn * Tn * Tn * Tn

        elseif (Tn .LE. 1210.0_dp) then

          Cp_T = 7802.0_dp
          Cp_T = Cp_T - 5.278_dp * Tn
          Cp_T = Cp_T - 3.676E-3_dp * Tn * Tn
          Cp_T = Cp_T + 1.388E-6_dp * Tn * Tn * Tn
          Cp_T = Cp_T + 1.031E-9_dp * Tn * Tn * Tn * Tn

        elseif  (Tn .LT. Tliq) then

          Cp_T = 703

        else

          ! Adapta��o de Suzuki et al.
          Cp_T = 720


        end if

      end function Cp_T

      function Cp_Eq(Tn)
        implicit none

        real(dp) :: Cp_Eq
        real(dp) :: Tn, Tliq, Tsol

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiquidFraction/ mode

        Cp_Eq = 0.0_dp

        Tliq = Tliquidus
        Tsol = Tsolidus

        if ((Tn .LE. Tliq) .AND. (Tn .GE. Tsol)) then

          if (mode .EQ. 1) then
            Cp_Eq = Cp_T(Tn) + LatentHeat() * df_dTLever(Tn)
          end if

          if (mode .EQ. 2) then
            Cp_Eq = Cp_T(Tn) + LatentHeat() * df_dTLin()
          end if

        else
          Cp_Eq = Cp_T(Tn)
        end if

      end function Cp_Eq


      function h_Total(h,Tf,Tp)

        implicit none

        real(dp) :: h_Total
        real(dp) :: h,Tf, Tp
        real(dp) :: sigma, eps, p1, p2

        sigma = 5.67E-8_dp
*        eps = emissiv(Tp)
        eps = 0.80_dp
        p1 = Tf + Tp
        p2 = Tf * Tf + Tp * Tp

        h_Total = h + sigma * eps * p1 * p2

      end function h_Total

      function emissiv(Tp)

        implicit none

        real(dp) :: emissiv
        real(dp) :: Tp
        real(dp) :: y1, y2, x1, x2, m

        if (Tp .LT. 590.0_dp) then
          emissiv = 0.69_dp
        elseif (Tp .LT. 755.0_dp) then
          y1 = 0.69_dp
          x1 = 590.0_dp
          y2 = 0.72_dp
          x2 = 755.0_dp
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 920.0_dp) then
          y1 = 0.72_dp
          x1 = 755.0_dp
          y2 = 0.76_dp
          x2 = 920.0_dp
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1090.0_dp) then
          y1 = 0.76_dp
          x1 = 920.0_dp
          y2 = 0.79_dp
          x2 = 1090.0_dp
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1255.0_dp) then
          y1 = 0.79_dp
          x1 = 1090.0_dp
          y2 = 0.82_dp
          x2 = 1255.0_dp
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        else
          emissiv = 0.82_dp
        end if


      end function emissiv

      function f_Liq_Lever(T)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: f_Liq_Lever
        real(dp) :: T, Tliq, Tsol

        real(dp) :: f, k

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        if (T .LT. Tsol) then
          f = 0.0_dp
        elseif (T .GT. Tliq) then
          f = 1.0_dp
        else
          k = (Tmelt - Tliq) / (Tmelt - Tsol)
          f = (1.0_dp) / (1.0_dp - k)
          f = f * (T - Tliq) / (T - Tmelt)
          f = 1.0_dp - f
        end if

        f_Liq_Lever = f

        if (debugmode .EQV. .TRUE.) then
          if (f_Liq_Lever .LT. 1.0_dp) then
*            print*,"f_Liq_Lever =", f_Liq_Lever
          end if
        end if

      end function f_Liq_Lever

      function df_dTLever(Tn)
        implicit none

        real(dp) :: df_dTLever
        real(dp) :: Tn, Tliq, Tsol, k, d

        Tliq = Tliquidus
        Tsol = Tsolidus
        k = (Tmelt - Tliq) / (Tmelt - Tsol)

        d = Tmelt - Tliq
        d = d / (1.0_dp - k)
        d = d / ((Tmelt - Tn) * (Tmelt - Tn) )

        df_dTLever = d

      end function

      function f_Linear(T)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: f_Linear
        real(dp) :: T, Tliq, Tsol
        real(dp) :: f

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        if (T .LT. Tsol) then
          f = 0.0_dp
        elseif (T .GT. Tliq) then
          f = 1.0_dp
        else
          f = (T - Tsol) / (Tliq - Tsol)
        end if

        f_Linear = f

        if (debugmode .EQV. .TRUE.) then
*          print*,"f_Linear =", f_Linear
        end if

      end function f_Linear

      function df_dTLin()
        implicit none

        real(dp) :: df_dTlin
        real(dp) :: Tliq, Tsol

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        df_dTlin = (1.0_dp) / (Tliq - Tsol)

      end function df_dTLin

      function Tliquidus()
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: Tliquidus

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        real(dp) :: T

        T = 78.0_dp * SteelChem(1)
        T = T + 4.9_dp * SteelChem(2)
        T = T + 7.6_dp* SteelChem(3)
        T = T + 34.4_dp * SteelChem(4)
        T = T + 38.0_dp * SteelChem(5)
        T = T + 3.6_dp * SteelChem(6)

        T = Tmelt - T

        if (debugmode .EQV. .TRUE.) then
          ! print*,"Tliquidus = ", T
        end if

        Tliquidus = T

      end function Tliquidus

      function Tsolidus()
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: Tsolidus

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        real(dp) :: T

        T = 415.5_dp * SteelChem(1)
        T = T + 6.8_dp * SteelChem(2)
        T = T + 12.3_dp* SteelChem(3)
        T = T + 124.5_dp * SteelChem(4)
        T = T + 183.9_dp * SteelChem(5)
        T = T + 4.1_dp * SteelChem(6)

        T = Tmelt - T

        if (debugmode .EQV. .TRUE.) then
          ! print*,"Tsolidus = ", T
        end if

        Tsolidus = T

      end function Tsolidus

      function Q_mould(time)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: Q_mould
        real(dp) :: time

        real(dp) :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        if (time .LE. tMould) then
          Q_mould = 1000.0_dp * (2680.0_dp - 220.0_dp * sqrt(time))
        else
          Q_mould = 0.0_dp
        end if

        Q_mould = - Q_mould

        if (debugmode .EQV. .TRUE.) then
*          print*,"Qmould = ", Q_mould
          ! read*
        end if

      end function Q_mould

      function h_SecCol(time,Temp)
        implicit none

        real(dp) :: h_SecCol
        real(dp) :: time, Temp

        real(dp) :: W1, W2, H3
        common /secCooling/ W1, W2, H3

        real(dp) :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: CumTimeM, CumTimeS1, CumTimeS2, CumTimeS3

        CumTimeM = tMould
        CumTimeS1 = CumTimeM + tS1
        CumTimeS2 = CumTimeS1 + tS2
        CumTimeS3 = CumTimeS2 + tS3

        h_SecCol = 0.0_dp

        if (time .LE. CumTimeM) then
          h_SecCol = 0.0_dp
        elseif (time .LE. CumTimeS1) then
          h_SecCol = 0.116
          h_SecCol = h_SecCol + 708.0_dp * (W1**0.75_dp) * (Temp**(-1.2_dp))
          h_SecCol = h_SecCol * 1000.0_dp
        elseif (time .LE. CumTimeS2) then
          h_SecCol = 0.116
          h_SecCol = h_SecCol + 708.0_dp * (W2**0.75_dp) * (Temp**(-1.2_dp))
          h_SecCol = h_SecCol * 1000.0_dp
        elseif (time .LE. CumTimeS3) then
          h_SecCol = H3
        end if

      end function h_SecCol

      function LatentHeat()
        implicit none

        real(dp) :: LatentHeat

        LatentHeat = 283000.0_dp

      end function LatentHeat
      end module Properties
