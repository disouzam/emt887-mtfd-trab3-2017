*******************************************************
* FILE: Properties.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      module Properties

            implicit none

            ! Melting temperature in Kelvin
            double precision, parameter :: Tmelt = 1809.15D0

      contains

      function k_Tn(Tn)


        implicit none

        double precision :: k_Tn
        double precision :: Tn

        ! Implementation according Pereira (2004)
        ! Modelamento matem�tico do escoamento turbulento, da transfer�ncia de calor e da solidifica��o no distribuidor
        ! e na m�quina de lingotamento cont�nuo
        ! Rodrigo Ottoni da Silva Pereira - PPGEM - UFMG - 2004 - Disserta��o de mestrado
        if (Tn .LT. 1000.15D0) then
          k_Tn = 64.1354D0 - 0.0427D0 * Tn
        else
          k_Tn = 16.9952D0 + 0.0115D0 * Tn
        end if

      end function k_Tn

      function k_Mushy(Tn)

        implicit none

        double precision :: k_Mushy
        double precision :: Tn
        double precision :: fLiquid
        double precision :: kS, C

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        double precision :: SteelChem(1:6)
        common /steelprop/ SteelChem

        k_Mushy = 0.0D0
        C = 5.0D0

        fLiquid = f_Liq(Tn)
        kS = k_Tn(Tn)

        k_Mushy = kS * (1.0D0 + (C - 1.0D0) * fLiquid * fLiquid)

      end function k_Mushy

      function k_Inter(Ti,Tip1,DelMinus,DelPlus)
        implicit none

        double precision :: k_Inter
        double precision :: Ti, Tip1,DelMinus,DelPlus
        double precision :: dummy,DelTotal

        DelTotal = DelMinus + DelPlus

        dummy = (DelMinus / DelTotal) * (1.0D0 / k_Mushy(Ti))
        dummy = dummy + (DelPlus / DelTotal) * (1.0D0 / k_Mushy(Tip1))

        dummy = 1.0D0 / dummy

        k_Inter = dummy

      end function k_Inter

      function Cp_T(Tn)

        implicit none

        double precision :: Cp_T
        double precision :: Tn
        double precision :: Tliq

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        double precision :: SteelChem(1:6)
        common /steelprop/ SteelChem

        ! Implementation according Pereira (2004)
        ! Modelamento matem�tico do escoamento turbulento, da transfer�ncia de calor e da solidifica��o no distribuidor
        ! e na m�quina de lingotamento cont�nuo
        ! Rodrigo Ottoni da Silva Pereira - PPGEM - UFMG - 2004 - Disserta��o de mestrado

        Tliq = Tliquidus()

        if (Tn .LE. 1020.0D0) then

          Cp_T = 2368.0D0
          Cp_T = Cp_T - 14.92D0 * Tn
          Cp_T = Cp_T + 4.107D-2 * Tn * Tn
          Cp_T = Cp_T - 4.696D-5 * Tn * Tn * Tn
          Cp_T = Cp_T + 1.953D-8 * Tn * Tn * Tn * Tn

        elseif (Tn .LE. 1210.0D0) then

          Cp_T = 7802.0D0
          Cp_T = Cp_T - 5.278D0 * Tn
          Cp_T = Cp_T - 3.676D-3 * Tn * Tn
          Cp_T = Cp_T + 1.388D-6 * Tn * Tn * Tn
          Cp_T = Cp_T + 1.031D-9 * Tn * Tn * Tn * Tn

        elseif  (Tn .LT. Tliq) then

          Cp_T = 703

        else

          ! Adapta��o de Suzuki et al.
          Cp_T = 720

        end if

      end function Cp_T

      function Cp_Eq(Tn)
        implicit none

        double precision :: Cp_Eq
        double precision :: Tn, Tliq, Tsol

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        Cp_Eq = 0.0D0

        Tliq = Tliquidus()
        Tsol = Tsolidus()

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

        double precision :: h_Total
        double precision :: h,Tf, Tp
        double precision :: sigma, eps, p1, p2

        sigma = 5.67D-8
        eps = 0.80D0
        p1 = Tf + Tp
        p2 = Tf * Tf + Tp * Tp

        h_Total = h + sigma * eps * p1 * p2

      end function h_Total

      function emissiv(Tp)

        implicit none

        double precision :: emissiv
        double precision :: Tp
        double precision :: y1, y2, x1, x2, m

        if (Tp .LT. 590.0D0) then
          emissiv = 0.69D0
        elseif (Tp .LT. 755.0D0) then
          y1 = 0.69D0
          x1 = 590.0D0
          y2 = 0.72D0
          x2 = 755.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 920.0D0) then
          y1 = 0.72D0
          x1 = 755.0D0
          y2 = 0.76D0
          x2 = 920.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1090.0D0) then
          y1 = 0.76D0
          x1 = 920.0D0
          y2 = 0.79D0
          x2 = 1090.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        elseif (Tp .LT. 1255.0D0) then
          y1 = 0.79D0
          x1 = 1090.0D0
          y2 = 0.82D0
          x2 = 1255.0D0
          m = (y2 - y1) / (x2 - x1)
          emissiv = m * (Tp - x1) + y1
        else
          emissiv = 0.82D0
        end if


      end function emissiv

      function f_Liq_Lever(T)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: f_Liq_Lever
        double precision :: T, Tliq, Tsol

        double precision :: f, k

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        if (T .LT. Tsol) then
          f = 0.0D0
        elseif (T .GT. Tliq) then
          f = 1.0D0
        else
          k = (Tmelt - Tliq) / (Tmelt - Tsol)
          f = (1.0D0) / (1.0D0 - k)
          f = f * (T - Tliq) / (T - Tmelt)
          f = 1.0D0 - f
        end if

        f_Liq_Lever = f

      end function f_Liq_Lever

      function df_dTLever(Tn)
        implicit none

        double precision :: df_dTLever
        double precision :: Tn, Tliq, Tsol, k, d

        Tliq = Tliquidus()
        Tsol = Tsolidus()
        k = (Tmelt - Tliq) / (Tmelt - Tsol)

        d = Tmelt - Tliq
        d = d / (1.0D0 - k)
        d = d / ((Tmelt - Tn) * (Tmelt - Tn) )

        df_dTLever = d

      end function

      function f_Linear(T)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: f_Linear
        double precision :: T, Tliq, Tsol
        double precision :: f

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        if (T .LT. Tsol) then
          f = 0.0D0
        elseif (T .GT. Tliq) then
          f = 1.0D0
        else
          f = (T - Tsol) / (Tliq - Tsol)
        end if

        f_Linear = f

      end function f_Linear

      function df_dTLin()
        implicit none

        double precision :: df_dTlin
        double precision :: Tliq, Tsol

        Tliq = Tliquidus()
        Tsol = Tsolidus()

        df_dTlin = (1.0D0) / (Tliq - Tsol)

      end function df_dTLin

      function Tliquidus()
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: Tliquidus

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        double precision :: SteelChem(1:6)
        common /steelprop/ SteelChem

        double precision :: T

        T = 78.0D0 * SteelChem(1)
        T = T + 4.9D0 * SteelChem(2)
        T = T + 7.6D0* SteelChem(3)
        T = T + 34.4D0 * SteelChem(4)
        T = T + 38.0D0 * SteelChem(5)
        T = T + 3.6D0 * SteelChem(6)

        T = Tmelt - T
        Tliquidus = T

      end function Tliquidus

      function Tsolidus()
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: Tsolidus

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        double precision :: SteelChem(1:6)
        common /steelprop/ SteelChem

        double precision :: T

        T = 415.5D0 * SteelChem(1)
        T = T + 6.8D0 * SteelChem(2)
        T = T + 12.3D0* SteelChem(3)
        T = T + 124.5D0 * SteelChem(4)
        T = T + 183.9D0 * SteelChem(5)
        T = T + 4.1D0 * SteelChem(6)

        T = Tmelt - T
        Tsolidus = T

      end function Tsolidus

      function Q_mould(time)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        double precision :: Q_mould
        double precision :: time

        double precision :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        if (time .LE. tMould) then
          Q_mould = 1000.0D0 * (2680.0D0 - 220.0D0 * sqrt(time))
        else
          Q_mould = 0.0D0
        end if

        Q_mould = - Q_mould

      end function Q_mould

      function h_SecCol(time,Temp)
        implicit none

        double precision :: h_SecCol
        double precision :: time, Temp

        double precision :: W1, W2, H3
        common /secCooling/ W1, W2, H3

        double precision :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        double precision :: CumTimeM, CumTimeS1, CumTimeS2, CumTimeS3

        CumTimeM = tMould
        CumTimeS1 = CumTimeM + tS1
        CumTimeS2 = CumTimeS1 + tS2
        CumTimeS3 = CumTimeS2 + tS3

        h_SecCol = 0.0D0

        if (time .LE. CumTimeM) then
          h_SecCol = 0.0D0
        elseif (time .LE. CumTimeS1) then
          h_SecCol = 0.116
          h_SecCol = h_SecCol + 708.0D0 * (W1**0.75D0) * (Temp**(-1.2D0))
          h_SecCol = h_SecCol * 1000.0D0
        elseif (time .LE. CumTimeS2) then
          h_SecCol = 0.116
          h_SecCol = h_SecCol + 708.0D0 * (W2**0.75D0) * (Temp**(-1.2D0))
          h_SecCol = h_SecCol * 1000.0D0
        elseif (time .LE. CumTimeS3) then
          h_SecCol = H3
        end if

      end function h_SecCol

      function LatentHeat()
        implicit none
        double precision :: LatentHeat

        LatentHeat = 283000.0D0

      end function LatentHeat


      function f_Liq(Tn)
        implicit none

        double precision :: f_Liq
        double precision :: Tn

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        f_Liq = 0.0D0

        if (mode .EQ. 1) then
          f_Liq = f_Liq_Lever(Tn)
        elseif (mode .EQ. 2) then
          f_Liq = f_Linear(Tn)
        end if

      end function

      end module Properties
