*******************************************************
* FILE: main.f
* Author: Dickson Alves de Souza
* Professor: Roberto Parreiras Tavares - UFMG
*******************************************************
      program Trab3

*Modules
        use CustomDouble
        use Coefficients
        use Geometry
        use Results
        use Solver
        use Properties

*Settings
        implicit none

*Variable declaration

        logical :: debugmode
        common /dbgMode/ debugmode

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        real(dp) :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: maxX, maxY, Xi(9), Yi(9), tsSavPos
        common /storeLimits/ maxX, maxY, Xi, Yi, tsSavPos

*       Del: array containing distance between neighboring nodes in all axis(to be calculated)
*       Index ranges from 1 to N(1) - 1 in X direction, from 1 to N(2) - 1 in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Delta: array containing distance between parallel interfaces of a volume element in all axis (to be calculated)
*       Index ranges from 1 to N(1)  in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Del(1:999,1:999,1:2),Delta(1:999,1:999,1:2)
        common /mesh/ Del, Delta

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        real(dp) :: mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf
        common /boundCondLen/ mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf

        real(dp) :: W1, W2, H3
        common /secCooling/ W1, W2, H3

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        real(dp) :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 50000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Specific heat
*          6: Thermal conductivity
*          7: Liquid fraction
        real(dp) :: TmRs(1:81,0:50000,1:7)
        common /transient/ TmRs

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        integer :: meth,relax_m, ADI
        real(dp) :: tol, alpha
        common /solMethod/ tol, alpha, meth, relax_m, ADI

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        real(dp) :: dens
        common /density/ dens
*       Tolerance value for iterative methods - it defines the stopping criteria
*       It is an absolute value, i.e. it is intented to be the maximum value
*         for the difference between temperatures at successive iterations
*       alpha: A constant to be used as a factor for relaxation in Gauss-Seidel and TDMA methods

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        mode = 2

        SteelChem(1) = 0.08_dp
        SteelChem(2) = 0.42_dp
        SteelChem(3) = 0.05_dp
        SteelChem(4) = 0.008_dp
        SteelChem(5) = 0.005_dp
        SteelChem(6) = 0.00_dp

*       In debug mode (=.TRUE.), some info is displayed in screen
        debugmode = .TRUE.

*       Header of program
        call printHeader()

* Instructions
        call readGeometry()

        print*
        print*
        if (debugmode .EQV. .TRUE.) then
          ! read*
        end if

        call readInitTemp()

        print*
        print*
        if (debugmode .EQV. .TRUE.) then
          ! read*
        end if

        call readBoundCond()

        print*
        print*
        if (debugmode .EQV. .TRUE.) then
          ! read*
        end if

        call readTimeInfo()

        print*
        print*
        if (debugmode .EQV. .TRUE.) then
          ! read*
        end if

        call readSolverSel()
        if (debugmode .EQV. .TRUE.) then
          ! read*
        end if

        Tnew(:,:) = Told(:,:)
        Tmax = ToldMax
        Tmin = ToldMin
        Tavg = ToldAvg
        tsSavPos = 0
        call SaveResults(0)
        call defineLimits()
        call storeTimeSteps(0.0_dp)

*        call readBoundCond()

        curTimeStep = 0
        totalIter = 0

        call solveStage()

        call saveTimeSteps()

        print*,"Simulation has finished!"
        print*,"Total iterations: ", totalIter
        read*

      end program Trab3



      subroutine readGeometry()
*       Read geometry parameters and generate mesh through subroutine Grid2DUni

        use CustomDouble
        use Geometry

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        real(dp) :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        print*,"-------------------------------------------------------------------------------"
        print*,"What are the DIMENSIONS of continuous casting strand cross-section(in meters)?"
        print*,"-------------------------------------------------------------------------------"

        print*,"Length in X direction (in meters):"
        if (debugmode .EQV. .TRUE.) then
          Lth(1) = 1.0_dp
          print*,"Lth(1) = ", Lth(1), " m"
        else
          read*, Lth(1)
        end if

        print*
        print*,"Length in Y direction (in meters):"
        if (debugmode .EQV. .TRUE.) then
          Lth(2) = 0.20_dp
          print*,"Lth(2) = ", Lth(2), " m"
        else
          read*, Lth(2)
        end if

        print*
        print*
        print*,"----------------------------------------------------------"
        print*,"What are the NUMBER OF NODES in each direction?"
        print*,"----------------------------------------------------------"

        print*,"Nodes in X direction: (max. 995 nodes)"
        if (debugmode .EQV. .TRUE.) then
          N(1) = 101
          print*,"N(1) = ", N(1), " nodes"
        else
          read*, N(1)
        end if

        print*
        print*,"Nodes in Y direction: (max. 995 nodes)"
        if (debugmode .EQV. .TRUE.) then
          N(2) = 21
          print*,"N(2) = ", N(2), " nodes"
        else
          read*, N(2)
        end if

*       Generate a 2D grid for solution
        call Grid2DUni(Lth,Np,Ni)

      end subroutine readGeometry




      subroutine readBoundCond()
        ! Read boundary conditions for the problem
        use CustomDouble

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode


102     format(' ', A, I3)

        real(dp) :: mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf
        common /boundCondLen/ mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf

        real(dp) :: tMould, tS1, tS2, tS3,CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: W1, W2, H3
        common /secCooling/ W1, W2, H3

        print*,"==================================================================================================================="
        print*,"                                        BOUNDARY CONDITIONS FOR CONTINUOUS CASTING                                 "
        print*,"==================================================================================================================="
        print*

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the MOULD LENGTH? (in meters)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, mouldLen
        else
          mouldLen = 0.8_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the LENGTH of FIRST ZONE of secondary cooling? (in meters)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, SecColLen1
        else
          SecColLen1 = 0.21_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the specific WATER FLOW RATE of FIRST ZONE of secondary cooling? (in liters/sq-m.C)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, W1
        else
          W1 = 5.56_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the LENGTH of SECOND ZONE of secondary cooling? (in meters)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, SecColLen2
        else
          SecColLen2 = 1.85_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the specific WATER FLOW RATE of SECOND ZONE of secondary cooling? (in liters/sq-m.C)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, W2
        else
          W2 = 0.83_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the LENGTH of THIRD ZONE of secondary cooling? (in meters)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, SecColLen3
        else
          SecColLen3 = 13.26_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the HEAT TRANSFER COEFFICIENT of THIRD ZONE of secondary cooling? (in W/sq-m.C)"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, H3
        else
          H3 = 15.0_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the CASTING SPEED (in meters/minute)?"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, CasSpd
          CasSpd = CasSpd / 60.0_dp
        else
          CasSpd = 1.0_dp / 60.0_dp
        end if

        print*,""
        print*,""

        print*,"-----------------------------------------------------------------------------------------------"
        print*,"What is the WATER TEMPERATURE of cooling sprays (in Kelvin)?"
        print*,"-----------------------------------------------------------------------------------------------"

        if (debugmode .EQV. .FALSE.) then
          read*, Tf
        else
          Tf = 298.0_dp
        end if

        tMould = mouldLen / CasSpd
        tS1 = SecColLen1 / CasSpd
        tS2 = SecColLen2 / CasSpd
        ts3 = SecColLen3 / CasSpd
        CasTime = tMould + tS1 + tS2 + tS3

      end subroutine readBoundCond




      subroutine printHeader()
        ! Print the header for the main program

        implicit none

        print*,"==============================================================="
        print*,"        Program for calculation of a 2D thermal profile        "
        print*,"          in cross section of a slab or billet in              "
        print*,"                     Continuous Casting                        "
        print*,"                                                               "
        print*,"  Heat Transfer is calculated considering only heat conduction "
        print*," and convection in liquid metal is taken into consideration by "
        print*," artificially increasing thermal conductivity                  "
        print*,"                                                               "
        print*," Heat released by solidification is modelled by                "
        print*," equivalent heat capacity method.                              "
        print*,"                                                               "
        print*," Boundary Conditions are:                                      "
        print*," 1) In the mold, a specified heat flow as a function of time   "
        print*," 2) In secondary cooling, a specified heat transfer coefficient"
        print*,"    and also radiactive heat transfer.                         "
        print*,"                                                               "
        print*," Velocity of strand is constant and equal to v sub z           "
        print*,"                                                               "
        print*,"==============================================================="
        print*," Author: Dickson Souza                                "
        print*,"                                                               "
        print*," Based on lectures by professor Roberto Parreiras Tavares      "
        print*,"                                                               "
        print*,"      and book Numerical Heat Transfer and Fluid Flow          "
        print*,"      by Suhas V. Patankar (1980)                              "
        print*,"                                                               "
        print*," Federal University of Minas Gerais                            "
        print*," November 18th, 2017                                            "
        print*,"==============================================================="
        print*,"                                                               "
        print*,"                                                               "
        print*,"                                                               "

      end subroutine printHeader




      subroutine readInitTemp()
        ! Read initial temperature

        use CustomDouble
        use Solver

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        real(dp) :: dummy

        print*,"---------------------------------------------------------------"
        print*,"Specify the INITIAL TEMPERATURE of liquid steel (in Kelvin):"
        print*,"---------------------------------------------------------------"

        if (debugmode .EQV. .TRUE.) then
          dummy = 1800.07_dp
          print*,"Initial liquid steel temperature = ", dummy, " K"
        else
          read*, dummy
        end if

        Told(:,:) = dummy

        ! Temperature outside domain was considered equal to their domain counterpart due to symmetry effects
        ! However there are new assignments on coefficient calculations that deals with this condition
        Told(0,:) = dummy
        Told(:,0) = dummy

        Told(N(1)+1,:) = 0.0_dp
        Told(:,N(2)+1) = 0.0_dp

        Imin = 0
        Jmin = 0
        call calcOldStats()

      end subroutine readInitTemp




      subroutine readTimeInfo()
        ! Read information about each stage to be simulated. Up to 10 stages are allowed
        use CustomDouble
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        real(dp) :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: CumTimeM, CumTimeS1, CumTimeS2, CumTimeS3

        real(dp) :: minTime, minTimeSteps

        CumTimeM = tMould
        CumTimeS1 = CumTimeM + tS1
        CumTimeS2 = CumTimeS1 + tS2
        CumTimeS3 = CumTimeS2 + tS3

        print*
        print*

101     format(' ', A, I3, A, I3)
201     format(' ', A, I3, A)
301     format(' ', A, I8)
401     format(' ', F12.0, A)
601     format(' ', A, F8.3, A)


        if (tMould .LT. tS1) then

          if (tMould .LT. tS2) then

            if (tmould .LT. tS3) then
              minTime = tMould
            else
              minTime = tS3
            end if

          else

            if (tS2 .LT. tS3) then
              minTime = tS2
            else
              minTime = tS3
            end if

          end if

        else

          if (tS1 .LT. tS2) then

            if (tS1 .LT. tS3) then
              minTime = tS1
            else
              minTime = tS3
            end if

          else

            if (tS2 .LT. tS3) then
              minTime = tS2
            else
              minTime = tS3
            end if

          end if

        end if

        minTimeSteps = 10 * (CasTime) / (minTime)

        simTime = CasTime
        print*,"---------------------------------------------------------------------------------------------------"
        print 601,"TOTAL simulation time (in seconds): ", simTime
        print*,"                                                                                                   "
        print*,"This time includes time spent in model and outside model "
        print*,"up to the end of THIRD zone of secondary cooling"
        print*,"---------------------------------------------------------------------------------------------------"

        print*,""
        print*,""

        print*,"---------------------------------------------------------------------------------------------------"
        print*,"Which one would you like to define: NUMBER OF STEPS or directly the TIME STEP?"
        print*,"---------------------------------------------------------------------------------------------------"
        print*,"Instructions:"
        print*,"Type 1 to define number of steps."
        print*,"Type 2 to define directly the time step."

        if (debugmode .EQV. .FALSE.) then
          read*, timeChoice
        else
          timeChoice = 2
        end if

        print*,""
        print*,""

        if (timeChoice .EQ. 1) then
          print*,"--------------------------------------------------------------"
          print 201,"What should be the NUMBER OF TIME STEPS?"
          print*,"--------------------------------------------------------------"

          timeStep = simTime
          do while (timeStep .GT. (minTime / 10.0_dp))

            if (debugmode .EQV. .FALSE.) then
              read*, num_steps
            else
              num_steps = 100
            end if

            timeStep = simTime / num_steps

            print*
            print 601, "Time step chosen is ", timeStep, " seconds."

            if (timeStep .GT. minTime) then
              print 601,"Time step should be lower than ", (minTime / 10.0_dp), " seconds."
              print *,"Please increase the number of time steps."
              print *,"The minimum number of time steps should be ", minTimeSteps
            end if

            print*
            print*

          end do

        end if

        if (timeChoice .EQ. 2) then
          print*,"--------------------------------------------------------------"
          print 201,"What should be the TIME STEP in seconds:"
          print*,"--------------------------------------------------------------"

          timeStep = minTime
          do while (timeStep .GT. (minTime / 10.0_dp))

            if (debugmode .EQV. .FALSE.) then
              read*, timeStep
            else
              timeStep = 0.0001_dp
            end if

            num_steps = ceiling(simTime / timeStep)

            print*
            print 401, num_steps," steps will be calculated."

            if (timeStep .GT. (minTime / 10.0_dp)) then
              print 601, "Time step should be lower than ", (minTime / 10.0_dp)
              print *,"Please choose a smaller time step."
            end if

            print*
            print*

          end do
        end if

        print*
        print*,"--------------------------------------------------------------------"
        print 201,"What should be the time step SAVING FREQUENCY?"
        print*,"---------------------------------------------------------------------"
        print*,"Instructions:"
        print*,"Type 1 to save each time step."
        print*,"Type 3 to save each third time step."
        print*,"Type n to save each n-th time step."
        print*

        if (debugmode .EQV. .FALSE.) then
          read*, TStepSavingPeriod
        else
          TStepSavingPeriod = 10000
        end if

        print*

        print 301,"Number of time steps to be saved: ", ceiling(num_steps/TStepSavingPeriod)

        print*
        print*
        print*
        print*

      end subroutine readTimeInfo




      subroutine readSolverSel()
      ! Read the choices about solver method, relaxation, use of ADI in TDMA method.

        use CustomDouble
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: meth,relax_m, ADI
        real(dp) :: tol, alpha
        common /solMethod/ tol, alpha, meth, relax_m, ADI

        print*,"========================================================"
        print*,"Select the SOLUTION METHOD:"
        print*,"========================================================"
        print*,"Instructions:"
        print*,"Type 1 to choose Jacobi method"
        print*,"Type 2 to choose Gauss-Seidel method"
        print*,"Type 3 to choose TDMA (tri-diagonal matrix algorithm)"
        print*
        if (debugmode .EQV. .TRUE.) then
          meth = 2
          print*,"meth = ", meth
        else
          read*, meth
        end if

        if (meth .EQ. 1) then
          print*,"-------------------------------------------"
          print*,"Jacobi method was selected as a solver."
          print*,"-------------------------------------------"
        end if
        if (meth .EQ. 2) then
          print*,"-------------------------------------------------"
          print*,"Gauss-Seidel method was selected as a solver."
          print*,"-------------------------------------------------"
        end if
        if (meth .EQ. 3) then
          print*,"----------------------------------------------------------------"
          print*,"TDMA (tri-diagonal matrix algorithm) was selected as a solver."
          print*,"---------------------------------------------------------------"
        end if


        if (meth .EQ. 3) then
          print*,"--------------------------------------------------------------------"
          print*,"Do you want to use ADI (alternating direction implicit) method?"
          print*,"--------------------------------------------------------------------"
          print*,"Instructions:"
          print*,"Type 1 for Yes"
          print*,"Type 2 for No"
          print*
          if (debugmode .EQV. .TRUE.) then
            ADI = 2
            print*,"ADI = ", ADI
          else
            read*, ADI
          end if

        end if

        print*
        print*
        print*,"---------------------------------------------"
        print*,"What TOLERANCE should be considered?"
        print*,"---------------------------------------------"
        print*,"Range: 1.0D-1 - 1.0D-8"
        if (debugmode .EQV. .TRUE.) then
          tol = 0.1_dp
          print*,"tol = ", tol
        else
          read*, tol

          if (tol .GT. 1.0D-1) then
            tol = 1.0D-1
          elseif (tol .LT. 1.0E-8_dp) then
            tol = 1.0E-8_dp
          end if
        end if

        print*
        print*
        print*,"---------------------------------------------"
        print*,"Do you want to use a RELAXATION FACTOR?"
        print*,"---------------------------------------------"
        print*,"Instructions:"
        print*,"Type 1 for Yes"
        print*,"Type 2 for No"


        if (debugmode .EQV. .TRUE.) then
          alpha = 1.0_dp

          if (alpha .EQ. 1.0_dp) then
            print*,"---------------------------------------------"
            print*,"No relaxation factor will be used"
            print*,"---------------------------------------------"
          end if
          print*,"alpha = ", alpha
        else
          read*, relax_m
          print*
          if (relax_m .EQ. 1) then
            print*,"---------------------------------------------------------------------"
            print*,"What RELAXATION FACTOR do you want to use?"
            print*,"Range: 0.01 - 2.00"
            print*,"Be cautious: the chosen value could cause divergence of solution."
            print*,"---------------------------------------------------------------------"
            print*
            read*,alpha

            if ((alpha .LT. 0.01) .OR. (alpha .GT. 2.0)) then
              print*,"You typed a wrong value for alpha. Default value (alpha = 1.0) will be considered, instead."
              alpha = 1.0_dp
            end if

          else
            alpha = 1.0_dp
          end if
        end if
      end subroutine readSolverSel




      subroutine solveStage()
      ! Routine responsible to call corresponding solver and to control the time evolution

        use CustomDouble
        use Solver

        implicit none

        integer :: iter, curTimeStep, totalIter
        common /control/ iter, curTimeStep, totalIter

        real(dp) :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        integer :: Nchange

        integer :: tsShLog
        common /convergence/ tsShLog

        logical :: debugmode
        common /dbgMode/ debugmode

        real(dp) :: mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf
        common /boundCondLen/ mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        real(dp) :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: CumTimeM, CumTimeS1, CumTimeS2, CumTimeS3

        real(dp) :: W1, W2, H3
        common /secCooling/ W1, W2, H3

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

*        real(dp) :: simTime, curTime, num_steps, timeStep
*        integer :: timeChoice, TStepSavingPeriod
*        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        integer :: meth,relax_m, ADI
        real(dp) :: tol, alpha
        common /solMethod/ tol, alpha, meth, relax_m, ADI

        real(dp) :: Delta_t
        logical :: timeStepChanged

        real(dp) :: startTime, endTime, nextEndTime


        character(LEN=50) :: name

160     format(' ', A, I6)
170     format(' ', A, F12.4, A)

        Nchange = 0
        Delta_t = 0

        if (timeChoice .EQ. 1) then
          Delta_t = simTime / num_steps
        elseif (timeChoice .EQ. 2) then
          Delta_t = timeStep
        end if

        CumTimeM = tMould
        CumTimeS1 = CumTimeM + tS1
        CumTimeS2 = CumTimeS1 + tS2
        CumTimeS3 = CumTimeS2 + tS3

        startTime = 0.0_dp
        endTime = CumTimeS3

        curTime = startTime + Delta_t
        timeStepChanged = .FALSE.

        tsShLog = 0

        do while (curTime .LE. endTime)

            print*
            print*
            print*
            print*
            print*,"--------------------------------------------------------------------------------------------------"
            print*
            print*
            print 170,"CALCULATING TEMPERATURES FOR TIME = ", curTime, " seconds."
            print*
            print*
            print*,"--------------------------------------------------------------------------------------------------"

*         Jacobi method
          if (meth .EQ. 1) then
            name = "Jacobi Method"
            call Jacobi(Delta_t,tol,alpha)
            if (debugmode .EQV. .TRUE.) then
*              call print_res2D(name,Tnew,Np,curTime)
            end if
            if ((totalIter - tsShLog) .LT. 50) then
                print 160,"Jacobi method has terminated - Iterations: ", iter
                print*,"#################################################################################################"
            end if
            if (debugmode .EQV. .TRUE.) then
            ! read*
            end if
          end if

*         Gauss-Seidel method
          if (meth .EQ. 2) then
            name = "Gauss-Seidel Method"
            call GaSe2D(Delta_t, tol,alpha)
            if (debugmode .EQV. .TRUE.) then
*              call print_res2D(name,Tnew,Np,curTime)
            end if

            if ((totalIter - tsShLog) .LT. 50) then
                print 160,"Gauss-Seidel method has terminated - Iterations: ", iter
                print*,"#################################################################################################"
            end if
            if (debugmode .EQV. .TRUE.) then
              ! read*
            end if
          end if

*         TDMA method
          if (meth .EQ. 3) then
            if (ADI .EQ. 1) then
              name = "TDMA with ADI"
              call TDMA2dADI(Delta_t,tol,alpha)
              if (debugmode .EQV. .TRUE.) then
*                call print_res2D(name,Tnew,Np,curTime)
              end if

              if ((totalIter - tsShLog) .LT. 50) then
                print 160,"TDMA method with ADI has terminated - Iterations: ", iter
                print*,"#################################################################################################"
              end if
              if (debugmode .EQV. .TRUE.) then
                ! read*
              end if
            elseif (ADI .EQ. 2) then
              name = "TDMA Simple (Non-ADI)"
              call TDMA2D(Delta_t,tol,alpha)
              if (debugmode .EQV. .TRUE.) then
*                call print_res2D(name,Tnew,Np,curTime)
              end if

              if ((totalIter - tsShLog) .LT. 50) then
                print 160,"TDMA method without ADI has terminated - Iterations: ", iter
                print*,"#################################################################################################"
              end if
              if (debugmode .EQV. .TRUE.) then
                ! read*
              end if
            end if
          end if

          curTimeStep = curTimeStep + 1

          nextEndTime = 0.0_dp

          if (curTime .LT. CumTimeM) then
            nextEndTime = CumTimeM
          elseif (curTime .LT. CumTimeS1) then
            nextEndTime = CumTimeS1
          elseif (curTime .LT. CumTimeS2) then
            nextEndTime = CumTimeS2
          elseif (curTime .LT. CumTimeS3) then
            nextEndTime = CumTimeS3
          end if


          if (mod(curtimeStep, TStepSavingPeriod) .EQ. 0) then
            call SaveResults(totalIter)
            call storeTimeSteps(curTime)
          elseif ((nextEndTime - curTime) .LT. 1E-6_dp) then
            call SaveResults(totalIter)
            call storeTimeSteps(curTime)
          end if


          if (timeStepChanged .EQV. .FALSE.) then

            if ((nextEndTime - curTime) .GT. 0.0001_dp) then
              if ((nextEndTime - curTime) .LT. Delta_t) then
                Delta_t = nextEndTime - curTime
                timeStepChanged = .TRUE.
              end if
            end if
          else
            Delta_t = timeStep
          end if

          if ((Nchange .EQ. 0) .AND. (timeStep .GT. 0.02))then
            if (((Tmin - Tliquidus()) .LT. 2.0_dp) .AND. (Tmin .GT. Tsolidus())) then
              timeStep = timeStep / 100.0_dp
              print*,timeStep
              TStepSavingPeriod = TStepSavingPeriod * 100
              print*,TStepSavingPeriod
              Delta_t = timeStep
              print*,Delta_t
              Nchange = 1
            end if
          end if

          if ((Nchange .EQ. 1) .AND. (timeStep .GT. 0.02)) then
            if (((Tmin - Tsolidus()) .LT. 2.0_dp) .AND. (Tmin .GT. Tsolidus())) then
              timeStep = timeStep * 100.0_dp
              print*,timeStep
              TStepSavingPeriod = TStepSavingPeriod / 100
              print*,TStepSavingPeriod
              Delta_t = timeStep
              print*,Delta_t
              Nchange = 2
            end if
          end if

          curTime = curTime + Delta_t

        end do

      end subroutine solveStage




      subroutine SaveResults(totalIter)
        use CustomDouble
        use Properties

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        integer :: totalIter

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        real(dp) :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        integer :: I, J

*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        real(dp) :: mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf
        common /boundCondLen/ mouldLen, SecColLen1, SecColLen2, SecColLen3, CasSpd, Tf

        real(dp) :: W1, W2, H3
        common /secCooling/ W1, W2, H3

        real(dp) :: tMould, tS1, tS2, tS3, CasTime
        common/timingCast/ tMould, tS1, tS2, tS3, CasTime

        real(dp) :: CumTimeM, CumTimeS1, CumTimeS2, CumTimeS3

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        real(dp) :: Tmax, Tmin, Tavg,ToldMax,ToldMin,ToldAvg
        integer :: Imin,Jmin
        common /stats/ Tmax,Tmin,Tavg, ToldMax,ToldMin,ToldAvg,Imin,Jmin

        real(dp) :: TmaxCh, TminCh, TavgCh

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        real(dp) :: liqFrac

        integer :: meth,relax_m, ADI
        real(dp) :: tol, alpha
        common /solMethod/ tol, alpha, meth, relax_m, ADI

        real(dp) :: Cp, k

        character(LEN = 100) :: filename

        CumTimeM = tMould
        CumTimeS1 = CumTimeM + tS1
        CumTimeS2 = CumTimeS1 + tS2
        CumTimeS3 = CumTimeS2 + tS3

        write(filename,"(A7,F8.2,A4)") "RESULTS", curTime, ".DAT"
        filename = trim(filename)

        print*,"Saving results for this time step in a text file..."

*       Writing results in RESULT.DAT file, stored in the same folder as the program
100     format(' ', 2I5,2F10.4,1F12.4,1F19.4,1F17.4,1F17.4)
200     format(' ', 2A5,2A10,1A12,1A19,1A17,1A17)
300     format(' ', A, F12.4)
400     format(' ', A, I12)
500     format(' ', A, I3, A, I3, A)
600     format(' ', A, F12.4, A, F12.4, A)
800     format(' ', A, F12.4, A)
        open(10,file=filename,status='UNKNOWN')

        write(10,*)"===================================================================================="
        write(10,*)"              Steady state solution for a 2D thermal profile                        "
        write(10,*)"                 in cross section of a slab or billet in                            "
        write(10,*)"                            Continuous Casting                                      "
        write(10,*)"===================================================================================="
        write(10,*)"                                                                                    "
        write(10,*)"  Heat Transfer is calculated considering only heat conduction                      "
        write(10,*)" and convection in liquid metal is taken into consideration by                      "
        write(10,*)" artificially increasing thermal conductivity                                       "
        write(10,*)"                                                                                    "
        write(10,*)" Heat exchanged due to solidification is modelled by                                "
        write(10,*)" equivalent heat capacity method.                                                   "
        write(10,*)"                                                                                    "
        write(10,*)" Boundary Conditions are:                                                           "
        write(10,*)" 1) In the mold, a specified heat flow as a function of time                        "
        write(10,*)" 2) In secondary cooling, a specified heat transfer coefficient                     "
        write(10,*)"    and also radiactive heat transfer.                                              "
        write(10,*)"                                                                                    "
        write(10,*)" Velocity of strand is constant and equal to v sub z                                "
        write(10,*)"                                                                                    "
        write(10,*)"===================================================================================="
        write(10,*)"Author: Dickson Alves de Souza                                                      "
        write(10,*)"                                                                                    "
        write(10,*)"Based on lectures by professor Roberto Parreiras Tavares                            "
        write(10,*)"                                                                                    "
        write(10,*)"and book Numerical Heat Transfer and Fluid Flow                                     "
        write(10,*)"by Suhas V. Patankar (1980)                                                         "
        write(10,*)"                                                                                    "
        write(10,*)"Federal University of Minas Gerais (UFMG)                                           "
        write(10,*)"November 18th, 2017                                                                 "
        write(10,*)"===================================================================================="
        write(10,*)"                             Geometry:                                              "
        write(10,*)"                                                                                    "
        write(10,300)"Lenght in X direction (meters):             ", Lth(1)
        write(10,300)"Lenght in Y direction (meters):             ", Lth(2)
        write(10,400)"Nodes in X direction:                       ", N(1)
        write(10,400)"Nodes in Y direction:                       ", N(2)
        write(10,*)"                                                                                    "
        write(10,300)"Mould Length (meters):                      ", mouldLen
        write(10,600)"Time at mould: ", tMould,"  - Cummulative time: ", CumTimeM," seconds"
        write(10,*)"                                                                                    "
        write(10,300)"Secondary Cooling - Zone 1 - Length (meters)", SecColLen1
        write(10,600)"Time at mould: ", tS1," seconds - Cummulative time: ", CumTimeS1," seconds"
        write(10,*)"                                                                                    "
        write(10,300)"Secondary Cooling - Zone 2 - Length (meters)", SecColLen2
        write(10,600)"Time at mould: ", tS2," seconds - Cummulative time: ", CumTimeS2," seconds"
        write(10,*)"                                                                                    "
        write(10,300)"Secondary Cooling - Zone 3 - Length (meters)", SecColLen3
        write(10,600)"Time at mould: ", tS3," seconds - Cummulative time: ", CumTimeS3," seconds"
        write(10,*)"                                                                                    "
        write(10,300)"Casting speed: ", CasSpd
        write(10,*)"===================================================================================="
        write(10,*)"                             Boundary Conditions:                                   "
        write(10,*)"                                                                                    "
        write(10,300)"Secondary Cooling - Zone 1 - Specific Water Flow Rate (liter/sq-m.s)", W1
        write(10,300)"Secondary Cooling - Zone 2 - Specific Water Flow Rate (liter/sq-m.s)", W2
        write(10,300)"Secondary Cooling - Zone 3 - Heat Transfer Coefficient (W/sq-m.s)", H3
        write(10,*)"===================================================================================="


        write(10,*)" "
        write(10,*)" "
        write(10,*)" "
        write(10,*)"===================================================================================="
        write(10,*)"                             Solver:                                      "
        if (meth .EQ. 1) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Jacobi method"
        elseif (meth .EQ. 2) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Gauss-Seidel method"
        elseif (meth .EQ. 3) then

          write(10,*) "SOLUTION of LINEAR SYSTEM: TDMA (tri-diagonal matrix algorithm)"

          if (ADI .EQ. 1) then
            write(10,*) "     ADI approach employed in solution."
          elseif (ADI .EQ. 2) then
            write(10,*) "     ADI approach NOT employed in solution."
          end if

        end if

        if (relax_m .EQ. 1) then
          write(10,*)"    Relaxation factor: ", alpha
        elseif (relax_m .EQ. 2) then
          write(10,*)"    No relaxation applied to solution."
        end if

        write(10,300)"Tolerance: ", tol
        write(10,*)""
        if (mode .EQ. 1) then
          write(10,*)"Lever rule was used to calculate liquid fraction."
        elseif (mode .EQ. 2) then
          write(10,*)"A linear relationship was used to calculate liquid fraction."
        end if
        write(10,*)"===================================================================================="

        write(10,*)" "
        write(10,*)" "
        write(10,*)" "

        TmaxCh = Tmax - ToldMax
        TminCh = Tmin - ToldMin
        TavgCh = Tavg - ToldAvg

        write(10,*)"==========================================================================================================="
        write(10,*)"                   Calculation Results                                              "
        write(10,800)"Current time:                                  ", curTime," seconds               "
        write(10,800)"Time Step:                                     ", timeStep," seconds               "
        write(10,400)"Iterations:                                    ", totalIter
        write(10,600)"Maximum Temperature :                          ", Tmax, " K - Change from last time step: ", TmaxCh ," K"
        write(10,600)"Minimum Temperature :                          ", Tmin, " K - Change from last time step: ", TminCh ," K"
        write(10,600)"Average(algebraic, non-weigthed) Temperature : ", Tavg, " K - Change from last time step: ", TavgCh ," K"
        write(10,800)"Liquidus Temperature:                          ", Tliquidus(), " K"
        write(10,800)"Solidus Temperature:                           ", Tsolidus(), " K"
        write(10,*)"==========================================================================================================="

        write(10,200)"I","J","X(m)","Y(m)","T(K)","Cp_Equiv(J/kg.K)", "K(W/sq-m.K)","Liquid_Fraction"

        write(10,*)"==========================================================================================================="

        do I = 1, N(1), 1
          do J = 1, N(2), 1
            Cp = Cp_Eq(Tnew(I,J))
            k = k_Mushy(Tnew(I,J))
            liqFrac = f_Liq(Tnew(I,J))
            write(10,100)I,J,Np(I,J,1),Np(I,J,2),Tnew(I,J),Cp,k,liqFrac
          end do
        end do

        endfile 10
        close(10,status='KEEP')
        print*,"Results were successfully saved. Next time step..."

      end subroutine SaveResults




      subroutine storeTimeSteps(curTime)
        use CustomDouble
        use Solver

        implicit none

        real(dp) :: curTime

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
        real(dp) :: TmRs(1:81,0:50000,1:7)
        common /transient/ TmRs

        ! SteelChem(1): Carbon
        ! SteelChem(2): Manganese
        ! SteelChem(3): Silicon
        ! SteelChem(4): Phosphorus
        ! SteelChem(5): Sulfur
        ! SteelChem(6): Aluminium
        real(dp) :: SteelChem(1:6)
        common /steelprop/ SteelChem

        ! Define the formula for liquid fraction calculation: lever rule or linear variation in mushy zone
        ! Mode = 1 > Lever rule
        ! Mode = 2 > Linear variation of liquid fraction
        integer :: mode
        common /LiqFracCalc/ mode

        real(dp) :: liqFrac, ShTh

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        integer :: I, J, Npos

        integer :: maxX, maxY, Xi(9), Yi(9), tsSavPos
        common /storeLimits/ maxX, maxY, Xi, Yi, tsSavPos
*       Np: array containing node positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) in X direction, from 1 to N(2) in Y direction.
*         The first two indices locate the node indexes  - I for x axis and J for y axis
*         The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)

*       Ni: array containing interface positions in X and Y axis (to be calculated)
*       Index ranges from 1 to N(1) + 1 in X direction, from 1 to N(2) + 1 in Y direction.
*           The first two indices locate the node indexes  - I for x axis and J for y axis
*           The third index indicates if its a distance in X direction (index = 1) or in Y direction (index = 2)
        real(dp) :: Np(1:999,1:999,1:2), Ni(1:999,1:999,1:2)
        common /coordinate/ Np, Ni

        real(dp) :: Tnew(0:999,0:999), Told(0:999,0:999)
        common /tempRes/ Tnew, Told

        Npos = 1

        do I = 1, maxX, 1
          do J = 1, maxY, 1

            ShTh = ShellThick(Xi(I),1)
            ShTh = ShellThick(Yi(J),2)

            TmRs(Npos,tsSavPos,1) = curTime
            TmRs(Npos,tsSavPos,2) = Np(Xi(I),Yi(J),1)
            TmRs(Npos,tsSavPos,3) = Np(Xi(I),Yi(J),2)
            TmRs(Npos,tsSavPos,4) = Tnew(Xi(I),Yi(J))
            TmRs(Npos,tsSavPos,5) = Cp_Eq(Tnew(Xi(I),Yi(J)))
            TmRs(Npos,tsSavPos,6) = k_Mushy(Tnew(Xi(I),Yi(J)))

            liqFrac = 0.0_dp

            if (mode .EQ. 1) then
              liqFrac = f_Liq_Lever(Tnew(I,J))
            elseif (mode .EQ. 2) then
              liqFrac = f_Linear(Tnew(I,J))
            end if

            TmRs(Npos,tsSavPos,7) = liqFrac

            Npos = Npos + 1
          end do
        end do

        tsSavPos = tsSavPos + 1

      end subroutine storeTimeSteps




      subroutine defineLimits()
        ! Define the nodes which will be saved for each step in a single file
        ! Up to 81 nodes can be saved to be inspected later regarding time evolution

        use CustomDouble
        implicit none

        integer :: rpt, I, J
        integer :: maxX, maxY, Xi(9), Yi(9), tsSavPos
        common /storeLimits/ maxX, maxY, Xi, Yi, tsSavPos

        integer :: N(1:2)
        common /gridsize/ N

        Xi(1) = 1
        Xi(2) = ceiling(1.0_dp * float(N(1)) / 8.0_dp)
        Xi(3) = ceiling(2.0_dp * float(N(1)) / 8.0_dp)
        Xi(4) = ceiling(3.0_dp * float(N(1)) / 8.0_dp)
        Xi(5) = ceiling(4.0_dp * float(N(1)) / 8.0_dp)
        Xi(6) = ceiling(5.0_dp * float(N(1)) / 8.0_dp)
        Xi(7) = ceiling(6.0_dp * float(N(1)) / 8.0_dp)
        Xi(8) = ceiling(7.0_dp * float(N(1)) / 8.0_dp)
        Xi(9) = N(1)

        maxX = 9
        rpt = 9

        do while (rpt .GT. 0)
          do I = 1, maxX - 2, 1

            if (I .EQ. 1) then
              rpt = 0
            end if

            if (Xi(I+1) .EQ. Xi(I)) then
              Xi(I+1) = Xi(I+2)
              rpt = rpt + 1
            end if

          end do

          if (rpt .GT. 0) then
            maxX = maxX - 1
          end if

        end do

        if (Xi(maxX - 1) .EQ. Xi(maxX)) then
          maxX = maxX - 1
        end if

        Yi(1) = 1
        Yi(2) = ceiling(1.0_dp * float(N(2)) / 8.0_dp)
        Yi(3) = ceiling(2.0_dp * float(N(2)) / 8.0_dp)
        Yi(4) = ceiling(3.0_dp * float(N(2)) / 8.0_dp)
        Yi(5) = ceiling(4.0_dp * float(N(2)) / 8.0_dp)
        Yi(6) = ceiling(5.0_dp * float(N(2)) / 8.0_dp)
        Yi(7) = ceiling(6.0_dp * float(N(2)) / 8.0_dp)
        Yi(8) = ceiling(7.0_dp * float(N(2)) / 8.0_dp)
        Yi(9) = N(2)

        maxY = 9
        rpt = 9

        do while (rpt .GT. 0)
          do J = 1, maxY - 2, 1

            if (J .EQ. 1) then
              rpt = 0
            end if

            if (Yi(J+1) .EQ. Yi(J)) then
              Yi(J+1) = Yi(J+2)
              rpt = rpt + 1
            end if

          end do

          if (rpt .GT. 0) then
            maxY = maxY - 1
          end if

        end do

        if (Yi(maxY - 1) .EQ. Yi(maxY)) then
          maxY = maxY - 1
        end if

      end subroutine defineLimits




      subroutine saveTimeSteps()
        use CustomDouble
        implicit none

        integer :: I, J, K, Npos

*       First index: 9 nodes x 9 nodes = 81 selected positions to be saved
*       Second index: up to 50000 time steps
*       Third index:
*          1: curTime
*          2: X position
*          3: Y position
*          4: Temperature
*          5: Specific equivalent heat
*          6: Thermal conductivity
*          7: Liquid fraction
        real(dp) :: TmRs(1:81,0:50000,1:7)
        common /transient/ TmRs

        integer :: maxX, maxY, Xi(9), Yi(9), tsSavPos
        common /storeLimits/ maxX, maxY, Xi, Yi, tsSavPos

*       Lth: Lenght along each direction
*          1,2 indices refer to x and y respectively
        real(dp) :: Lth(1:2)
        common /domainsize/ Lth

*       N: number of nodes in each direction
*          1,2 indices refer to x and y respectively
        integer :: N(1:2)
        common /gridsize/ N

        real(dp) :: simTime, curTime, num_steps, timeStep
        integer :: timeChoice, TStepSavingPeriod
        common /time/ simTime, curTime, num_steps, timeStep, TStepSavingPeriod, timeChoice

        character(LEN = 100) :: filename


120     format(' ', 6A12, 1A15)
220     format(' ', 7F12.4)
320     format(' ', A,I3)
420     format(' ', A9,I3,A3,I3,A1,F6.2,A6,A4)
520     format(' ', A, F12.4)

        write(filename,420) "TempEvol_", N(1), "_X_", N(2), "_", timeStep,"_s_CNT",".DAT"

        filename = trim(filename)

        open(20,file=filename,status='UNKNOWN')

        Npos = 1
        write(20,*)"========================================================================"
        write(20,*)"Mesh information:"
        write(20,320)"Nx = ", N(1)
        write(20,320)"Ny = ", N(2)
        write(20,*)" "
        write(20,*)"Slab / Billet dimensions:"
        write(20,520)"X(m) = ", Lth(1)
        write(20,520)"Y(m) = ", Lth(2)

        write(20,*)"                                                                        "
        write(20,*)"                                                                        "
        write(20,*)"========================================================================"

        do I = 1, maxX, 1
          do J = 1, maxY, 1
            write (20,320)"Npos = ", Npos
            write(20,120)"Time(s)","X(m)","Y(m)","Temp(K)","Cp(J/kg.K)","K(W/sq-m.K)", "Liquid fraction"

            do K = 0, tsSavPos, 1
              write(20,220)TmRs(Npos,K,1),TmRs(Npos,K,2),TmRs(Npos,K,3),TmRs(Npos,K,4),TmRs(Npos,K,5),TmRs(Npos,K,6), TmRs(Npos,K,7)
            end do

            Npos = Npos + 1

            write(20,*)"                                                                        "
            write(20,*)"                                                                        "
            write(20,*)"========================================================================"

          end do
        end do

        endfile 20
        close(20,status='KEEP')

      end subroutine saveTimeSteps
