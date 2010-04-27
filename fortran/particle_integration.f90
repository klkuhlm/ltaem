! $Id: particle_integration.f90,v 1.3 2008/09/07 15:57:10 kris Exp kris $
module particle_integrate
  implicit none
  private
  public :: rungekuttamerson, rungekutta, fwdEuler

contains

  subroutine rungekuttamerson(p,tee,coeff,lo,num)
    ! results are saved to Presult, which is shared via module with main program
    use constants, only : DP
    use inverse_laplace_transform, only : invlap => deHoog_invlap
    use element_specs, only : INValpha,INVtol,INVm,PARti, PARx, PARy, PARdt, PARtf, PARtol, PARmin, PARmaxStep
    use calc_shared_data, only  : Presult
    use calc_routines, only : velXCalc, velYCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: p(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    complex(DP), intent(in) :: coeff(1:,lo:,0:,1:)  !! coefficients of elements
    integer, intent(in) :: num  !! which particle is being integrated
    integer :: count, ilogc
    logical :: partEnd
    real(DP) :: xVelInit, yVelInit, xVelTrap, yVelTrap, xVelAB3, yVelAB3, xVelAB2, yVelAB2
    real(DP) :: xVelSF, yVelSF
    real(DP) :: FwdEulerX, FwdEulerY, TrapX, TrapY, halfABX, halfABY, fullABX, fullABY, SimpX, SimpY
    real(DP) :: error, pt, px, py, dt, tf, tmp, L
    real(DP), parameter :: safety = 0.9

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !        not using porosity correctly ?
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    ! read this in eventually: for now just set it here
    PARmin = 1.1D-8

    pt = PARti(num)
    px = PARx(num)
    py = PARy(num)
    dt = PARdt
    tf = PARtf(num)

    Presult(0,1,num) = pt
    Presult(0,2,num) = px
    Presult(0,3,num) = py
    count = 1

    partEnd = .false.
    write(*,'(A)')    '*******************************' 
    write(*,'(A,I3)') 'rkm integration, particle ',num
    write(*,'(A)')    '*******************************' 

    ! Runge-Kutta-Merson 4th-order adaptive integration scheme

    !! debugging
    write(77,'(A)') '#        pt            px          py          dt       ilogc     xV_0    &
         &    yV_0        x_1         y_1      ilogc     xV_1        yV_1        x_2         y_2 &
         &           xV_2        yV_2       x_3          y_3      ilogc     xV_3        yV_3   &
         &     x_4         y_4      ilogc     xV_4        yV_4        x_5        y_5           &
         &error        L '
    !! debugging

    rkm: do 
       if (partEnd .or. pt + dt >= tf) exit rkm

       !! debugging
       write(77,'(ES18.10,3(ES11.3E2,1X))',advance='no') pt,px,py,dt
       !!

       ! forward Euler 1/3-step  (predictor)
       ilogc = ceiling(log10(pt))

       xVelInit = invlap(INValpha,INVtol,pt,tee(ilogc), &
            & velXCalc(p(:,ilogc),px,py,coeff(:,ilogc,:,:)),INVm)
       yVelInit = invlap(INValpha,INVtol,pt,tee(ilogc), &
            & velYCalc(p(:,ilogc),px,py),INVm)

       FwdEulerX = px + dt/3.0_DP*xVelInit
       FwdEulerY = py + dt/3.0_DP*yVelInit

       !! debugging
       write(77,'(A,I3,1X,4(ES11.3E2,1X))',advance='no') '   ',&
            & ilogc,xVelInit,yvelinit,fwdeulerx,fwdeulery
       !!

       ! trapazoid rule 1/3-step (corrector)
       ilogc = ceiling(log10(pt + dt/3.0_DP))

       xVelTrap = invlap(INValpha,INVtol,pt + dt/3.0_DP,tee(ilogc), &
            & velXCalc(p(:,ilogc),FwdEulerX,FwdEulerY,coeff(:,ilogc,:,:)),INVm)
       yVelTrap = invlap(INValpha,INVtol,pt + dt/3.0_DP,tee(ilogc), &
            & velYCalc(p(:,ilogc),FwdEulerX,FwdEulerY),INVm)

       TrapX = px + dt/6.0_DP*(xVelInit + xVelTrap)
       TrapY = py + dt/6.0_DP*(yVelInit + yVelTrap)

       !! debugging
       write(77,'(A,I3,1X,4(ES11.3E2,1X))',advance='no') '   ',&
            & ilogc,xVelTrap,yveltrap,trapx,trapy
       !!

       ! Adams-Bashforth 1/2-step predictor 
       xVelAB3 = invlap(INValpha,INVtol,pt+dt/3.0_DP,tee(ilogc), &
            & velXCalc(p(:,ilogc),TrapX,TrapY,coeff(:,ilogc,:,:)),INVm)
       yVelAB3 = invlap(INValpha,INVtol,pt+dt/3.0_DP,tee(ilogc), &
            & velYCalc(p(:,ilogc),TrapX,TrapY),INVm)

       halfABX = px + dt/8.0_DP*(xVelInit + 3.0_DP*xVelAB3)
       halfABY = py + dt/8.0_DP*(yVelInit + 3.0_DP*yVelAB3)

       !! debugging
       write(77,'(A,4(ES11.3E2,1X))',advance='no') '   ',xVelAB3,yvelAB3,halfABx,halfABy
       !!

       ! full step Adams-Bashforth predictor
       ilogc = ceiling(log10(pt + dt/2.0_DP))

       xVelAB2 = invlap(INValpha,INVtol,pt+dt/2.0_DP,tee(ilogc), &
            & velXCalc(p(:,ilogc),halfABX,halfABY,coeff(:,ilogc,:,:)),INVm)
       yVelAB2 = invlap(INValpha,INVtol,pt+dt/2.0_DP,tee(ilogc), &
            & velYCalc(p(:,ilogc),halfABX,halfABY),INVm)

       fullABX = px + dt/2.0_DP*(xVelInit - 3.0_DP*xVelAB3 + 4.0_DP*xVelAB2)
       fullABY = py + dt/2.0_DP*(yVelInit - 3.0_DP*yVelAB3 + 4.0_DP*yVelAB2)

       !! debugging
       write(77,'(A,I3,1X,4(ES11.3E2,1X))',advance='no') '   ',&
            & ilogc,xVelAB2,yvelAB2,fullABx,fullABy
       !!

       ! full step Simpson's rule corrector
       ilogc = ceiling(log10(pt + dt))

       xVelSF = invlap(INValpha,INVtol,pt+dt,tee(ilogc), &
            & velXCalc(p(:,ilogc),halfABX,halfABY,coeff(:,ilogc,:,:)),INVm)
       yVelSF = invlap(INValpha,INVtol,pt+dt,tee(ilogc), &
            & velYCalc(p(:,ilogc),halfABX,halfABY),INVm)

       SimpX = px + dt/6.0_DP*(xVelInit + 4.0_DP*xvelAB2 + xvelSF)
       SimpY = py + dt/6.0_DP*(yVelInit + 4.0_DP*yvelAB2 + yvelSF)

       !! debugging
       write(77,'(A,I3,1X,4(ES11.3E2,1X))',advance='no') '   ',&
            & ilogc,xVelSF,yvelSF,Simpx,Simpy
       !!

       ! relative error
       error = max(abs((SimpX - FullABX)/SimpX), abs((SimpY - FullABY)/SimpY))

       ! magnitude of step taken
       L = sqrt((SimpX - px)**2 + (SimpY - py)**2)

       write(77,'(A,2(ES11.3E2,1X))') '  ',error,L

       ! only advance to next step if error level is acceptable
       ! _and_ resulting step is less than prescribed limit 
       ! _or_ we have gotten to the minimum step size (proceed anyways)
       if ((error < PARtol .and. L <= PARmaxStep) .or. dt <= PARmin) then
          pt = pt + dt
          px = SimpX
          py = SimpY
          Presult(count,1,num) = pt
          Presult(count,2,num) = px
          Presult(count,3,num) = py
          Presult(count-1,4,num) = xVelSF ! velocity associated with previous step
          Presult(count-1,5,num) = yVelSF
          count = count + 1;

          partEnd = sinkCheck(px,py,num)
          
          if (partEnd) then
             write(*,'(A,I3,A,ES12.6E2)') 'particle ',num,' entered a sink at t=',pt
          else 
!!$             write(*,'(A,ES12.6E2)') 'advance to t=',pt
             if (count == ubound(Presult,dim=1)) then
                ! if out of space in results array double the size
                call reallocate(ubound(Presult,dim=1)*2)
             end if
          end if
       else
!!$          write(*,'(2(A,ES12.6E2))') 're-do; error=',error, ' step=',L
       end if

       ! new step size based on error estimate: 
       ! Numerical Recipes, Press et al 1992, eqn 16.2.10, p 712
       if (.not. partEnd .and. pt <= PARtf(num)) then

          if (L > PARmaxStep) then
             ! cut time step to minimize distance
             ! integrated in one step

             tmp = dt/3.0

!             print *, 'reduce1:',dt,'->',tmp

          elseif (PARtol >= error) then
             ! enlarge dt to speed up integration when possible
             tmp = safety*dt*(abs(PARtol/error))**(0.2)

!             print *, 'enlarge:',dt,'->',tmp
          else
             ! reduce dt to increase accurace when needed
             tmp = safety*dt*(abs(PARtol/error))**(0.25)
!             print *, 'reduce2:',dt,'->',tmp
          end if

          dt = tmp

          if (dt < PARmin) then
             write(*,'(3(A,ES12.6E2))') 'min step; dt=',dt,' error=',error, ' pt=',pt
             dt = PARmin
          end if
          
       end if

    end do rkm

  end subroutine rungekuttamerson

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine rungekutta(p,tee,coeff,lo,num)
    use constants, only : DP, RTWO
    use inverse_laplace_transform, only : invlap => deHoog_invlap 
    use element_specs, only : INValpha,INVtol,INVm,PARti, PARx, PARy, PARdt, PARtf
    use calc_shared_data,only  : Presult
    use calc_routines, only : velXCalc, velYCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: coeff(1:,lo:,0:,1:)
    complex(DP), intent(in) :: p(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    integer, intent(in) :: num
    integer :: i, numdt, ilogc, np
    logical :: partEnd
    real(DP) :: xVelInit, yVelInit, xVelBkwdEuler, yVelBkwdEuler, xVelMidpt, yVelMidpt, xVelSimp, yVelSimp
    real(DP) :: FwdEulerX, FwdEulerY, BkwdEulerX, BkwdEulerY, MidptX, MidptY, SimpX, SimpY
    real(DP) :: pt, px, py, dt, tf
    complex(DP), dimension(1:size(p,1)) :: xp,yp

    np = size(p,1)
    pt = PARti(num)
    px = PARx(num)
    py = PARy(num)
    dt = PARdt
    tf = PARtf(num)

    Presult(0,1,num) = pt
    Presult(0,2,num) = px
    Presult(0,3,num) = py

    write(*,'(A)') '  '
    write(*,'(A)')    '*******************************' 
    write(*,'(A,I3)') 'rk integration, particle ',num
    write(*,'(A)')    '*******************************' 

    ! Runge-Kutta 4th order integration scheme (non-adaptive)
    numdt = ceiling((tf - pt)/dt)
    write(*,'(A,ES11.5,A,I6)') 'step size=',dt,' number steps needed=',numdt

    rk: do i = 1,numdt   

       if (mod(i,100) == 0) write(*,'(A,ES12.6E2)') 't=',pt

       ! see if particle will reach end this step
       if (pt + dt >= tf) then
          write(*,'(A,I3,A,ES12.6E2)') &
               & 'particle',num,' reached final time:',tf
          exit rk
       end if

       ! forward Euler 1/2-step  (predictor)
       ilogc = ceiling(log10(pt))

       xp(1:np) = velXCalc(p(:,ilogc),px,py,coeff(:,ilogc,:,:))
       yp(1:np) = velYCalc(p(:,ilogc),px,py)
       xVelInit = invlap(INValpha,INVtol,pt,tee(ilogc), xp(:) ,INVm)
       yVelInit = invlap(INValpha,INVtol,pt,tee(ilogc), yp(:) ,INVm)

       FwdEulerX = px + dt/RTWO*xVelInit
       FwdEulerY = py + dt/RTWO*yVelInit

       ! backward Euler 1/2-step (corrector)
       ilogc = ceiling(log10(pt + dt/RTWO))

       xp(1:np) = velXCalc(p(:,ilogc),FwdEulerX,FwdEulerY,coeff(:,ilogc,:,:))
       yp(1:np) = velYCalc(p(:,ilogc),FwdEulerX,FwdEulerY)
       xvelBkwdEuler = invlap(INValpha,INVtol,pt + dt/RTWO,tee(ilogc), xp(:) ,INVm)
       yvelBkwdEuler = invlap(INValpha,INVtol,pt + dt/RTWO,tee(ilogc), yp(:) ,INVm)

       BkwdEulerX = px + dt/RTWO*xvelBkwdEuler
       BkwdEulerY = py + dt/RTWO*yvelBkwdEuler

       ! midpoint rule full-step predictor
       ilogc = ceiling(log10(pt + dt))
      
       xp(1:np) = velXCalc(p(:,ilogc),BkwdEulerX,BkwdEulerY,coeff(:,ilogc,:,:))
       yp(1:np) = velYCalc(p(:,ilogc),BkwdEulerX,BkwdEulerY)
       xvelMidpt = invlap(INValpha,INVtol,pt + dt,tee(ilogc), xp(:) ,INVm)
       yvelMidpt = invlap(INValpha,INVtol,pt + dt,tee(ilogc), yp(:) ,INVm)

       MidptX = px + dt*xvelMidpt
       MidptY = py + dt*yvelMidpt

       ! Simpson's rule full-step corrector

       xp(1:np) = velXCalc(p(:,ilogc),MidptX,MidptY,coeff(:,ilogc,:,:))
       yp(1:np) = velYCalc(p(:,ilogc),MidptX,MidptY)
       xvelSimp = invlap(INValpha,INVtol,pt + dt,tee(ilogc), xp(:) ,INVm)
       yvelSimp = invlap(INValpha,INVtol,pt + dt,tee(ilogc), yp(:) ,INVm)

       SimpX = px + dt/6.0_DP*(xVelInit + RTWO*xvelBkwdEuler + RTWO*xvelMidpt + xvelSimp)
       SimpY = py + dt/6.0_DP*(yVelInit + RTWO*yvelBkwdEuler + RTWO*yvelMidpt + yvelSimp)

       pt = pt + dt
       px = SimpX
       py = SimpY
       
       Presult(i,1,num) = pt
       Presult(i,2,num) = px
       Presult(i,3,num) = py
       Presult(i-1,4,num) = xvelSimp
       Presult(i-1,5,num) = yvelSimp
       
       partEnd = sinkCheck(px,py,num)
       if (partEnd) then
          write(*,'(A,I3,A,ES12.6E2)') 'particle ',num,' entered a sink at t=',pt
          exit rk             
       end if      
    end do rk
  end subroutine rungekutta

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine fwdEuler(p,tee,coeff,lo,num)
    use constants, only : DP
    use inverse_laplace_transform, only : invlap => deHoog_invlap 
    use element_specs, only : INValpha,INVtol,INVm,PARti, PARx, PARy, PARdt, PARtf
    use calc_shared_data, only : Presult
    use calc_routines, only : velXCalc, velYCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: p(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    complex(DP), intent(in) :: coeff(1:,lo:,0:,1:)  ! coefficients for LT-AEM elements
    integer, intent(in) :: num  ! number of particles

    complex(DP), dimension(1:size(p,1)) :: xp,yp
    integer :: i, numdt, ilogc
    logical :: partEnd
    real(DP) :: pt, px, py, dt, xvel, yvel, tf

    pt = PARti(num)
    px = PARx(num)
    py = PARy(num)
    dt = PARdt
    tf = PARtf(num)

    Presult(0,1,num) = pt
    Presult(0,2,num) = px
    Presult(0,3,num) = py

    ! 1st order fwd Euler
    numdt = ceiling((tf - pt)/dt)
    
    write(*,'(A)')    '****************************************'
    write(*,'(A,I3)') 'forward Euler integration, particle ', num
    write(*,'(A)')    '****************************************'

    fe: do i = 1,numdt   
       if(mod(i,500) ==0 ) write(*,'(A,ES12.6E2)') 't=',pt

       ! see if particle will reach end this step
       if (pt + dt >= tf) then
          write(*,'(A,I3,A,ES12.6E2)') 'particle',num,' reached specified ending time:',tf
          exit fe
       end if

       ! full step forward Euler
       ilogc = ceiling(log10(pt))

       xp(:) = velXCalc(p(:,ilogc),px,py,coeff(:,ilogc,:,:))
       yp(:) = velYCalc(p(:,ilogc),px,py)

       xVel = invlap(INValpha,INVtol,pt,tee(ilogc),xp(:),INVm)
       yVel = invlap(INValpha,INVtol,pt,tee(ilogc),yp(:),INVm)

       px = px + dt*xvel
       py = py + dt*yvel

       pt = pt + dt

       Presult(i,1,num) = pt
       Presult(i,2,num) = px
       Presult(i,3,num) = py
       Presult(i-1,4,num) = xvel
       Presult(i-1,5,num) = yvel

       partEnd = sinkCheck(px,py,num)
       if (partEnd) then
          write(*,'(A,I3,A,ES12.6E2)') 'particle ',num,' entered a sink at t=',pt
          exit fe
       end if

    end do fe
  end subroutine fwdEuler

  !###########################################################################
  ! re-allocate Presult array to longer/shorter first dimension, saving contents
  ! Presult array is in a module

  subroutine reallocate(newub)
    use constants, only : DP
    use error_handler, only : subError
    use calc_shared_data, only : Presult
    use element_specs, only : PARnum
    implicit none

    integer, intent(in) :: newub
    integer :: oldub, ierr
    real(DP), dimension(0:ubound(Presult,dim=1),5,PARnum) :: tmp
    character(128) :: sn = 'allocate', csn = 'reallocate'

    oldub = ubound(Presult,dim=1)
    tmp(0:oldub,:,:) = Presult(0:oldub,:,:)

    deallocate(Presult)
    allocate(Presult(0:newub,5,PARnum), stat=ierr)
    if (ierr /= 0) call subError(sn,ierr,csn)

    Presult = 0.0_DP
    if (newub > oldub) then
       Presult(0:oldub,:,:) = tmp(0:oldub,:,:)
    else
       Presult(0:newub,:,:) = tmp(0:newub,:,:)
    end if

  end subroutine reallocate

  !###########################################################################
  ! check if particle has moved into a sink

  function sinkCheck(px,py,num) result(partEnd)
    use constants, only : DP
    use element_specs, only : WLx, WLy, WLr, WLnum, CIx, CIy, CIr, CInum, PARInclIn, CIibnd
    implicit none

    real(DP), dimension(WLnum) :: rw
    real(DP), intent(in) :: px,py
    integer, intent(in) :: num  ! index of particle
    logical :: partEnd
    integer :: i
    real(DP) :: ri

    partEnd = .false.

    ! did the particle enter a well?
    forall (i = 1:WLnum)
       rw(i) = sqrt((px - WLx(i))**2 + (py - WLy(i))**2)
    end forall
    if (any(rw <= WLr)) partEnd = .true.

    ! did the particle enter/leave a constant head/flux inclusion?            
    do i = 1,CInum     
       if (CIibnd(i) /= 0) then ! not a matching inclusion
          ri = sqrt((px - CIx(i))**2 + (py - CIy(i))**2)

          ! particle started outside a CH/CF inclusion
          if (.not. PARInclIn(num)) then
             if (ri <= CIr(i)) then ! it is now inside a CH/CF inclusion
                partEnd = .true.
             end if

          ! particle started inside a CH/CF inclusion
          else 
             if (ri >= CIr(i)) then ! it is now outside a CH/CF inclusion
                partEnd = .true. 
             end if
          end if
       end if
    end do

    ! dealing with area-sinks is more complex, will be dealt with later
  end function sinkCheck

end module particle_integrate

