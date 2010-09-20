module particle_integrate
  implicit none
  private
  public :: rungekuttamerson, rungekutta, fwdEuler

contains

  subroutine rungeKuttaMerson(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : invlap => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : velcalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: s(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    type(particle), intent(inout) :: p
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: sol
    type(domain), intent(in) :: dom
    type(element), intent(in) :: bg

    integer :: count, ilogc, los, his, ns
    logical :: partEnd
    complex(DP), dimension(size(s,dim=1),2) :: velp
    real(DP), dimension(2) :: vInit,vTrap,vAB3,vAB2,vSF
    real(DP), dimension(2) :: fwdEuler,Trap,halfAB,fullAB,Simp
    real(DP) :: error, pt, px, py, dt, tf, L
    real(DP), parameter :: SAFETY = 0.9

    !  TODO not using porosity correctly ?

    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    tf = p%tf

    p%r(0,1) = pt
    p%r(0,2) = px
    p%r(0,3) = py
    count = 1

    partEnd = .false.
    write(*,'(A)')    '*******************************' 
    write(*,'(A,I0)') 'rkm integration, particle ',p%id
    write(*,'(A)')    '*******************************' 

    ! Runge-Kutta-Merson 4th-order adaptive integration scheme
    rkm: do 
       if (partEnd .or. pt + dt >= tf) exit rkm

       ! forward Euler 1/3-step  (predictor)
       ilogc = ceiling(log10(pt))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(px,py,DP),s(:,ilogc),los,his,dom,c,e,bg)
       VInit(1:2) = invlap(pt,tee(ilogc),velp(:,1:2),sol%INVLT)

       FwdEuler(1:2) = [px,py] + dt/3.0*VInit(1:2)

       ! trapazoid rule 1/3-step (corrector)
       ilogc = ceiling(log10(pt + dt/3.0))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(FwdEuler(1),FwdEuler(2),DP),&
            & s(:,ilogc),los,his,dom,c,e,bg)
       VTrap(1:2) = invlap(pt + dt/3.0,tee(ilogc),velp(:,1:2),sol%INVLT)

       Trap(1:2) = [px,py] + dt/6.0*(VInit(1:2) + VTrap(1:2))

       ! Adams-Bashforth 1/2-step predictor 
       velp(1:ns,1:2) = velCalc(cmplx(Trap(1),Trap(2),DP),&
            & s(:,ilogc),los,his,dom,c,e,bg)
       VAB3(1:2) = invlap(pt+dt/3.0,tee(ilogc),velp(:,1:2),sol%INVLT)

       halfAB(1:2) = [px,py] + dt/8.0*(VInit(1:2) + 3.0*VAB3(1:2))

       ! full step Adams-Bashforth predictor
       ilogc = ceiling(log10(pt + dt/2.0))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(halfAB(1),halfAB(2),DP), &
            & s(:,ilogc),los,his,dom,c,e,bg)
       VAB2(1:2) = invlap(pt+dt/2.0,tee(ilogc),velp(:,1:2),sol%INVLT)

       fullAB(1:2) = [px,py] + dt/2.0* &
            & (VInit(1:2) - 3.0*VAB3(1:2) + 4.0*VAB2(1:2))

       ! full step Simpson's rule corrector
       ilogc = ceiling(log10(pt + dt))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(halfAB(1),halfAB(2),DP), &
            & s(:,ilogc),los,his,dom,c,e,bg)
       VSF = invlap(pt+dt,tee(ilogc), velp(:,1:2),sol%INVLT)

       Simp(1:2) = [px,py] + dt/6.0*(VInit + 4.0*vAB2 + vSF)

       ! relative error
       error = maxval(abs((Simp - fullAB)/Simp))

       ! magnitude of total step taken
       L = abs(cmplx(Simp(1),Simp(2),DP) - cmplx(px,py,DP))

       ! only advance to next step if error level is acceptable
       ! _and_ resulting step is less than prescribed limit 
       ! _or_ we have gotten to the minimum step size (proceed anyways)
       if ((error < p%tol .and. L <= p%maxL) .or. dt <= p%mindt) then
          pt = pt + dt
          px = Simp(1)
          py = Simp(2)
          p%r(count,1) = pt
          p%r(count,2) = px
          p%r(count,3) = py
          p%r(count-1,4) = VSF(1) ! velocity associated with previous step
          p%r(count-1,5) = VSF(2)
          count = count + 1

          partEnd = sinkCheck(px,py,c,e)
          
          if (partEnd) then
             write(*,'(A,I0,A,ES12.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          else 
             if (count == ubound(p%r,dim=1)) then
                ! if out of space in results array double the size
                call reallocate(p,ubound(p%r,dim=1)*2)
             end if
          end if
       end if

       ! new step size based on error estimate: 
       ! Numerical Recipes, Press et al 1992, eqn 16.2.10, p 712
       if (.not. partEnd .and. pt <= p%tf) then

          if (L > p%maxL) then
             ! cut time step to minimize distance
             ! integrated in one step
             dt = dt/3.0

          elseif (p%tol >= error) then
             ! enlarge dt to speed up integration when possible
             dt = SAFETY*dt*abs(p%tol/error)**0.200

          else
             ! reduce dt to increase accurace when needed
             dt = SAFETY*dt*abs(p%tol/error)**0.250

          end if

          if (dt <= p%mindt) then
             write(*,'(3(A,ES12.6E2))') 'min step; dt=',dt,' error=',error, ' pt=',pt
             dt = p%mindt
          end if
       end if
    end do rkm
    p%numt = count - 1
  end subroutine rungekuttamerson

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine rungekutta(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : invlap => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : velCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: s(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    type(particle), intent(inout) :: p
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: sol
    type(domain), intent(in) :: dom
    type(element), intent(in) :: bg

    integer :: i, numdt, ilogc, ns, los, his
    real(DP) :: pt, px, py, dt, tf
    complex(DP), dimension(1:size(s,dim=1),2) :: velp
    real(DP), dimension(2) :: FwdEuler,BkwdEuler,MidPt,Simp
    real(DP), dimension(2) :: vInit,vBkwdEuler,vMidpt,vSimp

    ! TODO not using porosity correctly?

    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    tf = p%tf

    ! initialize with starting position
    p%r(0,1) = pt
    p%r(0,2) = px
    p%r(0,3) = py

    write(*,'(A)') '  '
    write(*,'(A)')    '*******************************' 
    write(*,'(A,I0)') 'rk integration, particle ',p%id
    write(*,'(A)')    '*******************************' 

    ! Runge-Kutta 4th order integration scheme (non-adaptive)
    numdt = ceiling((tf - pt)/dt)
    write(*,'(A,ES11.5,A,I0)') 'step size=',dt,' number steps needed=',numdt

    rk: do i = 1,numdt   
       if (mod(i,20) == 0) write(*,'(I0,A,ES12.6E2)') i,' t=',pt

       ! see if particle will reach end this step
       if (pt + dt >= tf) then
          write(*,'(A,I3,A,ES12.6E2)') &
               & 'particle',p%id,' reached final time:',tf
          exit rk
       end if

       ! forward Euler 1/2-step  (predictor)
       ilogc = ceiling(log10(pt))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(px,py,DP),s(:,ilogc),los,his,dom,c,e,bg)
       vInit(1:2) = invlap(pt, tee(ilogc), velp(:,1:2), sol%INVLT)

       FwdEuler(1:2) = [px,py] + dt/2.0*VInit(1:2)

       ! backward Euler 1/2-step (corrector)
       ilogc = ceiling(log10(pt + dt/2.0))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(FwdEuler(1),FwdEuler(2),DP), &
            & s(:,ilogc),los,his,dom,c,e,bg)
       vBkwdEuler(1:2) = invlap(pt + dt/2.0, tee(ilogc), velp(:,1:2), sol%INVLT)

       BkwdEuler(1:2) = [px,py] + dt/2.0*vBkwdEuler(1:2)

       ! midpoint rule full-step predictor
       ilogc = ceiling(log10(pt + dt))
       los = (ilogc-lo)*ns + 1
       his = los + ns - 1
      
       velp(1:ns,1:2) = velCalc(cmplx(BkwdEuler(1),BkwdEuler(2),DP), &
            & s(:,ilogc),los,his,dom,c,e,bg)
       vMidpt(1:2) = invlap(pt + dt, tee(ilogc), velp(:,1:2), sol%INVLT)

       Midpt(1:2) = [px,py] + dt*vMidpt(1:2)

       ! Simpson's rule full-step corrector
       velp(1:ns,1:2) = velCalc(cmplx(Midpt(1),Midpt(2),DP), &
            & s(:,ilogc),los,his,dom,c,e,bg)
       vSimp(1:2) = invlap(pt + dt, tee(ilogc), velp(:,1:2), sol%INVLT)

       Simp(1:2) = [px,py] + dt/6.0* &
            & (vinit(1:2) + 2.0*vBkwdEuler(1:2) + 2.0*vMidpt(1:2) + vSimp(1:2))

       pt = pt + dt
       px = Simp(1)
       py = Simp(2)
       
       p%r(i,1) = pt
       p%r(i,2) = px
       p%r(i,3) = py
       p%r(i-1,4) = vSimp(1)
       p%r(i-1,5) = vSimp(2)
       
       if (sinkCheck(px,py,c,e)) then
          write(*,'(A,I0,A,ES12.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          exit rk             
       end if      
    end do rk
    p%numt = i-1
  end subroutine rungekutta

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine fwdEuler(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : invlap => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : velCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: s(1:,lo:) !! matrix of Laplace parameter np by nlogcycles
    type(particle), intent(inout) :: p
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: sol
    type(domain), intent(in) :: dom
    type(element), intent(in) :: bg

    complex(DP), dimension(1:size(s,1),2) :: velp
    integer :: i, numdt, ilogc, los, his, ns
    real(DP) :: pt, px, py, dt, tf
    real(DP), dimension(2) :: vel

    ! TODO not using porosity correctly?

    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    tf = p%tf

    ! initialize with starting position
    p%r(0,1) = pt
    p%r(0,2) = px
    p%r(0,3) = py

    ! 1st order fwd Euler
    numdt = ceiling((tf - pt)/dt)
    
    write(*,'(A)')    '****************************************'
    write(*,'(A,I0)') 'forward Euler integration, particle ', p%id
    write(*,'(A)')    '****************************************'

    fe: do i = 1,numdt   
       if( mod(i,100) ==0 ) write(*,'(I0,A,ES12.6E2)') i,' t=',pt

       ! see if particle will reach end this step
       if (pt + dt >= tf) then
          write(*,'(A,I0,A,ES12.6E2)') &
               &'particle',p%id,' reached specified ending time:',tf
          exit fe
       end if

       ! full step forward Euler
       ilogc = ceiling(log10(pt))

       los = (ilogc-lo)*ns + 1
       his = los + ns - 1

       velp(1:ns,1:2) = velCalc(cmplx(px,py,DP),s(:,ilogc),los,his,dom,c,e,bg)
       vel(1:2) = invlap(pt,tee(ilogc),velp(:,1:2),sol%INVLT)

       px = px + dt*vel(1)
       py = py + dt*vel(2)
       pt = pt + dt

       p%r(i,1) = pt
       p%r(i,2) = px
       p%r(i,3) = py
       p%r(i-1,4) = vel(1)
       p%r(i-1,5) = vel(2)

       if (sinkCheck(px,py,c,e)) then
          write(*,'(A,I0,A,ES12.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          exit fe
       end if
    end do fe
    p%numt = i - 1
          
  end subroutine fwdEuler

  !###########################################################################
  ! re-allocate particle%r array to longer first dimension

  subroutine reallocate(p,newub)
    use constants, only : DP
    use type_definitions, only : particle
    implicit none

    integer, intent(in) :: newub
    type(particle), intent(inout) :: p
    real(DP), allocatable :: tmp(:,:)

    intrinsic :: move_alloc

    allocate(tmp(0:newub,5))
    tmp(0:ubound(p%r,dim=1),1:5) = p%r(0:,1:5)
    
    ! see section 17.5.3 of Metcalf, Reid & Cohen
    call move_alloc(tmp,p%r)

  end subroutine reallocate

  !###########################################################################
  ! check if particle has moved into a sink

  function sinkCheck(px,py,c,e) result(partEnd)
    use constants, only : DP
    use type_definitions, only : ellipse, circle
    implicit none

    real(DP), intent(in) :: px,py
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    logical :: partEnd

    partEnd = .false.

    ! did the particle enter an ibnd==2 circle (well)?
    
    if (any(c%ibnd == 2) .and. any(abs(cmplx(px,py,DP) - cmplx(c%x,c%y,DP)) < c%r)) then
       partEnd = .true.
    end if
    
    ! TODO add check for flowing into ibnd==2 ellipse (line sink)

    ! TODO handle flowing into a constant head/flux element (from inside or
    ! TODO from outside, circles or ellipses

    ! TODO dealing with area-sinks is more complex, will be dealt with later
  end function sinkCheck

end module particle_integrate

