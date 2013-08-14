!
! Copyright (c) 2011,2012,2013 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

module particle_integrate
  implicit none
  private
  public :: rungekuttamerson, rungekutta, fwdEuler, analytic

contains

  subroutine rungeKuttaMerson(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : L => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : V => velcalc
    use utility, only : v2c
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

    integer :: count, lt, los, his, ns
    logical :: partEnd, done
    real(DP), dimension(2) :: vInit,vTrap,vAB3,vAB2,vSF
    real(DP), dimension(2) :: fwdEuler,Trap,halfAB,fullAB,Simp,x0
    real(DP) :: error, pt, px, py, dt, length, enddt
    real(DP), parameter :: SAFETY = 0.9

    !!complex(DP), dimension(21,2) :: tmp

    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    if (.not. p%forward) dt = -dt ! backwards particle tracking

    p%r(0,1:3) = [pt,px,py]
    count = 1

    partEnd = .false.
    write(*,'(A,I0)') '** rkm integration, particle ',p%id
    
    ! Runge-Kutta-Merson 4th-order adaptive integration scheme
    rkm: do
       ! did particle enter a sink in last step?
       if (partEnd) exit rkm
 
       ! will particle pass specified end time during this step?
       call trackDone(p%forward,pt,dt,p%tf,done,enddt) 
       if (done) then
          print *, 'p%forward,pt,dt,p%tf,done,enddt',p%forward,pt,dt,p%tf,done,enddt
          if (enddt < spacing(dt)) then
             ! close enough
             write(*,'(A,I0,2(A,ES13.6E2))') &
                  & 'particle ',p%id,' reached final time:',pt,'=',p%tf

             ! compute velocity at final location
             call getsrange(pt,lo,ns,los,his,lt)
             p%r(count-1,4:5) = L(pt,tee(lt), &
                  & V(x0,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
             
             exit rkm
          else
             dt = enddt
          end if
       end if
       
       x0 = [px,py]

       ! forward Euler 1/3-step  predictor
       call getsrange(pt,lo,ns,los,his,lt)
       vInit = L(pt,tee(lt), V(x0,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       FwdEuler(1:2) = x0(1:2) + dt/3.0*vInit(1:2)

       ! trapezoid rule 1/3-step corrector
       call getsrange(pt + dt/3.0,lo,ns,los,his,lt)
       vTrap = L(pt+dt/3.0,tee(lt), V(FwdEuler,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       Trap(1:2) = x0(1:2) + dt/6.0*(vInit(1:2) + vTrap(1:2))

       ! Adams-Bashforth 1/2-step predictor
       vAB3 = L(pt+dt/3.0,tee(lt), V(Trap,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       halfAB(1:2) = x0(1:2) + dt/8.0*(vInit(1:2) + 3*vAB3(1:2))

       ! full step Adams-Bashforth predictor
       call getsrange(pt + dt/2.0,lo,ns,los,his,lt)
       vAB2 = L(pt+dt/2.0,tee(lt), V(halfAB,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       fullAB(1:2) = x0(1:2) + dt/2.0*(vInit(1:2) - 3*vAB3(1:2) + 4*vAB2(1:2))

       ! full step Simpson's rule corrector
       call getsrange(pt + dt,lo,ns,los,his,lt)
       vSF = L(pt+dt,tee(lt), V(halfAB,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       Simp(1:2) = x0(1:2) + dt/6.0*(vInit(1:2) + 4*vAB2(1:2) + vSF(1:2))

       ! relative error (biggest of x- or y-direction)
       error = maxval(abs((Simp - fullAB)/Simp))

       ! magnitude of total step taken (straight line)
       length = abs(v2c(Simp) - v2c(x0))

       ! only advance to next step if error level is acceptable
       ! _and_ resulting step is less than prescribed limit
       ! _or_ we have gotten to the minimum step size (proceed anyways)
       if ((error < p%tol .and. length <= p%maxL) .or. abs(dt) <= p%mindt) then
          pt = pt + dt
          px = Simp(1)
          py = Simp(2)
          p%r(count,1:3) = [pt,px,py]
          ! velocity associated with previous step
          p%r(count-1,4:5) = vSF(1:2)
          count = count + 1

          partEnd = sinkCheck(px,py,c,e)

          if (partEnd) then
             write(*,'(A,I0,A,ES13.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          else
             if (count == ubound(p%r,dim=1)) then
                ! if out of space in results array double the size
                call reallocate(p,ubound(p%r,dim=1)*2)
             end if
          end if
       end if

       ! new step size based on error estimate:
       ! Numerical Recipes, Press et al 1992, eqn 16.2.10, p 712
       if ((.not. partEnd) .and. (.not. done)) then

          if (length > p%maxL) then
             ! cut time step to minimize distance
             ! integrated in one step
             dt = dt/3.0

          elseif (p%tol >= error) then
             ! enlarge dt to speed up integration when possible
             dt = SAFETY*dt*(p%tol/error)**2.0D-1

          else
             ! reduce dt to increase accuracy when needed
             dt = SAFETY*dt*(p%tol/error)**2.5D-1

          end if

          if (abs(dt) <= p%mindt) then
             write(*,'(3(A,ES13.6E2))') 'min step; dt=',dt,' error=',error, ' pt=',pt
             dt = p%mindt
             if (.not. p%forward) dt = -dt
          end if
       end if
    end do rkm
    p%numt = count - 1
  end subroutine rungekuttamerson

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine rungekutta(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : L => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : V => velCalc
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

    integer :: i, numdt, lt, ns, los, his
    real(DP) :: pt, px, py, dt
    real(DP), dimension(2) :: FwdEuler,BkwdEuler,MidPt,Simp,x0
    real(DP), dimension(2) :: vInit,vBkwdEuler,vMidpt,vSimp
    logical :: sink 

    sink = .false.
    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    if (.not. p%forward) dt = -dt ! backwards particle tracking

    ! initialize with starting position
    p%r(0,1:3) = [pt,px,py]

    write(*,'(A,I0)') '** rk integration, particle ',p%id

    ! Runge-Kutta 4th order integration scheme (non-adaptive)
    numdt = ceiling((p%tf - pt)/abs(dt))
    write(*,'(A,ES12.5,A,I0)') 'step size=',dt,' number steps needed=',numdt

    rk: do i = 1,numdt
       if (mod(i,100) == 0) write(*,'(I0,A,ES13.6E2)') i,' t=',pt

       x0 = [px,py]

       ! forward Euler 1/2-step  (predictor)
       call getsrange(pt,lo,ns,los,his,lt)
       vInit = L(pt,tee(lt), V(x0,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       FwdEuler(1:2) = x0(1:2) + dt/2.0*vInit(1:2)

       ! backward Euler 1/2-step (corrector)
       call getsrange(pt + dt/2.0,lo,ns,los,his,lt)
       vBkwdEuler = L(pt+dt/2.0,tee(lt), V(FwdEuler,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       BkwdEuler(1:2) = x0(1:2) + dt/2.0*vBkwdEuler(1:2)

       ! midpoint rule full-step predictor
       call getsrange(pt + dt,lo,ns,los,his,lt)
       vMidpt = L(pt+dt,tee(lt), V(BkwdEuler,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       Midpt(1:2) = x0(1:2) + dt*vMidpt(1:2)

       ! Simpson's rule full-step corrector
       vSimp = L(pt+dt,tee(lt), V(Midpt,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
       Simp(1:2) = x0(:) + dt/6.0*(vinit(:) + 2*vBkwdEuler(:) + 2*vMidpt(:) + vSimp(:))

       pt = pt + dt
       px = Simp(1)
       py = Simp(2)

       p%r(i,1:3) = [pt,px,py]
       p%r(i-1,4:5) = vSimp(1:2)

       if (sinkCheck(px,py,c,e)) then
          write(*,'(A,I0,A,ES13.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          sink = .true.
          exit rk
       end if
    end do rk

    if (.not. sink) then
       write(*,'(A,I0,2(A,ES13.6E2))') &
            & 'particle ',p%id,' reached final time:',pt,'=',p%tf
    
       ! compute velocity at final time
       call getsrange(pt,lo,ns,los,his,lt)
       p%r(numdt,4:5) = L(pt,tee(lt), &
            & V(x0,s(:,lt),los,his,dom,c,e,bg), sol%INVLT)
    end if
    
    p%numt = i-1
  end subroutine rungekutta

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine fwdEuler(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use inverse_laplace_transform, only : L => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : V => velCalc
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

    integer :: i, numdt, lt, los, his, ns
    real(DP) :: pt, px, py, dt
    real(DP), dimension(2) :: vel
    logical :: sink

    sink = .false.
    ns = size(s,dim=1)
    pt = p%ti
    px = p%x
    py = p%y
    dt = p%dt
    if (.not. p%forward) dt = -dt

    ! initialize with starting position
    p%r(0,1:3) = [pt,px,py]

    ! 1st order fwd Euler
    numdt = ceiling((p%tf - pt)/abs(dt))

    write(*,'(A,I0)') '** fwd Euler integration, particle ', p%id

    fe: do i = 1,numdt
       if(mod(i,500) == 0) write(*,'(I0,A,ES13.6E2)') i,' t=',pt

       ! full step forward Euler
       call getsrange(pt,lo,ns,los,his,lt)
       !!!print *, 'FE t>:',pt,' lo>:',lo,' ns>:',ns,' los<:',los,' his<:',his,' lt<:',lt

       vel = L(pt,tee(lt),V([px,py],s(:,lt),los,his,dom,c,e,bg),sol%INVLT)

       px = px + dt*vel(1)
       py = py + dt*vel(2)
       pt = pt + dt

       p%r(i,1:3) = [pt,px,py]
       p%r(i-1,4:5) = vel(1:2)

       if (sinkCheck(px,py,c,e)) then
          write(*,'(A,I0,A,ES13.6E2)') 'particle ',p%id,' entered a sink at t=',pt
          sink = .true.
          exit fe
       end if
    end do fe
    
    if (.not. sink) then
       write(*,'(A,I0,2(A,ES13.6E2))') &
            &'particle',p%id,' reached final time:',pt,'=',p%tf
       ! compute velocity at final time
       call getsrange(pt,lo,ns,los,his,lt)
       p%r(numdt,4:5) = L(pt,tee(lt),&
            & V([px,py],s(:,lt),los,his,dom,c,e,bg),sol%INVLT)
    end if
    
    p%numt = i - 1

  end subroutine fwdEuler

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine analytic(s,tee,c,e,bg,sol,dom,p,lo)
    use constants, only : DP
    use utility, only : outer
    use inverse_laplace_transform, only : L => deHoog_invlap
    use type_definitions, only : circle, ellipse, element, solution, particle, domain
    use calc_routines, only : V => velCalc
    implicit none

    integer, intent(in) :: lo
    real(DP), intent(in) :: tee(lo:)
    complex(DP), intent(in) :: s(1:,lo:) 
    type(particle), intent(inout) :: p
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: sol
    type(domain), intent(in) :: dom
    type(element), intent(in) :: bg

    integer :: i, lt, los, his, ns
    real(DP), dimension(3) :: t
    real(DP), dimension(2) :: loc, vel
    real(DP), dimension(2,3) :: ploc     ! location guesses
    ! f(x) and f(y) for root-finding
    complex(DP), dimension(size(s,1),2,3) :: fpv
    complex(DP), dimension(size(s,1),2) :: vv, fx3x2, fx3x1, fx2x1, fx3x2x1, w, r
    real(DP), parameter :: DTFRAC = 1.1, ITERTOL = 1.0E-3
    integer, parameter :: MAXITER = 100

    ns = size(s,dim=1)
    ! assume particles start at t=0?
    loc(1:2) = [p%x,p%y] ! starting position

    ! initialize with starting position
    p%r(0,1:3) = [0.0_DP,loc(:)]

    ! compute particle location at time t,
    ! given initial location at t=0, and desired time t=dt

    ! solution at t=ti + dt is a root-finding exercise in x and y
    ! use fwd euler to solve for initial guess where particle will be
    
    t = [p%ti + p%dt, p%dt/DTFRAC, p%dt*DTFRAC]

    call getsrange(t(1),lo,ns,los,his,lt)

    print *, 'analytical t>:',t(1),' lo>:',lo,' ns>:',ns,' los<:',los,' his<:',his,' lt<:',lt

    vv(:,1:2) = V(loc,s(:,lt),los,his,dom,c,e,bg)

    print *, 'analytical 002'  

    do i = 1, 3
       ! full step forward Euler to three slightly
       ! different times, re-using v(p) estimate, vv
       vel = L(t(i),tee(lt),vv(:,1:2),sol%INVLT)
       ! estimated coordinates
       ploc(1:2,i) = loc(:) + t(i)*vel(:)
       ! equation which is set to zero, and used for root-finding (p*x - x0) = V(x)
       fpv(1:ns,1:2,i) = V(ploc(:,i),s(:,lt),los,his,dom,c,e,bg) - &
            & (outer(s(:,lt),ploc(:,i)) - spread(loc(:),1,ns))
    end do
    
    print *, 'analytical 003'  

    ! Muller's algorithm for roots (Good for complex roots)
    ! compute roots for x & y, for each value of laplace parameter

    mul: do i = 1, MAXITER
       ! calculate divided differences
       fx3x2 = (fpv(:,:,2)-fpv(:,:,3))/spread(ploc(:,2)-ploc(:,3),1,ns)
       fx3x1 = (fpv(:,:,1)-fpv(:,:,3))/spread(ploc(:,1)-ploc(:,3),1,ns)
       fx3x1 = (fpv(:,:,1)-fpv(:,:,2))/spread(ploc(:,1)-ploc(:,2),1,ns)

       w = fx3x2 + fx3x1 - fx2x1
       fx3x2x1 = (fx2x1-fx3x2)/spread(ploc(:,1)-ploc(:,3),1,ns)

       if (any(abs(w) < epsilon(0.0) .and. &
             & abs(fx3x2x1) < epsilon(0.0))) then
          print *, 'cancelled with',t,ploc
          stop
       end if
       
       ploc(:,1) = ploc(:,2)
       fpv(:,:,1) = fpv(:,:,2)
       ploc(:,2) = ploc(:,3)
       fpv(:,:,2) = fpv(:,:,3)

       ! denominator should be as large as possible => choose sign
       r = sqrt(w**2 - 4.0*spread(ploc(:,3),1,ns)*fx3x2x1)
       where (abs(w - r) > abs(w + r))
          r = -r
       end where 
       ploc(:,3) = ploc(:,3) - 2.0*L(t(i),tee(lt),fpv(:,:,3)/(w+r),sol%INVLT)
       fpv(:,:,3) = V(ploc(:,3),s(:,lt),los,his,dom,c,e,bg)  - &
            & (outer(s(:,lt),ploc(:,3)) - spread(loc(:),1,ns))
       if (all(abs(ploc(:,3) - ploc(:,2)) < ITERTOL)) then
          ! compute inverse Laplace transform of optimized location
          vel = L(t(i),tee(lt),fpv(:,:,3),sol%INVLT)
          p%r(1,1) = p%dt
          p%r(1,2:3) = p%r(0,2:3) + p%dt*vel(:)
          exit mul
       end if
    end do mul

  end subroutine analytic


  !###########################################################################
  ! logical function for comparing tracking time to final time

  pure subroutine trackDone(fwd,tp,dt,te,done,newdt)
    use constants, only : DP
    logical, intent(in) :: fwd  ! is tracking forward or backwards?
    real(DP), intent(in) :: tp,dt,te ! particle and "end" time
    logical, intent(out) :: done
    real(DP), intent(out) :: newdt ! delta t to land on final time

    if (fwd) then
       ! forward particle tracking
       if (tp+dt >= te) then
          done = .true.
          newdt = te - tp
       else
          done = .false.
          newdt = 999.
       end if
    else
       ! backwards particle tracking
       if (tp-dt <= te) then
          done = .true.
          newdt = tp - te   
       else
          done = .false.
          newdt = 999.
       end if
    end if

  end subroutine trackDone


  !###########################################################################
  ! some common range-computing code

  subroutine getSRange(t,lo,ns,los,his,lt)
    use constants, only : DP
    real(DP), intent(in) :: t
    integer, intent(in) :: lo,ns
    integer, intent(out) :: los,his,lt

    lt = ceiling(log10(t))  ! what log-cycle does time fall into?
    los = (lt-lo)*ns + 1    ! low and high indices on p or s
    his = los + ns - 1

  end subroutine getSRange

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
  ! check if particle has moved into a sink (left flow domain)

  ! TODO: modify sinkCheck or create sourceCheck for reverse tracking
  function sinkCheck(px,py,c,e) result(partEnd)
    use constants, only : DP, EYE
    use utility, only : cacosh
    use type_definitions, only : ellipse, circle
    implicit none

    real(DP), intent(in) :: px,py
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    logical :: partEnd
    complex(DP) :: pz

    partEnd = .false.
    pz = cmplx(px,py,DP)

    ! TODO: try to check whether particle passed near a well during last time step
    ! TODO: then pass back a signal to take a smaller time step

    ! did the particle enter an ibnd==2 circle (well)?
    if (any(c%ibnd == 2 .and. abs(pz - c%z) <= c%r)) then
       partEnd = .true.
       goto 999
    end if

    ! did the particle enter an ibnd==2 ellipse (line sink)?
    if (any(e%ibnd == 2 .and. real(cacosh((pz-e%z)*exp(-EYE*e%theta)/e%f)) <= e%r)) then
       partEnd = .true.
       goto 999
    end if

    ! TODO handle flowing into a constant head/flux element (from inside or
    !      from outside, circles or ellipses

    ! TODO dealing with area-sinks is more complex, will be dealt with later

999 continue
  end function sinkCheck

end module particle_integrate

