! $Id: invlap_test_complex.f90,v 1.2 2007/11/07 20:06:18 kris Exp kris $

! program for testing de Hoog invlap method,
! but this test program is very general and would handle any method 
! in a module with the two functions called below

! main test program
program invlap_test
  use constants, only : DP, PI, RONE, RTWO
  use inverse_Laplace_Transform, only : invlap => deHoog_invlap , pvalues => deHoog_pvalues
  implicit none

  integer, parameter :: NFUN = 14
  character(50) :: filename
  real(DP), allocatable :: ft(:),tft(:),t(:)
  complex(DP), allocatable :: fp(:),p(:)
  integer :: i,j,nt,M
  real(DP) :: aa,Minput, alpha,tol, tee

  !! read parameters from input file
  open(unit=22,file='invlap_dehoog.in',status='old',action='read')
  read(22,*) tee, alpha, tol, Minput
  read(22,*) aa, nt
  read(22,*) filename

  M = nint(Minput)
  allocate(fp(0:2*M+1),p(0:2*M),t(1:nt),ft(1:nt),tft(1:nt))

  read(22,*) t(1:nt)
  close(22)

  !! compute necessary values of p
  p(0:2*M) = pvalues(tee, alpha, tol, M)

!!$  write(*,*) 't:',t
!!$  do j=0,N
!!$     write(*,*) j,p(j)
!!$  end do

  open(UNIT=20, FILE=filename, STATUS='REPLACE', ACTION='WRITE')  !! for plotting
  open(unit=77, file='invlap_norms.out', status='replace', action='write')  !! for automating

  ! output is in a gnuplot-friendly format 
  ! with each test function a different index

  write (20,*) '# inverse Laplace transform benchmark (deHoog method)'
  write (20,*) '# alpha:',alpha, 'tee:',tee, ' tol:',tol, ' M:',M, ' num_t:',nt

  do j=1,NFUN    
     write (20,*) '# function:', j 
     write (20,*) '#      t               true&
          &                inverted               diff'

     tft(1:nt) = foft(j,aa,t(1:nt))  !! true solution in time domain

     fp(0:2*M) = fofp(j,aa,p(0:2*M)) !! true solution in Laplace domain
     ft(1:nt) = invlap(alpha,tol,t(1:nt),tee,fp(0:2*M),M)  !! inverted time domain solution
     
     !! write L and Le norms of Davies and Martin (1979, eqns 2.1 & 2.2)
     write(77,888) j,sqrt(sum((ft(1:nt) - tft(1:nt))**2/30.0_DP)), &
                 & sqrt(sum((ft(1:nt) - tft(1:nt))**2*exp(-t(1:nt)))/&
                 & sum(exp(-t(1:nt))))

     do i=1,nt
        write(20,999) t(i), tft(i), ft(i), tft(i)-ft(i)
     end do
     write (20,'(//)')

  end do
  close(20)
  close(77)
888 format (i2,2(1x,ES22.14E3))
999 format (ES15.7E1,3(1x,ES22.14E3))

contains
  !! time domain functions
  function foft(idx,a,t)
    use constants, only : DP,RTWO,PI, RZERO,RONE
    use expint_approx

    integer, intent(in) :: idx
    real(DP), intent(in) :: a
    real(DP), intent(in), dimension(:) :: t
    real(DP), dimension(size(t)) :: foft

    real(dP), parameter :: w = 50.0_DP
    integer :: nt, i
    intrinsic :: derf,derfc,dbesj0
    nt = size(t)

    !! not vectorized, since built-in functions can't handle it
    do i=1,nt
       select case (idx)
       case(1)
          ! exponential * error function (A & S 29.3.40)
          foft(i) = exp(a**2*t(i))*derfc(a*sqrt(t(i)))
       case(2)
          ! cosine/sqrt(t), D&M f_2 
          foft(i) = cos(RTWO*sqrt(t(i)))/sqrt(PI*t(i))
       case(3)
          ! sine (A & S 29.3.22)
          foft(i) = t(i)/(RTWO*a)*sin(a*t(i))
       case(4)
          ! a step 'on' function
          if (t(i)<=a) then
             foft(i)=RZERO
          else
             foft(i)=RONE
          end if
       case(5)
          ! J_0, D&M f_1
          foft(i) = dbesj0(t(i))
       case(6)
          ! a stiff ODE
          foft(i) = (exp(-a*t(i))-exp(-10.0_DP*a*t(i)))
       case(7)
          ! a square wave of period a
          if (mod(ceiling(t(i)/(a/RTWO)),2) == 0) then
             foft(i)=RZERO
          else
             foft(i)=RONE
          end if
       case(8)
          ! a step 'off' function
          if (t(i)>=a) then
             foft(i)=RZERO
          else
             foft(i)=RONE
          end if
       case(9)
          ! exponential, D&M f_15
          foft(i) = RTWO*exp(-4.0_DP/t(i))/sqrt(PI*t(i)**3)
       case(10)
          ! constant
          foft(i) = RONE
       case(11)
          ! linear
          foft(i) = t(i)
       case(12)
          ! Theis solution, at r=1/a
          foft(i) = expint(1.0_DP/(a**2*4.0_DP*t(i)))/(4.0_DP*PI)
       case(13)
          ! "steep" error function step on at zero
          foft(i) = 0.5_DP*(1.0_DP + derf(t(i)*w))
       case(14)
          ! "steep" error function step on at a
          foft(i) = 0.5_DP*(1.0_DP + derf((t(i)-a)*w))
       case default
          stop 'incorrect index'   
       end select
    end do
  end function foft

  !! laplace domain functions
  function fofp(idx,a,p)
    use constants, only : DP,PI,RONE
    use complex_bessel, only :cbesk
    use double_complex_erf, only : dcerf
  
    integer, intent(in) :: idx
    real(DP), intent(in) :: a
    complex(DP), intent(in), dimension(:) :: p
    complex(DP), dimension(size(p)) :: fofp
    complex(DP), dimension(size(p)) :: ans, arg

    real(DP), dimension(2) :: erfc
    real(DP), parameter :: w = 50.0_DP
    integer :: np,j,nz,ierr
    np = size(p)

    !! these are vectorized wrt p
    select case(idx)
    case(1)
       fofp = RONE/(p+sqrt(p)*a)
    case(2)
       fofp = exp(-RONE/p)/sqrt(p)
    case(3)
       fofp = p/(p**2 + a**2)**2
    case(4)
       fofp = exp(-a*p)/p
    case(5)
       fofp = RONE/sqrt(p**2 + RONE)
    case(6)
       fofp = RONE/(p+a) - RONE/(p+10.0_DP*a)
    case(7)
       fofp = RONE/(p+p*exp(-a/2.0_DP*p))
    case(8)
       fofp = (RONE - exp(-a*p))/p
    case(9)
       fofp = exp(-4.0_DP*sqrt(p))
    case(10)
       fofp = RONE/p
    case(11)
       fofp = RONE/p**2
    case(12)  !! Theis soln at r=1/a
       do j=1,np
          call cbesk(sqrt(p(j))/a,0.0_DP,1,1,ans(j),nz,ierr)
          if (nz + ierr /= 0) print *, 'CBESK error',nz,ierr
       end do
       fofp = ans/(2.0_DP*PI*p)
    case(13)
       arg(1:np) = p(1:np)/(2.0_DP*w)
       do j=1,np
          call dcerf(1,[real(arg(j)),aimag(arg(j))],erfc(1:2))
          fofp(j) = 0.5_DP*(1.0_DP + exp(arg(j)**2)* &
               & cmplx(erfc(1),erfc(2),DP))/p(j)
       end do
    case(14)
       arg(1:np) = p(1:np)/(2.0_DP*w)
       do j=1,np
          call dcerf(1,[real(arg(j)),aimag(arg(j))],erfc(1:2))
          fofp(j) = 0.5_DP*(1.0_DP + exp(arg(j)**2)* &
               & cmplx(erfc(1),erfc(2),DP))/p(j) * exp(-a*p(j))
       end do
    case default
       stop 'incorrect index'   
    end select
  end function fofp

end program invlap_test

