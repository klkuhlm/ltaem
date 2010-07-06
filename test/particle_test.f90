program particle_test
  use constants, only : DP, PI, RTWO, RZERO
  use complex_bessel, only : cbesk

  implicit none

  integer, parameter :: MAXITER = 100
  real(DP), parameter :: ITERTOL = 1.0D-10
  integer, parameter :: NUMT = 1

  !! de Hoog method related variables
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2M+1 terms in fourier series approx (more terms -> smaller Gibbs effects)
  integer, parameter :: M = 8
  ! all singularities to the left of alpha (small but non-zero)
  ! desired accuracy (too small can cause overflow for some special functions??)  
  real(DP) :: DHalpha, tol, tee

  real(DP), dimension(NUMT) :: t

  real(DP) :: rinit
  real(DP) :: alpha, K, pump

  complex(DP) :: rlast, r, q
  complex(DP), dimension(2*M+1) :: ploc
  real(DP), dimension(NUMT) :: soln

  ! Laplace parameter values
  complex(DP), dimension(2*M+1) :: p
  complex(DP), dimension(1) :: K1
  
  integer :: i, j, ii, np, nz, ierr
  logical :: converge

  ! material parameters
  K = 1.0D0
  alpha = K/1.0D-3
  pump = 1.0D0  ! pumping rate

  ! time that particle location is wanted 
!!$  t = 1.0D1**[(real(i,DP)/real(NUMT,DP), i=0,NUMT-1)] * 0.5_DP
  t = 0.25_DP

  ! initial radial distance for particle
  rinit = 2.0D0

  DHalpha = RZERO
  tol = 1.0D-9
  np = 2*M+1

  open(unit=50,file='particle-test.out',status='replace',action='write')

  do ii = 1, NUMT

     ! parameter for de Hoog method
     tee = t(ii)*RTWO

     p(1:2*M+1) = cmplx(DHalpha - log(tol)/(RTWO*tee), real([(i,i=0,2*M)],DP)*PI/tee, DP)

     ! loop over values of p
     do j = 1, np
        print *, 'p:',j,p(j)

        q = sqrt(p(j)/alpha) 

        ! initialize search at starting location
!!$     if (ii == 1) then
        rlast = rinit/p(j)
!!$     else
!!$        rlast = soln(ii-1)
!!$     end if

        print *, 'initial p:', rlast
        converge = .false.
        iter: do i = 1, MAXITER

           call cbesk(rlast*q, 1.0D0, 1, 1, K1(1), nz, ierr)
           if (ierr /= 0) print *, 'Bessel function error:',ierr
           if (nz /=0) print *, 'Bessel function underflow:',nz

           r = K*pump*q/(RTWO*PI*p(j)**2)*K1(1) + rinit/p(j)

!!$           print *, 'iter:',i, r

           if (abs(r - rlast) < ITERTOL) then
              if(converge) then
                 ploc(j) = r
                 exit iter
              else
                 converge = .true.
              end if
           else
              rlast = r
              converge = .false.
           end if

           if(i == MAXITER) print *, 'failed to converge'

        end do iter

     end do

     ! do numerical inverse transform for each time
     soln(ii) = deHoog_invlap(DHalpha,tol,t(ii),tee,ploc(:),M,.false.)
     write(50,222) t(ii), soln(ii)

  end do

222 format(2(1X,ES22.14))


contains

  !! an implementation of the de Hoog method
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function deHoog_invLap(alpha,tol,t,tee,fp,M,smooth) result(ft)
    use constants, only : DP, PI, CZERO, CONE, EYE, HALF, RONE
    real(DP), intent(in) :: alpha  ! abcissa of convergence
    real(DP), intent(in) :: tol    ! desired accuracy
    real(DP), intent(in) :: t,tee  ! time and scaling factor
    integer, intent(in) :: M
    complex(DP), intent(in), dimension(0:2*M) :: fp
    logical, intent(in) :: smooth
    real(DP) :: ft ! output

    real(DP), dimension(0:2*M) :: lanczos
    real(DP) :: gamma
    integer :: r, rq, n, max, j
    complex(DP), dimension(0:2*M,0:M) :: e
    complex(DP), dimension(0:2*M,1:M) :: q
    complex(DP), dimension(0:2*M) :: d
    complex(DP), dimension(-1:2*M) :: A,B
    complex(DP) :: brem,rem,z
    
    ! there will be problems is fp(:)=0
    if(maxval(abs(fp)) > epsilon(1.0D0)) then

       lanczos(:) = RONE
       if(smooth) then
          do j=1,2*M
             lanczos(j) = sin(PI*real(j,DP)/real(2*M+1,DP))/(PI*real(j,DP)/real(2*M+1,DP))
          end do
       end if
       
       ! Re(p) -- this is the de Hoog parameter c
       gamma = alpha - log(tol)/(2.0_DP*tee)

       ! initialize Q-D table 
       e(0:2*M,0) = CZERO
       q(0,1) = fp(1)*lanczos(1)/(lanczos(0)*fp(0)*HALF) ! half first term
       q(1:2*M-1,1) = fp(2:2*M)*lanczos(2:2*M)/(fp(1:2*M-1)*lanczos(1:2*M-1))

       ! rhombus rule for filling in triangular Q-D table
       do r = 1,M
          ! start with e, column 1, 0:2*M-2
          max = 2*(M-r)
          e(0:max,r) = q(1:max+1,r) - q(0:max,r) + e(1:max+1,r-1)
          if (r /= M) then
             ! start with q, column 2, 0:2*M-3
             rq = r+1
             max = 2*(M-rq)+1
             q(0:max,rq) = q(1:max+1,rq-1) * e(1:max+1,rq-1) / e(0:max,rq-1)
          end if
       end do

       ! build up continued fraction coefficients
       d(0) = lanczos(0)*fp(0)*HALF ! half first term
       forall(r = 1:M)
          d(2*r-1) = -q(0,r) ! even terms
          d(2*r)   = -e(0,r) ! odd terms
       end forall

       ! seed A and B vectors for recurrence
       A(-1) = CZERO
       A(0) = d(0)
       B(-1:0) = CONE

       ! base of the power series
       z = exp(EYE*PI*t/tee)

       ! coefficients of Pade approximation
       ! using recurrence for all but last term
       do n = 1,2*M-1
          A(n) = A(n-1) + d(n)*A(n-2)*z
          B(n) = B(n-1) + d(n)*B(n-2)*z
       end do

       ! "improved remainder" to continued fraction
       brem = (CONE + (d(2*M-1) - d(2*M))*z)*HALF
       rem = -brem*(CONE - sqrt(CONE + d(2*M)*z/brem**2))

       ! last term of recurrence using new remainder
       A(2*M) = A(2*M-1) + rem*A(2*M-2)
       B(2*M) = B(2*M-1) + rem*B(2*M-2)

       ! diagonal Pade approximation
       ! F=A/B represents accelerated trapezoid rule
       ft =  exp(gamma*t)/tee * real(A(2*M)/B(2*M))

    else
       ft = RZERO
    end if

  end function deHoog_invLap


end program particle_test
