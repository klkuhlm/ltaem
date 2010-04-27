! $Id: invlap.f90,v 1.1 2007/03/15 00:17:36 kris Exp kris $
module inverse_Laplace_Transform
  use constants, only : DP
  implicit none

  private
  public :: deHoog_invlap, deHoog_pvalues
  
  interface deHoog_invlap
     module procedure deHoog_invlap_vect, deHoog_invlap_scal
  end interface
  
contains
  
  !! an implementation of the de Hoog method
  
  !! assumes proper f(p) have been computed for the p
  !! required for the vector of t passed to this function
  !! -- only one log-cycle of time should be passed at once --
  !! (no error checking done in this regard)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function deHoog_invLap_vect(alpha,tol,t,tee,fp,M) result(ft)
    use constants, only : DP, PI, CZERO, CONE, EYE, HALF, RONE, RZERO
    real(DP), intent(in) :: alpha  ! abcissa of convergence
    real(DP), intent(in) :: tol    ! desired accuracy
    real(DP), intent(in) :: tee    ! scaling factor
    real(DP), intent(in), dimension(:) :: t   ! vector of times
    integer, intent(in) :: M
    complex(DP), intent(in), dimension(0:2*M) :: fp
    real(DP), dimension(size(t)) :: ft ! output

    real(DP) :: gamma
    integer :: r, rq, n, max, nt
    complex(DP), dimension(0:2*M,0:M) :: e
    complex(DP), dimension(0:2*M,1:M) :: q
    complex(DP), dimension(0:2*M) :: d
    complex(DP), dimension(-1:2*M,size(t)) :: A,B
    complex(DP), dimension(size(t)) :: z,brem,rem
    
    nt = size(t)

    ! there will be problems is fp(:)==0
    if(maxval(abs(fp)) > tiny(1.0_DP)) then

       ! Re(p) -- this is the de Hoog parameter c
       gamma = alpha - log(tol)/(2.0_DP*tee)

       ! initialize Q-D table 
       e(0:2*M,0) = CZERO
       q(0,1) = fp(1)/(fp(0)*HALF) ! half first term
       q(1:2*M-1,1) = fp(2:2*M)/fp(1:2*M-1)

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
       d(0) = fp(0)*HALF ! half first term
       forall(r = 1:M)
          d(2*r-1) = -q(0,r) ! even terms
          d(2*r)   = -e(0,r) ! odd terms
       end forall

       ! seed A and B vectors for recurrence
       A(-1,1:nt) = CZERO
       A(0,1:nt) = d(0)
       B(-1:0,1:nt) = CONE

       ! base of the power series
       z(1:nt) = exp(EYE*PI*t(:)/tee)

       ! coefficients of Pade approximation
       ! using recurrence for all but last term
       do n = 1,2*M-1
          A(n,:) = A(n-1,:) + d(n)*A(n-2,:)*z(:)
          B(n,:) = B(n-1,:) + d(n)*B(n-2,:)*z(:)
       end do

       ! "improved remainder" to continued fraction
       brem(1:nt) = (CONE + (d(2*M-1) - d(2*M))*z(:))*HALF
       rem(1:nt) = -brem*(CONE - sqrt(CONE + d(2*M)*z(:)/brem**2))

       ! last term of recurrence using new remainder
       A(2*M,:) = A(2*M-1,:) + rem*A(2*M-2,:)
       B(2*M,:) = B(2*M-1,:) + rem*B(2*M-2,:)

       ! diagonal Pade approximation
       ! F=A/B represents accelerated trapezoid rule
       ft(1:nt) =  exp(gamma*t(:))/tee * real(A(2*M,:)/B(2*M,:))

    else  !! entire f(p) vector is zero
       ft = RZERO
       write(*,*) 'f(t) not computed at t=',t, ' because max|fp|=', maxval(abs(fp))

    end if

  end function deHoog_invLap_vect

  function deHoog_invLap_scal(alpha,tol,t,tee,fp,M) result(ft)
    use constants, only : DP
    real(DP), intent(in) :: alpha, tol, t, tee 
    integer, intent(in) :: M
    complex(DP), intent(in), dimension(0:2*M) :: fp
    real(DP) :: ft ! output
    
    ft = sum(deHoog_invLap_vect(alpha,tol,[t],tee,fp,M))

  end function deHoog_invLap_scal
  
  function deHoog_pvalues(tee,alpha,tol,M) result(p)
    use constants, only : DP, PI
    real(DP), intent(in) :: tee,alpha,tol
    integer, intent(in) :: M
    complex(DP), dimension(2*M+1) :: p

    integer :: i

    p(1:2*M+1) = cmplx(alpha-log(tol)/(2.0_DP*tee), PI*(/( i, i=0,2*M )/)/tee, DP)

  end function deHoog_pvalues
  

end module inverse_Laplace_Transform
