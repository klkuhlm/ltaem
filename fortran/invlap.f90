!
! Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

! this module implements the F.R. de Hoog, J.H. Knight, and A.N. Stokes
! numerical inverse Laplace transform algorithm.
! see "An improved method for numerical inversion of Laplace
!     transforms", SIAM J. Sci. Stat. Comp., 3, 357-366, 1982.

module inverse_Laplace_Transform
  implicit none

  private
  public :: deHoog_invlap, deHoog_pvalues

  interface deHoog_invlap
     module procedure deHoog_invlap_vect, deHoog_invlap_scalt, &
          & deHoog_invlap_vel_scalt, deHoog_invlap_vel
  end interface

contains

  !! an implementation of the de Hoog et al. method
  !! assumes proper f(p) have been computed for the p
  !! required for the vector of t passed to this function
  !! -- only one log-cycle of time should be passed at once --
  !! (no error checking done in this regard)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function deHoog_invLap_vect(t,tee,fp,lap) result(ft)
    use constants, only : DP, PI, EYE
    use type_definitions, only : INVLT

    real(DP), intent(in) :: tee              ! scaling factor (previously T=2*tmax, but potentially adjustable)
    real(DP), intent(in), dimension(:) :: t   ! vector of times
    type(INVLT), intent(in) :: lap            ! structure of inputs
    complex(DP), intent(in), dimension(0:2*lap%M) :: fp
    real(DP), dimension(size(t)) :: ft        ! output

    complex(DP), dimension(0:2*lap%M,0:lap%M) :: e
    complex(DP), dimension(0:2*lap%M,1:lap%M) :: q
    complex(DP), dimension(0:2*lap%M) :: d
    complex(DP), dimension(-1:2*lap%M,size(t)) :: A,B
    complex(DP), dimension(size(t)) :: z,brem,rem
    integer :: r, rq, n, max, nt, M
    real(DP) :: gamma

    M = lap%M
    nt = size(t)

    ! there will be problems is fp(:)==0, or any values are NaN
    if(maxval(abs(fp)) > tiny(1.0D0) .and. all(fp == fp)) then

       ! Re(p) -- this is the de Hoog parameter c
       gamma = lap%alpha - log(lap%tol)/(2.0*tee)

       ! initialize Q-D table
       e(0:2*M,0) = cmplx(0.0,0.0,DP)
       q(0,1) = fp(1)/(fp(0)/2.0) ! half first term
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
       d(0) = fp(0)/2.0 ! half first term
       forall(r = 1:M)
          d(2*r-1) = -q(0,r) ! even terms
          d(2*r)   = -e(0,r) ! odd terms
       end forall

       ! seed A and B vectors for recurrence
       A(-1,1:nt) = 0.0
       A(0,1:nt) = d(0)
       B(-1:0,1:nt) = 1.0

       ! base of the power series
       z(1:nt) = exp(EYE*PI*t(:)/tee)

       ! coefficients of Pade approximation
       ! using recurrence for all but last term
       do n = 1,2*M-1
          A(n,:) = A(n-1,:) + d(n)*A(n-2,:)*z(:)
          B(n,:) = B(n-1,:) + d(n)*B(n-2,:)*z(:)
       end do

       ! "improved remainder" to continued fraction
       brem(1:nt) = (1.0 + (d(2*M-1) - d(2*M))*z(:))/2.0
       rem(1:nt) = -brem*(1.0 - sqrt(1.0 + d(2*M)*z(:)/brem**2))

       ! last term of recurrence using new remainder
       A(2*M,:) = A(2*M-1,:) + rem*A(2*M-2,:)
       B(2*M,:) = B(2*M-1,:) + rem*B(2*M-2,:)

       ! diagonal Pade approximation
       ! F=A/B represents accelerated trapezoid rule
       ft(1:nt) =  exp(gamma*t(:))/tee * real(A(2*M,:)/B(2*M,:))

    else  !! entire f(p) vector is zero
       ft = 0.0
!!       write(*,*) 'f(t) not computed: t=',t, &
!!            & ' because max|fp|=', maxval(abs(fp)), ' any(NaN) ? ->',any(fp/=fp)
    end if
  end function deHoog_invLap_vect

  function deHoog_invLap_scalt(t,tee,fp,lap) result(ft)
    use constants, only : DP
    use type_definitions, only : INVLT
    real(DP), intent(in) ::  t, tee
    type(INVLT), intent(in) :: lap
    complex(DP), intent(in), dimension(0:2*lap%M) :: fp
    real(DP) :: ft ! output

    ft = sum(deHoog_invLap_vect([t],tee,fp,lap))
  end function deHoog_invLap_scalt

  function deHoog_invLap_vel(t,tee,fp,lap) result(ft)
    ! for use with velocity results, where second dimension is 2 (x,y)
    use constants, only : DP
    use type_definitions, only : INVLT
    real(DP), dimension(:), intent(in) :: t
    real(DP), intent(in) ::  tee
    complex(DP), intent(in), dimension(:,:) :: fp
    type(INVLT), intent(in) :: lap
    real(DP), dimension(size(t),2) :: ft ! output

    ft(:,1) = deHoog_invLap_vect(t,tee,fp(:,1),lap) ! x vel
    ft(:,2) = deHoog_invLap_vect(t,tee,fp(:,2),lap) ! y vel
  end function deHoog_invLap_vel

  function deHoog_invLap_vel_scalt(t,tee,fp,lap) result(ft)
    ! for use with velocity results, where second dimension is 2 (x,y)
    use constants, only : DP
    use type_definitions, only : INVLT
    real(DP), intent(in) :: t
    real(DP), intent(in) ::  tee
    complex(DP), intent(in), dimension(:,:) :: fp
    type(INVLT), intent(in) :: lap
    real(DP), dimension(2) :: ft ! output

    ft(1) = sum(deHoog_invLap_vect([t],tee,fp(:,1),lap)) ! x vel
    ft(2) = sum(deHoog_invLap_vect([t],tee,fp(:,2),lap)) ! y vel
  end function deHoog_invLap_vel_scalt

  pure function deHoog_pvalues(tee,lap) result(p)
    use constants, only : DP, PI
    use type_definitions, only : INVLT
    type(INVLT), intent(in) :: lap
    real(DP), intent(in) :: tee
    complex(DP), dimension(2*lap%M+1) :: p
    integer :: i

    forall (i=0:2*lap%M)
       p(i+1) = cmplx(lap%alpha - log(lap%tol)/(2.0*tee), PI*i/tee, DP)
    end forall
  end function deHoog_pvalues
end module inverse_Laplace_Transform
