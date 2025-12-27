!
! Copyright (c) 2011-2025 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

! this module contains the basic functions defining the head or flux
! effects of a elliptical element.

module elliptical_elements
  implicit none

  private
  public :: ellipse_match, ellipse_calc, ellipse_deriv, line

  interface ellipse_match
     module procedure ellipse_match_self, ellipse_match_other
  end interface

contains
  function ellipse_match_self(e,p,i,debug) result(r)
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit
    use constants, only : DP
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use type_definitions, only : ellipse, match_result, element
    use mathieu_functions, only : ce, se, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    use utility, only : ynot
    implicit none

    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: i ! indicates which value of p (global state)
    logical, intent(in) :: debug
    type(match_result) :: r

    integer :: j, N, M, loM, hiM, nrows, ncols
    ! mod. radial Mathieu function (K{e,o} or I{e,o})
    complex(DP), allocatable :: RMn(:,:), dRMn(:,:)
    complex(DP), dimension(1:e%M, 0:e%N-1, 0:1) :: cemat
    complex(DP), dimension(1:e%M, 1:e%N-1, 0:1) :: semat
    integer, dimension(0:e%N-1) :: vi

    N = e%N
    M = e%M
    do concurrent (j = 0:N-1)
      vi(j) = j
    end do

    if (e%ibnd == 0) then
       nrows = 2*M
       ncols = 4*N-2
       loM = M+1
       hiM = 2*M
    elseif (e%ibnd == 2) then
       ! specified flux line source has no unknowns
       ! only appears on RHS of other elements
       nrows = 0
       ncols = 0
       loM = 1
       hiM = M
    else
       if (e%calcin .and. (e%ibnd == -1 .or. e%ibnd == 1)) then
          nrows = M
          ncols = 4*N-2
          loM = 1
          hiM = M
          ! TODO: calculating inside a specified head/flux ellipse is broken
          write(stderr,*) 'WARNING: e%ibnd=={-1,1} and calcin are BROKEN?'
          write(stderr,*) 'UNIMPLEMENTED: e%ibnd=={-1,1} and calcin'
          stop 777
       else
          ! ibnd = {-1,1}, but .not. calcin
          nrows = M
          ncols = 2*N-1
          loM = 1
          hiM = M
       end if
    end if

    if (debug) then
       write(stdout,'(A,I0,A,5(I0,1X))') 'ELLIPSE_MATCH_SELF parent: ',&
            & e%parent%id,' (ibnd,nrows,ncols,loM,hiM): ',&
            & e%ibnd,nrows,ncols,loM,hiM
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows))

    ! outside angular modified Mathieu functions (last dimension inside/outside)
    cemat(1:M,0:N-1,0) = ce(e%parent%mat(i), vi(0:N-1), e%Pcm(:))
    semat(1:M,1:N-1,0) = se(e%parent%mat(i), vi(1:N-1), e%Pcm(:))
    if (e%ibnd == 0 .or. (e%calcin .and. (e%ibnd == 1 .or. e%ibnd == -1))) then
       ! inside
       cemat(1:M,0:N-1,1) = ce(e%mat(i), vi(0:N-1), e%Pcm(:))
       semat(1:M,1:N-1,1) = se(e%mat(i), vi(1:N-1), e%Pcm(:))
    end if

    ! setup LHS
    ! matching or specified total head
    if (e%ibnd == 0 .or. e%ibnd == -1) then
       r%LHS(1:M,1:N) =       cemat(:,0:N-1,0)/e%parent%K ! a_n head
       r%LHS(1:M,N+1:2*N-1) = semat(:,1:N-1,0)/e%parent%K ! b_n head

       if (e%ibnd == 0 .or. (e%ibnd == -1 .and. e%calcin)) then
          r%LHS(1:M,2*N:3*N-1) = -cemat(:,0:N-1,1)/e%K ! c_n head
          r%LHS(1:M,3*N:4*N-2) = -semat(:,1:N-1,1)/e%K ! d_n head
       end if
    end if

    ! matching or specified flux (no ibnd==2 case for ellipses)
    if (e%ibnd == 0 .or. e%ibnd == +1) then
       allocate(RMn(0:N-1,0:1),dRMn(0:N-1,0:1))

       ! radial functions last dimension is even/odd
       RMn(0:N-1,0) =   Ke(e%parent%mat(i), vi(0:N-1), e%r) ! even fn
       RMn(1:N-1,1) =   Ko(e%parent%mat(i), vi(1:N-1), e%r) ! odd fn
       RMn(0,1) = 0.0_DP
       dRMn(0:N-1,0) = dKe(e%parent%mat(i), vi(0:N-1), e%r) ! even deriv
       dRMn(1:N-1,1) = dKo(e%parent%mat(i), vi(1:N-1), e%r) ! odd deriv
       dRMn(0,1) = 0.0_DP

       r%LHS(loM:hiM,1:N) =       spread(dRMn(0:N-1,0) / &
            & RMn(0:N-1,0),1,M) * cemat(:,0:N-1,0) ! a_n flux
       r%LHS(loM:hiM,N+1:2*N-1) = spread(dRMn(1:N-1,1) / &
            & RMn(1:N-1,1),1,M) * semat(:,1:N-1,0) ! b_n flu

       if (e%ibnd == 0 .or. (e%ibnd == 1 .and. e%calcin)) then
          RMn(0:N-1,0) =   Ie(e%mat(i), vi(0:N-1), e%r)
          RMn(1:N-1,1) =   Io(e%mat(i), vi(1:N-1), e%r)
          RMn(0,1) = 0.0_DP
          dRMn(0:N-1,0) = dIe(e%mat(i), vi(0:N-1), e%r)
          dRMn(1:N-1,1) = dIo(e%mat(i), vi(1:N-1), e%r)
          dRMn(0,1) = 0.0_DP

          r%LHS(loM:hiM,2*N:3*N-1) = -spread(dRMn(0:N-1,0) / &
               & RMn(0:N-1,0), 1,M) * cemat(:,0:N-1,1) ! c_n flux
          r%LHS(loM:hiM,3*N:4*N-2) = -spread(dRMn(1:N-1,1) / &
               & RMn(1:N-1,1), 1,M) * semat(:,1:N-1,1) ! d_n flux
       end if
       deallocate(RMn,dRMn)
    end if

    ! setup RHS
    select case(e%ibnd)
    case(-1)
       ! put specified head on RHS
       r%RHS(1:M) = timef(p,e%time,.false.)*e%bdryQ
    case(0)
       ! put constant area source term effects on RHS
       ! optional 3rd argument -> kappa**2
       r%RHS(1:M) = -timef(p,e%time,.true.) * e%areaQ * e%Ss / &
            & kappa(p,e%parent,.true.)
       r%RHS(M+1:2*M) = 0.0_DP ! area source has no flux effects
    case(1)
       ! put specified flux effects on RHS
       r%RHS(1:M) = timef(p,e%time,.false.)*e%bdryQ/ynot(e%r,e%f)
    case(2)
       continue ! no ellipse wellbore storage
    end select
  end function ellipse_match_self

  function ellipse_match_other(src,trg,dom,p,idx,debug) result(r)
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit
    use constants, only : DP, CZERO
    use type_definitions, only : ellipse, element, domain
    use type_definitions, only : matching, match_result
    use mathieu_functions, only : ce, se, dce, dse
    use mathieu_functions, only : Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    use utility, only : rotate_vel_mat
    implicit none

    type(ellipse), intent(in) :: src ! source ellipse
    type(matching), intent(in) :: trg ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    logical, intent(in) :: debug
    type(match_result) :: r

    integer :: j, s, t, N, M, loN, hiN, loM, hiM, nrows, ncols
    complex(DP), allocatable :: cemat(:,:), semat(:,:), dcemat(:,:), dsemat(:,:)
    integer, dimension(0:src%N-1) :: vi
    complex(DP), allocatable :: RMn(:,:,:), dRMn(:,:,:), RMn0(:,:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:)
    complex(DP), allocatable :: dPot_dX(:,:), dPot_dY(:,:)
    real(DP), allocatable :: hsq(:,:)
    real(DP) :: K, factor

    N = src%N ! number of coefficients in the source elliptical element
    t = trg%id
    s = src%id
    do concurrent (j = 0:N-1)
      vi(j) = j
    end do

    if (debug) then
       loM = -999
       hiM = -888
       loN = -777
       hiN = -666
    end if

    M = trg%M
    ! target element determines number of rows
    if (trg%ibnd == 0) then
       nrows = 2*M
       loM = M+1
       hiM = 2*M
    elseif (trg%ibnd == 2) then
       nrows = 0
    else
       nrows = M
       loM = 1
       hiM = M
    end if

    ! source element determines number of columns
    if (src%ibnd == 0) then
       ncols = 4*N-2
    elseif (src%ibnd == 2) then
       ! line source (ibound==2) later re-allocated to zero size LHS
       ncols = 2*N-1
    elseif (src%calcin .and. (src%ibnd == -1 .or. src%ibnd == 1)) then
       ! specified flux/head w/ calc inside
       ncols = 4*N-2
    else
       ! specified flux/head w/ no calc inside
       ncols = 2*N-1
    end if

    if (debug) then
       write(stdout,'(A,3(I0,1X),A,4(I0,1X))') &
            & 'ELLIPSE_MATCH_OTHER (parent,s,t): ',&
            & src%parent%id,s,t,' (nrows,ncols,loM,hiM): ',nrows,ncols,loM,hiM
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows))
    r%LHS = CZERO
    r%RHS = CZERO

    if (nrows > 0) then ! is target matching?

       if (dom%inclBg(s,t) .or. dom%InclIn(s,t) .or. dom%InclIn(t,s)) then

          allocate(RMn(M,0:N-1,0:1), RMn0(0:N-1,0:1))
          allocate(cemat(M,0:N-1), semat(M,N-1))

          ! setup LHS
          ! for matching or specified total head target elements
          if (trg%ibnd == 0 .or. trg%ibnd == -1) then

             if (dom%inclBg(s,t) .or. dom%inclin(t,s)) then
                ! can the target element see the outside of the source element?
                ! use exterior angular and radial modified Mathieu functions
                cemat(1:M,0:N-1) = ce(src%parent%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                semat(1:M,1:N-1) = se(src%parent%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                RMn(1:M,0:N-1,0) = Ke(src%parent%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                RMn(1:M,1:N-1,1) = Ko(src%parent%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                RMn(1:M,0,1) = 0.0_DP
                RMn0(0:N-1,0) = Ke(src%parent%mat(idx), vi(0:N-1), src%r)
                RMn0(1:N-1,1) = Ko(src%parent%mat(idx), vi(1:N-1), src%r)
                RMn0(0,1) = 0.0_DP
                ! odd functions not needed for line source,
                ! computed anyway (?)
                K = src%parent%K

                loN = 1
                hiN = 2*N-1

                ! head effects due to ellipse on outside other element
                r%LHS(1:M,loN:loN+N-1) = RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1)/K ! a_n
                r%LHS(1:M,loN+N:hiN)   = RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,1:N-1)/K ! b_n

             else
                ! can target element see the inside of the source element?
                ! i.e., is the source element the parent?
                ! use interior angular and radial modified Mathieu functions
                cemat(1:M,0:N-1) = ce(src%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                semat(1:M,1:N-1) = se(src%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                RMn(1:M,0:N-1,0) = Ie(src%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                RMn(1:M,1:N-1,1) = Io(src%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                RMn(1:M,0,1) = 0.0_DP
                RMn0(0:N-1,0) = Ie(src%mat(idx), vi(0:N-1), src%r)
                RMn0(1:N-1,1) = Io(src%mat(idx), vi(1:N-1), src%r)
                RMn0(0,1) = 0.0_DP
                K = src%K

                if (src%ibnd == 0) then
                   ! is source the inside of matching element?
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! is target inside of specified head/flux element? (no other part)
                   loN = 1
                   hiN = 2*N-1
                end if

                ! head effects due to ellipse on inside other element
                r%LHS(1:M,loN:loN+N-1) = -RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1)/K ! c_n
                r%LHS(1:M,loN+N:hiN)   = -RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,0:N-1)/K ! d_n
             end if

             if (src%ibnd == 2 .and. (dom%inclBg(s,t) .or. dom%inclIn(t,s))) then
                ! save head effects of line source onto RHS
                if (dom%inclBg(s,t)) then
                   ! source well in background of target element
                   factor = -1.0_DP
                else
                   ! source well inside target element
                   factor = 1.0_DP
                end if

                ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ! TODO: copied from circular elements -- check/fix
                ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                ! specified flux (finite-radius well no storage)
                ! save head effects of well onto RHS
                r%RHS(1:M) = factor*line(src,p,idx)*r%LHS(1:M,1)
                r%LHS(1:M,1) = 0.0_DP ! LHS matrix re-allocated to zero size below
             end if

          end if

          ! for matching, specified total flux, or specified elemental flux target element
          if (trg%ibnd == 0 .or. trg%ibnd == +1 .or. trg%ibnd == +2) then
             allocate(dRMn(M,0:N-1,0:1), dcemat(M,0:N-1), dsemat(M,N-1), &
                  & dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), dPot_dX(M,2*N-1),dPot_dY(M,2*N-1))

             ! flux effects of source ellipse on target element
             if (dom%inclBg(s,t) .or. dom%inclIn(t,s)) then
                ! use exterior angular and radial modified Mathieu functions
                if (.not. trg%ibnd == 0) then
                   cemat(1:M,0:N-1) = ce(src%parent%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                   RMn(1:M,0:N-1,0) = Ke(src%parent%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                   RMn0(0:N-1,0)    = Ke(src%parent%mat(idx), vi(0:N-1), src%r)
                   semat(1:M,1:N-1) = se(src%parent%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                   RMn(1:M,1:N-1,1) = Ko(src%parent%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                   RMn(1:M,0,1) = 0.0_DP
                   RMn0(1:N-1,1) = Ko(src%parent%mat(idx), vi(1:N-1), src%r)
                   RMn0(0,1) = 0.0_DP
                   K = src%parent%K
                else
                   K = -999  ! compiler complained this might be used uninitialized
                end if
                dcemat(1:M,0:N-1) = dce(src%parent%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                dRMn(1:M,0:N-1,0) = dKe(src%parent%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                dsemat(1:M,1:N-1) = dse(src%parent%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                dRMn(1:M,1:N-1,1) = dKo(src%parent%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                dRMn(1:M,0,1) = 0.0_DP

                loN = 1
                hiN = 2*N-1

             else
                ! use interior angular and radial modified Mathieu functions
                if (.not. trg%ibnd == 0) then
                   cemat(1:M,0:N-1) = ce(src%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                   semat(1:M,1:N-1) = se(src%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                   RMn(1:M,0:N-1,0) = Ie(src%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                   RMn(1:M,1:N-1,1) = Io(src%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                   RMn(1:M,0,1) = 0.0_DP
                   RMn0(0:N-1,0) = -Ie(src%mat(idx), vi(0:N-1), src%r) ! apply in/out-side sign here
                   RMn0(1:N-1,1) = -Io(src%mat(idx), vi(1:N-1), src%r)
                   RMn0(0,1) = 0.0_DP
                end if
                dcemat(1:M,0:N-1) = dce(src%mat(idx), vi(0:N-1), src%G(t)%Pgm(:))
                dsemat(1:M,1:N-1) = dse(src%mat(idx), vi(1:N-1), src%G(t)%Pgm(:))
                dRMn(1:M,0:N-1,0) = dIe(src%mat(idx), vi(0:N-1), src%G(t)%Rgm(:))
                dRMn(1:M,1:N-1,1) = dIo(src%mat(idx), vi(1:N-1), src%G(t)%Rgm(:))
                dRMn(1:M,0,1) = 0.0_DP
                K = src%K

                if (src%ibnd == 0) then
                   ! inside matching element (last half)
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! inside specified flux element (only part)
                   loN = 1
                   hiN = 2*N-1
                end if

             end if

             ! derivative with respect to radius of source element
             dPot_dR(1:M,1:N) =       dRMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1)
             dPot_dR(1:M,N+1:2*N-1) = dRMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,1:N-1)

             ! derivative with respect to angle of source element
             dPot_dP(1:M,1:N) =       RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*dcemat(:,0:N-1)
             dPot_dP(1:M,N+1:2*N-1) = RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*dsemat(:,1:N-1)

             ! project these from elliptical onto Cartesian coordinates
             allocate(hsq(size(src%G(t)%Rgm),2*N-1))

             ! squared metric factor -- less a common f
             hsq = spread(src%f*0.5_DP*(cosh(2.0_DP*src%G(t)%Rgm) - cos(2.0_DP*src%G(t)%Pgm)),2,2*N-1)
             dPot_dX = (dPot_dR*spread(sinh(src%G(t)%Rgm)*cos(src%G(t)%Pgm),2,2*N-1) - &
                      & dPot_dP*spread(cosh(src%G(t)%Rgm)*sin(src%G(t)%Pgm),2,2*N-1))/hsq
             dPot_dY = (dPot_dR*spread(cosh(src%G(t)%Rgm)*sin(src%G(t)%Pgm),2,2*N-1) + &
                      & dPot_dP*spread(sinh(src%G(t)%Rgm)*cos(src%G(t)%Pgm),2,2*N-1))/hsq

             ! rotate to compensate for potentially arbitrary source ellipse
             call rotate_vel_mat(dPot_dX,dPot_dY,src%theta)

             deallocate(hsq)

             ! project from Cartesian to radial coordinate of target element
             if (trg%id <= dom%num(1)) then
                ! other element is a circle
                if (trg%ibnd == 2) then
                   ! other element is a well with wellbore storage (Type III BC)
                   r%LHS(1:M,loN:loN+N-1) = RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1)/K
                   r%LHS(1:M,loN+N:hiN) =   RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,1:N-1)/K

                   ! head effects of source ellipse
                   r%LHS(1:M,loN:hiN) = -(trg%r*p/trg%parent%T)*r%LHS(:,loN:hiN)

                   ! radial flux effects of element
                   r%LHS(1:M,loN:hiN) = r%LHS + (2.0_DP + trg%r**2*trg%dskin*p/trg%parent%T)* &
                        & (dPot_dX*spread(cos(trg%Pcm),2,2*N-1) + &
                        &  dPot_dY*spread(sin(trg%Pcm),2,2*N-1))
                else
                   ! other element is a standard circle
                   r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(cos(trg%Pcm),2,2*N-1) + &
                        & dPot_dY*spread(sin(trg%Pcm),2,2*N-1)
                end if
             else
                ! rotate to compensate for arbitrary target ellipse
                call rotate_vel_mat(dPot_dX,dPot_dY,-trg%theta)

                ! other element is a different ellipse
                r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(trg%f*sinh(trg%r)*cos(trg%Pcm(:)),2,2*N-1) + &
                     & dPot_dY*spread(trg%f*cosh(trg%r)*sin(trg%Pcm(:)),2,2*N-1)

             end if
             deallocate(dRMn,dcemat,dsemat,dPot_dR,dPot_dP,dPot_dX,dPot_dY)
          end if
          deallocate(RMn,RMn0,cemat,semat)

       end if
    end if
    if (src%ibnd == 2) then
       if (nrows > 0) then
          ! sum line source effects and move to RHS, re-setting LHS to 0 unknowns
          ! only uses even-order even coefficients (~1/4 of 2N-1)
          r%RHS(1:hiM) = -sum(spread(line(src,p,idx),1,hiM)*r%LHS(1:hiM,1:N:2))
       end if
       deallocate(r%LHS)
       allocate(r%LHS(nrows,0))
    end if

  end function ellipse_match_other

  function ellipse_calc(p,e,lo,hi,Rgp,Pgp,inside) result(H)
    use constants, only : DP
    use type_definitions, only : ellipse
    use time_mod, only : timef
    use mathieu_functions, only : Ke,Ko, Ie,Io, ce,se

    complex(DP), dimension(:), intent(in) :: p
    type(ellipse), intent(in) :: e
    integer, intent(in) :: lo,hi
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(size(p,1)) :: H

    integer, dimension(0:e%N-1) :: vi
    complex(DP), dimension(size(p,1),0:e%N-1) :: aa,bb
    complex(DP), dimension(size(p,1),0:e%N-1,0:1) :: RMRgp,RMR0,AM
    integer :: n0, np, i, j, N

    N = e%N
    np = size(p,1)
    do concurrent (i = 0:N-1)
      vi(i) = i
    end do

    if (inside) then
       if (e%match) then
          n0 = 2*N ! inside of matching ellipse
       else
          n0 = 1   ! inside of specified boundary ellipse
       end if
       do i = 1,np
          j = lo+i-1
          RMRgp(i,0:N-1,0) = Ie(e%mat(j),vi(0:N-1),Rgp)
          RMRgp(i,1:N-1,1) = Io(e%mat(j),vi(1:N-1),Rgp)

          RMR0(i,0:N-1,0) =  Ie(e%mat(j),vi(0:N-1),e%r)
          RMR0(i,1:N-1,1) =  Io(e%mat(j),vi(1:N-1),e%r)

          AM(i,0:N-1,0) =    ce(e%mat(j),vi(0:N-1),Pgp)
          AM(i,1:N-1,1) =    se(e%mat(j),vi(1:N-1),Pgp)
       end do
    else
       n0 = 1
       do i = 1,np
          j = lo+i-1
          RMRgp(i,0:N-1,0) = Ke(e%parent%mat(j),vi(0:N-1),Rgp)
          RMRgp(i,1:N-1,1) = Ko(e%parent%mat(j),vi(1:N-1),Rgp)

          RMR0(i,0:N-1,0) =  Ke(e%parent%mat(j),vi(0:N-1),e%r)
          RMR0(i,1:N-1,1) =  Ko(e%parent%mat(j),vi(1:N-1),e%r)

          AM(i,0:N-1,0) =    ce(e%parent%mat(j),vi(0:N-1),Pgp)
          AM(i,1:N-1,1) =    se(e%parent%mat(j),vi(1:N-1),Pgp)
       end do
    end if

    aa(1:np,0:N-1) = e%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,1:N-1) = e%coeff(lo:hi,n0+N:n0+2*N-2)
    H(1:np) = sum(RMRgp(1:np,0:N-1,0)/RMR0(1:np,0:N-1,0)*aa(1:np,0:N-1)*AM(1:np,0:N-1,0), 2) + &
            & sum(RMRgp(1:np,1:N-1,1)/RMR0(1:np,1:N-1,1)*bb(1:np,1:N-1)*AM(1:np,1:N-1,1), 2)

  end function ellipse_calc

  function ellipse_deriv(p,e,lo,hi,Rgp,Pgp,inside) result(dH)
    use constants, only : DP
    use type_definitions, only : ellipse
    use time_mod, only : timef
    use mathieu_functions, only : Ke,Ko, Ie,Io, ce,se, dKe,dKo, dIe,dIo, dce,dse

    complex(DP), dimension(:), intent(in) :: p
    type(ellipse), intent(in) :: e
    integer, intent(in) :: lo,hi
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(size(p,1),2) :: dH ! dPot_dEta, dPot_dPsi

    integer, dimension(0:e%N-1) :: vi
    complex(DP), dimension(size(p,1),0:e%N-1) :: aa,bb
    complex(DP), dimension(size(p,1),0:e%N-1,0:1) :: RMRgp,RMR0,AM,dRMRgp,dAM
    integer :: n0, np, i, j, N

    N = e%N
    np = size(p,1)
    do concurrent (i = 0:N-1)
      vi(i) = i
    end do

    if (inside) then
       if (e%match) then
          n0 = 2*N ! inside of matching ellipse
       else
          n0 = 1   ! inside of specified boundary ellipse
       end if
       do i = 1,np
          j = lo+i-1
          ! radial Mathieu functions -> to calculation point
          RMRgp(i,0:N-1,0) =   Ie(e%mat(j),vi(0:N-1),Rgp)
          RMRgp(i,1:N-1,1) =   Io(e%mat(j),vi(1:N-1),Rgp)
          dRMRgp(i,0:N-1,0) = dIe(e%mat(j),vi(0:N-1),Rgp)
          dRMRgp(i,1:N-1,1) = dIo(e%mat(j),vi(1:N-1),Rgp)
          ! radial Mathieu functions -> to boundary of ellipse
          RMR0(i,0:N-1,0) =    Ie(e%mat(j),vi(0:N-1),e%r)
          RMR0(i,1:N-1,1) =    Io(e%mat(j),vi(1:N-1),e%r)
          ! angular Mathieu functions -> to calculation point
          AM(i,0:N-1,0) =      ce(e%mat(j),vi(0:N-1),Pgp)
          AM(i,1:N-1,1) =      se(e%mat(j),vi(1:N-1),Pgp)
          dAM(i,0:N-1,0) =    dce(e%mat(j),vi(0:N-1),Pgp)
          dAM(i,1:N-1,1) =    dse(e%mat(j),vi(1:N-1),Pgp)
       end do
    else
       n0 = 1
       do i = 1,np
          j = lo+i-1
          RMRgp(i,0:N-1,0) =   Ke(e%parent%mat(j),vi(0:N-1),Rgp)
          RMRgp(i,1:N-1,1) =   Ko(e%parent%mat(j),vi(1:N-1),Rgp)
          dRMRgp(i,0:N-1,0) = dKe(e%parent%mat(j),vi(0:N-1),Rgp)
          dRMRgp(i,1:N-1,1) = dKo(e%parent%mat(j),vi(1:N-1),Rgp)

          RMR0(i,0:N-1,0) =    Ke(e%parent%mat(j),vi(0:N-1),e%r)
          RMR0(i,1:N-1,1) =    Ko(e%parent%mat(j),vi(1:N-1),e%r)

          AM(i,0:N-1,0) =      ce(e%parent%mat(j),vi(0:N-1),Pgp)
          AM(i,1:N-1,1) =      se(e%parent%mat(j),vi(1:N-1),Pgp)
          dAM(i,0:N-1,0) =    dce(e%parent%mat(j),vi(0:N-1),Pgp)
          dAM(i,1:N-1,1) =    dse(e%parent%mat(j),vi(1:N-1),Pgp)
       end do
    end if

    aa(1:np,0:N-1) = e%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,1:N-1) = e%coeff(lo:hi,n0+N:n0+2*N-2)
    dH(1:np,1) = sum(dRMRgp(1:np,0:N-1,0)/RMR0(1:np,0:N-1,0)*aa(1:np,0:N-1) *AM(1:np,0:N-1,0), 2) + &
               & sum(dRMRgp(1:np,1:N-1,1)/RMR0(1:np,1:N-1,1)*bb(1:np,1:N-1) *AM(1:np,1:N-1,1), 2)
    dH(1:np,2) = sum(RMRgp(1:np,0:N-1,0) /RMR0(1:np,0:N-1,0)*aa(1:np,0:N-1)*dAM(1:np,0:N-1,0), 2) + &
               & sum(RMRgp(1:np,1:N-1,1) /RMR0(1:np,1:N-1,1)*bb(1:np,1:N-1)*dAM(1:np,1:N-1,1), 2)

  end function ellipse_deriv

  function line(e,p,idx) result(a2n)
    ! this function returns the coefficients for a specified-flux line source
    use constants, only : DP, TWOPI
    use time_mod, only : timef
    use type_definitions, only : ellipse
    use mathieu_functions, only : Ke,dKe

    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    ! only even coefficients of even order
    complex(DP), dimension(ceiling(e%N*0.5_DP)) :: a2n
    real(DP), dimension(0:max(e%ms-1,e%N)) :: vs
    real(DP), dimension(1:e%ms,ceiling(e%N*0.5_DP)) :: arg
    integer, dimension(0:max(e%ms-1,e%N)) :: vi
    integer :: i, N, MS, nmax

    N = e%N
    MS = e%ms
    nmax = ceiling(e%N*0.5_DP)
    do concurrent (i = 0:max(MS,N)-1)
      vi(i) = i ! integer vector
      if (mod(vi(i),2) == 0) then
        vs = 1.0_DP
      else
        vs = -1.0_DP
      end if
    end do

    arg(1:MS,1:nmax) = spread(vs(0:MS-1)/real(1.0_DP-(2.0_DP*vi(0:MS-1))**2,DP),2,nmax)

    ! factor of 4 different from Kuhlman & Neuman (J. Eng. Mathematics) paper
    a2n(1:nmax) = timef(p,e%time,.false.)*e%bdryQ/TWOPI* &
            & (-vs(0:N-1:2))*sum(arg*conjg(e%parent%mat(idx)%A(1:MS,0:nmax-1,0)),dim=1)

  end function line

end module elliptical_elements
