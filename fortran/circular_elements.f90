!
! Copyright (c) 2011-2014 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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
! effects of a circular element.

module circular_elements
  implicit none

  private
  public :: circle_match, circle_calc, circle_deriv, well

  interface circle_match
     module procedure circle_match_self, circle_match_other
  end interface

contains
  function circle_match_self(c,p,debug) result(r)
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use utility, only : outer
    use type_definitions, only : circle, match_result
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    logical, intent(in) :: debug
    type(match_result) :: r

    integer :: j, N, M, loM, hiM, nrows, ncols
    complex(DP), allocatable :: Bn(:), dBn(:) ! mod. Bessel function (K or I)
    complex(DP) :: kap
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi

    N = c%N
    M = c%M
    forall (j = 0:N-1) vi(j) = real(j,DP)
    
    if (c%ibnd == 0) then
       ! matching
       nrows = 2*M
       ncols = 4*N-2
       ! loM:hiM is index range for flux matching portion, beyond head-matching part
       loM = M+1
       hiM = 2*M
    elseif (c%ibnd == 2) then
       ! simple well has no unknowns
       ! only appears on RHS of other elements
       loM = 1
       hiM = M
       if (c%storIn) then
          nrows = M ! one unknown, one? matching point
          ncols = 1
       else
          nrows = 0
          ncols = 0
       end if
    else
       if (c%calcin .and. (c%ibnd == -1 .or. c%ibnd == 1)) then
          ! if computing inside a specified head/flux element
          ! then need to compute solutions on both sides (even if
          ! only computing on the inside -- just to be safe)
          nrows = M
          ncols = 4*N-2
          ! only flux or head matching, so loM:him = 1:M
          loM = 1
          hiM = M
       else
          ! other cases (no calc inside)
          nrows = M
          ncols = 2*N-1
          ! here only head or flux matching, so loM:him = 1:M
          loM = 1
          hiM = M
       end if
    end if

    if (debug) then
       print '(2(A,I0),A,5(I0,1X))', 'CIRCLE_MATCH_SELF ',c%id,', parent: ',&
            & c%parent%id,' (ibnd,nrows,ncols,loM,hiM): ',&
            & c%ibnd,nrows,ncols,loM,hiM
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows))

    cmat(1:M,0:N-1) = cos(outer(c%Pcm(:),vi(0:N-1)))
    smat(1:M,1:N-1) = sin(outer(c%Pcm(:),vi(1:N-1)))

    ! setup LHS
    ! matching or specified total head
    if (c%ibnd == 0 .or. c%ibnd == -1) then

       r%LHS(1:M,1:N) =       cmat(:,0:N-1)/c%parent%K ! a_n head
       r%LHS(1:M,N+1:2*N-1) = smat(:,1:N-1)/c%parent%K ! b_n head

       if (c%ibnd == 0 .or. (c%ibnd == -1 .and. c%calcin)) then
          r%LHS(1:M,2*N:3*N-1) = -cmat(:,0:N-1)/c%K ! c_n head
          r%LHS(1:M,3*N:4*N-2) = -smat(:,1:N-1)/c%K ! d_n head
       end if
    end if

    ! matching or specified total flux
    if (c%ibnd == 0 .or. c%ibnd == 1 .or. (c%ibnd == 2 .and. c%storIn)) then
       allocate(Bn(0:N-1), dBn(0:N-1))
       kap = kappa(p,c%parent)
       call dBK(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
       dBn(0:N-1) = kap*dBn(0:N-1)

       r%LHS(loM:hiM,1:N) =       spread(dBn(0:N-1)/Bn(0:N-1),1,M)*cmat(:,0:N-1) ! a_n flux
       r%LHS(loM:hiM,N+1:2*N-1) = spread(dBn(1:N-1)/Bn(1:N-1),1,M)*smat(:,1:N-1) ! b_n flux

       if (c%ibnd == 0 .or. (c%ibnd == 1 .and. c%calcin)) then
          kap = kappa(p,c%element)
          call dBI(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
          dBn(0:N-1) = kap*dBn(0:N-1)

          r%LHS(loM:hiM,2*N:3*N-1) = -spread(dBn(0:N-1)/Bn(0:N-1),1,M)*cmat(:,0:N-1) ! c_n flux
          r%LHS(loM:hiM,3*N:4*N-2) = -spread(dBn(1:N-1)/Bn(1:N-1),1,M)*smat(:,1:N-1) ! d_n flux
       end if
       deallocate(Bn,dBn)
    end if

    ! setup RHS
    select case(c%ibnd)
    case(-1)
       ! put specified head on RHS
       r%RHS(1:M) = timef(p,c%time,.false.)*c%bdryQ
    case(0)
       ! put constant area source term effects (from inside the element) on RHS
       ! TODO : handle area source in background
       r%RHS(1:M) = -timef(p,c%time,.true.)*c%areaQ*c%Ss/kappa(p,c%element)**2
       r%RHS(M+1:2*M) = 0.0 ! constant area source has no flux effects
    case(1)
       ! put specified flux effects on RHS
       ! TODO : check addition of aquifer thickness to denominator
       r%RHS(1:M) = timef(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*c%b)
    case(2)
       if (c%StorIn) then
          ! effects of wellbore storage and skin on finite-radius well
          ! effects of other elements on this one show up in off-diagonals
          r%LHS(1:M,1) = storwell(c,p)*r%LHS(1:M,1)
          r%RHS(1:M) = timef(p,c%time,.false.)*c%bdryQ/(PI*c%r*c%parent%T)
       else
          continue ! no wellbore storage; nothing to do, since matrix is zero-sized
       end if
    end select
  end function circle_match_self

  function circle_match_other(src,trg,dom,p,debug) result(r)
    use constants, only : DP
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use utility, only : outer, rotate_vel_mat
    use type_definitions, only : circle, domain, matching, match_result
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(circle), intent(in) :: src   ! source circle
    type(matching), intent(in) :: trg ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    logical, intent(in) :: debug
    type(match_result) :: r

    integer :: j, s, t, N, M, loM, hiM, loN, hiN, nrows, ncols
    real(DP) :: K, factor
    real(DP), allocatable :: cmat(:,:), smat(:,:)
    real(DP), dimension(0:src%N-1) :: vi
    complex(DP), allocatable :: Bn(:,:), dBn(:,:), Bn0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    complex(DP) :: kap

    N = src%N ! number of coefficients in the source circular element
    t = trg%id
    s = src%id
    forall (j = 0:N-1) vi(j) = real(j,DP)

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
       if (trg%storIn) then
          nrows = M
          loM = 1
          hiM = M
       else
          nrows = 0
       end if
    else
       nrows = M
       loM = 1
       hiM = M
    end if

    ! source element determines number of columns
    if (src%ibnd == 0) then
       ncols = 4*N-2
    elseif (src%ibnd == 2) then
       if(src%storIn) then
          ncols = 1
       else
          ! compute effects due to well for RHS
          ncols = 1  ! reset to zero at end of routine
       end if
    elseif (src%calcin .and. (src%ibnd == -1 .or. src%ibnd == 1)) then
       ! specified flux/head with calc inside
       ncols = 4*N-2
    else
       ! specified flux/head w/o calc inside
       ncols = 2*N-1
    end if

    if (debug) then
       print '(A,3(I0,1X),A,4(I0,1X))', 'CIRCLE_MATCH_OTHER (parent,s,t): ',&
            & src%parent%id,s,t,' (nrows,ncols,loM,hiM): ',nrows,ncols,loM,hiM
    end if
    
    allocate(r%LHS(nrows,ncols), r%RHS(nrows))
    r%LHS = cmplx(0,0,DP)
    r%RHS = cmplx(0,0,DP)

    if (nrows > 0) then

       if (dom%inclBg(s,t) .or. dom%InclIn(s,t) .or. dom%InclIn(t,s)) then  

          allocate(Bn(M,0:N-1),Bn0(0:N-1),cmat(M,0:N-1),smat(M,N-1))

          cmat(1:M,0:N-1) = cos(outer(src%G(t)%Pgm(:),vi(0:N-1)))
          smat(1:M,1:N-1) = sin(outer(src%G(t)%Pgm(:),vi(1:N-1)))

          ! setup LHS
          ! for matching or specified total head target elements
          if (trg%ibnd == 0 .or. trg%ibnd == -1) then

             if (dom%inclBg(s,t) .or. dom%inclIn(t,s)) then   
                ! can the target element "see" the outside of the source element?
                ! use exterior Bessel functions (Kn)

                kap = kappa(p,src%parent)
                Bn(1:M,0:N-1) = bK(kap*src%G(t)%Rgm(:),N)
                Bn0(0:N-1) =    bK(kap*src%r,N)

                ! outside is always in 1:2*N-1 slot
                loN = 1
                hiN = 2*N-1

                if (debug) then
                   print '(A,2(I0,1X))', 'CIRCLE_MATCH_OTHER HEAD OUTSIDE: (loN,hiN): ', loN,hiN
                end if
                
                ! head effects on target element
                r%LHS(1:M,loN:loN+N-1) = Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1)/src%parent%K ! a_n
                r%LHS(1:M,loN+N:hiN)   = Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/src%parent%K ! b_n

             else
                ! can target element "see" the inside of the source element?
                ! i.e., is the source element the parent?
                ! use interior Bessel functions (In)

                kap = kappa(p,src%element)
                Bn(1:M,0:N-1) = bI(kap*src%G(t)%Rgm(:),N)
                Bn0(0:N-1) =    bI(kap*src%r,N)

                if (src%ibnd == 0) then
                   ! is source the inside of matching element?
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! is source inside of specified head/flux element? (no other previous part)
                   loN = 1
                   hiN = 2*N-1
                end if

                if (debug) then
                   print '(A,2(I0,1X))', 'CIRCLE_MATCH_OTHER HEAD INSIDE: (loN,hiN): ', loN,hiN
                end if

                ! head effects on other element
                r%LHS(1:M,loN:loN+N-1) = -Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1)/src%K ! c_n
                r%LHS(1:M,loN+N:hiN)   = -Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/src%K ! d_n
             end if

             if (src%ibnd == 2 .and. (dom%inclBg(s,t) .or. dom%inclIn(t,s))) then
                if (src%StorIn) then

                   ! wellbore storage and skin from finite-radius well
                   r%LHS(1:M,1) = storwell(src,p)*r%LHS(1:M,1)
                   r%RHS(1:M) = 0.0
                else
                   ! save head effects of well onto RHS
                   if (dom%inclBg(s,t)) then
                      ! source well in background of target element
                      factor = -1.0_DP
                   else
                      ! source well inside target element
                      factor = 1.0_DP
                   end if

                   ! specified flux (finite-radius well no storage)
                   ! save head effects of well onto RHS
                   r%RHS(1:M) = factor*well(src,p)*r%LHS(1:M,1)
                   r%LHS(1:M,1) = 0.0 ! LHS matrix re-allocated to zero size below
                end if
             end if
          end if

          ! for matching, specified total flux, or well with wellbore storage target element
          if (trg%ibnd == 0 .or. trg%ibnd == +1 .or. (trg%ibnd == +2 .and. trg%storIn)) then
             allocate(dBn(M,0:N-1), dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), &
                  & dPot_dX(M,2*N-1), dPot_dY(M,2*N-1))

             ! flux effects of source circle on target element
             if (dom%inclBg(s,t) .or. dom%inclIn(t,s)) then 
                ! use exterior Bessel functions (Kn)

                kap = kappa(p,src%parent)
                call dBK(kap*src%G(t)%Rgm(:),N,Bn(:,0:N-1),dBn(:,0:N-1))
                dBn(1:M,0:N-1) = kap*dBn(:,0:N-1)
                Bn0(0:N-1) = bK(kap*src%r,N)
                K = src%parent%K

                ! outside part is always 1:2*N-1
                loN = 1
                hiN = 2*N-1

                if (debug) then
                   print '(A,2(I0,1X))', 'CIRCLE_MATCH_OTHER FLUX OUTSIDE: (loN,hiN): ', loN,hiN
                end if

             else
                ! use interior Bessel functions (In)
                kap = kappa(p,src%element)
                call dBI(kap*src%G(t)%Rgm(:),N,Bn(:,0:N-1),dBn(:,0:N-1))
                dBn(1:M,0:N-1) = -kap*dBn(:,0:N-1) ! apply in/out-side sign here
                Bn0(0:N-1) = bI(kap*src%r,N)
                K = src%K

                if (src%ibnd == 0) then
                   ! inside of matching element
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! inside of specified flux element
                   loN = 1
                   hiN = 2*N-1
                end if
                
                if (debug) then
                   print '(A,2(I0,1X))', 'CIRCLE_MATCH_OTHER FLUX INSIDE: (loN,hiN): ', loN,hiN
                end if
             end if

             ! derivative wrt radius of source element
             dPot_dR(1:M,1:N) =       dBn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1)
             dPot_dR(1:M,N+1:2*N-1) = dBn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)

             if (trg%ibnd == 2 .and. dom%inclBg(s,t)) then
                ! wellbore storage and skin from finite-radius well
                dPot_dR(1:M,1) = storwell(src,p)*dPot_dR(:,1)
                dPot_dP(1:M,1:2*N-1) = 0.0 ! wells have angular symmetry
             else
                ! derivative wrt angle of source element for more general circular elements
                dPot_dP(1:M,1) = 0.0
                dPot_dP(1:M,2:N) =      -Bn(:,1:N-1)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*smat(:,1:N-1)
                dPot_dP(1:M,N+1:2*N-1) = Bn(:,1:N-1)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*cmat(:,1:N-1)
             end if

             ! project these from cylindrical onto Cartesian coordinates
             dPot_dX(1:M,1:2*N-1) = dPot_dR*spread(cos(src%G(t)%Pgm),2,2*N-1) - &
                                  & dPot_dP*spread(sin(src%G(t)%Pgm)/src%G(t)%Rgm,2,2*N-1)
             dPot_dY(1:M,1:2*N-1) = dPot_dR*spread(sin(src%G(t)%Pgm),2,2*N-1) + &
                                  & dPot_dP*spread(cos(src%G(t)%Pgm)/src%G(t)%Rgm,2,2*N-1)

             ! project from Cartesian to "radial" coordinate of target element
             if (trg%id <= dom%num(1)) then
                ! other element is a circle
                if (trg%ibnd == 2) then
                   ! other element is a well with wellbore storage (Type III BC)
                   ! need head effects too
                   r%LHS(1:M,loN:loN+N-1) = Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat/K
                   r%LHS(1:M,loN+N:hiN) =   Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/K

                   ! head effects of element
                   r%LHS(1:M,loN:hiN) = -(trg%r*p/trg%parent%T)*r%LHS(1:M,loN:hiN)

                   ! radial flux effects of element
                   r%LHS(1:M,loN:hiN) = (r%LHS + (2.0 + trg%r**2*trg%dskin*p/trg%parent%T)* &
                        & (dPot_dX*spread(cos(trg%Pcm),2,2*N-1) + &
                        &  dPot_dY*spread(sin(trg%Pcm),2,2*N-1)))
                else
                   ! other element is a 'normal' circular element without wellbore storage
                   r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(cos(trg%Pcm),2,2*N-1) + &
                                          & dPot_dY*spread(sin(trg%Pcm),2,2*N-1)
                end if
             else
                ! rotate to allow for arbitrary oriented ellipse
                call rotate_vel_mat(dPot_dX,dPot_dY,-trg%theta)

                ! other element is an ellipse
                r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(trg%f*sinh(trg%r)*cos(trg%Pcm(:)),2,2*N-1) + &
                                       & dPot_dY*spread(trg%f*cosh(trg%r)*sin(trg%Pcm(:)),2,2*N-1)
             end if
             deallocate(dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY)
          end if
          deallocate(Bn,Bn0,cmat,smat)
       end if
    endif

    if (src%ibnd == 2 .and. (.not. src%storin)) then
       if (nrows > 0) then

          ! save flux effects of well onto RHS
          if (dom%inclBg(s,t)) then
             ! source well in background of target element
             factor = -1.0_DP
          else
             ! source well inside target element
             factor = 1.0_DP
          end if

          r%RHS(loM:hiM) = factor*well(src,p)*r%LHS(loM:hiM,1)
       end if

       if (trg%ibnd == 0) then
          hiM = 2*M
       else
          hiM = M
       end if

       deallocate(r%LHS)
       allocate(r%LHS(nrows,0))
    end if

  end function circle_match_other

  function well(c,p) result(a0)
    ! this function returns the a_0 coefficient for a simple "well"
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use type_definitions, only : circle
    use bessel_functions, only : bK
    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    complex(DP) :: a0, kap
    complex(DP), dimension(0:1) ::Kn

    kap = kappa(p,c%parent)
    Kn(0:1) = bK(kap*c%r,2)
    ! TODO: should this have a factor of "b" (formation thickness) in the denominator?
    a0 = Kn(0)*timef(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn(1)*kap)
  end function well

  function storwell(c,p) result(a0)
    ! this function returns the a_0 coefficient for a
    ! well with wellbore storage and skin
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use type_definitions, only : circle
    use bessel_functions, only : bK

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    complex(DP) :: a0, kap
    complex(DP), dimension(0:1) :: Kn

    kap = kappa(p,c%parent)
    Kn(0:1) = bK(kap*c%r,2)
    ! TODO : should this have a factor of "b" in the denominator?
    ! TODO : off by a factor of 1/kappa?
    a0 = -Kn(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
               & (Kn(0)*c%r*p)/(2.0*PI*c%r*kap*Kn(1)*c%parent%T))
  end function storwell

  function circle_calc(p,c,lo,hi,Rgp,Pgp,inside) result(H)
    use constants, only : DP
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use type_definitions, only : circle
    use bessel_functions, only : bK, bI

    complex(DP), dimension(:), intent(in) :: p
    type(circle), intent(in) :: c
    integer, intent(in) :: lo,hi
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(size(p,1)) :: H

    real(DP), dimension(0:c%N-1) :: vr
    complex(DP), dimension(size(p,1),0:c%N-1) :: aa,bb,BRgp,BR0
    integer :: n0, np, i, N
    complex(DP), dimension(size(p,1)) :: kap

    N = c%N
    np = size(p,1)
    forall (i = 0:N-1) vr(i) = real(i,DP)

    if (inside) then
       if (c%ibnd == 0) then
          n0 = 2*N ! inside of matching circle
       else
          n0 = 1   ! inside of specified {head,flux} boundary circle
       end if
       kap(1:np) = kappa(p(:),c%element)
       BRgp(1:np,0:N-1) = bI(Rgp*kap(:),N)
       BR0(1:np,0:N-1) =  bI(c%r*kap(:),N)

    else
       n0 = 1
       kap(1:np) = kappa(p(:),c%parent)
       BRgp(1:np,0:N-1) = bK(Rgp*kap(:),N)
       BR0(1:np,0:N-1) =  bK(c%r*kap(:),N)
    end if

    ! if inside
    ! c_n for specified 1:N,      for matching 2N:3N-1
    ! d_n for specified N+1:2N-1, for matching 3N:4N-2

    ! if outside >> a_n is 1:N, b_n is N+1:2N-1

    aa(1:np,0:N-1) = c%coeff(lo:hi,n0:n0+N-1) ! a_n or c_n
    bb(1:np,0) = 0.0 ! make odd/even conformable
    bb(1:np,1:N-1) = c%coeff(lo:hi,n0+N:n0+2*N-2) ! b_n or d_n

    H(1:np) = sum(BRgp(1:np,0:N-1)/BR0(1:np,0:N-1)* &
         & ( aa(1:np,0:N-1)*spread(cos(vr(0:N-1)*Pgp),1,np) + &
         &   bb(1:np,0:N-1)*spread(sin(vr(0:N-1)*Pgp),1,np) ),dim=2)

  end function circle_calc

  function circle_deriv(p,c,lo,hi,Rgp,Pgp,inside) result(dH)
    use constants, only : DP
    use kappa_mod, only : kappa
    use time_mod, only : timef
    use type_definitions, only : circle
    use bessel_functions, only : bK, bI, dbk, dbi

    complex(DP), dimension(:), intent(in) :: p
    type(circle), intent(in) :: c
    integer, intent(in) :: lo,hi
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(size(p,1),2) :: dH ! dPot_dR, dPot_dTheta

    real(DP), dimension(0:c%N-1) :: vr
    complex(DP), dimension(size(p,1),0:c%N-1) :: aa,bb,BRgp,BR0,dBRgp
    integer :: n0, np, i, N
    complex(DP), dimension(size(p,1)) :: kap

    N = c%N
    np = size(p,1)

    forall (i = 0:N-1) vr(i) = real(i,DP)

    if (inside) then
       if (c%match) then
          n0 = 2*N ! inside of matching circle
       else
          n0 = 1   ! inside of specified boundary circle
       end if
       kap(1:np) = kappa(p(:),c%element)
       call dBI(Rgp*kap(:),N,BRgp(1:np,0:N-1),dBRgp(1:np,0:N-1))
       dBRgp(1:np,0:N-1) = spread(kap(:),2,N)*dBRgp(:,:)
       BR0(1:np,0:N-1) = bI(c%r*kap(:),N)
    else
       n0 = 1
       kap(1:np) = kappa(p(:),c%parent)
       call dBK(Rgp*kap(:),N,BRgp(1:np,0:N-1),dBRgp(1:np,0:N-1))
       dBRgp(1:np,0:N-1) = spread(kap(1:np),2,N)*dBRgp(:,:)
       BR0(1:np,0:N-1) =  bK(c%r*kap(:),N)
    end if

    aa(1:np,0:N-1) = c%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,0) = 0.0 ! make odd/even conformable
    bb(1:np,1:N-1) = c%coeff(lo:hi,n0+N:n0+2*N-2)

    ! dPot_dR
    dH(1:np,1) = sum(dBRgp(1:np,0:N-1)/BR0(1:np,0:N-1)* &
         & ( aa(1:np,0:N-1)*spread(cos(vr(0:N-1)*Pgp),1,np) + &
         &   bb(1:np,0:N-1)*spread(sin(vr(0:N-1)*Pgp),1,np) ),dim=2)

    ! dPot_dTheta
    dH(1:np,2) = sum(BRgp(1:np,0:N-1)/BR0(1:np,0:N-1)*spread(vr(0:N-1),1,np)* &
         & ( bb(1:np,0:N-1)*spread(cos(vr(0:N-1)*Pgp),1,np) - &
         &   aa(1:np,0:N-1)*spread(sin(vr(0:N-1)*Pgp),1,np) ),dim=2)

  end function circle_deriv
end module circular_elements

