! this module contains the basic functions defining the head or flux 
! effects of a circular element.

module circular_elements
  implicit none

  private
  public :: circle_match, circle_calc, circle_deriv

  interface circle_match
     module procedure circle_match_self, circle_match_other
  end interface
    
contains
  function circle_match_self(c,p) result(r)
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use time_mod, only : time
    use utility, only : outerprod
    use type_definitions, only : circle, match_result
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    type(match_result) :: r

    integer :: j, N, M, loM, hiM, ierr, nrows, ncols
    complex(DP), allocatable :: Bn(:), dBn(:) ! mod. bessel function (K or I)
    complex(DP) :: kap
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi

    N = c%N
    M = c%M
    vi = real([(j,j=0,N-1)],DP)

    if (c%ibnd == 0) then
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
          nrows = M
          ncols = 1
       else
          nrows = 0
          ncols = 0
       end if
    else
       nrows = M
       ncols = 2*N-1
       ! here only flux matching, no head matching, so loM:him = 1:M
       loM = 1
       hiM = M
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows), stat=ierr)
    if (ierr /= 0) then
       stop 'circular_elements.f90:circle_match_self() error allocating: r%LHS, r%RHS'
    end if

    cmat = cos(outerprod(c%Pcm(1:M), vi(0:N-1)))
    smat = sin(outerprod(c%Pcm(1:M), vi(1:N-1)))

    ! setup LHS
    ! matching or specified total head
    if (c%ibnd == 0 .or. c%ibnd == -1) then
       r%LHS(1:M,1:N) =       cmat/c%parent%K ! a_n head
       r%LHS(1:M,N+1:2*N-1) = smat/c%parent%K ! b_n head

       if (c%ibnd == 0 .or. (c%ibnd == -1 .and. c%calcin)) then
          r%LHS(1:M,2*N:3*N-1) = -cmat/c%K ! c_n head
          r%LHS(1:M,3*N:4*N-2) = -smat/c%K ! d_n head
       end if
    end if
    
    ! matching or specified total flux
    if (c%ibnd == 0 .or. c%ibnd == +1 .or. (c%ibnd == 2 .and. c%storIn)) then
       allocate(Bn(0:N-1),dBn(0:N-1), stat=ierr)
       if (ierr /= 0) stop 'circular_elements.f90 error allocating: Bn,dBn'
       kap = kappa(p,c%parent) 
       call dBK(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
       dBn(0:N-1) = kap*dBn(0:N-1)

       r%LHS(loM:hiM,1:N) =       spread(dBn(0:N-1)/Bn(0:N-1), 1,M)*cmat ! a_n flux
       r%LHS(loM:hiM,N+1:2*N-1) = spread(dBn(1:N-1)/Bn(1:N-1), 1,M)*smat ! b_n flux
       
       if (c%ibnd == 0 .or. (c%ibnd == 1 .and. c%calcin)) then
          kap = kappa(p,c%element)
          call dBI(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
          dBn(0:N-1) = kap*dBn(0:N-1)
          
          r%LHS(loM:hiM,2*N:3*N-1) = spread(dBn(0:N-1)/Bn(0:N-1), 1,M)*cmat ! c_n flux
          r%LHS(loM:hiM,3*N:4*N-2) = spread(dBn(1:N-1)/Bn(1:N-1), 1,M)*smat ! d_n flux
       end if
       deallocate(Bn,dBn, stat=ierr)
       if (ierr /= 0) stop 'circular_elements.f90 error deallocating: Bn,dBn'
    end if
    
    ! setup RHS
    select case(c%ibnd)
    case(-1)
       ! put specified head on RHS
       r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ
    case(0)
       ! put constant area source term effects on RHS
       r%RHS(1:M) = -time(p,c%time,.true.)*c%areaQ*c%Ss/kappa(p,c%element)**2
       r%RHS(M+1:2*M) = 0.0 ! area source has no flux effects
    case(1)
       ! put specified flux effects on RHS
       r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r)
    case(2)
       if (c%StorIn) then
          ! effects of wellbore storage and skin on finite-radius well
          ! effects of other elements on this one show up in off-diagonals
          r%LHS(1:,1) = storwell(c,p)*r%LHS(1:M,1)
          r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(PI*c%r*c%parent%T)
       end if
    end select
  end function circle_match_self

  function circle_match_other(c,el,dom,p) result(r)
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use time_mod, only : time
    use utility, only : outerprod
    use type_definitions, only : circle, domain, matching, match_result
    use bessel_functions, only : bK, bI, dbk, dbi
    implicit none

    type(circle), intent(in) :: c ! source circle
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    type(match_result) :: r

    integer :: j, src, targ, N, M, loM, hiM, loN, hiN, ierr, nrows, ncols
    real(DP) :: K
    real(DP), allocatable :: cmat(:,:), smat(:,:)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: Bn(:,:), dBn(:,:), Bn0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    complex(DP) :: kap

    N = c%N ! number of coefficients in the source circular element
    targ = el%id
    src = c%id
    vi = real([(j,j=0,N-1)],DP)

    M = el%M
    ! target element determines number of rows
    if (el%ibnd == 0) then
       nrows = 2*M
       loM = M+1
       hiM = 2*M
    elseif (el%ibnd == 2 .and. .not. el%storIn) then
       nrows = 0
    else
       nrows = M
       loM = 1
       hiM = M
    end if
    ! source element determines number of columns
    if (c%ibnd == 0) then
       ncols = 4*N-2
    elseif (c%ibnd == 2 .and. .not. c%storIn) then
       ncols = 1  ! reset to zero at end of routine
    else
       ncols = 2*N-1
    end if
    
    allocate(r%LHS(nrows,ncols), r%RHS(nrows), stat=ierr)
    if (ierr /= 0) stop 'circular_elements.f90:circle_match_other() error allocating:&
            &r%LHS, r%RHS'
    r%LHS = 0.0
    r%RHS = 0.0

    if (nrows > 0) then
       if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then

          allocate(Bn(1:M,0:N-1),Bn0(0:N-1),cmat(1:M,0:N-1),smat(1:M,0:N-1), stat=ierr)
          if (ierr /= 0) stop 'circular_elements.f90:circle_match_other() error allocating:&
               &Bn, Bn0, cmat, smat'

          cmat(1:M,0:N-1) = cos(outerprod(c%G(targ)%Pgm(1:M), vi(0:N-1)))
          smat(1:M,0:N-1) = sin(outerprod(c%G(targ)%Pgm(1:M), vi(0:N-1)))

          ! setup LHS 
          ! $$$$$$$$$$ head effects of source (c) on target (el) $$$$$$$$$$
          ! for matching or specified total head target elements
          if (el%ibnd == 0 .or. el%ibnd == -1) then

             if (dom%inclBg(src,targ)) then
                ! can the target element "see" the outside of the source element?
                ! use exterior Bessel functions (Kn)
                kap = kappa(p,c%parent)
                Bn(1:M,0:N-1) = bK(kap*c%G(targ)%Rgm(1:M),N)
                Bn0(0:N-1) =    bK(kap*c%r,N)

                ! head effects on outside of other element
                r%LHS(1:M,1:N) =       Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1)/c%parent%K ! a_n
                if (N > 1) then
                   r%LHS(1:M,N+1:2*N-1) = Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/c%parent%K ! b_n
                end if

                ! outside is always in 1:2*N-1 slot
                loN = 1
                hiN = 2*N-1

             else
                ! can target element "see" the inside of the source element?
                ! i.e., is the source element the parent?
                ! use interior Bessel functions (In)
                kap = kappa(p,c%element)
                Bn(1:M,0:N-1) = bI(kap*c%G(targ)%Rgm(1:M),N)
                Bn0(0:N-1) =    bI(kap*c%r,N)

                if (c%ibnd == 0) then
                   ! is source the inside of matching element?
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! is source inside of specified head/flux element? (no other part)
                   loN = 1
                   hiN = 2*N-1
                end if

                ! head effects on other element
                r%LHS(1:M,loN:loN+N-1) = Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1)/c%K ! c_n
                if (N > 1) then
                   r%LHS(1:M,loN+N:hiN)   = Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/c%K ! d_n
                end if
             end if

             if (c%ibnd == 2 .and. dom%inclBg(src,targ)) then
                if (c%StorIn) then

                   ! wellbore storage and skin from finite-radius well
                   r%LHS(1:M,1) = storwell(c,p)*r%LHS(1:M,1)
                   r%RHS(1:M) = 0.0
                else
                   ! specified flux (finite-radius well no storage)
                   r%RHS(1:M) = well(c,p)*r%LHS(1:M,1)
                   ! no LHS for specified simple well (deallocate away at end of routine)
                end if
             end if
          end if
          ! $$$$$$$$$$ flux effects of source (c) on target (el) $$$$$$$$$$
          ! for matching, specified total flux, or well with wellbore storage target element
          if (el%ibnd == 0 .or. el%ibnd == +1 .or. (el%ibnd == +2 .and. el%storIn)) then
             allocate(dBn(M,0:N-1), dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), &
                  & dPot_dX(M,2*N-1), dPot_dY(M,2*N-1), stat=ierr)
             if (ierr /= 0) stop 'circular_elements.f90 error allocating:&
                  & dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY'

             ! flux effects of source circle on target element
             if (dom%inclBg(src,targ)) then
                ! use exterior Bessel functions (Kn)
                kap = kappa(p,c%parent) 
                call dBK(kap*c%G(targ)%Rgm(1:M),N,Bn(1:M,0:N-1),dBn(1:M,0:N-1))
                dBn(1:M,0:N-1) = kap*dBn(1:M,0:N-1)
                Bn0(0:N-1) = bK(kap*c%r,N)
                K = c%parent%K

                ! outside part is always 1:2*N-1
                loN = 1
                hiN = 2*N-1

             else
                ! use interior Bessel functions (In)
                kap = kappa(p,c%element)
                call dBI(kap*c%G(targ)%Rgm(1:M),N,Bn(1:M,0:N-1),dBn(1:M,0:N-1))
                dBn(1:M,0:N-1) = kap*dBn(1:M,0:N-1)
                Bn0(0:N-1) = bI(kap*c%r,N)
                K = c%K

                if (c%ibnd == 0) then
                   ! inside of matching element
                   loN = 2*N
                   hiN = 4*N-2
                else
                   ! inside of specified flux element
                   loN = 1
                   hiN = 2*N-1
                end if
             end if
             ! derivative wrt radius of source element
             dPot_dR(1:M,1:N) =       dBn(1:M,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(1:M,0:N-1)
             if (N>1) then
                dPot_dR(1:M,N+1:2*N-1) = dBn(1:M,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(1:M,1:N-1)
             end if

             if (el%ibnd == 2 .and. dom%inclBg(src,targ)) then
                if (c%StorIn) then
                   ! wellbore storage and skin from finite-radius well
                   dPot_dR(1:M,1) = storwell(c,p)*dPot_dR(1:M,1)
                else
                   ! specified flux (finite-radius well no storage)
                   dPot_dR(1:M,1) = well(c,p)*dPot_dR(1:M,1)
                end if

                dPot_dP(1:M,1:2*N-1) = 0.0 ! wells have angular symmetry
             else
                ! derivative wrt angle of source element for more general circular elements
                dPot_dP(1:M,1:N) =      -Bn(:,0:N-1)*spread(vi(0:N-1)/Bn0(0:N-1),1,M)*smat(:,0:N-1)
                if (N>1) then
                   dPot_dP(1:M,N+1:2*N-1) = Bn(:,1:N-1)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*cmat(:,1:N-1)
                end if

             end if
             ! project these from cylindrical onto Cartesian coordinates
             dPot_dX(1:M,1:2*N-1) = dPot_dR*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
                  & dPot_dP*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
             dPot_dY(1:M,1:2*N-1) = dPot_dR*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
                  & dPot_dP*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
             ! project from Cartesian to "radial" coordinate of target element
             if (el%id <= dom%num(1)) then
                ! other element is a circle
                if (el%ibnd == 2) then
                   if (el%StorIn) then
                      ! other element is a well with wellbore storage (Type III BC)
                      ! need head effects too
                      r%LHS(1:M,loN:loN+N-1) = Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat/K
                      r%LHS(1:M,loN+N:hiN) =   Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1)/K

                      ! head effects of element
                      r%LHS(1:M,loN:hiN) = -(el%r*p/el%parent%T)*r%LHS(1:M,loN:hiN)

                      ! radial flux effects of element
                      r%LHS(1:M,loN:hiN) = r%LHS + (2.0 + el%r**2*el%dskin*p/el%parent%T)* &
                           & (dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                           &  dPot_dY*spread(sin(el%Pcm),2,2*N-1))
                   end if
                else
                   ! other element is a 'normal' circular element without wellbore storage
                   r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                        & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
                end if
             else
                ! other element is an ellipse
                r%LHS(loM:hiM,loN:hiN) = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                     & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)
             end if

             deallocate(dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY, stat=ierr)
             if (ierr /= 0) then
                stop 'circular_elements.f90 error deallocating:&
                     & dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY'
             end if

          end if
          deallocate(Bn,Bn0,cmat,smat, stat=ierr)
          if (ierr /= 0) then
             stop 'circular_elements.f90 error deallocating: Bn,Bn0,cmat,smat'
          end if

       end if
    end if
       
    if (c%ibnd == 2 .and. (.not. c%storin)) then
       if (el%ibnd == 0) then
          hiM = 2*M
       else
          hiM = M
       end if
       deallocate(r%LHS,stat=ierr)
       if (ierr /= 0) then
          stop 'circular_elements:circle_match_self() error deallocating r%LHS'
       end if
       
       allocate(r%LHS(hiM,0),stat=ierr)
       if (ierr /= 0) then
          stop 'circular_elements:circle_match_self() error re-allocating r%LHS'
       end if
       
    end if

  end function circle_match_other

  function well(c,p) result(a0)
    ! this function returns the a_0 coefficient for a simple "well"
    use constants, only : DP, PI
    use kappa_mod, only : kappa
    use time_mod, only : time
    use type_definitions, only : circle
    use bessel_functions, only : bK
    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    complex(DP) :: a0
    complex(DP), dimension(0:1) ::Kn
    
    Kn(0:1) = bK(kappa(p,c%parent)*c%r,2)
    a0 = Kn(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn(1))
    
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
    a0 = -Kn(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
               & (Kn(0)*c%r*p)/(2.0*PI*c%r*kap*Kn(1)*c%parent%T))    
  end function storwell

  function circle_calc(p,c,lo,hi,Rgp,Pgp,inside) result(H)
    use constants, only : DP
    use kappa_mod, only : kappa
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

#ifdef DEBUG
    print *, 'circle_calc: p:',p,' c:',c%id,' lo:',lo,' hi:',hi,' Rgp:',Rgp,' Pgp:',Pgp,' inside:',inside
#endif

    N = c%N
    np = size(p,1)
    vr(0:N-1) = real([(i,i=0,N-1)],DP)

    if (inside) then
       if (c%match) then
          n0 = 2*N ! inside of matching circle
       else
          n0 = 1   ! inside of specified boundary circle
       end if
       kap(1:np) = kappa(p(1:np),c%element)
       BRgp(1:np,0:N-1) = bK(Rgp*kap(1:np),N)
       BR0(1:np,0:N-1) =  bK(c%r*kap(1:np),N)
    else
       n0 = 1
       kap(1:np) = kappa(p(1:np),c%parent)
       BRgp(1:np,0:N-1) = bI(Rgp*kap(1:np),N)
       BR0(1:np,0:N-1) =  bI(c%r*kap(1:np),N)
    end if
    aa(1:np,0:N-1) = c%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,0) = 0.0 ! insert zero to make odd/even same shape
    bb(1:np,1:N-1) = c%coeff(lo:hi,n0+N:n0+2*N-2)
    H(1:np) = sum(BRgp(1:np,0:N-1)/BR0(1:np,0:N-1)* &
         & ( aa(1:np,0:N-1)*spread(cos(vr(0:N-1)*Pgp),1,np) + &
         &   bb(1:np,0:N-1)*spread(sin(vr(0:N-1)*Pgp),1,np) ),dim=2)
  end function circle_calc

  function circle_deriv(p,c,lo,hi,Rgp,Pgp,inside) result(dH)
    use constants, only : DP
    use kappa_mod, only : kappa
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

#ifdef DEBUG
    print *, 'circle_deriv: p:',p,' c:',c%id,' lo:',lo,' hi:',hi,' Rgp:',Rgp,' Pgp:',Pgp,' inside:',inside
#endif

    N = c%N
    np = size(p,1)
    vr(0:N-1) = real([(i,i=0,N-1)],DP)

    if (inside) then
       if (c%match) then
          n0 = 2*N ! inside of matching circle
       else
          n0 = 1   ! inside of specified boundary circle
       end if
       kap(1:np) = kappa(p(:),c%element)
       call dBK(Rgp*kap(1:np),N,BRgp(1:np,0:N-1),dBRgp(1:np,0:N-1))
       dBRgp(1:np,0:N-1) = spread(kap(1:np),dim=2,ncopies=N)*dBRgp(1:np,0:N-1)
       BR0(1:np,0:N-1) = bK(c%r*kap(1:np),N)
    else
       n0 = 1
       kap(1:np) = kappa(p(:),c%parent)
       call dBI(Rgp*kap(:),N,BRgp(1:np,0:N-1),dBRgp(1:np,0:N-1))
       dBRgp(1:np,0:N-1) = spread(kap(1:np),dim=2,ncopies=N)*dBRgp(1:np,0:N-1)
       BR0(1:np,0:N-1) =  bI(c%r*kap(1:np),N)
    end if
    aa(1:np,0:N-1) = c%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,0) = 0.0 ! insert zero to make odd/even same shape
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

