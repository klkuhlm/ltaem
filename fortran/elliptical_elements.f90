! this module contains the basic functions defining the head or flux 
! effects of a elliptical element.

module elliptical_elements
  use constants, only : DP, PI
  use kappa_mod
  use time_mod

  implicit none
  private

  interface ellipse_match
     module procedure ellipse_match_self, ellipse_match_other
  end interface

contains
  function ellipse_match_self(e,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, match_result
    use mathieu_functions, only : ce, se, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    implicit none

    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, inten(in) :: idx ! indicates which value of p (global state)
    type(match_result) :: r

    integer :: j, N, M, lo, hi, nmax
    complex(DP), allocatable :: RMn(:,:), dRMn(:,:) ! mod. radial Mathieu function (K{e,o} or I{e,o})
    complex(DP) :: kap
    real(DP) :: cemat(1:e%M,0:e%N-1,0:1), semat(1:e%M,1:e%N-1,0:1)
    real(DP), dimension(0:e%N-1) :: vi
    real(DP), allocatable :: vs(:,:)

    N = e%N; M = e%M
    vi = real([(j,j=0,N-1)],DP)

    if (e%ibnd == 0) then
       allocate(r%LHS(2*M,4*N-2), r%RHS(2*M))
       lo = M+1; hi = 2*M
    elseif (e%calcin) then
       allocate(r%LHS(M,4*N-2), r%RHS(M))
       lo = 1; hi = M
    else
       allocate(r%LHS(M,2*N-1), r%RHS(M))
       lo = 1;  hi = M
    end if

    ! outside ang. mod. Mathieu fcns
    cemat(:,0:N-1,0) = ce(e%parent%mat(idx), vi(0:N-1), e%Pcm(1:M))
    semat(:,1:N-1,0) = se(e%parent%mat(idx), vi(1:N-1), e%Pcm(1:M))
    if (e%ibnd == 0 .or. (e%calcin .and. (e%ibnd == 1 .or. e%ibnd == -1))) then
       ! inside
       cemat(:,0:N-1,1) = ce(e%mat(idx), vi(0:N-1), e%Pcm(1:M))
       semat(:,1:N-1,1) = se(e%mat(idx), vi(1:N-1), e%Pcm(1:M))
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
    
    ! matching or specified total flux
    if (e%ibnd == 0 .or. e%ibnd == +1 .or. e%ibnd == 2) then
       allocate(RMn(0:N-1,0:1),dRMn(0:N-1,0:1))
       RMn(0:N-1,0) =   Ke(e%parent%mat(idx), vi(0:N-1), e%r) ! even fn
       RMn(1:N-1,1) =   Ko(e%parent%mat(idx), vi(1:N-1), e%r) ! odd fn
       dRMn(0:N-1,0) = dKe(e%parent%mat(idx), vi(0:N-1), e%r) ! even deriv
       dRMn(1:N-1,1) = dKo(e%parent%mat(idx), vi(1:N-1), e%r) ! odd deriv

       r%LHS(lo:hi,1:N) =       spread(dRMn(0:N-1)/RMn(0:N-1), 1,M)*cemat ! a_n flux
       r%LHS(lo:hi,N+1:2*N-1) = spread(dRMn(1:N-1)/RMn(1:N-1), 1,M)*semat ! b_n flux
       
       if (e%ibnd == 0 .or. (e%ibnd == 1 .and. e%calcin)) then
          RMn(0:N-1,0) =   Ie(e%mat(idx), vi(0:N-1), e%r)
          RMn(1:N-1,1) =   Io(e%mat(idx), vi(1:N-1), e%r)
          dRMn(0:N-1,0) = dIe(e%mat(idx), vi(0:N-1), e%r)
          dRMn(1:N-1,1) = dIo(e%mat(idx), vi(1:N-1), e%r)

          r%LHS(lo:hi,2*N:3*N-1) = spread(dRMn(0:N-1)/RMn(0:N-1), 1,M)*cemat ! c_n flux
          r%LHS(lo:hi,3*N:4*N-2) = spread(dRMn(1:N-1)/RMn(1:N-1), 1,M)*semat ! d_n flux
       end if
       deallocate(RMn,dRMn)
    end if
    
    ! setup RHS
    select case(e%ibnd)
    case(-1)
       ! put specified head on RHS
       r%RHS(1:M) = time(p,e%time,.false.)*e%bdryQ
    case(0)
       ! put constant area source term effects on RHS
       r%RHS(1:M) = -time(p,e%time,.true.)*e%areaQ*e%Ss/kappa(p,e%parent)**2
       r%RHS(M+1:2*M) = 0.0 ! area source has no flux effects
    case(1)
       ! put specified flux effects on RHS
       r%RHS(1:M) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI*e%r)
    case(2)
       ! specified flux (line-source of specified strength)
       ! even coefficients are computed analytically (odd are zero by symmetry)
       allocate(vs(0:e%ms-1,2))
       vs(:,1) = real([(j,j=0,e%ms-1)],DP) ! vector of integers
       vs(:,2) = (-1.0_DP)**vs(0:e%ms-1,1) ! vector of signs
       nmax = size(vs(0:N-1,2),1)
       r%RHS(1:M) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI)* &
            & Ke(e%parent%mat(idx), vi(0:N-1:2), e%r)/dKe(e%parent%mat(idx), vi(0:N-1:2), e%r)&
            & sum(spread(-vs(0:N-1:2,2)*sum(&
            & spread(vs(0:e%ms,2),2,nmax)*conjg(e%mat(idx)%A(:,0:nmax,0))/ &
            & (1.0_DP - vs(0:e%ms,2)**2),dim=1),dim=1,ncopies=M)*r%LHS(1:M,0:2*N-1:2),dim=2)
       deallocate(vs)
    end select
  end function ellipse_match_self

  function ellipse_match_other(e,el,dom,p) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, domain, matching, match_result
    use bessel_functions, only : bK, bI
    implicit none

    type(ellipse), intent(in) :: e ! source ellipse
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    type(match_result) :: r 

    integer :: j, src, targ, N, M, lo, hi
    real(DP), allocatable :: cmat(:,:), smat(:,:)
    real(DP), dimension(0:e%N-1) :: vi
    complex(DP), allocatable :: Bn(:,:), dBn(:,:), Bn0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    complex(DP) :: kap

    N = e%N ! number of coefficients in the source elliptical element
    targ = el%id; src = e%id
    vi = real([(j,j=0,N-1)],DP)

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then

       M = el%M
       if (el%ibnd == 0) then
          allocate(r%LHS(2*M,2*N-1), r%RHS(2*M))
          lo = M+1; hi = 2*M
       else
          allocate(r%LHS(M,2*N-1), r%RHS(M))
          lo = 1;  hi = M
       end if

       cmat = cos(outerprod(e%G(targ)%Pgm(1:M), vi(0:N-1)))
       smat = sin(outerprod(e%G(targ)%Pgm(1:M), vi(1:N-1)))
       allocate(Bn(M,0:N-1),Bn0(0:N-1))

       ! setup LHS 
       ! for matching or specified total head target elements
       if (el%ibnd == 0 .or. el%ibnd == -1) then

          if (dom%inclBg(src,targ)) then
             ! can the target element "see" the outside of the source element?
             ! use exterior Bessel functions (Kn)
             kap = kappa(p,e%parent)
             Bn(0:N-1,1:M) = bK(kap*e%G(targ)%Rgm(1:M),N)
             Bn0(0:N-1) =    bK(kap*e%r,N)
          else
             ! can target element "see" the inside of the source element?
             ! i.e., is the source element the parent?
             ! use interior Bessel functions (In)
             kap = kappa(p,e%element)
             Bn(0:N-1,1:M) = bI(kap*e%G(targ)%Rgm(1:M),N)
             Bn0(0:N-1) =    bI(kap*e%r,N)
          end if
          
          ! head effects on other element
          r%LHS(1:M,1:N) =       Bn(0:N-1,:)/spread(Bn0(0:N-1),1,M)*cmat/e%parent%K ! a_n || c_n
          r%LHS(1:M,N+1:2*N-1) = Bn(1:N-1,:)/spread(Bn0(1:N-1),1,M)*smat/e%parent%K ! b_n || d_n
          
          if (e%ibnd == 2 .and. dom%inclBg(src,targ)) then
             if (e%StorIn) then
                ! wellbore storage and skin from finite-radius well
                r%LHS(1:M,1) = -Bn0(0)*((2.0 + e%r**2*e%dskin*p/e%parent%T)/(2.0*PI*e%r) + &
                     & (Bn0(0)*e%r*p)/(2.0*PI*e%r*kap*Bn0(1)*e%parent%T))*r%LHS(1:M,1)
                r%RHS(1:M) = 0.0
             else
                ! specified flux (finite-radius well no storage)
                r%RHS(1:M) = Bn0(0)*time(p,e%time,.false.)*e%bdryQ/(2.0*PI*e%r*Bn0(1))*r%LHS(1:M,1)
                r%LHS(1:M,1) = 0.0
             end if
          end if
       end if
          
       ! for matching, specified total flux, or specified elemental flux target element
       if (el%ibnd == 0 .or. el%ibnd == +1 .or. el%ibnd == +2) then
          allocate(dBn(M,0:N-1), dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), &
               & dPot_dX(M,2*N-1),dPot_dY(M,2*N-1))

          ! flux effects of source circle on target element
          if (dom%inclBg(src,targ)) then
             ! use exterior Bessel functions (Kn)
             kap = kappa(p,e%parent) 
             call bKD(kap*e%G(targ)%Rgm(1:M),N+1,Bn,dBn)
             dBn = kap*dBn
             Bn0(0:N-1) = bK(kap*e%r,N)
          else
             ! use interior Bessel functions (In)
             kap = kappa(p,e%element)
             call bId(kap*e%G(targ)%Rgm(1:M),N+1,Bn,dBn)
             dBn = kap*dBn
             Bn0(0:N-1) = bI(kap*e%r,N)
          end if

          ! derivative wrt radius of source element
          dPot_dR(1:M,1:N) =       dBn(1:M,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat
          dPot_dR(1:M,N+1:2*N-1) = dBn(1:M,1:N-1)/spread(Bn0(1:N-1),1,M)*smat

          if (el%ibnd == 2 .and. dom%inclBg(src,targ)) then
             if (e%StorIn) then
                ! wellbore storage and skin from finite-radius well
                dPot_dR(1:M,1) = -Bn0(0)*((2.0 + e%r**2*e%dskin*p/e%parent%T)/(2.0*PI*e%r) + &
                     & (Bn0(0)*e%r*p)/(2.0*PI*e%r*kap*Bn0(1)*e%parent%T))*dPot_dR(1:M,1)
             else
                ! specified flux (finite-radius well no storage)
                dPot_dR(1:M,1) = Bn0(0)*time(p,e%time,.false.)*e%bdryQ/(2.0*PI*e%r*Bn0(1))*dPot_dR(1:M,1)
             end if
             dPot_dP(:,:) = 0.0 ! wells are radially-symmetric; no angular deriv. contribution
          else
             ! derivative wrt angle of source element for more general elliptical elements
             dPot_dP(1:M,1:N) =      -Bn(0:N-1,:)*spread(vi(0:N-1)/Bn0(0:N-1),1,M)*smat
             dPot_dP(1:M,N+1:2*N-1) = Bn(1:N-1,:)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*cmat
          end if

          ! project these from cylindrical onto Cartesian coordinates
          dPot_dX = dPot_dR*spread(cos(e%G(targ)%Pgm),2,2*N-1) - &
               & dPot_dP*spread(sin(e%G(targ)%Pgm)/e%G(targ)%Rgm,2,2*N-1)
          dPot_dY = dPot_dR*spread(sin(e%G(targ)%Pgm),2,2*N-1) + &
               & dPot_dP*spread(cos(e%G(targ)%Pgm)/e%G(targ)%Rgm,2,2*N-1)
          
          ! project from Cartesian to "radial" coordinate of target element
          if (el%id <= dom%num(1)) then
             ! other element is a circle
             if (el%ibnd == 2) then
                if (el%StorIn) then
                   ! other element is a well with wellbore storage (Type III BC)
                   r%LHS(1:M,1:N) =       Bn(0:N-1,:)/spread(Bn0(0:N-1),1,M)*smat/e%parent%K
                   r%LHS(1:M,N+1:2*N-1) = Bn(1:N-1,:)/spread(Bn0(1:N-1),1,M)*cmat/e%parent%K
                   
                   ! head effects of element
                   r%LHS(1:M,:) = -(el%r*p/el%parent%T)*r%LHS(1:M,:)
                   
                   ! radial flux effects of element
                   r%LHS(1:M,:) = r%LHS + (2.0 + el%r**2*el%dskin*p/el%parent%T)* &
                        & (dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                        &  dPot_dY*spread(sin(el%Pcm),2,2*N-1))
                end if
             else
                ! other element is a circle without wellbore storage
                r%LHS(lo:hi,:) = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                               & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
             end if
          else
             ! other element is an ellipse
             r%LHS(lo:hi,:) = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                            & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)
          end if
          deallocate(dBn,dPot_dR,dPot_dP,dPot_dX,dPot_dY)
       end if
       deallocate(Bn,Bn0)
    end if
  end function circle_match_other

end module elliptical_elements

