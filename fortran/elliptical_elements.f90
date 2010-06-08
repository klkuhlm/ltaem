! this module contains the basic functions defining the head or flux 
! effects of a circular element.

module circular_elements
  use constants, only : DP, PI
  use kappa_mod
  use time_mod

  implicit none
  private

  interface circle_match
     module procedure circle_match_self, circle_match_other
  end interface

contains
  function circle_match_self(c,p) result(r)
    use utility, only : outerprod
    use type_definitions, only : circle, match_result
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    type(match_result) :: r

    integer :: j, N, M, lo, hi
    complex(DP), allocatable :: Bn(:), dBn(:) ! mod. bessel function (K or I)
    complex(DP) :: kap
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi

    N = c%N; M = c%M
    vi = real([(j,j=0,N-1)],DP)

    if (c%ibnd == 0) then
       allocate(r%LHS(2*M,4*N-2), r%RHS(2*M))
       lo = M+1; hi = 2*M
    elseif (c%calcin) then
       allocate(r%LHS(M,4*N-2), r%RHS(M))
       lo = 1; hi = M
    else
       allocate(r%LHS(M,2*N-1), r%RHS(M))
       lo = 1;  hi = M
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
    if (c%ibnd == 0 .or. c%ibnd == +1 .or. c%ibnd == 2) then
       allocate(Bn(0:N-1),dBn(0:N-1))
       kap = kappa(p,c%parent) 
       call bKD(kap*c%r,N,Bn,dBn)
       dBn = kap*dBn

       r%LHS(lo:hi,1:N) =       spread(dBn(0:N-1)/Bn(0:N-1), 1,M)*cmat ! a_n flux
       r%LHS(lo:hi,N+1:2*N-1) = spread(dBn(1:N-1)/Bn(1:N-1), 1,M)*smat ! b_n flux
       
       if (c%ibnd == 0 .or. (c%ibnd == 1 .and. c%calcin)) then
          kap = kappa(p,c%element)
          call bID(kap*c%r,N,Bn,dBn)
          dBn = kap*dBn
          
          r%LHS(lo:hi,2*N:3*N-1) = spread(dBn(0:N-1)/Bn(0:N-1), 1,M)*cmat ! c_n flux
          r%LHS(lo:hi,3*N:4*N-2) = spread(dBn(1:N-1)/Bn(1:N-1), 1,M)*smat ! d_n flux
       end if
       deallocate(Bn,dBn)
    end if
    
    ! setup RHS
    select case(c%ibnd)
    case(-1)
       ! put specified head on RHS
       r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ
    case(0)
       ! put constant area source term effects on RHS
       r%RHS(1:M) = -time(p,c%time,.true.)*c%areaQ*c%Ss/kappa(p,c%parent)**2
       r%RHS(M+1:2*M) = 0.0 ! area source has no flux effects
    case(1)
       ! put specified flux effects on RHS
       r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r)
    case(2)
       allocate(Bn(0:1))
       Bn(0:1) = bK(kap*c%r,2)
          
       if (c%StorIn) then
          ! effects of wellbore storage and skin on finite-radius
          ! well, where a_0 is computed (generally depends on
          ! other elements, too; these show up in off-diagonal sub-matrices)
          r%LHS(1:M,1) = -Bn(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
               & (Bn(0)*c%r*p)/(2.0*PI*c%r*kap*Bn(1)*c%parent%T))*r%LHS(1:M,1)
          r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(PI*c%r*c%parent%T)
       else
          ! specified flux (finite-radius well no storage)
          ! a_0 coefficient is computed analytically
          r%RHS(1:M) = Bn(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Bn(1))*r%LHS(1:M,1)
          r%LHS(1:M,1) = 0.0
       end if
       deallocate(Bn)
    end select
  end function circle_match_self

  function circle_match_other(c,el,dom,p) result(r)
    use utility, only : outerprod
    use type_definitions, only : circle, domain, matching, match_result
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c ! source circle
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    type(match_result) :: r 

    integer :: j, src, targ, N, M
    real(DP), allocatable :: cmat(:,:), smat(:,:)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: Bn(:,:), dBn(:,:), Bn0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:
    complex(DP) :: kap

    N = c%N ! number of coefficients in the source circular element
    targ = el%id; src = c%id
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

       cmat = cos(outerprod(c%G(targ)%Pgm(1:M), vi(0:N-1)))
       smat = sin(outerprod(c%G(targ)%Pgm(1:M), vi(1:N-1)))
       allocate(Bn(M,0:N-1),Bn0(0:N-1))

       ! setup LHS 
       ! for matching or specified total head target elements
       if (el%ibnd == 0 .or. el%ibnd == -1) then

          if (dom%inclBg(src,targ)) then
             ! can the target element "see" the outside of the source element?
             ! use exterior Bessel functions (Kn)
             kap = kappa(p,c%parent)
             Bn(0:N-1,1:M) = bK(kap*c%G(targ)%Rgm(1:M),N)
             Bn0(0:N-1) =    bK(kap*c%r,N)
          else
             ! can target element "see" the inside of the source element?
             ! i.e., is the source element the parent?
             ! use interior Bessel functions (In)
             kap = kappa(p,c%element)
             Bn(0:N-1,1:M) = bI(kap*c%G(targ)%Rgm(1:M),N)
             Bn0(0:N-1) =    bI(kap*c%r,N)
          end if
          
          ! head effects on other element
          r%LHS(1:M,1:N) =       Bn(0:N-1,:)/spread(Bn0(0:N-1),1,M)*cmat/c%parent%K ! a_n || c_n
          r%LHS(1:M,N+1:2*N-1) = Bn(1:N-1,:)/spread(Bn0(1:N-1),1,M)*smat/c%parent%K ! b_n || d_n
          
          if (c%ibnd == 2 .and. dom%inclBg(src,targ)) then
             if (c%StorIn) then
                ! wellbore storage and skin from finite-radius well
                r%LHS(1:M,1) = -Bn0(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                     & (Bn0(0)*c%r*p)/(2.0*PI*c%r*kap*Bn0(1)*c%parent%T))*r%LHS(1:M,1)
                r%RHS(1:M) = 0.0
             else
                ! specified flux (finite-radius well no storage)
                r%RHS(1:M) = Bn0(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Bn0(1))*r%LHS(1:M,1)
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
             kap = kappa(p,c%parent) 
             call bKD(kap*c%G(targ)%Rgm(1:M),N+1,Bn,dBn)
             dBn = kap*dBn
             Bn0(0:N-1) = bK(kap*c%r,N)
          else
             ! use interior Bessel functions (In)
             kap = kappa(p,c%element)
             call bId(kap*c%G(targ)%Rgm(1:M),N+1,Bn,dBn)
             dBn = kap*dBn
             Bn0(0:N-1) = bI(kap*c%r,N)
          end if

          ! derivative wrt radius of source element
          dPot_dR(1:M,1:N) =       dBn(1:M,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat
          dPot_dR(1:M,N+1:2*N-1) = dBn(1:M,1:N-1)/spread(Bn0(1:N-1),1,M)*smat

          if (el%ibnd == 2 .and. dom%inclBg(src,targ)) then
             if (c%StorIn) then
                ! wellbore storage and skin from finite-radius well
                dPot_dR(1:M,1) = -Bn0(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                     & (Bn0(0)*c%r*p)/(2.0*PI*c%r*kap*Bn0(1)*c%parent%T))*dPot_dR(1:M,1)
             else
                ! specified flux (finite-radius well no storage)
                dPot_dR(1:M,1) = Bn0(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Bn0(1))*dPot_dR(1:M,1)
             end if
             dPot_dP(:,:) = 0.0 ! wells are radially-symmetric; no angular deriv. contribution
          else
             ! derivative wrt angle of source element for more general circular elements
             dPot_dP(1:M,1:N) =      -Bn(0:N-1,:)*spread(vi(0:N-1)/Bn0(0:N-1),1,M)*smat
             dPot_dP(1:M,N+1:2*N-1) = Bn(1:N-1,:)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*cmat
          end if

          ! project these from cylindrical onto Cartesian coordinates
          dPot_dX = dPot_dR*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
               & dPot_dP*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
          dPot_dY = dPot_dR*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
               & dPot_dP*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
          
          ! project from Cartesian to "radial" coordinate of target element
          if (el%id <= dom%num(1)) then
             ! other element is a circle
             if (el%ibnd == 2) then
                if (el%StorIn) then
                   ! other element is a well with wellbore storage (Type III BC)
                   r%LHS(1:M,1:N) =       Bn(0:N-1,:)/spread(Bn0(0:N-1),1,M)*smat/c%parent%K
                   r%LHS(1:M,N+1:2*N-1) = Bn(1:N-1,:)/spread(Bn0(1:N-1),1,M)*cmat/c%parent%K
                   
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

end module circular_elements

