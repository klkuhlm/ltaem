! this module contains the basic functions defining the head or flux 
! effects of a circular or elliptical element, as well as the time
! behaviors, and the routine for computing general kappa for MHE.

module elements
  use constants, only : DP, PI
  use kappa_mod
  use time_mod

  implicit none

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface circle_head
     module procedure circle_match_head_self, &
          & circle_match_head_other, circle_head_calc
  end interface
  interface circle_flux
     module procedure circle_match_flux_self, &
          & circle_match_flux_other, circle_flux_calc
  end interface

  interface ellipse_head
     module procedure ellipse_match_head_self, &
          & ellipse_match_head_other, ellipse_head_calc
  end interface
  interface ellipse_flux
     module procedure ellipse_match_flux_self, &
          & ellipse_match_flux_other, ellipse_flux_calc
  end interface

contains
  function circle_match_head_self(c,p) result(r)
    use utility, only : outerprod
    use type_definitions, only : circle, match_result
    implicit none

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    type(match_result) :: r

    integer :: j, N, M
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi

    N = c%N; M = c%M
    vi = real([(j,j=0,N-1)],DP)

    ! LHS dim=2 is 4N-2 for matching, 2N-1 for spec. total head
    allocate(r%LHS(M,(2*N-1)*(2-abs(c%ibnd))), r%RHS(M))
    cmat = cos(outerprod(c%Pcm(1:M), vi(0:N-1)))
    smat = sin(outerprod(c%Pcm(1:M), vi(1:N-1)))

    ! setup LHS
    ! matching or specified head (always first M rows); no dependence on p
    ! ibnd==-2 would be here, but it doesn't actually make physical sense
    if (c%ibnd==0 .or. c%ibnd==-1) then

       r%LHS(1:M,1:N) =       cmat/c%parent%K
       r%LHS(1:M,N+1:2*N-1) = smat/c%parent%K
       
       ! setup RHS
       select case(c%ibnd)
       case(-1)
          ! put specified head on RHS
          r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ
       case(0)
          ! put constant area source term effects on RHS
          r%RHS(1:M) = -time(p,c%time,.true.)*c%areaQ*c%Ss/kappa(p,c%parent)**2
       end select

       if (c%ibnd==0 .or. c%calcin) then
          r%LHS(1:M,2*N:3*N-1) = -cmat/c%K
          r%LHS(1:M,3*N:4*N-2) = -smat/c%K
       end if
    else
       r%LHS = -huge(1.0)
       r%RHS =  huge(1.0)
    end if
  end function circle_match_head_self

  function circle_match_flux_self(c,p) result(r)
    use utility, only : outerprod
    use type_definitions, only : circle, match_result
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    complex(DP), intent(in) :: p
    type(match_result) :: r

    complex(DP), dimension(c%M,2*c%N-1) :: tmp
    integer :: j, N, M
    complex(DP), allocatable :: Kn(:), dKn(:), In(:), dIn(:)
    complex(DP) :: kap
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi

    N = c%N; M = c%M
    vi = real([(j,j=0,N-1)],DP)

    allocate(r%LHS(M,(2*N-1)*(2-abs(c%ibnd))), r%RHS(M))
    cmat = cos(outerprod(c%Pcm(1:M), vi(0:N-1)))
    smat = sin(outerprod(c%Pcm(1:M), vi(1:N-1)))

    ! matching (second M) or specified flux (first M); depends on p
    if (c%ibnd==0 .or. c%ibnd==1 .or. c%ibnd==2) then
       allocate(Kn(0:N),dKn(0:N))
       kap = kappa(p,c%parent) 
       call bKD(kap*c%r,N+1,Kn,dKn)
       dKn = kap*dKn

       tmp(1:M,1:N) =       spread(dKn(0:N-1)/Kn(0:N-1), 1,M)*cmat/c%parent%K
       tmp(1:M,N+1:2*N-1) = spread(dKn(1:N-1)/Kn(1:N-1), 1,M)*smat/c%parent%K
       deallocate(Kn,dKn)

       select case(c%ibnd)
       case(2)
          allocate(Kn(0:1))
          Kn(0:1) = bK(kap*c%r,2)
          
          if (c%StorIn) then
             ! effects of wellbore storage and skin on finite-radius
             ! well, where a_0 is computed (generally depends on
             ! other elements, too; these show up in off-diagonal sub-matrices)
             r%LHS(1:M,1) = -Kn(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                  & (Kn(0)*c%r*p)/(2.0*PI*c%r*kap*Kn(1)*c%parent%T))*tmp(1:M,1)
             r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(PI*c%r*c%parent%T)
          else
             ! specified flux (finite-radius well no storage)
             ! a_0 coefficient is computed analytically
             r%LHS(1:M,1) = 0.0
             r%RHS(1:M) = Kn(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn(1))*tmp(1:M,1)
          end if
          deallocate(Kn)
       case(1)
          ! put specified flux effects on RHS
          r%LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          r%RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r)
       case(0)
          ! no area source term effects on flux matchinig
          r%LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          r%RHS(1:M) = 0.0 
       end select
    
       if (c%ibnd==0 .or. c%calcin) then
          allocate(In(0:N),dIn(0:N))
          kap = kappa(p,c%element)
          call bID(kap*c%r,N+1,In,dIn)
          dIn = kap*dIn
          
          r%LHS(1:M,2*N:3*N-1) = spread(dIn(0:N-1)/In(0:N-1), 1,M)*cmat/c%K
          r%LHS(1:M,3*N:4*N-2) = spread(dIn(1:N-1)/In(1:N-1), 1,M)*smat/c%K
          deallocate(dIn,In)
       end if
    else
       r%LHS = -huge(1.0)
       r%RHS =  huge(1.0)
    end if
  end function circle_match_flux_self

  function circle_match_head_other(c,el,dom,p) result(r)
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
    complex(DP), allocatable :: Kn(:,:), Kn0(:), In(:,:), In0(:)
    complex(DP), allocatable :: tmp(:,:)
    complex(DP) :: kap

    N = c%N ! number of coefficients in the source circular element
    targ = el%id; src = c%id
    vi = real([(j,j=0,N-1)],DP)

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then
       ! common stuff between both branches below
       M = el%M
       allocate(r%LHS(M,(2*N-1)*(2-abs(el%ibnd))), r%RHS(M), &
            & cmat(M,0:N-1), smat(M,1:N-1))
       cmat = cos(outerprod(c%G(targ)%Pgm(1:M), vi(0:N-1)))
       smat = sin(outerprod(c%G(targ)%Pgm(1:M), vi(1:N-1)))

       ! can the target element "see" the outside of the source element?
       if (dom%inclBg(src,targ)) then
          allocate(tmp(M,2*N-1))

          ! setup LHS (head effects due to source element)
          ! for constant head (-1), or matching (0)
          if (el%ibnd==0 .or. el%ibnd==-1) then
             allocate(Kn(M,0:N-1),Kn0(0:N-1))
             kap = kappa(p,c%parent)
             Kn(0:N-1,1:M) = bK(kap*c%G(targ)%Rgm(1:M),N)
             Kn0(0:N-1) = bK(kap*c%r,N)

             ! head effects on other element
             tmp(1:M,1:N) =       Kn(0:N-1,:)/spread(Kn0(0:N-1),1,M)*cmat/c%parent%K
             tmp(1:M,N+1:2*N-1) = Kn(1:N-1,:)/spread(Kn0(1:N-1),1,M)*smat/c%parent%K
             deallocate(Kn,Kn0)

             select case (c%ibnd)
             case(2)
                allocate(Kn0(0:1))
                kap = kappa(p,c%parent)
                Kn0(0:1) = bK(kap*c%r,2)

                if (c%StorIn) then
                   ! wellbore storage and skin from finite-radius well
                   r%LHS(1:M,1) = -Kn0(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                        & (Kn0(0)*c%r*p)/(2.0*PI*c%r*kap*Kn0(1)*c%parent%T))*tmp(1:M,1)
                   r%RHS(1:M) = 0.0
                else
                   ! specified flux (finite-radius well no storage)
                   R%LHS(1:M,1) = 0.0
                   r%RHS(1:M) = Kn0(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn0(1))*tmp(1:M,1)
                end if
             case(-1,0,1)
                r%LHS(1:M,1:2*N-1) = tmp(1:M,1:2*N-1)
             end select
          end if
       else
          ! can target element "see" the inside of the source element?
          !  i.e., is the source element the parent?

          if (el%ibnd==0 .or. el%ibnd==-1) then
             allocate(In(M,0:N-1),In0(0:N-1))
             kap = kappa(p,c%element)
             In(0:N-1,1:M) = bI(kap*c%G(targ)%Rgm(1:M),N)
             In0(0:N-1) = bI(kap*c%r,N)

             ! head effects on other element
             r%LHS(1:M,1:N) =       In(0:N-1,:)/spread(In0(0:N-1),1,M)*cmat/c%K
             r%LHS(1:M,N+1:2*N-1) = In(1:N-1,:)/spread(In0(1:N-1),1,M)*smat/c%K
             deallocate(In,In0)
          end if
       end if
    end if
  end function circle_match_head_other

  function circle_match_flux_other(c,el,dom,p) result(r)
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
    complex(DP), allocatable :: Kn(:,:), Kn0(:), In(:,:), In0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    complex(DP) :: kap

    N = c%N ! number of coefficients in the source circular element
    targ = el%id; src = c%id
    vi = real([(j,j=0,N-1)],DP)

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then
       ! common stuff between both branches below
       M = el%M
       allocate(r%LHS(M,(2*N-1)*(2-abs(el%ibnd))), r%RHS(M), &
            & cmat(M,0:N-1), smat(M,1:N-1))
       cmat = cos(outerprod(c%G(targ)%Pgm(1:M), vi(0:N-1)))
       smat = sin(outerprod(c%G(targ)%Pgm(1:M), vi(1:N-1)))

       ! can the target element "see" the outside of the source element?
       if (dom%inclBg(src,targ)) then
          allocate(tmp(M,2*N-1),2)

          ! setup LHS (flux effects due to source element)
          ! for constant head (-1), or matching (0)
          if (el%ibnd==0 .or. el%ibnd==-1) then

             ! flux effects of source circle on target element
             allocate(Kn(M,0:N),dKn(M,0:N),Kn0(0:N-1), &
                  & dPot_dR(M,2*N-1),dPot_dP(M,2*N-1),dPot_dX(M,2*N-1),dPot_dY(M,2*N-1))
             kap = kappa(p,c%parent) 
             call bKD(kap*c%G(targ)%Rgm(1:M),N+1,Kn,dKn)
             Kn0(0:N-1) = bK(kap*c%r,N)
             dKn = kap*dKn

             ! derivative wrt radius of source element             
             dPot_dR(:,1:N) =       dKn(0:N-1)/spread(Kn(0:N-1),1,M)*cmat/c%parent%K
             dPot_dR(:,N+1:2*N-1) = dKn(1:N-1)/spread(Kn(1:N-1),1,M)*smat/c%parent%K
             deallocate(dKn)

             ! derivative wrt angle of source element
             dPot_dP(:,1:N) =      -Kn(0:N-1,:)*spread(vi(0:N-1)/Kn0(0:N-1),1,M)*smat/c%parent%K
             dPot_dP(:,N+1:2*N-1) = Kn(1:N-1,:)*spread(vi(1:N-1)/Kn0(1:N-1),1,M)*cmat/c%parent%K
             deallocate(Kn)

             ! project these from cylindrical onto Cartesian coordinates
             dPot_dX = dPot_dR*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
                     & dPot_dP*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

             dPot_dY = dPot_dR*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
                     & dPot_dP*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

             ! project from Cartesian to "radial" coordinate of target element
             if (el%id <= dom%num(1)) then
                ! other element is a circle
                r%LHS = dPot_dX*spread() + dPot_dY*spread()

             else
                ! other element is an ellipse
                r%LHS = dPot_dX*spread(cos()*sinh()*f**2*(cosh(2*) - cos(2*))) + &
                      & dPot_dY*spread()

             end if
             

             deallocate(Kn,Kn0)

             select case (c%ibnd)
             case(2)
                allocate(Kn0(0:1))
                kap = kappa(p,c%parent)
                Kn0(0:1) = bK(kap*c%r,2)

                if (c%StorIn) then
                   ! wellbore storage and skin from finite-radius well
                   r%LHS(1:M,1) = -Kn0(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                        & (Kn0(0)*c%r*p)/(2.0*PI*c%r*kap*Kn0(1)*c%parent%T))*tmp(1:M,1)
                   r%RHS(1:M) = 0.0
                else
                   ! specified flux (finite-radius well no storage)
                   R%LHS(1:M,1) = 0.0
                   r%RHS(1:M) = Kn0(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn0(1))*tmp(1:M,1)
                end if
             case(-1,0,1)
                r%LHS(1:M,1:2*N-1) = tmp(1:M,1:2*N-1)
             end select
          end if
       else
          ! can target element "see" the inside of the source element?
          !  i.e., is the source element the parent?

          if (el%ibnd==0 .or. el%ibnd==-1) then
             allocate(In(M,0:N-1),In0(0:N-1))
             kap = kappa(p,c%element)
             In(0:N-1,1:M) = transpose(bI(kap*c%G(targ)%Rgm(1:M),N))
             In0(0:N-1) = bI(kap*c%r,N)

             ! head effects on other element
             r%LHS(1:M,1:N) =       In(0:N-1,:)/spread(In0(0:N-1),1,M)*cmat/c%K
             r%LHS(1:M,N+1:2*N-1) = In(1:N-1,:)/spread(In0(1:N-1),1,M)*smat/c%K
             deallocate(In,In0)
          end if
       end if
    end if

  end function circle_match_flux_other

  

end module elements

