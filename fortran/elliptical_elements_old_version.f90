! this module contains the basic functions defining the head or flux 
! effects of an elliptical element.

module elliptical_elements
  use constants, only : DP, PI
  use kappa_mod
  use time_mod

  implicit none

  private
  public :: ellipse_head, ellipse_flux

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface ellipse_head
     module procedure ellipse_match_head_self, ellipse_match_head_other
  end interface
  interface ellipse_flux
     module procedure ellipse_match_flux_self, ellipse_match_flux_other
  end interface

contains
  function ellipse_match_head_self(e,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, match_result
    use mathieu_functions, only : mathieu_init, ce, se
    implicit none

    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx  
    type(match_result) :: r

    integer :: j, N, M
    real(DP) :: cemat(1:e%M,0:e%N-1), semat(1:e%M,1:e%N-1)
    real(DP), dimension(0:e%N-1) :: vi

    N = e%N; M = e%M
    vi = real([(j,j=0,N-1)],DP)

    ! LHS dim=2 is 4N-2 for matching, 2N-1 for spec. total head
    allocate(r%LHS(M,(2*N-1)*(2-abs(e%ibnd))), r%RHS(M))

    ! setup LHS
    ! matching or specified head (always first M rows); no dependence on p
    ! ibnd==-2 would be here, but it doesn't actually make physical sense
    if (e%ibnd==0 .or. e%ibnd==-1) then

       cemat = ce(e%parent%mat(idx), vi(0:N-1), e%Pcm(1:M))
       semat = se(e%parent%mat(idx), vi(1:N-1), e%Pcm(1:M))

       r%LHS(1:M,1:N) =       cemat/e%parent%K
       r%LHS(1:M,N+1:2*N-1) = semat/e%parent%K
       
       ! setup RHS
       select case(e%ibnd)
       case(-1)
          ! put specified head on RHS
          r%RHS(1:M) = time(p,e%time,.false.)*e%bdryQ
       case(0)
          ! put constant area source term effects on RHS
          r%RHS(1:M) = -time(p,e%time,.true.)*e%areaQ*e%Ss/kappa(p,e%parent)**2
       end select

       if (e%ibnd==0 .or. e%calcin) then
          r%LHS(1:M,2*N:3*N-1) = -cemat/e%K
          r%LHS(1:M,3*N:4*N-2) = -semat/e%K
       end if
    else
       r%LHS = -huge(1.0)
       r%RHS =  huge(1.0)
    end if
  end function ellipse_match_head_self

  function ellipse_match_flux_self(e,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, match_result
    use mathieu_functions, only : ce, se, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    implicit none

    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    type(match_result) :: r

    integer :: j, N, M
    complex(DP), allocatable :: mR(:,:), dmR(:,:)
    complex(DP) :: kap
    real(DP) :: cmat(1:e%M,0:e%N-1), smat(1:e%M,1:e%N-1)
    real(DP), dimension(0:e%N-1) :: vi

    N = e%N; M = e%M
    vi = real([(j,j=0,N-1)],DP)

    allocate(r%LHS(M,(2*N-1)*(2-abs(e%ibnd))), r%RHS(M))

    ! matching (second M) or specified flux (first M); depends on p
    if (e%ibnd==0 .or. e%ibnd==1 .or. e%ibnd==2) then
       cemat = ce(e%parent%mat(idx), vi(0:N-1), e%Pcm(1:M)) ! outside
       semat = se(e%parent%mat(idx), vi(1:N-1), e%Pcm(1:M))

       allocate(mR(0:N-1,0:1),dmR(0:N-1,0:1))
       mR(0:N-1,0) =   Ke(e%parent%mat(idx), vi(0:N-1), e%r) ! even fn
       mR(1:N-1,1) =   Ko(e%parent%mat(idx), vi(1:N-1), e%r) ! odd fn
       dmR(0:N-1,0) = dKe(e%parent%mat(idx), vi(0:N-1), e%r) ! even deriv
       dmR(1:N-1,1) = dKo(e%parent%mat(idx), vi(1:N-1), e%r) ! odd deriv

       r%LHS(1:M,1:N) =       spread(dmR(0:N-1,0)/mR(0:N-1,0), 1,M)*cemat/e%parent%K
       r%LHS(1:M,N+1:2*N-1) = spread(dmR(1:N-1,1)/mR(1:N-1,1), 1,M)*semat/e%parent%K
       deallocate(mR,dmR)

       select case(e%ibnd)
       case(2)
          ! specified flux (line source/sink)
          ! a_n & b_n coefficients computed analytically
          r%RHS(1:M) = Kn(0)*time(p,e%time,.false.)*e%bdryQ/(2.0*PI*e%r*Kn(1))*r%LHS(1:M,1)
          r%LHS(1:M,1) = 0.0
       case(1)
          ! put specified flux effects on RHS
          r%LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          r%RHS(1:M) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI*e%r)
       case(0)
          ! no area source term effects on flux matchinig
          r%LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          r%RHS(1:M) = 0.0 
       end select
    
       if (e%ibnd==0 .or. e%calcin) then
          cemat = ce(e%mat(idx), vi(0:N-1), e%Pcm(1:M)) ! inside
          semat = se(e%mat(idx), vi(1:N-1), e%Pcm(1:M))

          allocate(mR(0:N-1,0:1),dmR(0:N-1,0:1))
          mR(0:N-1,0) =   Ie(e%mat(idx), vi(0:N-1), e%r)
          mR(1:N-1,1) =   Io(e%mat(idx), vi(1:N-1), e%r)
          dmR(0:N-1,0) = dIe(e%mat(idx), vi(0:N-1), e%r)
          dmR(1:N-1,1) = dIo(e%mat(idx), vi(1:N-1), e%r)
          
          r%LHS(1:M,2*N:3*N-1) = spread(dmR(0:N-1,0)/mR(0:N-1,0), 1,M)*cemat/e%K
          r%LHS(1:M,3*N:4*N-2) = spread(dmR(1:N-1,1)/mR(1:N-1,1), 1,M)*semat/e%K
          deallocate(dmR,mR)
       end if
    else
       r%LHS = -huge(1.0)
       r%RHS =  huge(1.0)
    end if
  end function ellipse_match_flux_self

  function ellipse_match_head_other(e,el,dom,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, domain, matching, match_result
    use mathieu_functions, only : ce, se, Ke, Ko, Ie, Io
    implicit none

    type(ellipse), intent(in) :: e ! source ellipse
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    type(match_result) :: r 

    integer :: j, src, targ, N, M
    real(DP), allocatable :: cemat(:,:), semat(:,:)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: mR(:,:,:), mR0(:,:)

    N = c%N ! number of coefficients in the source circular element
    targ = el%id; src = c%id
    vi = real([(j,j=0,N-1)],DP)

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then
       ! common stuff between both branches below
       M = el%M
       allocate(r%LHS(M,(2*N-1)*(2-abs(el%ibnd))), r%RHS(M), &
            & cemat(M,0:N-1), semat(M,1:N-1))

       ! can the target element "see" the outside of the source element?
       if (dom%inclBg(src,targ)) then
          cemat = ce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
          semat = se(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))

          ! setup LHS (head effects due to source element)
          ! for constant head (-1), or matching (0)
          if (el%ibnd==0 .or. el%ibnd==-1) then
             allocate(mR(M,0:N-1,0:1),mR0(0:N-1,0:1))
             mR(1:M,0:N-1,0) = Ke(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             mR(1:M,1:N-1,1) = Ko(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             mR0(0:N-1,0) = Ke(e%parent%mat(idx), vi(0:N-1), e%r)
             mR0(1:N-1,1) = Ko(e%parent%mat(idx), vi(1:N-1), e%r)

             ! head effects on other element
             r%LHS(1:M,1:N) =       Kn(0:N-1,:)/spread(Kn0(0:N-1),1,M)*cmat/c%parent%K
             r%LHS(1:M,N+1:2*N-1) = Kn(1:N-1,:)/spread(Kn0(1:N-1),1,M)*smat/c%parent%K
             deallocate(Kn,Kn0)

             select case (c%ibnd)
             case(2)
                allocate(Kn0(0:1))
                kap = kappa(p,c%parent)
                Kn0(0:1) = bK(kap*c%r,2)

                if (c%StorIn) then
                   ! wellbore storage and skin from finite-radius well
                   r%LHS(1:M,1) = -Kn0(0)*((2.0 + c%r**2*c%dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                        & (Kn0(0)*c%r*p)/(2.0*PI*c%r*kap*Kn0(1)*c%parent%T))*r%LHS(1:M,1)
                   r%RHS(1:M) = 0.0
                else
                   ! specified flux (finite-radius well no storage)
                   r%RHS(1:M) = Kn0(0)*time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn0(1))*r%LHS(1:M,1)
                   r%LHS(1:M,1) = 0.0
                end if
             case(-1,0,1)
                ! do nothing, results already in LHS
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
  end function ellipse_match_head_other

  function ellipse_match_flux_other(e,el,dom,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, domain, matching, match_result
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(ellipse), intent(in) :: e ! source ellipse
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    type(match_result) :: r

    integer :: j, src, targ, N, M
    real(DP), allocatable :: cemat(:,:), semat(:,:)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: mR(:,:,:), dmR(:,:,:), mR0(:,:)
    complex(DP), allocatable :: dPot_dEta(:,:), dPot_dPsi(:,:), dPot_dX(:,:), dPot_dY(:,:)

    N = c%N ! number of coefficients in the source circular element
    targ = el%id; src = c%id
    vi = real([(j,j=0,N-1)],DP)

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then
       ! common stuff between both branches below
       M = el%M
       allocate(r%LHS(M,(2*N-1)*(2-abs(el%ibnd))), r%RHS(M), &
            & cemat(M,0:N-1), semat(M,1:N-1))
       cemat = ce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
       semat = se(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))

       ! can the target element "see" the outside of the source element?
       if (dom%inclBg(src,targ)) then

          ! setup LHS (flux effects due to source element)
          ! for constant head (-1), or matching (0)
          if (el%ibnd==0 .or. el%ibnd==-1) then

             ! flux effects of source ellipse on target element
             allocate(mR(M,0:N,0:1), dmR(M,0:N,0:1), mR0(0:N-1,0:1), &
                  & dPot_dEta(M,2*N-1), dPot_dPsi(M,2*N-1),&
                  & dPot_dX(M,2*N-1), dPot_dY(M,2*N-1))
             mR(1:M,0:N-1,0) =   Ke(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             mR(1:M,1:N-1,1) =   Ko(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             dmR(1:M,0:N-1,0) = dKe(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             dmR(1:M,1:N-1,1) = dKo(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             mR0(0:N-1,0) =      Ke(e%parent%mat(idx), vi(0:N-1), e%r)
             mR0(1:N-1,1) =      Ko(e%parent%mat(idx), vi(1:N-1), e%r)

             ! derivative wrt radius of source element             
             dPot_dEta(:,1:N) =       dmR(0:N-1,:,0)/spread(mR0(0:N-1,0),1,M)*cemat/c%parent%K
             dPot_dEta(:,N+1:2*N-1) = dmR(1:N-1,:,1)/spread(mR0(1:N-1,1),1,M)*semat/c%parent%K
             deallocate(dKn)

             select case (c%ibnd)
             case(2)
                ! radially-symmetric, angular derivative zero by default
                dPot_dPsi(:,:) = 0.0
             case(-1,0,1)
                ! derivative wrt angle of source element for more general circular elements
                dPot_dPsi(:,1:N) =       mR(0:N-1,:,0)/spread(mR0(0:N-1,0),1,M)* &
                     & dce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))/c%parent%K
                dPot_dPsi(:,N+1:2*N-1) = mR(1:N-1,:,1)/spread(mR0(1:N-1,1),1,M)* &
                     & dse(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))/c%parent%K
             end select

             ! project these from elliptical onto Cartesian coordinates
!!$             dPot_dX = dPot_dEta*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
!!$                     & dPot_dPsi*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
!!$
!!$             dPot_dY = dPot_dEta*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
!!$                     & dPot_dPsi*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

             ! project from Cartesian to "radial" coordinate of target element
             if (el%id <= dom%num(1)) then
                ! other element is a circle
                select case(el%ibnd)
                case(2)
                   if (el%StorIn) then
                      ! other element is a well with wellbore storage (Type III BC)
                      r%LHS(:,1:N) =       Kn(0:N-1,:)/spread(Kn0(0:N-1),1,M)*smat/c%parent%K
                      r%LHS(:,N+1:2*N-1) = Kn(1:N-1,:)/spread(Kn0(1:N-1),1,M)*cmat/c%parent%K

                      ! head effects of element
                      r%LHS = -(el%r*p/el%parent%T)*r%LHS

                      ! radial flux effects of element
                      r%LHS = r%LHS + (2.0 + el%r**2*el%dskin*p/el%parent%T)* &
                           & (dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                           &  dPot_dY*spread(sin(el%Pcm),2,2*N-1))
                   end if
                case(-1,0)
                   r%LHS = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                         & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
                end select
  
             else
                ! other element is an ellipse
                r%LHS = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                      & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)
             end if
             deallocate(Kn,Kn0)
          end if
       else
          ! can target element "see" the inside of the source element?
          !  i.e., is the source element the parent?

          if (el%ibnd==0 .or. el%ibnd==-1) then
             allocate(In(M,0:N-1),dIn(M,0:N),In0(0:N-1), &
                  & dPot_dR(M,2*N-1),dPot_dP(M,2*N-1),dPot_dX(M,2*N-1),dPot_dY(M,2*N-1))
             kap = kappa(p,c%element)
             call bId(kap*c%G(targ)%Rgm(1:M),N+1,In,dIn)
             In0(0:N-1) = bI(kap*c%r,N)
             dIn = kap*dIn

             ! derivative wrt radius of source element
             dPot_dR(:,1:N) =       dIn(0:N-1,:)/spread(In0(0:N-1),1,M)*cmat/c%K
             dPot_dR(:,N+1:2*N-1) = dIn(1:N-1,:)/spread(In0(1:N-1),1,M)*smat/c%K
             deallocate(dIn)

             ! derivative wrt angle of source element 
             dPot_dP(:,1:N) =      -In(0:N-1,:)*spread(vi(0:N-1)/In0(0:N-1),1,M)*smat/c%K
             dPot_dP(:,N+1:2*N-1) = In(1:N-1,:)*spread(vi(1:N-1)/In0(1:N-1),1,M)*cmat/c%K 

             ! project these from cylindrical onto Cartesian coordinates
             dPot_dX = dPot_dR*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
                     & dPot_dP*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

             dPot_dY = dPot_dR*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
                     & dPot_dP*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

             ! project from Cartesian to "radial" coordinate of target element
             if (el%id <= dom%num(1)) then
                ! other element is a circle
                select case(el%ibnd)
                case(2)
                   if (el%StorIn) then
                      ! other element is a well with wellbore storage (Type III BC)
                      r%LHS(:,1:N) =       In(0:N-1,:)/spread(In0(0:N-1),1,M)*smat/c%K
                      r%LHS(:,N+1:2*N-1) = In(1:N-1,:)/spread(In0(1:N-1),1,M)*cmat/c%K

                      ! head effects of element
                      r%LHS = -(el%r*p/el%T)*r%LHS

                      ! radial flux effects of element
                      r%LHS = r%LHS + (2.0 + el%r**2*el%dskin*p/el%T)* &
                           & (dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                           &  dPot_dY*spread(sin(el%Pcm),2,2*N-1))
                   end if
                case(-1,0)
                   r%LHS = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                         & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
                end select
  
             else
                ! other element is an ellipse
                r%LHS = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                      & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)
             end if
             deallocate(In,In0)
          end if
       end if
    end if

  end function ellipse_match_flux_other

end module elliptical_elements

