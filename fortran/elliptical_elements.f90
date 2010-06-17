! this module contains the basic functions defining the head or flux 
! effects of a elliptical element.

module elliptical_elements
  use constants, only : DP, PI
  use kappa_mod
  use time_mod

  implicit none
  private
  public :: ellipse_match

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
    integer, intent(in) :: idx ! indicates which value of p (global state)
    type(match_result) :: r

    integer :: j, N, M, lo, hi
    complex(DP), allocatable :: RMn(:,:), dRMn(:,:) ! mod. radial Mathieu function (K{e,o} or I{e,o})
    complex(DP) :: cemat(1:e%M,0:e%N-1,0:1), semat(1:e%M,1:e%N-1,0:1)
    integer, dimension(0:e%N-1) :: vi

    N = e%N; M = e%M
    vi = [(j,j=0,N-1)]

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

       r%LHS(lo:hi,1:N) =       spread(dRMn(0:N-1,0)/RMn(0:N-1,0), 1,M)*cemat(:,:,0) ! a_n flux
       r%LHS(lo:hi,N+1:2*N-1) = spread(dRMn(1:N-1,1)/RMn(1:N-1,1), 1,M)*semat(:,:,1) ! b_n flux
       
       if (e%ibnd == 0 .or. (e%ibnd == 1 .and. e%calcin)) then
          RMn(0:N-1,0) =   Ie(e%mat(idx), vi(0:N-1), e%r)
          RMn(1:N-1,1) =   Io(e%mat(idx), vi(1:N-1), e%r)
          dRMn(0:N-1,0) = dIe(e%mat(idx), vi(0:N-1), e%r)
          dRMn(1:N-1,1) = dIo(e%mat(idx), vi(1:N-1), e%r)

          r%LHS(lo:hi,2*N:3*N-1) = spread(dRMn(0:N-1,0)/RMn(0:N-1,0), 1,M)*cemat(:,:,0) ! c_n flux
          r%LHS(lo:hi,3*N:4*N-2) = spread(dRMn(1:N-1,1)/RMn(1:N-1,1), 1,M)*semat(:,:,1) ! d_n flux
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
       r%RHS(1:M) = sum(spread(line(e,p,idx),dim=1,ncopies=M)*r%LHS(1:M,1:N:2),dim=2)
       deallocate(r%LHS)
       allocate(r%LHS(1:M,0))
    end select
  end function ellipse_match_self

  function ellipse_match_other(e,el,dom,p,idx) result(r)
    use utility, only : outerprod
    use type_definitions, only : ellipse, domain, matching, match_result
    use mathieu_functions, only : ce, se, dce, dse, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    implicit none

    type(ellipse), intent(in) :: e ! source ellipse
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    type(match_result) :: r 

    integer :: j, src, targ, N, M, lo, hi, nmax
    complex(DP), allocatable :: cemat(:,:), semat(:,:), dcemat(:,:), dsemat(:,:) 
    integer, dimension(0:e%N-1) :: vi
    complex(DP), allocatable :: RMn(:,:,:), dRMn(:,:,:), RMn0(:,:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    real(DP), allocatable :: hsq(:,:)
    real(DP) :: K

    N = e%N ! number of coefficients in the source elliptical element
    targ = el%id; src = e%id
    vi = [(j,j=0,N-1)]

    if (dom%inclBg(src,targ) .or. dom%InclIn(src,targ)) then

       M = el%M
       if (el%ibnd == 0) then
          allocate(r%LHS(2*M,2*N-1), r%RHS(2*M))
          lo = M+1; hi = 2*M
       else
          allocate(r%LHS(M,2*N-1), r%RHS(M))
          lo = 1;  hi = M
       end if
       
       allocate(RMn(M,0:N-1,0:1), RMn0(0:N-1,0:1), cemat(1:M,0:N-1), semat(1:M,1:N-1))

       ! setup LHS 
       ! for matching or specified total head target elements
       if (el%ibnd == 0 .or. el%ibnd == -1) then

          if (dom%inclBg(src,targ)) then
             ! can the target element "see" the outside of the source element?
             ! use exterior angular and radial modified Mathieu functions
             cemat(:,0:N-1) =   ce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
             RMn(0:N-1,1:M,0) = Ke(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             RMn0(0:N-1,0) =    Ke(e%parent%mat(idx), vi(0:N-1), e%r)
             if (.not. e%ibnd == 2) then
                ! odd functions not needed for line sources
                semat(:,1:N-1) =   se(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
                RMn(1:N-1,1:M,1) = Ko(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
                RMn0(1:N-1,1) =    Ko(e%parent%mat(idx), vi(1:N-1), e%r)
                nmax = N
             end if
             K = e%parent%K
             nmax = 2*N-1
          else
             ! can target element "see" the inside of the source element?
             ! i.e., is the source element the parent?
             ! use interior angular and radial modified Mathieu functions
             cemat(:,0:N-1) =   ce(e%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
             semat(:,1:N-1) =   se(e%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
             RMn(0:N-1,1:M,0) = Ie(e%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             RMn(1:N-1,1:M,1) = Io(e%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             RMn0(0:N-1,0) =    Ie(e%mat(idx), vi(0:N-1), e%r)
             RMn0(1:N-1,1) =    Io(e%mat(idx), vi(1:N-1), e%r)
             K = e%K
             nmax = 2*N-1
          end if
          
          ! head effects on other element
          r%LHS(1:M,1:N) =       RMn(0:N-1,:,0)/spread(RMn0(0:N-1,0),1,M)*cemat/K ! a_n || c_n
          if (.not. e%ibnd == 2) then
             r%LHS(1:M,N+1:2*N-1) = RMn(1:N-1,:,1)/spread(RMn0(1:N-1,1),1,M)*semat/K ! b_n || d_n
          end if          
       end if
          
       ! for matching, specified total flux, or specified elemental flux target element
       if (el%ibnd == 0 .or. el%ibnd == +1 .or. el%ibnd == +2) then
          allocate(dRMn(M,0:N-1,0:1), dcemat(1:M,0:N-1), dsemat(1:M,1:N-1), &
               & dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), dPot_dX(M,2*N-1),dPot_dY(M,2*N-1))

          ! flux effects of source ellpise on target element
          if (dom%inclBg(src,targ)) then
             ! use exterior angular and radial modified mathieu functions
             if (.not. el%ibnd == 0) then
                cemat(:,0:N-1) =   ce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
                RMn(0:N-1,1:M,0) = Ke(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
                RMn0(0:N-1,0) =    Ke(e%parent%mat(idx), vi(0:N-1), e%r)
                if (.not. e%ibnd == 2) then
                   ! line sources don't need odd functions
                   semat(:,1:N-1) =   se(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
                   RMn(1:N-1,1:M,1) = Ko(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
                   RMn0(1:N-1,1) =    Ko(e%parent%mat(idx), vi(1:N-1), e%r)
                end if
                K = e%parent%K
             end if
             dcemat(:,0:N-1) =   dce(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
             dRMn(0:N-1,1:M,0) = dKe(e%parent%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             if (.not. e%ibnd == 2) then
                dsemat(:,1:N-1) =   dse(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
                dRMn(1:N-1,1:M,1) = dKo(e%parent%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             end if
          else
             ! use interior angular and radial modified Mathieu functions
             if (.not. el%ibnd == 0) then
                cemat(:,0:N-1) =   ce(e%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
                semat(:,1:N-1) =   se(e%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
                RMn(0:N-1,1:M,0) = Ie(e%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
                RMn(1:N-1,1:M,1) = Io(e%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
                RMn0(0:N-1,0) =    Ie(e%mat(idx), vi(0:N-1), e%r)
                RMn0(1:N-1,1) =    Io(e%mat(idx), vi(1:N-1), e%r)
             end if
             dcemat(:,0:N-1) =   dce(e%mat(idx), vi(0:N-1), e%G(targ)%Pgm(1:M))
             dsemat(:,1:N-1) =   dse(e%mat(idx), vi(1:N-1), e%G(targ)%Pgm(1:M))
             dRMn(0:N-1,1:M,0) = dIe(e%mat(idx), vi(0:N-1), e%G(targ)%Rgm(1:M))
             dRMn(1:N-1,1:M,1) = dIo(e%mat(idx), vi(1:N-1), e%G(targ)%Rgm(1:M))
             K = e%K
          end if

          ! derivative wrt radius of source element
          dPot_dR(1:M,1:N) =          dRMn(1:M,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat
          if (.not. e%ibnd == 2) then
             dPot_dR(1:M,N+1:2*N-1) = dRMn(1:M,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat
          end if

          ! derivative wrt angle of source element 
          dPot_dP(1:M,1:N) =          RMn(0:N-1,:,0)/spread(RMn0(0:N-1,0),1,M)*dcemat
          if (.not. e%ibnd == 2) then
             dPot_dP(1:M,N+1:2*N-1) = RMn(1:N-1,:,1)/spread(RMn0(1:N-1,1),1,M)*dsemat
          end if

          ! project these from elliptical onto Cartesian coordinates
          allocate(hsq(size(e%G(targ)%Rgm),2*N-1))
          ! squared metric factor -- less a common f
          hsq = spread(e%f/2.0*(cosh(2.0*e%G(targ)%Rgm) - cos(2.0*e%G(targ)%Pgm)),2,2*N-1) 
          dPot_dX = (dPot_dR*spread(sinh(e%G(targ)%Rgm)*cos(e%G(targ)%Pgm),2,2*N-1) - &
                   & dPot_dP*spread(cosh(e%G(targ)%Rgm)*sin(e%G(targ)%Pgm),2,2*N-1))/hsq
          dPot_dY = (dPot_dR*spread(cosh(e%G(targ)%Rgm)*sin(e%G(targ)%Pgm),2,2*N-1) + &
                   & dPot_dP*spread(sinh(e%G(targ)%Rgm)*cos(e%G(targ)%Pgm),2,2*N-1))/hsq
          deallocate(hsq)

          ! project from Cartesian to "radial" coordinate of target element
          if (el%id <= dom%num(1)) then
             ! other element is a circle
             if (el%ibnd == 2) then
                ! other element is a specified elemental flux circle
                if (el%StorIn) then
                   ! other element is a well with wellbore storage (Type III BC)
                   r%LHS(1:M,1:N) = RMn(0:N-1,:,0)/spread(RMn0(0:N-1,0),1,M)*semat/K
                   if (.not. e%ibnd == 2) then
                      r%LHS(1:M,N+1:2*N-1) = RMn(1:N-1,:,1)/spread(RMn0(1:N-1,1),1,M)*cemat/K
                   end if
                   
                   ! head effects of source ellipse
                   r%LHS(1:M,1:nmax) = -(el%r*p/el%parent%T)*r%LHS(1:M,1:nmax)
                   
                   ! radial flux effects of element
                   r%LHS(1:M,1:nmax) = r%LHS + (2.0 + el%r**2*el%dskin*p/el%parent%T)* &
                        & (dPot_dX*spread(cos(el%Pcm),2,nmax) + &
                        &  dPot_dY*spread(sin(el%Pcm),2,nmax))
                else
                   continue
                   ! wells without wellbore storage are specified and don't
                   ! need to have head/flux computed along their boundary
                end if
             else
                ! other element is a standard circle 
                r%LHS(lo:hi,1:nmax) = dPot_dX*spread(cos(el%Pcm),2,nmax) + &
                                    & dPot_dY*spread(sin(el%Pcm),2,nmax)
             end if
          else
             ! other element is a different ellipse
             r%LHS(lo:hi,:) = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,nmax) + &
                            & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,nmax)

          end if
          deallocate(dRMn,dPot_dR,dPot_dP,dPot_dX,dPot_dY)
       end if
       deallocate(RMn,RMn0)
    end if
    
    if (e%ibnd == 2) then
       ! sum line source effects and move to RHS, re-setting LHS to 0 unknowns
       r%RHS(1:hi) = sum(spread(line(e,p,idx),1,hi)*r%LHS(1:hi,1:N:2))
       deallocate(r%LHS)
       allocate(r%LHS(1:hi,0))
    end if
    
  end function ellipse_match_other
  
  function line(e,p,idx) result(a2n)
    ! this function returns the coefficients for a specified-flux line source
    use type_definitions, only : ellipse
    use mathieu_functions, only : Ke,dKe
    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    complex(DP), dimension(ceiling(e%N/2.0)) :: a2n ! only even coefficients of even order
    real(DP), dimension(0:e%ms-1) :: vs
    real(DP), dimension(1:e%ms,ceiling(e%N/2.0)) :: arg
    integer, dimension(0:e%ms-1) :: vi
    integer :: i, N, MS, nmax

    N = e%N; MS = e%ms
    nmax = ceiling(e%N/2.0)
    vi = [(i,i=0,MS-1)]  ! integer vector
    vs = (-1.0_DP)**vi(0:MS-1)    ! sign vector
    arg = spread(vs(0:MS-1)/real(1 - (2*vi(0:MS-1))**2,DP),dim=2,ncopies=nmax)
    
    ! factor of 4 different from Kuhlman&Neuman paper
    ! include Radial/dRadial MF here to balance with those in general solution
    a2n(1:nmax) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI)* &
            & Ke(e%parent%mat(idx), vi(0:N-1:2), e%r) / dKe(e%parent%mat(idx), vi(0:N-1:2), e%r)* &
            & (-vs(0:N-1:2))*sum(arg(:,:)*conjg(e%mat(idx)%A(1:MS,0:nmax-1,0)),dim=1)

  end function line

  function ellipse_calc(p,e,lo,hi,Rgp,Pgp,inside) result(H)
    use type_definitions, only : ellipse
    use mathieu_functions, only : Ke,Ko, Ie,Io, ce,se

    complex(DP), dimension(:), intent(in) :: p
    type(ellipse), intent(in) :: e
    integer, intent(in) :: lo,hi 
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(size(p,1)) :: H

    integer, dimension(e%N) :: vi
    complex(DP), dimension(size(p,1),e%N) :: aa,bb
    complex(DP), dimension(size(p,1),e%N,0:1) :: RMRgp,RMR0,AM
    integer :: n0, np, i, j, N

    N = e%N
    np = size(p,1)
    vi = [(i,i=0,N-1)]

    if (inside) then
       if (e%match) then
          n0 = 2*N ! inside of matching ellipse
       else
          n0 = 1   ! inside of specified boundary ellipse
       end if
       do i=1,np
          j = lo+i-1
          RMRgp(j,0:N-1,0) = Ie(e%mat(j),vi(0:N-1),Rgp)
          RMRgp(j,1:N-1,1) = Io(e%mat(j),vi(1:N-1),Rgp)
          RMR0(j,0:N-1,0) =  Ie(e%mat(j),vi(0:N-1),e%r)
          RMR0(j,1:N-1,1) =  Io(e%mat(j),vi(1:N-1),e%r)
          AM(j,0:N-1,0) =    ce(e%mat(j),vi(0:N-1),Pgp)
          AM(j,1:N-1,1) =    se(e%mat(j),vi(0:N-1),Pgp)
       end do
    else
       n0 = 1
       do i=1,np
          j = lo+i-1
          RMRgp(j,0:N-1,0) = Ke(e%parent%mat(j),vi(0:N-1),Rgp)
          RMRgp(j,1:N-1,1) = Ko(e%parent%mat(j),vi(1:N-1),Rgp)
          RMR0(j,0:N-1,0) =  Ke(e%parent%mat(j),vi(0:N-1),e%r)
          RMR0(j,1:N-1,1) =  Ko(e%parent%mat(j),vi(1:N-1),e%r)
          AM(j,0:N-1,0) =    ce(e%parent%mat(j),vi(0:N-1),Pgp)
          AM(j,1:N-1,1) =    se(e%parent%mat(j),vi(0:N-1),Pgp)
       end do
    end if
    aa(1:np,1:N) = e%coeff(lo:hi,n0:n0+N-1)
    bb(1:np,2:N) = e%coeff(lo:hi,n0+N:n0+2*N-2)
    H(1:np) = sum(RMRgp(:,0:N-1,0)/RMR0(:,0:N-1,0)*aa(:,1:N)*AM(:,0:N-1,0), 2) + &
            & sum(RMRgp(:,1:N-1,1)/RMR0(:,1:N-1,1)*bb(:,2:N)*AM(:,1:N-1,1), 2)
  end function ellipse_calc
end module elliptical_elements

