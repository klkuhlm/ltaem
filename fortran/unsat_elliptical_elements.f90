! this module contains the basic functions defining the head or flux 
! effects of a elliptical element.

module unsat_elliptical_elements
  implicit none

  private
  public :: ellipse_match, ellipse_calc, ellipse_deriv, line

  interface ellipse_match
     module procedure ellipse_match_self, ellipse_match_other
  end interface

contains
  function ellipse_match_self(e) result(r)
    use constants, only : DP, PI
    use unsat_type_definitions, only : ellipse, match_result, element
    use mathieu_functions, only : ce, se, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    use utility, only : ynot 
    implicit none

    type(ellipse), target, intent(in) :: e
    type(match_result) :: r
    type(element), pointer :: f => null()  ! shortcut for parent of element 

    integer :: j, N, M, loM, hiM, ierr, nrows, ncols
    ! mod. radial Mathieu function (K{e,o} or I{e,o})
    complex(DP), allocatable :: RMn(:,:), dRMn(:,:)
    complex(DP), dimension(1:e%M, 0:e%N-1, 0:1) :: cemat
    complex(DP), dimension(1:e%M, 1:e%N-1, 0:1) :: semat
    integer, dimension(0:e%N-1) :: vi
    complex(DP), allocatable :: H(:,:),dH(:,:)

    N = e%N
    M = e%M
    vi(0:N-1) = [(j,j=0,N-1)]

    if (e%ibnd == 0) then
       nrows = 2*M
       ncols = 4*N-2
       loM = M+1
       hiM = 2*M
    elseif (e%ibnd == 2) then
       ! specified flux line source has no unknowns
       ! only appears on RHS of other elements
       nrows = 0
       ncols = 2*N-1
       loM = 1
       hiM = M
    else
       nrows = M
       ncols = 2*N-1
       loM = 1
       hiM = M
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows), stat=ierr)
    if (ierr /= 0) stop 'elliptical_elements.f90:ellipse_match_self()&
         & error allocating: r%LHS, r%RHS'

    if (e%ibnd /= 2) then
       f => e%parent

       ! outside ang. mod. Mathieu fcns (last dimension is inside/outside)
       cemat(1:M,0:N-1,0) = transpose(ce(f%mat, vi(0:N-1), e%Pcm(:)))
       semat(1:M,1:N-1,0) = transpose(se(f%mat, vi(1:N-1), e%Pcm(:)))
       if (e%ibnd == 0) then
          ! inside
          cemat(1:M,0:N-1,1) = transpose(ce(e%mat, vi(0:N-1), e%Pcm(:)))
          semat(1:M,1:N-1,1) = transpose(se(e%mat, vi(1:N-1), e%Pcm(:)))
       end if

       allocate(h(M,ncols),dH(M,ncols),stat=ierr)
       if (ierr /= 0) stop 'unsat_elliptical_elements.f90 error allocating: h,dh'

       ! setup LHS
       H(:,1:N) =       cemat(1:M,0:N-1,0) ! a_n head
       H(:,N+1:2*N-1) = semat(1:M,1:N-1,0) ! b_n head
       if (e%ibnd == 0) then
          H(:,2*N:3*N-1) = -cemat(1:M,0:N-1,1) ! c_n head
          H(:,3*N:4*N-2) = -semat(1:M,1:N-1,1) ! d_n head
       end if
       
       ! matching or specified total head
       if (e%ibnd == 0 .or. e%ibnd == -1) then
          r%LHS(1:M,1:2*N-1) = H(:,1:2*N-1)* &
               & spread(exp(f%alpha*e%f/2.0*sinh(e%r)*sin(e%Pcm(1:M))),2,2*N-1)

          if (e%ibnd == 0) then
             r%LHS(1:M,2*N:4*N-2) = H(:,2*N:4*N-2)*&
                  & spread(exp(e%alpha*e%f/2.0*sinh(e%r)*sin(e%Pcm(1:M))),2,2*N-1)
          end if
       end if

       ! matching or specified flux
       if (e%ibnd == 0 .or. e%ibnd == +1) then
          allocate(RMn(0:N-1,0:1),dRMn(0:N-1,0:1), stat=ierr)
          if (ierr /= 0) stop 'elliptical_elemen ts.f90 error allocating: RMn, dRMn'

          ! radial functions last dimension is even/odd
          RMn(0:N-1,0) =   Ke(f%mat, vi(0:N-1), e%r) ! even fn
          RMn(1:N-1,1) =   Ko(f%mat, vi(1:N-1), e%r) ! odd fn
          RMn(0,1) = 0.0
          dRMn(0:N-1,0) = dKe(f%mat, vi(0:N-1), e%r) ! even deriv
          dRMn(1:N-1,1) = dKo(f%mat, vi(1:N-1), e%r) ! odd deriv
          dRMn(0,1) = 0.0

          dH(:,1:N) =       spread(dRMn(0:N-1,0)/RMn(0:N-1,0),1,M)*cemat(:,0:N-1,0) ! a_n flux
          dH(:,N+1:2*N-1) = spread(dRMn(1:N-1,1)/RMn(1:N-1,1),1,M)*semat(:,1:N-1,0) ! b_n flux

          r%LHS(loM:hiM,1:2*N-1) = (dH(:,1:2*N-1) + H(:,1:2*N-1)*&
               & spread(e%f/2.0*cosh(e%r)*sin(e%Pcm(:))*f%alpha,2,2*N-1))*&
               & spread(exp(e%f/2.0*sinh(e%r)*sin(e%Pcm(:))*f%alpha),2,2*N-1)

          if (e%ibnd == 0) then
             RMn(0:N-1,0) =   Ie(e%mat, vi(0:N-1), e%r)
             RMn(1:N-1,1) =   Io(e%mat, vi(1:N-1), e%r)
             RMn(0,1) = 0.0
             dRMn(0:N-1,0) = dIe(e%mat, vi(0:N-1), e%r)
             dRMn(1:N-1,1) = dIo(e%mat, vi(1:N-1), e%r)
             dRMn(0,1) = 0.0

             dH(:,1:N) =       -spread(dRMn(0:N-1,0)/RMn(0:N-1,0),1,M)*cemat(:,0:N-1,1) ! c_n flux
             dH(:,N+1:2*N-1) = -spread(dRMn(1:N-1,1)/RMn(1:N-1,1),1,M)*semat(:,1:N-1,1) ! d_n flux

             r%LHS(loM:hiM,2*N:4*N-2) = (dH(:,1:2*N-1) + H(:,1:2*N-1)*&
               & spread(e%f/2.0*cosh(e%r)*sin(e%Pcm(:))*e%alpha,2,2*N-1))*&
               & spread(exp(e%f/2.0*sinh(e%r)*sin(e%Pcm(:))*e%alpha),2,2*N-1)

          end if
          deallocate(RMn,dRMn,stat=ierr)
          if (ierr /= 0) stop 'elliptical_elements.f90 error deallocating: RMn,dRMn'
       end if

       ! setup RHS
       select case(e%ibnd)
       case(-1)
          ! put specified head on RHS
          r%RHS(1:M) = e%bdryQ
       case(0)
          ! put constant area source term effects on RHS
          ! TODO: handle unsaturated area sources???
          r%RHS(1:2*M) = 0.0_DP  ! no effects
       case(1)
          ! put specified flux effects on RHS
          r%RHS(1:M) = e%bdryQ/ynot(e%r,e%f)
       end select

       r%RHS(1:M) = r%RHS(:)*exp(e%f/2.0*sinh(e%r)*sin(e%Pcm(:))*f%alpha)

       f => null()
    end if
  end function ellipse_match_self

  function ellipse_match_other(e,el,dom) result(r)
    use constants, only : DP
    use unsat_type_definitions, only : ellipse, domain, matching, match_result, element
    use mathieu_functions, only : ce, se, dce, dse, Ke, Ko, dKe, dKo, Ie, Io, dIe, dIo
    use utility, only : rotate_vel_mat
    implicit none

    type(ellipse), target, intent(in) :: e ! source ellipse
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    type(match_result) :: r
    type(element), pointer :: f => null()

    integer :: j, src, targ, N, M, loM, hiM, ierr, nrows, ncols
    complex(DP), allocatable :: cemat(:,:), semat(:,:), dcemat(:,:), dsemat(:,:) 
    integer, dimension(0:e%N-1) :: vi
    complex(DP), allocatable :: RMn(:,:,:), dRMn(:,:,:), RMn0(:,:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    real(DP), allocatable :: hsq(:,:)

    complex(DP), allocatable :: H(:,:), dH(:,:)

    N = e%N ! number of coefficients in the source elliptical element
    targ = el%id
    src = e%id
    vi(0:N-1) = [(j,j=0,N-1)]

    M = el%M
    ! target element determines number of rows
    if (el%ibnd == 0) then
       nrows = 2*M
       loM = M+1
       hiM = 2*M
    elseif (el%ibnd == 2) then
       nrows = 0
    else
       nrows = M
       loM = 1
       hiM = M
    end if

    ! source element determines number of columns
    if (e%ibnd == 0) then
       ncols = 4*N-2
    else
       ncols = 2*N-1
       ! line sources (ibound==2) are later re-allocated to zero size LHS
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows), stat=ierr)
    if (ierr /= 0) stop 'elliptical_elements.f90:ellipse_match_other() error allocating&
         &r%LHS, r%RHS'
    r%LHS = 0.0
    r%RHS = 0.0

    if (nrows > 0) then
       f => e%parent

       allocate(RMn(1:M,0:N-1,0:1), RMn0(0:N-1,0:1), cemat(1:M,0:N-1), semat(1:M,1:N-1), &
            & H(M,2*N-1), dH(M,2*N-1), stat=ierr)
       if (ierr /= 0) stop 'elliptical_elements.f90:ellipse_match_other() error allocating:&
            &RMn, RMn0, cemat, semat'

       ! can the target element "see" the outside of the source element?
       ! use exterior angular and radial modified Mathieu functions
       cemat(1:M,0:N-1) = transpose(ce(f%mat, vi(0:N-1), e%G(targ)%Pgm(:)))
       semat(1:M,1:N-1) = transpose(se(f%mat, vi(1:N-1), e%G(targ)%Pgm(:)))
       RMn(1:M,0:N-1,0) = transpose(Ke(f%mat, vi(0:N-1), e%G(targ)%Rgm(:)))
       RMn(1:M,1:N-1,1) = transpose(Ko(f%mat, vi(1:N-1), e%G(targ)%Rgm(:)))
       RMn(1:M,0,1) = 0.0
       RMn0(0:N-1,0) = Ke(e%parent%mat, vi(0:N-1), e%r)
       RMn0(1:N-1,1) = Ko(e%parent%mat, vi(1:N-1), e%r)
       RMn0(0,1) = 0.0
       ! odd functions not needed for line source, but computed anyway
       
       H(:,1:N) =       RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1) ! a_n
       H(:,N+1:2*N-1) = RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,1:N-1) ! b_n

       ! setup LHS 
       ! for matching or specified total head target elements
       if (el%ibnd == 0 .or. el%ibnd == -1) then

          ! head effects due to ellipse on outside other element
          r%LHS(1:M,1:2*N-1) = H(:,1:2*N-1)* &
               & spread(exp(f%alpha*e%f/2.0*sinh(e%G(targ)%Rgm(:))*&
               & sin(e%G(targ)%Pgm(:))),2,2*N-1)
       elseif (e%ibnd == 2) then        
          ! case handled at end of routine
       end if

       ! for matching, specified total flux, or specified elemental flux target element
       if (el%ibnd == 0 .or. el%ibnd == +1) then
          allocate(dRMn(M,0:N-1,0:1), dcemat(1:M,0:N-1), dsemat(1:M,1:N-1), &
               & dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), dPot_dX(M,2*N-1),dPot_dY(M,2*N-1), stat=ierr)
          if (ierr /= 0) stop 'elliptical_elements.f90 error allocating:&
               & dRMn, dcemat, dsemat, dPot_dR, dPot_dP, dPot_dX, dPot_dY'

          ! flux effects of source ellpise on target element
          ! use exterior angular and radial modified mathieu functions
          dcemat(1:M,0:N-1) = transpose(dce(f%mat, vi(0:N-1), e%G(targ)%Pgm(:)))
          dRMn(1:M,0:N-1,0) = transpose(dKe(f%mat, vi(0:N-1), e%G(targ)%Rgm(:)))
          dsemat(1:M,1:N-1) = transpose(dse(f%mat, vi(1:N-1), e%G(targ)%Pgm(:)))
          dRMn(1:M,1:N-1,1) = transpose(dKo(f%mat, vi(1:N-1), e%G(targ)%Rgm(:)))
          dRMn(1:M,0,1) = 0.0
          
          ! derivative wrt radius of source element
          dH(1:M,1:N) =       dRMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*cemat(:,0:N-1)
          dH(1:M,N+1:2*N-1) = dRMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*semat(:,1:N-1)

          dPot_dR(:,1:2*N-1) = (dH(:,1:2*N-1) + H(:,1:2*N-1)*spread(&
               & f%alpha*e%f/2.0*cosh(e%G(targ)%Rgm(:))*sin(e%G(targ)%Pgm(:)),2,2*N-1))* &
               & spread(exp(f%alpha*e%f/2.0*sinh(e%G(targ)%Rgm(:))*&
               & sin(e%G(targ)%Pgm(:))),2,2*N-1)

          ! derivative wrt angle of source element 
          dH(1:M,1:N) =       RMn(:,0:N-1,0)/spread(RMn0(0:N-1,0),1,M)*dcemat(:,0:N-1)
          dH(1:M,N+1:2*N-1) = RMn(:,1:N-1,1)/spread(RMn0(1:N-1,1),1,M)*dsemat(:,1:N-1)

          dPot_dP(:,1:2*N-1) = (dH(:,1:2*N-1) + H(:,1:2*N-1)*spread(&
               & f%alpha*e%f/2.0*sinh(e%G(targ)%Rgm(:))*cos(e%G(targ)%Pgm(:)),2,2*N-1))* &
               & spread(exp(f%alpha*e%f/2.0*sinh(e%G(targ)%Rgm(:))*&
               & sin(e%G(targ)%Pgm(:))),2,2*N-1)

          ! project these from elliptical onto Cartesian coordinates
          allocate(hsq(size(e%G(targ)%Rgm),2*N-1), stat=ierr)
          if (ierr /= 0) stop 'elliptical_elements.f90 error allocating: hsq'

          ! squared metric factor -- less a common f
          hsq = spread(e%f/2.0*(cosh(2.0*e%G(targ)%Rgm) - cos(2.0*e%G(targ)%Pgm)),2,2*N-1) 
          dPot_dX = (dPot_dR*spread(sinh(e%G(targ)%Rgm)*cos(e%G(targ)%Pgm),2,2*N-1) - &
                   & dPot_dP*spread(cosh(e%G(targ)%Rgm)*sin(e%G(targ)%Pgm),2,2*N-1))/hsq
          dPot_dY = (dPot_dR*spread(cosh(e%G(targ)%Rgm)*sin(e%G(targ)%Pgm),2,2*N-1) + &
                   & dPot_dP*spread(sinh(e%G(targ)%Rgm)*cos(e%G(targ)%Pgm),2,2*N-1))/hsq

          ! rotate to compensate for potentially arbitrary source ellipse
          call rotate_vel_mat(dPot_dX,dPot_dY,e%theta)

          deallocate(hsq,stat=ierr)
          if (ierr /= 0) stop 'elliptical_elements.f90 error deallocating: hsq'

          ! project from Cartesian to "radial" coordinate of target element
          if (el%id <= dom%num(1)) then
             ! other element is a standard circle 
             r%LHS(loM:hiM,1:2*N-1) = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                                    & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
          else
             ! rotate to compensate for arbitrary target ellipse
             call rotate_vel_mat(dPot_dX,dPot_dY,-el%theta)

             ! other element is a different ellipse
             r%LHS(loM:hiM,1:2*N-1) = dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                                    & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)

          end if
          deallocate(dRMn,dcemat,dsemat,dPot_dR,dPot_dP,dPot_dX,dPot_dY,stat=ierr)
          if (ierr /= 0) stop 'elliptical_elements.f90, error deallocating:&
               & dRMn, dcemat, dsemat, dPot_dR, dPot_dP, dPot_dX, dPot_dY'

       end if
       deallocate(RMn,RMn0,cemat,semat,stat=ierr)
       if (ierr /= 0) stop 'elliptical_elements.f90:elliptical_match_other(), &
            &error deallocating: RMn, RMn0, cemat, semat'

       if (e%ibnd == 2) then
          ! sum line source effects and move to RHS, re-setting LHS to 0 unknowns
          ! only uses even-order even coefficients (~1/4 of 2N-1)
          r%RHS(1:hiM) = -sum(spread(line(e),1,hiM)*r%LHS(1:hiM,1:N:2))
          deallocate(r%LHS,stat=ierr)

          if (ierr /= 0) then
             stop 'elliptical_elements.f90 error deallocating: r%LHS'
          end if

          allocate(r%LHS(1:hiM,0),stat=ierr)
          if (ierr /= 0) then
             stop 'elliptical_elements.f90 error re-allocating: r%LHS'
          end if

       end if
       f => null()
    end if

  end function ellipse_match_other

  function line(e) result(a2n)
    ! this function returns the coefficients for a specified-flux line source
    use constants, only : DP, PI
    use unsat_type_definitions, only : ellipse
    use mathieu_functions, only : Ke,dKe

    type(ellipse), intent(in) :: e
    complex(DP), dimension(ceiling(e%N/2.0)) :: a2n ! only even coefficients of even order
    real(DP), dimension(0:e%ms-1) :: vs
    real(DP), dimension(1:e%ms,ceiling(e%N/2.0)) :: arg
    integer, dimension(0:e%ms-1) :: vi
    integer :: i, N, MS, nmax

    N = e%N 
    MS = e%ms
    nmax = ceiling(e%N/2.0)
    vi(0:MS-1) = [(i,i=0,MS-1)]  ! integer vector
    vs = -1.0 ! sign vector
    where (mod(vi,2)==0) vs = 1.0

    arg(1:MS,1:nmax) = spread(vs(0:MS-1)/real(1-(2*vi(0:MS-1))**2,DP),2,nmax)
    
    ! factor of 4 different from Kuhlman&Neuman paper
    ! include Radial/dRadial MF here to balance with those in general solution
    a2n(1:nmax) = e%bdryQ/(2.0*PI)* &
            & Ke(e%parent%mat, vi(0:N-1:2), e%r) / dKe(e%parent%mat, vi(0:N-1:2), e%r)* &
            & (-vs(0:N-1:2))*sum(arg(1:MS,1:nmax)*conjg(e%parent%mat%A(1:MS,0:nmax-1,0)),dim=1)

  end function line
  
  function ellipse_calc(e,Rgp,Pgp,inside) result(H)
    use constants, only : DP
    use unsat_type_definitions, only : ellipse
    use mathieu_functions, only : Ke,Ko, Ie,Io, ce,se

    type(ellipse), intent(in) :: e
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP) :: H

    integer, dimension(0:e%N-1) :: vi
    complex(DP), dimension(0:e%N-1) :: aa,bb
    complex(DP), dimension(0:e%N-1,0:1) :: RMRgp,RMR0,AM
    integer :: n0, i, N

    N = e%N
    vi(0:N-1) = [(i,i=0,N-1)]

    if (inside) then
       n0 = 2*N ! inside of matching ellipse

       RMRgp(0:N-1,0) = Ie(e%mat,vi(0:N-1),Rgp)
       RMRgp(1:N-1,1) = Io(e%mat,vi(1:N-1),Rgp)
       
       RMR0(0:N-1,0) =  Ie(e%mat,vi(0:N-1),e%r)
       RMR0(1:N-1,1) =  Io(e%mat,vi(1:N-1),e%r)

       AM(0:N-1,0) =    ce(e%mat,vi(0:N-1),Pgp)
       AM(1:N-1,1) =    se(e%mat,vi(1:N-1),Pgp)
    else
       n0 = 1
       RMRgp(0:N-1,0) = Ke(e%parent%mat,vi(0:N-1),Rgp)
       RMRgp(1:N-1,1) = Ko(e%parent%mat,vi(1:N-1),Rgp)

       RMR0(0:N-1,0) =  Ke(e%parent%mat,vi(0:N-1),e%r)
       RMR0(1:N-1,1) =  Ko(e%parent%mat,vi(1:N-1),e%r)

       AM(0:N-1,0) =    ce(e%parent%mat,vi(0:N-1),Pgp)
       AM(1:N-1,1) =    se(e%parent%mat,vi(1:N-1),Pgp)
    end if

    aa(0:N-1) = e%coeff(n0:n0+N-1)
    bb(1:N-1) = e%coeff(n0+N:n0+2*N-2)

    H = sum(RMRgp(0:N-1,0)/RMR0(0:N-1,0)*aa(0:N-1)*AM(0:N-1,0)) + &
      & sum(RMRgp(1:N-1,1)/RMR0(1:N-1,1)*bb(1:N-1)*AM(1:N-1,1))

  end function ellipse_calc

  function ellipse_deriv(e,Rgp,Pgp,inside) result(dH)
    use constants, only : DP
    use unsat_type_definitions, only : ellipse
    use mathieu_functions, only : Ke,Ko, Ie,Io, ce,se, dKe,dKo, dIe,dIo, dce,dse

    type(ellipse), intent(in) :: e
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(2) :: dH ! dPot_dEta, dPot_dPsi
    complex(DP) :: H, dEta, dPsi

    integer, dimension(0:e%N-1) :: vi
    complex(DP), dimension(0:e%N-1) :: aa,bb
    complex(DP), dimension(0:e%N-1,0:1) :: RMRgp,RMR0,AM,dRMRgp,dAM
    integer :: n0, i, N
    complex(DP) :: alpha
    
    N = e%N
    vi(0:N-1) = [(i,i=0,N-1)]

    if (inside) then
       n0 = 2*N ! inside of matching ellipse

       ! radial Mathieu functions -> to calc point
       RMRgp(0:N-1,0) =   Ie(e%mat,vi(0:N-1),Rgp)
       RMRgp(1:N-1,1) =   Io(e%mat,vi(1:N-1),Rgp)
       dRMRgp(0:N-1,0) = dIe(e%mat,vi(0:N-1),Rgp)
       dRMRgp(1:N-1,1) = dIo(e%mat,vi(1:N-1),Rgp)
       ! radial Mathieu functions -> to bdry of ellipse
       RMR0(0:N-1,0) =    Ie(e%mat,vi(0:N-1),e%r)
       RMR0(1:N-1,1) =    Io(e%mat,vi(1:N-1),e%r)
       ! angular mathieu funcitons -> to calc point
       AM(0:N-1,0) =      ce(e%mat,vi(0:N-1),Pgp)
       AM(1:N-1,1) =      se(e%mat,vi(1:N-1),Pgp)
       dAM(0:N-1,0) =    dce(e%mat,vi(0:N-1),Pgp)
       dAM(1:N-1,1) =    dse(e%mat,vi(1:N-1),Pgp)
       alpha = e%alpha
    else
       n0 = 1
       RMRgp(0:N-1,0) =   Ke(e%parent%mat,vi(0:N-1),Rgp)
       RMRgp(1:N-1,1) =   Ko(e%parent%mat,vi(1:N-1),Rgp)
       dRMRgp(0:N-1,0) = dKe(e%parent%mat,vi(0:N-1),Rgp)
       dRMRgp(1:N-1,1) = dKo(e%parent%mat,vi(1:N-1),Rgp)

       RMR0(0:N-1,0) =    Ke(e%parent%mat,vi(0:N-1),e%r)
       RMR0(1:N-1,1) =    Ko(e%parent%mat,vi(1:N-1),e%r)

       AM(0:N-1,0) =      ce(e%parent%mat,vi(0:N-1),Pgp)
       AM(1:N-1,1) =      se(e%parent%mat,vi(1:N-1),Pgp)
       dAM(0:N-1,0) =    dce(e%parent%mat,vi(0:N-1),Pgp)
       dAM(1:N-1,1) =    dse(e%parent%mat,vi(1:N-1),Pgp)
       alpha = e%parent%alpha
    end if

    aa(0:N-1) = e%coeff(n0:n0+N-1)
    bb(1:N-1) = e%coeff(n0+N:n0+2*N-2)

    H = sum(RMRgp(0:N-1,0)/RMR0(0:N-1,0)*aa(0:N-1)*AM(0:N-1,0)) + &
      & sum(RMRgp(1:N-1,1)/RMR0(1:N-1,1)*bb(1:N-1)*AM(1:N-1,1))

    dEta = sum(dRMRgp(0:N-1,0)/RMR0(0:N-1,0)*aa(0:N-1)*AM(0:N-1,0)) + &
         & sum(dRMRgp(1:N-1,1)/RMR0(1:N-1,1)*bb(1:N-1)*AM(1:N-1,1))

    dPsi = sum(RMRgp(0:N-1,0)/RMR0(0:N-1,0)*aa(0:N-1)*dAM(0:N-1,0)) + &
         & sum(RMRgp(1:N-1,1)/RMR0(1:N-1,1)*bb(1:N-1)*dAM(1:N-1,1))

    dH(1) = (dEta + H*e%f/2.0*cosh(Rgp)*sin(Pgp)*alpha)
    dH(2) = (dEta + H*e%f/2.0*sinh(Rgp)*cos(Pgp)*alpha)

  end function ellipse_deriv
end module unsat_elliptical_elements

