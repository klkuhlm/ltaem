! this module contains the basic functions defining the head or flux 
! effects of a circular element.

module unsat_circular_elements
  implicit none

  private
  public :: circle_match, circle_calc, circle_deriv, well

  interface circle_match
     module procedure circle_match_self, circle_match_other
  end interface
    
contains
  function circle_match_self(c) result(r)
    use constants, only : DP, PI
    use utility, only : outer
    use unsat_type_definitions, only : circle, match_result, element
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(circle), target, intent(in) :: c
    type(match_result) :: r
    type(element), pointer :: f => null()

    integer :: j, N, M, loM, hiM, ierr, nrows, ncols
    complex(DP), allocatable :: Bn(:), dBn(:) ! mod. bessel function (K or I)
    complex(DP) :: kap
    real(DP) :: cmat(1:c%M,0:c%N-1), smat(1:c%M,1:c%N-1)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: H(:,:), dH(:,:)

    N = c%N
    M = c%M
    vi(0:N-1) = real([(j,j=0,N-1)],DP)

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
       nrows = 0
       ncols = 0
    else
       nrows = M
       ncols = 2*N-1
       ! here only flux matching, no head matching, so loM:him = 1:M
       loM = 1
       hiM = M
    end if

    allocate(r%LHS(nrows,ncols), r%RHS(nrows), stat=ierr)
    if (ierr /= 0) stop 'circular_elements.f90:circle_match_self() error &
         &allocating: r%LHS, r%RHS'

    if (c%ibnd /= 2) then
       f => c%parent

       cmat(1:M,0:N-1) = cos(outer(c%Pcm(1:M),vi(0:N-1)))
       smat(1:M,1:N-1) = sin(outer(c%Pcm(1:M),vi(1:N-1)))

       allocate(H(M,ncols),dH(M,ncols),stat=ierr)
       if (ierr /= 0) stop 'unsat_circular_elements.f90 error allocating: h,dh'

       ! setup LHS

       ! head needed for  both head (Type I) and flux (Type III) matching
       H(:,1:N) =       cmat(:,0:N-1) ! a_n head
       H(:,N+1:2*N-1) = smat(:,1:N-1) ! b_n head
       if (c%ibnd == 0) then
          H(:,2*N:3*N-1) = -cmat(:,0:N-1) ! c_n head
          H(:,3*N:4*N-2) = -smat(:,1:N-1) ! d_n head
       end if

       ! matching or specified total head
       if (c%ibnd == 0 .or. c%ibnd == -1) then
          r%LHS(1:M,1:2*N-1) = H(:,1:2*N-1)* &
               & spread(exp(-f%alpha*c%r*sin(c%Pcm(1:M))/2.0),2,2*N-1)
          if (c%ibnd == 0) then
             r%LHS(1:M,2*N:4*N-2) = H(:,2*N:4*N-2)* &
                  & spread(exp(-c%alpha*c%r*sin(c%Pcm(1:M))/2.0),2,2*N-1)
          end if
       end if

       ! matching or specified total flux
       if (c%ibnd == 0 .or. c%ibnd == +1 .or. c%ibnd == 2) then
          allocate(Bn(0:N-1), dBn(0:N-1), stat=ierr)
          if (ierr /= 0) stop 'circular_elements.f90 error allocating: Bn,dBn'
          kap = (f%alpha/2.0)**2 
          call dBK(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
          dBn(0:N-1) = kap*dBn(0:N-1)

          dH(:,1:N) =       spread(dBn(0:N-1)/Bn(0:N-1),1,M)*cmat(:,0:N-1) ! a_n flux
          dH(:,N+1:2*N-1) = spread(dBn(1:N-1)/Bn(1:N-1),1,M)*smat(:,1:N-1) ! b_n flux

          ! flux matching becomes type-III BC due to exponential transformation
          r%LHS(loM:hiM,1:2*N-1) = (dH(:,1:2*N-1) - H(:,1:2*N-1)*&
               & spread(sin(c%Pcm(1:M))*f%alpha/2.0,2,2*N-1))*&
               & spread(exp(-f%alpha*c%r*sin(c%Pcm(1:M))/2.0),2,2*N-1)

          if (c%ibnd == 0) then
             kap = (c%alpha/2.0)**2
             call dBI(kap*c%r,N,Bn(0:N-1),dBn(0:N-1))
             dBn(0:N-1) = kap*dBn(0:N-1)

             dH(:,2*N:3*N-1) = -spread(dBn(0:N-1)/Bn(0:N-1),1,M)*cmat(:,0:N-1) ! c_n flux
             dH(:,3*N:4*N-2) = -spread(dBn(1:N-1)/Bn(1:N-1),1,M)*smat(:,1:N-1) ! d_n flux

             r%LHS(loM:hiM,1:2*N-1) = (dH(:,1:2*N-1) - H(:,1:2*N-1)*&
                  & spread(sin(c%Pcm(1:M))*c%alpha/2.0,2,2*N-1))*&
                  & spread(exp(-c%alpha*c%r*sin(c%Pcm(1:M))/2.0),2,2*N-1)

          end if
          deallocate(Bn,dBn, stat=ierr)
          if (ierr /= 0) stop 'circular_elements.f90 error deallocating: Bn,dBn'

       end if

       ! setup RHS
       select case(c%ibnd)
       case(-1)
          ! put specified head out _outside_ of element on RHS
          r%RHS(1:M) = c%bdryQ
       case(0)
          ! TODO : handle unsaturated area sources???
          r%RHS(1:2*M) = 0.0 ! currently no effects
       case(1)
          ! put specified flux on _outside_ of element on RHS
          ! TODO : check addition of aquifer thickness to denominator
          r%RHS(1:M) = c%bdryQ/(2.0*PI*c%r)
       end select

       ! factor to convert from pressure (little psi) to helmholtz variable (big psi)
       r%RHS(1:M) = r%RHS(:)*exp(-f%alpha/2.0*c%r*sin(c%Pcm(1:M)))

       f => null()
    end if
  end function circle_match_self

  function circle_match_other(c,el,dom) result(r)
    use constants, only : DP, PI
    use utility, only : outer, rotate_vel_mat
    use unsat_type_definitions, only : circle, domain, matching, match_result, element
    use bessel_functions, only : bK, bI, dbK, dbI
    implicit none

    type(circle), target, intent(in) :: c ! source circle
    type(matching), intent(in) :: el ! target element (circle or ellipse)
    type(domain), intent(in) :: dom
    type(match_result) :: r
    type(element), pointer :: f => null()

    integer :: j, src, targ, N, M, loM, hiM, ierr, nrows, ncols
    real(DP), allocatable :: cmat(:,:), smat(:,:)
    real(DP), dimension(0:c%N-1) :: vi
    complex(DP), allocatable :: Bn(:,:), dBn(:,:), Bn0(:)
    complex(DP), allocatable :: dPot_dR(:,:), dPot_dP(:,:), dPot_dX(:,:), dPot_dY(:,:)
    complex(DP) :: kap
    complex(DP), allocatable :: H(:,:),dH(:,:)

    N = c%N ! number of coefficients in the source circular element
    targ = el%id
    src = c%id
    vi(0:N-1) = real([(j,j=0,N-1)],DP)

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
    if (c%ibnd == 0) then
       ncols = 4*N-2
    elseif (c%ibnd == 2) then
       ! compute effects due to well for RHS
       ncols = 1  ! reset to zero at end of routine
    else
       ncols = 2*N-1
    end if
    
    allocate(r%LHS(nrows,ncols), r%RHS(nrows), H(M,2*N-1), dH(M,2*N-1), stat=ierr)
    if (ierr /= 0) stop 'circular_elements.f90:circle_match_other() error allocating:&
            &r%LHS, r%RHS'
    r%LHS = 0.0
    r%RHS = 0.0

    if (nrows > 0) then
       f => c%parent

       allocate(Bn(1:M,0:N-1),Bn0(0:N-1),cmat(1:M,0:N-1),smat(1:M,1:N-1), stat=ierr)
       if (ierr /= 0) stop 'circular_elements.f90:circle_match_other() error allocating:&
            &Bn, Bn0, cmat, smat'

       cmat(1:M,0:N-1) = cos(outer(c%G(targ)%Pgm(1:M), vi(0:N-1)))
       smat(1:M,1:N-1) = sin(outer(c%G(targ)%Pgm(1:M), vi(1:N-1)))

       ! setup LHS 
       kap = (f%alpha/2.0)**2
       Bn(1:M,0:N-1) = bK(kap*c%G(targ)%Rgm(1:M),N)
       Bn0(0:N-1) =    bK(kap*c%r,N)

       H(:,1:N) =       Bn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(:,0:N-1) ! a_n
       H(:,N+1:2*N-1) = Bn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(:,1:N-1) ! b_n

       ! $$$$$$$$$$ head effects of source (c) on target (el) $$$$$$$$$$
       ! for matching or specified total head target elements
       if (el%ibnd == 0 .or. el%ibnd == -1) then

          ! head effects on outside of other element
          r%LHS(1:M,1:2*N-1) = H(:,1:2*N-1)* &
               & spread(exp(-f%alpha*c%G(targ)%Rgm(1:M)*&
               & sin(c%G(targ)%Pgm(1:M))/2.0),2,2*N-1)

       elseif (c%ibnd == 2) then
          ! specified flux (finite-radius well no storage)
          ! save head effects of well onto RHS
          r%RHS(1:M) = -well(c)*r%LHS(1:M,1)
          r%LHS(1:M,1) = 0.0 ! LHS matrix re-allocated to zero size below
       end if

       ! $$$$$$$$$$ flux effects of source (c) on target (el) $$$$$$$$$$
       ! for matching, specified total flux
       if (el%ibnd == 0 .or. el%ibnd == +1) then
          allocate(dBn(M,0:N-1), dPot_dR(M,2*N-1), dPot_dP(M,2*N-1), &
               & dPot_dX(M,2*N-1), dPot_dY(M,2*N-1), stat=ierr)
          if (ierr /= 0) stop 'circular_elements.f90 error allocating:&
               & dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY'

          ! flux effects of source circle on target element
          ! use exterior Bessel functions (Kn)
          kap = (f%alpha/2.0)**2
          call dBK(kap*c%G(targ)%Rgm(1:M),N,Bn(1:M,0:N-1),dBn(1:M,0:N-1))
          dBn(1:M,0:N-1) = kap*dBn(1:M,0:N-1)
          Bn0(0:N-1) = bK(kap*c%r,N)

          ! derivative wrt radius of source element
          dH(:,1:N) =       dBn(:,0:N-1)/spread(Bn0(0:N-1),1,M)*cmat(1:M,0:N-1)
          dH(:,N+1:2*N-1) = dBn(:,1:N-1)/spread(Bn0(1:N-1),1,M)*smat(1:M,1:N-1)

          dPot_dR(1:M,1:2*N-1) = (dH(:,1:2*N-1) - H(:,1:2*N-1)*&
               & spread(sin(c%G(targ)%Pgm(:))*f%alpha/2.0,2,2*N-1))*&
               & spread(exp(-f%alpha*c%G(targ)%Rgm(:)*&
               & sin(c%G(targ)%Pgm(:))/2.0),2,2*N-1)    

          ! derivative wrt angle of source element
          dH(:,1) = 0.0
          dH(:,2:N) =      -Bn(:,1:N-1)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*smat(:,1:N-1)
          dH(:,N+1:2*N-1) = Bn(:,1:N-1)*spread(vi(1:N-1)/Bn0(1:N-1),1,M)*cmat(:,1:N-1)

          dPot_dP(:,1:2*N-1) = (dH(:,1:2*N-1) - H(:,1:2*N-1)*spread(&
               & c%G(targ)%Rgm(:)*cos(c%G(targ)%Pgm(:))*f%alpha/2.0,2,2*N-1))*&
               & spread(exp(-f%alpha*c%G(targ)%Rgm(:)*&
               & sin(c%G(targ)%Pgm(:))/2.0),2,2*N-1)

          ! project these from cylindrical onto Cartesian coordinates
          dPot_dX(1:M,1:2*N-1) = dPot_dR*spread(cos(c%G(targ)%Pgm),2,2*N-1) - &
                               & dPot_dP*spread(sin(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)
          dPot_dY(1:M,1:2*N-1) = dPot_dR*spread(sin(c%G(targ)%Pgm),2,2*N-1) + &
                               & dPot_dP*spread(cos(c%G(targ)%Pgm)/c%G(targ)%Rgm,2,2*N-1)

          ! project from Cartesian to "radial" coordinate of target element
          if (el%id <= dom%num(1)) then
             ! other element is a  circular element without wellbore storage
             r%LHS(loM:hiM,1:2*N-1) = dPot_dX*spread(cos(el%Pcm),2,2*N-1) + &
                  & dPot_dY*spread(sin(el%Pcm),2,2*N-1)
          else
             ! rotate to allow for arbitrary oriented ellipse
             call rotate_vel_mat(dPot_dX,dPot_dY,-el%theta)

             ! other element is an ellipse
             r%LHS(loM:hiM,1:2*N-1) = &
                  & dPot_dX*spread(el%f*sinh(el%r)*cos(el%Pcm(1:M)),2,2*N-1) + &
                  & dPot_dY*spread(el%f*cosh(el%r)*sin(el%Pcm(1:M)),2,2*N-1)
          end if

          deallocate(dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY, stat=ierr)
          if (ierr /= 0) stop 'circular_elements.f90 error deallocating:&
               & dBn, dPot_dR, dPot_dP, dPot_dX, dPot_dY'

       end if
       deallocate(Bn,Bn0,cmat,smat, stat=ierr)
       if (ierr /= 0) stop 'circular_elements.f90 error deallocating: Bn,Bn0,cmat,smat'

       f => null()
    end if
       
    if (c%ibnd == 2) then
       
       ! source on _outside_ of target element
       r%RHS(loM:hiM) = -well(c)*r%LHS(loM:hiM,1)

       if (el%ibnd == 0) then
          hiM = 2*M
       else
          hiM = M
       end if

       deallocate(r%LHS,stat=ierr)
       if (ierr /= 0) stop 'circular_elements:circle_match_self() error deallocating r%LHS'
       allocate(r%LHS(hiM,0),stat=ierr)
       if (ierr /= 0) stop 'circular_elements:circle_match_self() error re-allocating r%LHS'
    end if

  end function circle_match_other

  function well(c) result(a0)
    ! this function returns the a_0 coefficient for a simple "well"
    use constants, only : DP, PI
    use unsat_type_definitions, only : circle
    use bessel_functions, only : bK
    type(circle), intent(in) :: c
    complex(DP) :: a0
    complex(DP), dimension(0:1) ::Kn
    
    Kn(0:1) = bK((c%parent%alpha/(2.0,0.0))**2 *c%r,2)
    a0 = Kn(0)*c%bdryQ/(2.0*PI*c%r*Kn(1))
  end function well
  
  function circle_calc(c,Rgp,Pgp,inside) result(H)
    use constants, only : DP
    use unsat_type_definitions, only : circle
    use bessel_functions, only : bK, bI

    type(circle), intent(in) :: c
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP) :: H

    real(DP), dimension(0:c%N-1) :: vr
    complex(DP), dimension(0:c%N-1) :: aa,bb,BRgp,BR0
    integer :: n0, i, N
    complex(DP) :: kap

    N = c%N
    vr(0:N-1) = real([(i,i=0,N-1)],DP)

    if (inside) then
       n0 = 2*N ! inside of matching circle 
       kap = (c%alpha/2.0)**2
       BRgp(0:N-1) = bI(Rgp*kap,N)
       BR0(0:N-1) =  bI(c%r*kap,N)
       
    else
       n0 = 1 
       kap = (c%parent%alpha/2.0)**2
       BRgp(0:N-1) = bK(Rgp*kap,N)
       BR0(0:N-1) =  bK(c%r*kap,N)
    end if

    ! if inside >> a_n is 2N:3N-1, b_n is 3N:4N-2
    ! if outside >> a_n is 1:N, b_n is N+1:2N-1 

    aa(0:N-1) = c%coeff(n0:n0+N-1) ! a_n or c_n
    bb(0) = 0.0 ! make odd/even conformable
    bb(1:N-1) = c%coeff(n0+N:n0+2*N-2) ! b_n or d_n

    H = sum(BRgp(0:N-1)/BR0(0:N-1)* &
         & ( aa(0:N-1)*cos(vr(0:N-1)*Pgp) + &
         &   bb(0:N-1)*sin(vr(0:N-1)*Pgp) ))

  end function circle_calc

  function circle_deriv(c,Rgp,Pgp,inside) result(dH)
    use constants, only : DP
    use unsat_type_definitions, only : circle
    use bessel_functions, only : bK, bI, dbk, dbi

    type(circle), intent(in) :: c
    real(DP), intent(in) :: Rgp, Pgp
    logical, intent(in) :: inside
    complex(DP), dimension(2) :: dH ! dPot_dR, dPot_dTheta
    complex(DP) :: H, dR, dTh

    real(DP), dimension(0:c%N-1) :: vr
    complex(DP), dimension(0:c%N-1) :: aa,bb,BRgp,BR0,dBRgp
    integer :: n0, i, N
    complex(DP) :: kap

    N = c%N
    vr(0:N-1) = real([(i,i=0,N-1)],DP)

    if (inside) then
       n0 = 2*N ! inside of matching circle
       kap = (c%alpha/2.0)**2
       call dBI(Rgp*kap,N,BRgp(0:N-1),dBRgp(0:N-1))
       dBRgp(0:N-1) = kap*dBRgp(:)
       BR0(0:N-1) = bI(c%r*kap,N)
    else
       n0 = 1
       kap = (c%parent%alpha/2.0)**2
       call dBK(Rgp*kap,N,BRgp(0:N-1),dBRgp(0:N-1))
       dBRgp(0:N-1) = kap*dBRgp(:)
       BR0(0:N-1) =  bK(c%r*kap,N)
    end if

    aa(0:N-1) = c%coeff(n0:n0+N-1)
    bb(0) = 0.0 ! make odd/even conformable
    bb(1:N-1) = c%coeff(n0+N:n0+2*N-2)

    H = circle_calc(c,Rgp,Pgp,inside)

    ! dPot_dR
    dR = sum(dBRgp(0:N-1)/BR0(0:N-1)* &  
         & ( aa(0:N-1)*cos(vr(0:N-1)*Pgp) + &
         &   bb(0:N-1)*sin(vr(0:N-1)*Pgp) ))

    ! dPot_dTheta
    dTh = sum(BRgp(0:N-1)/BR0(0:N-1)*vr(0:N-1)* &  
         & ( bb(0:N-1)*cos(vr(0:N-1)*Pgp) - &
         &   aa(0:N-1)*sin(vr(0:N-1)*Pgp) ))    

    dH(1) = (dR - H*sin(Pgp)*c%alpha/2.0)*&
         & exp(-c%alpha*Rgp*sin(Pgp)/2.0)

    dH(2) = (dTh - H*Rgp*cos(Pgp)*c%alpha/2.0)*&
         & exp(-c%alpha*Rgp*sin(Pgp)/2.0)

  end function circle_deriv
end module unsat_circular_elements

