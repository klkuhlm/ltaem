! this module contains the routines that compute the head
! or flux at a calculation location, assuming the coefficients
! for each element have already been computed.

module unsat_calc_routines
  implicit none

  private
  public :: headCalc, velCalc

contains

  !##################################################
  function headCalc(Z,dom,c,e,bg) result(H)
    use unsat_type_definitions, only : element, domain, circle, ellipse
    use constants, only : DP
    use unsat_circular_elements, only : circle_calc
    use unsat_elliptical_elements, only : ellipse_calc

#ifdef DEBUG
    use utility, only : ccosh
    use constants, only : EYE, PI
#endif

    complex(DP), intent(in) :: Z  ! location for calculation (complex coordinates)
    type(domain), intent(in) :: dom
    type(circle),  target, dimension(:), intent(in) :: c
    type(ellipse), target, dimension(:), intent(in) :: e
    type(element), intent(in) :: bg
    complex(DP) :: H  

    type(element), pointer :: elin => null()
    real(DP), dimension(sum(dom%num)) :: Rgp, Pgp
    integer :: nc, ne, ntot, j, in

#ifdef DEBUG
    complex(DP) :: Ztmp
#endif

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = nc+ne

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(Z,c,e,dom,Rgp,Pgp,in) 
   
#ifdef DEBUG
    if (any(Rgp < epsilon(0.0) .or. Pgp < -PI .or. Pgp > PI)) then
       print *, 'headCalc error in results returned from calcLocation Rgp:',Rgp,' Pgp:',Pgp
       stop
    end if
#endif

    H = 0.0

#ifdef DEBUG
    ! units opened in ltaem_main.f90
    ! which element is each calculation point located in?
    write(303,*) real(Z),aimag(Z),in
    ! write cartesian vectors pointing from center of each element to
    ! calc point for visual checking in gnuplot
    do j = 1,nc
       write(404,*) c(j)%x,c(j)%y,Rgp(j)*cos(Pgp(j)),Rgp(j)*sin(Pgp(j))
    end do
    do j = 1,ne
       Ztmp = e(j)%f*ccosh(cmplx(Rgp(j+nc),Pgp(j+nc),DP))*exp(EYE*e(j)%theta)
       write(404,*) e(j)%x,e(j)%y,real(Ztmp),aimag(Ztmp)
    end do
    write(404,'(/)')
#endif

    ! TODO there should be a way to combine the two branches of this 
    ! if statement into a more general single branch.

!!$    print '(A,I0)', 'in ==',in
    !##################################################
    !! calculation point is outside all elements (in background)
    if (in == 0) then
       do j = 1,nc
          H = H + circle_calc(c(j),Rgp(j),Pgp(j),.false.)
       end do       
       do j = nc+1,ntot
          H = H + ellipse_calc(e(j-nc),Rgp(j),Pgp(j),.false.)
       end do
       ! z is "down" coordinate in problem
       ! z = -y (y is imaginary part of Z)
       H = H*exp(-bg%alpha*aimag(Z)/2.0) ! convert to head
       
!!$       H = H + bg%qz0/bg%alpha

    !##################################################
    !! calculation point is inside an element (not background)
    else
       if (in <= nc) then
          elin => c(in)%element
          ! calculation point is inside (or on bdry of) a circular element
          if (c(in)%calcin) then
             H = circle_calc(c(in),Rgp(in),Pgp(in),.true.)
          else
             H = 0.0 ! inside an element labeled no-calc
          end if
       else 
          elin => e(in-nc)%element
          ! calculation point is inside (or on bdry of) an elliptical element
          if (e(in-nc)%calcin) then
             H = ellipse_calc(e(in-nc),Rgp(in),Pgp(in),.true.)
          else
             H = 0.0 ! inside an element labeled no-calc
          end if
       end if
              
       ! apply potential source term on inside of element
       ! TODO : handle unsaturated source terms???

       H = H*exp(-elin%alpha*aimag(Z)/2.0) ! convert to head
       elin => null()
    end if
  end function headCalc

  !##################################################
  function velCalc(Z,dom,c,e,bg) result(v)
   
    use unsat_type_definitions, only : element, domain, circle, ellipse
    use constants, only : DP
    use unsat_circular_elements, only : circle_deriv
    use unsat_elliptical_elements, only : ellipse_deriv
    use utility, only : rotate_vel
#ifdef DEBUG
    use constants, only : PI
#endif

    complex(DP), intent(in) :: Z
    type(domain), intent(in) :: dom
    type(circle),  target, dimension(:), intent(in) :: c
    type(ellipse), target, dimension(:), intent(in) :: e
    type(element), intent(in) :: bg
    complex(DP), dimension(2) :: v
    complex(DP), dimension(1,2) :: dH

    type(element), pointer :: elin
    real(DP) :: hsq
    real(DP), dimension(sum(dom%num)) :: Rgp, Pgp
    integer :: nc, ne, ntot, j, in

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = nc+ne

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(Z,c,e,dom,Rgp,Pgp,in) 

#ifdef DEBUG
    if (any(Rgp < epsilon(0.0) .or. Pgp < -PI .or. Pgp > PI)) then
       print *, 'velCalc error in results returned from calcLocation Rgp:',Rgp,' Pgp:',Pgp
       stop
    end if
#endif

    v(1:2) = 0.0

    ! TODO
    ! there should be a way to combine the two branches of this 
    ! if statement into a more general single branch.

    !##################################################
    !! calculation point is outside all inclusions
    if (in == 0) then
       do j = 1,nc
          dH(1,1:2) = circle_deriv(c(j),Rgp(j),Pgp(j),.false.)
          ! project onto X and Y
          v(1) = v(1) + cos(Pgp(j))*dH(1,1) - sin(Pgp(j))*dH(1,2)/Rgp(j) ! x
          v(2) = v(2) + sin(Pgp(j))*dH(1,1) + cos(Pgp(j))*dH(1,2)/Rgp(j) ! y
       end do
       
       do j = nc+1,ntot
          ! eta and psi components of gradient
          dH(1,1:2) = ellipse_deriv(e(j-nc),Rgp(j),Pgp(j),.false.)
          dH(:,1:2) = rotate_vel(dH(:,1:2),e(j-nc)%theta)

          hsq = e(j-nc)%f/2.0_DP*(cosh(2.0_DP*Rgp(j)) - cos(2.0_DP*Pgp(j)))
          v(1) = v(1) + (sinh(Rgp(j))*cos(Pgp(j))*dH(1,1) - &
                      &  cosh(Rgp(j))*sin(Pgp(j))*dH(1,2))/hsq
          v(2) = v(2) + (cosh(Rgp(j))*sin(Pgp(j))*dH(1,1) + &
                      &  sinh(Rgp(j))*cos(Pgp(j))*dH(1,2))/hsq
       end do
       ! return velocity (gradient points uphill)
       ! convert to head
       v(1:2) = -v(1:2)*exp(-bg%alpha*aimag(Z)/2.0)
       
    !##################################################
    else
       if (in <= nc) then
          elin => c(in)%element
          !! calculation point is inside (or on bdry of) a circular element
          if (c(in)%calcin) then
             dH(1,1:2) = circle_deriv(c(in),Rgp(in),Pgp(in),.true.)
             v(1) = cos(Pgp(in))*dH(1,1) - sin(Pgp(in))*dH(1,2)/Rgp(in)
             v(2) = sin(Pgp(in))*dH(1,1) + cos(Pgp(in))*dH(1,2)/Rgp(in)
          end if
       else 
          elin => e(in-nc)%element
          ! calculation point is inside or on the boundary of an elliptical element
          if (e(in-nc)%calcin) then
             dH(1,1:2) = ellipse_deriv(e(in-nc),Rgp(in),Pgp(in),.true.)
             dH(:,1:2) = rotate_vel(dH(:,1:2),e(in-nc)%theta)
             
             hsq = e(in-nc)%f/2.0_DP*(cosh(2.0_DP*Rgp(in)) - cos(2.0_DP*Pgp(in)))
             v(1) = (sinh(Rgp(in))*cos(Pgp(in))*dH(1,1) - &
                  &  cosh(Rgp(in))*sin(Pgp(in))*dH(1,2))/hsq
             v(2) = (cosh(Rgp(in))*sin(Pgp(in))*dH(1,1) + &
                  &  sinh(Rgp(in))*cos(Pgp(in))*dH(1,2))/hsq

          end if
       end if
       
       v(1:2) = -v(1:2)*exp(-elin%alpha*aimag(Z)/2.0)
       elin => null()

    end if
  end function velCalc

  !##################################################
  subroutine calcLocation(Z,c,e,dom,Rgp,Pgp,inside)
    use constants, only : DP, EYE
    use unsat_type_definitions, only : circle, ellipse, domain
    use utility, only : ccosh, cacosh

    complex(DP), intent(in) :: Z
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(domain), intent(in) :: dom
    real(DP), dimension(sum(dom%num)), intent(out) :: Rgp, Pgp
    integer, intent(out) :: inside

    complex(DP), dimension(sum(dom%num)) :: Zgp
    complex(DP), dimension(dom%num(2)) :: Ztmp
    integer, dimension(sum(dom%num)) :: inout

    integer :: nc, ne, ntot, j, k

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = nc+ne

    ! determine if observation point is inside an inclusion
    inout(1:ntot) = 0
    k = 0

    ! components of vector from center of circle to observation point
    Zgp(1:nc) = Z - cmplx(c(:)%x,c(:)%y,DP)
    Rgp(1:nc) = abs(Zgp(1:nc))
    Pgp(1:nc) = atan2(aimag(Zgp(1:nc)), real(Zgp(1:nc)))
    do j = 1, nc
       if (Rgp(j) <= c(j)%r) then    ! inside or on boundary
          k = k+1
          inout(k) = j
       end if
    end do

    ! components of vector from center of ellipse to observation point
    Zgp(nc+1:ntot) = Z - cmplx(e(:)%x,e(:)%y,DP) 
    Ztmp(1:ne) = cacosh( Zgp(nc+1:ntot)*exp(-EYE*e(:)%theta)/e(:)%f )
    Rgp(nc+1:ntot) =  real(Ztmp(1:ne))
    Pgp(nc+1:ntot) = aimag(Ztmp(1:ne))
    do j = 1, ne
       if (Rgp(j+nc) <= e(j)%r) then    ! inside or on boundary
          k = k+1
          inout(k) = j+nc
       end if
    end do

    if (ntot > 0) then
       select case (k)
       case (0) ! inside nothing (background)
          inside = 0
       case (1) ! inside one element
          inside = inout(1)
       case default
          write(*,*) 'unsat CALCLOCATION ERROR: nested inclusions not allowed', dom%InclIn
          stop 
       end select
    else
       ! no matching elements
       inside = 0
       Rgp = -999.
       Pgp = -999.
    end if

    ! move observation points inside ibnd==2 elements 
    where (Rgp(1:nc) < c%r .and. c%ibnd == 2)
       Rgp(1:nc) = c%r
    end where
    where (Rgp(nc+1:ntot) < e%r .and. e%ibnd == 2)
       Rgp(nc+1:ntot) = e%r
    end where

  end subroutine  calcLocation
  
end module unsat_calc_routines
