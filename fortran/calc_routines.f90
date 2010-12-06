! this module contains the routines that compute the head
! or flux at a calculation location, assuming the coefficients
! for each element have already been computed.

module calc_routines
  implicit none

  private
  public :: headCalc, velCalc

contains

  !##################################################
  function headCalc(Z,p,lo,hi,dom,c,e,bg) result(H)
    use type_definitions, only : element, domain, circle, ellipse
    use constants, only : DP
    use circular_elements, only : circle_calc
    use elliptical_elements, only : ellipse_calc
    use kappa_mod, only : kappa
    use time_mod, only : time
#ifdef DEBUG
    use utility, only : ccosh
    use constants, only : EYE, PI
#endif

    complex(DP), intent(in) :: Z  ! location for calculation (complex coordinates)
    complex(DP), dimension(:), intent(in) :: p  ! vector of Laplace parameters
    integer, intent(in) :: lo,hi  ! lo,hi bounds of p relative to overall s
    type(domain), intent(in) :: dom
    type(circle),  target, dimension(:), intent(in) :: c
    type(ellipse), target, dimension(:), intent(in) :: e
    type(element), intent(in) :: bg
    complex(DP), dimension(size(p,1)) :: H  

    real(DP), dimension(sum(dom%num)) :: Rgp, Pgp
    type(element), pointer :: elin => null()
    integer :: nc, ne, ntot, np, j
    integer :: in, oth
#ifdef DEBUG
    complex(DP) :: Ztmp
#endif

    call check_np(p,lo,hi)
    nc = dom%num(1)
    ne = dom%num(2)
    ntot = nc+ne
    np = size(p,dim=1)

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(Z,c,e,dom,Rgp,Pgp,in) 
   
#ifdef DEBUG
    if (any(Rgp < spacing(0.0) .or. Pgp < -PI .or. Pgp > PI)) then
       print *, 'headCalc error in results returned from calcLocation Rgp:',Rgp,' Pgp:',Pgp
       stop
    end if
#endif

    H(:) = 0.0

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
          if (dom%InclUp(j) == 0) then  ! circle is also in background
             H(1:np) = H(1:np) + &
                  & circle_calc(p,c(j),lo,hi,Rgp(j),Pgp(j),.false.)
          end if
       end do
       
       do j = nc+1,ntot
          if (dom%InclUp(j) == 0) then  ! ellipse is also in background
             H(1:np) = H(1:np) + &
                  & ellipse_calc(p,e(j-nc),lo,hi,Rgp(j),Pgp(j),.false.)
          end if
       end do
       H(1:np) = H(1:np)/bg%k ! convert to head

    !##################################################
    !! calculation point is inside an element (not background)
    else
       if (in <= nc) then
          ! calculation point is inside (or on bdry of) a circular element
          if (c(in)%calcin) then
             H(1:np) = H(1:np) + &
                  & circle_calc(p,c(in),lo,hi,Rgp(in),Pgp(in),.true.)
          end if
       else 
          ! calculation point is inside (or on bdry of) an elliptical element
          if (e(in-nc)%calcin) then
             H(1:np) = H(1:np) + &
                  & ellipse_calc(p,e(in-nc),lo,hi,Rgp(in),Pgp(in),.true.)
          end if
       end if

       ! effects of any other elements which may be inside same circle/ellipse too
       do oth = 1,nc
          ! oth element is a circle
          if (dom%InclIn(in,oth)) then
             H(1:np) = H(1:np) + &
                  & circle_calc(p,c(oth),lo,hi,Rgp(oth),Pgp(oth),.false.)
          end if
       end do
       do oth = nc+1,ntot
          ! other element is an ellipse
          if (dom%InclIn(in,oth)) then
             H(1:np) = H(1:np) + &
                  & ellipse_calc(p,e(oth-nc),lo,hi,Rgp(oth),Pgp(oth),.false.)
          end if
       end do
       
       if (in <= nc) then
          elin => c(in)%element     ! circle 
       else
          elin => e(in-nc)%element  ! ellipse
       end if
       
       ! apply potential source term on inside of element
       H(1:np) = H(:) - elin%areaQ*elin%Ss*time(p,elin%time,.true.)/kappa(p,elin)**2
       H(1:np) = H(:)/elin%K ! convert to head

       elin => null()

    end if
  end function headCalc

  !##################################################
  function velCalc(Z,p,lo,hi,dom,c,e,bg) result(v)
   
    use type_definitions, only : element, domain, circle, ellipse
    use constants, only : DP
    use circular_elements, only : circle_deriv
    use elliptical_elements, only : ellipse_deriv
    use kappa_mod, only : kappa
    use time_mod, only : time
    use utility, only : rotate_vel
#ifdef DEBUG
    use constants, only : PI
#endif

    complex(DP), intent(in) :: Z
    complex(DP), dimension(:), intent(in) :: p
    integer, intent(in) :: lo,hi  ! lo,hi bounds of p relative to overall s
    type(domain), intent(in) :: dom
    type(circle),  target, dimension(:), intent(in) :: c
    type(ellipse), target, dimension(:), intent(in) :: e
    type(element), intent(in) :: bg
    complex(DP), dimension(size(p,1),2) :: v, dH

    real(DP) :: hsq
    real(DP), dimension(sum(dom%num)) :: Rgp, Pgp
    type(element), pointer :: elin => null()
    integer :: in, oth
    integer :: np, nc, ne, ntot, j

    call check_np(p,lo,hi)
    nc = dom%num(1)
    ne = dom%num(2)
    ntot = nc+ne
    np = size(p,dim=1)

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(Z,c,e,dom,Rgp,Pgp,in) 

#ifdef DEBUG
    if (any(Rgp < spacing(0.0) .or. Pgp < -PI .or. Pgp > PI)) then
       print *, 'velCalc error in results returned from calcLocation Rgp:',Rgp,' Pgp:',Pgp
       stop
    end if
#endif

    v(1:np,1:2) = 0.0

    ! TODO
    ! there should be a way to combine the two branches of this 
    ! if statement into a more general single branch.

    !##################################################
    !! calculation point is outside all inclusions
    if (in == 0) then
       do j = 1,nc
          if (dom%InclUp(j) == 0) then  ! circle is also in background
             dH(1:np,1:2) = circle_deriv(p,c(j),lo,hi,Rgp(j),Pgp(j),.false.)
             ! project onto X and Y
             v(1:np,1) = v(:,1) + &
                  & cos(Pgp(j))*dH(:,1) - sin(Pgp(j))*dH(:,2)/Rgp(j) ! x
             v(1:np,2) = v(:,2) + &
                  & sin(Pgp(j))*dH(:,1) + cos(Pgp(j))*dH(:,2)/Rgp(j) ! y
          end if
       end do
       
       do j = nc+1,ntot
          if (dom%InclUp(j) == 0) then  ! ellipse is also in background
             
             ! eta and psi components of gradient
             dH(1:np,1:2) = ellipse_deriv(p,e(j-nc),lo,hi,Rgp(j),Pgp(j),.false.)
             dH(1:np,1:2) = rotate_vel(dH(:,1:2),e(j-nc)%theta)

             hsq = e(j-nc)%f/2.0_DP*(cosh(2.0_DP*Rgp(j)) - cos(2.0_DP*Pgp(j)))
             v(1:np,1) = v(:,1) + (sinh(Rgp(j))*cos(Pgp(j))*dH(:,1) - &
                                &  cosh(Rgp(j))*sin(Pgp(j))*dH(:,2))/hsq
             v(1:np,2) = v(:,2) + (cosh(Rgp(j))*sin(Pgp(j))*dH(:,1) + &
                                &  sinh(Rgp(j))*cos(Pgp(j))*dH(:,2))/hsq
          end if
       end do
       v(1:np,1:2) = -v(1:np,1:2)/bg%por ! return velocity (gradient points uphill)
       
    !##################################################
    else
       if (in <= nc) then
          !! calculation point is inside (or on bdry of) a circular element
          if (c(in)%calcin) then
             dH(1:np,1:2) = circle_deriv(p,c(in),lo,hi,Rgp(in),Pgp(in),.true.)
             v(1:np,1) = v(:,1) + &
                  & cos(Pgp(in))*dH(:,1) - sin(Pgp(in))*dH(:,2)/Rgp(in)
             v(1:np,2) = v(:,2) + &
                  & sin(Pgp(in))*dH(:,1) + cos(Pgp(in))*dH(:,2)/Rgp(in)
          end if
       else 
          ! calculation point is inside or on the boundary of an elliptical element
          if (e(in-nc)%calcin) then
             dH(1:np,1:2) = ellipse_deriv(p,e(in-nc),lo,hi,Rgp(in),Pgp(in),.true.)
             dH(1:np,1:2) = rotate_vel(dH(:,1:2),e(in-nc)%theta)
             
             hsq = e(in-nc)%f/2.0_DP*(cosh(2.0_DP*Rgp(in)) - cos(2.0_DP*Pgp(in)))
             v(1:np,1) = v(:,1) + (sinh(Rgp(in))*cos(Pgp(in))*dH(:,1) - &
                                &  cosh(Rgp(in))*sin(Pgp(in))*dH(:,2))/hsq
             v(1:np,2) = v(:,2) + (cosh(Rgp(in))*sin(Pgp(in))*dH(:,1) + &
                                &  sinh(Rgp(in))*cos(Pgp(in))*dH(:,2))/hsq

          end if
       end if
       ! effects of any other elements which may be inside same circle/ellipse too
       do oth = 1,nc
          ! other element is a circle
          if (dom%InclIn(in,oth)) then
             dH(1:np,1:2) = circle_deriv(p,c(oth),lo,hi,Rgp(oth),Pgp(oth),.false.)
             v(1:np,1) = v(:,1) + &
                  & cos(Pgp(oth))*dH(:,1) - sin(Pgp(oth))*dH(:,2)/Rgp(oth)
             v(1:np,2) = v(:,2) + &
                  & sin(Pgp(oth))*dH(:,1) + cos(Pgp(oth))*dH(:,2)/Rgp(oth)
          end if
       end do
       do oth = nc+1,ntot
          ! other element is an ellipse
          if (dom%InclIn(in,oth)) then
             dH(1:np,1:2) = ellipse_deriv(p,e(oth-nc),lo,hi,Rgp(oth),Pgp(oth),.false.)
             dH(1:np,1:2) = rotate_vel(dH(:,1:2),e(oth-nc)%theta) 

             hsq = e(oth-nc)%f/2.0_DP*(cosh(2.0_DP*Rgp(oth)) - cos(2.0_DP*Pgp(oth)))
             v(1:np,1) = v(:,1) + (sinh(Rgp(oth))*cos(Pgp(oth))*dH(:,1) - &
                                 & cosh(Rgp(oth))*sin(Pgp(oth))*dH(:,2))/hsq
             v(1:np,2) = v(:,2) + (cosh(Rgp(oth))*sin(Pgp(oth))*dH(:,1) + &
                                &  sinh(Rgp(oth))*cos(Pgp(oth))*dH(:,2))/hsq
          end if
       end do
       
       if (in <= nc) then
          elin => c(in)%element     ! circle 
       else
          elin => e(in-nc)%element  ! ellipse
       end if

       ! area source has no flux effects, since it is a constant WRT space
       v(1:np,1:2) = -v(1:np,1:2)/elin%por ! gradient points uphill

       elin => null()

    end if
  end function velCalc

  !##################################################
  subroutine calcLocation(Z,c,e,dom,Rgp,Pgp,inside)
    use constants, only : DP, EYE
    use type_definitions, only : circle, ellipse, domain
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

    integer :: nc, ne, ntot, j, first, second, third, k

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
       if (Rgp(j) < c(j)%r) then    ! inside (not including on boundary)
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
       if (Rgp(j+nc) < e(j)%r) then    ! inside (not including on boundary)
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
       case (2) ! inside two elements
          if (dom%InclUp(inout(1)) == 0) then
             inside = inout(2)
          else
             inside = inout(1)
          end if
       case (3) ! inside three elements
          do first = 1,3
             do second = 1,3
                if (second /= first) then
                   do third = 1,3
                      if ((third /= first) .and. (third /= second)) then
                         if (dom%InclIn(inout(first),inout(second))  .and. &
                              & dom%InclIn(inout(third),inout(first))) then 
                            inside = inout(second)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       case default
          write(*,*) 'CALCLOCATION ERROR: more than triply-nested inclusions', dom%InclIn
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

  subroutine check_np(p,lo,hi)
    use constants, only : DP
    complex(DP), dimension(:), intent(in) :: p
    integer, intent(in) :: lo,hi
    if (hi-lo+1 /= size(p,1)) then
       write(*,'(A,3(1X,I0))') 'ERROR: lo,hi do not match dimensions of p',lo,hi,size(p,1)
       stop
    end if
  end subroutine check_np
  
end module calc_routines
