! some general? functions collected in one place

module utility
  implicit none

  private
  public :: diag, logspace, linspace, outer, ccosh, cacosh, ynot, rotate_vel, rotate_vel_mat

  interface diag
     module procedure diagonal_z, diagonal_d
  end interface
  
  interface logspace
     module procedure logspace_int, logspace_r, logspace_d
  end interface
  
  interface outer
     module procedure outerprod_r, outerprod_d, outerprod_z, outerprod_dz, outerprod_zd
  end interface
  
  interface outerdiff
     module procedure outerdiff_r
  end interface outerdiff
  
contains

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  pure function outerprod_d(a,b) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: a,b
    real(DP), dimension(size(a),size(b)) :: c
    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  end function outerprod_d

  pure function outerprod_r(a,b) result(c)
    real, intent(in), dimension(:) :: a,b
    real, dimension(size(a),size(b)) :: c
    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  end function outerprod_r

  pure function outerprod_z(a,b) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: a,b
    complex(DP), dimension(size(a),size(b)) :: c
    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  end function outerprod_z

  pure function outerprod_dz(ca,db) result(c)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: ca
    real(DP), intent(in), dimension(:) :: db
    complex(DP), dimension(size(ca),size(db)) :: c
    c = spread(ca,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(ca))
  end function outerprod_dz

  pure function outerprod_zd(da,cb) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da
    complex(DP), intent(in), dimension(:) :: cb
    complex(DP), dimension(size(da),size(cb)) :: c
    c = spread(da,dim=2,ncopies=size(cb))*spread(cb,dim=1,ncopies=size(da))
  end function outerprod_zd

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
  pure function outerdiff_r(da,db) result(dc) 
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da, db
    real(DP), dimension(size(da),size(db)) :: dc
    c = spread(da,dim=2,ncopies=size(da))-spread(db,dim=1,ncopies=size(db))
  end function outerdiff_r
  
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
  pure elemental function ccosh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    x = real(z)
    y = aimag(z)
    f = cmplx(cosh(x)*cos(y), sinh(x)*sin(y),DP)
  end function ccosh

  pure elemental function cacosh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f
    
    ! compared ~10^-15 against acosh() in Matlab
    ! branch cut is left of +1, along x axis
    ! -pi <= aimag(f) <= +pi
    if(real(z) >= 0.0) then
       f = log(z + sqrt(z**2 - 1.0))
    else
       f = -log(z + sqrt(z**2 - 1.0))
    end if
  end function cacosh

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  pure function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: rnum, range, sgn

    rnum = real(num - 1,DP)
    range = abs(hi - lo) 
    sgn = sign(1.0_DP,hi-lo) ! if lo > high, count backwards
    forall (i=0:num-1) v(i+1) = lo + sgn*real(i,DP)*range/rnum
  end function linspace

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  pure function logspace_int(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)
  end function logspace_int

  pure function logspace_r(lo,hi,num) result(v)
    use constants, only : DP
    real, intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    v = logspace_int(int(lo),int(hi),num)
  end function logspace_r

  pure function logspace_d(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    v = logspace_int(int(lo),int(hi),num)
  end function logspace_d

   ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   pure function diagonal_z(A,n) result(d)
     use constants, only : DP
     integer, intent(in) :: n
     complex(DP), intent(in), dimension(:,:) :: A
     complex(DP), dimension(size(A,dim=1)-abs(n)) :: d
     integer :: i
     if (n >= 0) then
        ! diagonal in lower triangle (+) or main diag (0)
        forall (i=1:size(A,dim=1)-n) d(i) = A(i+n,i)
     else
        ! diagonal in upper triangle (-)
        forall (i=1:size(A,dim=1)+n) d(i) = A(i,i-n)
     end if
   end function diagonal_z
   
   pure function diagonal_d(A,n) result(d)
     use constants, only : DP
     integer, intent(in) :: n
     real(DP), intent(in), dimension(:,:) :: A
     real(DP), dimension(size(A,dim=1)-abs(n)) :: d
     integer :: i
     if (n >= 0) then
        ! diagonal in lower triangle (+) or main diag (0)
        forall (i=1:size(A,dim=1)-n) d(i) = A(i+n,i)
     else
        ! diagonal in upper triangle (-)
        forall (i=1:size(A,dim=1)+n) d(i) = A(i,i-n)
     end if
   end function diagonal_d

  ! return the approximate circumference of an ellipse
  ! given the radius in elliptical coordinates, and 
  ! the semi-focal distance (YNOT formula)
  elemental function ynot(eta,f) result(P)
    use constants, only : DP, LN2, LNPIOV2
    real(DP), intent(in) :: eta,f
    real(DP) :: P, y

    ! semi-major/minor length a:=f*cosh/sinh(psi)
    y = LN2/LNPIOV2
    P = 4.0_DP*((f*cosh(eta))**y + (f*sinh(eta))**y)**(1.0_DP/y)

  end function ynot

  function rotate_vel(v,theta) result(w)
    use constants, only : DP
    complex(DP), dimension(:,:), intent(in) :: v
    real(DP), intent(in) :: theta

    complex(DP), dimension(size(v,dim=1),2) :: w
    real(DP), dimension(2,2) :: rot 
    integer :: i

    ! can't use complex math (i.e., exp(-i theta)) because "x" and "y"
    ! components are themselves complex here
    rot(1,1) = cos(theta)
    rot(2,1) = sin(theta)
    rot(2,2) = rot(1,1)
    rot(1,2) = -rot(2,1)
    
    forall (i = 1:size(v,dim=1))
       w(i,1:2) = matmul(rot(1:2,1:2),v(i,1:2))
    end forall
    
  end function rotate_vel

  subroutine rotate_vel_mat(u,v,theta)
    use constants, only : DP
    complex(DP), dimension(:,:), intent(inout) :: u,v
    real(DP), intent(in) :: theta
    complex(DP), dimension(size(u,dim=2),2) :: yin, yout
    integer :: i

    do i = 1,size(u,dim=1)
       yin(:,1) = u(i,:)
       yin(:,2) = v(i,:)
       yout(:,1:2) = rotate_vel(yin,theta)
       u(i,:) = yout(:,1)
       v(i,:) = yout(:,2)
    end do
    
  end subroutine rotate_vel_mat

end module utility
