!
! Copyright (c) 2011-2025 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

! some general? functions collected in one place

module utility
  implicit none

  private
  public :: v2c, diag, logspace, linspace, outer, ynot, rotate_vel, rotate_vel_mat
  public :: cos_recurrence, sin_recurrence, cosh, acosh, tanh

  interface diag
     module procedure diagonal_z, diagonal_d
  end interface

  interface logspace
     module procedure logspace_int, logspace_r, logspace_d
  end interface

  interface outer
     module procedure outerprod_r, outerprod_d, outerprod_z, outerprod_dz, outerprod_zd
  end interface

  ! overload intrinsic hyperbolic functions for complex argument
  interface cosh
     module procedure ccosh
  end interface cosh

  interface acosh
     module procedure cacosh
  end interface acosh

  interface tanh
     module procedure ctanh
  end interface tanh

contains

  function cos_recurrence(x,n) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: x
    integer, intent(in) :: n
    real(DP), dimension(size(x),0:n-1) :: c
    integer :: i

    c(:,0) = 1.0_DP
    c(:,1) = cos(x)
    do i = 2,N-2
       c(:,i) = 2.0_DP*c(:,1)*c(:,i-1) - c(:,i-2)
    end do
  end function cos_recurrence

  function sin_recurrence(x,n) result(s)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: x
    integer, intent(in) :: n
    real(DP), dimension(size(x),1:n-1) :: s
    real(DP), dimension(size(x)) :: c
    integer :: i

    s(:,1) = sin(x)
    s(:,2) = sin(2.0_DP*x)
    c = cos(x)
    do i = 3,N-2
       s(:,i) = 2.0_DP*c*s(:,i-1) - s(:,i-2)
    end do
  end function sin_recurrence

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  pure function v2c(v) result(z)
    use constants, only : DP
    real(DP), intent(in), dimension(2) :: v
    complex(DP) :: z
    z = cmplx(v(1),v(2),DP)
  end function v2c

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
    complex(DP), dimension(size(db,1)) :: cb
    complex(DP), dimension(size(ca),size(db)) :: c
    cb = cmplx(db,kind=DP)
    c = spread(ca,dim=2,ncopies=size(db))*spread(cb,dim=1,ncopies=size(ca))
  end function outerprod_dz

  pure function outerprod_zd(da,cb) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da
    complex(DP), dimension(size(da,1)) :: ca
    complex(DP), intent(in), dimension(:) :: cb
    complex(DP), dimension(size(da),size(cb)) :: c
    ca = cmplx(da,kind=DP)
    c = spread(ca,dim=2,ncopies=size(cb))*spread(cb,dim=1,ncopies=size(da))
  end function outerprod_zd

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  elemental function ccosh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    x = real(z)
    y = aimag(z)
    f = cmplx(cosh(x)*cos(y), sinh(x)*sin(y),DP)
  end function ccosh

  elemental function ctanh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    x = real(z)
    y = aimag(z)
    f = cmplx(tanh(2.0_DP*x)/(1.0_DP+cos(2.0_DP*y)/cosh(2.0_DP*x)), &
         & sin(2.0_DP*y)/(cosh(2.0_DP*x)+cos(2.0_DP*y)),DP)
  end function ctanh

  elemental function cacosh(z) result(f)
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: f

    ! compared ~10^-15 against acosh() in Matlab
    ! branch cut is left of +1, along x axis
    ! -pi <= aimag(f) <= +pi
    if(real(z) >= 0.0_DP) then
       f = log(z + sqrt(z**2 - 1.0_DP))
    else
       f = -log(z + sqrt(z**2 - 1.0_DP))
    end if
  end function cacosh

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  pure function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: dx

    if (num == 1) then
       v = [lo]
    else
       dx = (hi-lo)/(num-1)
       do concurrent (i=1:num)
         v(i) = lo + (i-1)*dx
       end do
    end if
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
       do concurrent (i=1:size(A,dim=1)-n)
         d(i) = A(i+n,i)
       end do
     else
        ! diagonal in upper triangle (-)
       do concurrent (i=1:size(A,dim=1)+n)
         d(i) = A(i,i-n)
       end do
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
       do concurrent (i=1:size(A,dim=1)-n)
         d(i) = A(i+n,i)
       end do
     else
        ! diagonal in upper triangle (-)
       do concurrent (i=1:size(A,dim=1)+n)
         d(i) = A(i,i-n)
       end do
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
    interface
       ! Level-1 BLAS routine for real rotation of complex vector
       ! rotation is in opposite sense, so use -theta below
       subroutine ZDROT(N,CX,INCX,CY,INCY,C,S)
         integer, intent(in) :: n,incx,incy
         real(8), intent(in) :: c,s
         complex(8), dimension(n), intent(inout) :: cx,cy
       end subroutine ZDROT
    end interface

    complex(DP), dimension(:,:), intent(in) :: v
    real(DP), intent(in) :: theta
    complex(DP), dimension(size(v,dim=1),2) :: w

    w = v
    call zdrot(size(v,dim=1),w(:,1),1,w(:,2),1,cos(-theta),sin(-theta))

  end function rotate_vel

  subroutine rotate_vel_mat(u,v,theta)
    use constants, only : DP
    interface
       subroutine ZDROT(N,CX,INCX,CY,INCY,C,S)
         integer, intent(in) :: n,incx,incy
         real(8), intent(in) :: c,s
         complex(8), dimension(n), intent(inout) :: cx,cy
       end subroutine ZDROT
    end interface

    complex(DP), dimension(:,:), intent(inout) :: u,v
    complex(DP), dimension(product(shape(u))) :: uu,vv
    real(DP), intent(in) :: theta
    integer :: n,m,nel

    n = size(u,dim=1)
    m = size(u,dim=2)
    nel = n*m

    uu = reshape(u,[nel])
    vv = reshape(v,[nel])
    call zdrot(nel,uu,1,vv,1,cos(-theta),sin(-theta))
    u = reshape(uu,[n,m])
    v = reshape(vv,[n,m])

  end subroutine rotate_vel_mat
end module utility
