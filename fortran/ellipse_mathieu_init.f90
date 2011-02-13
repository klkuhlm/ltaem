!
! Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module ellipse_mathieu_init
  implicit none

  private
  public :: ellipse_init

contains
  subroutine ellipse_init(e,bg,s)
    use constants, only : DP
    use mathieu_functions, only : mathieu_init
    use type_definitions, only : ellipse, element
    use kappa_mod, only : kappa 

    type(ellipse), dimension(:), intent(inout) :: e
    type(element), intent(inout) :: bg
    complex(DP), dimension(:,:), intent(in) :: s
    complex(DP), dimension(product(shape(s))) :: p
    
    integer :: tnp, i, j
    integer, dimension(size(p,1)) :: dim
    complex(DP), dimension(size(p,1)) :: kap

    tnp = size(s)
    p = reshape(s,[tnp])

    ! allocate background
    allocate(bg%mat(tnp))
    kap = kappa(p,bg)

    ! initialize background for each value of p
    dim(:) = shirts(maxval(e(:)%N), kap)

    ! TODO: make a decision about matrix size here and
    ! either bail out if too large, or potentially use an 
    ! asymptotic expansion (needs investigation)

    write(*,'(A)',advance='no') 'bg: q-MS '
    do i = 1, tnp
       bg%mat(i) = mathieu_init(kap(i), MM=max(bg%ms,dim(i)), CUTOFF=bg%cutoff)
       write(*,'(3(I0,1X))',advance='no') i,bg%MS,dim(i)
    end do
    write(*,'(/)')

    ! allocate/initialize each element for each value of p
    do j = 1, size(e)
       if (e(j)%ibnd == 0 .or. e(j)%calcin) then
          allocate(e(j)%mat(tnp))

          kap = kappa(p,e(j)%element)
          dim(:) = shirts(e(j)%N, kap)

          write(*,'(A,I0,1X)',advance='no') 'el',j,': q-MS '
          do i = 1, tnp
             e(j)%mat(i) = mathieu_init(kap(i), MM=max(e(j)%MS,dim(i)), CUTOFF=e(j)%cutoff)
             write(*,'(3(I0,1X))',advance='no') i,e(j)%MS,dim(i)
          end do
          write(*,'(/)')
       end if
    end do
  end subroutine ellipse_init
  
  elemental function shirts(n,q) result(dim)
    use constants, only : DP
    integer, intent(in) :: n
    complex(DP), intent(in) :: q
    integer :: dim
    real(DP) :: C,D
    
    ! estimate required matrix size to achieve accuracy ~ 1.0E-12,
    ! based on rational approximation due to  Shirts, R.B., 1993.
    ! "The Computation of Eigenvalues and Solutions of Mathieu's
    ! Differential Equation for Noninteger Order", ACM TOMS 19(3) pp377-390.

    C = (8.46_DP +   0.444_DP*n)/(1.0_DP + 0.085_DP*n)
    D = (0.240_DP + 0.0214_DP*n)/(1.0_DP + 0.059_DP*n)
    dim = int(n + 3.0_DP + C*abs(q)**D)

  end function shirts
end module ellipse_mathieu_init
