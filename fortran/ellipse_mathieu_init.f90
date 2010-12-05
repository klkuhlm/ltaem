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
    ! asymptotic expansion

    do i = 1, tnp
       if (bg%cutoff < 0.0) then
          bg%mat(i) = mathieu_init(kap(i), MM=max(bg%ms,dim(i)))
       else  ! only include if read from file (default = -999.)
          bg%mat(i) = mathieu_init(kap(i), MM=max(bg%ms,dim(i)), CUTOFF=bg%cutoff)
       end if
    end do

    ! allocate/initialize each element for each value of p
    do j = 1, size(e)
       if (e(j)%ibnd == 0 .or. e(j)%calcin) then
          allocate(e(j)%mat(tnp))

          kap = kappa(p,e(j)%element)
          dim(:) = shirts(e(j)%N, kap)

          do i = 1, tnp
             ! TODO adjustable cutoff not handled except in background 
             ! either add here or remove above.
             e(j)%mat(i) = mathieu_init(kap(i), MM=max(e(j)%MS,dim(i)))
          end do
       end if
    end do
  end subroutine ellipse_init
  
  elemental function shirts(n,q) result(dim)
    use constants, only : DP
    integer, intent(in) :: n
    complex(DP), intent(in) :: q
    integer :: dim
    real(DP) :: C,D
    
    ! estimate required matrix size, based on rational approximation due to 
    ! Shirts, R.B., 1993. "The Computation of Eigenvalues and Solutions of Mathieu's
    ! Differential Equation for Noninteger Order", ACM TOMS 19(3) pp377-390.

    C = (8.46 +   0.444*n)/(1.0_DP + 0.085*n)
    D = (0.240 + 0.0214*n)/(1.0_DP + 0.059*n)
    dim = int(n + 3.0_DP + C*abs(q)**D)

  end function shirts
end module ellipse_mathieu_init
