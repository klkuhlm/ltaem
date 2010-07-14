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
    complex(DP), dimension(size(p,1)) :: kap

    tnp = product(shape(s))
    p = reshape(s,[tnp])

    ! allocate background
    allocate(bg%mat(tnp))
    kap = kappa(p,bg)

    ! initialize for each value of p
    do i = 1, tnp
       bg%mat(i) = mathieu_init(kap(i),MM=bg%ms)
    end do

    ! allocate/initialize each element for each value of p
    do j = 1, size(e,1)
       if (e(j)%ibnd == 0 .or. e(j)%calcin) then
          allocate(e(j)%mat(tnp))
          kap = kappa(p,e(j)%element)

          do i = 1, tnp
             e(j)%mat(i) = mathieu_init(kap(i),MM=e(j)%MS)
          end do
       end if
    end do
  end subroutine ellipse_init
end module ellipse_mathieu_init
