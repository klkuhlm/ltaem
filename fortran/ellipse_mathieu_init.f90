module ellipse_mathieu_init
  implicit none

contains
  subroutine ellipse_init(e,bg,p)
    use constants, only : DP
    use mathieu_functions, only : mathieu_init
    use type_definitions, only : ellipse, element
    use kappa_mod, only : kappa 

    type(ellipse), dimension(:), intent(inout) :: e
    type(element), intent(inout) :: bg
    complex(DP), dimension(:), intent(in) :: p

    integer :: np, i, j
    complex(DP), dimension(size(p,1)) :: kap

#ifdef DEBUG
    print *, 'ellipse_init: e:',e%id,' bg:',bg%id,' p:',p
#endif

    np = size(p,1)

    ! allocate background
    allocate(bg%mat(np))
    kap = kappa(p,bg)

    ! initialize for each value of p
    do i = 1, np
       bg%mat(i) = mathieu_init(kap(i))
    end do

    ! allocate/initialize each element for each value of p
    do j = 1, size(e,1)
       if (e(j)%ibnd == 0 .or. e(j)%calcin) then
          allocate(e(j)%mat(np))
          kap = kappa(p,e(j)%element)
          do i = 1, np
             e(j)%mat(i) = mathieu_init(kap(i))
          end do
       end if
    end do
  end subroutine ellipse_init
end module ellipse_mathieu_init
