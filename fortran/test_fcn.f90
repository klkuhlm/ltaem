function imagefcn(p) result(img)

  use constants, only : DP, PI,CONE
  use invlap_shared_data
  implicit none

  complex(DP), dimension(:), intent(in) :: p
  complex(DP), dimension(size(p,1)) :: img
  
  integer :: np
  
  np = size(p,1)

  select case(index)
  case(1)
     img(1:np) = CONE/(p+sqrt(p)*a)
  case(2)
     img(1:np) = exp(-a/p)/sqrt(p)
  case(3)
     img(1:np) = p/(p**2 + a**2)**2
  case(4)
     img(1:np) = exp(-a*p)/p
  case(5)
     img(1:np) = exp(-a/p)/p
  case(6)
     img(1:np) = CONE/(p+a) - CONE/(p+10._DP*a)
  case(7)
     img(1:np) = CONE/(p+p*exp(-a/2.0_DP*p))
  end select
        
end function imagefcn

  
