module kappa_mod
  implicit none

  private
  public :: kappa

  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains

  function kappa_pVect(p,el) result(q)
    use constants, only : DP, PI
    use type_definitions, only : element

    complex(DP), intent(in), dimension(:) :: p
    type(element), intent(in) :: el
    complex(DP), dimension(size(p)) :: q

    integer :: np
    complex(DP), dimension(size(p)) ::  kap2, exp2z

    np = size(p)
    if(el%leakFlag /= 0) then
       kap2(1:np) = sqrt(p(:)*el%aquitardSs/el%aquitardK)
       exp2z(1:np) = exp(-2.0*kap2(:)*el%aquitardb)
    end if
    
    !! leaky-ness
    !! ##############################
    select case(el%leakFlag)
    case(0)
       !! no leaky layer, standard definition
       q(:) = p(1:np)/el%alpha
    case(1)
       !! case I, no-drawdown condition at top of aquitard
       q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
            & (1.0 + exp2z(:))/(1.0 - exp2z(:))
    case(2)
       !! case II, no-flow condition at top of aquitard
       q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
            & (1.0 - exp2z(:))/(1.0 + exp2z(:))
    case(3)
       !! aquitard thickness -> infinity
       q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)
    case default
       stop 'ERROR: incorrect value for leakFlag parameter -> (1,2,3)'
    end select

    !! unconfined-ness (if confined do nothing)
    !! ##############################
    if(el%unconfinedFlag) then

       ! results of integrating Neuman 1972 unconfined solution from 0 -> b in z
       ! Kz parameter is integrated away?
       q(1:np) = q(:) + el%Sy*p(:)/(el%b*el%K)
       
    end if
    
    !! sources are only additive under the square root
    q = sqrt(q)

  end function kappa_pVect
  
  !! scalar version useful in matching
  function kappa_pscal(p,el) result(q)
    use constants, only : DP
    use type_definitions, only : element
    complex(DP), intent(in) :: p
    type(element), intent(in) :: el
    complex(DP) :: q

    !! sum away singleton first dimension
    q = sum(kappa_pVect([p],el))

  end function kappa_pscal
end module kappa_mod
