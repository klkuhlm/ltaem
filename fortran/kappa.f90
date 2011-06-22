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

module kappa_mod
  implicit none

  private
  public :: kappa

  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains

  function kappa_pVect(p,el) result(q)
    use constants, only : DP
    use type_definitions, only : element

    complex(DP), intent(in), dimension(:) :: p
    type(element), intent(in) :: el ! circle, ellipse, or background
    complex(DP), dimension(size(p)) :: q

    integer :: np
    complex(DP), allocatable ::  kap2(:), exp2z(:)

    np = size(p)
    if (el%leakFlag /= 0) then
       allocate(kap2(np),exp2z(np))
       kap2(1:np) = sqrt(p(:)*el%aquitardSs/el%aquitardK)
       exp2z(1:np) = exp(-2.0*kap2(:)*el%aquitardb)
    end if
    
    !! leaky-ness
    !! ##############################
    select case(el%leakFlag)
    case(0)
       !! no leaky layer, standard definition
       q(1:np) = p(:)/el%alpha
    case(1)
       !! case I, no-drawdown condition at top of aquitard
       q(1:np) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
            & (1.0 + exp2z(:))/(1.0 - exp2z(:))
    case(2)
       !! case II, no-flow condition at top of aquitard
       q(1:np) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
            & (1.0 - exp2z(:))/(1.0 + exp2z(:))
    case(3)
       !! aquitard thickness -> infinity
       q(1:np) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)
    end select

    !! unconfined-ness (if confined do nothing)
    !! ##############################
    if(el%unconfinedFlag) then

       ! results of integrating Neuman 1972 unconfined solution from 0 -> b in z
!!$       q(1:np) = q(:) + el%Sy*p(:)/(el%b*el%K)
       q(1:np) = q(:) + el%Kz*p(:)/(el%K*el%b*(el%Kz/el%Sy + p(:)))
    end if
    
    !! sources are only additive under the square root
    q = sqrt(q)

    if (el%leakFlag /= 0) deallocate(kap2,exp2z)
   
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
