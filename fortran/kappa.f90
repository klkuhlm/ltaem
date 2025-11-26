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

module kappa_mod
  implicit none

  private
  public :: kappa

  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains

  function kappa_pVect(p,el,sq_arg) result(q)
    use constants, only : DP, PISQ
    use type_definitions, only : element
    use utility, only : tanh

    complex(DP), intent(in), dimension(:) :: p
    type(element), intent(in) :: el ! circle, ellipse, or background
    complex(DP), dimension(size(p)) :: q
    logical, intent(in), optional :: sq_arg
    logical :: squared

    integer :: np, i
    complex(DP), allocatable ::  kap2(:), exp2z(:)

    real(DP) :: omega, beta
    real(DP), allocatable :: a(:), pdf(:)

    if (.not. present(sq_arg)) then
      squared = .false.
    else
      squared = sq_arg
    end if

    np = size(p)

    !! standard confined definition
    q(1:np) = p(:)/el%alpha

    if (el%wave) then
      !! wave solution (another time derivative)
      q(1:np) = q(:) * p
    end if

    !! leaky
    !! ##############################
    if (el%leakFlag /= 0) then
       allocate(kap2(np),exp2z(np))
       kap2(1:np) = sqrt(p(:)*el%aquitardSs/el%aquitardK)
       exp2z(1:np) = exp(-2.0_DP*kap2(:)*el%aquitardb)

       select case(el%leakFlag)
       case(1)
          !! case I, no-drawdown condition at top of aquitard
          q(1:np) = q + kap2(:)*el%aquitardK/(el%b*el%K)*&
               & (1.0_DP + exp2z(:))/(1.0_DP - exp2z(:))
       case(2)
          !! case II, no-flow condition at top of aquitard
          q(1:np) = q + kap2(:)*el%aquitardK/(el%b*el%K)*&
               & (1.0_DP - exp2z(:))/(1.0_DP + exp2z(:))
       case(3)
          !! aquitard thickness -> infinity
          q(1:np) = q + kap2(:)*el%aquitardK/(el%b*el%K)
       end select

       deallocate(kap2,exp2z)
    end if

    !! unconfined
    !! ##############################
    if(el%unconfinedFlag) then
       ! integrated (0->b dz) Neuman (1972, WRR) solution
       q(1:np) = q + el%Kz*p(:)/(el%K*el%b*(el%Kz/el%Sy + p(:)))
    end if

    !! dual porosity
    !! ##############################
    if(el%dualPorosityFlag) then

       if (el%multiporosityDiffusion == 0) then
          ! solution adapted from Dougherty and Babu (1984, WRR)
          q(1:np) = q + &
               & (1.0_DP - el%lambda/(el%matrixSs*p(:) + el%lambda))*el%lambda/el%K
       else
          ! multiporosity diffusion
          omega = el%Ss/(el%Ss + el%matrixSs) ! fraction of total storage in fractures
          beta = (1.0_DP - omega)/omega ! "capacity ratio"
          allocate(a(el%NDiffterms),pdf(el%NDiffterms))

          if (el%NDiffTerms > 0) then
             select case(el%multiporosityDiffusion)
             case(1)
                do concurrent (i=1:el%NDiffterms)
                  ! Warren-Root lambda = kappa/((1-omega)*LD**2)
                  a(i) = real((2*i-1)**2,DP) *PISQ*el%lambda*0.25_DP
                  pdf(i) = 8.0_DP*beta/(real((2*i-1)**2,DP) *PISQ)
                end do
             case(2)
                stop 'cylindrical multiporosity matrix diffusion not impelemented yet'
             case(3)
                do concurrent (i=1:el%NDiffterms)
                  a(i) = real(i**2,DP) *PISQ*el%lambda
                  pdf(i) = 6.0_DP*beta/(real(i**2,DP) *PISQ)
                end do
             case default
                stop 'invalid multiporosity matrix diffusion index'
             end select

             ! TODO: this is not an additive term (FIX?)
             q = p*omega*(1.0_DP + sum(spread(el%lambda*pdf*a,2,np)/ &
                  & (spread(p,1,el%NDiffterms) + spread(a,2,np)),dim=1))

          elseif (el%NDiffTerms == 0) then
             ! analytical solutions for diffusion into matrix
             ! activated for case where Nterms == 0
             select case(el%multiporosityDiffusion)
             case(1)
                q = p*(omega + sqrt(el%lambda*(1.0_DP-omega)/(3.0_DP*p))*&
                     & tanh(sqrt(3.0_DP*(1.0_DP-omega)*p/el%lambda)))
             case(2)
                stop 'cylindrical multiporosity matrix diffusion not impelemented yet'
             case(3)
                q = p*(omega + (sqrt(15.0_DP*(1.0_DP-omega*p)/el%lambda)/ &
                     & tanh(sqrt(15.0_DP*(1.0_DP-omega*p)/el%lambda)) - 1.0_DP)/(5.0_DP*p))
             case default
                stop 'invalid multiporosity matrix diffusion index'
             end select
          end if
          deallocate(a,pdf)
       end if
    end if

    !! sources additive under square root
    if (.not. squared) then
      q = sqrt(q)
    end if

  end function kappa_pVect

  !! scalar version useful in matching
  function kappa_pscal(p,el,sq_arg) result(q)
    use constants, only : DP
    use type_definitions, only : element
    complex(DP), intent(in) :: p
    type(element), intent(in) :: el
    logical, optional, intent(in) :: sq_arg
    complex(DP) :: q
    logical :: squared

    if (.not. present(sq_arg)) then
      squared = .false.
    else
      squared = sq_arg
    end if

    q = sum(kappa_pVect([p],el,squared))

  end function kappa_pscal

end module kappa_mod
