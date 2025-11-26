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

module time_mod
  implicit none

  private
  public :: timef

  interface timef
     module procedure Time_pScal, Time_pVect
  end interface

contains

  ! ##################################################
  ! general time behavior function, for all elements
   pure function Time_pvect(p,t,area) result(mult)
    use constants, only: DP
    use type_definitions, only : time
    use utility, only : outer
    implicit none

    type(time), intent(in) :: t
    complex(DP), dimension(:), intent(in) :: p
    logical, intent(in) :: area
    complex(DP), dimension(size(p,1)) :: mult

    integer :: n, np, flag
    real(DP), allocatable :: ti(:), y(:), dy(:), W(:), par(:)
    real(DP), allocatable :: denom(:), numer(:)
    real(DP) :: tf

    np = size(p,1)

    ! consolidate area and boundary time functions
    if (area) then
       flag = t%areaTime
       allocate(par(size(t%ATPar,dim=1)))
       par(:) = t%AtPar(:)
    else
       flag = t%BdryTime
       allocate(par(size(t%BTPar,dim=1)))
       par(:) = t%BTPar(:)
    end if

    select case(flag)
    case(1)
       ! step on at time=par1
       mult(1:np) = exp(-par(1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult(1:np) = (exp(-par(1)*p) - exp(-par(2)*p))/p
    case(3)
       ! instantaneous at t=par1
       mult(1:np) = exp(-par(1)*p)
    case(4)
       ! "step test": increasing by integer multiples of Q each
       ! integer multiple of par1 time, off at par2
       mult(1:np) = 1.0_DP/(p - p*exp(-par(1)*p)) * &
                      & (1.0_DP - exp(-par(2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-par(2)*p)/(p + p*exp(-par(1)*p))
    case(6)
       ! cosine wave  - cos(at)
       ! shifted to start at t=par2
       mult(1:np) = exp(-par(2)*p)*p/(p**2 + par(1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-par(2)*p) / p**2 * &
            &(exp(par(1)*p) - exp(par(1)*p))/ &
            &(exp(par(1)*p) + exp(par(1)*p))
    case(8)
       ! full square wave (only +, then -), period 2*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-par(2)*p)* &
            &  (1.0_DP - exp(-par(1)*p*0.5_DP))/ &
            & ((1.0_DP + exp(-par(1)*p*0.5_DP))*p)
    case(9)
       ! sine wave shifted to start at t=par2
       mult(1:np) = exp(-par(2)*p)*par(1)/(p**2 + par(1)**2)
    case(-100:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = -flag
       allocate(ti(n),y(0:n),dy(n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       y(0) = 0.0_DP
       y(1:n) = par(n+2:2*n+1)
       dy = y(1:n) - y(0:n-1)

       mult(1:np) = (sum(spread(dy(:),2,np)*exp(-outer(ti(1:n),p(:))),1) - &
            & sum(dy(:))*exp(-tf*p(:)))/p(1:np)

       deallocate(ti,y,dy)

    case(:-101)
       !! piecewise linear pumping rate with n steps, from ti(1) to tf
       !! no jumps in value (no vertical slopes)
       n = -flag - 100
       allocate(ti(n),W(0:n+1),y(1:n+1),denom(1:n),numer(1:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       y(1:n) = par(n+2:2*n+1)

       ! compute slope between each pair of points
       W(0) = 0.0_DP
       W(n+1) = 0.0_DP
       y(n+1) = 0.0_DP
       denom = [ti(2:n),tf] - ti(1:n)
       numer = y(2:n+1) - y(1:n)
       where (abs(denom) < epsilon(abs(numer)))
          ! TODO: probably a better way to handle this
          denom = denom + spacing(denom)
       end where
       W(1:n) = numer/denom ! rise/run

       mult(1:np) = (sum(spread(W(1:n) - W(0:n-1),2,np)*&
            & exp(-outer(ti(1:n),p(1:np))),dim=1) - &
            & sum(W(1:n) - W(0:n-1))*exp(-tf*p(:)))/p(:)**2

       deallocate(ti,W,y,denom)

    end select
    deallocate(par)

  end function Time_pvect

  ! wrapper to allow calling time routine with scalar p
  pure function Time_pscal(p,t,area) result(mult)
    use constants, only : DP
    use type_definitions, only : time

    type(time), intent(in) :: t
    complex(DP), intent(in) :: p
    logical, intent(in) :: area
    complex(DP) :: mult

    mult = sum(Time_pVect([p],t,area))

  end function Time_pscal
end module time_mod

