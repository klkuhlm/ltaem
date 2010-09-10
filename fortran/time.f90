module time_mod
  implicit none

  private 
  public :: time

  interface time
     module procedure Time_pScal, Time_pVect
  end interface

contains

  ! ##################################################
  ! general time behavior function, for all elements
   function Time_pvect(p,t,area) result(mult)
    use constants, only: DP
    use type_definitions, only : time
    use utility, only : outer
    implicit none

    type(time), intent(in) :: t
    complex(DP), dimension(:), intent(in) :: p
    logical, intent(in) :: area
    complex(DP), dimension(size(p,1)) :: mult

    integer :: n, np, flag
    real(DP), allocatable :: ti(:), q(:), par(:)
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
       mult(1:np) = exp(-par(1)*p)/p - exp(-par(2)*p)/p
    case(3)
       ! instantaneous at t=par1
       mult(1:np) = exp(-par(1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult(1:np) = 1.0/(p - p*exp(-par(1)*p)) * &
                      & (1.0 - exp(-par(2)*p))/p
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
            &  (1.0 - exp(-par(1)*p/2.0))/ &
            & ((1.0 + exp(-par(1)*p/2.0))*p)
    case(:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = -flag
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       Q(0) = 0.0
       Q(1:n) = par(n+2:2*n+1)
       
       mult(1:np) = (sum(spread(Q(1:n) - Q(0:n-1),2,np)*&
            & exp(-outer(ti(1:n),p(1:np))),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p(:)))/p(:)

       deallocate(ti,Q)

    case default
       write(*,'(A,I0)') 'Time_pvect: error in case for time behavior ',flag
       stop
    end select

    deallocate(par)

  end function Time_pvect

  ! wrapper to allow calling time routine with scalar p
  function Time_pscal(p,t,area) result(mult)
    use constants, only : DP
    use type_definitions, only : time

    type(time), intent(in) :: t
    complex(DP), intent(in) :: p
    logical, intent(in) :: area
    complex(DP) :: mult

    mult = sum(Time_pVect([p],t,area))

  end function Time_pscal
end module time_mod
