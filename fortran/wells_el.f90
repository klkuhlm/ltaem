! $Id: wells_el.f90,v 1.1 2007/07/30 05:58:04 kris Exp kris $
module wells_el
  implicit none
  private
  public  :: WellHead, WellFlux, ElTimeArea, ElTimeBdry !this needs a different home eventually

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface WellHead
     module procedure wellhead_rvect
     module procedure wellhead_pvect
  end interface

  interface WellFlux
     module procedure wellflux_rvect
     module procedure wellflux_pvect
  end interface

  interface WellTime
     module procedure welltime_pscal
     module procedure welltime_pvect
  end interface
  
  interface ElTimeArea
     module procedure eltimearea_pscal
     module procedure eltimearea_pvect
  end interface

  interface ElTimeBdry
     module procedure eltimebdry_pscal
     module procedure eltimebdry_pvect
  end interface


contains

  !##################################################
  ! head effects of instantaneous point source
  ! radius is a vector, p and q are scalars (matching)

   function WellHead_rvect(id,p,q,r) result(Lhead)
    use constants, only: DP, TWOPI, CZERO
    use element_specs, only : Wlr, Wlq
    use complex_bessel, only : cbesk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), intent(in) :: p, q
    real(DP), dimension(:), intent(in) :: r
    complex(DP), dimension(size(r,1)) :: Lhead
    complex(DP), dimension(size(r,1)) :: bessk0
    complex(DP),dimension(1) :: bessk1
    integer :: nR, ierr,nz,i

    nR = size(r,1)

    ! calculate head response due to instataneous point source
    call cbesk(q*WLr(id),1.0_DP,1,1,bessk1(1),ierr,nz)
    do i=1,nR
       call cbesk(q*r(i),0.0_DP,1,1,bessk0(i),ierr,nz)
    end do
    
    Lhead(1:nR) = WLq(id)*bessk0(1:nR)*WellTime(id,p)/(TWOPI*q*Wlr(id)*bessk1(1))
  end function WellHead_rvect

  !##################################################
  ! radius is a scalar, p and q are vectors (calculation)

   function WellHead_pvect(id,p,q,r) result(Lhead)
    use constants, only: DP, TWOPI, CZERO
    use element_specs, only: Wlq, Wlr
    use complex_bessel, only : cbesk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p, q
    real(DP), intent(in) :: r
    complex(DP), dimension(size(p,1)) :: Lhead
    complex(DP), dimension(size(p,1)) :: bessk0, bessk1
    integer ::  nP,ierr,nz,i

    nP = size(p,1)
    
    do i=1,nP
       call cbesk(q(i)*WLr(id),1.0_DP,1,1,bessk1(i),ierr,nz)
       call cbesk(q(i)*r,      0.0_DP,1,1,bessk0(i),ierr,nz)
    end do

    Lhead(1:nP) = WLq(id)*bessk0(1:nP)*WellTime(id,p)/(TWOPI*q*Wlr(id)*bessk1(1:nP))
  end function WellHead_pvect

  !##################################################
  ! flux effects of instantaneous point source
  ! radius is a vector, p and q are scalars (matching)

   function WellFlux_rvect(id,p,q,r) result(Lflux)
    use constants, only: DP, TWOPI, CZERO
    use element_specs, only : Wlq, Wlr
    use complex_bessel, only : cbesk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), intent(in) :: p, q
    real(DP), dimension(:), intent(in) :: r
    complex(DP), dimension(size(r,1)) :: Lflux
    complex(DP), dimension(size(r,1)) :: bessk1
    complex(DP), dimension(1) :: bessk1r0
    integer :: nR,ierr,nz,i

    nR = size(r,1)

    call cbesk(q*WLr(id),1.0_DP,1,1,bessk1r0(1),ierr,nz)
    do i=1,nR
       call cbesk(q*r(i),1.0_DP,1,1,bessk1(i),ierr,nz)
    end do

    Lflux(1:nR) = -WLq(id)*bessk1(1:nR)*WellTime(id,p)/&
         & (TWOPI*Wlr(id)*bessk1r0(1))
  end function WellFlux_rvect

  !##################################################
  ! radius is a scalar, p and q are vectors (calculation)

   function WellFlux_pvect(id,p,q,r) result(Lflux)
    use constants, only: DP, TWOPI, RONE, CZERO
    use element_specs, only : Wlr, Wlq
    use complex_bessel, only: cbesk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p, q
    real(DP), intent(in) :: r
    complex(DP), dimension(size(p,1)) :: Lflux
    ! two values of K1 for rw and rwm (one-based vector)
    complex(DP), dimension(size(p,1)) :: bessk1, bessk1r0
    integer :: nP,i,ierr,nz

    nP = size(p,1)
    
    do i=1,nP
       call cbesk(q(i)*r,      1.0_DP,1,1,bessk1(i),ierr,nz)
       call cbesk(q(i)*WLr(id),1.0_DP,1,1,bessk1r0(i),ierr,nz)
    end do

    Lflux(1:nP) = -WLq(id)*bessk1(1:nP)/(TWOPI*Wlr(id)*bessk1r0(1:nP))*WellTime(id,p)
  end function WellFlux_pvect

  !##################################################
  ! this function provides the time-dependent portion
  ! of the well response. Solution comes from convolving 
  ! this with the instantan. point source (* in Lap domain)

   function WellTime_pscal(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : WLtime, WLtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), intent(in) :: p
    complex(DP) :: mult
    
    select case(WLtime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-WLtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-WLtpar(id,1)*p)/p * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = RONE
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-WLtpar(id,1)*p)) * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-WLtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + WLtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(WLtpar(id,1)*p) - exp(WLtpar(id,1)*p))/ &
            &(exp(WLtpar(id,1)*p) + exp(WLtpar(id,1)*p))
    end select
  end function WellTime_pscal

  ! ##################################################
  ! version for p a vector (code is the same)
   function WellTime_pvect(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only: WLtime, WLtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    
    select case(WLtime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-WLtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-WLtpar(id,1)*p)/p * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = cmplx(RONE,0.0_DP)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-WLtpar(id,1)*p)) * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-WLtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + WLtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(WLtpar(id,1)*p) - exp(WLtpar(id,1)*p))/ &
            &(exp(WLtpar(id,1)*p) + exp(WLtpar(id,1)*p))
    end select
  end function WellTime_pvect

  !##################################################
  ! this function provides the time-dependent portion
  ! of a circular inclusion's AREA FLUX response. 
  ! solution comes from convolving 
  ! this with the instantan. solution (* in Lap domain)

   function ElTimeArea_pscal(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : EIAreaTime, EIAtpar
    implicit none

    integer, intent(in) :: id
    complex(DP),  intent(in) :: p
    complex(DP) :: mult
    
    select case(EIAreaTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-EIAtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-EIAtpar(id,1)*p)/p * &
            & (RONE - exp(-EIAtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = cmplx(RONE,0.0_DP)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-EIAtpar(id,1)*p)) * &
            & (RONE - exp(-EIAtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-EIAtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + EIAtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(EIAtpar(id,1)*p) - exp(EIAtpar(id,1)*p))/ &
            &(exp(EIAtpar(id,1)*p) + exp(EIAtpar(id,1)*p))
    end select
  end function ElTimeArea_pscal

  function ElTimeArea_pvect(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : EIAreaTime, EIAtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    
    select case(EIAreaTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-EIAtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-EIAtpar(id,1)*p)/p * &
            & (RONE - exp(-EIAtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = cmplx(RONE,0.0_DP)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-EIAtpar(id,1)*p)) * &
            & (RONE - exp(-EIAtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-EIAtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + EIAtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(EIAtpar(id,1)*p) - exp(EIAtpar(id,1)*p))/ &
            &(exp(EIAtpar(id,1)*p) + exp(EIAtpar(id,1)*p))
    end select
  end function ElTimeArea_pvect

  !##################################################
  ! this function provides the time-dependent portion
  ! of a circular inclusion's BOUNDARY HEAD/FLUX response. 
  ! solution comes from convolving 
  ! this with the instantan. solution (* in Lap domain)

   function ElTimeBdry_pscal(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : EIBdryTime, EIBtpar
    implicit none

    integer, intent(in) :: id
    complex(DP),  intent(in) :: p
    complex(DP) :: mult
    
    select case(EIBdryTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-EIBtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-EIBtpar(id,1)*p)/p * &
            & (RONE - exp(-EIBtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = cmplx(RONE,0.0_DP)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-EIBtpar(id,1)*p)) * &
            & (RONE - exp(-EIBtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-EIBtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + EIBtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(EIBtpar(id,1)*p) - exp(EIBtpar(id,1)*p))/ &
            &(exp(EIBtpar(id,1)*p) + exp(EIBtpar(id,1)*p))
    end select
  end function ElTimeBdry_pscal

  function ElTimeBdry_pvect(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : EIBdryTime,EIBtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    
    select case(EIBdryTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-EIBtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-EIBtpar(id,1)*p)/p * &
            & (RONE - exp(-EIBtpar(id,2)*p))/p
    case(3)
       ! instantaneous
       mult = cmplx(RONE,0.0_DP)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE / (p - p*exp(-EIBtpar(id,1)*p)) * &
            & (RONE - exp(-EIBtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = RONE /(p + p*exp(-EIBtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = RONE /(p**2 + EIBtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = RONE / p**2 * &
            &(exp(EIBtpar(id,1)*p) - exp(EIBtpar(id,1)*p))/ &
            &(exp(EIBtpar(id,1)*p) + exp(EIBtpar(id,1)*p))
    end select
  end function ElTimeBdry_pvect

end module wells_el
