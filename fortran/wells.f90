! $Id: wells.f90,v 1.7 2008/12/10 02:47:51 kris Exp kris $
module wells
  implicit none
  private
  public  :: WellHead, WellFlux, CircTimeArea, CircTimeBdry !this needs a different home eventually

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
  
  interface CircTimeArea
     module procedure circtimearea_pscal
     module procedure circtimearea_pvect
  end interface

  interface CircTimeBdry
     module procedure circtimebdry_pscal
     module procedure circtimebdry_pvect
  end interface


contains

  !##################################################
  ! does some error checking and calculates Bessel functions 
  ! needed for head effects of instantaneous point source
  
  ! radius is a vector, p and kappa are scalars (matching)

   function WellHead_rvect(id,p,kappa,r) result(Lhead)
    use constants, only: DP, TWOPI, CZERO, PI
    use element_specs, only : Wlr, Wlq, WLstor, WLdskin, kv, CIWellUp
    use bessel_functions, only : besselk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), intent(in) :: p, kappa
    real(DP), dimension(:), intent(in) :: r
    complex(DP), dimension(size(r,1)) :: Lhead, aw,bessk0
    complex(DP) :: bessk1
    integer :: nR
    real(DP) :: k, cw

    nR = size(r,1)
    Lhead = CZERO
    
    bessk1 = besselk(kappa*WLr(id),1)
    bessk0(1:nR) = besselk(kappa*r(:),0)

    if(.not.WLstor(id)) then
       aw(1:nR) = WLq(id)*WellTime(id,p)
    else
       k = kv(CIWellUp(id))
       ! assume casing radius = well radius
       ! assume wellbore storage = area
       cw = PI*WLr(id)**2  
       
!!$       aw(1:nR) = WLq(id)*WellTime(id,p)/(1.0_DP + p/(TWOPI*k)*&
!!$              & (WLDSkin(id)*cw +  besselk(kappa*WLr(id),0)/&
!!$              & (bessk1*kappa*WLr(id))))
       aw(1:nR) = WLq(id)*WellTime(id,p)/(1.0_DP + cw*p/(TWOPI*k)*&
              & (WLDSkin(id) +  besselk(kappa*WLr(id),0)/&
              & (bessk1*kappa*WLr(id))))
    end if

    ! potential due to finite radius well
    where (r > WLr(id))
       ! solution in the aquifer
       Lhead(1:nR) = aw(:)*bessk0(:)/(TWOPI*kappa*Wlr(id)*bessk1)
    elsewhere
       ! solution in the pumping well
       Lhead(1:nR) = aw(:)*(bessk0(:)/(kappa*bessk1) + &
            & WLDSkin(id))/(TWOPI*WLr(id))
    end where
    

  end function WellHead_rvect

  !##################################################
  ! radius is a scalar, p and kappa are vectors (calculation)

   function WellHead_pvect(id,p,kappa,r) result(Lhead)
    use constants, only: DP, TWOPI, PI, CZERO
    use element_specs, only: Wlq, Wlr, Wlstor, WLdskin, kv, CIWellUp
    use bessel_functions, only : besselk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p, kappa
    real(DP), intent(in) :: r
    complex(DP), dimension(size(p,1)) :: Lhead, aw, bessk0, bessk1
    integer ::  nP
    real(DP) :: k, cw

    nP = size(p,1)
    Lhead = CZERO
    
    bessk0(1:nP) = besselk(kappa(:)*r,0)
    bessk1(1:nP) = besselk(kappa(:)*WLr(id),1)

    if (.not. Wlstor(id)) then
       aw(1:nP) = WLq(id)*WellTime(id,p(:))
    else
       k = kv(ciwellup(id))
       cw = PI*WLr(id)**2       

       aw(1:nP) = WLq(id)*WellTime(id,p(:))/(1.0_DP + cw*p(:)/(TWOPI*k)*&
            &(WLDSkin(id) + besselk(kappa(:)*WLr(id),0)/&
            & (bessk1(:)*kappa(:)*WLr(id)))) 
    end if

    if (r > WLr(id)) then
       Lhead(1:nP) = aw(:)*bessk0(:)/(TWOPI*kappa(:)*Wlr(id)*bessk1(:))
    else
       Lhead(1:nP) = aw(:)*(bessk0(:)/(kappa(:)*bessk1(:)) + &
            & WLDSkin(id))/(TWOPI*WLr(id))
    end if
    

  end function WellHead_pvect

  !##################################################
  ! does some error checking and calculates Bessel functions 
  ! needed for flux effects of instantaneous point source
  
  ! radius is a vector, p and kappa are scalars (matching)

   function WellFlux_rvect(id,p,kappa,r) result(Lflux)
    use constants, only: DP, TWOPI, PI, CZERO
    use element_specs, only : Wlq, Wlr, Wlstor, WLdskin, kv, CIWellUp
    use bessel_functions, only : besselk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), intent(in) :: p, kappa
    real(DP), dimension(:), intent(in) :: r
    complex(DP), dimension(size(r,1)) :: Lflux, aw, bessk1, bessk0r0
    complex(DP) :: bessk1r0
    integer :: nR
    real(DP) :: k, cw

    nR = size(r,1)
    Lflux = CZERO

    bessk1r0 = besselk(kappa*WLr(id),1)
    bessk1(1:nR) = besselk(kappa*r(:),1)

    if(.not. WLstor(id)) then
       aw(1:nR) = WLq(id)*WellTime(id,p)

    else
       k = kv(ciwellup(id))
       cw = PI*WLr(id)**2 
       bessk0r0 = besselk(kappa*WLr(id),0)
       
       aw(1:nR) = WLq(id)*WellTime(id,p)/(1.0_DP + cw*p/(TWOPI*k)*&
            &(WLDSkin(id) + bessk0r0/(bessk1*kappa*WLr(id)))) 
    end if

    ! skin does not effect flux into well like it does head?
    Lflux(1:nR) = -aw(:)*bessk1(:)/(TWOPI*Wlr(id)*bessk1r0)   

  end function WellFlux_rvect

  !##################################################
  ! radius is a scalar, p and kappa are vectors (calculation)

   function WellFlux_pvect(id,p,kappa,r) result(Lflux)
    use constants, only: DP, TWOPI, RONE, PI, CZERO
    use element_specs, only : Wlr, Wlq, WLstor, WLdskin, kv, CIWellUp
    use bessel_functions, only : besselk
    implicit none
    
    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p, kappa
    real(DP), intent(in) :: r
    complex(DP), dimension(size(p,1)) :: Lflux, aw, bessk1, bessk1r0
    integer :: nP
    real(DP) :: k, cw

    nP = size(p,1)
    Lflux =  CZERO
    
    bessk1(1:nP) = besselk(kappa(1:nP)*r,1)
    bessk1r0(1:nP) = besselk(kappa(1:nP)*WLr(id),1)

    if(.not. WLstor(id)) then    
       aw(1:nP) = WLq(id)*WellTime(id,p(:))

    else
       k = kv(ciwellup(id))
       cw = PI*WLr(id)**2 

       aw(1:nP) = WLq(id)*WellTime(id,p(:))/(1.0_DP + cw*p(:)/(TWOPI*k)*&
            &(WLDSkin(id) + besselk(kappa(:)*WLr(id),0)/&
            &(bessk1(:)*kappa(:)*WLr(id))))            
    end if
    
    Lflux(1:nP) = -aw(:)*bessk1(:)/(TWOPI*Wlr(id)*bessk1r0(:))

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
    real(DP), allocatable :: ti(:),q(:)
    integer :: n
    real(DP) :: tf
    complex(DP) :: mult
    
    select case(WLtime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-WLtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-WLtpar(id,1)*p)/p - &
            & exp(-WLtpar(id,2)*p)/p
    case(3)
       ! instantaneous at time=par1
       mult = exp(-WLtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE/(p - p*exp(-WLtpar(id,1)*p)) * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult = exp(-WLtpar(id,2)*p)/(p + p*exp(-WLtpar(id,1)*p))
    case(6)
       ! cosine wave  - cos(at)
       ! shifted to start at t=par2
       mult = exp(-WLtpar(id,2)*p)*p/(p**2 + WLtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult = exp(-WLtpar(id,2)*p) / p**2 * &
            &(exp(WLtpar(id,1)*p) - exp(WLtpar(id,1)*p))/ &
            &(exp(WLtpar(id,1)*p) + exp(WLtpar(id,1)*p))
    case(8)
       ! full square wave (only +, then -), period 2*par1
       ! shifted to start at t=par2
       mult = exp(-WLtpar(id,2)*p)* &
            &  (RONE - exp(-WLtpar(id,1)*p/2.0))/ &
            & ((RONE + exp(-WLtpar(id,1)*p/2.0))*p)
    case(:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = abs(WLtime(id))
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = WLTpar(id,1:n); tf = WLTpar(id,n+1)
       Q(1:n) =  WLTpar(id,n+2:2*n+1); Q(0) = 0.0_DP
       
       mult = (sum((Q(1:n) - Q(0:n-1))*exp(-ti(1:n)*p)) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p))/p

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
    integer :: n, np
    real(DP), allocatable :: ti(:),q(:)
    real(DP) :: tf
    
    np = size(p,1)

    select case(WLtime(id))
    case(1)
       ! step on at time=par1
       mult(1:np) = exp(-WLtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult(1:np) = exp(-WLtpar(id,1)*p)/p - &
            & exp(-WLtpar(id,2)*p)/p
    case(3)
       ! instantaneous at t=par1
       mult(1:np) = exp(-WLtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult(1:np) = RONE/(p - p*exp(-WLtpar(id,1)*p)) * &
            & (RONE - exp(-WLtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-WLtpar(id,2)*p)/(p + p*exp(-WLtpar(id,1)*p))
    case(6)
       ! cosine wave  - cos(at)
       ! shifted to start at t=par2
       mult(1:np) = exp(-WLtpar(id,2)*p)*p/(p**2 + WLtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-WLtpar(id,2)*p) / p**2 * &
            &(exp(WLtpar(id,1)*p) - exp(WLtpar(id,1)*p))/ &
            &(exp(WLtpar(id,1)*p) + exp(WLtpar(id,1)*p))
    case(8)
       ! full square wave (only +, then -), period 2*par1
       ! shifted to start at t=par2
       mult(1:np) = exp(-WLtpar(id,2)*p)* &
            &  (RONE - exp(-WLtpar(id,1)*p/2.0))/ &
            & ((RONE + exp(-WLtpar(id,1)*p/2.0))*p)
    case(:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = abs(WLtime(id))
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = WLTpar(id,1:n); tf = WLTpar(id,n+1)
       Q(1:n) =  WLTpar(id,n+2:2*n+1); Q(0) = 0.0_DP
       
       mult(1:np) = (sum(spread(Q(1:n) - Q(0:n-1),2,np)*&
            & exp(-spread(ti(1:n),2,np)*spread(p(1:np),1,n)),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p(:)))/p(:)

    end select
  end function WellTime_pvect

  !##################################################
  ! this function provides the time-dependent portion
  ! of a circular inclusion's AREA FLUX response. 
  ! solution comes from convolving 
  ! this with the instantan. solution (* in Lap domain)

   function CircTimeArea_pscal(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : CIAreaTime, CIAtpar
    implicit none

    integer, intent(in) :: id
    complex(DP),  intent(in) :: p
    complex(DP) :: mult
    
    select case(CIAreaTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-CIAtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-CIAtpar(id,1)*p)/p - &
            & exp(-CIAtpar(id,2)*p)/p
    case(3)
       ! instantaneous at t=par1
       mult = exp(-CIAtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE/(p - p*exp(-CIAtpar(id,1)*p)) * &
            & (RONE - exp(-CIAtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult = exp(-CIAtpar(id,2)*p)/(p + p*exp(-CIAtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       ! shifted to start at t=par2
       mult = exp(-CIAtpar(id,2)*p)/(p**2 + CIAtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult = exp(-CIAtpar(id,2)*p)/p**2 * &
            &(exp(CIAtpar(id,1)*p) - exp(CIAtpar(id,1)*p))/ &
            &(exp(CIAtpar(id,1)*p) + exp(CIAtpar(id,1)*p))
    end select
  end function CircTimeArea_pscal

  function CircTimeArea_pvect(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : CIAreaTime, CIAtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    
    select case(CIAreaTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-CIAtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-CIAtpar(id,1)*p)/p - &
            & exp(-CIAtpar(id,2)*p)/p
    case(3)
       ! instantaneous at t=par1
       mult = exp(-CIAtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE/(p - p*exp(-CIAtpar(id,1)*p)) * &
            & (RONE - exp(-CIAtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = exp(-CIAtpar(id,2)*p)/(p + p*exp(-CIAtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = exp(-CIAtpar(id,2)*p)/(p**2 + CIAtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = exp(-CIAtpar(id,2)*p)/ p**2 * &
            &(exp(CIAtpar(id,1)*p) - exp(CIAtpar(id,1)*p))/ &
            &(exp(CIAtpar(id,1)*p) + exp(CIAtpar(id,1)*p))
    end select
  end function CircTimeArea_pvect

  !##################################################
  ! this function provides the time-dependent portion
  ! of a circular inclusion's BOUNDARY HEAD/FLUX response. 
  ! solution comes from convolving 
  ! this with the instantan. solution (* in Lap domain)

   function CircTimeBdry_pscal(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : CIBdryTime, CIBtpar
    implicit none

    integer, intent(in) :: id
    complex(DP),  intent(in) :: p
    complex(DP) :: mult
    
    select case(CIBdryTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-CIBtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-CIBtpar(id,1)*p)/p - &
            & exp(-CIBtpar(id,2)*p)/p
    case(3)
       ! instantaneous
       mult = exp(-CIBtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE/(p - p*exp(-CIBtpar(id,1)*p)) * &
            & (RONE - exp(-CIBtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = exp(-CIBtpar(id,2)*p)/(p + p*exp(-CIBtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = exp(-CIBtpar(id,2)*p)/(p**2 + CIBtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = exp(-CIBtpar(id,2)*p)/ p**2 * &
            &(exp(CIBtpar(id,1)*p) - exp(CIBtpar(id,1)*p))/ &
            &(exp(CIBtpar(id,1)*p) + exp(CIBtpar(id,1)*p))
    end select
  end function CircTimeBdry_pscal

  function CircTimeBdry_pvect(id,p) result(mult)
    use constants, only: DP, RONE
    use element_specs, only : CIBdryTime,CIBtpar
    implicit none

    integer, intent(in) :: id
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    
    select case(CIBdryTime(id))
    case(1)
       ! step on at time=par1
       mult = exp(-CIBtpar(id,1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult = exp(-CIBtpar(id,1)*p)/p - &
            & exp(-CIBtpar(id,2)*p)/p
    case(3)
       ! instantaneous
       mult = exp(-CIBtpar(id,1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult = RONE/(p - p*exp(-CIBtpar(id,1)*p)) * &
            & (RONE - exp(-CIBtpar(id,2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       mult = exp(-CIBtpar(id,1)*p)/(p + p*exp(-CIBtpar(id,1)*p))
    case(6)
       ! sine wave  - sin(par1*t)/par1
       mult = exp(-CIBtpar(id,1)*p)/(p**2 + CIBtpar(id,1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       mult = exp(-CIBtpar(id,1)*p)/p**2 * &
            &(exp(CIBtpar(id,1)*p) - exp(CIBtpar(id,1)*p))/ &
            &(exp(CIBtpar(id,1)*p) + exp(CIBtpar(id,1)*p))
    end select
  end function CircTimeBdry_pvect

end module wells
