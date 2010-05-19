! this module contains the basic functions defining the head or flux 
! effects of a circular or elliptical element, as well as the time
! behaviors, and the routine for computing general kappa for MHE.

module elements
  implicit none
  private

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface CircleHead
     module procedure Circle_Head_match_self, Circle_Head_match_other, Circle_Head_calc
  end interface
  interface CircleFlux
     module procedure CircleFlux_match, CircleFlux_calc
  end interface

  interface EllipseHead
     module procedure EllipseHead_match, EllipseHead_calc
  end interface
  interface EllipseFlux
     module procedure EllipseFlux_match, EllipseFlux_calc
  end interface

  interface Time
     module procedure Time_pScal, Time_pVect
  end interface
  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains

  function circle_head_match_self(c,p,dom) result(res)
    use constants, only : DP, PI
    use element_specs, only : circle, domain
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    type(domain), intent(in) :: dom
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(c%M,c%n+1,size(p,1)) :: res    
    integer :: nP

    cop = cos(outer_prod())

    res(1:c%M,)

  end function circle_head_match_self
  

  function circle_head_match_other(c,r,p,dom,in) result(res)
    use constants, only : DP, PI
    use element_specs, only : circle, domain
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    type(domain), intent(in) :: dom
    complex(DP), dimension(:), intent(in) :: p
    logical, intent(in) :: in ! calculating inside element?

    ! size = number of matching locations on target element
    real(DP), dimension(:), intent(in) :: r
    ! second dimension is either LHS + RHS 
    ! if Theis well specified flux source, n=0,
    ! and solution is the RHS vector for each p value (known),
    ! otherwise solution is LHS, a 2D matrix for each value of p (unknown)
    complex(DP), dimension(size(r,1),c%n+1,size(p,1)) :: res
    
    complex(DP), allocatable :: besk(:,:,:), besi(:,:,:), kap(:)
    integer :: nR,nP

    nR = size(r,1)
    nP = size(p,1)

    if (c%n <= 1) then
       if (.not. c%storIn) then
          ! finite-radius "Theis" well (n==0) with specified strength
          ! well given in van Everdingen & Hurst, no wellbore storage.

          allocate(besk(nR,nP,0:1),kap(nP))
          kap(1:nP) = kappa(p,c%matching)
          besk = bK(outerprod(kap(1:nP)*r(1:nR)),2)

          ! result is LHS vector (unknown strength) for each value of p
          res(1:nR,1,1:nP) = spread(Time(p,c%time,.false.)/kap(1:nP),dim=1,ncopies=nR)*&
               & besk(:,:,0)/(2.0*PI*c%r*besk(:,:,1))
          deallocate(besk,kap)
       
          if (c%n == 0) then
             ! result is RHS vector (known strength) for each value of p
             res(:,1,:) = C%spec*res(:,1,:)
          end if
       else
          ! finite-radius well with wellbore storage
          
       end if
    end if
    
    
    if (something) then
    else
       if (in) then
          


       else

       end if
    end if
    

  end function circle_head_match_other

  !##################################################
  ! does some error checking and calculates Bessel functions 
  ! needed for head effects of instantaneous point source
  ! radius is a vector, p and kappa are scalars (matching)

   function CircleHead_rvect(id,p,kappa,r) result(Lhead)
    use constants, only: DP, PI
    use element_specs, only : circle
    
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

  ! ##################################################
  ! general time behavior function, for all elements
   function Time_pvect(p,t,area) result(mult)
    use constants, only: DP
    use element_specs, only : time
    implicit none

    type(time), intent(in) :: t
    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p,1)) :: mult
    integer :: n, np
    real(DP), allocatable :: ti(:),q(:)
    real(DP) :: tf
    
    np = size(p,1)
    if (area) then
       time = t%areaTime
       allocate(par(size(t%ATPar)))
       par = t%AtPar
    else
       time = t%BdryTime
       allocate(par(size(t%BTPar)))
       par = t%BTPar
    end if

    select case(time)
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
       n = abs(time)
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = par(1:n)
       tf = par(n+1)
       Q(0) = 0.0
       Q(1:n) = par(n+2:2*n+1)
       
       mult(1:np) = (sum(spread(Q(1:n) - Q(0:n-1),2,np)*&
            & exp(-spread(ti(1:n),2,np)*spread(p(1:np),1,n)),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p(:)))/p(:)
    end select
  end function Time_pvect

  ! wrapper to allow calling time routine with scalar p
  function Time_pscal(p,t,area) result(mult)
    type(time), intent(in) :: t
    complex(DP), intent(in) :: p
    complex(DP) :: mult

    mult = sum(Time_pVect([p],t,area))

  end function Time_pscal

  function kappa_pVect(p,m) result(q)
    use constants, only : DP, PI
    use element_specs, only : matching

    integer, parameter :: NTERMS = 200, MAXITER = 200
    
    complex(DP), intent(in), dimension(:) :: p
    complex(DP), dimension(size(p),0:CInum) :: q

    integer :: i, ni, np
    complex(DP), dimension(size(p),0:CInum) :: kap2
    complex(DP), dimension(size(p)) :: exp2z

    !! boulton stuff
!!$    real(DP), dimension(NTERMS) :: guess, gamma
!!$    real(DP), dimension(0:CInum) :: sigma
!!$    complex(DP), dimension(size(p)) :: kernel
!!$    real(DP) :: x, delta
!!$    integer :: k, kk
!!$    logical, save :: first = .true.
    real(DP) :: boulton

    np = size(p)
    ni = CInum
!!$    sigma(0:ni) = sqrt(sv/Syv)

    do i=0,ni
       !! leaky-ness
       !! ##############################
       if(leakv(i) == 0) then
          !! no leaky layer, standard definition
          q(1:np,i) = p(1:np)/av(i)
       else
          kap2(1:np,i) = sqrt(p(:)/a2v(i))
          exp2z(1:np) = exp(-2.0_DP*kap2(:,i)*b2v(i))

          if(leakv(i) == 1) then
             !! case I, no-drawdown condition at top of aquitard
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))*&
                  & (1.0_DP + exp2z(:))/(1.0_DP - exp2z(:))
          elseif(leakv(i) == 2) then
             !! case II, no-flow condition at top of aquitard
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))*&
                  & (1.0_DP - exp2z(:))/(1.0_DP + exp2z(:))
          elseif(leakv(i) == 3) then
             !! aquitard thickness -> infinity
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))
          else
             stop 'ERROR: incorrect value for CIAquitardLeak parameter -> (1,2,3)'
          end if
       end if
       
       !! unconfined-ness 
       !! ##############################
       if(unconfv(i) == 0) then
          !! do nothing, q already computed above
       else
          !! Boulton unconfined source (Herrera infinite sum Kernel)
          !! guess is halfway between asymptotes of cot()
!!$
!!$          if (first) then
!!$             if(.not. allocated(root)) then
!!$                allocate(root(NTERMS,0:CInum))
!!$             end if
!!$             
!!$             !! roots are not a function of p (just sigma)- only compute once
!!$             guess(2:NTERMS) = PI*(real((/(k, k=1,NTERMS-1)/)) + 0.5_DP)/sigma(i)
!!$             guess(1) = 1.7D0
!!$
!!$             !! first root is hard to find with NR, 
!!$             !! use TS approximation for tangent and re-arrange
!!$             x = guess(1)
!!$             NR1: do kk = 1,MAXITER
!!$                delta = (x + (sigma(i) - 1.0_DP/sigma(i))*(x*sigma(i) + (x*sigma(i))**3/3.0_DP + &
!!$                     & 2.0_DP*(sigma(i)*x)**5/15.0_DP) + 17.0_DP*(x*sigma(i))**7/315.0_DP)/ &
!!$                     & (1.0_DP - (1.0_DP/sigma(i) - sigma(i))*(sigma(i) + x**2*sigma(i)**3 + &
!!$                     & 2.0_DP*x**4*sigma(i)**5/3.0_DP + 17.0_DP*x**6*sigma(i)**7/45.0_DP))
!!$                x = x - delta
!!$                if (abs(delta) <= 1.0D-10) then
!!$                   root(1,i) = x
!!$                   exit NR1
!!$                end if
!!$                if(kk == MAXITER) print *, '1 failed to converge'
!!$             end do NR1
!!$
!!$             do k = 2, NTERMS
!!$                x = guess(k)
!!$                NR: do kk = 1,MAXITER
!!$                   delta = (1.0_DP/tan(x*sigma(i)) + (sigma(i) - 1.0_DP/sigma(i))/x)/&
!!$                        & (sigma(i)/(sin(sigma(i)*x)**2) + (sigma(i) + 1.0_DP/sigma(i))/x**2)
!!$                   x = x + delta
!!$                   if (abs(delta) <= spacing(x)*10.0) then
!!$                      root(k,i) = x
!!$                      exit NR
!!$                   end if
!!$                   if(kk == MAXITER) print *, k,'failed to converge'
!!$                end do NR
!!$             end do
!!!!!$             if (i == ni) first = .false.
!!!!!$          end if
!!$
!!$          gamma(1:NTERMS) = Kzv(i)*root(1:NTERMS,i)**2/(bgb*Syv(i))
!!$
!!$          kernel(1:np) = 2.0_DP*sum(spread(gamma(1:NTERMS),2,np)/&
!!$               & ((spread(root(1:NTERMS,i)**2,2,np) - 1.0_DP + sigma(i)**2)* &
!!$               & (spread(p(1:np),1,NTERMS) + spread(gamma(1:NTERMS),2,np))),dim=1)
!!$          
!!$
!!$          q(1:np,i) = q(1:np,i) + p(1:np)*Syv(i)/Kv(i)*kernel

          ! scrap Herrera's infinite sum for Boulton's original
          ! rough-n-ready alpha, with a semi-physical expression for it
          boulton = 3.0_DP*Kzv(i)/(bgb*Syv(i))
          q(1:np,i) = q(:,i) + Syv(i)*p(1:np)*boulton/(Kv(i)*(boulton + p(1:np)))
       end if
    end do
    
    !! sources are only additive under the square root
    q = sqrt(q);

  end function compute_CIleaky_qv
  
  !! scalar version useful in matching
  function compute_CIleaky_qs(p) result(q)
    use constants, only : DP
    use element_specs, only : CInum

    complex(DP), intent(in) :: p
    complex(DP), dimension(0:CInum) :: q

    !! sum away singleton first dimension
    kappa = sum(compute_CIleaky_qv( [p],m))

  end function compute_CIleaky_qs
  

end module elements
