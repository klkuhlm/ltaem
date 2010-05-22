! this module contains the basic functions defining the head or flux 
! effects of a circular or elliptical element, as well as the time
! behaviors, and the routine for computing general kappa for MHE.

module elements
  implicit none
  private

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface Circle_Effect
     module procedure Circle_match_self, Circle_match_other, Circle_calc
  end interface

  interface Ellipse_Effect
     module procedure Ellipse_match_self, Ellipse_match_other, Ellipse_calc
  end interface

  interface Time
     module procedure Time_pScal, Time_pVect
  end interface
  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains

  subroutine circle_match_self(c,p,dom,LHS,RHS)
    use constants, only : DP, PI
    use utilities, only : outer_prod
    use element_specs, only : circle, domain
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p

    ! ibnd={-2,-1,0,1,2}: spec {el,tot} flux, match, spec {tot,el} head
    
    ! LHS dim=1 is 2M for matching, M for spec. total head/flux,
    !   and 0 for specified elemental head/flux
    ! LHS dim=2 is 4N-2 for matching, 2N-1 for spec. total head/flux
    !   and 0 for specified elemental head/flux
    complex(DP), intent(out), &
         & dimension(c%M*(2-abs(c%ibnd)), &
         &       (2*N-1)*(2-abs(c%ibnd))) :: LHS
    complex(DP) :: dimension(c%M,2*N-1) :: tmp

    ! RHS is 2M for matching, M for spec. total head/flux,
    !   and M for spec. elemental head/flux
    complex(DP), intent(out), &
         & dimension(c%M*(2-count([c%ibnd]/=0))) :: RHS

    integer :: nP, j, N, M, lo, hi
    real(DP), dimension(0:c%N) :: vi, k(0:1)
    complex(DP), allocatable :: Kn(:), dKn(:), In(:), dIn(:)
    complex(DP) :: kap

    real(DP) :: cmat(1:M,0:N-1), smat(1:M,1:N-1)

    N = c%N; M = c%M
    forall(j=0,N) vi(j)=real(j,DP)

    k(0) = c%K   ! K (hyd. cond.) inside
    k(1) = c%parent%K ! K of parent

    cmat = cos(outer_prod(c%Pcm(1:M),vi(0:N-1)))
    smat = sin(outer_prod(c%Pcm(1:M),vi(1:N-1)))
    RHS = 0.0

    ! setup LHS
    ! matching or specified head (always first M rows); no dependence on p
    if (c%ibnd==0 .or. c%ibnd==-1) then

       LHS(1:M,1:N) =       cmat/c%parent%K
       LHS(1:M,N+1:2*N-1) = smat/c%parent%K
       
       ! assume always calculate outside an element.  
       if (c%calcin) then
          LHS(1:M,2*N:3*N-1) = -cmat/c%K
          LHS(1:M,3*N:4*N-2) = -smat/c%K
       end if
    end if

    ! matching (second M) or specified flux (first M); depends on p
    if (c%ibnd==0 .or. c%ibnd==1 .or. c%ibnd==2) then
       lo = M+1-c%ibnd*M; hi = 2*M-c%ibnd*M
       allocate(Kn(0:N),dKn(0:N))
       kap = kappa(p,c%parent) 
       call bKD(kap*c%r,N+1,Kn,dKn)

       tmp(1:M,1:N) =       &
            & spread(kap*dKn(0:N-1)/Kn(0:N-1), 1,M) * cmat/c%parent%K
       tmp(1:M,N+1:2*N-1) = &
            & spread(kap*dKn(1:N-1)/Kn(1:N-1), 1,M) * smat/c%parent%K
       deallocate(Kn,dKn)

       if (c%ibnd==2) then
          allocate(Kn(0:1))
          Kn(0:1) = bK(kap*c%r,2)
          
          if (.not. c%StorIn) then
             ! specified flux (finite-radius well no storage)
             ! a_0 coefficient is computed analytically
             LHS(1:M,1) = 0.0
             RHS(1:M) = time(p,c%time)*c%bdryQ/(2.0*PI*c%r*Kn(1)) * &
                  & tmp(1:M,1)
          else
             ! effects of wellbore storage and skin on finite-radius
             ! well, where a_0 is computed (generally depends on
             ! other elements, too)
             LHS(1:M,1) = -((2.0 + c%r**2*dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                  & (Kn(0)*c%r*p)/(2.0*PI*c%r*kap*Kn(1)*c%parent%T))
             RHS(1:M) = time(p,c%time)*c%bdryQ/(PI*c%r*c%parent%T)
          end if
          deallocate(Kn)
       else
          LHS(lo:hi,1:N) = tmp(M+1:2*M,1:2*N-1)
       end if
    
       if (c%calcin) then
          allocate(In(0:N),dIn(0:N))
          kap = kappa(p,c%element)
          call bID(kap*c%r,N+1,In,dIn)
          
          LHS(lo:hi,2*N:3*N-1) = &
               & spread(kap*dIn(0:N-1)/In(0:N-1), 1,M) * cmat/c%K
          LHS(lo:hi,3*N:4*N-2) = &
               & spread(kap*dIn(1:N-1)/In(1:N-1), 1,M) * smat/c%K
          deallocate(dIn,In)

       end if
    end if
  
    if(c%ibnd == 0) then
       ! area source
       RHS(1:M) = 
    end if
    
  end subroutine circle_match_self

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
       ti(1:n) = par(1:n); tf = par(n+1)
       Q(0) = 0.0; Q(1:n) = par(n+2:2*n+1)
       
       mult(1:np) = (sum(spread(Q(1:n) - Q(0:n-1),2,np)*&
            & exp(-spread(ti(1:n),2,np)*spread(p(1:np),1,n)),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p(:)))/p(:)
    case default
       write(*,'(A,I0)') 'Time_pvect: error in case for time behavior ',time
       stop
    end select
  end function Time_pvect

  ! wrapper to allow calling time routine with scalar p
  function Time_pscal(p,t) result(mult)
    type(time), intent(in) :: t
    complex(DP), intent(in) :: p
    complex(DP) :: mult

    mult = sum(Time_pVect([p],t,area))

  end function Time_pscal

  function kappa_pVect(p,el) result(q)
    use constants, only : DP, PI
    use element_specs, only : element

    integer, parameter :: NTERMS = 200, MAXITER = 200
    
    complex(DP), intent(in), dimension(:) :: p
    type(element), intent(in) :: el
    complex(DP), dimension(size(p)) :: q

    integer :: i, ni, np
    complex(DP), dimension(size(p),0:CInum) :: kap2
    complex(DP), dimension(size(p)) :: exp2z
    complex(DP) :: boulton

    np = size(p)
    ni = CInum

    !! leaky-ness
    !! ##############################
    if(el%leakFlag == 0) then
       !! no leaky layer, standard definition
       q(:) = p(1:np)/el%alpha
    else
       kap2(1:np) = sqrt(p(:)*el%aquitardSs/el%aquitardK)
       exp2z(1:np) = exp(-2.0*kap2(:)*el%aquitardb)
       
       if(el%leakFlag == 1) then
          !! case I, no-drawdown condition at top of aquitard
          q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
               & (1.0 + exp2z(:))/(1.0 - exp2z(:))
       elseif(el%leakFlag == 2) then
          !! case II, no-flow condition at top of aquitard
          q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)*&
               & (1.0 - exp2z(:))/(1.0 + exp2z(:))
       elseif(el%leakFlag == 3) then
          !! aquitard thickness -> infinity
          q(:) = p(:)/el%alpha + kap2(:)*el%aquitardK/(el%b*el%K)
       else
          stop 'ERROR: incorrect value for leakFlag parameter -> (1,2,3)'
       end if
    end if
    
    !! unconfined-ness 
    !! ##############################
    if(el%unconfinedFlag) then
       !! do nothing, q already computed above
    else
       !! Boulton unconfined source (Herrera infinite sum Kernel)
       !! guess is halfway between asymptotes of cot()
       
       ! scrap Herrera's infinite sum for Boulton's original
       ! rough-n-ready alpha, with a semi-physical expression for it
       boulton = 3.0*el%Kz/(el%b*el%Sy)
       q(1:np) = q(:) + el%Sy*p(:)*boulton/(el%K*(boulton + p(:)))
    end if
    
    !! sources are only additive under the square root
    q = sqrt(q)

  end function kappa_pVect
  
  !! scalar version useful in matching
  function kappa_pscal(p,el) result(q)
    use constants, only : DP
    use element_specs, only : element

    complex(DP), intent(in) :: p
    type(element), intent(in) :: el
    complex(DP) :: q

    !! sum away singleton first dimension
    q = sum(kappa_pVect([p],el))

  end function kappa_pscal

end module elements
