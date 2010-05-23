! this module contains the basic functions defining the head or flux 
! effects of a circular or elliptical element, as well as the time
! behaviors, and the routine for computing general kappa for MHE.

module elements
  implicit none
  private

  ! the four basic functions are overloaded for either 
  ! p being a vector (during inversion) or a scalar (during matching)

  interface circle_head
     module procedure circle_match_head_self, &
          & circle_match_head_other, circle_head_calc
  end interface
  interface circle_flux
     module procedure circle_match_flux_self, &
          & circle_match_flux_other, circle_flux_calc
  end interface

  interface ellipse_head
     module procedure ellipse_match_head_self, &
          & ellipse_match_head_other, ellipse_head_calc
  end interface
  interface ellipse_flux
     module procedure ellipse_match_flux_self, &
          & ellipse_match_flux_other, ellipse_flux_calc
  end interface

  interface time
     module procedure Time_pScal, Time_pVect
  end interface
  interface kappa
     module procedure  kappa_pVect, kappa_pScal
  end interface

contains
  subroutine circle_match_head_self(c,p,dom,LHS,RHS)
    use constants, only : DP, PI
    use utilities, only : outer_prod
    use type_definitions, only : circle, domain
    implicit none

    type(circle), intent(in) :: c
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p

    ! LHS dim=2 is 4N-2 for matching, 2N-1 for spec. total head
    complex(DP), intent(out), &
         & dimension(c%M,(2*N-1)*(2-abs(c%ibnd))) :: LHS
    complex(DP), intent(out), dimension(c%M) :: RHS

    complex(DP) :: dimension(c%M,2*N-1) :: tmp
    integer :: j, N, M
    complex(DP) :: kap
    real(DP) :: cmat(1:M,0:N-1), smat(1:M,1:N-1)

    N = c%N; M = c%M
    cmat = cos(outer_prod(c%Pcm(1:M), real([(j,j=0,N-1)],DP)))
    smat = sin(outer_prod(c%Pcm(1:M), real([(j,j=1,N-1)],DP)))

    ! setup LHS
    ! matching or specified head (always first M rows); no dependence on p
    ! ibnd==-2 would be here, but it doesn't actually make physical sense
    if (c%ibnd==0 .or. c%ibnd==-1) then

       LHS(1:M,1:N) =       cmat/c%parent%K
       LHS(1:M,N+1:2*N-1) = smat/c%parent%K
       
       ! setup RHS
       select case(c%ibnd)
       case(-1)
          ! put specified head on RHS
          RHS(1:M) = time(p,c%time,.false.)*c%bdryQ
       case(0)
          ! put constant area source term effects on RHS
          RHS(1:M) = -time(p,c%time,.true.)*c%areaQ*c%Ss/kappa(p,c%parent)**2
       end select

       if (c%ibnd==0 .or. c%calcin) then
          LHS(1:M,2*N:3*N-1) = -cmat/c%K
          LHS(1:M,3*N:4*N-2) = -smat/c%K
       end if
    else
       LHS = 0.0
       RHS = 0.0
    end if
  end subroutine circle_match_head_self

  subroutine circle_match_flux_self(c,p,dom,LHS,RHS)
    use constants, only : DP, PI
    use utilities, only : outer_prod
    use type_definitions, only : circle, domain
    use bessel_functions, only : bK, bI
    implicit none

    type(circle), intent(in) :: c
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p

    complex(DP), intent(out), &
         & dimension(c%M,(2*N-1)*(2-abs(c%ibnd))) :: LHS
    complex(DP), intent(out), dimension(c%M) :: RHS

    complex(DP) :: dimension(c%M,2*N-1) :: tmp
    integer :: j, N, M
    complex(DP), allocatable :: Kn(:), dKn(:), In(:), dIn(:)
    complex(DP) :: kap
    real(DP) :: cmat(1:M,0:N-1), smat(1:M,1:N-1)

    N = c%N; M = c%M
    cmat = cos(outer_prod(c%Pcm(1:M), real([(j,j=0,N-1)],DP)))
    smat = sin(outer_prod(c%Pcm(1:M), real([(j,j=1,N-1)],DP)))

    ! matching (second M) or specified flux (first M); depends on p
    if (c%ibnd==0 .or. c%ibnd==1 .or. c%ibnd==2) then
       allocate(Kn(0:N),dKn(0:N))
       kap = kappa(p,c%parent) 
       call bKD(kap*c%r,N+1,Kn,dKn)

       tmp(1:M,1:N) =       spread(kap*dKn(0:N-1)/Kn(0:N-1), 1,M)*cmat/c%parent%K
       tmp(1:M,N+1:2*N-1) = spread(kap*dKn(1:N-1)/Kn(1:N-1), 1,M)*smat/c%parent%K
       deallocate(Kn,dKn)

       select case(c%ibnd)
       case(2)
          allocate(Kn(0:1))
          Kn(0:1) = bK(kap*c%r,2)
          
          if (c%StorIn) then
             ! effects of wellbore storage and skin on finite-radius
             ! well, where a_0 is computed (generally depends on
             ! other elements, too; these show up in off-diagonal sub-matrices)
             LHS(1:M,1) = -((2.0 + c%r**2*dskin*p/c%parent%T)/(2.0*PI*c%r) + &
                  & (Kn(0)*c%r*p)/(2.0*PI*c%r*kap*Kn(1)*c%parent%T))
             RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(PI*c%r*c%parent%T)
          else
             ! specified flux (finite-radius well no storage)
             ! a_0 coefficient is computed analytically
             LHS(1:M,1) = 0.0
             RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r*Kn(1))*tmp(1:M,1)
          end if
          deallocate(Kn)
       case(1)
          ! put specified flux effects on RHS
          LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          RHS(1:M) = time(p,c%time,.false.)*c%bdryQ/(2.0*PI*c%r)
       case(0)
          ! no area source term effects on flux matchinig
          LHS(1:M,1:N) = tmp(M+1:2*M,1:2*N-1)
          RHS(1:M) = 0.0 
       end if
    
       if (c%ibnd==0 .or. c%calcin) then
          allocate(In(0:N),dIn(0:N))
          kap = kappa(p,c%element)
          call bID(kap*c%r,N+1,In,dIn)
          
          LHS(1:M,2*N:3*N-1) = spread(kap*dIn(0:N-1)/In(0:N-1), 1,M)*cmat/c%K
          LHS(1:M,3*N:4*N-2) = spread(kap*dIn(1:N-1)/In(1:N-1), 1,M)*smat/c%K
          deallocate(dIn,In)
       end if
    else
       LHS = 0.0
       RHS = 0.0
    end if
  end subroutine circle_match_flux_self


  ! ##################################################
  ! general time behavior function, for all elements
   function Time_pvect(p,t,area) result(mult)
    use constants, only: DP
    use type_definitions, only : time
    implicit none

    type(time), intent(in) :: t
    complex(DP), dimension(:), intent(in) :: p
    logical, intent(in) :: area
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
    use type_definitions, only : element

    complex(DP), intent(in), dimension(:) :: p
    type(element), intent(in) :: el
    complex(DP), dimension(size(p)) :: q

    integer :: i, ni, np
    complex(DP), dimension(size(p)) ::  kap2, exp2z
    complex(DP) :: boulton

    np = size(p); ni = CInum

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

    !! integrate neuman 72 solution and include here instead

    !! unconfined-ness (if confined do nothing)
    !! ##############################
    if(el%unconfinedFlag) then
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
    use type_definitions, only : element

    complex(DP), intent(in) :: p
    type(element), intent(in) :: el
    complex(DP) :: q

    !! sum away singleton first dimension
    q = sum(kappa_pVect([p],el))

  end function kappa_pscal

end module elements
