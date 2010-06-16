

module calc_routines
  use constants, only : DP
  implicit none

  private
  public :: headCalc, velxCalc, velyCalc

  ! assuming velx to be called before vely - much of calculation will be saved
  ! in these module-wide variables

  complex(DP), save, allocatable :: CIdPot_dAng(:,:,:), CIdPot_dR(:,:,:)

contains

  !##################################################
  function headCalc(Z,p,lo,hi,sol,dom,c,e,bg) result(H)
    use type_definitions, only : solution, domain, circle, ellipse
    use constants, only : DP, PI
    use bessel_functions, only : bK, bI
    use mathieu_functions, only : Ke, Ko, Ie, Io, ce, se
    use kappa_mod
    use time_mod

    complex(DP), intent(in) :: Z  ! location for calculation (complex coordinates)
    complex(DP), dimension(:), intent(in) :: p  ! vector of Laplace parameters
    integer, intent(in) :: lo,hi  ! lower and upper bounds of p relative to overall s
    type(solution), intent(inout) :: sol 
    type(domain), intent(in) :: dom
    type(circle), dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(element), intent(in) :: bg

    complex(DP), dimension(size(p,1)) :: H  
    complex(DP), dimension(size(p,1),0:CIn) :: PotInclPar
    complex(DP), dimension(size(p,1),0:CIn, CInum) :: PotInclBg, PotInclIn
    complex(DP), dimension(size(p,1),WLnum) :: PotWell
    real(DP), dimension(CInum) :: Rcp, Pcp
    real(DP), dimension(WLnum) :: Rwp, Pwp
    complex(DP), allocatable :: aa(:,:),bb(:,:),cc(:,:),dd(:,:) ! local coefficients
    complex(DP), allocatable :: RMRgp(:,:,:), RMR0(:,:,:), AM(:,:,:) ! radial / angular MF

    real(DP), allocatable :: vr(:)
    integer, allocatable :: vi(:)
    integer :: nc, ne, ntot, np, N, j, i
    
    nc = dom%num(1); ne = dom%num(2)
    ntot = nc + ne
    np = size(p,1)
    if (hi-lo+1 /= np) then
       write(*,'(A,3(1X,I0))') 'ERROR: lo,hi do not match dimensions of p',lo,hi,size(p,1)
       stop
    end if
    
    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(Z,e,c,dom,Rgp,Pgp,inside) 
    
    !##################################################
    !! calculation point is outside all elements (in background)
    if (inside == 0) then
       H = 0.0
       do j = 1,nc
          if (dom%InclUp(j) == 0) then  ! circle is also in background
             N = c(j)%N
             kap = kappa(p(:),c%parent%element)
             allocate(vr(0:N),aa(np,N),bb(np,N))
             vr = real([(i,i=0,N-1)],DP)
             aa(1:np,1:N) = c(j)%A(lo:hi,1:N)
             bb(1:np,1) = 0.0 ! insert zero to make odd/even same shape
             bb(1:np,2:N) = c(j)%A(lo:hi,N+1:2*N-1)

             H(1:np) = H + sum(bK(Rgp(j)*kap,N) / bK(c(j)%r*kap,N)* &
                  & ( aa(1:N)*spread(cos(vr(0:N-1)*Pgp(j)),1,np) + &
                  &   bb(1:N)*spread(sin(vr(0:N-1)*Pgp(j)),1,np) ),dim=2)
             deallocate(vr,aa,bb)
          end if
       end do
       
       do j = 1,ne
          if (dom%InclUp(nc+j) == 0) then  ! ellipse is also in background
             allocate(vi(0:N),RMRgp(np,N,0:1),RMR0(np,N,0:1), &
                  & AM(np,N,0:1),aa(np,N),bb(np,N))
             vi = [(i,i=0,N-1)]
             aa(1:np,1:N) = e(j)%A(lo:hi,1:N)
             bb(1:np,2:N) = e(j)%A(lo:hi,N+1:2*N-1)

             do i=1,np
                RMRgp(i,0:N-1,0) = Ke(bg%mat(lo+i-1),vi(0:N-1),Rgp(j))
                RMRgp(i,1:N-1,1) = Ko(bg%mat(lo+i-1),vi(1:N-1),Rgp(j))
                 RMR0(i,0:N-1,0) = Ke(bg%mat(lo+i-1),vi(0:N-1),e(j)%r)
                 RMR0(i,1:N-1,1) = Ko(bg%mat(lo+i-1),vi(1:N-1),e(j)%r)
                   AM(i,0:N-1,0) = ce(bg%mat(lo+i-1),vi(0:N-1),Pgp(j))
                   AM(i,1:N-1,1) = se(bg%mat(lo+i-1),vi(0:N-1),Pgp(j))
             end do

             H(1:np) = H + &
                  & sum(RMRgp(:,0:N-1,0)/RMR0(:,0:N-1,0)*aa(:,1:N)*AM(:,0:N-1,0), 2) + &
                  & sum(RMRgp(:,1:N-1,1)/RMR0(:,1:N-1,1)*bb(:,2:N)*AM(:,1:N-1,1), 2)
             deallocate(vi,RMRgp,RMR0,AM,aa,bb)
          end if
       end do

       ! result for invlap (head, not discharge potential)
       H(1:np) = H(:)/bg%k

    !##################################################
    !! calculation point is inside (or on bdry of) an inclusion
    else
       PotInclPar = CZERO; PotInclIn = CZERO
       if (CIcalcin(inside)) then
          
          besiRcp(1:np,0:N) = besseli(Rcp(inside)*q(1:np,inside),0,N+1)
          besiR0(1:np,0:N) = besseli(CIr(inside)*q(1:np,inside),0,N+1)
          
          ! "parent" inclusion, which calculation point is within

          PotInclPar(1:np,0:N) = besiRcp(1:np,0:N)/besiR0(1:np,0:N)* &
               & ( c(1:np,0:N,inside)*spread(cos(rk(0:N)*Pcp(inside)),1,np) + &
               &   d(1:np,0:N,inside)*spread(sin(rk(0:N)*Pcp(inside)),1,np) )

          ! effects of any inclusions which may be inside same parent too
          do other = 1,ni
             if (CIInclIn(inside,other)) then
                beskRcp(1:np,0:N) = besselk(Rcp(other)*q(1:np,inside),0,N+1)
                beskR0(1:np,0:N) = besselk(CIr(other)*q(1:np,inside),0,N+1)

                PotInclIn(1:np,0:N,other) = beskRcp(1:np,0:N)/beskR0(1:np,0:N)* &
                     & ( a(1:np,0:N,other)*spread(cos(rk(0:N)*Pcp(other)),1,np) + &
                     &   b(1:np,0:N,other)*spread(sin(rk(0:N)*Pcp(other)),1,np) )
             end if
          end do
       end if  !if CIcalcin

       ! result for invlap (head, not discharge potential)
       H(1:np) = (sum(PotWell,dim=2) + sum(sum(PotInclIn,dim=3),dim=2) + sum(PotInclPar,dim=2) + &
            & CIarea(inside)*circTimeArea(inside,p(:))*sv(inside)/q(:,inside)**2)/kv(inside)
    end if ! inout
  end function headCalc

  !##################################################
  function velXCalc(p,calcX,calcY,calcCoeff) result(vx)

    use element_specs, only : CIn,CInum,WLnum,porv,CIWellIn,CIInclUp,CICalcIn,CIInclIn,CIr
    use constants, only : DP, PI, CZERO, CTWO
    use bessel_functions
    use wells, only : wellFlux
    use leaky_q

    complex(DP), dimension(:), intent(in) :: p
    real(DP), intent(in) :: calcX,calcY
    complex(DP), dimension(1:,0:,1:), intent(in) :: calcCoeff

    complex(DP), dimension(size(p,1)) :: vx
    complex(DP), dimension(size(p,1),0:CIn, 1:CInum) :: a,b,c,d
    complex(DP), dimension(size(p,1),0:CInum) :: q
    complex(DP), dimension(size(p,1)) :: FluxInclPar
    complex(DP), dimension(size(p,1),CInum) :: FluxInclBg, FluxInclIn
    complex(DP), dimension(size(p,1),WLnum) :: FluxWell
    complex(DP), dimension(size(p,1),-1:CIn+1) :: beskRcp, besiRcp
    complex(DP), dimension(size(p,1),0:CIn) :: beskR0, besiR0
    complex(DP), dimension(size(p,1),0:CIn) :: cosnPcp, sinnPcp
    real(DP), dimension(CInum) :: Rcp, Pcp
    real(DP), dimension(WLnum) :: Rwp, Pwp

    integer :: ni, nw, incl, inside, other, well, i
    integer :: np, N
    real(DP), dimension(0:CIn) :: rk

    ni = CINum; nw = WLnum
    np = size(p,1);

    ! number of terms in sin/bessel funciton expansion 
    N = CIn;

    forall (i = 0:N)
       rk(i) = real(i,DP)
    end forall

    q(1:np,0:ni) = compute_CIleaky_q(p(1:np))

    b(:,0,:) = CZERO; d(:,0,:) = CZERO;! zero out unused terms

    ! coefficients passed through module
    a(1:np, 0:N, 1:ni) = calcCoeff(1:np, 0:N, 1:ni)
    b(1:np, 1:N, 1:ni) = calcCoeff(1:np, N+1:2*N, 1:ni)  ! b_0 not included
    c(1:np, 0:N, 1:ni) = calcCoeff(1:np, 2*N+1:3*N+1, 1:ni)
    d(1:np, 1:N, 1:ni) = calcCoeff(1:np, 3*N+2:4*N+1, 1:ni)

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(CalcX,CalcY,inside,Rcp,Pcp,Rwp,Pwp)

    ! effects of wells inside or outside inclusions
    FluxWell = CZERO
    do well = 1,nw
       if (CIWellIn(inside,well)) then ! if well is in background
          FluxWell(1:np,well) = cos(Pwp(well))*WellFlux(well,p(1:np),q(1:np,inside),Rwp(well))
       end if
    end do

    !##################################################
    !! calculation point is outside all inclusions
    if (inside == 0) then

       if (.not. allocated(CIdPot_dAng)) then
          allocate(CIdPot_dAng(1:np,1:ni,0:N),CIdPot_dR(1:np,1:ni,0:N))
       end if

       ! zero out
       CIdPot_dAng = CZERO; CIdPot_dR = CZERO;  FluxInclBg = CZERO

       do incl = 1,ni
          if (CIInclUp(incl) == 0) then  ! inclusion is in background

             beskRcp(1:np,0:N+1) = besselk(Rcp(incl)*q(:,0),0,N+2)
             beskR0(1:np,0:N) = besselk(CIr(incl)*q(:,0),0,N+1)
             beskRcp(:,-1) = beskRcp(:,+1)             

             cosnPcp(1:np,0:N) = spread(cos(rk(0:N)*Pcp(incl)),dim=1,ncopies=np)
             sinnPcp(1:np,0:N) = spread(sin(rk(0:N)*Pcp(incl)),dim=1,ncopies=np)

             ! compute derivatives in local coordinates of circular element
             CIdPot_dR(1:np,incl,0:N) = -spread(q(:,0),dim=2,ncopies=N+1)* &
                  & (beskRcp(:,-1:N-1) + beskRcp(:,1:N+1))/(CTWO*beskR0(:,0:N))* &
                  & ( a(:,0:N,incl)*cosnPcp(:,0:N) + b(:,0:N,incl)*sinnPcp(:,0:N) )

             CIdPot_dAng(1:np,incl,0:N) = spread(rk(0:N),1,np)*beskRcp(:,0:N)/beskR0(:,0:N)* &
                  & ( b(:,0:N,incl)*cosnPcp(:,0:N) - a(:,0:N,incl)*sinnPcp(:,0:N) )
             
             ! project this flux into x-direction due to background inclusions
             FluxInclBg(1:np,incl) =  cos(Pcp(incl))*sum(CIdPot_dR(:,incl,:),dim=2) - & 
                  & sin(Pcp(incl))*sum(CIdPot_dAng(:,incl,:),dim=2)/Rcp(incl)

          end if
       end do
       
       ! result for invlap velocity in x direction
       vx(1:np) = -(sum(FluxInclBg,dim=2) + sum(FluxWell,dim=2))/porv(0)

    !##################################################
    !! calculation point is inside (or on bdry of) an inclusion
    else
       if (allocated(CIdPot_dAng)) deallocate(CIdPot_dAng,CIdPot_dR)
       allocate(CIdPot_dAng(1:np,1:ni,0:N),CIdPot_dR(1:np,1:ni,0:N))

       ! zero out
       CIdPot_dAng = CZERO; CIdPot_dR = CZERO;
       FluxInclIn = CZERO;  FluxInclPar = CZERO

       if (CIcalcin(inside)) then
          besiRcp(1:np,0:N+1) = besseli(Rcp(inside)*q(:,inside),0,N+2)
          besiR0(1:np,0:N) = besseli(CIr(inside)*q(:,inside),0,N+1)
          besiRcp(:,-1) = besiRcp(:,+1)
          
          cosnPcp(1:np,0:N) = spread(cos(rk(0:N)*Pcp(inside)),dim=1,ncopies=np)
          sinnPcp(1:np,0:N) = spread(sin(rk(0:N)*Pcp(inside)),dim=1,ncopies=np)

          CIdPot_dR(1:np,inside,0:N) = spread(q(:,inside),dim=2,ncopies=N+1)* &
               (besiRcp(:,-1:N-1) + besiRcp(:,1:N+1))/(CTWO*besiR0(:,0:N))* &
               ( c(:,0:N,inside)*cosnPcp(:,0:N) + d(:,0:N,inside)*sinnPcp(:,0:N) )
          CIdPot_dAng(1:np,inside,0:N) = spread(rk(0:N),1,np)*besiRcp(:,0:N)/besiR0(:,0:N)* &
               ( d(:,0:N,inside)*cosnPcp(:,0:N) - c(:,0:N,inside)*sinnPcp(:,0:N) )

          FluxInclPar(1:np) = cos(Pcp(inside))*sum(CIdPot_dR(1:np,inside,0:N),dim=2) - &
               & sin(Pcp(inside))*sum(CIdPot_dAng(1:np,inside,0:N),dim=2)/Rcp(inside)

          ! effects of any inclusions, which may be inside parent too
          do other = 1,ni
             if (CIInclIn(inside,other)) then
                beskRcp(1:np,0:N+1) = besselk(Rcp(other)*q(:,inside),0,N+2)
                beskR0(1:np,0:N) = besselk(CIr(other)*q(:,inside),0,N+1)
                beskRcp(:,-1) = beskRcp(:,+1)

                cosnPcp(1:np,0:N) = spread(cos(rk(0:N)*Pcp(other)),dim=1,ncopies=np)
                sinnPcp(1:np,0:N) = spread(sin(rk(0:N)*Pcp(other)),dim=1,ncopies=np)

                CIdPot_dR(1:np,other,0:N) =  -spread(q(:,inside),dim=2,ncopies=N+1)* &
                     (beskRcp(:,-1:N-1) + beskRcp(:,1:N+1))/(CTWO*beskR0(:,0:N))* &
                     ( a(:,0:N,other)*cosnPcp(:,0:N) + b(:,0:N,other)*sinnPcp(:,0:N) )
                CIdPot_dAng(1:np,other,0:N) = spread(rk(0:N),1,np)*beskRcp(:,0:N)/beskR0(:,0:N)* &
                     ( b(:,0:N,other)*cosnPcp(:,0:N) - a(:,0:N,other)*sinnPcp(:,0:N) )
                
                FluxInclIn(1:np,other) = cos(Pcp(other))*sum(CIdPot_dR(:,other,:),dim=2) &
                     & - sin(Pcp(other))*sum(CIdPot_dAng(:,other,:),dim=2)/Rcp(other)
             end if
          end do

       end if  !if CIcalcin

       ! result for invlap
       vx(1:np) = -(sum(FluxWell,dim=2) + sum(FluxInclIn,dim=2) + FluxInclPar)/porv(inside)
    end if
  end function velXCalc 

  !##################################################
  function velYCalc(p,calcX,calcY) result(vy)
    
    ! most of the calculations are actually done in the velXCalc() function 
    ! which must be called before this one
    use element_specs, only : CInum,WLnum,porv,CIWellIn,CIInclUp,CICalcIn,CIInclIn
    use constants, only : DP, PI, CZERO
    use wells, only : wellFlux
    use leaky_q

    complex(DP), dimension(:), intent(in) :: p
    real(DP), intent(in) :: calcX,calcY

    complex(DP), dimension(size(p,1)) :: vy
    complex(DP), dimension(size(p,1),0:CInum) :: q
    complex(DP), dimension(size(p,1),1:CInum) :: FluxInclBg,FluxInclIn
    complex(DP), dimension(size(p,1),1:WLnum) :: FluxWell
    complex(DP), dimension(size(p,1)) :: FluxInclPar
    real(DP), dimension(CInum) :: Rcp, Pcp
    real(DP), dimension(WLnum) :: Rwp, Pwp

    integer :: ni, nw, incl, inside, other, well
    integer :: np

    ni = CINum; nw = WLnum
    np = size(p,1);

    q(1:np,0:ni) = compute_CIleaky_q(p(1:np))

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(CalcX,CalcY,inside,Rcp,Pcp,Rwp,Pwp)

    ! well flux inside or outside
    FluxWell = CZERO
    do well = 1,nw
       if (CIWellIn(inside,well)) then
          FluxWell(1:np,well) = sin(Pwp(well))*WellFlux(well,p(:),q(:,inside),Rwp(well))
       end if
    end do

    !##################################################
    !! calculation point is outside all inclusions
    if (inside == 0) then
       FluxInclBg = CZERO
       
       forall (incl = 1:ni, CIInclUp(incl) == 0)
          
          ! flux in y-direction due to background inclusions
          FluxInclBg(1:np,incl) = sin(Pcp(incl))*sum(CIdPot_dR(:,incl,:),dim=2) + &
               & cos(Pcp(incl))*sum(CIdPot_dAng(:,incl,:),dim=2)/Rcp(incl)
       end forall

       ! result for invlap: velocity in y direction
       vy(1:np) = -(sum(FluxInclBg,dim=2) + sum(FluxWell,dim=2))/porv(0)

    !##################################################
    !! calculation point is inside (or on bdry of) an inclusion
    else
       FluxInclIn = CZERO;  FluxInclPar = CZERO

       if (CIcalcin(inside)) then

          FluxInclPar(1:np) = sin(Pcp(inside))*sum(CIdPot_dR(:,inside,:),dim=2) + &
               & cos(Pcp(inside))*sum(CIdPot_dAng(:,inside,:),dim=2)/Rcp(inside)

          ! effects of any inclusions, which may be inside parent too
          forall (other = 1:ni, CIInclIn(inside,other))
             FluxInclIn(1:np,other) = sin(Pcp(other))*sum(CIdPot_dR(:,other,:),dim=2) &
                     & + cos(Pcp(other))*sum(CIdPot_dAng(:,other,:),dim=2)/Rcp(other)
          end forall
       end if  !if I.calcin

       ! velocity in y-direction for invlap
       vy(1:np) = -(sum(FluxWell,dim=2) + sum(FluxInclIn,dim=2) + FluxInclPar(:))/porv(inside)
    end if
  end function velYCalc

  !##################################################
  subroutine calcLocation(Z,e,c,dom,Rgp,Pgp,inside)
    use constants, only : DP, EYE
    use type_definitions, only : circle, ellipse, domain
    use utility, only : ccosh, cacosh

    complex(DP), intent(in) :: Z
    type(circle), dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    real(DP), dimension(sum(dom%num)), intent(out) :: Rgp, Pgp
    integer, intent(out) :: inside

    complex(DP), dimension(sum(dom%num)) :: Zgp
    integer, dimension(sum(dom%num)) :: inout

    complex(DP) :: tmp
    integer :: nc, ne, ntot, j, first, second, third, k

    nc = dom%num(1); ne = dom%num(2)
    ntot = nc + ne

    ! determine if observation point is inside an inclusion
    inout(1:ntot) = 0
    k = 0;

    do j = 1, nc
       ! components of vector from center of circle to observation point
       Zgp(j) = Z - cmplx(c(j)%x,c(j)%y,DP)
       Rgp(j) = abs(Zgp(j))
       Pgp(j) = atan2(aimag(Zgp(j)), real(Zgp(j)))
       if (Rgp(j) <= c(j)%r) then    ! inside or on boundary
          k = k + 1;
          inout(k) = j;
       end if
    end do

    do j = nc+1, ntot
       ! components of vector from center of ellipse to observation point
       Zgp(j) = Z - cmplx(e(j-nc)%x,e(j-nc)%y,DP)
       tmp = cacosh(Zgp(j))*exp(-EYE*e(j-nc)%theta)/e(j-nc)%f
       Rgp(j) = real(tmp)
       Pgp(j) = aimag(tmp)
       if (Rgp(j) <= e(j-nc)%r) then    ! inside or on boundary
          k = k + 1;
          inout(k) = j;
       end if
    end do

    if (ntot > 0) then
       select case (k)
       case (0) ! inside nothing (background)
          inside = 0
       case (1) ! inside one element
          inside = inout(1)
       case (2) ! inside two elements
          if (dom%InclUp(inout(1)) == 0) then
             inside = inout(2);
          else
             inside = inout(1);
          end if
       case (3) ! inside three elements
          do first = 1,3
             do second = 1,3
                if (second /= first) then
                   do third = 1,3
                      if ((third /= first) .and. (third /= second)) then
                         if (dom%InclIn(inout(first),inout(second))  .and. &
                              & dom%InclIn(inout(third),inout(first))) then 
                            inside = inout(second)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       case default
          write(*,*) 'CALCLOCATION ERROR: more than triply-nested inclusions', dom%InclIn
          stop 
       end select
    else
       inside = 0
       Xcp = 0.0; Ycp = 0.0
       Rcp = 0.0; Pcp = 0.0
    end if

    ! move observation points inside ibnd==2 elements
    do j = 1,nc 
       if (Rgp(j) < c(j)%r) then
          Rgp(j) = c(j)%r
       end if
    end do
    do j = 1,ne 
       if (Rgp(nc+j) < e(j)%r) then
          Rgp(nc+j) = e(j)%r
       end if
    end do
  end subroutine  calcLocation
end module calc_routines
