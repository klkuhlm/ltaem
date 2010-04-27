! $Id: calc_routines.f90,v 1.6 2008/12/10 02:46:59 kris Exp kris $
module calc_routines
  use constants, only : DP
  implicit none

  private
  public :: headCalc, velxCalc, velyCalc, freecalcmem

  ! assuming velx to be called before vely - much of calculation will be saved
  ! in these module-wide variables

  complex(DP), save, allocatable :: CIdPot_dAng(:,:,:), CIdPot_dR(:,:,:)

contains

  !##################################################
  function headCalc(p,calcX,calcY,calcCoeff) result(H)
    use element_specs, only : CIn,CInum,WLnum,kv,sv,CIWellIn,CIInclUp,&
         & CIr,CICalcIn,CIInclIn,CIArea,sv
    use constants, only : DP, PI, CZERO, RZERO
    use bessel_functions
    use wells, only : wellHead, circTimeArea
    use leaky_q

    complex(DP), dimension(:), intent(in) :: p
    real(DP), intent(in) :: calcX,calcY
    complex(DP), dimension(1:,0:,1:), intent(in) :: calcCoeff

    complex(DP), dimension(size(p,1)) :: H  !! f(p)
    complex(DP), dimension(size(p,1), 0:CIn, CInum) :: a,b,c,d
    complex(DP), dimension(size(p,1),0:CInum) :: q
    complex(DP), dimension(size(p,1),0:CIn) :: PotInclPar
    complex(DP), dimension(size(p,1),0:CIn, CInum) :: PotInclBg, PotInclIn
    complex(DP), dimension(size(p,1),WLnum) :: PotWell
    complex(DP), dimension(size(p,1),0:CIn) :: beskRcp, beskR0, besiRcp, besiR0
    real(DP), dimension(CInum) :: Rcp, Pcp
    real(DP), dimension(WLnum) :: Rwp, Pwp
    
    integer :: ni, nw, incl, inside, other, well
    integer :: np, N, i
    Real(DP), dimension(0:CIn) :: rk


    ni = CINum; nw = WLnum
    np = size(p,1)

    ! number of terms in sin/bessel funciton expansion 
    N = CIn;

    forall (i = 0:N)
       rk(i) = real(i,DP)
    end forall
        
    q(1:np,0:ni) = compute_CIleaky_q(p(1:np))

    b(:,0,:) = CZERO; d(:,0,:) = CZERO;
    
    ! coefficients passed through module
    a(1:np, 0:N, 1:ni) = calcCoeff(1:np, 0:N, 1:ni)
    b(1:np, 1:N, 1:ni) = calcCoeff(1:np, N+1:2*N, 1:ni)      ! b_0 not included
    c(1:np, 0:N, 1:ni) = calcCoeff(1:np, 2*N+1:3*N+1, 1:ni)
    d(1:np, 1:N, 1:ni) = calcCoeff(1:np, 3*N+2:4*N+1, 1:ni)

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(CalcX,CalcY,inside,Rcp,Pcp,Rwp,Pwp) !Pwp not used here

    ! effects of wells inside or outside inclusion
    PotWell = CZERO
    do well = 1,nw 
       if (CIWellIn(inside,well)) then
          PotWell(1:np,well) = WellHead(well,p(:),q(:,inside),Rwp(well))
       end if
    end do
    
    !##################################################
    !! calculation point is outside all inclusions
    if (inside == 0) then
       PotInclBg = CZERO
       do incl = 1,ni
          if (CIInclUp(incl) == 0) then  ! if inclusion is in background

             beskR0(1:np,0:N) = besselk(CIr(incl)*q(:,0),0,N+1)
             beskRcp(1:np,0:N) = besselk(Rcp(incl)*q(:,0),0,N+1)

             PotInclBg(1:np,0:N,incl) = beskRcp(:,0:N)/beskR0(:,0:N)* &
                  & ( a(:,0:N,incl)*spread(cos(rk(0:N)*Pcp(incl)),1,np) + &
                  &   b(:,0:N,incl)*spread(sin(rk(0:N)*Pcp(incl)),1,np) )
          end if
       end do ! for ni

       ! result for invlap (head, not discharge potential)

       H(1:np) = (sum(sum(PotInclBg,dim=3),dim=2) + sum(PotWell,dim=2))/kv(0)

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
  subroutine calcLocation(CalcX,CalcY,inside,Rcp,Pcp,Rwp,Pwp)
    use constants, only : DP
    use element_specs, only : CInum, CIx, CIy, CIr, CIInclUp, CIInclIn, WLnum, &
         & WLx, WLy, CIWellIn,  WLr

    real(DP), intent(in) :: CalcX, CalcY
    integer, intent(out) :: inside
    real(DP), dimension(1:CInum), intent(out) :: Rcp, Pcp
    real(DP), dimension(1:WLnum), intent(out) :: Rwp, Pwp
    real(DP), dimension(1:CInum) :: Xcp, Ycp
    real(DP), dimension(1:WLnum) :: Xwp, Ywp
    integer :: incl, well, first, second, third, ni, nw, k
    integer, dimension(CInum) :: inout

    ni = CInum; nw = WLnum

    !! this routine should be replaced or share code with the similar routine
    !! in circular_geometry.f90 which is more complicated, but much more general than this

    ! determine if observation point is inside an inclusion
    inout(1:ni) = 0
    k = 0;

    do incl = 1,ni
       ! components of vector from center of inclusion to observation point
       Xcp(incl) = CalcX - CIx(incl)
       Ycp(incl) = CalcY - CIy(incl)
       Rcp(incl) = sqrt(Xcp(incl)**2 + Ycp(incl)**2)
       Pcp(incl) = atan2(Ycp(incl), Xcp(incl))
       if (Rcp(incl) <= CIr(incl)) then    ! inside or on boundary
          k = k + 1;
          inout(k) = incl;
       end if
    end do

    if (ni > 0) then
       select case (k)
       case (0) !background
          inside = 0
       case (1) !one inclusion
          inside = inout(1)
       case (2) !two inclusions 
          if (CIInclUp(inout(1)) == 0) then
             inside = inout(2);
          else
             inside = inout(1);
          end if
       case (3) !three inclusions
          do first = 1,3
             do second = 1,3
                if (second /= first) then
                   do third = 1,3
                      if ((third /= first) .and. (third /= second)) then
                         if (CIInclIn(inout(first),inout(second))  .and. &
                              & CIInclIn(inout(third),inout(first))) then 
                            inside = inout(second)
                         end if
                      end if
                   end do
                end if
             end do
          end do
       case default
          stop 'ERROR: more than triply-nested inclusions'
       end select
    else
       inside = 0
       Xcp = 0.0; Ycp = 0.0
       Rcp = 0.0; Pcp = 0.0
    end if

    ! determine geometry to pertinant wells from observation point
    forall (well = 1:nw, CIWellIn(inside,well))
       Xwp(well) = CalcX - WLx(well)
       Ywp(well) = CalcY - WLy(well)
       Rwp(well) = sqrt(Xwp(well)**2 + Ywp(well)**2)
       Pwp(well) = atan2(Ywp(well),Xwp(well))
    end forall

    ! move observation points inside wells to edge of wells
    do well = 1,nw 
       if (Rwp(well) < WLr(well)) then
          Rwp(well) = WLr(well)
       end if
    end do

  end subroutine  calcLocation

 ! free memory shared between velx and vely calcs
  subroutine freecalcmem()
    deallocate(CIdPot_dAng, CIdPot_dR)
  end subroutine freecalcmem

end module calc_routines
