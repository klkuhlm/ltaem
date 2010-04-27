! $Id: calc_ellipse_routines.f90,v 1.9 2007/08/28 23:15:04 kris Exp kris $
module calc_ellipse_routines
  use constants, only : DP
  implicit none

  private
  public :: headCalc, velCalc

contains

  !##################################################
  function headCalc(p,calcX,calcY,calcCoeff,AA,BB,qq) result(H)
    use element_specs, only : EIn,EInum,WLnum,av,kv,sv,EIWellIn,EIInclUp,&
         & EIeta,EICalcIn,EIInclIn,EIArea,EIms
    use constants, only : DP, PI, CZERO, RZERO
    use shared_mathieu, only : Amat => A, Bmat => B,q
    use mathieu_functions
    use wells_el, only : wellHead, ElTimeArea

    complex(DP), dimension(:), intent(in) :: p
    real(DP), intent(in) :: calcX,calcY
    complex(DP), dimension(1:,0:,1:), intent(in) :: calcCoeff
    complex(DP), intent(in), dimension(1:EIms,0:EIms-1,1:2,2*EInum,size(p,1)) :: AA,BB
    complex(DP), intent(in), dimension(size(p,1),1:2*EInum) :: qq

    complex(DP), dimension(size(p,1)) :: H  !! f(p)
    complex(DP), dimension(size(p,1), 0:EIn, EInum) :: a,c
    complex(DP), dimension(size(p,1), 1:EIn, EInum) :: b,d
    complex(DP), dimension(size(p,1),0:EInum) :: kap
    complex(DP), dimension(size(p,1)) :: PotInclPar
    complex(DP), dimension(size(p,1), EInum) :: PotInclBg, PotInclIn
    complex(DP), dimension(size(p,1),WLnum) :: PotWell
    complex(DP), dimension(size(p,1),0:EIn) :: Ke, Ke0, Ko, Ko0
    complex(DP), dimension(size(p,1),0:EIn) :: Ie, Ie0, Io, Io0
    complex(DP), dimension(size(p,1),0:EIn) :: ce
    complex(DP), dimension(size(p,1),1:EIn) :: se
    real(DP), dimension(EInum) :: Eep, Pep
    real(DP), dimension(WLnum) :: Rwp, Pwp
    
    integer :: ni, nw, el, inside, other, well
    integer :: np, N, i, k
    integer, dimension(0:EIn) :: vi

    ni = EINum; nw = WLnum
    np = size(p,1)

    ! number of terms in sin/bessel funciton expansion 
    N = EIn;

    forall (i = 0:N) vi(i) = i
    forall (i = 0:ni) kap(1:np,i) = sqrt(p(1:np)/av(i))

    ! coefficients passed through module
    a(1:np, 0:N, 1:ni) = calcCoeff(1:np, 0:N, 1:ni)
    b(1:np, 1:N, 1:ni) = calcCoeff(1:np, N+1:2*N, 1:ni)      ! b_0 not included
    c(1:np, 0:N, 1:ni) = calcCoeff(1:np, 2*N+1:3*N+1, 1:ni)
    d(1:np, 1:N, 1:ni) = calcCoeff(1:np, 3*N+2:4*N+1, 1:ni)

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(CalcX,CalcY,inside,Eep,Pep,Rwp,Pwp) !Pwp not used here

    ! effects of wells inside or outside inclusion
    PotWell(1:np,1:nw) = CZERO
    do well = 1,nw 
       if (EIWellIn(inside,well)) then
          PotWell(1:np,well) = WellHead(well,p(:),kap(:,inside),Rwp(well))
       end if
    end do
    
!!#ifdef DEBUG    
    write(*,'(I1)',advance='no')inside
!!#endif
    
    !##################################################
    !! calculation point is outside all elements (background == 0)
    if (inside == 0) then
       PotInclBg = CZERO
       do el = 1,ni
          if (EIInclUp(el) == 0) then  ! if inclusion is in background
             
             !! compute external Mathieu functions for all values of p
             do k=1,np
                Amat = AA(:,:,:,ni+el,k); Bmat = BB(:,:,:,ni+el,k)
                q = qq(k,ni+el)
                call mmatKeKo(vi(0:N),Eep(el),   Ke(k,0:N), Ko(k,0:N))
                call mmatKeKo(vi(0:N),EIeta(el),Ke0(k,0:N),Ko0(k,0:N))
                
                ce(k,0:N) = mmatce(vi(0:N),Pep(el))
                se(k,1:N) = mmatse(vi(1:N),Pep(el))
             end do

             PotInclBg(1:np,el) = &
                  & sum(Ke(:,0:N)/Ke0(:,0:N)*a(:,0:N,el)*ce(:,0:N),dim=2) + &
                  & sum(Ko(:,1:N)/Ko0(:,1:N)*b(:,1:N,el)*se(:,1:N),dim=2)
          end if
       end do ! for ni
       
       ! result for invlap (head, not discharge potential)
       !! should have area source term in background too, eventually.
       H(1:np) = (sum(PotInclBg,dim=2) + sum(PotWell,dim=2))/kv(0)
       
    !##################################################
    !! calculation point is inside (or on bdry of) an inclusion
    else
       PotInclPar = CZERO; PotInclIn = CZERO
       if (EIcalcin(inside)) then
          
          !! compute internal Mathieu functions for all values of p
          do k=1,np
             Amat = AA(:,:,:,inside,k); Bmat = BB(:,:,:,inside,k)
             q = qq(k,inside)
             call mmatIeIo(vi(0:N),Eep(inside),   Ie(k,0:N), Io(k,0:N))
             call mmatIeIo(vi(0:N),EIeta(inside),Ie0(k,0:N),Io0(k,0:N))
             
             ce(k,0:N) = mmatce(vi(0:N),Pep(inside))
             se(k,1:N) = mmatse(vi(1:N),Pep(inside))
          end do
          
          ! element immediately surrounding calc point
          PotInclPar(1:np) = &
               & sum(Ie(:,0:N)/Ie0(:,0:N)*c(:,0:N,inside)*ce(:,0:N),dim=2) + &
               & sum(Io(:,1:N)/Io0(:,1:N)*d(:,1:N,inside)*se(:,1:N),dim=2)
          
          ! effects of any active elements that may be inside parent too
          do other = 1,ni
             if (EIInclIn(inside,other)) then

                do k=1,np
                   Amat = AA(:,:,:,ni+other,k); Bmat = BB(:,:,:,ni+other,k)
                   q = qq(k,ni+other)
                   call mmatKeKo(vi(0:N),Eep(other),   Ke(k,0:N), Ko(k,0:N))
                   call mmatKeKo(vi(0:N),EIeta(other),Ke0(k,0:N),Ko0(k,0:N))
                   
                   ce(k,0:N) = mmatce(vi(0:N),Pep(other))
                   se(k,1:N) = mmatse(vi(1:N),Pep(other))
                end do

                PotInclIn(1:np,other) = &
                     & sum(Ke(:,0:N)/Ke0(:,0:N)*a(:,0:N,other)*ce(:,0:N),dim=2) + &
                     & sum(Ko(:,1:N)/Ko0(:,1:N)*b(:,1:N,other)*se(:,1:N),dim=2)
             end if
          end do
       end if  !if EIcalcin
       
       ! result for invlap (head, not discharge potential)
       H(1:np) = (sum(PotWell,dim=2) + sum(PotInclIn,dim=2) + PotInclPar + &
            & EIarea(inside)*ElTimeArea(inside,p(:))* &
            & sv(inside)/kap(:,inside)**2)/kv(inside)
    end if ! inout

  end function headCalc

  !##################################################
  function velCalc(p,calcX,calcY,calcCoeff,AA,BB,qq) result(vel)

    use element_specs, only : EIn,EInum,WLnum,av,porv,EIWellIn,EIInclUp,&
         & EICalcIn,EIInclIn,EIeta,EIf,EIms
    use constants, only : DP, PI, CZERO, CTWO
    use shared_mathieu, only : Amat => A, Bmat => B, q
    use mathieu_functions
    use wells_el, only : wellFlux

    complex(DP), dimension(:), intent(in) :: p
    real(DP), intent(in) :: calcX,calcY
    complex(DP), dimension(1:,0:,1:), intent(in) :: calcCoeff
    complex(DP), intent(in), dimension(1:EIms,0:EIms-1,1:2,2*EInum,size(p,1)) :: AA,BB
    complex(DP), intent(in), dimension(size(p,1),2*EInum) :: qq

    complex(DP), dimension(size(p,1)) :: EIdPot_dPsi, EIdPot_dEta
    complex(DP), dimension(size(p,1),2) :: vel
    complex(DP), dimension(size(p,1),0:EInum) :: kap
    complex(DP), dimension(size(p,1),0:EIn, 1:EInum) :: a,c
    complex(DP), dimension(size(p,1),1:EIn, 1:EInum) :: b,d
    complex(DP), dimension(size(p,1),2) :: FluxInclPar
    complex(DP), dimension(size(p,1),EInum,2) :: FluxInclBg, FluxInclIn
    complex(DP), dimension(size(p,1),WLnum,2) :: FluxWell
    complex(DP), dimension(size(p,1),0:EIn) :: Ke, Ko, Ie, Io, Ke0, Ko0, DKe, &
         & DKo, Ie0, Io0, DIe, DIo
    complex(DP), dimension(size(p,1),0:EIn) :: ce, Dce
    complex(DP), dimension(size(p,1),1:EIn) :: se, Dse
    real(DP), dimension(EInum) :: Eep, Pep
    real(DP), dimension(WLnum) :: Rwp, Pwp

    integer :: ni, nw, inside, other, well, i
    integer :: np, N, el, k
    integer, dimension(0:EIn) :: vi
    real(DP), dimension(EInum) :: hsq

    ni = EINum; nw = WLnum
    np = size(p,1);

    ! number of terms in sin/bessel function expansion 
    N = EIn;

    forall (i = 0:N) vi(i) = i
    forall (i = 0:ni) kap(1:np,i) = sqrt(p(1:np)/av(i))

    ! coefficients passed through module
    a(1:np, 0:N, 1:ni) = calcCoeff(1:np, 0:N, 1:ni)
    b(1:np, 1:N, 1:ni) = calcCoeff(1:np, N+1:2*N, 1:ni)  ! b_0 not included
    c(1:np, 0:N, 1:ni) = calcCoeff(1:np, 2*N+1:3*N+1, 1:ni)
    d(1:np, 1:N, 1:ni) = calcCoeff(1:np, 3*N+2:4*N+1, 1:ni)

!!$    write(*,*) 'calc a:',maxval(abs(a(:,0:N,1)))
!!$    write(*,*) 'calc b:',maxval(abs(b(:,1:N,1)))
!!$    write(*,*) 'calc c:',maxval(abs(c(:,0:N,1)))
!!$    write(*,*) 'calc d:',maxval(abs(d(:,1:N,1)))

    ! determine which inclusion this point is in, and related geometry
    call CalcLocation(CalcX,CalcY,inside,Eep,Pep,Rwp,Pwp)

    !! metric coefficient
    hsq(1:ni) = EIf(1:ni)**2/2.0_DP*(cosh(2.0_DP*Eep(1:ni)) - cos(2.0_DP*Pep(1:ni)))

    ! effects of wells inside or outside inclusions
    FluxWell = CZERO
    do well = 1,nw
       if (EIWellIn(inside,well)) then ! if well is in background
          FluxWell(:,well,1) = cos(Pwp(well))*WellFlux(well,p(:),kap(:,inside),Rwp(well))
          FluxWell(:,well,2) = sin(Pwp(well))*WellFlux(well,p(:),kap(:,inside),Rwp(well))
       end if
    end do

    !##################################################
    !! calculation point is outside all inclusions
    if (inside == 0) then

       ! zero out
       FluxInclBg = CZERO

       do el = 1,ni
          if (EIInclUp(el) == 0) then  ! inclusion is in background

             !! compute external Mathieu functions for all values of p
             do k=1,np
                Amat = AA(:,:,:,ni+el,k); Bmat = BB(:,:,:,ni+el,k)
                q = qq(k,ni+el)
                call mmatDKeDKo(vi(0:N),Eep(el),DKe(k,0:N),DKo(k,0:N))
                call mmatKeKo(vi(0:N),Eep(el),   Ke(k,0:N), Ko(k,0:N)) 
                call mmatKeKo(vi(0:N),EIeta(el),Ke0(k,0:N),Ko0(k,0:N)) 
                
                ce(k,0:N) = mmatce(vi(0:N),Pep(el))
                se(k,1:N) = mmatse(vi(1:N),Pep(el))
                Dce(k,0:N) = mmatDce(vi(0:N),Pep(el))
                Dse(k,1:N) = mmatDse(vi(1:N),Pep(el))
             end do

             EIdPot_dEta(1:np) = &
                  & sum(DKe(:,0:N)/Ke0(:,0:N)*a(:,0:N,el)*ce(:,0:N),dim=2) + &
                  & sum(DKo(:,1:N)/Ko0(:,1:N)*b(:,1:N,el)*se(:,1:N),dim=2)
             EIdPot_dPsi(1:np) = &
                  & sum(Ke(:,0:N)/Ke0(:,0:N)*a(:,0:N,el)*Dce(:,0:N),dim=2) + &
                  & sum(Ko(:,1:N)/Ko0(:,1:N)*b(:,1:N,el)*Dse(:,1:N),dim=2)
             
             ! flux in x-direction due to background inclusions            
             FluxInclBg(1:np,el,1) = EIf(el)/hsq(el)*(&
                  & sinh(Eep(el))*cos(Pep(el))*EIdPot_dEta(:) -&
                  & cosh(Eep(el))*sin(Pep(el))*EIdPot_dPsi(:))
             ! y-flux
             FluxInclBg(1:np,el,2) = EIf(el)/hsq(el)*(&
                  & cosh(Eep(el))*sin(Pep(el))*EIdPot_dEta(:) +&
                  & sinh(Eep(el))*cos(Pep(el))*EIdPot_dPsi(:))
          end if
       end do

#ifdef DEBUG       
       if (Eep(1) < 0.75) then
          write(*,*) 'E:',Eep(1),' P:',Pep(1)
          write(*,*) 'OUT xFIBG:',FluxInclBg(1:np,1,1)
          write(*,*) 'OUT yFIBG:',FluxInclBg(1:np,1,2)
          write(*,*) 'OUT xFW:',sum(FluxWell(1:np,1:nw,1),dim=2)
          write(*,*) 'OUT yFW:',sum(FluxWell(1:np,1:nw,2),dim=2)
       end if
#endif

       ! result for invlap velocity
       vel(1:np,1:2) = (sum(FluxInclBg,dim=2) + sum(FluxWell,dim=2))/porv(0)

    !##################################################
    !! calculation point is inside (or on bdry of) an inclusion
    else

       ! zero out
       FluxInclIn = CZERO;  FluxInclPar = CZERO

#ifdef DEBUG
       write(*,*) 'X',calcX,' Y',calcY,' h2',hsq(inside),' E',Eep(inside),' P',Pep(inside)
#endif

       if (EIcalcin(inside)) then

          !! compute internal Mathieu functions for all values of p
          do k=1,np
             Amat = AA(:,:,:,inside,k); Bmat = BB(:,:,:,inside,k)
             q = qq(k,inside)
             call mmatDIeDIo(vi(0:N),Eep(inside),DIe(k,0:N),DIo(k,0:N))
             call mmatIeIo(vi(0:N),Eep(inside),   Ie(k,0:N), Io(k,0:N))
             call mmatIeIo(vi(0:N),EIeta(inside),Ie0(k,0:N),Io0(k,0:N))
             
             ce(k,0:N) = mmatce(vi(0:N),Pep(inside))
             se(k,1:N) = mmatse(vi(1:N),Pep(inside))
             Dce(k,0:N) = mmatDce(vi(0:N),Pep(inside))
             Dse(k,1:N) = mmatDse(vi(1:N),Pep(inside))
          end do

          EIdPot_dEta(1:np) = &
               & sum(DIe(:,0:N)/Ie0(:,0:N)*c(:,0:N,inside)*ce(:,0:N),dim=2) + &
               & sum(DIo(:,1:N)/Io0(:,1:N)*d(:,1:N,inside)*se(:,1:N),dim=2)
          EIdPot_dPsi(1:np) = &
               & sum(Ie(:,0:N)/Ie0(:,0:N)*c(:,0:N,inside)*Dce(:,0:N),dim=2) + &
               & sum(Io(:,1:N)/Io0(:,1:N)*d(:,1:N,inside)*Dse(:,1:N),dim=2)

          FluxInclPar(1:np,1) = EIf(inside)/hsq(inside)*(&
               & sinh(Eep(inside))*cos(Pep(inside))*EIdPot_dEta(:) -&
               & cosh(Eep(inside))*sin(Pep(inside))*EIdPot_dPsi(:))

          FluxInclPar(1:np,2) = EIf(inside)/hsq(inside)*(&
               & cosh(Eep(inside))*sin(Pep(inside))*EIdPot_dEta(:) +&
               & sinh(Eep(inside))*cos(Pep(inside))*EIdPot_dPsi(:))

          ! effects of any inclusions, which may be inside parent too
          do other = 1,ni
             if (EIInclIn(inside,other)) then
                
                !! compute external Mathieu functions for all values of p
                do k=1,np
                   Amat = AA(:,:,:,ni+other,k); Bmat = BB(:,:,:,ni+other,k)
                   q = qq(k,ni+other)
                   call mmatDKeDKo(vi(0:N),Eep(other),DKe(k,0:N),DKo(k,0:N))
                   call mmatKeKo(vi(0:N),Eep(other),   Ke(k,0:N), Ko(k,0:N)) 
                   call mmatKeKo(vi(0:N),EIeta(other),Ke0(k,0:N),Ko0(k,0:N)) 
                   
                   ce(k,0:N) = mmatce(vi(0:N),Pep(other))
                   se(k,1:N) = mmatse(vi(1:N),Pep(other))
                   Dce(k,0:N) = mmatDce(vi(0:N),Pep(other))
                   Dse(k,1:N) = mmatDse(vi(1:N),Pep(other))
                end do

                EIdPot_dEta(1:np) = &
                     & sum(DKe(:,0:N)/Ke0(:,0:N)*a(:,0:N,other)*ce(:,0:N),dim=2) + &
                     & sum(DKo(:,1:N)/Ko0(:,1:N)*b(:,1:N,other)*se(:,1:N),dim=2)
                EIdPot_dPsi(1:np) = &
                     & sum(Ke(:,0:N)/Ke0(:,0:N)*a(:,0:N,other)*Dce(:,0:N),dim=2) + &
                     & sum(Ko(:,1:N)/Ko0(:,1:N)*b(:,1:N,other)*Dse(:,1:N),dim=2)
                
                FluxInclIn(1:np,other,1) = EIf(other)/hsq(other)*(&
                     & sinh(Eep(other))*cos(Pep(other))*EIdPot_dEta(:) -&
                     & cosh(Eep(other))*sin(Pep(other))*EIdPot_dPsi(:))
                FluxInclIn(1:np,other,2) = EIf(other)/hsq(other)*(&
                     & cosh(Eep(other))*sin(Pep(other))*EIdPot_dEta(:) + &
                     & sinh(Eep(other))*cos(Pep(other))*EIdPot_dPsi(:))
             end if
          end do

       end if  !if EIcalcin

#ifdef DEBUG
       write(*,*) 'IN dEta:',EidPot_dEta(1:np)
       write(*,*) 'IN dPsi:',EIdPot_dPsi(1:np)
       write(*,*) 'IN xFIP:',FluxInclPar(1:np,1)
       write(*,*) 'IN yFIP:',FluxInclPar(1:np,2)
#endif

       ! result for invlap
       vel(1:np,1:2) = (sum(FluxWell(1:np,1:nw,1:2),dim=2) + &
               & sum(FluxInclIn(1:np,1:ni,1:2),dim=2) + FluxInclPar(1:np,1:2))/porv(inside)
       
    end if
  end function velCalc

  !##################################################
  subroutine calcLocation(CalcX,CalcY,inside,Eep,Pep,Rwp,Pwp)
    use constants, only : DP, EYE
    use element_specs, only : EInum, EIx, EIy, EItheta, EIeta, EIf, &
         & EIInclUp, EIInclIn, WLnum, WLx, WLy, WLr, EIWellIn

    real(DP), intent(in) :: CalcX, CalcY
    integer, intent(out) :: inside
    real(DP), dimension(1:EInum), intent(out) :: Eep, Pep
    real(DP), dimension(1:WLnum), intent(out) :: Rwp, Pwp
    real(DP), dimension(1:WLnum) :: Xwp, Ywp
    integer :: incl, well, first, second, third, ni, nw, k
    integer, dimension(EInum) :: inout
    complex(DP), dimension(1:EInum) :: calcW, Z0
    complex(DP) :: calcZ

    ni = EInum; nw = WLnum

    ! determine if observation point is inside an inclusion
    inout(1:ni) = 0
    k = 0;

    Z0(1:EInum) = cmplx(EIx(1:EInum),EIy(1:EInum),DP)  ! centers of ellipses
    calcZ = cmplx(CalcX,CalcY,DP)

    calcW(1:EInum) = acosh(exp(-EYE*EItheta(:))*(calcZ - Z0(:))/EIf(:))
    Eep(1:Einum) = real(calcW(1:EInum))
    Pep(1:Einum) = aimag(calcW(1:EInum))

    do incl = 1,ni
       if (Eep(incl) <= EIeta(incl)) then    ! inside or on boundary
          k = k + 1;
          inout(k) = incl;
       end if
    end do

    select case (k)
    case (0) !background
       inside = 0
    case (1) !one inclusion
       inside = inout(1)
    case (2) !two inclusions 
       if (EIInclUp(inout(1)) == 0) then
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
                      if (EIInclIn(inout(first),inout(second))  .and. &
                           & EIInclIn(inout(third),inout(first))) then 
                         inside = inout(second)
                      end if
                   end if
                end do
             end if
          end do
       end do
    case default
       stop 'ERROR: more than triply-nested inclusions, &
            & general case not implemente yet'
    end select

    ! determine geometry to pertinant wells from observation point
    forall (well = 1:nw, EIWellIn(inside,well))
       Xwp(well) = CalcX - WLx(well)
       Ywp(well) = CalcY - WLy(well)
       Rwp(well) = abs(cmplx(Xwp(well),Ywp(well),DP))
       Pwp(well) = atan2(Ywp(well),Xwp(well))
    end forall

    ! move observation points inside wells to edge of wells
    do well = 1,nw
       if(Rwp(well) < WLr(well)) Rwp(well) = WLr(well)
    end do

  end subroutine  calcLocation

   ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  elemental function acosh(z) result(f)
    use constants, only : DP, RONE

    complex(DP), intent(in) :: z
    complex(DP) :: f

    if(real(z) > 0.0) then
       f = log(z + sqrt(z**2 - RONE))
    else
       f = log(z - sqrt(z**2 - RONE))
    end if

  end function acosh

end module calc_ellipse_routines
