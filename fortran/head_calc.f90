! $Id: head_calc.f90,v 1.1 2005/12/05 19:17:53 kris Exp kris $
module calc_routines
  implicit none

  private
  public :: headCalc, velxCalc, velyCalc

contains

function headCalc(p) result(H)
  use constants, only : DP, PI, CZERO
  use calc_shared_data, only : calcPoint, calcCoeff
  use complex_Bessel, only : cbesk
  use wells, only : wellHead

  implicit none

  complex(DP), dimension(:), intent(in) :: p
  complex(DP), dimension(size(p,1)) :: H
  complex(DP), dimension(INVm*2+1, 0:CIn, CInum) :: a,b,c,d
  complex(DP), dimension(INVm*2+1,0:CInum) :: q
  complex(DP), dimension(INVm*2+1) :: circ, PotInclBg, PotWellBg, PotInclPar, PotInclIn, PotWellIn
  complex(DP), dimension(INVm*2+1:0:CIn) :: beskRcp, beskR0, besiRcp, besiR0
  real(DP), dimension(CInum) :: Xcp, Ycp, Rcp, Pcp
  integer, dimension(CInum) :: inout
  integer :: ni, nw, incl,j,k, first, second, third, inside
  integer :: nz, ierr
  real(DP) :: rk
  
  ! get values from struct
  ni = ICNum nw = WLnum

  Xop = obs(1); 
  Yop = obs(2);

  Xow = Q.Xw; Yow = Q.Yw;    

  ! in this subroutine most everything is a column vector with each row 
  ! being for a different value of p (Laplace parameter).
  Pn = INVm*2+1;

  ! number of terms in sin/bessel funciton expansion 
  N = CIn;

  ! normal setup
  q(1:Pn,0) = sqrt( p(1:Pn)*/av(0) )
  do incl = 1,ni
     q(1:Pn,incl) = sqrt( p(1:Pn)/av(incl) )
  enddo

  ! coefficients passed through module
  a(1:Pn, 0:N, 1:ni) = calcCoeff(1:Pn, 0:N, 1:ni)
  b(1:Pn, 1:N, 1:ni) = calcCoeff(1:Pn, N+1:2*N, 1:ni)  ! b_0 not included
  c(1:Pn, 0:N, 1:ni) = calcCoeff(1:Pn, 2*N+1:3*N+1, 1:ni)
  d(1:Pn, 1:N, 1:ni) = calcCoeff(1:Pn, 3*N+2:4*N+1, 1:ni)  ! d_0 not included

  ! determine if observation point is inside an inclusion
  inout(1:ni) = 0
  k = 0;
  do incl = 1,ni
     ! components of vector from center of inclusion to observation point
     Xcp(incl) = Xop - ICx(incl)
     Ycp(incl) = Yop - ICy(incl)
     Rcp(incl) = sqrt(Xcp(incl)**2 + Ycp(incl)**2)
     Pcp(incl) = atan2(Ycp(incl), Xcp(incl))
     if (Rcp(incl) <= ICr(incl)) then    ! inside or on boundary
        k = k + 1;
        inout(k) = incl;
     endif
  enddo
  select case (k)
  case (0)
     inside = 0
  case (1)
     inside = inout(1)
  case (2)
     if (CIInclUp(inout(1)) == 0) then
        inside = inout(2);
     else
        inside = inout(1);
     end if
  case (3)
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
     stop 'AYE! too many nested inclusions'
  end select

  !##################################################
  !! calculation point is outside all inclusions
  if (inside == 0) then
     PotInclBg(1:Pn) = CZERO
     do incl = 1,ni
        if (CIInclUp(incl) == 0) then  ! if inclusion is in background
           circ(1:Pn) = CZERO
           do j = 1,Pn
              call cbesk(Rcp(incl)*q(j,0), 0.0_DP, 1, N+1, beskRcp(j,0:N), nz, ierr)
              if (nz /= 0 .or. ierr /= 0) write (*,*) 'error1 on return of cbesk in head_calc'
              call cbesk(ICr(incl)*q(j,0), 0.0_DP, 1, N+1, beskR0(j,0:N), nz, ierr)
              if (nz /= 0 .or. ierr /= 0) write (*,*) 'error2 on return of cbesk in head_calc'
           end do
           do  k = 0,N
              rk = real(k,DP)
              circ(1:Pn) = circ(:) + beskRcp(:,k)/beskR0(:,k)* &
                   & (a(:,k,incl)*cos(rk*Pcp(incl)) + &
                   &  b(:,k,incl)*sin(rk*Pcp(incl)))
           end do
           PotInclBg(1:Pn) = PotInclBg(:) + circ(:)
        end if
     end do ! for ni
  
     PotWellBg = CZERO
     do well = 1,nw
        if (CIWellUp(well) == 0) then ! if well is in background
           Rwp = sqrt((Xop - WLx(well))**2 + (Yop - WLy(well))**2)
           if (Rwp < WLr(well)) Rwp = WLr(well)
           PotWellBg(1:Pn) = PotWellBg(:) + WellHead(well,p(1:Pn),q(1:Pn,0),Rwp)
        end if
     end do
     ! result for invlap (head, not discharge potential)
     H(1:Pn) = (PotInclBg(:) + PotWellBg(:))/kv(0) 
     
  !##################################################
  !! calculation point is inside (or on bdry of) an inclusion
  else
     PotInclPar = CZERO; PotInclIn = CZERO;  PotWellIn = CZERO
     if (CIcalcin(inside)) then
  
        do j = 1,Pn
           call cbesi(Rcp(inside)*q(j,inside), 0.0_DP, 1, N+1, besiRcp(j,0:N), nz, ierr)
           if (nz /= 0 .or. ierr /= 0) write (*,*) 'error1 on return of cbesi in head_calc'
           call cbesi(ICr(inside)*q(j,inside), 0.0_DP, 1, N+1, besiR0(j,0:N), nz, ierr)
           if (nz /= 0 .or. ierr /= 0) write (*,*) 'error2 on return of cbesi in head_calc'
        end do

        ! "parent" inclusion, which point is within
        do k = 0,N         ! D_0 = zero
           rk = real(k,DP)
           PotInclPar(1:Pn) = PotInclPar(:) + besiRcp(:,k)/besiR0(:,k)* &
                & ( c(:,k,inside)*cos(rk*Pcp(inside)) + &
                &   d(:,k,inside)*sin(rk*Pcp(inside)) )
        end do
    
        ! constant flux over area of inclusion
        PotInclPar(1:Pn) = PotInclPar(:) + circHead(inside,p(1:Pn))/q(:,inside)**2
    
        ! effects of any inclusions, which may be inside parent too
        do other = 1,ni
           if (CIInclIn(inside,other)) then
              do j = 1,Pn
                 call cbesk(Rcp(other)*q(j,inside), 0.0_DP, 1, N+1, beskRcp(j,0:N), nz, ierr)
                 if (nz /= 0 .or. ierr /= 0) write (*,*) 'error1 on return of cbesk in head_calc'
                 call cbesk(ICr(other)*q(j,inside), 0.0_DP, 1, N+1, beskR0(j,0:N), nz, ierr)
                 if (nz /= 0 .or. ierr /= 0) write (*,*) 'error2 on return of cbesk in head_calc'
              enddo
              circ = CZERO
              do k = 0,N
                 rk = real(k,DP)
                 circ(1:Pn) = circ(:) + beskRcp(:,k)/beskR0(:,k)* &
                      & ( a(:,k,other)*cos(rk*Pcp(other)) + &
                      &   b(:,k,other)*sin(rk*Pcp(other)))
              enddo
              PotInclIn(1:Pn) = PotInclIn(:) + circ(:)
           endif
        enddo
    
        ! effects of wells inside inclusion
        do well = 1,nw
           if (CIWellIn(inside,well)) then
              Rwp = sqrt((Xop - WLx(well))^2 + (Yop - WLy(well))^2);
              if (Rwp < Q.rw(well)) then
                 Rwp = Q.rw(well)
              end if
              PotWellIn(1:Pn,1) = PotWellIn(:,1) + WellHead(inside,p(1:Pn),q(1:Pn,inside),Rwp)
           end if
        end do
     end if  !if I.calcin

     ! result for invlap (head, not discharge potential)
     H(1:Pn) = (PotWellIn(:) + PotInclIn(:) + PotInclPar(:))/kv(inside)
  end if ! inout


end function headCalc
end module calc_routines
