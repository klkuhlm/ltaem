! $Id: matching.f90,v 1.8 2008/12/10 02:47:22 kris Exp kris $

module matching_old
  implicit none

  private 
  public :: CircInverse_old, MatchIter !!, freeMatchMem

contains

  !##################################################
  ! calculate inverse exploiting the  orthogonallity of the columns of Am
  ! this subroutine called once for each value of p before beginning iteration
  subroutine CircInverse_old(p,Gm)
    use constants, only: DP,PI,TWOPI,CZERO,CONE,CTWO,RONE,RZERO
    use element_specs, only : CIn,CIm,CInum,kv,CImatch,CIInclUp, CIibnd,CIr,CICalcIn
    use bessel_functions ! my wrapper for Amos complex Bessel fcn subroutine
    use shared_matching_data, only : CIPcm
    use matrix_inverse
    use leaky_q

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! rows of Gm are orthonormal eigenvectors, multiplication with b gives LS-solution
    complex(DP), dimension(1:4*CIn+2, 1:2*CIm, 1:CInum), intent(out) :: Gm

    ! INTERNAL VARIABLES
      ! norm (length squared) of columns of Am
    real(DP), dimension(1:2*Cin+1) :: Diag
      ! orthogonal eigenvectors from separation of variables
    complex(DP), dimension(1:2*CIm, 1:4*CIn+2, 1:CInum) :: Am
      ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
      ! vectors of Bessel functions to minimize calls to Bessel fcn subroutine
    complex(DP), dimension(-1:CIn+1) :: bessk, bessi
      ! cos/sin(theta) -- where theta is spaced evenly along circumference
    real(DP), dimension(1:CIm,0:CIn) :: cosnCIPcm, sinnCIPcm
      ! counters
    integer :: ni, inc, par, k, M, N
    real(DP), dimension(0:CIn) :: rk

    M = CIm; N = CIn; ni = CInum;

    q(0:ni) = compute_CIleaky_q(p)

    sinnCIPcm(:,0) = RZERO;
    Am = CZERO; Gm = CZERO
    
    rk(0:N) = real((/ (k, k=0,N) /),DP)
    

    !******************** Am matrix setup ********************
    ! calcualte coefficient matrix (Am) which depends on the geometry 
    ! of the problem and p (calcualte once for each inclusion)

    do inc = 1,ni
       if (CImatch(inc)) then
          par = CIInclUp(inc)

          select case (CIibnd(inc))
          !^^^^^^^^^^^^^^^^^^^^ head & flux matching ^^^^^^^^^^^^^^^^^^^^
          case (0)

             ! pre-compute bessels (during matching all radial distances are Rcm)
             bessk(0:N+1) = besselk(CIr(inc)*q(par),0,N+2)
             bessi(0:N+1) = besseli(CIr(inc)*q(inc),0,N+2)
             bessk(-1) = bessk(1); bessi(-1) = bessi(1) !since BF sub only for index >= 0

             cosnCIPcm(1:M,0:N) = cos(outer_prod(CIPcm(1:M),rk(0:N)))
             sinnCIPcm(1:M,1:N) = sin(outer_prod(CIPcm(1:M),rk(1:N)))

             ! head (PHI/k) matching
             Am(1:M, 1:N+1, inc)    =   -cosnCIPcm(1:M,0:N)/kv(par) ! a_n
             Am(1:M, N+2:2*N+1, inc)  = -sinnCIPcm(1:M,1:N)/kv(par) ! b_n
             Am(1:M, 2*N+2:3*N+2, inc) = cosnCIPcm(1:M,0:N)/kv(inc) ! c_n
             Am(1:M, 3*N+3:4*N+2, inc) = sinnCIPcm(1:M,1:N)/kv(inc) ! d_n

                ! flux (dPHI/dr) matching
             Am(M+1:2*M, 1:N+1, inc) = q(par)*&
                  & spread((bessk(-1:N-1)+bessk(1:N+1))/(CTWO*bessk(0:N)),1,M)*&
                  & cosnCIPcm(1:M,0:N) ! a_n
             Am(M+1:2*M, N+2:2*N+1, inc) = q(par)*&
                  & spread((bessk(0 :N-1)+bessk(2:N+1))/(CTWO*bessk(1:N)),1,M)*&
                  & sinnCIPcm(1:M,1:N) ! b_n
             
             Am(M+1:2*M, 2*N+2:3*N+2, inc) = q(inc)*&
                  & spread((bessi(-1:N-1)+bessi(1:N+1))/(CTWO*bessi(0:N)),1,M)*&
                  & cosnCIPcm(1:M,0:N) ! c_n
             Am(M+1:2*M, 3*N+3:4*N+2, inc) = q(inc)*&
                  & spread((bessi(0 :N-1)+bessi(2:N+1))/(CTWO*bessi(1:N)),1,M)*&
                  & sinnCIPcm(1:M,1:N) ! d_n

             ! matching uses the general LU-decomposition matrix inversion method
             Gm(1:4*N+2,1:2*M,inc) = matmulc2(inverse(matmulc1(&
                  & Am(1:2*M,1:4*N+2,inc),Am(1:2*M,1:4*N+2,inc))), &
                  & Am(1:2*M,1:4*N+2,inc))

          !^^^^^^^^^^^^^^^^^ specified head ^^^^^^^^^^^^^^^^^^^^
          case (-1)
!!$             print *, 'inv spec head'

             ! head (PHI/k) matching
             Am(1:M, 1:N+1, inc)     = cos(outer_prod(CIPcm(1:M),rk(0:N))) ! a_n
             Am(1:M, N+2:2*N+1, inc) = sin(outer_prod(CIPcm(1:M),rk(1:N))) ! b_n
             
             ! due to orthogonallity of columns, ATA is just a diagonal
             do k = 1,2*N+1
                Diag(k) = length_sq(Am(1:M,k,inc))
             end do

             ! inverse of diagonal is trivial; multiplication w/ AT is simple too
             Gm(1:2*N+1,1:M,inc) = spread(RONE/Diag(1:2*N+1),dim=2,ncopies=M)*&
                  & transpose(conjg(Am(1:M,1:2*N+1,inc)))

             !>>>>> calculate inside the specified head inclusion <<<<<
             if (CIcalcin(inc)) then
                print *, 'inside inv spec head'

                ! head (PHI/k) matching
                Am(1:M, 2*N+2:3*N+2, inc) = cos(outer_prod(CIPcm(1:M),rk(0:N))) ! c_n
                Am(1:M, 3*N+3:4*N+2, inc) = sin(outer_prod(CIPcm(1:M),rk(1:N))) ! d_n

                do k = 1,2*N+1
                   Diag(k) = length_sq(Am(1:M,k+2*N+1,inc))
                end do

                Gm(2*N+2:4*N+2,1:M,inc) = spread(RONE/Diag(1:2*N+1),dim=2,ncopies=M)*&
                     & transpose(conjg(Am(1:M,2*N+2:4*N+2,inc)))

             end if ! if calcin

          !^^^^^^^^^^^^^^^^^^^^ specified flux ^^^^^^^^^^^^^^^^^^^^
          case (+1)
!!$             print *, 'inv spec flux'

             bessk(0:N+1) = besselk(CIr(inc)*q(par),0,N+2)
             bessk(-1) = bessk(1)

             Am(M+1:2*M, 1:N+1, inc) =     -q(par)*&
                  &spread((bessk(-1:N-1) + bessk(1:N+1))/(CTWO*bessk(0:N)),1,M)*&
                  & cos(outer_prod(CIPcm(1:M),rk(0:N))) ! a_n

             Am(M+1:2*M, N+2:2*N+1, inc) = -q(par)*&
                  &spread((bessk(0:N-1) + bessk(2:N+1))/(CTWO*bessk(1:N)),1,M)*&
                  & sin(outer_prod(CIPcm(1:M),rk(1:N))) ! b_n
             
             do k = 1,2*N+1
                Diag(k) = length_sq(Am(M+1:2*M,k,inc))
             end do

             Gm(1:2*N+1,M+1:2*M,inc) = spread(RONE/Diag(1:2*N+1),dim=2,ncopies=M)* &
                  & transpose(conjg(Am(M+1:2*M,1:2*N+1,inc)))

             !>>>>> calculate inside the specified flux inclusion <<<<<
             if (CIcalcin(inc)) then
!!$                print *, 'inside inv spec flux'

                bessi(0:N+1) = besseli(CIr(inc)*q(inc),0,N+2)
                bessi(-1) = bessi(1)

                Am(M+1:2*M, 2*N+2:3*N+2, inc) = q(inc)*&
                     & spread((bessi(-1:N-1) + bessi(1:N+1))/(CTWO*bessi(0:N)),1,M)*&
                     & cos(outer_prod(CIPcm(1:M),rk(0:N))) ! c_n
                Am(M+1:2*M, 3*N+3:4*N+2, inc) = q(inc)*&
                     & spread((bessi(0:N-1) + bessi(2:N+1))/(CTWO*bessi(1:N)),1,M)*&
                     & sin(outer_prod(CIPcm(1:M),rk(1:N))) ! d_n

                do k = 1,2*N+1
                   Diag(k) = length_sq(Am(M+1:2*M,k+2*N+1,inc))
                end do

                Gm(2*N+2:4*N+2,M+1:2*M,inc) = spread(RONE/Diag(1:2*N+1),dim=2,ncopies=M)*&
                     & transpose(conjg(Am(M+1:2*M,2+N+2:4*N+2,inc)))

             end if
          end select  ! CIibnd
       end if
    end do
  end subroutine CircInverse_old
  
  !##################################################
  ! calcuale effects of wells (passive elements)
  ! this subroutine called once for each value of p before beginning iteration
  subroutine wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    use shared_matching_data, only : CIPcm, CIRwm, CIPwm
    use constants, only: DP, CZERO
    use element_specs, only : CIm,CInum,WLnum,CImatch,CIInclUp,CIWellIn,CIWellBg
    use wells ! generalized well-effect functions
    use leaky_q

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! potential and flux due to wells inside and in background of each circular element
    complex(DP), dimension(CIm,CInum), intent(out) :: PotWellIn,FluxWellIn,&
         &PotWellBg,FluxWellBg

    ! INTERNAL VARIABLES
      ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(CIm) :: dPot_dRw, dPot_dX, dPot_dY
    integer :: M, well, inc, par, ni, nw

    M = CIm; ni = CInum; nw = WLnum;
    q = compute_CIleaky_q(p)

    PotWellIn = CZERO; FluxWellIn = CZERO;
    PotWellBg = CZERO; FluxWellBg = CZERO;

    !******************** well effects (head and flux) ********************
    do inc = 1,ni
       if (CImatch(inc)) then
          par = CIInclUp(inc)

          do well=1,nw
             if(CIWellIn(inc,well)) then

                !^^^^^^^^^^^^^^^^^^^^ well inside this inclusion ^^^^^^^^^^^^^^^^^^^^
                ! potential due to wells
                PotWellIn(1:M,inc) = PotWellIn(1:M,inc) + &
                     & wellHead(well,p,q(inc),CIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,q(inc),CIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(CIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local polar coords and add to total
                FluxWellIn(1:M,inc) = FluxWellIn(1:M,inc) + &
                     & dPot_dX(1:M)*cos(CIPcm(1:M)) + dPot_dY(1:M)*sin(CIPcm(1:M))

             elseif(CIWellBg(inc,well)) then

                !^^^^^^^^^^^^^^^^^^^^ well in bg of this inclusion ^^^^^^^^^^^^^^^^^^^^
                ! potential due to wells
                PotWellBg(1:M,inc) = PotWellBg(1:M,inc) + &
                     & wellHead(well,p,q(par),CIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,q(par),CIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(CIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local polar coords and add to total
                FluxWellBg(1:M,inc) = FluxWellBg(1:M,inc) + &
                     & dPot_dX(1:M)*cos(CIPcm(1:M)) + dPot_dY(1:M)*sin(CIPcm(1:M))
             end if
          end do
       end if
    end do
  end subroutine wellEffects


!!$  !##################################################
!!$  ! calcuale effects of all elements (passive and active elements)
!!$  ! on the wells where wellbore storage must be accounted for
!!$
!!$  subroutine OnStor(p,a,b,c,d,PotOnStor,FluxOnStor)
!!$    use constants, only: DP, CZERO, CTWO
!!$    use element_specs, only : CIwellup, CIWellBg, CIWellIn, WLstor, WLnum, CInum, &
!!$         & WLx, WLy, CIn, CIx, CIy, CIr, CIInclIn
!!$    use wells 
!!$    use bessel_functions
!!$    use leaky_q
!!$
!!$    ! ARGUMENTS EXTERNAL TO SUBROUTINE
!!$      ! scalar Laplace parameter
!!$    complex(DP), intent(in) :: p
!!$    complex(DP), intent(in), dimension(0:CIn,CInum) :: a,b,c,d 
!!$      ! potential and flux due to wells inside and in background of each circular element
!!$    complex(DP), dimension(WLnum), intent(out) :: PotOnStor,FluxOnStor
!!$
!!$    ! INTERNAL VARIABLES
!!$      ! sqrt(p/alpha)
!!$    complex(DP), dimension(0:CInum) :: q
!!$    complex(DP), dimension(WLnum) :: PotWellOnStor, FluxWellOnStor
!!$    complex(DP), dimension(WLnum) :: PotInclBg, FluxInclBg
!!$    complex(DP) :: dPot_dR, dPot_dAng
!!$    
!!$    integer :: par, nw, ni, targW, srcW, src, w, N
!!$    real(DP), dimension(WLnum,WLnum) :: WLXw,WLYw,WLRw
!!$    real(DP), dimension(CInum,WLnum) :: CIXcp,CIYcp,CIRcp,CIPcp
!!$    Real(DP), dimension(0:CIn) :: rk
!!$
!!$    complex(DP), dimension(0:CIn) :: beskR0,besiR0
!!$    complex(DP), dimension(-1:CIn+1) :: beskRcp,besiRcp
!!$    real(DP), dimension(0:CIn) :: cosPC, sinPC
!!$
!!$    nw = WLnum; ni = CInum; N = CIn
!!$    q = compute_CIleaky_q(p)
!!$
!!$    PotWellonStor = CZERO; FluxWellonStor = CZERO;
!!$
!!$    !! well-on-stor portion
!!$    !! ((((((((((((((((((((((((((((((((((((((((((((((((((
!!$
!!$    ! x,y,r & theta from well1 to well2
!!$    forall (targW = 1:nw, srcW = 1:nw, CIWellUp(targW) == CIWellUp(srcW))
!!$       WLXw(targW,srcW) = WLx(targW) - WLx(srcW)
!!$       WLYw(targW,srcW) = WLy(targW) - WLy(srcW)
!!$       WLRw(targW,srcW) = sqrt(WLXw(targW,srcW)**2 + WLYw(targW,srcW)**2)
!!$       ! angle not needed
!!$    end forall
!!$
!!$    do targW = 1,nw
!!$       if (WLstor(targW)) then
!!$          par = CIWellUp(targW)  
!!$
!!$          do srcW = 1,nw
!!$             if(targW /= srcW .and. CIWellUp(srcW) == par) then
!!$                ! potential 
!!$                
!!$                print *, srcW,p,q(par),[WLRw(targW,srcW)]
!!$
!!$                PotWellonStor(targW) = PotWellonStor(targW) + &
!!$                     & sum(wellHead(srcW,p,q(par),[WLRw(targW,srcW)]))
!!$
!!$                ! radial flux
!!$                dPot_dR = sum(wellFlux(srcW,p,q(par),[WLRw(targW,srcW)]))
!!$
!!$                ! project
!!$                FluxWellOnStor(targW) = FluxWellOnStor(targW) + &
!!$                     &  dPot_dR*(WLXw(targW,srcW)*WLx(targW) + &
!!$                     &           WLYw(targW,srcW)*WLy(targW))               
!!$
!!$             end if
!!$          end do
!!$       end if
!!$    end do
!!$
!!$    !! circ-on-stor portion
!!$    !! ((((((((((((((((((((((((((((((((((((((((((((((((((
!!$
!!$    forall (w = 0:N)  rk(w) = real(w,DP)
!!$    PotInclBg = CZERO;  FluxInclBg = CZERO
!!$
!!$    forall(src = 1:ni, w = 1:nw, CIWellIn(src,w) .or. CIWellBg(src,w))
!!$       ! components of vector from center of inclusion to observation point
!!$       CIXcp(src,w) = WLx(w) - CIx(src)
!!$       CIYcp(src,w) = WLy(w) - CIy(src)
!!$       CIRcp(src,w) = sqrt(CIXcp(src,w)**2 + CIYcp(src,w)**2)
!!$       CIPcp(src,w) = atan2(CIYcp(src,w), CIXcp(src,w))
!!$    end forall
!!$
!!$    write(*,*) 'a-d',maxval(abs(a)),maxval(abs(b)),maxval(abs(c)), maxval(abs(d))
!!$
!!$    ! loop over wellbore storage wells
!!$    do w = 1,nw
!!$       if (WLstor(w)) then
!!$          par = CIWellUp(w)
!!$
!!$          ! this storage well inside a circle; 0 is "background"
!!$          if (par > 0) then
!!$             print *, 'd01',par,w
!!$             cosPC(0:N) = cos(rk(0:N)*CIPcp(par,w))
!!$             sinPC(0:N) = sin(rk(0:N)*CIPcp(par,w))
!!$
!!$             besiR0(0:N) =    besseli(CIr(par)*q(par),0,N+1)
!!$             besiRcp(0:N+1) = besseli(CIRcp(par,w)*q(par),0,N+2)
!!$             besiRcp(-1) = besiRcp(+1)
!!$
!!$             ! parent potential 
!!$             PotInclBg(w) = PotInclBg(w) + sum(besiRcp(0:N)/besiR0(0:N)* &
!!$                  & ( c(0:N,par)*cosPC(0:N) + d(0:N,par)*sinPC(0:N) ))
!!$
!!$             ! parent flux
!!$             dPot_dR = sum(q(par)*(besiRcp(-1:N-1) + &
!!$                  & besiRcp(1:N+1))/(CTWO*besiR0(0:N))* &
!!$                  ( c(0:N,par)*cosPC(0:N) + d(0:N,par)*sinPC(0:N) ))
!!$
!!$             dPot_dAng = sum(rk(0:N)*besiRcp(0:N)/besiR0(0:N)* &
!!$                  ( d(0:N,par)*cosPC(0:N) - c(0:N,par)*sinPC(0:N) ))
!!$
!!$             ! project flux
!!$             FluxInclBg(w) = FluxInclBg(w) + dPot_dR*(CIXcp(par,w)*WLx(w) + &
!!$                  & CIYcp(par,w)*WLy(w)) + dPot_dAng/CIRcp(par,w)*&
!!$                  & (CIXcp(par,w)*WLy(w) - CIYcp(par,w)*WLx(w))
!!$          end if
!!$
!!$          ! effects of any inclusions which may be in same parent
!!$          ! (including special case of background)
!!$          do src = 1,ni
!!$             if (CIInclIn(par,src)) then
!!$                print *, 'd02',par,w,src
!!$                cosPC(0:N) = cos(rk(0:N)*CIPcp(src,w))
!!$                sinPC(0:N) = sin(rk(0:N)*CIPcp(src,w))
!!$
!!$                beskR0(0:N) =    besselk(CIr(src)*q(par),0,N+1)
!!$                beskRcp(0:N+1) = besselk(CIRcp(src,w)*q(par),0,N+2)
!!$                beskRcp(-1) = beskRcp(+1)  
!!$
!!$                ! potential
!!$                PotInclBg(w) = PotInclBg(w) + sum(beskRcp(0:N)/beskR0(0:N)* &
!!$                     & ( a(0:N,src)*cosPC(0:N) + b(0:N,src)*sinPC(0:N) ))
!!$
!!$                ! flux
!!$                dPot_dR =  -sum(q(par)*(beskRcp(-1:N-1) + &
!!$                     & beskRcp(1:N+1))/(CTWO*beskR0(0:N))* &
!!$                     ( a(0:N,src)*cosPC(0:N) + b(0:N,src)*sinPC(0:N) ))
!!$
!!$                dPot_dAng = sum(rk(0:N)*beskRcp(0:N)/beskR0(0:N)* &
!!$                     ( b(0:N,src)*cosPC(0:N) - a(0:N,src)*sinPC(0:N) ))
!!$
!!$                ! project flux
!!$                FluxInclBg(w) = FluxInclBg(w) + &
!!$                     & dPot_dR*(CIXcp(src,w)*WLx(w) + &
!!$                     & CIYcp(src,w)*WLy(w)) + dPot_dAng/CIRcp(src,w)*&
!!$                     & (CIXcp(src,w)*WLy(w) - CIYcp(src,w)*WLx(w))
!!$             end if
!!$          end do
!!$       end if
!!$    end do
!!$    PotOnStor(1:nw) =  PotWellOnStor +  PotInclBg 
!!$    FluxOnStor(1:nw) = FluxWellOnStor  + FluxInclBg
!!$
!!$  end subroutine OnStor

  !##############################################################################
  ! this iterative subroutine calculates the generalized Fourier coefficients
  ! for the active circular elements which depend non-linearly on each other
  subroutine matchIter(Gm,coeff,p,matchtol)
    use shared_matching_data, only : CIPcm,CIRgm,CIPgm
    use constants, only : DP, TWOPI, CZERO, RZERO, RONE, PI
    use circ_passive_elements, only : CircFluxElementEffects, CircHeadElementEffects
    use wells
    use element_specs, only : Cin,CInum,CIm, CIMatchOmega, kv,sv,WLnum, &
         & CIr,CIInclUp,CIInclBg,CIInclIn,CIInclIn,CIibnd,CIarea,CICalcIn,CIspec, &
         & WLstor, WLnum, WLstorCoeff !!, CIWellUp, WLr, WLdskin
    use bessel_functions
    use leaky_q

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! generalized Fourier coefficients (previous value passed in as first guess)
    complex(DP), dimension(0:4*CIn+1,CInum), intent(inout) :: coeff
      ! from CircInverse(); "weight" matrix in LS-solution
    complex(DP), dimension(1:4*CIn+2, 1:2*CIm, 1:CInum), intent(in) :: Gm
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! iteration tolerance for SOR/Gauss-Seidel iteration
    real(DP), intent(in) :: matchtol

    ! INTERNAL VARIABLES
      ! "data" vector in LS-solution (effects of other elements on current)
    complex(DP), dimension(2*CIm, CInum) :: bm
      ! "solution" vector in LS-solution (current and previous iteration)
    complex(DP), dimension(4*CIn+2) :: xm, xm1
      ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
      ! vectors of Bessel functions to minimize calls to Bessel fcn subroutine
    complex(DP), dimension(CIm,0:CIn,CInum) :: besskRcm, bessiRcm
    complex(DP), dimension(CIm,-1:CIn+1,CInum,CInum) :: besskRgm, bessiRgm
      ! same data as x, but in common name (previous iteration for convergence test)
    complex(DP), dimension(0:CIn,CInum) :: a,a1,b,b1,c,c1,d,d1
    complex(DP), dimension(CIm,0:CIn,CInum) :: amv,bmv,cmv,dmv
      ! results from subroutines for wells and passive circular elements
    complex(DP), dimension(CIm,CInum) :: PotWellIn,FluxWellIn,PotWellBg,FluxWellBg
!!$    complex(DP), dimension(Wlnum) :: PotonStor,FluxOnStor
    complex(DP), dimension(CIm,CInum) :: PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(CIm) ::  dPot_dX, dPot_dY
    complex(DP), dimension(CIm,0:CIn) :: dPot_dRg, dPot_dPg
      ! potential due to other circular elements (extra dimension for vectorizing)
    complex(DP), dimension(CIm,0:CIn) :: PotInclBg, PotInclIn, PotInclPar
      ! flux due to other circular elements
    complex(DP), dimension(CIm) :: PotInclParFlux, FluxInclPar, FluxInclBg, FluxInclIn
      ! vector of trig functions to minimize calls to built-in functions
    complex(DP), dimension(CIm,0:CIn) :: cosnPgm, sinnPgm
      ! SOR parameter (imag part = 0)
    complex(DP) :: om
!!$    complex(DP), dimension(WLnum) :: WLStorCoeff1
    integer :: N, M, inc, oth, ni, iteration, par
    integer :: k, loM,hiM,loN,hiN
    real(DP) :: residual
!!$    real(DP), dimension(WLnum) :: rc
    real(DP), dimension(0:CIn) :: rk

    N = CIn; M = CIm; ni = CInum
    om = cmplx(CImatchOmega,0.0_DP,DP)

    q(0:ni) = compute_CIleaky_q(p)


    ! initialize
    besskRgm = CZERO; bessiRgm = CZERO
    PotCElmIn = CZERO; FluxCElmIn = CZERO; PotCElmBg = CZERO; FluxCElmBg = CZERO
    b(0,:) = CZERO; d(0,:) = CZERO; xm1 = CZERO
    iteration = 0; residual = huge(1.0)
    
    ! first guess from previous direct solution
     a(0:N,1:ni) = coeff(0:N,1:ni)      ! are these actually needed anymore?
     b(1:N,1:ni) = coeff(N+1:2*N,1:ni)  ! eventually re-write to only use _mv versions
     c(0:N,1:ni) = coeff(2*N+1:3*N+1,1:ni)
     d(1:N,1:ni) = coeff(3*N+2:4*N+1,1:ni)

!!$    a(0:N,1:ni) = CZERO; b(1:N,1:ni) = CZERO; ! temporarily disabling old guesses
!!$    c(0:N,1:ni) = CZERO; d(1:N,1:ni) = CZERO;

     amv(1:M,0:N,1:ni) = spread(a(0:N,1:ni),dim=1,ncopies=M)
     bmv(1:M,0:N,1:ni) = spread(b(0:N,1:ni),dim=1,ncopies=M)
     cmv(1:M,0:N,1:ni) = spread(c(0:N,1:ni),dim=1,ncopies=M)
     dmv(1:M,0:N,1:ni) = spread(d(0:N,1:ni),dim=1,ncopies=M)

!!$     if (any(WLStor)) then
!!$        WLStorCoeff = CZERO
!!$        WLStorCoeff1 = CZERO
!!$     end if

    rk(0:N) = real((/ (k, k=0,N) /),DP)

    if (WLnum > 0) then
       call wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    end if

    ! Calculate non-matching circular element effects and save their coefficients
    if (any(CIibnd(1:ni) == -2)) then
       call circHeadElementEffects(p,PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg,coeff)
    end if
    if (any(CIibnd(1:ni) == +2)) then
       call circFluxElementEffects(p,PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg,coeff)
    end if

    !******************** Successive Over-Relaxation iteration loop ********************
    SOR: do
!!$       print *, 'iter'
       if (residual <= matchtol ) then
          write (*,'(1x,i4,a6,ES10.4)') iteration,' norm=',residual
          exit SOR
       end if
  
       if (iteration > 0 .and. mod(iteration,5)==0) &
            & write(*,'(1x,i4,a6,ES10.4)') iteration,' norm=',residual

       iteration = iteration + 1
       a1 = a; b1 = b; c1 = c; d1 = d
      
       THIS: do inc = 1,ni

!!$          print *, 'this',inc
          par = CIInclUp(inc) ! parent of this or any other bg inclusions
          PotInclBg = CZERO;  FluxInclBg = CZERO
          PotInclIn = CZERO;  FluxInclIn = CZERO
    
          !^^^^^^^^^^^^^^^^^^^^ loop over other inclusions ^^^^^^^^^^^^^^^^^^^^
          OTHER: do oth = 1,ni

!!$             print *, 'other',oth  
             ! @@@@@ other in background of current @@@@@
             if (CIInclBg(oth,inc)) then
                
                ! pre-compute Bessel-K
                if (iteration == 1) then 
                   besskRcm(1:M,0:N,oth) = spread(besselk(CIr(oth)*q(par),0,N+1),1,M)
                   besskRgm(1:M,0:N+1,oth,inc) = besselk(CIRgm(1:M,oth,inc)*q(par),0,N+2)
                   besskRgm(1:M,-1,oth,inc) = besskRgm(1:M,+1,oth,inc)
                end if

                ! pre-compute sin/cos
                cosnPgm(1:M,1:N) = cos(outer_prod(CIPgm(1:M,oth,inc),rk(1:N)))
                cosnPgm(1:M,0) = RONE
                sinnPgm(1:M,1:N) = sin(outer_prod(CIPgm(1:M,oth,inc),rk(1:N)))
                sinnPgm(1:M,0) = RZERO

                !##########################################
                ! head effects on outside of this inclusion
                !##########################################
                if (CIibnd(inc) /= +1) then ! not specified flux 
                   ! potential effects of other inlcusions
                   ! accumulate effects over all other inclusions
                   PotInclBg(1:M,0:N) = PotInclBg(1:M,0:N) + &
                        & besskRgm(1:M,0:N,oth,inc)/besskRcm(1:M,0:N,oth)* &
                        & (amv(1:M,0:N,oth)*cosnPgm(1:M,0:N) + bmv(1:M,0:N,oth)*sinnPgm(1:M,0:N))
                end if

                !##########################################
                ! flux effects on outside of this inclusion
                !##########################################
                if (CIibnd(inc) /= -1) then ! not constant head 

                   ! flux; deriv wrt _radius_ of other inclusion
                   dPot_dRg(1:M,0:N) = - q(par)/besskRcm(1:M,0:N,oth)*&
                        &(besskRgm(1:M,-1:N-1,oth,inc) + besskRgm(1:M,1:N+1,oth,inc))* &
                        &( amv(1:M,0:N,oth)*cosnPgm(1:M,0:N) + bmv(1:M,0:N,oth)*sinnPgm(1:M,0:N))   

                   ! flux; deriv wrt to _angle_ of other inclusion
                   dPot_dPg(1:M,0:N) = spread(rk(0:N),1,M)*&
                        & besskRgm(1:M,0:N,oth,inc)/ besskRcm(1:M,0:N,oth)*&
                        &(-amv(1:M,0:N,oth)*sinnPgm(1:M,0:N) + bmv(1:M,0:N,oth)*cosnPgm(1:M,0:N))
       
                   ! project onto general Cartesian coords
                   dPot_dX(:) = &
                       & cos(CIPgm(:,oth,inc))*                 sum(dPot_dRg(:,:),dim=2) - &
                       & sin(CIPgm(:,oth,inc))/CIRgm(:,oth,inc)*sum(dPot_dPg(:,:),dim=2)
                   dPot_dY(:) = &
                       & sin(CIPgm(:,oth,inc))*                 sum(dPot_dRg(:,:),dim=2) + &
                         cos(CIPgm(:,oth,inc))/CIRgm(:,oth,inc)*sum(dPot_dPg(:,:),dim=2)
          
                   ! project onto local radial coords (add to running total)
                   FluxInclBg(:) = FluxInclBg(:) + dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
                end if
        
             !@@@@ other inside current @@@@@
             elseif (CIInclIn(inc,oth) .and. CIcalcin(inc)) then 
                if (iteration == 1) then
                   besskRcm(1:M,0:N,oth) = spread(besselk(CIr(oth)*q(inc),0,N+1),1,M)
                   besskRgm(1:M,0:N+1,oth,inc) = besselk(CIRgm(:,oth,inc)*q(inc),0,N+2)
                   besskRgm(1:M,-1,oth,inc) = besskRgm(1:M,+1,oth,inc)
                end if

                cosnPgm(1:M,1:N) = cos(outer_prod(CIPgm(1:M,oth,inc),rk(1:N)))
                cosnPgm(1:M,0) = RONE
                sinnPgm(1:M,1:N) = sin(outer_prod(CIPgm(1:M,oth,inc),rk(1:N)))
                sinnPgm(1:M,0) = RZERO

                !#########################################
                ! head effects on inside of this inclusion
                !######################################### 
                if (CIibnd(inc) /= +1) then ! not constant flux
!!$                   print *, 'db5'
                   ! potential effects of other inlcusions
                   PotInclIn(1:M,0:N) = PotInclIn(1:M,0:N) + &
                        & besskRgm(1:M,0:N,oth,inc)/besskRcm(1:M,0:N,oth)* &
                        & (amv(1:M,0:N,oth)*cosnPgm(1:M,0:N) + &
                        &  bmv(1:M,0:N,oth)*sinnPgm(1:M,0:N))
                end if

                !#########################################
                ! flux effects on inside of this inclusion
                !#########################################
                if (CIibnd(inc) /= -1) then ! not constant head 

                   ! flux; deriv wrt radius of other inclusion
                   dPot_dRg(:,0:N) = - q(inc)/besskRcm(1:M,0:N,oth)* &
                        & (besskRgm(1:M,-1:N-1,oth,inc) + besskRgm(1:M,1:N+1,oth,inc)) *&
                        & (amv(:,0:N,oth)*cosnPgm(:,0:N) + bmv(:,0:N,oth)*sinnPgm(:,0:N))  

                   ! flux; deriv wrt to angle of other inclusion
                   dPot_dPg(:,0:N) = spread(rk(0:N),1,M)*&
                        & besskRgm(1:M,0:N,oth,inc)/besskRcm(1:M,0:N,oth)*&
                        &(-amv(:,0:N,oth)*sinnPgm(:,0:N) + bmv(:,0:N,oth)*cosnPgm(:,0:N))

                   ! project onto general Cartesian coords
                   dPot_dX(:) = &
                       & cos(CIPgm(:,oth,inc))*                 sum(dPot_dRg,dim=2) - &
                       & sin(CIPgm(:,oth,inc))/CIRgm(:,oth,inc)*sum(dPot_dPg,dim=2)
                   dPot_dY(:) = &
                       & sin(CIPgm(:,oth,inc))*                 sum(dPot_dRg,dim=2) + &
                       & cos(CIPgm(:,oth,inc))/CIRgm(:,oth,inc)*sum(dPot_dPg,dim=2)

                   ! project onto local radial coords
                   FluxInclIn(:) = FluxInclIn(:) + dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
                end if
             end if !if CIInclIn
          end do OTHER
    
          !^^^^^^^^^^ effects of parental inclusion ^^^^^^^^^^
          PotInclPar = CZERO; FluxInclPar = CZERO; PotInclParFlux = CZERO
          if (CIInclUp(inc) /= 0) then 

             ! pre-compute Bessel-I
             if (iteration == 1) then
                bessiRcm(1:M,0:N,par) = spread(besseli(CIr(par)*q(par),0,N+1),1,M)
                bessiRgm(1:M,0:N+1,par,inc) = besseli(CIRgm(:,par,inc)*q(par),0,N+2)
                bessiRgm(1:M,-1,par,inc) = bessiRgm(:,+1,par,inc)
             end if

             cosnPgm(1:M,1:N) = cos(outer_prod(CIPgm(1:M,par,inc),rk(1:N)))
             cosnPgm(1:M,0) = RONE
             sinnPgm(1:M,1:N) = sin(outer_prod(CIPgm(1:M,par,inc),rk(1:N)))
             sinnPgm(1:M,0) = RZERO

             !##################################
             ! parental head effects on outside
             !##################################
             if (CIibnd(inc) /= +1) then ! not specified flux
                ! potential effects of parent inclusion
                PotInclPar(1:M,0:N) = &
                     & bessiRgm(1:M,0:N,par,inc)/bessiRcm(1:M,0:N,par)* &
                     & (cmv(1:M,0:N,par)*cosnPgm(1:M,0:N) + &
                     &  dmv(1:M,0:N,par)*sinnPgm(1:M,0:N))
             end if

             ! constant flux over area of inclusion [[ need to generalize this ]]
             PotInclParFlux(1:M) = CircTimeArea(par,p)*CIarea(par)*sv(par)/q(par)**2

             !##################################
             ! parental flux effects on outside
             !##################################
             if (CIibnd(inc) /= -1) then ! not specified head

                ! flux; deriv wrt radius of parent inclusion (sign difference
                ! with analog for K Bessel functions)
                dPot_dRg(:,0:N) = q(par)/ bessiRcm(1:M,0:N,par)* &
                     & (bessiRgm(:,-1:N-1,par,inc) + bessiRgm(:,1:N+1,par,inc))* &
                     & (cmv(:,0:N,par)*cosnPgm(:,:) + dmv(:,0:N,par)*sinnPgm(:,:))

                ! flux; deriv wrt to angle of parent inclusion
                dPot_dPg(:,0:N) = spread(rk(0:N),1,M)*&
                     & bessiRgm(:,0:N,par,inc)/bessiRcm(1:M,0:N,par)*&
                     & (-cmv(:,0:N,par)*sinnPgm(:,:) + dmv(:,0:N,par)*cosnPgm(:,:))

                ! project onto general Cartesian coords
                dPot_dX(:) =&
                    & cos(CIPgm(:,par,inc))*                 sum(dPot_dRg,dim=2) - &
                    & sin(CIPgm(:,par,inc))/CIRgm(:,par,inc)*sum(dPot_dPg,dim=2)
                dPot_dY(:) = &
                    & sin(CIPgm(:,par,inc))*                 sum(dPot_dRg,dim=2) + &
                    & cos(CIPgm(:,par,inc))/CIRgm(:,par,inc)*sum(dPot_dPg,dim=2)
        
                ! project onto local radial coords
                FluxInclPar(1:M) = dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
             end if
          end if

          !^^^^^^^^^^^^^^^^^^^^ setup RHS (bm) ^^^^^^^^^^^^^^^^^^^^
          !$$$$$ head and flux matching $$$$$
          select case (CIibnd(inc))
          case (0)
             
             ! this assumes you _always_ want to calculate inside a
             ! matching inclusion - inside and outside are calcualted together

             bm(1:M,inc) = &
                  & (PotWellBg(:,inc) + sum(PotInclBg(:,0:N),dim=2) + PotCElmBg(:,inc) + &
                  & sum(PotInclPar,dim=2) + PotInclParFlux(:))/kv(par) - &
                  & (PotWellIn(:,inc) + sum(PotInclIn(:,0:N),dim=2) + PotCElmIn(:,inc) + &
                  & CircTimeArea(inc,p)*CIarea(inc)*sv(inc)/q(inc)**2)/kv(inc)
             bm(M+1:2*M,inc) = &
                  & FluxWellBg(:,inc) + FluxInclBg(:) + FluxCElmBg(:,inc) + &
                  & FluxInclPar(:) - FluxWellIn(:,inc) - FluxInclIn(:) - &
                  & FluxCElmIn(:,inc)

             ! calculate least squares using generalized inverse 
             xm(1:4*N+2) = matvec_sor(om,Gm(1:4*N+2,1:2*M,inc),bm(1:2*M,inc),xm1(1:4*N+2))

             ! put results into common named arrays
             a(0:N,inc) = xm(1    :N+1)    ! N+1 terms
             b(1:N,inc) = xm(N+2  :2*N+1)  ! N terms
             c(0:N,inc) = xm(2*N+2:3*N+2)
             d(1:N,inc) = xm(3*N+3:4*N+2)
      
          !$$$$$ specified head $$$$$
          case(-1)

             !<<<< outside specified head inclusion >>>>
             ! currently assumes specified value inside = outside value
             bm(1:M,inc) = CIspec(inc)*kv(par)*CircTimeBdry(inc,p) - &
                  & PotWellBg(:,inc) - sum(PotInclBg(1:M,0:N),dim=2) - PotCElmBg(:,inc) - &
                  &  sum(PotInclPar,dim=2) - PotInclParFlux(:)
   
             loM = 1;  hiM = M;
             loN = 1;  hiN = 2*N+1;
             xm(loN:hiN) = matvec_sor(om,Gm(loN:hiN,loM:hiM,inc),bm(loM:hiM,inc),xm1(loN:hiN))
      
             a(0:N,inc) = xm(loN:N+1) ! N+1 terms
             b(1:N,inc) = xm(N+2:hiN) ! N terms
      
             !>>>>> inside specified head inclusion <<<<<
             ! currently assumes specified value inside = outside value
             if (CIcalcin(inc)) then
                bm(1:M,inc) = CIspec(inc)*kv(inc)*CircTimeBdry(inc,p) - &
                     & PotWellIn(:,inc) - sum(PotInclIn(1:M,0:N),dim=2) - PotCElmIn(:,inc) - &
                     & CircTimeArea(inc,p)*CIarea(inc)*sv(inc)/q(inc)**2
        
                loM = 1;      hiM = M;
                loN = 2*N+2;  hiN = 4*N+2;
                xm(loN:hiN) =  matvec_sor(om,Gm(loN:hiN,loM:hiM,inc),bm(loM:hiM,inc),xm1(loN:hiN))
                
                c(0:N,inc) = xm(loN  :3*N+2) ! N+1 terms
                d(1:N,inc) = xm(3*N+3:hiN)   ! N terms
             end if
      
          !$$$$$ specified flux $$$$$
          case(+1)

             !<<<< outside specified flux inclusion >>>>
             bm(M+1:2*M,inc) = CIspec(inc)*CircTimeBdry(inc,p)/(TWOPI*CIr(inc))  - &
                  & FluxWellBg(:,inc) - FluxInclBg(:) - FluxCElmBg(:,inc) - FluxInclPar(:)
   
             loM = M+1;  hiM = 2*M;
             loN = 1;    hiN = 2*N+1;
             xm(loN:hiN) =  matvec_sor(om,Gm(loN:hiN,loM:hiM,inc),bm(loM:hiM,inc),xm1(loN:hiN))
      
             a(0:N,inc) = xm(loN:N+1) ! N+1 terms
             b(1:N,inc) = xm(N+2:hiN) ! N terms
      
             !>>>>> inside specified flux inclusion <<<<<
             if (CIcalcin(inc)) then
                bm(M+1:2*M,inc) = CIspec(inc)*CircTimeBdry(inc,p)/(TWOPI*CIr(inc)) - &
                     & FluxWellIn(:,inc) - FluxInclIn(:) - FluxCElmIn(:,inc)
        
                loM = M+1;     hiM = 2*M;
                loN = 2*N+2;   hiN = 4*N+2;
                xm(loN:hiN) =  matvec_sor(om,Gm(loN:hiN,loM:hiM,inc),bm(loM:hiM,inc),xm1(loN:hiN))
        
                c(0:N,inc) = xm(loN  :3*N+2) ! N+1 terms
                d(1:N,inc) = xm(3*N+3:hiN)   ! N terms
             end if
          end select

          ! "spread" arrays for vectorizing wrt M and N simultaneously
          amv(1:M,0:N,inc) = spread(a(0:N,inc),dim=1,ncopies=M)
          bmv(1:M,0:N,inc) = spread(b(0:N,inc),dim=1,ncopies=M)
          cmv(1:M,0:N,inc) = spread(c(0:N,inc),dim=1,ncopies=M)
          dmv(1:M,0:N,inc) = spread(d(0:N,inc),dim=1,ncopies=M)

       end do THIS

       !! update solution for wells with wellbore storage
       !! effects of other elements on wells sporting wellbore storage and skin
       if(any(WLstor)) then
          
!!$          WLstorCoeff1 = WLstorCoeff

          !! effects on wells and circular elements on wells with wellbore storage
!!$          call OnStor(p,a,b,c,d,PotOnStor,FluxOnStor)

          ! casing radius assumed equal to wellscreen radius for now
!!$          rc(1:WLnum) = WLr(1:WLnum)

          !! compute coefficients for wellbore storage wells
!!$          do k=1,WLnum
!!$             if (WLstor(k)) then
!!$                WLstorCoeff(k) = (PI*WLr(k)*FluxOnStor(k)*&
!!$                     & (2.0_DP*kv(CIWellUp(k))/(p*rc(k)**2) + WLDSkin(k)) - PotOnStor(k)*PI)
!!$             end if
!!$          end do
          WLstorCoeff(:) = (0.0_DP,0.0_DP)
          
!!$          print *, WLstorCoeff1 ,'<-old, new->',WLstorCoeff

          !! re-compute well effects on circular elements
          !! using new coefficient for wellbore storage
          call wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
       end if       

       xm1 = xm ! last iteration's solution - for SOR

       ! calcualte max absolute change from last iteration
       residual = max(maxval(abs(a-a1)),maxval(abs(b-b1)),&
                    & maxval(abs(c-c1)),maxval(abs(d-d1))) !!,&
!!$                    & maxval(abs(WLstorcoeff-WLstorcoeff1)))

    end do SOR ! while (SOR loop)

    ! pack results into coeff matrix to pass back to calling routine
    coeff(0:N,1:ni)         = a(0:N,1:ni)
    coeff(N+1:2*N,1:ni)     = b(1:N,1:ni)
    coeff(2*N+1:3*N+1,1:ni) = c(0:N,1:ni)
    coeff(3*N+2:4*N+1,1:ni) = d(1:N,1:ni)
    
  end subroutine matchIter

!!$  !########################################
!!$  ! free memory related to points along circular inclusions 
!!$  ! (not needed after coeff are calculated)
!!$  subroutine freematchmem()
!!$    use shared_matching_data
!!$    deallocate(CIXcm, CIYcm, CIXom, CIYom, CIRwm, CIPwm, CIRgm, CIPgm, CIPcm)
!!$  end subroutine freematchmem

!##################################################
! functions for hiding BLAS routines, allowing switching
!##################################################

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! matrix-matrix multiply, conjg transpose of first matrix
  function matmulc1(a,b) result(c)
    use constants, only : DP, CZERO, CONE
    complex(DP), dimension(:,:), intent(in) :: A,B
    complex(DP), dimension(size(A,dim=2),size(B,dim=2)) :: C

    ! A and B are both rectangular; A'*B=C is square

    ! BLAS version
    integer :: n,k
    external zgemm
    n = size(B,dim=2); k = size(B,dim=1) ! m=n (left out)

    call zgemm('Conjg Trans','No Trans',n,n,k,CONE,A,k,B,k,CZERO,C,n)
!!$
!!$    ! f90 version
!!$    c = matmul(transpose(conjg(a)),b)

  end function matmulc1
  
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! matrix-matrix multiply, conjg transpose of second matrix
  function matmulc2(a,b) result(c)
    use constants, only : DP, CZERO, CONE
    complex(DP), dimension(:,:), intent(in) :: A,B
    complex(DP), dimension(size(A,dim=1),size(B,dim=1)) :: C

    ! A is square, B is rectangular; A*B'=C is same shape as B

    ! BLAS version
    integer :: m,n
    external zgemm
    m = size(A,dim=1); n = size(B,dim=1) ! k=m (left out)

    call zgemm('No Trans','Conjg Trans',m,n,m,CONE,A,m,B,n,CZERO,C,m)
!!$
!!$    ! f90 version
!!$    c = matmul(a,transpose(conjg(b)))

  end function matmulc2


  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! length squared (norm, N) of complex vector
  ! used in construction of inv(A'*A)
  function length_sq(x) result(z)
    use constants, only : DP
      ! complex vector (column of Am)
    complex(DP), dimension(:), intent(in) :: x
      ! length squared
    real(DP) :: z

   ! BLAS version
    complex(8) :: zdotc 
    external zdotc
    z = real(zdotc(size(x),x,1,x,1)) ! optimized BLAS version

!!$    ! f90 version
!!$    z = real(dot_product(x,x))

  end function length_sq

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! matrix-vector multiplication
  ! this operation performed _each iteration_
  function matvec_sor(omega,G,b,oldx) result(x)
    use constants, only : DP, CONE
    
      ! generalized inverse => inv(A'*A)*A'
    complex(DP), dimension(:,:), intent(in) :: G 
      ! right hand side "data" vector
    complex(DP), dimension(:), intent(in) :: b
      ! previous iteration's solution
    complex(DP), dimension(:), intent(in) :: oldx
      ! SOR parameter (0< omega <= 2)
    complex(DP), intent(in) :: omega
      ! solution vector (gen. Fourier coefficients)
    complex(DP), dimension(size(G,dim=1)) :: x

    ! BLAS version
    integer :: m,n
    external zgemv
    m = size(G,dim=1)
    n = size(G,dim=2)
    x = oldx ! BLAS works on x 'in-place'
    call zgemv('No transpose',m,n,omega,G,m,b,1,CONE-omega,x,1) 

!!$    ! f90 version
!!$    x = (CONE-omega)*oldx + omega*matmul(G,b) 

  end function matvec_sor

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! outer product, for vectorizing operations 
  ! in two directions simultaneously
  function outer_prod(a,b) result(x)
    use constants, only : DP
    
    real(DP), dimension(:), intent(in) :: a
    real(DP), dimension(:), intent(in) :: b
    real(DP), dimension(size(a),size(b)) :: x
    
    x = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

  end function outer_prod

end module matching_old
