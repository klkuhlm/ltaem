! $Id: matching_ellipse_only.f90,v 1.17 2007/09/19 21:13:45 kris Exp kris $

module matching_ellipse_only
  implicit none

  private 
  public :: EllipseInverse_old, MatchIter

contains

  !##################################################
  ! calculate inverse exploiting the  orthogonallity of the columns of Am
  ! this subroutine called once for each value of p before beginning iteration
  subroutine EllipseInverse_old(p,Am,AA,BB,qq)
    use constants, only: DP,PI,TWOPI,CZERO,CONE,CTWO,RONE,RZERO
    use element_specs, only : EIn,EIm,EInum,kv,EImatch,EIInclUp, EIibnd,EIeta,EICalcIn, EIms
    use mathieu_functions
    use shared_mathieu, only : A,B,q
    use shared_matching_el_data, only : EIPcm

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! geometry of this specific problem, in matrix form
    complex(DP), intent(out), dimension(1:2*EIm, 1:4*EIn+2, 1:EInum) :: Am

      ! eigenvectors needed for Mathieu function routines, Mathieu parameter
    complex(DP), dimension(1:EIms,0:EIms-1,1:2,2*EInum), intent(in) :: AA, BB
    complex(DP), dimension(2*EInum), intent(in) :: qq

    ! INTERNAL VARIABLES
      ! vectors of radial Mathieu functions,
      ! Ke & Ko <-> ceout & seout;  Ie & Io <-> cein & sein
    complex(DP), dimension(0:EIn) :: Ke, Ko, Ie, Io
    complex(DP), dimension(0:EIn) :: DKe, DKo, DIe, DIo
      ! angular MF, different Mathieu parameter inside and out
    complex(DP), dimension(1:EIm,0:EIn) :: cein, ceout
    complex(DP), dimension(1:EIm,1:EIn) :: sein, seout
      ! counters
    integer :: ni, inc, par, k, M, N
    integer, dimension(0:EIn) :: vi

    M = EIm; N = EIn; ni = EInum;

    Am = CZERO
    forall (k=0:N) vi(k)=k
    
    !******************** Am matrix setup ********************
    ! calcualte coefficient matrix (Am) which depends on the geometry 
    ! of the problem and p (calcualte once for each inclusion)

    do inc = 1,ni
       if (EImatch(inc)) then
          par = EIInclUp(inc)

          select case (EIibnd(inc))
          !^^^^^^^^^^^^^^^^^^^^ head & flux matching ^^^^^^^^^^^^^^^^^^^^
          case (0)
             ! pre-compute Mathieu functions
             
             !! mathieu functions related to exterior problem
             A = AA(:,:,:,ni+inc);  B = BB(:,:,:,ni+inc)
             q = qq(ni+inc)
             call   mmatKeKo(vi(0:N),EIeta(inc),Ke,Ko)
             call mmatDKeDKo(vi(0:N),EIeta(inc),DKe,DKo)
             ceout(1:M,0) = mmatce(0,EIPcm(1:M))
             do k=1,N
                ceout(1:M,k) = mmatce(k,EIPcm(1:M))
                seout(1:M,k) = mmatse(k,EIPcm(1:M))
             end do
             
             !! mathieu functions related to interior problem
             A = AA(:,:,:,inc);  B = BB(:,:,:,inc)
             q = qq(inc)
             call   mmatIeIo(vi(0:N),EIeta(inc),Ie,Io)
             call mmatDIeDIo(vi(0:N),EIeta(inc),DIe,DIo)
             cein(1:M,0) = mmatce(0,EIPcm(1:M))
             do k=1,N
                cein(1:M,k) = mmatce(k,EIPcm(1:M))
                sein(1:M,k) = mmatse(k,EIPcm(1:M))
             end do

             ! head matching
             Am(1:M, 1:N+1, inc)    =    ceout(1:M,0:N)/kv(par) ! a_n
             Am(1:M, N+2:2*N+1, inc)  =  seout(1:M,1:N)/kv(par) ! b_n
             Am(1:M, 2*N+2:3*N+2, inc) = -cein(1:M,0:N)/kv(inc) ! c_n
             Am(1:M, 3*N+3:4*N+2, inc) = -sein(1:M,1:N)/kv(inc) ! d_n

             ! flux matching
             Am(M+1:2*M, 1:N+1, inc) = ceout(1:M,0:N)* &
                  & spread(DKe(0:N)/Ke(0:N),1,M) ! a_n
             Am(M+1:2*M, N+2:2*N+1, inc) = seout(1:M,1:N)* &
                  & spread(DKo(1:N)/Ko(1:N),1,M) ! b_n
             
             Am(M+1:2*M, 2*N+2:3*N+2, inc) = -cein(1:M,0:N)* &
                  & spread(DIe(0:N)/Ie(0:N),1,M) ! c_n
             Am(M+1:2*M, 3*N+3:4*N+2, inc) = -sein(1:M,1:N)* &
                  & spread(DIo(1:N)/Io(1:N),1,M) ! d_n

          !^^^^^^^^^^^^^^^^^ specified head ^^^^^^^^^^^^^^^^^^^^
          case (-1)
             print *, 'specified head', inc
             A = AA(:,:,:,ni+inc);  B = BB(:,:,:,ni+inc)
             q = qq(ni+inc)
             
             ! head matching
             Am(1:M, 1, inc) = mmatce(0,EIPcm(1:M))
             do k=1,N
                Am(1:M, k+1, inc)   =  mmatce(k,EIPcm(1:M)) ! a_n
                Am(1:M, N+1+k, inc) =  mmatse(k,EIPcm(1:M)) ! b_n
             end do
             
             !>>>>> calculate inside the specified head inclusion <<<<<
             if (EIcalcin(inc)) then
                print *, 'inside sp h',inc
                A = AA(:,:,:,inc);  B = BB(:,:,:,inc)
                q = qq(inc)

                ! head matching
                Am(1:M, 2*N+2, inc) = mmatce(0,EIPcm(1:M))
                do k=1,N
                   Am(1:M, 2*N+2+k, inc) = mmatce(k,EIPcm(1:M)) ! c_n
                   Am(1:M, 3*N+2+k, inc) = mmatse(k,EIPcm(1:M)) ! d_n
                end do
             end if ! if calcin

          !^^^^^^^^^^^^^^^^^^^^ specified flux ^^^^^^^^^^^^^^^^^^^^
          case (+1)
             print *, 'specified flux',inc
             A = AA(:,:,:,ni+inc);  B = BB(:,:,:,ni+inc)
             q = qq(ni+inc)
             
             call   mmatKeKo(vi(0:N),EIeta(inc),Ke,Ko)
             call mmatDKeDKo(vi(0:N),EIeta(inc),DKe,DKo)
             
             ceout(1:M,0) = mmatce(0,EIPcm(1:M))
             do k=1,N
                ceout(1:M,k) = mmatce(k,EIPcm(1:M))
                seout(1:M,k) = mmatse(k,EIPcm(1:M))
             end do
             
             Am(M+1:2*M, 1:N+1, inc) = &
                  & ceout(1:M,0:N)*spread(DKe(0:N)/Ke(0:N),1,M) ! a_n

             Am(M+1:2*M, N+2:2*N+1, inc) = &
                  & seout(1:M,1:N)*spread(DKo(1:N)/Ko(1:N),1,M) ! b_n
             
             !>>>>> calculate inside the specified flux inclusion <<<<<
             if (EIcalcin(inc)) then
                print *, 'inside sp q',inc
                A = AA(:,:,:,inc);  B = BB(:,:,:,inc)
                q = qq(inc)
             
                call mmatIeIo(vi(0:N),EIeta(inc),Ie,Io)
                call mmatDIeDIo(vi(0:N),EIeta(inc),DIe,DIo)
             
                cein(1:M,0) = mmatce(0,EIPcm(1:M))
                do k=1,N
                   cein(1:M,k) = mmatce(k,EIPcm(1:M))
                   sein(1:M,k) = mmatse(k,EIPcm(1:M))
                end do

                Am(M+1:2*M, 2*N+2:3*N+2, inc) = &
                     & cein(1:M,0:N)*spread(DIe(0:N)/Ie(0:N),1,M) ! c_n
                Am(M+1:2*M, 3*N+3:4*N+2, inc) = &
                     & sein(1:M,1:N)*spread(DIo(1:N)/Io(1:N),1,M) ! d_n

             end if
          end select  ! EIibnd
       end if
    end do
  end subroutine EllipseInverse_old
  
  !##################################################
  ! calcuale effects of wells (passive elements)
  ! this subroutine called once for each value of p before beginning iteration
  subroutine wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    use shared_matching_el_data, only : EIPcm, EIRwm, EIPwm
    use constants, only: DP, CZERO
    use element_specs, only : EIm,EInum,WLnum,av,EImatch,EIInclUp,EIWellIn,EIWellBg, &
         & EIeta, EIf
    use wells_el ! generalized well-effect functions

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! potential and flux due to wells inside and in background of each circular element
    complex(DP), dimension(EIm,EInum), intent(out) :: PotWellIn,FluxWellIn,PotWellBg,FluxWellBg

    ! INTERNAL VARIABLES
      ! sqrt(p/alpha)
    complex(DP), dimension(0:EInum) :: kap
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(EIm) :: dPot_dRw, dPot_dX, dPot_dY
    integer :: M, well, inc, par, ni, nw

    M = EIm; ni = EInum; nw = WLnum;
    kap(0:ni) = sqrt(p/av(0:ni))

    PotWellIn = CZERO; FluxWellIn = CZERO;
    PotWellBg = CZERO; FluxWellBg = CZERO;

    !******************** well effects (head and flux) ********************
    do inc = 1,ni
       if (EImatch(inc)) then
          par = EIInclUp(inc);
          do well = 1,nw 
             !^^^^^^^^^^^^^^^^^^^^ well inside this inclusion ^^^^^^^^^^^^^^^^^^^^
             if (EIWellIn(inc,well)) then
                ! potential due to wells
                PotWellIn(1:M,inc) = PotWellIn(1:M,inc) + &
                     &  wellHead(well,p,kap(inc),EIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,kap(inc),EIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(EIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(EIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local elliptical coords and add to total
                FluxWellIn(1:M,inc) = FluxWellIn(1:M,inc) + &
                     &  (dPot_dX(1:M)*sinh(EIeta(inc))*cos(EIPcm(1:M)) &
                     & + dPot_dY(1:M)*cosh(EIeta(inc))*sin(EIPcm(1:M)))*EIf(inc)

             !^^^^^^^^^^^^^^^^^^^^ well in bg of this inclusion ^^^^^^^^^^^^^^^^^^^^
             elseif (EIWellBg(inc,well)) then
                ! potential due to wells
                PotWellBg(1:M,inc) = PotWellBg(1:M,inc) + &
                     &  wellHead(well,p,kap(par),EIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,kap(par),EIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(EIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(EIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local elliptical radial coordinate (eta) and add to total
                FluxWellBg(1:M,inc) = FluxWellBg(1:M,inc) + &
                     &  (dPot_dX(1:M)*sinh(EIeta(inc))*cos(EIPcm(1:M)) &
                     & + dPot_dY(1:M)*cosh(EIeta(inc))*sin(EIPcm(1:M)))*EIf(inc)
             else
                ! if well is inside a different inclusion, do nothing
             end if
          end do
       end if
    end do
  end subroutine wellEffects

  !##############################################################################
  ! this iterative subroutine calculates the generalized Fourier coefficients
  ! for the active elliptical elements which depend on each other
  subroutine matchIter(Am,coeff,p,matchtol,AA,BB,qq)
    use shared_matching_el_data, only : EIPcm,EIEgm,EIPgm
    use constants, only : DP, TWOPI, CZERO, RZERO, RONE
    use wells_el
    use element_specs, only : EIn,EInum,EIm, EIn, EIms,av,kv,sv,WLnum,EIf, &
         & EIeta,EIInclUp,EIInclBg,EIInclIn,EIInclIn,EIibnd,EIarea,EICalcIn,EIspec
    use mathieu_functions
    use shared_mathieu, only : Amat => A, Bmat => B,q

    character(276) :: fmt

    !! least-squares for potentially degenerate case (using QR decomposition, rather
    !! than simpler, but less-stable normal equations)
    interface
       subroutine zgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
            & WORK, LWORK, RWORK, INFO)
         integer, intent(in) ::  M, N, NRHS, LWORK, LDA, LDB
         integer, intent(out) ::  RANK, INFO
         real(kind=8), intent(in) :: RCOND
         real(kind=8), intent(inout) ::   RWORK(5*min(M,N))
         real(kind=8), intent(out) :: S(min(M,N))
         complex(kind=8), intent(inout) ::  A(LDA,N), B(LDB,NRHS), WORK(max(1,LWORK))
       end subroutine zgelss
    end interface

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! generalized Fourier coefficients (previous value passed in as first guess)
    complex(DP), dimension(0:4*EIn+1,EInum), intent(out) :: coeff
      ! ~ orthogonal eigenfunctions used for LS solution
    complex(DP), dimension(1:2*EIm, 1:4*EIn+2, 1:EInum), intent(in) :: Am
      ! scalar Laplace parameter
    complex(DP), intent(in) :: p
      ! iteration tolerance for SOR/Gauss-Seidel iteration
    real(DP), intent(in) :: matchtol
    complex(DP), intent(in), dimension(1:EIms,0:EIms-1,1:2,2*EInum) :: AA,BB
    complex(DP), intent(in), dimension(2*EInum) :: qq

    ! INTERNAL VARIABLES
    complex(DP), dimension(1:2*EIm, 1:4*EIn+2, 1:EInum) :: Am2 !! temporary copy

      ! "data" vector in LS-solution (effects of other elements on current)
    complex(DP), dimension(2*EIm, EInum) :: bm
      ! sqrt(p/alpha)
    complex(DP), dimension(0:EInum) :: kap
      ! vectors of Mathieu functions, filled first iteration
    complex(DP), dimension(EIm,0:EIn,EInum) :: KeEcm, KoEcm, IeEcm, IoEcm
    complex(DP), dimension(EIm,0:EIn,EInum,EInum) :: KeEgm, KoEgm, IeEgm, IoEgm
    complex(DP), dimension(EIm,0:EIn,EInum,EInum) :: DKeEgm, DKoEgm, DIeEgm, DIoEgm
      ! same data as x, but in common name (previous iteration for convergence test)
    complex(DP), dimension(0:EIn,EInum) :: a,a1,c,c1
    complex(DP), dimension(1:EIn,EInum) :: b,b1,d,d1
    complex(DP), dimension(EIm,0:EIn,EInum) :: amv,cmv
    complex(DP), dimension(EIm,1:EIn,EInum) :: bmv,dmv
      ! results from subroutines for wells
    complex(DP), dimension(EIm,EInum) :: PotWellIn,FluxWellIn,PotWellBg,FluxWellBg
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(EIm) ::  dPot_dX, dPot_dY
    complex(DP), dimension(EIm) :: dPot_dEg, dPot_dPg
      ! potential due to other circular elements
    complex(DP), dimension(EIm) :: PotInclBg, PotInclIn, PotInclPar
      ! flux due to other circular elements
    complex(DP), dimension(EIm) :: PotInclParFlux, FluxInclPar, FluxInclBg, FluxInclIn
      ! vector of angular Mathieu functions
    complex(DP), dimension(EIm,0:EIn,EInum,EInum) :: cein, ceout, Dceout, Dcein
    complex(DP), dimension(EIm,1:EIn,EInum,EInum) :: sein, seout, Dseout, Dsein

    !! zgelss things
    !! only keep eigenvectors with singular values > 1% max sing val.
    !! I guess this tolerance is arbitrary, but it only seems to kick in when
    !! near a double point in q
    real(DP) :: rcond = 1.0D-2  
    integer :: rank, info
    integer, parameter :: LWORK = 410
    complex(DP), dimension(LWORK) :: work
    real(DP) :: rwork(5*(4*EIn+2)), s(4*EIn+2)

    real(DP), dimension(EIm,EInum) :: hsq

    integer :: N, M,  ni, iteration, par
    integer :: k, loM,hiM,loN,hiN,tg,sc, numM, numN
    real(DP) :: residual
    integer, dimension(0:EIn) :: vi

    N = EIn; M = EIm; ni = EInum
    kap(0:ni) = sqrt(p/av(0:ni))

    ! initialize
    iteration = 0 ; residual = huge(1.0)
    
    amv(1:M,0:N,1:ni) = CZERO; bmv(1:M,1:N,1:ni) = CZERO
    cmv(1:M,0:N,1:ni) = CZERO; dmv(1:M,1:N,1:ni) = CZERO
    
    forall(k=0:N) vi(k) = k
    
    !! compute metric coefficient for points on boundary of ellipses
    do k=1,ni
       hsq(1:M,k) = EIf(k)**2/2.0_DP*(cosh(2.0_DP*EIeta(k)) - cos(2.0_DP*EIPcm(1:M)))
    end do
    
    if (WLnum > 0) then
       call wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    end if
    
    KeEcm = CZERO; KoEcm = CZERO; KeEgm = CZERO; KoEgm = CZERO;
    ceout = CZERO; seout = CZERO; Dceout = CZERO; Dseout = CZERO;
    IeEcm = CZERO; IoEcm = CZERO; IeEgm = CZERO; IoEgm = CZERO;
    cein = CZERO; sein = CZERO; Dcein = CZERO; Dsein = CZERO;

    !******************** fixed-point iteration loop ********************
    FP: do 
       if (residual <= matchtol) then
          write (*,'(1x,i4,a6,ES10.4)') iteration,' norm=',residual
          exit FP
       end if
  
       if (iteration > 0 .and. mod(iteration,5)==0) then
          write(*,'(1x,i4,a6,ES10.4)') iteration,' norm=',residual
       end if

       iteration = iteration + 1
       a1 = a; b1 = b; c1 = c; d1 = d
       Am2 = Am
      
       TARGET: do tg = 1,ni
!!$          if(iteration == 1) print *, 'target:',tg

          par = EIInclUp(tg) ! parent of this or any other bg inclusions
          PotInclBg = CZERO;  FluxInclBg = CZERO
          PotInclIn = CZERO;  FluxInclIn = CZERO
    
          !^^^^^^^^^^^^^^^^^^^^ loop over other elliptical elements ^^^^^^^^^^^^^^^^^^^^
          SOURCE: do sc = 1,ni
  
             ! @@@@@ source in background of target @@@@@
             if (EIInclBg(sc,tg)) then
!!$                if(iteration == 1) print *, 'BG source:',sc
        
                ! pre-compute exterior Mathieu functions
                if (iteration == 1) then 
!!$                   print *, 'compute exterior MF'
                   Amat = AA(:,:,:,ni+sc); Bmat = BB(:,:,:,ni+sc); q = qq(ni+sc)

                   call mmatKeKo(vi(0:N),EIeta(sc),KeEcm(1,0:N,sc),KoEcm(1,0:N,sc))
                   KeEcm(2:M,0:N,sc) = spread(KeEcm(1,0:N,sc),1,M-1)
                   KoEcm(2:M,1:N,sc) = spread(KoEcm(1,1:N,sc),1,M-1)
                   do k=1,M
                      call mmatKeKo(vi(0:N),EIEgm(k,sc,tg),&
                           & KeEgm(k,0:N,sc,tg),KoEgm(k,0:N,sc,tg))
                   end do

                   ceout(1:M,0,sc,tg) = mmatce(0,EIPgm(1:M,sc,tg))
                   do k=1,N
                      ceout(1:M,k,sc,tg) = mmatce(k,EIPgm(1:M,sc,tg))
                      seout(1:M,k,sc,tg) = mmatse(k,EIPgm(1:M,sc,tg))
                   end do
                end if

                !##########################################
                ! head effects on outside of target inclusion
                !##########################################

                if (EIibnd(tg) /= +1) then ! not specified flux 
!!$                   if(iteration == 1)  print *, 'head effects of',sc,' on',tg
                   ! potential effects of other inlcusions
                   ! accumulate effects over all other inclusions
                   PotInclBg(1:M) = PotInclBg(1:M) + sum(&
                        & KeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*ceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(KoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)*&
                        & bmv(1:M,1:N,sc)*seout(1:M,1:N,sc,tg),dim=2)
                end if

                !##########################################
                ! flux effects on outside of target inclusion
                !##########################################
 
                if (EIibnd(tg) /= -1) then ! not constant head 
                   if(iteration == 1) then
!!$                      print *, 'flux effects of',sc,' on',tg
!!$                      print *, 'compute exterior DMF'
                      Amat = AA(:,:,:,ni+sc); Bmat = BB(:,:,:,ni+sc); q = qq(ni+sc)
                      do k=1,M
                         call mmatDKeDKo(vi(0:N),EIEgm(k,sc,tg),&
                              & DKeEgm(k,0:N,sc,tg),DKoEgm(k,0:N,sc,tg))
                      end do

                      Dceout(1:M,0,sc,tg) = mmatDce(0,EIPgm(1:M,sc,tg))
                      do k=1,N
                         Dceout(1:M,k,sc,tg) = mmatDce(k,EIPgm(1:M,sc,tg))
                         Dseout(1:M,k,sc,tg) = mmatDse(k,EIPgm(1:M,sc,tg))
                      end do
                   end if

                   ! deriv wrt _radius_ of source inclusion
                   dPot_dEg(1:M) = sum(DKeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*ceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(DKoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)* &
                        & bmv(1:M,1:N,sc)*seout(1:M,1:N,sc,tg),dim=2)

                   ! deriv wrt to _angle_ of source inclusion
                   dPot_dPg(1:M) = sum(KeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*Dceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(KoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)*&
                        & bmv(1:M,1:N,sc)*Dseout(1:M,1:N,sc,tg),dim=2)
       
                   ! project onto general Cartesian coords
                   dPot_dX(:) = EIf(sc)/hsq(1:M,sc) *(&
                        & sinh(EIEgm(:,sc,tg))*cos(EIPgm(:,sc,tg))*dPot_dEg(:) -&
                        & cosh(EIEgm(:,sc,tg))*sin(EIPgm(:,sc,tg))*dPot_dPg(:))

                   dPot_dY(:) = EIf(sc)/hsq(1:M,sc) *(&
                        & cosh(EIEgm(:,sc,tg))*sin(EIPgm(:,sc,tg))*dPot_dEg(:) +&
                        & sinh(EIEgm(:,sc,tg))*cos(EIPgm(:,sc,tg))*dPot_dPg(:))
          
                   ! project onto target eta coord (add to running total)
                   FluxInclBg(:) = FluxInclBg(:) + EIf(tg)*(&
                        & dPot_dX(:)*sinh(EIeta(tg))*cos(EIPcm(:)) + &
                        & dPot_dY(:)*cosh(EIeta(tg))*sin(EIPcm(:)))
                end if
        
             !@@@@ source inside target @@@@@
             elseif (EIInclIn(tg,sc) .and. EIcalcin(tg)) then 
!!$                print *, 'inside source:',sc

                ! pre-compute exterior Mathieu functions
                if (iteration == 1) then 
!!$                   print *, 'compute exterior MF'
                   Amat = AA(:,:,:,ni+sc); Bmat = BB(:,:,:,ni+sc); q = qq(ni+sc)

                   call mmatKeKo(vi(0:N),EIeta(sc),KeEcm(1,0:N,sc),KoEcm(1,0:N,sc))
                   KeEcm(2:M,0:N,sc) = spread(KeEcm(1,0:N,sc),dim=1,ncopies=M-1)
                   KoEcm(2:M,0:N,sc) = spread(KoEcm(1,0:N,sc),dim=1,ncopies=M-1)
                   do k=1,M
                      call mmatKeKo(vi(0:N),EIEgm(k,sc,tg),&
                           & KeEgm(k,0:N,sc,tg),KoEgm(k,0:N,sc,tg))
                   end do

                   ceout(1:M,0,sc,tg) = mmatce(0,EIPgm(1:M,sc,tg))
                   do k=1,N
                      ceout(1:M,k,sc,tg) = mmatce(k,EIPgm(1:M,sc,tg))
                      seout(1:M,k,sc,tg) = mmatse(k,EIPgm(1:M,sc,tg))
                   end do
                end if

                !#########################################
                ! head effects on inside of this inclusion
                !######################################### 
                if (EIibnd(tg) /= +1) then ! not constant flux
!!$                   print *, 'head effects of',sc,' on',tg
                   ! potential effects of other inlcusions
                   PotInclIn(1:M) = PotInclIn(1:M) + sum(&
                        & KeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*ceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(KoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)*&
                        & bmv(1:M,1:N,sc)*seout(1:M,1:N,sc,tg),dim=2)
                end if

                !#########################################
                ! flux effects on inside of this inclusion
                !#########################################
                if (EIibnd(tg) /= -1) then ! not constant head  
!!$                   print *, 'flux effects of',sc,' on',tg

                   if(iteration == 1) then
!!$                      print *, 'compute exterior MF'
                      Amat = AA(:,:,:,ni+sc); Bmat = BB(:,:,:,ni+sc); q = qq(ni+sc)
                      do k=1,M
                         call mmatDKeDKo(vi(0:N),EIEgm(k,sc,tg),&
                              & DKeEgm(k,0:N,sc,tg),DKoEgm(k,0:N,sc,tg))
                      end do
                      
                      Dceout(1:M,0,sc,tg) = mmatDce(0,EIPgm(1:M,sc,tg))
                      do k=1,N
                         Dceout(1:M,k,sc,tg) = mmatDce(k,EIPgm(1:M,sc,tg))
                         Dseout(1:M,k,sc,tg) = mmatDse(k,EIPgm(1:M,sc,tg))
                      end do
                   end if

                   ! flux; deriv wrt radius of other inclusion
                   dPot_dEg(:) = sum(DKeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*ceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(DKoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)* &
                        & bmv(1:M,1:N,sc)*seout(1:M,1:N,sc,tg),dim=2)  

                   ! flux; deriv wrt to angle of other inclusion
                   dPot_dPg(:) = sum(KeEgm(1:M,0:N,sc,tg)/KeEcm(1:M,0:N,sc)* &
                        & amv(1:M,0:N,sc)*Dceout(1:M,0:N,sc,tg),dim=2) + &
                        & sum(KoEgm(1:M,1:N,sc,tg)/KoEcm(1:M,1:N,sc)*&
                        & bmv(1:M,1:N,sc)*Dseout(1:M,1:N,sc,tg),dim=2)

                   ! project onto general Cartesian coords
                   dPot_dX(:) = EIf(sc)/hsq(1:M,sc) *(&
                        & sinh(EIEgm(:,sc,tg))*cos(EIPgm(:,sc,tg))*dPot_dEg(:) -&
                        & cosh(EIEgm(:,sc,tg))*sin(EIPgm(:,sc,tg))*dPot_dPg(:))
                   dPot_dY(:) = EIf(sc)/hsq(1:M,sc) *(&
                        & cosh(EIEgm(:,sc,tg))*sin(EIPgm(:,sc,tg))*dPot_dEg(:) +&
                        & sinh(EIEgm(:,sc,tg))*cos(EIPgm(:,sc,tg))*dPot_dPg(:))

                   ! project onto local radial coords
                   FluxInclIn(:) = FluxInclIn(:) + EIf(tg)*(&
                        & dPot_dX(:)*sinh(EIeta(tg))*cos(EIPcm(:)) + &
                        & dPot_dY(:)*cosh(EIeta(tg))*sin(EIPcm(:)))
                end if
             end if !if EIInclIn
          end do SOURCE

    
          !^^^^^^^^^^ effects of parental inclusion ^^^^^^^^^^
          PotInclPar = CZERO; FluxInclPar = CZERO; PotInclParFlux = CZERO
          if (EIInclUp(tg) /= 0) then 
!!$             print *, 'parental source:',par

             ! pre-compute interior Mathieu functions
             if (iteration == 1) then ! shared stuff
!!$                print *, 'compute interior MF' 
                Amat = AA(:,:,:,par); Bmat = BB(:,:,:,par); q = qq(par)
                
                call mmatIeIo(vi(0:N),EIeta(par),IeEcm(1,0:N,par),IoEcm(1,0:N,par))
                IeEcm(2:M,0:N,par) = spread(IeEcm(1,0:N,par),dim=1,ncopies=M-1)
                IoEcm(2:M,0:N,par) = spread(IoEcm(1,0:N,par),dim=1,ncopies=M-1)
                do k=1,M
                   call mmatIeIo(vi(0:N),EIEgm(k,par,tg),&
                        & IeEgm(k,0:N,par,tg),IoEgm(k,0:N,par,tg))
                end do
                
                cein(1:M,0,par,tg) = mmatce(0,EIPgm(1:M,par,tg))
                do k=1,N
                   cein(1:M,k,par,tg) = mmatce(k,EIPgm(1:M,par,tg))
                   sein(1:M,k,par,tg) = mmatse(k,EIPgm(1:M,par,tg))
                end do
             end if

             !##################################
             ! parental head effects on outside of target
             !##################################
             if (EIibnd(tg) /= +1) then ! not specified flux
!!$                print *, 'head effects of parent',par
                ! potential effects of parent inclusion
                PotInclPar(1:M) = sum(&
                     & IeEgm(1:M,0:N,par,tg)/IeEcm(1:M,0:N,par)* &
                     & cmv(1:M,0:N,par)*cein(1:M,0:N,par,tg),dim=2) + &
                     & sum(IoEgm(1:M,1:N,par,tg)/IoEcm(1:M,1:N,par)*&
                     & dmv(1:M,1:N,par)*sein(1:M,1:N,par,tg),dim=2)
             end if

             ! constant flux over area of inclusion [[ need to generalize this ]]
             PotInclParFlux(1:M) = ElTimeArea(par,p)*EIarea(par)*sv(par)/kap(par)**2

             !##################################
             ! parental flux effects on outside
             !##################################
             if (EIibnd(tg) /= -1) then ! not specified head
!!$                print *, 'flux effects of parent',par

                if(iteration == 1) then
!!$                   print *, 'compute interior MF'
                   Amat = AA(:,:,:,par); Bmat = BB(:,:,:,par); q = qq(par)
                   do k=1,M
                      call mmatDIeDIo(vi(0:N),EIEgm(k,par,tg),&
                           & DIeEgm(k,0:N,par,tg),DIoEgm(k,0:N,par,tg))
                   end do
                   
                   Dcein(1:M,0,par,tg) = mmatDce(0,EIPgm(1:M,par,tg))
                   do k=1,N
                      Dcein(1:M,k,par,tg) = mmatDce(k,EIPgm(1:M,par,tg))
                      Dsein(1:M,k,par,tg) = mmatDse(k,EIPgm(1:M,par,tg))
                   end do
                end if

                ! deriv wrt radius of parent inclusion
                dPot_dEg(:) = sum(DIeEgm(1:M,0:N,par,tg)/IeEcm(1:M,0:N,par)* &
                        & cmv(1:M,0:N,par)*cein(1:M,0:N,par,tg),dim=2) + &
                        & sum(DIoEgm(1:M,1:N,par,tg)/IoEcm(1:M,1:N,par)* &
                        & dmv(1:M,1:N,par)*sein(1:M,1:N,par,tg),dim=2)

                ! flux; deriv wrt to angle of parent inclusion
                dPot_dPg(:) = sum(IeEgm(1:M,0:N,par,tg)/IeEcm(1:M,0:N,par)* &
                        & cmv(1:M,0:N,par)*Dcein(1:M,0:N,par,tg),dim=2) + &
                        & sum(IoEgm(1:M,1:N,par,tg)/IoEcm(1:M,1:N,par)*&
                        & dmv(1:M,1:N,par)*Dsein(1:M,1:N,par,tg),dim=2)

                ! project onto general Cartesian coords
                dPot_dX(:) = EIf(par)/hsq(1:M,par) *(&
                        & sinh(EIEgm(:,par,tg))*cos(EIPgm(:,par,tg))*dPot_dEg(:) -&
                        & cosh(EIEgm(:,par,tg))*sin(EIPgm(:,par,tg))*dPot_dPg(:))
                dPot_dY(:) = EIf(par)/hsq(1:M,par) *(&
                        & cosh(EIEgm(:,par,tg))*sin(EIPgm(:,par,tg))*dPot_dEg(:) +&
                        & sinh(EIEgm(:,par,tg))*cos(EIPgm(:,par,tg))*dPot_dPg(:))
        
                ! project onto target eta coord
                FluxInclPar(1:M) = FluxInclIn(:) + EIf(tg)*(&
                        & dPot_dX(:)*sinh(EIeta(tg))*cos(EIPcm(:)) + &
                        & dPot_dY(:)*cosh(EIeta(tg))*sin(EIPcm(:)))
             end if
          end if

          !^^^^^^^^^^^^^^^^^^^^ setup RHS (bm) ^^^^^^^^^^^^^^^^^^^^
          !$$$$$ head and flux matching $$$$$
          
          select case (EIibnd(tg))
          case (0)
             
             ! this assumes you _always_ want to calculate inside a
             ! matching inclusion - inside and outside are calcualted together
             !! can't think of a reason why you would do otherwise

             bm(1:M,tg) = &
                  & (PotWellIn(:,tg) + PotInclIn(:) +  &
                  &  ElTimeArea(tg,p)*EIarea(tg)*sv(tg)/kap(tg)**2)/kv(tg) &
                  &  - (PotWellBg(:,tg) + PotInclBg +  &
                  &  PotInclPar + PotInclParFlux(:))/kv(par) 
             bm(M+1:2*M,tg) = &
                  & FluxWellIn(:,tg) + FluxInclIn(:) - (FluxWellBg(:,tg) + &
                  & FluxInclBg(:) + FluxInclPar(:))

             !! compute least-squares solution using QR decomposition
             call zgelss(2*M,4*N+2,1,Am2(:,:,tg),2*M,bm(:,tg),2*M,s,&
                  & rcond, rank, work, LWORK, rwork, info)
             if(info /= 0) print *, 'ZGELSS: error', info
             
             
             if (rank < 4*N+2 .and. iteration == 1) then
                fmt = '(A,I4,A,I4,A,    (1X,F5.3))'
                write(fmt(14:17),'(I4.4)') 4*N+2
                write(*,fmt) 'rank:',rank, ' of possible:',4*N+2, ' s(:)/s(1):', s/s(1)
             end if

             ! put results into common-named arrays
             a(0:N,tg) = bm(1    :N+1,tg)    ! N+1 terms
             b(1:N,tg) = bm(N+2  :2*N+1,tg)  ! N terms
             c(0:N,tg) = bm(2*N+2:3*N+2,tg)
             d(1:N,tg) = bm(3*N+3:4*N+2,tg)
      
          !$$$$$ specified head $$$$$
          case(-1)

             !<<<< outside specified head inclusion >>>>
             ! currently assumes specified value inside = outside value
             bm(1:M,tg) = EIspec(tg)*kv(par)*ElTimeBdry(tg,p) - &
                  & (PotWellBg(:,tg) + PotInclBg(1:M) + &
                  &  PotInclPar + PotInclParFlux(:))
   
             loM = 1;  hiM = M;  numM = hiM-loM+1
             loN = 1;  hiN = 2*N+1; numN = hiN-loN+1

             call zgelss(numM,numN,1,Am2(loM:hiM,loN:hiN,tg),numM,&
                  & bm(loM:hiM,tg),numM,s(loN:hiN), &
                  & rcond, rank, work, LWORK, rwork(1:5*hiN), info)
             if(info /= 0) print *, 'ZGELSS: error', info
      
             a(0:N,tg) = bm(loN:N+1,tg) ! N+1 terms
             b(1:N,tg) = bm(N+2:hiN,tg) ! N terms
      
             !>>>>> inside specified head inclusion <<<<<
             ! currently assumes specified value inside = outside value
             if (EIcalcin(tg)) then
                bm(1:M,tg) = EIspec(tg)*kv(tg)*ElTimeBdry(tg,p) - &
                     & (PotWellIn(:,tg) + PotInclIn(1:M) +  &
                     & ElTimeArea(tg,p)*EIarea(tg)*sv(tg)/kap(tg)**2)
        
                loM = 1;      hiM = M;     numM = hiM-loM+1
                loN = 2*N+2;  hiN = 4*N+2; numN = hiN-loN+1

                call zgelss(numM,numN,1,Am2(loM:hiM,loN:hiN,tg),numM,&
                     & bm(loM:hiM,tg),numM,s(loN:hiN), &
                     & rcond, rank, work, LWORK, rwork(1:5*hiN), info)
                if(info /= 0) print *, 'ZGELSS: error', info
                
                c(0:N,tg) = bm(loN  :3*N+2,tg) ! N+1 terms
                d(1:N,tg) = bm(3*N+3:hiN,tg)   ! N terms
             end if
      
          !$$$$$ specified flux $$$$$
          case(+1)

             !! need a routine for the circumference of an ellipse, to replace 4f below

             !<<<< outside specified flux inclusion >>>>
             bm(M+1:2*M,tg) = EIspec(tg)*ElTimeBdry(tg,p)/(4.0_DP*EIf(tg)) - &
                  & (FluxWellBg(:,tg) + FluxInclBg(:) + FluxInclPar(:))
   
             loM = M+1;  hiM = 2*M;   numM = hiM-loM+1
             loN = 1;    hiN = 2*N+1; numN = hiN-loN+1

             call zgelss(numM,numN,1,Am2(loM:hiM,loN:hiN,tg),numM,&
                  & bm(loM:hiM,tg),numM,s(loN:hiN), &
                  & rcond, rank, work, LWORK, rwork(1:5*hiN), info)
             if(info /= 0) print *, 'ZGELSS: error', info
      
             a(0:N,tg) = bm(loN:N+1,tg) ! N+1 terms
             b(1:N,tg) = bm(N+2:hiN,tg) ! N terms
      
             !>>>>> inside specified flux inclusion <<<<<
             if (EIcalcin(tg)) then
                bm(M+1:2*M,tg) = EIspec(tg)*ElTimeBdry(tg,p)/(4.0_DP*EIf(tg)) - &
                     & (FluxWellIn(:,tg) + FluxInclIn(:))
        
                loM = M+1;     hiM = 2*M;   numM = hiM-loM+1
                loN = 2*N+2;   hiN = 4*N+2; numN = hiN-loN+1

                call zgelss(numM,numN,1,Am2(loM:hiM,loN:hiN,tg),numM,&
                     & bm(loM:hiM,tg),numM,s(loN:hiN), &
                     & rcond, rank, work, LWORK, rwork(1:5*hiN), info)
                if(info /= 0) print *, 'ZGELSS: error', info
        
                c(0:N,tg) = bm(loN  :3*N+2,tg) ! N+1 terms
                d(1:N,tg) = bm(3*N+3:hiN,tg)   ! N terms
             end if
          end select

          ! "spread" arrays for vectorizing wrt M and N simultaneously
          amv(1:M,0:N,tg) = spread(a(0:N,tg),dim=1,ncopies=M)
          bmv(1:M,1:N,tg) = spread(b(1:N,tg),dim=1,ncopies=M)
          cmv(1:M,0:N,tg) = spread(c(0:N,tg),dim=1,ncopies=M)
          dmv(1:M,1:N,tg) = spread(d(1:N,tg),dim=1,ncopies=M)

       end do TARGET
       
       ! calcualte max absolute change from last iteration
       residual = max(abs(maxval(abs(a-a1))),abs(maxval(abs(b-b1))),&
                    & abs(maxval(abs(c-c1))),abs(maxval(abs(d-d1))))

    end do FP 

    ! pack results into coeff matrix to pass back to calling routine
    coeff(0:N,1:ni)         = a(0:N,1:ni)
    coeff(N+1:2*N,1:ni)     = b(1:N,1:ni)
    coeff(2*N+1:3*N+1,1:ni) = c(0:N,1:ni)
    coeff(3*N+2:4*N+1,1:ni) = d(1:N,1:ni)

!!$    write(*,*) 'a:',maxval(abs(a(0:N,1))),a(0:N,1)
!!$    write(*,*) 'b:',maxval(abs(b(1:N,1))),b(1:N,1)
!!$    write(*,*) 'c:',maxval(abs(c(0:N,1))),c(0:N,1)
!!$    write(*,*) 'd:',maxval(abs(d(1:N,1))),d(1:N,1)
    
  end subroutine matchIter
end module matching_ellipse_only
