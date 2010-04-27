! $Id: ellipse_matching.f90,v 1.3 2007/09/25 00:49:09 kris Exp kris $
module ellipse_match
  implicit none

  private
  public :: mathieulinehead, pointhead

contains

  !! determine strengths for elliptical elements (no circles for now)
  function ellipse_only_matching(p,calcX,calcY, LINalpha, LINEalpha, LINq, LINsf, &
         & LINeta0, calc_head, calc_velx, LINk, LINgamma, &
         & LINEK, match, flux,n,ms) result(f)
    use mathieu_functions_norm, only : mmatce, mmatse, mmatDce, mmatDse, mmatKeKo, mmatDKeDKo, &
         & mmatDIeDIo, mmatIeIo
    use shared_mathieu, only : A,B,q
    use mcn_matrix_method, only : mcn_eigenvalues
    use complex_bessel, only : cbesk
    use constants, only : DP, TWOPI, PI, CZERO, RONE

    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p)) :: f

    ! elliptical coords, semi-focal length
    real(DP) :: psi, eta, rw1 !!, theta
    complex(DP) :: W, well, calcpt, phiDr, phiDeta, phiDpsi, phiWell

    !! things which were in calcshareddata module
    logical, intent(in) :: flux, match
    ! parameters of the line source
    real(DP), intent(in) :: LINq, LINalpha, LINEalpha, LINk, &
         & LINsf, LINeta0, LINgamma, LINEk
    logical, intent(in) :: calc_head, calc_velx
    real(DP), intent(in) :: calcX, calcY

    integer, intent(in) :: n,ms

    ! mathieu parameter (McLachlan definiton)
    complex(DP), dimension(size(p),2) :: qq
    integer :: nP, i, k, j, nz, ierr

    ! m: number of matching points
    ! n: number of terms in generalized fourier series
    ! ms: size of "infinite" matrix 
    !   used to estimate eigen{values,vectors} (calculated below)
    integer, save :: M
    integer, save, allocatable :: vi(:)
    real(DP), save, allocatable :: sig(:), denom(:)

    ! coefficient matrix used in least-squares process for finding x (strengths of gFc)
    complex(DP), allocatable :: Am(:,:)
    complex(DP), allocatable :: bm(:)        ! RHS vector
    real(DP), allocatable :: Diag(:)     ! norm vector

    complex(DP), allocatable :: Zw(:), Rcw(:) !!, Pw(:)
    real(DP), allocatable :: rw(:) ! distance to well at 1,1

    ! for saving results of radial mf routines (since they are subroutines)
    complex(DP), allocatable :: Ko(:), DKo(:), Ko0(:), Ke(:), DKe(:), &
         & Ke0(:), Io(:), Ie(:), Io0(:), Ie0(:), DIe(:), DIo(:), DKe0(:), DKo0(:), &
         & DKein(:), DKoin(:), Ke0in(:), Ko0in(:), Kein(:), Koin(:), &
         & Je(:), DJe0(:), Ye(:), DYe0(:)
    complex(DP) :: kw(0:1), qb
    complex(DP) , allocatable :: matcein(:,:), matceout(:,:), matsein(:,:), matseout(:,:)

    ! vectors of a_n and b_n for each value of p
    complex(DP), allocatable, save :: x(:,:)

    ! is this the first time through this subroutine?
    logical, save :: first = .true. 

    ! Mathieu coefficients (eigenvectors) for each value of p
    complex(DP), allocatable, save :: AA(:,:,:,:,:), BB(:,:,:,:,:)
    real(DP) :: metric ! scale factors for coordinate system
    real(DP), allocatable, save :: psiM(:) ! vector of matching locations

    well = cmplx(1.0,1.0,DP)
    nP = size(p)

    ! location of calculation point in complex Cartesian coords
    calcpt = cmplx(calcX,calcY,DP)

    ! location in complex elliptical coords
    W = cacosh(calcpt/LINsf)

    ! radial and angular elliptical components
    eta = real(W)  ! radial
    psi = aimag(W) ! azimuthal

    ! mathieu parameter (McLachlan definition)
    ! NB: negative removed from definition !!
    ! (mathieu function re-defined accordingly)
    qq(1:nP,1) = -LINsf**2*p(1:np)/(LINalpha*4.0_DP)  !! outside q
    if (match) then
       qq(1:nP,2) = -LINsf**2*p(1:np)/(LINEalpha*4.0_DP)  !! inside q
    else
       qq(1:nP,2) = cmplx(0.0,0.0,DP)
    end if
    
    allocate(Ko(0:N), DKo(0:N), Ko0(0:N), Ke(0:N), DKe(0:N), Ke0(0:N), &
         & Io(0:N), Ie(0:N), Io0(0:N), Ie0(0:N), DIe(0:N), DIo(0:N), &
         & DKe0(0:N), DKo0(0:N))

    allocate(DKein(0:N), DKoin(0:N), Ke0in(0:N), Ko0in(0:N), &
         & Kein(0:N), Koin(0:N), Je(0:N), DJe0(0:N), Ye(0:N), DYe0(0:N))
    
    if(first) then
       M = 2*(2*N+1)
       allocate(A(MS,0:MS-1,2),B(MS,0:MS-1,2),AA(MS,0:MS-1,2,nP,2),BB(MS,0:MS-1,2,nP,2))
       print *, '=======================================>'
       
       if(match) then
          print *, 'matching with well at (1,1) of strength:', LINq
          print *, 'inside K:',LINEk,' outside K:',LINk
          print *, 'inside alpha:', LINEalpha, ' outside alpha:', LINalpha
          print *, 'area flux over element:',LINgamma
          
          allocate(Am(1:2*M,0:4*N+1),bm(1:2*M),x(0:4*N+1,1:nP))
       elseif(flux) then
          print *, 'specified flux element; q=',LINq
       else
          print *, 'specified head element; h=',LINq
          allocate(Am(1:M,0:2*N),bm(1:M),Diag(0:2*N),x(0:4*N+1,1:nP))
       end if
       print *, 'N:',N,' M: ',M
       print *, 'matrix size:',MS
       
       allocate(psiM(1:M), rw(1:M), Zw(1:M), Rcw(1:M), vi(0:N))
       allocate(matcein(1:M,0:N),matceout(1:M,0:N),&
            & matsein(1:M,1:N),matseout(1:M,1:N))
       vi(0:N) = (/ (i, i=0,N) /)

       ! matching points along circumference of ellipse
       psiM(1:M) = linspace(-PI,+PI,M)
       allocate(sig(0:MS),denom(0:MS))
       sig(0:MS) = real((/ ((-1)**i, i=0,MS) /),DP)
       denom(0:MS) = real((/ (1 - (2*i)**2 ,i=0,MS) /),DP)

       ! for checking size of q vs size of matrix
       print *, 'largest term in Matheiu parameter vector:',maxval(abs(qq))
    else
       !! variables re-used in later calls, but only one term
       allocate(rw(1),Rcw(1)) !! ,Pw(1)
    end if

    do i = 1,np
       if (first) then
          write(*,'(A)',advance='no') '*'
          ! calcualte eigenvectors for a single value of p (McLachlan norm)
          call mcn_eigenvalues(qq(i,1),AA(1:MS,0:MS-1,1:2,i,1),&
               & BB(1:MS,0:MS-1,1:2,i,1),1) 

          ! compute similar properties for inside, if matching
          if (match) then
             call mcn_eigenvalues(qq(i,2),AA(1:MS,0:MS-1,1:2,i,2),&
               & BB(1:MS,0:MS-1,1:2,i,2),1)
          else
             AA(:,:,1:2,1:np,2) = cmplx(0.0,0.0,DP)
             BB(:,:,1:2,1:np,2) = cmplx(0.0,0.0,DP)
          end if
          
          ! save Mathieu coefficients and corresponding q+ into module variables
          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
          q = qq(i,1)

          if(match .or. (.not. flux))then
             Am(:,:) = cmplx(0.0,0.0,DP)
             bm(:) = cmplx(0.0,0.0,DP)
          end if
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          if ((.not. flux) .and. (.not. match)) then ! specified head element

             ! normalizing terms cancel out radial Mathieu functions

             Am(1:M,0) = mmatce(0,psiM(1:M))/LINk   ! a_0
             do k=1,N
                Am(1:M,k) =   mmatce(k,psiM(1:M))/LINk ! a_n
                Am(1:M,N+k) = mmatse(k,psiM(1:M))/LINk ! b_n
             end do

             ! RHS vector for constant head line
             bm(1:M) = LINq/p(i)

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
          elseif (match) then ! matching element with a well in the background

             !!@@@@@@@@@@ pre-compute outside stuff
             call mmatDKeDKo(vi(0:N),LINeta0,DKe,DKo)
             call   mmatKeKo(vi(0:N),LINeta0,Ke0,Ko0)

             matceout(1:M,0) = mmatce(0,psiM(1:M))
             do k = 1,N
                matceout(1:M,k) = mmatce(k,psiM(1:M))
                matseout(1:M,k) = mmatse(k,psiM(1:M))
             end do

             !!@@@@@@@@@@ pre-compute inside stuff                          
             A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,2)
             B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,2)
             q = qq(i,2)

             call mmatDIeDIo(vi(0:N),LINeta0,DIe,DIo)
             call   mmatIeIo(vi(0:N),LINeta0,Ie0,Io0)
             call mmatDKeDKo(vi(0:N),LINeta0,DKein,DKoin)  !! need K fcns inside too?
             call   mmatKeKo(vi(0:N),LINeta0,Ke0in,Ko0in)

             matcein(1:M,0) = mmatce(0,psiM(1:M))
             do k = 1,N
                matcein(1:M,k) = mmatce(k,psiM(1:M))
                matsein(1:M,k) = mmatse(k,psiM(1:M))
             end do
             
             ! head matching
             Am(1:M,0) =     matceout(1:M,0)/LINk    ! a_0
             Am(1:M,2*N+1) = -matcein(1:M,0)/LINEk   ! c_0

             ! flux matching
             Am(M+1:2*M,0) =     matceout(1:M,0)*DKe(0)/Ke0(0) ! a_0
             Am(M+1:2*M,2*N+1) = -matcein(1:M,0)*&
                  &(DIe(0)/Ie0(0) + DKein(0)/Ke0in(0)) ! c_0

             do k = 1,N
                ! head matching
                Am(1:M,k) =   matceout(1:M,k)/LINk ! a_n
                Am(1:M,N+k) = matseout(1:M,k)/LINk ! b_n

                Am(1:M,2*N+k+1) = -matcein(1:M,k)/LINEk ! c_n
                Am(1:M,3*N+k+1) = -matsein(1:M,k)/LINEk ! d_n

                ! flux matching
                Am(M+1:2*M,k) =   matceout(1:M,k)*DKe(k)/Ke0(k) ! a_n
                Am(M+1:2*M,N+k) = matseout(1:M,k)*DKo(k)/Ko0(k) ! b_n

                Am(M+1:2*M,2*N+k+1) = -matcein(1:M,k)*&
                     &(DIe(k)/Ie0(k) + DKein(k)/Ke0in(k)) ! c_n
                Am(M+1:2*M,3*N+k+1) = -matsein(1:M,k)*&
                     &(DIo(k)/Io0(k) + DKoin(k)/Ko0in(k)) ! d_n
             end do

             ! RHS vector: effect of well at 1,1
             !-------------------------------------------

             ! cartesian coordinates of points along circumference of ellipse
             Zw(1:M) = LINsf*ccosh(cmplx(spread(LINeta0,1,M),psiM(1:M),DP))

             open(unit=99,file='ellipse.bdry')
             do j=1,M
                write(99,*) real(Zw(j)),aimag(Zw(j))
             end do
             write(99,*)  real(Zw(1)),aimag(Zw(1)) !! connect it back up
             close(99)
             
             ! cartesian components of vector from well to pts on bndry of ellipse
             Rcw(1:M) = Zw(1:M) - well

             ! distance from well to pt on circ of ellipse
             rw(1:M) = abs(Rcw(1:M))

             ! q associated with Bessel functions
             qb = sqrt(p(i)/LINalpha)

             do j=1,M
                call cbesk(rw(j)*qb,0.0_DP,1,2,kw(0:1),nz,ierr)
                if (nz + ierr /= 0) print *, 'bessel function error.'

                !! well fcn is pos, becomes neg when moved to RHS
                phiWell = LINq*kw(0)/(TWOPI*p(i)*LINk)

                !! deriv of well fcn is neg, becomes pos when moved to RHS
                phiDr =  -LINq*qb*kw(1)/(TWOPI*p(i)*LINk)

                !! area source is negative, stays neg on RHS
                bm(j) = -phiWell/LINk - LINgamma*LINEalpha/(p(i)*LINEk)

                ! flux effects
                bm(M+j) = -phiDr* &
                     & LINsf*(real(Rcw(j))/rw(j)*sinh(LINeta0)*cos(psiM(j)) + &
                     & aimag(Rcw(j))/rw(j)*cosh(LINeta0)*sin(psiM(j)))
             end do
             
          end if

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          ! solve for coefficients in specified head case

          if ((.not. match) .and. (.not. flux))  then

             ! discrete version of norm calculation
             forall(k = 0:2*N)
                Diag(k) = real(dot_product(Am(1:M,k),Am(1:M,k)))
             end forall

             ! least-squares solution for coefficients
             x(0:2*N,i) = matmul(spread(RONE/Diag(0:2*N),dim=2,ncopies=M)* &
                  & transpose(conjg(Am(1:M,0:2*N))),bm(1:M))

          elseif (match) then
             ! solve normal equations for matching case
             x(0:4*N+1,i) = matmul(matmul(inverse(matmul(&
                  & transpose(conjg(Am(1:2*M,0:4*N+1))),Am(1:2*M,0:4*N+1))),&
                  & transpose(conjg(Am(1:2*M,0:4*N+1)))),bm(1:2*M))
          endif
       endif

       A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
       B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
       q = qq(i,1)

       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       ! use x(:,i) to calculate results
       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       
       phiDr = cmplx(0.0,0.0,DP)
       phiDeta = cmplx(0.0,0.0,DP)

       if (calc_head .and. eta>=LINeta0 .and.(match.or.(.not.flux))) then
          ! head effects of a matching, or specified head element OUTSIDE

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
          q = qq(i,1)

          call mmatKeKo(vi(0:N),eta,    Ke, Ko)
          call mmatKeKo(vi(0:N),LINeta0,Ke0,Ko0)

          F(i) = x(0,i)*mmatce(0,psi)*Ke(0)/Ke0(0) ! a_0
          do k=1,N
             F(i) = F(i) + &
                  & x(k,i)*  mmatce(k,psi)*Ke(k)/Ke0(k) + &
                  & x(N+k,i)*mmatse(k,psi)*Ko(k)/Ko0(k)
          end do
          
          ! add effect of pumping well at 1,1
          if (match) then

             rw1 = abs(calcpt - well)
             qb = sqrt(p(i)/LINalpha)

             call cbesk(rw1*qb, 0.0_DP,1,1,kw(0),nz,ierr)
             if (nz + ierr /= 0) print *, 'bessel function error'

             F(i) = F(i) + LINq*kw(0)/(TWOPI*p(i)*LINk)
          end if

          !! convert solution to head
          F(i) = F(i)/LINk  

       elseif ((.not.calc_head) .and. eta>=LINeta0 .and.(match.or.(.not.flux))) then
          ! normal flux effects of a matching or head element OUTSIDE

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
          q = qq(i,1)

          call mmatDKeDKo(vi(0:N),eta,    DKe,DKo)
          call   mmatKeKo(vi(0:N),LINeta0,Ke0,Ko0)

          phiDeta = x(0,i)*mmatce(0,psi)*DKe(0)/Ke0(0)
          do k=1,N 
             phiDeta = phiDeta + &
                  & x(k,i)*mmatce(k,psi)*DKe(k)/Ke0(k) + &
                  & x(k+N,i)*mmatse(k,psi)*DKo(k)/Ko0(k)
          end do

          call mmatKeKo(vi(0:N),eta,Ke,Ko)

          phiDpsi = x(0,i)*mmatDce(0,psi)*Ke(0)/Ke0(0)
          do k=1,N 
             phiDpsi = phiDpsi + &
                  & x(k,i)*mmatDce(k,psi)*Ke(k)/Ke0(k) + &
                  & x(k+N,i)*mmatDse(k,psi)*Ko(k)/Ko0(k)
          end do

          if (match) then
             rw1 = abs(calcpt - well)
             qb = sqrt(p(i)/LINalpha)

             call cbesk(rw1*qb,1.0_DP,1,1,kw(1),nz,ierr)
             if (nz + ierr /= 0) print *, 'bessel function error'
             phiDr = -LINq*qb*kw(1)/(TWOPI*p(i)*LINk)
          end if

       elseif (match .and. calc_head .and. eta<LINeta0) then
          ! head effects of a matching element INSIDE

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,2)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,2)
          q = qq(i,2)

          call mmatIeIo(vi(0:N),eta,    Ie, Io)
          call mmatIeIo(vi(0:N),LINeta0,Ie0,Io0)
          call mmatKeKo(vi(0:N),eta,    Kein,Koin)
          call mmatKeKo(vi(0:N),LINeta0,Ke0in,Ko0in)

          F(i) = x(2*N+1,i)*mmatce(0,psi)*&
               & (Ie(0)/Ie0(0) + Kein(0)/Ke0in(0))
          do k=1,N
             F(i) = F(i) + &
                  & x(2*N+1+k,i)*mmatce(k,psi)*&
                  & (Ie(k)/Ie0(k) + Kein(k)/Ke0in(k)) + &
                  & x(3*N+1+k,i)*mmatse(k,psi)*&
                  & (Io(k)/Io0(k) + Koin(k)/Ko0in(k))
          end do

          !! area source effects & DisPot -> head
          F(i) = (F(i) - LINEalpha*LINgamma/p(i))/LINEk

       elseif (match .and. (.not.calc_head) .and. eta<LINeta0) then
          ! normal flux effects of a matching element INSIDE

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,2)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,2)
          q = qq(i,2)

          call mmatDIeDIo(vi(0:N),eta,    DIe,DIo)
          call   mmatIeIo(vi(0:N),LINeta0,Ie0,Io0)
          call mmatDKeDKo(vi(0:N),eta,    DKein,DKoin)
          call   mmatKeKo(vi(0:N),LINeta0,Ke0in,Ko0in)

          phiDeta = x(2*N+1,i)*mmatce(0,psi)* &
               & (DIe(0)/Ie0(0) + DKein(0)/Ke0in(0))
          do k=1,N
             phiDeta = phiDeta + &
                  & x(2*N+1+k,i)*mmatce(k,psi)*&
                  &(DIe(k)/Ie0(k) + DKein(k)/Ke0in(k)) + &
                  & x(3*N+1+k,i)*mmatse(k,psi)*&
                  &(DIo(k)/Io0(k) + DKoin(k)/Ko0in(k)) 
          end do

          call mmatIeIo(vi(0:N),eta,Ie,Io)
          call mmatKeKo(vi(0:N),eta,Kein,Koin)

          phiDpsi = x(2*N+1,i)*mmatDce(0,psi)*&
               &(Ie(0)/Ie0(0) + Kein(0)/Ke0in(0))
          do k=1,N
             phiDpsi = phiDpsi + &
                  & x(2*N+1+k,i)*mmatDce(k,psi)*&
                  & (Ie(k)/Ie0(k) + Kein(k)/Ke0in(k))&
                  & x(3*N+1+k,i)*mmatDse(k,psi)*&
                  & (Io(k)/Io0(k) + Koin(k)/Ko0in(k))
          end do

       elseif((.not.match) .and. flux .and. calc_head) then

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
          q = qq(i,1)

          call mmatJeJo(vi(0:N),eta,Je,Ko)
          call mmatDJeDJo(vi(0:N),0.0_DP,DJe0,DKo0)
          call mmatYeYo(vi(0:N),eta,Ye,Ko)
          call mmatDYeDYo(vi(0:N),0.0_DP,DYe0,DKo0)

          ! head effects for a specified flux line element
!!$          call mmatDKeDKo(vi(0:N),0.0_DP,DKe0,DKo0)
!!$          call   mmatKeKo(vi(0:N),eta,Ke,Ko)

          F(i) = cmplx(0.0,0.0,DP)
          do k=0,floor(N/2.0)
             F(i) = F(i) + linQ/PI*sig(k)* &
                  & sum(sig(0:MS-1)*conjg(A(1:MS,k,1))/denom(0:MS-1))* &
                  & mmatce(2*k,psi)*(Je(2*k)/DJe0(2*k) + (0.0_DP,1.0_DP)*&
                  & Ye(2*k)/DYe0(2*k)) 
          end do
          
          !! apply time behavior (constant) and convert DP -> h
          F(i) = F(i)/(LINk*p(i))

       elseif((.not.match) .and. flux .and. (.not.calc_head)) then
          ! flux effects for a specified flux line element

          A(1:MS,0:MS-1,1:2) = AA(1:MS,0:MS-1,1:2,i,1)
          B(1:MS,0:MS-1,1:2) = BB(1:MS,0:MS-1,1:2,i,1)
          q = qq(i,1)

          call mmatDKeDKo(vi(0:N),0.0_DP,DKe0,DKo0)
          call mmatKeKo(vi(0:N),0.0_DP,Ke0,Ko0)

          phiDeta = cmplx(0.0,0.0,DP)
          do k=0,floor(N/2.0)
             phiDeta = phiDeta + linQ/PI*sig(k)* &
                  & sum(sig(0:MS-1)*conjg(A(1:MS,k,1))/denom(0:MS-1))* &
                  & mmatce(2*k,psi)*Ke(2*k)/DKe0(2*k) 
          end do

          call mmatKeKo(vi(0:N),eta,Ke,Ko)

          phiDpsi = cmplx(0.0,0.0,DP)
          do k=0,floor(N/2.0)
             phiDpsi = phiDpsi + linQ/PI*sig(k)* &
                  & sum(sig(0:MS-1)*conjg(A(1:MS,k,1))/denom(0:MS-1))* &
                  & mmatDce(2*k,psi)*Ke(2*k)/DKe0(2*k) 
          end do

          !! apply time behavior
          phiDeta = phiDeta/p(i)
          phiDpsi = phiDpsi/p(i)
       end if

       if (.not. calc_head) then
          !! project elliptical and polar onto cartesian x or y
          metric = LINsf**2*(0.5_DP*(cosh(2.0_DP*eta) - cos(2.0_DP*psi)))

          if (calc_velx) then
             F(i) = LINsf/metric*(&
                  & phiDeta*sinh(eta)*cos(psi) - phiDpsi*cosh(eta)*sin(psi)) + &
                  & phiDr*real(calcpt-well)/abs(calcpt-well)
             
          else !! y-flux
             F(i) = LINsf/metric*(&
                  & phiDeta*cosh(eta)*sin(psi) + phiDpsi*sinh(eta)*cos(psi)) + &
                  & phiDr*aimag(calcpt-well)/abs(calcpt-well)
          end if
       end if
       
    end do
    if(first) write(*,'(/A)',advance='no') 'continuing->'

    first = .false.
  end function ellipse_only_matching

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  elemental function cacosh(z) result(f)
    use constants, only : DP, RONE

    complex(DP), intent(in) :: z
    complex(DP) :: f

    if(real(z) >= 0.0) then
       f = log(z + sqrt(z**2 - RONE))
    else
       f = -log(z + sqrt(z**2 - RONE))
    end if

  end function cacosh

  elemental function ccosh(z) result(f)
    use constants, only : DP

    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    
    x = real(z)
    y = aimag(z)

    f = cmplx(cosh(x)*cos(y),sinh(x)*sin(y),DP)

  end function ccosh
  

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v

    integer :: i
    real(DP) :: rnum, range

    if(lo >= hi) stop "LINSPACE: lower bound must be less than upper bound."

    rnum = real(num - 1,DP)
    range = hi - lo

    v = (/ ( lo + real(i,DP)*range/rnum, i=0,num-1) /)

  end function linspace

  function inverse(AI) result(INV)
    use constants, only : DP
    
    ! LAPACK LU decomposition 
    INTERFACE 
       SUBROUTINE ZGETRF(M,N,A,LDA,IPIV,INFO)
         INTEGER, intent(in) :: LDA, M, N
         COMPLEX(KIND=8), intent(inout) :: A(M,N)
         INTEGER, intent(inout) :: IPIV(MIN(M,N))
         INTEGER, intent(inout) :: INFO
       END SUBROUTINE ZGETRF
    END INTERFACE

    ! LAPACK inverse calculation from results of LU
    INTERFACE 
       SUBROUTINE ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
         INTEGER, intent(in) :: LDA, N
         COMPLEX(KIND=8), intent(inout) :: A(LDA,N)
         INTEGER, intent(inout) :: IPIV(N)
         INTEGER, intent(in) :: LWORK
         COMPLEX(KIND=8), intent(inout) :: WORK(LWORK)
         INTEGER, intent(inout) :: INFO
       END SUBROUTINE ZGETRI
    END INTERFACE

    integer :: n, ierr
    !! see comment in mcn_eigenvalues for optimum lwork sizes
    integer, parameter :: LWORK = 410
    complex(DP), dimension(:,:), intent(in) :: ai
    complex(DP), dimension(size(ai,1),size(ai,1)) :: inv
    integer, dimension(size(ai,1)) :: indx
    complex(DP), dimension(LWORK) :: work

    indx = 0
    n = size(ai,1)
    inv = ai    

    call zgetrf(n,n,inv,n,indx,ierr)
    if (ierr /= 0) write(*,*) 'error returned from ZGETRF',ierr

    call zgetri(n,inv,n,indx,work,LWORK,ierr)
    if (ierr /= 0) write(*,*) 'error returned from ZGETRI',ierr
    
  end function inverse
  

end module ellipse_match
