! $Id: mathieu_line.f90,v 1.24 2010/04/08 14:46:45 klkuhlm Exp klkuhlm $
module mathieu_line
  use constants, only : DP
  use mathieu_functions, only : mathieu
  implicit none

  type, public :: line_external
     ! everything in this library is public, and set in the driving program

     ! parameters of the line source / aquifer
     real(DP) :: alpha, Ealpha  ! alpha = K/Ss (E means inside ellipse)
     real(DP) :: Ss, ESs, k, Ek  ! specific storage and hydraulic conductivity (horizontal)
     real(DP) :: Q  ! flowrate of a specified-flux line source
     real(DP) :: sf, eta0  ! semi-focal length, elliptical radius
     real(DP) :: gamma  ! area-distributed constant-strength source term
     integer :: N, MS ! number of matching points along circumference and size of inf matrix
     
     ! problem type flags set in the main 
     logical :: flux, match, cartesian, calc_head, calc_velx, first
     complex(DP) :: WW ! complex elliptical coordinates to solve for

     ! leaky/unconfined aquifer stuff
     integer :: leak  ! leaky type (1=CH bdry, 2=no-flow bdry, 3=thickness->infinity)
     logical :: unconfined
     real(DP):: k2,ss2,b2  ! hydraulic conductivity (horizontal), specific storage & thickness of aquitard
     real(DP) :: sy,kz,b ! specific yield, vertical K, main aquifer thickness

     ! well location
     complex(DP) :: well

  end type line_external
  
contains

  function mathieulinehead(p,mf,calcx,calcy,lin) result(f)
    use mathieu_functions
    use complex_bessel, only : cbesk
    use constants, only : DP, TWOPI, PI, PIOV2, CZERO, RONE
    use utility, only : cacosh, outerprod, linspace

    interface  ! LAPACK least-squares routine
       subroutine zgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
            & WORK, LWORK, RWORK, INFO)
         INTEGER, intent(in) ::  M, N, NRHS, LWORK, LDA, LDB
         integer, intent(out) ::  RANK, INFO
         real(kind=8), intent(in) :: RCOND
         real(kind=8), intent(inout) :: RWORK(5*min(M,N))
         real(kind=8), intent(out) :: S(min(M,N))
         COMPLEX(kind=8), intent(inout) ::  A(LDA,N), B(LDB,NRHS), WORK(max(1,LWORK))
       end subroutine zgelss
    end interface

    real(DP), intent(in) :: calcx,calcy
    type(line_external), intent(in) :: lin
    complex(DP), dimension(:), intent(in) :: p    
    type(mathieu), intent(inout), dimension(size(p),2) :: mf
    complex(DP), dimension(size(p)) :: f

    ! coefficient matrix used in least-squares process for finding x (strengths of gFc)
    complex(DP), allocatable :: Am(:,:), bm(:)
    ! vectors of a_n and b_n for each value of p
    complex(DP), allocatable, save :: x(:,:)
    integer, allocatable :: vi(:) ! a vector of integers, counting up from zero

    ! kappa^2 in modified Helmholtz equation, including leaky, unconfined & wave effects, etc.
    complex(DP), allocatable, save :: kappasq(:)
          
    real(DP) :: metricsq ! scale factor for elliptical coordinates
    real(DP), allocatable :: psiM(:) ! vector or matching locations along circumference

    ! elliptical coords, semi-focal length
    real(DP) :: psi, eta, rw1
    complex(DP) :: W, calcpt, phiDr, phiDeta, phiDpsi, phiWell

    !! zgelss things
    real(DP) :: rcond = -1.0_DP
    integer :: rank, info, lwork
    complex(DP), allocatable :: work(:)
    real(DP), allocatable :: rwork(:), s(:)

    integer :: nP, M, i, j, nz, ierr, MS, N, nmax
    real(DP), save, allocatable :: sig(:), denom(:)

    complex(DP), allocatable :: Zw(:), Rcw(:)
    complex(DP), allocatable, save :: qq(:,:)
    real(DP), allocatable :: rw(:) ! distance to well at 1,1
    complex(DP) :: kw(0:1), qb

!!$    write(*,'(6(A,1L))') 'flux:',lin%flux,' match:',lin%match, &
!!$         & ' cart:',lin%cartesian,' c_head:',lin%calc_head,&
!!$         & ' c_velx:',lin%calc_velx,' first:',lin%first

    nP = size(p)
    N = lin%N
    M = 2*(2*N+1)
    lwork = 4*M
    allocate(work(lwork),stat=ierr)
    if ( ierr /= 0 ) then
       write(*,'(A,I0)') 'error allocating: work ',lwork
       stop 222
    end if
    
    ! location of calculation point in complex Cartesian coords
    if (lin%cartesian) then
       calcpt = cmplx(calcx,calcy,DP)

       ! location in complex elliptical coords
       W = cacosh(calcpt/lin%sf)
    else
       ! complex definition of elliptical coords
       W = lin%WW

       ! complex definition of Cartesian coords
       calcpt = cmplx(lin%sf*cosh(real(W))*cos(aimag(W)), &
                   &  lin%sf*sinh(real(W))*sin(aimag(W)))
    end if
    
    eta = real(W)  ! radial elliptical coordinate
    psi = aimag(W) ! azimuthal elliptical coordinate

    if (lin%match) then
!!$       write(*,'(A)') 'matching'
       if (.not. allocated(qq)) then
          allocate(qq(nP,2),stat=ierr)
          if (ierr /= 0) then
             write(*,'(A,(I0,1X))') 'error allocating: qq ',np
             stop 100
          end if
       end if
       qq(1:nP,1) = p(1:np)*lin%sf**2/(lin%alpha*4.0_DP)   !! outside q
       qq(1:nP,2) = p(1:np)*lin%sf**2/(lin%Ealpha*4.0_DP)  !! inside q
    else
       if (lin%first) then
          if(.not. allocated(kappasq)) then
             allocate(kappasq(np),qq(nP,2),stat=ierr)
             if (ierr /= 0) then
                write(*,'(A,I0)') 'error allocating: kappasq ',np
                stop 101
             end if
          end if
          ! incorporate leaky and unconfined effects for line source here
          print *, 'compute leaky parameters'
          kappasq(1:nP) = compute_ELleaky_qv(p,lin)
       end if

       ! mathieu parameter (p/alpha come from leaky/unconfined subroutine)
       qq(1:nP,1) = kappasq(1:nP)*lin%sf**2/4.0_DP  !! outside q
       qq(1:nP,2) = cmplx(-huge(1.),huge(1.),DP)    !! invalid
    end if
    
    if(lin%first) then
       if(.not. allocated(mf(1,1)%A)) then
          do i=1,nP
             write(*,666) 'qq(',i,',1) = (',real(qq(i,1)),',',&
                  & aimag(qq(i,1)),'), M=',lin%MS,ceiling(abs(2.0*qq(i,1))), &
                  & max(lin%MS,16,ceiling(abs(2.0*qq(i,1))))
             ! allocate matrices and compute mcn and A & B matrices
             mf(i,1) = mathieu_init(qq(i,1), &
                  & MM=max(lin%MS,16,ceiling(abs(2.0*qq(i,1)))),CUTOFF=1.0D-9)
             if (lin%match) then
                write(*,666) 'qq(',i,',2) = (',real(qq(i,2)),',',&
                     & aimag(qq(i,2)),'), M=',lin%MS,ceiling(abs(2.0*qq(i,2))),&
                     & max(lin%MS,16,ceiling(abs(2.0*qq(i,2))))
                mf(i,2) = mathieu_init(qq(i,2),  &
                     & MM=max(lin%MS,16,ceiling(abs(2.0*qq(i,2)))),CUTOFF=1.0D-9)
             end if
          end do
          if (lin%match) then
             MS = maxval(mf(:,1:2)%M)
          else
             MS = maxval(mf(:,1:1)%M)
          end if
       end if
666    format(A,I0,2(A,ES10.3E2),A,3(I0,1X))
       print *, 'c'
       if(lin%match) then
          if(.not. allocated(Am)) then
             allocate(Am(1:2*M,0:4*N+1),bm(1:2*M),stat=ierr)
             if (ierr /= 0) then
                write(*,'(A,2(I0,1X))') 'error allocating: Am,bm ',2*M,4*N+1
                stop 102
             end if
          end if
          if(.not. allocated(x)) then
             allocate(x(0:4*N+1,1:nP),stat=ierr)
             if (ierr /= 0) then
                write(*,'(A,2(I0,1X))') 'error allocating: x ',4*N,np
                stop 103
             end if
          end if
       elseif(lin%flux) then
          write(*,'(A,ES12.5,A)') 'specified flux element; q=',lin%q
       else
           write(*,'(A,ES12.5,A)') 'specified head element; h=',lin%q
          if(.not. allocated(Am)) then
             allocate(Am(1:M,0:2*N),bm(1:M),x(0:4*N+1,1:nP),stat=ierr)
             if (ierr /= 0) then
                write(*,'(A,3(I0,1X))') 'error allocating Am,bm,x',M,2*N,nP
                stop 104
             end if
          end if
       end if
       
       if(.not. allocated(psiM)) then
          allocate(psiM(1:M), rw(1:M), Zw(1:M), Rcw(1:M), vi(0:N),stat=ierr)
          if (ierr /= 0) then
             write(*,'(A,2(I0,1X))') 'error allocating: psiM,rw,Zw,Rcw,vi ',M,N
             stop 105
          end if
       end if
       vi(0:N) = [(i, i=0,N)]

       ! matching points along circumference of ellipse
       psiM(1:M) = linspace(-PI,+PI,M)
       if(.not. allocated(sig)) then
          allocate(sig(0:MS-1),denom(0:MS-1),stat=ierr)
          if (ierr /= 0) then
             write(*,'(A,I0)') 'error allocating: sig,denom ',MS-1
             stop 106
          end if
       end if

       ! related to the infinite matrices used in MF calculation
       sig(0:MS-1) =   real([((-1)**i, i=0,MS-1)],DP)
       denom(0:MS-1) = real([ (1 - (2*i)**2 ,i=0,MS-1) ],DP)

       ! for checking size of q vs size of matrix
    else
       !! variables re-used in later calls, but only one term
       allocate(rw(1),Rcw(1),vi(0:N),stat=ierr)
       if (ierr /= 0) then
          write(*,'(A,(I0,1X))') 'error allocating: rw,Rcw,vi ',N
          stop 107
       end if
       vi(0:N) = [(i, i=0,N)]
    end if

    do i = 1,np
!!$       write(*,'(I0,1X)',advance='no') i
       if (lin%first) then
          
          write(*,'(A)',advance='no') '*';

          !! specified flux element
          if (.not. lin%match .and. lin%flux) then
          end if
          
          if(lin%match .or. (.not. lin%flux))then
             Am(:,:) = cmplx(-huge(1.),+huge(1.),DP)
             bm(:) =   cmplx(+huge(1.),-huge(1.),DP)
          end if

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          if ((.not. lin%flux) .and. (.not. lin%match)) then ! specified head element

             ! normalizing terms cancel out radial Mathieu functions
             Am(1:M,0:N) =     transpose(ce(mf(i,1),vi(0:N),psiM(1:M)))/lin%k ! a_0:n
             Am(1:M,N+1:2*N) = transpose(se(mf(i,1),vi(1:N),psiM(1:M)))/lin%k ! b_1:n

             ! RHS vector for constant head line
             bm(1:M) = lin%q/p(i)

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
          elseif (lin%match) then ! matching element with a well in the background

             Am(1:M,0:N) =     transpose(ce(mf(i,1),vi(0:N),psiM(1:M)))/lin%k ! a_n
             Am(1:M,N+1:2*N) = transpose(se(mf(i,1),vi(1:N),psiM(1:M)))/lin%k ! b_n
             Am(1:M,2*N+1:3*N+1) = transpose(-ce(mf(i,2),vi(0:N),psiM(1:M)))/lin%Ek ! c_n
             Am(1:M,3*N+2:4*N+1) = transpose(-se(mf(i,2),vi(1:N),psiM(1:M)))/lin%Ek ! d_n

             ! flux matching
             Am(M+1:2*M,0:N) =     transpose(ce(mf(i,1),vi(0:N),psiM(1:M)))* &
                  & spread(DKe(mf(i,1),vi(0:N),lin%eta0)/Ke(mf(i,1),vi(0:N),lin%eta0),1,M) ! a_n
             Am(M+1:2*M,N+1:2*N) = transpose(se(mf(i,1),vi(1:N),psiM(1:M)))* &
                  & spread(DKo(mf(i,1),vi(1:N),lin%eta0)/Ko(mf(i,1),vi(1:N),lin%eta0),1,M) ! b_n

             Am(M+1:2*M,2*N+1:3*N+1) = transpose(-ce(mf(i,2),vi(0:N),psiM(1:M)))* &
                  & spread(DIe(mf(i,2),vi(0:N),lin%eta0)/Ie(mf(i,2),vi(0:N),lin%eta0),1,M) ! c_n
             Am(M+1:2*M,3*N+2:4*N+1) = transpose(-se(mf(i,2),vi(1:N),psiM(1:M)))* &
                  & spread(DIo(mf(i,2),vi(1:N),lin%eta0)/Io(mf(i,2),vi(1:N),lin%eta0),1,M) ! d_n

             ! RHS vector: effect of well at 1,1
             !-------------------------------------------
             ! cartesian coordinates of points along circumference of ellipse
             Zw(1:M) = lin%sf*cacosh(cmplx(spread(lin%eta0,1,M),psiM(1:M),DP))

             open(unit=99,file='ellipse.bdry')
             do j=1,M
                write(99,*) real(Zw(j)),aimag(Zw(j))
             end do
             write(99,*)  real(Zw(1)),aimag(Zw(1)) !! connect it back up
             close(99)

             ! cartesian components of vector from well to pts on bndry of ellipse
             Rcw(1:M) = Zw(1:M) - lin%well

             ! distance from well to pt on circ of ellipse
             rw(1:M) = abs(Rcw(1:M))

             ! q associated with Bessel functions
             qb = sqrt(p(i)/lin%alpha)

             do j=1,M
                call cbesk(rw(j)*qb,0.0_DP,1,2,kw(0:1),nz,ierr)
                if (nz + ierr /= 0) write(*,'(A,(I0,1X))') 'bessel function error ',ierr

                !! well fcn is pos, becomes neg when moved to RHS
                phiWell = lin%q*kw(0)/(TWOPI*p(i)*lin%k)

                !! deriv of well fcn is neg, becomes pos when moved to RHS
                phiDr =  -lin%q*qb*kw(1)/(TWOPI*p(i)*lin%k)

                !! area source is negative, stays neg on RHS
                bm(j) = -phiWell/lin%k - lin%gamma*lin%Ealpha/(p(i)*lin%Ek)

                ! flux effects
                bm(M+j) = -phiDr* &
                     & lin%sf*(real(Rcw(j))/rw(j)*sinh(lin%eta0)*cos(psiM(j)) + &
                     & aimag(Rcw(j))/rw(j)*cosh(lin%eta0)*sin(psiM(j)))
             end do
             
          end if

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          ! solve for coefficients in specified head case

          if ((.not. lin%match) .and. (.not. lin%flux))  then

             if(.not. allocated(rwork)) then
                allocate(rwork(5*min(M,2*N+1)),s(min(M,2*N+1)),stat=ierr)
                if (ierr /= 0) then
                   write(*,*) 'error allocating: rwork,s ',M,N
                   stop 108
                end if
             end if
             
             call zgelss(M,2*N+1,1,Am(1:M,0:2*N),M,bm(1:M),M,s,rcond,rank, &
                  & work, LWORK, RWORK, info)

             if(info /= 0) write(*,'(A,(I0,1X))') 'ZGELSS: error ', info
             ! copy results into x
              x(0:2*N,i) = bm(1:2*N+1)

          elseif (lin%match) then
             if(.not. allocated(rwork)) then
                allocate(rwork(2*max(2*M,4*N+2)+1),s(min(2*M,4*N+2)),stat=ierr)
                if (ierr /= 0) then
                   write(*,'(A,2(I0,1X))') 'error allocating: rwork,s ',M,N
                end if
             end if
             
             call zgelss(2*M,4*N+2,1,Am(1:2*M,0:4*N+1),2*M,bm(1:2*M),2*M,s,&
                  & rcond,rank, work, LWORK, RWORK, info)

             if(info /= 0) write(*,'(A,(I0,1X))') 'ZGELSS: error ', info
             ! copy results into x
              x(0:4*N+1,i) = bm(1:4*N+2)
           end if
        end if

       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       ! use x(:,i) to calculate results
       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       
       phiDr = cmplx(0.0,0.0,DP)
       phiDeta = cmplx(0.0,0.0,DP)
       phiDpsi = cmplx(0.0,0.0,DP)

       if (lin%calc_head .and. eta>=lin%eta0 .and. &
            & (lin%match .or. (.not.lin%flux))) then
          ! head effects of a matching, or specified head element OUTSIDE

          F(i) = sum(x(0:N,i)* ce(mf(i,1),vi(0:N),psi)* &
               &   Ke(mf(i,1),vi(0:N),eta)/ Ke(mf(i,1),vi(0:N),lin%eta0)) + &
               & sum(x(N+1:2*N,i)* se(mf(i,1),vi(1:N),psi)* &
               &   Ko(mf(i,1),vi(1:N),eta)/ Ko(mf(i,1),vi(1:N),lin%eta0))

          ! add effect of pumping well at 1,1
          if (lin%match) then

             rw1 = abs(calcpt - lin%well)
             qb = sqrt(p(i)/lin%alpha)

             call cbesk(rw1*qb, 0.0_DP,1,1,kw(0),nz,ierr)
             if (nz + ierr /= 0) write(*,'(A,(I0,1X))') 'bessel function K0 error ',ierr

             F(i) = F(i) + lin%q*kw(0)/(TWOPI*p(i)*lin%k)
          end if

          !! convert solution to head
          F(i) = F(i)/lin%k  

       elseif ((.not. lin%calc_head) .and. eta>=lin%eta0 .and. &
            & (lin%match .or. (.not. lin%flux))) then
          ! normal flux effects of a matching or head element OUTSIDE

          phiDeta = sum(x(0:N,i)*ce(mf(i,1),vi(0:N),psi)* &
               & DKe(mf(i,1),vi(0:N),eta)/Ke(mf(i,1),vi(0:N),lin%eta0)) + &
               & sum(x(N+1:2*N,i)*se(mf(i,1),vi(1:N),psi)* &
               & DKo(mf(i,1),vi(1:N),eta)/Ko(mf(i,1),vi(1:N),lin%eta0))

          phiDpsi = sum(x(0:N,i)*Dce(mf(i,1),vi(0:N),psi)* &
               & Ke(mf(i,1),vi(0:N),eta)/Ke(mf(i,1),vi(0:N),lin%eta0)) + &
               & sum(x(N+1:2*N,i)*Dse(mf(i,1),vi(1:N),psi)* &
               & Ko(mf(i,1),vi(1:N),eta)/Ko(mf(i,1),vi(1:N),lin%eta0))

          if (lin%match) then
             rw1 = abs(calcpt - lin%well)
             qb = sqrt(p(i)/lin%alpha)

             call cbesk(rw1*qb,1.0_DP,1,1,kw(1),nz,ierr)
             if (nz + ierr /= 0) write(*,'(A,(I0,1X))') 'bessel function K1 error ',ierr
             phiDr = -lin%q*qb*kw(1)/(TWOPI*p(i)*lin%k)
          end if

       elseif (lin%match .and. lin%calc_head .and. eta<lin%eta0) then
          ! head effects of a matching element INSIDE

          F(i) = sum(x(2*N+1:3*N+1,i)*ce(mf(i,2),vi(0:N),psi)* &
               & Ie(mf(i,2),vi(0:N),eta)/Ie(mf(i,2),vi(0:N),lin%eta0)) + &
               & sum(x(3*N+2:4*N+1,i)*se(mf(i,2),vi(1:N),psi)* &
               & Ie(mf(i,2),vi(1:N),eta)/Ie(mf(i,2),vi(1:N),lin%eta0))
 
          !! area source effects & DisPot -> head
          F(i) = (F(i) - lin%Ealpha*lin%gamma/p(i))/lin%Ek

       elseif (lin%match .and. (.not.lin%calc_head) .and. eta<lin%eta0) then
          ! normal flux effects of a matching element INSIDE

          phiDeta = sum(x(2*N+1:3*N+1,i)*ce(mf(i,2),vi(0:N),psi)* &
               & DIe(mf(i,2),vi(0:N),eta)/Ie(mf(i,2),vi(0:N),lin%eta0)) + &
               & sum(x(3*N+2:4*N+1,i)*se(mf(i,2),vi(1:N),psi)* &
               & DIo(mf(i,2),vi(1:N),eta)/Io(mf(i,2),vi(1:N),lin%eta0))

          phiDpsi = sum(x(2*N+1:3*N+1,i)*Dce(mf(i,2),vi(0:N),psi)* &
               & Ie(mf(i,2),vi(0:N),eta)/Ie(mf(i,2),vi(0:N),lin%eta0)) + &
               & sum(x(3*N+2:4*N+1,i)*Dse(mf(i,2),vi(1:N),psi)* &
               & Io(mf(i,2),vi(1:N),eta)/Io(mf(i,2),vi(1:N),lin%eta0))

       elseif((.not. lin%match) .and. lin%flux .and. lin%calc_head) then
          ! head effects for a specified flux line element          
          nmax = floor(N/2.0)-1

          !! removed factor of 4 from numerator
          !! apply time behavior (constant) and convert DP -> h in denominator

          F(i) = -lin%Q/(lin%k*p(i)*PI)*sum(sig(0:nmax)*sum(&
               & spread(sig(0:mf(i,1)%M-1)/denom(0:mf(i,1)%M-1),2,nmax+1)* &
               & conjg(mf(i,1)%A(:,0:nmax,0)),dim=1) * ce(mf(i,1),2*vi(0:nmax),psi)* &
               & Ke(mf(i,1),2*vi(0:nmax),eta)/DKe(mf(i,1),2*vi(0:nmax),0.0_DP))

       elseif((.not. lin%match) .and. lin%flux .and. (.not.lin%calc_head)) then
          ! flux effects for a specified flux line element
          nmax = floor(N/2.0)-1

          ! time behavior (constant) in denominator
          ! changed 4 in numerator to 2
          phiDeta = lin%Q*2.0_DP/(PI*p(i))*sum(sig(0:nmax)*sum(&
               & spread(sig(0:mf(i,1)%M-1)/denom(0:mf(i,1)%M-1),2,nmax+1)* &
               & conjg(mf(i,1)%A(:,0:nmax,0)),dim=1) * ce(mf(i,1),2*vi(0:nmax),psi)* &
               & DKe(mf(i,1),2*vi(0:nmax),eta)/DKe(mf(i,1),2*vi(0:nmax),0.0_DP))
       end if

       if (.not. lin%calc_head) then
          !! project elliptical and polar onto cartesian x or y
          metricsq = lin%sf**2*(cosh(2.0_DP*eta) - cos(2.0_DP*psi))
          
          if (lin%cartesian) then
             if (lin%calc_velx) then
                !! x-flux
                F(i) = lin%sf/metricsq*(&
                     & phiDeta*sinh(eta)*cos(psi) - phiDpsi*cosh(eta)*sin(psi)) + &
                     & phiDr*real(calcpt-lin%well)/abs(calcpt-lin%well)
                
             else !! y-flux
                F(i) = lin%sf/metricsq*(&
                     & phiDeta*cosh(eta)*sin(psi) + phiDpsi*sinh(eta)*cos(psi)) + &
                  & phiDr*aimag(calcpt-lin%well)/abs(calcpt-lin%well)
             end if
          else
             F(i) = phiDeta  !! just return elliptic normal flux
          end if
       end if
    
!!$       if (i==nP) write(*,'(/)',advance='no')
    end do
!!$    if (allocated(work)) deallocate(work)
  end function mathieulinehead
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! an approximation to the constant line source using a sum of point sources
  function pointhead(p,calcX,calcY,lin,pts) result(PH)
    use constants, only : DP, CZERO, TWOPI, PI
    use complex_bessel, only : cbesk
    use utility, only : outerprod
    implicit none

    complex(DP), intent(in), dimension(:) :: p
    real(DP), intent(in) :: calcX,calcY
    type(line_external), intent(in) :: lin
    complex(DP), dimension(size(p)) :: PH
    complex(DP), dimension(size(p)) :: q ! not the Mathieu function q
    integer, intent(in) :: pts
    complex(DP), dimension(size(p),pts) :: k0, k1
    real(DP), dimension(1:pts) :: r,xi
    real(DP) :: begin, delxi
    integer :: i, j, np, nz, ierr

    if (pts <= 0) then
       write(*,'(A,I0)') '# points for point approximation to line must be >0 ',pts
       stop 222
    end if
    
    np = size(p)
    q(1:np) = sqrt(p(1:np)/lin%alpha)

    begin = -lin%sf   !! start at left end
    delxi = 2.0_DP*lin%sf/real(PTS-1,DP)  !! step size

    PH = CZERO

    !! xi is 'arc length' in line integral along line source
    xi(1:pts) = begin + delxi*[(real(i-1,DP),i=1,pts)]
    r(1:pts) = abs(cmplx(calcX-xi(1:pts),calcY,DP))

    do i=1,pts
       !! compute required bessel functions for each value of p
       do j=1, np
          if (lin%calc_head) then
             call cbesk(r(i)*q(j),0.0_DP,1,1,k0(j,i),nz,ierr)
             if (nz /=0 .or. ierr /= 0) then
                write(*,'(A,2(I0,1X))') 'K0 bessel function error ', nz, ierr
             end if
          else
             call cbesk(r(i)*q(j),1.0_DP,1,1,k1(j,i),nz,ierr)
             if (nz /=0 .or. ierr /= 0) then
                write(*,'(A,2(I0,1X))') 'K1 bessel function error ', nz, ierr
             end if
          end if
       end do
    end do
    
    !! is everything off by a sign here?

    if (lin%calc_head) then ! head
       PH(1:np) = sum(k0(1:np,1:pts),dim=2)*&
            & lin%q/(PI*real(PTS,DP)*p(1:np)*lin%k)
       
    elseif (lin%calc_velx) then ! x flux
       PH(1:np) = sum(k1(1:np,1:pts)*&
            & outerprod(lin%q*q(1:np)/(PI*real(PTS,DP)*p(1:np)*lin%k),&
            & (calcX-xi(1:pts))/r(1:pts)),dim=2)
       
    else ! y flux
       PH(1:np) = sum(k1(1:np,1:pts)*&
            & outerprod(q(1:np)*lin%q/(PI*real(PTS,DP)*p(1:np)*lin%k),&
            & calcY/r(1:pts)),dim=2)
    end if

  end function pointhead

  !! this now produces kappa^2, rather than kappa, because that is more useful
  function compute_ELleaky_qv(p,lin) result(q)
    use constants, only : DP, PI
    integer, parameter :: NTERMS = 300, MAXITER = 2000
    complex(DP), intent(in), dimension(:) :: p
    type(line_external), intent(in) :: lin
    complex(DP), dimension(size(p)) :: q  !! actually q (or kappa) SQUARED

    integer :: np
    complex(DP), dimension(size(p)) :: kap2
    complex(DP), dimension(size(p)) :: exp2z

    real(DP) :: sigma, alpha, alpha2
    real(DP), dimension(NTERMS) :: guess, root, gamma
    complex(DP), dimension(size(p)) :: kernel
    real(DP) :: x, delta
    integer :: k, kk

    np = size(p)
    sigma = sqrt(lin%Ss/lin%Sy)
    alpha = lin%k/lin%Ss
    alpha2 = lin%k2/lin%Ss2

    !! leaky-ness
    !! ##############################
    if(lin%leak == 0) then
       !! no leaky layer, standard definition
       q(1:np) = p(1:np)/alpha
       print *, 'non-leaky'
    else
       print *, 'leaky type',lin%leak
       kap2(1:np) = sqrt(p(:)/alpha2)
       exp2z(1:np) = exp(-2.0_DP*kap2(:)*lin%b2)

       if(lin%leak == 1) then
          !! case I, no-drawdown condition at top of aquitard
          q(:) = p(:)/alpha + kap2(:)*lin%k2/(lin%b*lin%k)*&
               & (1.0_DP + exp2z(:))/(1.0_DP - exp2z(:))
       elseif(lin%leak == 2) then
          !! case II, no-flow condition at top of aquitard
          q(:) = p(:)/alpha + kap2(:)*lin%k2/(lin%b*lin%k)*&
               & (1.0_DP - exp2z(:))/(1.0_DP + exp2z(:))
       elseif(lin%leak == 3) then
          !! aquitard thickness -> infinity
          q(:) = p(:)/alpha + kap2(:)*lin%k2/(lin%b*lin%k)
       else
          stop 'ERROR: incorrect value for leaky flag'
       end if
    end if

    !! unconfined-ness 
    !! ##############################
    if( .not. lin%unconfined ) then
       print *, 'confined'
       !! do nothing, q already computed above
    else
       print *, 'unconfined'
       !! Boulton unconfined source (Herrera infinite sum Kernel)
       !! guess is halfway between asymptotes of cot()

       !! should be able to solve for Neuman solution, integrated across the thickness...

       !! kernel is not a function of p
       guess(2:NTERMS) = PI*(real([(k, k=1,NTERMS-1)],DP) + 0.5_DP)/sigma
       guess(1) = 1.7D0

       !! first root is hard to find with NR, 
       !! use TS approximation for tangent and re-arrange
       x = guess(1)
       NR1: do kk = 1,MAXITER
          delta = (x + (sigma - 1.0_DP/sigma)*(x*sigma + (x*sigma)**3/3.0_DP + &
               & 2.0_DP*(sigma*x)**5/15.0_DP) + 17.0_DP*(x*sigma)**7/315.0_DP)/ &
               & (1.0_DP - (1.0_DP/sigma - sigma)*(sigma + x**2*sigma**3 + &
               & 2.0_DP*x**4*sigma**5/3.0_DP + 17.0_DP*x**6*sigma**7/45.0_DP))
          x = x - delta
          if (abs(delta) <= 1.0D-10) then
             root(1) = x
             exit NR1
          end if
          if(kk == MAXITER) write(*,'(A)') 'failed to converge for first root'
       end do NR1

       do k = 2, NTERMS
          x = guess(k)
          NR: do kk = 1,MAXITER
             delta = (1.0_DP/tan(x*sigma) + (sigma - 1.0_DP/sigma)/x)/&
                  & (sigma/(sin(sigma*x)**2) + (sigma + 1.0_DP/sigma)/x**2)
             x = x + delta
             if (abs(delta) <= spacing(x)*10.0) then
                root(k) = x
                exit NR
             end if
             if(kk == MAXITER) write(*,'(I0,A)') k,'th root failed to converge'
          end do NR
       end do

       gamma(1:NTERMS) = lin%Kz*root(1:NTERMS)**2/(lin%b*lin%Sy)
       kernel(1:np) = 2.0_DP*sum(spread(gamma(1:NTERMS),2,np)/&
            & ((spread(root(1:NTERMS)**2,2,np) - 1.0_DP + sigma**2)* &
            & (spread(p(1:np),1,NTERMS) + spread(gamma(1:NTERMS),2,np))),dim=1)

       q(1:np) = q(1:np) + p(1:np)*lin%Sy/lin%k*kernel
    end if

  end function compute_ELleaky_qv
  
end module mathieu_line
