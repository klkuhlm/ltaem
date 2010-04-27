! $Id: arc_test.f90,v 1.8 2007/07/17 23:16:25 kris Exp kris $
module arc_test
  implicit none

  private
  public :: arc_head, pointhead

contains

  function arc_head(p,calcX,calcY,alpha, arcR1, arcTH1, arcTH2, calc_head, calc_velx, &
       & Wellq,arcQ,kBG, passive, arcX1, arcY1) result(F)
    use complex_bessel, only : cbesk, cbesi
    use constants, only : DP, TWOPI, PI, CZERO, RONE

!!$    character(3) :: chnp

    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p)) :: F

    complex(DP) :: well, calcpt, arcCent, dpotra, dpotrw

    !! things which were in calcshareddata module
    logical, intent(in) :: passive
    ! parameters of the line source
    real(DP), intent(in) :: alpha, kBG, arcR1, arcTH1, arcTH2, arcX1, arcY1,WellQ,arcQ
    logical, intent(in) :: calc_head, calc_velx
    real(DP), intent(in) :: calcX, calcY

    real(DP) :: r, theta
    complex(DP), dimension(0:1) :: kw
    complex(DP), dimension(size(p)) :: kappa
    integer :: nP, i, j, nz, ierr !! ,k

    integer, save :: N, M
    real(DP), save, allocatable :: vi(:)

    ! coefficient matrix used in least-squares process for finding x (strengths of gFc)
    complex(DP), allocatable :: Am(:,:), bm(:)
    complex(DP), save, allocatable :: x(:,:)

    complex(DP), allocatable :: besk(:), Dbesk(:), besk0(:), &
         & besi(:), Dbesi(:),besi0(:)

    complex(DP), allocatable :: Zw(:), Rcw(:), dpotrwv(:) !! Pw(:), 
    real(DP), allocatable :: rw(:) ! distance to well at 1,1

    ! is this the first time through this subroutine?
    logical, save :: first = .true. 
    real(DP), allocatable, save :: phiM(:,:) ! vector of matching locations

    well = cmplx(1.0_DP,1.0_DP,DP)
    nP = size(p)
    arcCent = cmplx(arcX1,arcY1,DP)
    kappa(1:np) = sqrt(p(1:np)/alpha)
        
    if(first) then
       N = 30
       M = (2*N+1)
       if(.not. passive) then
          print *, 'active arc element with well at (1,1) of strength:', Wellq
          print *, 'arc radius:',arcR1, ' angular range:',arcTH1,arcTH2
          print *, 'background K:',kBG
          print *, 'background alpha:', alpha
       else
          print *, 'passive specified flux arc element; q=',ARCQ
          print *, 'arc radius:',arcR1, ' angular range:',arcTH1,arcTH2
          print *, 'background K:',kBG
          print *, 'background alpha:', alpha
       end if
       allocate(Am(1:3*M,0:4*N+1),bm(1:3*M),x(0:4*N+1,1:np))
       print *, 'N:',N,' M: ',M
       
       allocate(phiM(2,1:M), rw(1:2*M), Zw(1:2*M), &
            & Rcw(1:2*M), vi(0:N), dpotrwv(1:2*M)) !! , Pw(1:2*M)
       vi(0:N) = real((/(i, i=0,N)/),DP)

       ! matching points along arc (crowded near ends)
       ! M points in each interval (th1 -> th2 is interval 1, 
       ! the rest is second interval
!!$       phiM(1:2,1:M) = chebspace(arcTH1,arcTH2,M)
       phiM(1,1:M) = linspace(arcTH1,arcTH2,M)
       phiM(2,1:M) = linspace_no_ends(arcTH2,arcTH2+(arcTH2 - arcTH1),M)
       where(phiM(2,1:M) > PI)
          phiM(2,1:M) = phiM(2,1:M) - TWOPI
       end where
       

    else
       !! variables re-used in later calls, but only one term
       allocate(rw(1),Rcw(1)) !!,Pw(1)
    end if

    allocate(besk(0:N), Dbesk(0:N), besi(0:N), Dbesi(0:N), besk0(0:N),besi0(0:N))
    
    do i = 1,np
       if (first) then
          write(*,'(A)',advance='no') '*'

          !! there are three 'panels' to the Am matrix 
          !! for an active arc element
          
          !! 1) flux matching associated with specified flux portion
          !! 2) head and 3) flux associated with matching portion
          
          !! same aquifer properties 'inside' and 'outside' the arc

          call cbesk(arcR1*kappa(i),0.0_DP,1,N+1,besk0(0:N),nz,ierr)
          if(nz + ierr /= 0) print *, 'K bessel function error'
          call cbesi(arcR1*kappa(i),0.0_DP,1,N+1,besi0(0:N),nz,ierr)
          if(nz + ierr /= 0) print *, 'I bessel function error'

          Dbesk(0:N) = besk_deriv(besk0(0:N),arcR1*kappa(i),kappa(i))
          Dbesi(0:N) = besi_deriv(besi0(0:N),arcR1*kappa(i),kappa(i))

          ! a_n: (outside on specified flux portion)
          Am(1:M,0:N) = cos(outerprod(phiM(1,1:M),vi(0:N)))*&
               & spread(Dbesk(0:N)/besk0(0:N),1,M)

          ! b_n: 
          Am(1:M,N+1:2*N) = sin(outerprod(phiM(1,1:M),vi(1:N)))*&
               & spread(Dbesk(1:N)/besk0(1:N),1,M)

          ! c_n: (inside)
          Am(1:M,2*N+1:3*N+1) = -cos(outerprod(phiM(1,1:M),vi(0:N)))*&
               & spread(Dbesi(0:N)/besi0(0:N),1,M)

          ! d_n: 
          Am(1:M,3*N+2:4*N+1) = -sin(outerprod(phiM(1,1:M),vi(1:N)))*&
               & spread(Dbesi(1:N)/besi0(1:N),1,M)

          !##############################
          ! a_n: (outside on head matching portion)
          Am(M+1:2*M,0:N) = cos(outerprod(phiM(2,1:M),vi(0:N)))

          ! b_n: 
          Am(M+1:2*M,N+1:2*N) = sin(outerprod(phiM(2,1:M),vi(1:N)))

          ! c_n: (inside)
          Am(M+1:2*M,2*N+1:3*N+1) = -cos(outerprod(phiM(2,1:M),vi(0:N)))

          ! d_n: 
          Am(M+1:2*M,3*N+2:4*N+1) = -sin(outerprod(phiM(2,1:M),vi(1:N)))

          !##############################
          ! a_n: (outside on normal flux matching portion)
          Am(2*M+1:3*M,0:N) = cos(outerprod(phiM(2,1:M),vi(0:N)))*&
               & spread(Dbesk(0:N)/besk0(0:N),1,M)

          ! b_n: 
          Am(2*M+1:3*M,N+1:2*N) = sin(outerprod(phiM(2,1:M),vi(1:N)))*&
               & spread(Dbesk(1:N)/besk0(1:N),1,M)

          ! c_n: (inside)
          Am(2*M+1:3*M,2*N+1:3*N+1) = -cos(outerprod(phiM(2,1:M),vi(0:N)))*&
               & spread(Dbesi(0:N)/besi0(0:N),1,M)

          ! d_n: 
          Am(2*M+1:3*M,3*N+2:4*N+1) = -sin(outerprod(phiM(2,1:M),vi(1:N)))*&
               & spread(Dbesi(1:N)/besi0(1:N),1,M)       

          if(.not. passive) then
             
             ! RHS vector: effect of well at 1,1
             !-------------------------------------------

             ! cartesian coordinates of matching points along arc
             Zw(1:M) =     arcr1*cmplx(cos(phiM(1,1:M)),sin(phiM(1,1:M))) + arcCent
             Zw(M+1:2*M) = arcr1*cmplx(cos(phiM(2,1:M)),sin(phiM(2,1:M))) + arcCent

             open(unit=99,file='arc.bdry',action='write',status='replace')
             write(99,'(A)') '# active no-flux boundary portion of circular boundary'
             do j=1,M  
                write(99,*) real(Zw(j)),aimag(Zw(j))
             end do
             write(99,'(//A)') '# matching portion of circular boundary'
             do j=1,M
                write(99,*) real(Zw(M+j)),aimag(Zw(M+j))
             end do
             close(99)
       
             ! cartesian components of vector to well from pt on circ. of ellipse
             Rcw(1:2*M) = well - Zw(1:2*M)

             ! distance from well to pt on arc
             rw(1:2*M) = abs(Rcw(1:2*M))

!!$             ! angle associated with radial vector pointing from arc
!!$             Pw(1:2*M) = atan2(aimag(Rcw(1:2*M)),real(Rcw(1:2*M)))

             !! flux effects on specified flux arc
             do j=1,M
                call cbesk(rw(j)*kappa(i),1.0_DP,1,1,kw(1),nz,ierr)
                if (nz + ierr /= 0) print *, 'Well bessel function error'

                dpotrwv(j) = WellQ*kappa(i)*kw(1)/(TWOPI*p(i)*Kbg)
             end do

             do j=M+1,2*M
                call cbesk(rw(j)*kappa(i),0.0_DP,1,2,kw(0:1),nz,ierr)
                if (nz + ierr /= 0) print *, 'Well bessel function error'

                !! head effects of well on opening
                bm(j) =  -WellQ*kw(0)/(TWOPI*p(i)*Kbg)
                !! flux effects of well on opening
                dpotrwv(j) = WellQ*kappa(i)*kw(1)/(TWOPI*p(i)*Kbg)
             end do

             !! project flux from one radial coordinate to another
             bm(1:M) = dpotrwv(1:M) * &
                  &   (real(Rcw(1:M))/rw(1:M)*cos(phiM(1,1:M)) &
                  & + aimag(Rcw(1:M))/rw(1:M)*sin(phiM(1,1:M)))
             bm(2*M+1:3*M) = dpotrwv(M+1:2*M) * &
                  &   (real(Rcw(M+1:2*M))/rw(M+1:2*M)*cos(phiM(2,1:M)) &
                  & + aimag(Rcw(M+1:2*M))/rw(M+1:2*M)*sin(phiM(2,1:M)))

          else !! passive
             bm(1:M) = arcQ/(arcr1*abs(arcTH1 - arcTH2)*p(i))
             bm(M+1:3*M) = cmplx(0.0,0.0,DP)
          end if

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          !! solve for coefficients, active element

          x(0:4*N+1,i) = matmul(inverse(matmul(conjg(transpose(&
               & Am(1:3*M,0:4*N+1))), Am(1:3*M,0:4*N+1))),&
               & matmul(conjg(transpose(Am(1:3*M,0:4*N+1))),bm(1:3*M)))
       endif


       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       ! use x(:,i) to calculate results
       !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

       calcpt = cmplx(calcX,calcY,DP)
       r = abs(calcpt - arccent)
       theta = atan2(aimag(calcpt - arccent),real(calcpt - arccent))
       dpotrw = cmplx(0.0,0.0,DP)

       !! consolidate BF calls, since same K/Ss everywhere
       if(r >= arcR1) then
          call cbesk(r*kappa(i),0.0_DP,1,N+1,besk(0:N),nz,ierr)
          call cbesk(arcR1*kappa(i),0.0_DP,1,N+1,besk0(0:N),nz,ierr)
          if(nz + ierr /= 0) print *, 'K bessel function error'          
       else
          call cbesi(r*kappa(i),0.0_DP,1,N+1,besi(0:N),nz,ierr)
          call cbesi(arcR1*kappa(i),0.0_DP,1,N+1,besi0(0:N),nz,ierr)
          if(nz + ierr /= 0) print *, 'I bessel function error'
       end if
       
       if (calc_head) then
          if(r >= arcR1) then

             ! head effects outside specified flux arc element
             F(i) = besk(0)/besk0(0)*x(0,i)
             F(i) = F(i) + sum(besk(1:N)/besk0(1:N)*&
                  & (x(1:N,i)*    cos(vi(1:N)*theta) + &
                  &  x(N+1:2*N,i)*sin(vi(1:N)*theta)))
             
             if(.not. passive) then
                !! include head effects of well
                rw(1) = abs(calcpt - well)
                
                call cbesk(rw(1)*kappa(i),0.0_DP,1,1,kw(0),nz,ierr)
                if (nz + ierr /= 0) print *, 'well bessel function error.'
                
                F(i) = F(i) + WellQ*kw(0)/(TWOPI*p(i)*Kbg)
             end if
          
             F(i) = F(i)/kbg  !! return head, not discharge potential

          else

             !! inside effects  specified flux arc element
             F(i) = besi(0)/besi0(0)*x(2*N+1,i)
             F(i) = F(i) + sum(besi(1:N)/besi0(1:N)*&
                  & (x(2*N+2:3*N+1,i)*cos(vi(1:N)*theta) + &
                  &  x(3*N+2:4*N+1,i)*sin(vi(1:N)*theta)))
             
             F(i) = F(i)/kbg
             !! no background or particular solutions inside arc

          end if
       else  !! calculat flux effects
          if(r >= arcR1) then
             Dbesk(0:N) = besk_deriv(besk(0:N),r*kappa(i),kappa(i))

             ! normal flux outside
             dpotra = Dbesk(0)/besk0(0)*x(0,i)
             dpotra = dpotra + sum(Dbesk(1:N)/besk0(1:N)*&
                  & (x(1:N,i)*    cos(vi(1:N)*theta) + &
                  &  x(N+1:2*N,i)*sin(vi(1:N)*theta)))

             if(.not. passive) then
                !! flux effects of well
                rw(1) = abs(calcpt - well)
!!$                pw(1) = atan2(aimag(calcpt-well),real(calcpt-well))
                
                call cbesk(rw(1)*kappa(i),1.0_DP,1,1,kw(1),nz,ierr)
                if (nz + ierr /= 0) print *, 'well bessel function error.'
                
                dpotrw = -WellQ*kw(1)*kappa(i)/(TWOPI*p(i)*Kbg)
             end if

          else
             Dbesi(0:N) = besi_deriv(besi(0:N),r*kappa(i),kappa(i))
             
             !! normal flux inside
             dpotra = Dbesi(0)/besi0(0)*x(2*N+1,i)
             dpotra = dpotra + sum(Dbesi(1:N)/besi0(1:N)*&
                  & (x(2*N+2:3*N+1,i)*cos(vi(1:N)*theta) + &
                  &  x(3*N+2:4*N+1,i)*sin(vi(1:N)*theta)))
          end if
          
          if(calc_velx) then
             F(i) = dpotra*real(calcpt - arccent)/r + &
                  & dpotrw*real(calcpt - well)/rw(1)
          else
             F(i) = dpotra*aimag(calcpt - arccent)/r + &
                  & dpotrw*aimag(calcpt - well)/rw(1)
          end if
          
       end if
    end do

    first = .false.
  end function arc_head

  ! an approximation to the passive arc source using a sum of point sources
  function pointhead(p,calcX,calcY,alpha, kbg, arcq, arcR, arcX, arcy, &
       & arct1, arct2, calc_head, calc_velx) result(PH)
    use constants, only : DP, CZERO, TWOPI
    use complex_bessel, only : cbesk
    implicit none

    complex(DP), intent(in), dimension(:) :: p
    complex(DP), dimension(size(p)) :: PH
    complex(DP), dimension(size(p)) :: q, k0, k1
    integer, parameter :: PTS = 30 !! number of points to approximate line
    real(DP) :: r, deltheta, thetarange
    integer :: i, j, np, nz, ierr
    real(DP), dimension(PTS) :: x,y

    ! parameters of the line source
    real(DP), intent(in) :: alpha,calcX, calcY, kbg, arcr, arcx, arcy, &
         & arct1, arct2, arcQ
    logical, intent(in) :: calc_head, calc_velx

    np = size(p)
    q(1:np) = sqrt(p(:)/alpha)

    thetarange = arct2 - arct1
    deltheta = thetarange/real(PTS-1,DP)

    x(1:PTS) = arcX + arcR*cos(arct1 + deltheta*real((/(i,i=0,PTS-1)/),DP))
    y(1:PTS) = arcY + arcR*sin(arct1 + deltheta*real((/(i,i=0,PTS-1)/),DP))

    open(unit=55,file='arc_pt.bdry')

    PH = CZERO

    do i=1,pts
       write(55,*) x(i),y(i)
       !! radial distance from calc point to point on arc
       r = sqrt((calcX-x(i))**2 + (calcY-y(i))**2)
       
       !! compute required bessel functions for each value of p
       do j=1, np
          if (calc_head) then
             call cbesk(r*q(j),0.0_DP,1,1,k0(j),nz,ierr)
             if (nz /=0 .or. ierr /= 0) then
                print *, 'K0 bessel function error', nz, ierr
             end if
          else
             call cbesk(r*q(j),1.0_DP,1,1,k1(j),nz,ierr)
             if (nz /=0 .or. ierr /= 0) then
                print *, 'K1 bessel function error', nz, ierr
             end if
          end if
       end do
       
       
       if (calc_head) then ! head
          PH(1:np) = PH + k0(1:np)/kbg
       elseif (calc_velx) then ! x flux
          PH(1:np) = PH - q(1:np)*k1(1:np)*(calcX-x(i))/r
       else ! y flux
          PH(1:np) = PH - q(1:np)*k1(1:np)*(calcY-y(i))/r
       end if
       
    end do
    close(55)

    ! factor out common stuff to apply here
    PH = PH*arcq/(TWOPI*real(PTS,DP)*p(1:np)*kbg)
    
  end function pointhead
  
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  !! takes a vector of K bessel functions, (order 0 to nk)
  !! computes a vector of K' using recurrence
  function besk_deriv(bk,arg,Darg) result(Dbk)
    use constants, only : DP
    complex(DP), intent(in), dimension(0:) :: bk
    complex(DP), intent(in) :: arg, Darg
    complex(DP), dimension(0:size(bk)-1) :: Dbk

    integer :: nk

    nk = size(bk)-1

    Dbk(0) = -bk(1)*Darg
    Dbk(1:nk-1) = -0.5_DP*(bk(0:nk-2) + bk(2:nk))*Darg
    Dbk(nk) = -(bk(nk-1) + real(nk,DP)/arg*bk(nk))*Darg
    
  end function besk_deriv
 
   function besi_deriv(bi,arg,Darg) result(Dbi)
    use constants, only : DP
    complex(DP), intent(in), dimension(0:) :: bi
    complex(DP), intent(in) :: arg, Darg
    complex(DP), dimension(0:size(bi)-1) :: Dbi

    integer :: ni

    ni = size(bi)-1

    Dbi(0) = bi(1)*Darg
    Dbi(1:ni-1) = 0.5_DP*(bi(0:ni-2) + bi(2:ni))*Darg
    Dbi(ni) = (bi(ni-1) - real(ni,DP)/arg*bi(ni))*Darg
    
  end function besi_deriv

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v

    integer :: i
    real(DP) :: rnum, range

    if(lo >= hi) stop "LINSPACE: lower bound must be &
         &less than upper bound."

    rnum = real(num - 1,DP)
    range = hi - lo

    v = (/ ( lo + real(i,DP)*range/rnum, i=0,num-1) /)

  end function linspace

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  function linspace_no_ends(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v

    integer :: i
    real(DP) :: rnum, range

    if(lo >= hi) stop "LINSPACE_NO_ENDS: lower bound &
         &must be less than upper bound."

    rnum = real(num + 1,DP)
    range = hi - lo

    v = (/ ( lo + real(1+i,DP)*range/rnum, i=0,num-1) /)

  end function linspace_no_ends

!!$  !! given a low and high theta for the first interval, give
!!$  !! two vectors of points distributed in a chebyshev-style 
!!$  !! spacing over each interval (but covering all of the -pi to +pi 
!!$  !! interval (taking wrapping into account)
!!$  function chebspace(th1,th2,num) result(v)
!!$    use constants, only : DP, PI, TWOPI
!!$    real(DP), intent(in) ::  th1, th2
!!$    integer, intent(in) :: num
!!$    real(DP), dimension(2,num) :: v
!!$
!!$    real(DP) :: range1, range2
!!$    real(DP), dimension(num) :: cheb1,cheb2
!!$
!!$    cheb1(1:num) = 0.5_DP*(1.0_DP - cos(linspace(0.0_DP,PI,num)))
!!$    cheb2(1:num) = 0.5_DP*(1.0_DP - cos(linspace_no_ends(0.0_DP,PI,num)))
!!$
!!$    if(th1 < th2) then
!!$       !! +/-pi cut in second interval
!!$       range1 = th2-th1
!!$       v(1,:) = th1 + range1*cheb1(:)
!!$       
!!$       range2 = TWOPI - range1
!!$       !! th2 is now th1
!!$       v(2,:) = th2 + range2*cheb2(:)
!!$
!!$       !! wrap interval around
!!$       where (v(2,:) > PI)
!!$          v(2,:) = v(2,:) - TWOPI
!!$       end where
!!$    else
!!$       !! +/-pi cut in first interval
!!$       range2 = th1-th2
!!$       v(2,:) = th2 + range2*cheb1(:)
!!$       
!!$       range1 = TWOPI - range2
!!$       !! th2 is now th1
!!$       v(1,:) = th1 + range1*cheb2(:)
!!$
!!$       !! wrap interval around
!!$       where (v(1,:) > PI)
!!$          v(1,:) = v(1,:) - TWOPI
!!$       end where
!!$    end if
!!$  end function chebspace

  function outerprod(a,b) result(c)
    use constants, only : DP
    real(DP), dimension(:), intent(in) :: a,b
    real(DP), dimension(size(a),size(b)) :: c

    c = spread(a,dim=2,ncopies=size(b))*&
         & spread(b,dim=1,ncopies=size(a))

  end function outerprod
  


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
  

end module arc_test
