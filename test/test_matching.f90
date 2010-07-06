! $Id: test_matching.f90,v 1.1 2006/08/08 20:29:36 kris Exp kris $
module new_matching
  implicit none

  private
  public :: new_match

contains

  function new_match(p) result(F)
    use calc_shared_data, only : calcX, calcY, LINk, LINalpha, LINSs
    use complex_bessel, only : cbesk, cbesi
    use constants, only : DP, TWOPI, PI, CZERO, RONE, SQRT2OV2

    complex(DP), dimension(:), intent(in) :: p
    complex(DP), dimension(size(p)) :: F

    complex(DP), dimension(size(p)) :: qout, qin
    integer :: nP, i, k, nz, ierr

    ! dimensions of stuff
    integer, parameter :: N = 8, M = 3*(2*N+1)
    real(DP), dimension(0:N), parameter :: vi = real([(i, i=0,N)],DP)

    complex(DP), dimension(0:N+1) :: besk, besi, besk0, besi0
    complex(DP), save, allocatable :: a(:,:),b(:,:),c(:,:),d(:,:)
    complex(DP), save, dimension(0:N) :: vn
    
    complex(DP), dimension(2*M,0:2*N) :: Am
    complex(DP), dimension(2*M) :: bm        ! RHS vector
    real(DP), dimension(0:2*N) :: Diag     ! norm vector

    real(DP) :: rc, rw, tc
    complex(DP), dimension(1) :: well

    ! vectors of a_n and b_n for each value of p
    complex(DP), allocatable, save :: x(:,:)

    complex(DP), dimension(M,0:1) :: beskwell ! order zero and one k bessel functions
    real(DP), dimension(M) :: phiM, r

    ! is this the first time through this subroutine?
    logical, save :: first = .true. 

    nP = size(p)

    qout(1:nP) = sqrt(p(1:nP)/LINalpha)
    qin(1:np)  = sqrt(p(1:np)*LINSs)         ! k inside circle = 1

    ! matching points along circumference of ellipse
    if (first) then
       phiM(1:M) = linspace(0.0_DP,2.0_DP*PI,M)
       allocate(a(0:2*N,np),b(0:2*N,np),c(0:2*N,np),d(0:2*N,np),x(0:2*N,np))
    end if

    do i = 1,np
       if (first) then

          Am(1:M,0:N) =     cos(outer_prod(phiM(1:M),vi(0:N)))
          Am(1:M,N+1:2*N) = sin(outer_prod(phiM(1:M),vi(1:N)))
          Am(M+1:2*M,:) = Am(1:M,:)

          r(1:M) = sqrt((2.0_DP - cos(phiM(:)))**2 + (2.0_DP - sin(phiM(:)))**2)

          do k = 1,M
             call cbesk(r(k)*qout(i),0.0_DP,1,2,beskwell(k,0:1),nz,ierr)
             if (nz /=0 .or. ierr /= 0) print *, 'bessel function error', nz, ierr
          end do

          ! head effects of background well
          bm(1:M) = -beskwell(1:M,0)/(TWOPI*LINk*p(i))

          ! flux effect of background well
          bm(M+1:2*M) = +beskwell(1:M,1)*qout(i)*SQRT2OV2/(TWOPI*LINk*p(i)) &
               & *(cos(phiM(1:M)) + sin(phiM(1:M)))

          ! discrete version of norm calculation
          forall(k = 0:2*N)
             Diag(k) = real(dot_product(Am(1:2*M,k),Am(1:2*M,k)))
          end forall

          ! this is similar to a discrete fourier transform
          ! least-squares solution for coefficients
          x(0:2*N,i) = matmul(spread(RONE/Diag(0:2*N),dim=2,ncopies=2*M)* &
               & transpose(conjg(Am(1:2*M,0:2*N))),bm(1:2*M))

          call cbesk(qout(i),0.0_DP,1,N+2,besk(0:N+1),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel K function error', nz, ierr
          call cbesi(qin(i),0.0_DP,1,N+2,besi(0:N+1),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel I function error', nz, ierr
          
          ! k == 0
          vn(0) = -qout(i)*besk(1)/(qin(i)*besi(1))*besi(0)/besk(0)
          a(0,i) = x(0,i)/(1.0_DP/LINk - vn(0))
          c(0,i) = x(0,i)/(1.0_DP/(vn(0)*LINk) - 1.0_DP)
          b(0,i) = CZERO; d(0,i) = CZERO

          do k = 1,N
             vn(k) = -qout(i)/qin(i)*(besk(k-1) + besk(k+1))/(besi(k-1) + besi(k+1))&
                  & *besi(k)/besk(k)
             a(k,i) = x(k,i)/(1.0_DP/LINk - vn(k))
             b(k,i) = x(N+k,i)/(1.0_DP/LINk - vn(k))
             c(k,i) = x(k,i)/(1.0_DP/(vn(k)*LINk) - 1.0_DP)
             d(k,i) = x(N+k,i)/(1.0_DP/(vn(k)*LINk) - 1.0_DP)
          end do
       end if

       ! distance to orgin-centered circle
       rc = sqrt(calcX**2 + calcY**2)
       tc = atan2(calcY,calcX)

       ! outside circle
       if(rc > 1.0_DP) then

          ! distance to pumping well
          rw = sqrt((2.0_DP - calcX)**2 + (2.0_DP - calcY)**2)

          call cbesk(rw*qout(i),0.0_DP,1,1,well(1),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel K function error', nz, ierr
          F(i) = well(1)/(TWOPI*p(i))

          call cbesk(qout(i)*rc,0.0_DP,1,N+1,besk(0:N),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel K function error', nz, ierr
          call cbesk(qout(i),0.0_DP,1,N+1,besk0(0:N),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel K function error', nz, ierr

          F(i) = F(i) + a(0,i)*besk(0)/besk0(0)

          do k=1,N
             F(i) = F(i) + &
                  & besk(k)/besk0(k)*(a(k,i)*cos(tc) + b(k,i)*sin(tc))
          end do


       ! inside (or on circumference of) circle
       else
          call cbesi(qin(i)*rc,0.0_DP,1,N+1,besi(0:N),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel I function error', nz, ierr
          call cbesi(qin(i),0.0_DP,1,N+1,besi0(0:N),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel I function error', nz, ierr
          
          F(i) = c(0,i)*besi(0)/besi0(0)

          do k=1,N
             F(i) = F(i) + &
                  & besi(k)/besi0(k)*(c(k,i)*cos(tc) + d(k,i)*sin(tc))
          end do

       end if
    end do

    first = .false.
  end function new_match

 
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

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  function logspace(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v

    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)

  end function logspace

  function outer_prod(a,b) result(c)
    use constants, only : DP
    real(DP), dimension(:), intent(in) :: a,b
    real(DP), dimension(size(a),size(b)) :: c

    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

  end function outer_prod
  
end module new_matching
