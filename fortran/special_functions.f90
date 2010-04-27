! $Id: special_functions.f90,v 1.9 2008/03/30 20:24:37 kris Exp kris $

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! there are FOUR freaking modules in this file
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

module error_handler
  implicit none
  
  private
  public :: fileError, subError
  
  contains

  !##################################################
  ! error handler for file opening errors -> kill program
  subroutine fileError(fname,error,callsub,flag)

    ! fname: file which threw error trying to open
    ! error: integer error returned from open() call
    ! callsub: string identifying which subroutine/ what location error happened
    ! flag:  =0 -> write problem;  <>0 -> read problem
    character(len=128), intent(in) :: fname, callsub
    integer, intent(in) :: error, flag
    
    if (flag == 0) then
       print *, trim(callsub),': error',error,'opening file ',trim(fname),' for writing'
    else
       print *, trim(callsub),': error',error,'opening file ',trim(fname),' for reading'
    end if
    stop 'quitting due to file error'

  end subroutine fileError

  !##################################################
  ! error handler for external subroutine errors (bessel functions, matrix decomposition)
  subroutine subError(sub,error,callsub)

    ! sub: name of subroutine which passed back error code
    ! error: integer error returned by subroutine
    ! callsub: string identifying which subroutine/ what location error happened
    character(len=128), intent(in) :: sub, callsub
    integer, intent(in) :: error
    
    print *, trim(callsub),': error',error,' returned from ',trim(sub)
    stop 'quitting due to subroutine error'
    
  end subroutine subError
end module error_handler

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! an interface to cbesk and cbesi subroutines with
! wrappers to handle passing vectors/matricies as argument

module bessel_functions
  use Complex_Bessel, only : cbesk, cbesi
  use Error_Handler, only : subError
  use constants, only : DP
  
  implicit none

  private
  public :: besselk, besseli


  interface besselk
     module procedure besk_zsingle, besk_zvect, besk_zvect_nsingle, besk_zscal
  end interface

  interface besseli
     module procedure besi_zsingle, besi_zvect, besi_zscal
  end interface

  contains

    ! K Bessel function for scalar argument / vector of N
    function besk_zscal(z,first,num) result(besk)
      complex(DP), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(first:first+num-1) :: besk
      character(128) :: sn = 'besk_zscal', kn = 'cbesk'
      integer :: ierr1, ierr2

      call cbesk(z, real(first,DP), 1, num, besk(:), ierr1, ierr2)
      ! either 0 or 3 are acceptable return codes
      if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
         call subError(sn,ierr2,kn)
      end if
      

    end function besk_zscal

    ! K Bessel function for vector argument / vector of N
    function besk_zvect(z,first,num) result(besk)
      complex(DP), dimension(:), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(size(z,1),first:first+num-1) :: besk
      character(128) :: sn = 'besk_zvect', kn = 'cbesk'
      integer :: ierr1, ierr2, i
      
      do i = 1, size(z,1)
         call cbesk(z(i), real(first,DP), 1, num, besk(i,:), ierr1, ierr2)
         if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
            call subError(sn,ierr2,kn) 
         end if
         
      end do

    end function besk_zvect

    ! K Bessel function for single argument, single value of N
    function besk_zsingle(z,n) result(besk)
      complex(DP), intent(in) :: z
      integer, intent(in) :: n
      complex(DP) :: besk
      complex(DP), dimension(1) :: temp
      character(128) :: sn = 'besk_zsingle', kn = 'cbesk'
      integer :: ierr1, ierr2
      
      call cbesk(z, real(n,DP), 1, 1, temp, ierr1, ierr2)
      if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
         call subError(sn,ierr2,kn)
      end if
      
      besk = temp(1) ! convert to scalar
      
    end function besk_zsingle

    ! K Bessel function for vector argument / one N
    function besk_zvect_nsingle(z,n) result(besk)
      complex(DP), dimension(:), intent(in) :: z
      integer, intent(in) :: n
      complex(DP), dimension(size(z,1)) :: besk
      character(128) :: sn = 'besk_zvect', kn = 'cbesk'
      integer :: ierr1, ierr2, i

      do i = 1, size(z,1)
         call cbesk(z(i), real(n,DP), 1, 1, besk(i), ierr1, ierr2)
         if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
            call subError(sn,ierr2,kn)
         end if
      end do

    end function besk_zvect_nsingle

    ! I Bessel function for scalar argument / vector of N
    function besi_zscal(z,first,num) result(besi)
      complex(DP), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(first:first+num-1) :: besi
      character(128) :: sn = 'besi_zscal', kn = 'cbesi'
      integer :: ierr1, ierr2

      call cbesi(z, real(first,DP), 1, num, besi(:), ierr1, ierr2)
      if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
         call subError(sn,ierr2,kn)
      end if
      
    end function besi_zscal

    ! I Bessel function for vector argument / vector of N
    function besi_zvect(z,first,num) result(besi)
      complex(DP), dimension(:), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(size(z,1), first:first+num-1) :: besi
      character(128) :: sn = 'besi_zvect', kn = 'cbessi'
      integer :: ierr1, ierr2, i

      do i = 1, size(z,1)
         call cbesi(z(i), real(first,DP), 1, num, besi(i,:), ierr1, ierr2)
         if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
            call subError(sn,ierr2,kn)
         end if
      end do

    end function besi_zvect

    ! K Bessel function for single argument, single value of N
    function besi_zsingle(z,n) result(besi)
      complex(DP), intent(in) :: z
      integer, intent(in) :: n
      complex(DP) :: besi
      complex(DP), dimension(1) :: temp
      character(128) :: sn = 'besi_zsingle', kn = 'cbesi'
      integer :: ierr1, ierr2
      
      call cbesi(z, real(n,DP), 1, 1, temp, ierr1, ierr2)
      if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
         call subError(sn,ierr2,kn)
      end if
            
      besi = temp(1) ! convert to scalar
      
    end function besi_zsingle

end module bessel_functions


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

module matrix_inverse
  use constants, only : DP
  implicit none
  private
  public :: inverse !, decomp, backsub, ludcmp
contains


  !##################################################
  ! using double-precision complex routines from LAPACK 3.3
  function inverse(ai) result(inv)

    ! interfaces to LAPACK generated by ifort -gen-interfaces
    ! LAPACK LU decomposition 
    INTERFACE 
       SUBROUTINE ZGETRF(M,N,A,LDA,IPIV,INFO)
         INTEGER, intent(in) :: LDA, M, N
         COMPLEX(KIND=8), intent(inout) :: A(LDA,*)
         INTEGER, intent(inout) :: IPIV(*)
         INTEGER, intent(inout) :: INFO
       END SUBROUTINE ZGETRF
    END INTERFACE

    ! LAPACK inverse calculation from results of LU
    INTERFACE 
       SUBROUTINE ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
         INTEGER, intent(in) :: LDA, N
         COMPLEX(KIND=8), intent(inout) :: A(LDA,*)
         INTEGER, intent(inout) :: IPIV(*)
         COMPLEX(KIND=8), intent(inout) :: WORK(*)
         INTEGER, intent(in) :: LWORK
         INTEGER, intent(inout) :: INFO
       END SUBROUTINE ZGETRI
    END INTERFACE

    integer :: n, ierr
    integer, parameter :: LWORK = 6400
    complex(DP), dimension(:,:), intent(in) :: ai
    complex(DP), dimension(size(ai,1),size(ai,1)) :: inv
    integer, dimension(size(ai,1)) :: indx
    complex(DP), dimension(LWORK) :: work

    indx = 0
    n = size(ai,1)
    inv = ai    

    call zgetrf(n,n,inv,n,indx,ierr)
    if (ierr /= 0) write (*,*) 'error returned from ZGETRF'

    call zgetri(n,inv,n,indx,work,LWORK,ierr)
    if (ierr /= 0) write (*,*) 'error returned from ZGETRI'
    
  end function inverse

end module matrix_inverse

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

module leaky_q
  use constants, only : DP
  implicit none

  interface compute_CIleaky_q
     module procedure  compute_CIleaky_qv, compute_CIleaky_qs
  end interface

!!$  real(DP), save, allocatable :: root(:,:)

  private
  public :: compute_CIleaky_q

contains

  !! vector version useful during calc phase
  function compute_CIleaky_qv(p) result(q)
    use constants, only : DP, PI
    use element_specs, only : CInum, a2v,b2v,kv,av,leakv,k2v,b2v,sv,syv,kzv,unconfv,bgb

    integer, parameter :: NTERMS = 200, MAXITER = 200
    
    complex(DP), intent(in), dimension(:) :: p
    complex(DP), dimension(size(p),0:CInum) :: q

    integer :: i, ni, np
    complex(DP), dimension(size(p),0:CInum) :: kap2
    complex(DP), dimension(size(p)) :: exp2z

    !! boulton stuff
!!$    real(DP), dimension(NTERMS) :: guess, gamma
!!$    real(DP), dimension(0:CInum) :: sigma
!!$    complex(DP), dimension(size(p)) :: kernel
!!$    real(DP) :: x, delta
!!$    integer :: k, kk
!!$    logical, save :: first = .true.
    real(DP) :: boulton

    np = size(p)
    ni = CInum
!!$    sigma(0:ni) = sqrt(sv/Syv)

    do i=0,ni
       !! leaky-ness
       !! ##############################
       if(leakv(i) == 0) then
          !! no leaky layer, standard definition
          q(1:np,i) = p(1:np)/av(i)
       else
          kap2(1:np,i) = sqrt(p(:)/a2v(i))
          exp2z(1:np) = exp(-2.0_DP*kap2(:,i)*b2v(i))

          if(leakv(i) == 1) then
             !! case I, no-drawdown condition at top of aquitard
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))*&
                  & (1.0_DP + exp2z(:))/(1.0_DP - exp2z(:))
          elseif(leakv(i) == 2) then
             !! case II, no-flow condition at top of aquitard
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))*&
                  & (1.0_DP - exp2z(:))/(1.0_DP + exp2z(:))
          elseif(leakv(i) == 3) then
             !! aquitard thickness -> infinity
             q(:,i) = p(:)/av(i) + kap2(:,i)*k2v(i)/(bgb*kv(i))
          else
             stop 'ERROR: incorrect value for CIAquitardLeak parameter -> (1,2,3)'
          end if
       end if
       
       !! unconfined-ness 
       !! ##############################
       if(unconfv(i) == 0) then
          !! do nothing, q already computed above
       else
          !! Boulton unconfined source (Herrera infinite sum Kernel)
          !! guess is halfway between asymptotes of cot()
!!$
!!$          if (first) then
!!$             if(.not. allocated(root)) then
!!$                allocate(root(NTERMS,0:CInum))
!!$             end if
!!$             
!!$             !! roots are not a function of p (just sigma)- only compute once
!!$             guess(2:NTERMS) = PI*(real((/(k, k=1,NTERMS-1)/)) + 0.5_DP)/sigma(i)
!!$             guess(1) = 1.7D0
!!$
!!$             !! first root is hard to find with NR, 
!!$             !! use TS approximation for tangent and re-arrange
!!$             x = guess(1)
!!$             NR1: do kk = 1,MAXITER
!!$                delta = (x + (sigma(i) - 1.0_DP/sigma(i))*(x*sigma(i) + (x*sigma(i))**3/3.0_DP + &
!!$                     & 2.0_DP*(sigma(i)*x)**5/15.0_DP) + 17.0_DP*(x*sigma(i))**7/315.0_DP)/ &
!!$                     & (1.0_DP - (1.0_DP/sigma(i) - sigma(i))*(sigma(i) + x**2*sigma(i)**3 + &
!!$                     & 2.0_DP*x**4*sigma(i)**5/3.0_DP + 17.0_DP*x**6*sigma(i)**7/45.0_DP))
!!$                x = x - delta
!!$                if (abs(delta) <= 1.0D-10) then
!!$                   root(1,i) = x
!!$                   exit NR1
!!$                end if
!!$                if(kk == MAXITER) print *, '1 failed to converge'
!!$             end do NR1
!!$
!!$             do k = 2, NTERMS
!!$                x = guess(k)
!!$                NR: do kk = 1,MAXITER
!!$                   delta = (1.0_DP/tan(x*sigma(i)) + (sigma(i) - 1.0_DP/sigma(i))/x)/&
!!$                        & (sigma(i)/(sin(sigma(i)*x)**2) + (sigma(i) + 1.0_DP/sigma(i))/x**2)
!!$                   x = x + delta
!!$                   if (abs(delta) <= spacing(x)*10.0) then
!!$                      root(k,i) = x
!!$                      exit NR
!!$                   end if
!!$                   if(kk == MAXITER) print *, k,'failed to converge'
!!$                end do NR
!!$             end do
!!!!!$             if (i == ni) first = .false.
!!!!!$          end if
!!$
!!$          gamma(1:NTERMS) = Kzv(i)*root(1:NTERMS,i)**2/(bgb*Syv(i))
!!$
!!$          kernel(1:np) = 2.0_DP*sum(spread(gamma(1:NTERMS),2,np)/&
!!$               & ((spread(root(1:NTERMS,i)**2,2,np) - 1.0_DP + sigma(i)**2)* &
!!$               & (spread(p(1:np),1,NTERMS) + spread(gamma(1:NTERMS),2,np))),dim=1)
!!$          
!!$
!!$          q(1:np,i) = q(1:np,i) + p(1:np)*Syv(i)/Kv(i)*kernel

          ! scrap Herrera's infinite sum for Boulton's original
          ! rough-n-ready alpha, with a semi-physical expression for it
          boulton = 3.0_DP*Kzv(i)/(bgb*Syv(i))
          q(1:np,i) = q(:,i) + Syv(i)*p(1:np)*boulton/(Kv(i)*(boulton + p(1:np)))
       end if
    end do
    
    !! sources are only additive under the square root
    q = sqrt(q);

  end function compute_CIleaky_qv
  
  !! scalar version useful in matching
  function compute_CIleaky_qs(p) result(q)
    use constants, only : DP
    use element_specs, only : CInum

    complex(DP), intent(in) :: p
    complex(DP), dimension(0:CInum) :: q

    !! sum away singleton first dimension
    q(0:CInum) = sum(compute_CIleaky_qv( (/p/)),dim=1)

  end function compute_CIleaky_qs
  

end module leaky_q
  
