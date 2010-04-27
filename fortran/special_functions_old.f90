! -*-f90-*- $Id: special_functions.f90,v 1.8.1.3 2006/07/11 15:37:13 kris Exp kris $

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! there are THREE modules in this file
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
      complex(DP), dimension(first:num) :: besk
      character(128) :: sn = 'besk_zscal', kn = 'cbesk'
      integer :: ierr1, ierr2

      call cbesk(z, real(first,DP), 1, num, besk(first:num), ierr1, ierr2)
      if (ierr2 /= 0) call subError(sn,ierr2,kn)

    end function besk_zscal

    ! K Bessel function for vector argument / vector of N
    function besk_zvect(z,first,num) result(besk)
      complex(DP), dimension(:), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(size(z,1),1:num) :: besk
      character(128) :: sn = 'besk_zvect', kn = 'cbesk'
      integer :: ierr1, ierr2, i
      
      do i = 1, size(z,1)
         call cbesk(z(i), real(first,DP), 1, num, besk(i,1:num), ierr1, ierr2)
         if (ierr2 /= 0) call subError(sn,ierr2,kn) 
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
      if (ierr2 /= 0) call subError(sn,ierr2,kn)
      
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
         if (ierr2 /= 0)  call subError(sn,ierr2,kn)
      end do

    end function besk_zvect_nsingle

    ! I Bessel function for scalar argument / vector of N
    function besi_zscal(z,first,num) result(besi)
      complex(DP), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(1:num) :: besi
      character(128) :: sn = 'besi_zscal', kn = 'cbesi'
      integer :: ierr1, ierr2

      call cbesi(z, real(first,DP), 1, num, besi(1:num), ierr1, ierr2)
      if (ierr2 /= 0) call subError(sn,ierr2,kn)

    end function besi_zscal

    ! I Bessel function for vector argument / vector of N
    function besi_zvect(z,first,num) result(besi)
      complex(DP), dimension(:), intent(in) :: z
      integer, intent(in) :: first, num
      complex(DP), dimension(size(z,1), 1:num) :: besi
      character(128) :: sn = 'besi_zvect', kn = 'cbessi'
      integer :: ierr1, ierr2, i

      do i = 1, size(z,1)
         call cbesi(z(i), real(first,DP), 1, num, besi(i,1:num), ierr1, ierr2)
         if (ierr2 /= 0) call subError(sn,ierr2,kn)
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
      if (ierr2 /= 0) call subError(sn,ierr2,kn)
      
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

