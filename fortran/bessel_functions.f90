module bessel_functions

  ! this module is a wrapper for the complex Bessel funcitons
  ! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995

  use Complex_Bessel, only : cbesk, cbesi
  use constants, only : DP

  implicit none

  private
  public :: bK, bI, bId, bIk

  interface bK
     module procedure besk_zsingle, besk_zvect, besk_zvect_nsingle, besk_zscal
  end interface bK

  interface bI
     module procedure besi_zsingle, besi_zvect, besi_zscal
  end interface bI

contains

  ! K Bessel function for scalar argument / vector of N
  function besk_zscal(z,first,num) result(besk)
    complex(DP), intent(in) :: z
    integer, intent(in) :: first, num
    complex(DP), dimension(first:first+num-1) :: besk
    integer :: ierr1, ierr2

    call cbesk(z, real(first,DP), 1, num, besk(:), ierr1, ierr2)
    ! either 0 or 3 are acceptable return codes
    if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
       write(*,'(A,2(1X,I0))') 'besk_zscal error',ierr1,ierr2
       stop 
    end if
  end function besk_zscal

  ! K Bessel function for vector argument / vector of N
  function besk_zvect(z,first,num) result(besk)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: first, num
    complex(DP), dimension(size(z,1),first:first+num-1) :: besk
    integer :: ierr1, ierr2, i

    do i = 1, size(z,1)
       call cbesk(z(i), real(first,DP), 1, num, besk(i,:), ierr1, ierr2)
       if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
          write(*,'(A,2(1X,I0))') 'besk_zvect error',ierr1,ierr2
          stop 
       end if
    end do
  end function besk_zvect

  ! K Bessel function for single argument, single value of N
  function besk_zsingle(z,n) result(besk)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP) :: besk
    complex(DP), dimension(1) :: temp
    integer :: ierr1, ierr2

    call cbesk(z, real(n,DP), 1, 1, temp, ierr1, ierr2)
    if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
       write(*,'(A,2(1X,I0))') 'besk_zsingle error',ierr1,ierr2
       stop          
    end if
    besk = temp(1) ! convert to scalar
  end function besk_zsingle

  ! K Bessel function for vector argument / one N
  function besk_zvect_nsingle(z,n) result(besk)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), dimension(size(z,1)) :: besk
    integer :: ierr1, ierr2, i

    do i = 1, size(z,1)
       call cbesk(z(i), real(n,DP), 1, 1, besk(i), ierr1, ierr2)
       if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
          write(*,'(A,2(1X,I0))') 'besk_zvect_nsingle error',ierr1,ierr2
          stop 
       end if
    end do
  end function besk_zvect_nsingle

  ! I Bessel function for scalar argument / vector of N
  function besi_zscal(z,first,num) result(besi)
    complex(DP), intent(in) :: z
    integer, intent(in) :: first, num
    complex(DP), dimension(first:first+num-1) :: besi
    integer :: ierr1, ierr2

    call cbesi(z, real(first,DP), 1, num, besi(:), ierr1, ierr2)
    if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
       write(*,'(A,2(1X,I0))') 'besi_zscal error',ierr1,ierr2
       stop 
    end if
  end function besi_zscal

  ! I Bessel function for vector argument / vector of N
  function besi_zvect(z,first,num) result(besi)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: first, num
    complex(DP), dimension(size(z,1), first:first+num-1) :: besi
    integer :: ierr1, ierr2, i

    do i = 1, size(z,1)
       call cbesi(z(i), real(first,DP), 1, num, besi(i,:), ierr1, ierr2)
       if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
          write(*,'(A,2(1X,I0))') 'besi_zvect error',ierr1,ierr2
          stop 
       end if
    end do
  end function besi_zvect

  ! K Bessel function for single argument, single value of N
  function besi_zsingle(z,n) result(besi)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP) :: besi
    complex(DP), dimension(1) :: temp
    integer :: ierr1, ierr2

    call cbesi(z, real(n,DP), 1, 1, temp, ierr1, ierr2)
    if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
       write(*,'(A,2(1X,I0))') 'besi_zsingle error',ierr1,ierr2
       stop 
    end if
    besi = temp(1) ! convert to scalar
  end function besi_zsingle

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine bId(z,n,I,ID)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID

    I(0:n-1) = bI(z,n)
    ID(1:n-2) = 0.5_DP*(I(0:n-3) + I(2:n-1))   ! middle
    ID(0) = I(1)                               ! low end
    ID(n-1) = I(n-2) - real(n-1,DP)/z*I(n-1)   ! high end
  end subroutine BID

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine bkd(z,n,K,KD)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD

    K(0:n-1) = bK(z,n)
    KD(1:n-2) = -0.5_DP*(K(0:n-3) + K(2:n-1))    ! middle
    KD(0) = - K(1)                               ! low end
    KD(n-1) = -(K(n-2) + real(n-1,DP)/z*K(n-1))  ! high end
  end subroutine Bkd

end module bessel_functions

