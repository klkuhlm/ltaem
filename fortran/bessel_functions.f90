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

  ! K Bessel function for 2D matrix argument / vector of N
  function besk_matz(z,num) result(besk)
    complex(DP), dimension(:,:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,1),size(z,2),0:num-1) :: besk
    integer :: ierr1, ierr2, i,j

    do i = 1, size(z,1)
       do j = 1, size(z,2)
          call cbesk(z(i,j), real(first,DP), 1, num, besk(i,j,:), ierr1, ierr2)
          ! either 0 or 3 are acceptable return codes
          if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
             write(*,'(A,4(1X,I0))') 'besk_matz error',ierr1,ierr2,i,j
             stop 
          end if
       end do
    end do
  end function besk_matz

  ! K Bessel function for vector argument / vector of N
  function besk_vectz(z,num) result(besk)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,1),0:num-1) :: besk
    besk = sum(besk_matz(spread(z,dim=2,ncopies=1),num),dim=2)
  end function besk_zvect

  ! K Bessel function for scalar argument / vector of N
  function besk_zscal(z,num) result(besk)
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: besk
    besk = sum(sum(besk_matz(spread([z],dim=2,ncopies=1),num),dim=2),dim=1)
  end function besk_zscal

  ! I Bessel function for 2D matrix argument / vector of N
  function besi_matz(z,num) result(besi)
    complex(DP), dimension(:,:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,1),size(z,2),0:num-1) :: besi
    integer :: ierr1, ierr2, i,j

    do i = 1, size(z,1)
       do j = 1, size(z,2)
          call cbesi(z(i,j), real(first,DP), 1, num, besi(i,j,:), ierr1, ierr2)
          ! either 0 or 3 are acceptable return codes
          if ((ierr2 >= 1 .and. ierr2 <= 2) .or. ierr2 >= 4) then
             write(*,'(A,4(1X,I0))') 'besi_matz error',ierr1,ierr2,i,j
             stop 
          end if
       end do
    end do
  end function besi_matz

  ! I Bessel function for vector argument / vector of N
  function besi_vectz(z,num) result(besi)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,1),0:num-1) :: besi
    besi = sum(besi_matz(spread(z,dim=2,ncopies=1),num),dim=2)
  end function besi_zvect

  ! I Bessel function for scalar argument / vector of N
  function besi_zscal(z,num) result(besi)
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: besi
    besi = sum(sum(besi_matz(spread([z],dim=2,ncopies=1),num),dim=2),dim=1)
  end function besi_zscal

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
  end subroutine bId

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine bKd(z,n,K,KD)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD

    K(0:n-1) = bK(z,n)
    KD(1:n-2) = -0.5_DP*(K(0:n-3) + K(2:n-1))    ! middle
    KD(0) = - K(1)                               ! low end
    KD(n-1) = -(K(n-2) + real(n-1,DP)/z*K(n-1))  ! high end
  end subroutine bKd

end module bessel_functions

