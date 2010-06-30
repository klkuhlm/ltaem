module bessel_functions

  ! this module is a wrapper for the complex Bessel funcitons
  ! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995

  use Complex_Bessel, only : cbesk, cbesi
  use constants, only : DP
  implicit none

  private
  public :: bK, bI, dbI, dbK

  interface bK
     module procedure besk_zscal, besk_vectz
  end interface bK
  interface bI
     module procedure besi_zscal, besi_vectz
  end interface bI

  interface dbK
     module procedure beskd_zscal, beskd_zvect
  end interface dbK
  interface dbI
     module procedure besid_zscal, besid_zvect
  end interface dbI

contains

  ! K Bessel function for vector argument / vector of N
  function besk_vectz(z,num) result(besk)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,dim=1),0:num-1) :: besk
    integer :: numzero, ierr, i
    integer, parameter :: kode = 1

    do i = 1, size(z,dim=1)
       call cbesk(z(i), 0.0_DP, kode, num, besk(i,0:num-1), numzero, ierr)
       ! either 0 or 3 are acceptable return codes
       if ((ierr >= 1 .and. ierr <= 2) .or. ierr >= 4) then
          write(*,'(A,3(1X,I0))') 'besk_vectz error',numzero,ierr,i
          stop 222
       end if
    end do
  end function besk_vectz

  ! K Bessel function for scalar argument / vector of N
  function besk_zscal(z,num) result(besk)
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: besk
    besk = sum(besk_vectz([z],num),dim=1)
  end function besk_zscal

  ! I Bessel function for vector argument / vector of N
  function besi_vectz(z,num) result(besi)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,dim=1),0:num-1) :: besi
    integer :: numzero, ierr, i
    integer, parameter :: kode = 1

    do i = 1, size(z,dim=1)
       call cbesi(z(i), 0.0_DP, kode, num, besi(i,0:num-1), numzero, ierr)
       ! either 0 or 3 are acceptable return codes
       if ((ierr >= 1 .and. ierr <= 2) .or. ierr >= 4) then
          write(*,'(A,3(1X,I0))') 'besi_vectz error',numzero,ierr,i
          stop 223
       end if
    end do
  end function besi_vectz

  ! I Bessel function for scalar argument / vector of N
  function besi_zscal(z,num) result(besi)
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: besi
    besi = sum(besi_vectz([z],num),dim=1)
  end function besi_zscal

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! these routines do not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine besId_zvect(z,n,I,ID)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,dim=1),0:n-1) :: I, ID
    complex(DP), dimension(size(z,dim=1),0:max(2,n)-1) :: Itmp
    
    Itmp(:,0:max(n,2)-1) = bI(z,max(2,n))
    ID(:,0) = I(:,1)   ! low end
    if (n >= 2) then
       I(:,0:n-1) = Itmp(:,0:max(n,2)-1)
       ID(:,n-1) = I(:,n-2) - (n-1)/z(:)*I(:,n-1) ! high end
       if (n >= 3) then
          ID(:,1:n-2) = 0.5_DP*(I(:,0:n-3) + I(:,2:n-1)) ! middle
       end if
    else
       I(:,0) = Itmp(:,0)
    end if
  end subroutine besId_zvect

  subroutine besId_zscal(z,n,I,ID)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID
    complex(DP), dimension(1,0:n-1) :: tI, tId
    call besId_zvect([z],n,tI(1,0:n-1),tId(1,0:n-1))
     I(0:n-1) =  tI(1,0:n-1)
    ID(0:n-1) = tId(1,0:n-1)
  end subroutine besId_zscal

  subroutine besKd_zvect(z,n,K,KD)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,dim=1),0:n-1) :: K, KD
    complex(DP), dimension(size(z,dim=1),0:max(n,2)-1) :: Ktmp
    Ktmp(:,0:max(n,2)-1) = bK(z,max(n,2))
    KD(:,0) = -Ktmp(:,1)  ! low end
    if (n >= 2) then
       K(:,0:n-1) = Ktmp(:,0:max(n,2)-1)
       KD(:,n-1) = -(K(:,n-2) + (n-1)/z*K(:,n-1)) ! high end
       if (n >= 3) then
          KD(:,1:n-2) = -0.5_DP*(K(:,0:n-3) + K(:,2:n-1)) ! middle
       end if
    else
       K(:,0) = Ktmp(:,0)
    end if
  end subroutine besKd_zvect

  subroutine besKd_zscal(z,n,K,KD)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD
    complex(DP), dimension(1,0:n-1) :: tK, tKd
    call besKd_zvect([z],n,tK(1,0:n-1),tKd(1,0:n-1))
     K(0:n-1) =  tK(1,0:n-1)
    KD(0:n-1) = tKd(1,0:n-1)
  end subroutine besKd_zscal

end module bessel_functions

