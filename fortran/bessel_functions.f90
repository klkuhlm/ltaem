module bessel_functions

  ! this module is a wrapper for the complex Bessel funcitons
  ! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995

  use Complex_Bessel, only : cbesk, cbesi
  use constants, only : DP
  implicit none

  private
  public :: bK, bI, dbI, dbK

  interface bK
     module procedure besk_zscal, besk_vectz, besk_matz
  end interface bK
  interface bI
     module procedure besi_zscal, besi_vectz, besi_matz
  end interface bI

  interface dbK
     module procedure beskd_zscal, beskd_zvect, beskd_zmat
  end interface dbK
  interface dbI
     module procedure besid_zscal, besid_zvect, besid_zmat
  end interface dbI

contains

  ! K Bessel function for 2D matrix argument / vector of N
  function besk_matz(z,num) result(besk)
    complex(DP), dimension(:,:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,1),size(z,2),0:num-1) :: besk
    integer :: nz, ierr, i,j
    integer, parameter :: kode = 1

    do i = 1, size(z,1)
       do j = 1, size(z,2)
          call cbesk(z(i,j), 0.0_DP, kode, num, besk(i,j,0:num-1), nz, ierr)
          ! either 0 or 3 are acceptable return codes
          if ((ierr >= 1 .and. ierr <= 2) .or. ierr >= 4) then
             write(*,'(A,4(1X,I0))') 'besk_matz error',nz,ierr,i,j
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
  end function besk_vectz

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
    integer :: nz, ierr, i,j
    integer, parameter :: kode = 1

    do i = 1, size(z,1)
       do j = 1, size(z,2)
          call cbesi(z(i,j), 0.0_DP, kode, num, besi(i,j,0:num-1), nz, ierr)
          ! either 0 or 3 are acceptable return codes
          if ((ierr >= 1 .and. ierr <= 2) .or. ierr >= 4) then
             write(*,'(A,4(1X,I0))') 'besi_matz error',nz,ierr,i,j
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
  end function besi_vectz

  ! I Bessel function for scalar argument / vector of N
  function besi_zscal(z,num) result(besi)
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: besi
    besi = sum(sum(besi_matz(spread([z],dim=2,ncopies=1),num),dim=2),dim=1)
  end function besi_zscal

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! these routines do not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine besId_zmat(z,n,I,ID)
    complex(DP), dimension(:,:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,1),size(z,2),0:n-1) :: I, ID
    complex(DP), dimension(size(z,1),size(z,2),0:max(2,n)-1) :: Itmp
    
    Itmp(:,:,0:max(n,2)-1) = bI(z,max(2,n))
    ID(:,:,0) = I(:,:,1)   ! low end
    if (n >= 2) then
       I = Itmp
       ID(:,:,n-1) = I(:,:,n-2) - (n-1)/z*I(:,:,n-1) ! high end
       if (n >= 3) then
          ID(:,:,1:n-2) = 0.5_DP*(I(:,:,0:n-3) + I(:,:,2:n-1)) ! middle
       end if
    else
       I(:,:,0) = Itmp(:,:,0)
    end if
  end subroutine besId_zmat

  subroutine besId_zvect(z,n,I,ID)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,1),0:n-1) :: I, ID
    complex(DP), dimension(size(z,1),1,0:n-1) :: tI, tId
    call besId_zmat(spread(z,2,1),n,tI,tId)
    I = tI(:,1,:)
    ID = tId(:,1,:)
  end subroutine besId_zvect

  subroutine besId_zscal(z,n,I,ID)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID
    complex(DP), dimension(1,1,0:n-1) :: tI, tId
    call besId_zmat(spread([z],2,1),n,tI,tId)
    I = tI(1,1,:)
    ID = tId(1,1,:)
  end subroutine besId_zscal

  subroutine besKd_zmat(z,n,K,KD)
    complex(DP), dimension(:,:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,1),size(z,2),0:n-1) :: K, KD
    complex(DP), dimension(size(z,1),size(z,2),0:max(n,2)-1) :: Ktmp
    Ktmp(:,:,0:max(n,2)-1) = bK(z,max(n,2))
    KD(:,:,0) = -Ktmp(:,:,1)  ! low end
    if (n >= 2) then
       K = Ktmp
       KD(:,:,n-1) = -(K(:,:,n-2) + (n-1)/z*K(:,:,n-1)) ! high end
       if (n >= 3) then
          KD(:,:,1:n-2) = -0.5_DP*(K(:,:,0:n-3) + K(:,:,2:n-1)) ! middle
       end if
    else
       K(:,:,0) = Ktmp(:,:,0)
    end if
  end subroutine besKd_zmat

  subroutine besKd_zvect(z,n,K,KD)
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,1),0:n-1) :: K, KD
    complex(DP), dimension(size(z,1),1,0:n-1) :: tK, tKd
    call besKd_zmat(spread(z,2,1),n,tK,tKd)
    K = tK(:,1,:)
    KD = tKd(:,1,:)
  end subroutine besKd_zvect

  subroutine besKd_zscal(z,n,K,KD)
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD
    complex(DP), dimension(1,1,0:n-1) :: tK, tKd
    call besKd_zmat(spread([z],2,1),n,tK,tKd)
    K = tK(1,1,:)
    KD = tKd(1,1,:)
  end subroutine besKd_zscal

end module bessel_functions

