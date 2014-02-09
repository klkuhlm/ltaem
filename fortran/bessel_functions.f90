!
! Copyright (c) 2011-2014 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

module bessel_functions

  ! this module is a wrapper for the complex Bessel functions
  ! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995

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
  function besk_vectz(z,num) result(K)
    use constants, only : DP
    use complex_bessel, only : cbesk
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,dim=1),0:num-1) :: K
    complex(DP), dimension(0:num-1) :: tmp
    integer :: numzero, ierr, j

    do j = 1,size(z,dim=1)
       call cbesk(z=z(j), fnu=0.0_DP, kode=1, n=num, cy=tmp(0:num-1), nz=numzero, ierr=ierr)
       ! either 0 or 3 are acceptable return codes
       if (.not.(ierr == 0 .or. ierr == 3)) then
          write(stderr,*) 'CBESK_VECTZ error (numzero=',numzero,', ierr=',ierr,', j=',j,&
               & ', num=',num,') z(j)=',z(j),'z',z
          call abort() ! to dump for checking backtrace in gdb
          stop 222
       end if
       K(j,0:num-1) = tmp(0:num-1)
    end do
  end function besk_vectz

  ! K Bessel function for scalar argument / vector of N
  function besk_zscal(z,num) result(K)
    use constants, only : DP
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: K
    K = sum( besk_vectz([z],num), dim=1)
  end function besk_zscal

  ! I Bessel function for vector argument / vector of N
  function besi_vectz(z,num) result(I)
    use constants, only : DP
    use complex_bessel, only : cbesi
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(size(z,dim=1),0:num-1) :: I
    complex(DP), dimension(0:num-1) :: tmp
    integer :: numzero, ierr, j

    do j = 1,size(z,dim=1)
       call cbesi(z=z(j), fnu=0.0_DP, kode=1, n=num, cy=tmp(0:num-1), nz=numzero, ierr=ierr)
       ! either 0 or 3 are acceptable return codes
       if (.not.(ierr == 0 .or. ierr == 3)) then
          write(stderr,'(A,3(1X,I0))') 'CBESI_VECTZ error (numzero=',numzero,', ierr=',ierr,&
               &', j=',j, ', num=',num,') z(j)=',z(j),'z',z
          call abort()
          stop 223
       end if
       I(j,0:num-1) = tmp(0:num-1)
    end do
  end function besi_vectz

  ! I Bessel function for scalar argument / vector of N
  function besi_zscal(z,num) result(I)
    use constants, only : DP
    complex(DP), intent(in) :: z
    integer, intent(in) :: num
    complex(DP), dimension(0:num-1) :: I
    I = sum( besi_vectz([z],num), dim=1)
  end function besi_zscal

  ! use recurrence relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! these routines do not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine besId_zvect(z,n,I,ID)
    use constants, only : DP
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,dim=1),0:n-1) :: I, ID
    complex(DP), dimension(size(z,dim=1),0:max(2,n)) :: Itmp
    integer :: nz, mn
    nz = size(z,dim=1)
    mn = max(n,2)

    Itmp(1:nz,0:mn) = besi_vectz(z,mn+1)
    ID(1:nz,0) = Itmp(1:nz,1)   ! low end
    if (n >= 2) then
       I(1:nz,0:n-1) = Itmp(1:nz,0:n-1)
       ! since I(0) is finite, wrote this to not use 1/z form
       ! but it does require computing I of one higher order 
       ID(1:nz,n-1) = 0.5_DP*(I(1:nz,n-2) + Itmp(1:nz,n)) ! high end
       if (n >= 3) then
          ID(1:nz,1:n-2) = 0.5_DP*(I(1:nz,0:n-3) + I(1:nz,2:n-1)) ! middle
       end if
    else
       I(1:nz,0) = Itmp(1:nz,0)
    end if
  end subroutine besId_zvect

  subroutine besId_zscal(z,n,I,ID)
    use constants, only : DP
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID
    complex(DP), dimension(1,0:n-1) :: tI, tId
    call besId_zvect([z],n,tI(1,0:n-1),tId(1,0:n-1))
     I(0:n-1) =  tI(1,0:n-1)
    ID(0:n-1) = tId(1,0:n-1)
  end subroutine besId_zscal

  subroutine besKd_zvect(z,n,K,KD)
    use constants, only : DP
    complex(DP), dimension(:), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(size(z,dim=1),0:n-1) :: K, KD
    complex(DP), dimension(size(z,dim=1),0:max(n,2)) :: Ktmp
    integer :: nz, mn
    nz = size(z,dim=1)
    mn = max(n,2)

    Ktmp(1:nz,0:mn) = besk_vectz(z,mn+1)
    KD(1:nz,0) = -Ktmp(1:nz,1)  ! low end (always used)
    if (n >= 2) then
       K(1:nz,0:n-1) = Ktmp(1:nz,0:n-1)
       ! even though K(0) is typically infinite, wrote to not use 1/z
       ! but it does require computing K of one higher order
       KD(1:nz,n-1) = -0.5_DP*(K(1:nz,n-2) + Ktmp(1:nz,n)) ! high end
       if (n >= 3) then
          KD(1:nz,1:n-2) = -0.5_DP*(K(1:nz,0:n-3) + K(1:nz,2:n-1)) ! middle
       end if
    else
       ! only one order requested (n=1)
       K(1:nz,0) = Ktmp(1:nz,0)
    end if
  end subroutine besKd_zvect

  subroutine besKd_zscal(z,n,K,KD)
    use constants, only : DP
    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD
    complex(DP), dimension(1,0:n-1) :: tK, tKd
    call besKd_zvect([z],n,tK(1,0:n-1),tKd(1,0:n-1))
     K(0:n-1) =  tK(1,0:n-1)
    KD(0:n-1) = tKd(1,0:n-1)
  end subroutine besKd_zscal

end module bessel_functions

