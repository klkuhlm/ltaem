Module Complex_Bessel

  ! provide interface to arb modified bessel function routines
  ! using interface like amos routines

  use constants, only : DP
  implicit none

  private
  public :: cbesk, cbesi

  ! interfaces to c wrappers
  
  interface
     function arb_K(nu,z,scaled) bind(c,name="arb_K") result(K)
       use, intrinsic :: iso_c_binding, only : C_DOUBLE_COMPLEX, C_DOUBLE, C_INT
       complex(C_DOUBLE_COMPLEX), intent(in), value :: z
       real(C_DOUBLE), intent(in), value :: nu
       integer(C_INT), intent(in), value :: scaled
       complex(C_DOUBLE_COMPLEX) :: K
     end function arb_K
  end interface
  interface
     function arb_I(nu,z,scaled) bind(c,name="arb_I") result(I)
       use, intrinsic :: iso_c_binding, only : C_DOUBLE_COMPLEX, C_DOUBLE, C_INT
       complex(C_DOUBLE_COMPLEX), intent(in), value :: z
       real(C_DOUBLE), intent(in), value :: nu
       integer(C_INT), intent(in), value :: scaled
       complex(C_DOUBLE_COMPLEX) :: I
     end function arb_I
  end interface
  contains

    ! use same call signature as amos routines    
    subroutine cbesk(z, fnu, kode, n, cy, nz, ierr)
      COMPLEX (dp), INTENT(IN)   :: z      ! argument
      REAL (dp), INTENT(IN)      :: fnu    ! lowest order requested
      INTEGER, INTENT(IN)        :: kode   ! unscaled =1, scaled /= 1
      INTEGER, INTENT(IN)        :: n      ! number orders requested
      COMPLEX (dp), INTENT(OUT)  :: cy(n)  ! result
      INTEGER, INTENT(OUT)       :: nz     ! number good values
      INTEGER, INTENT(OUT)       :: ierr   ! error code

      !!real(DP), dimension(n) :: rcy, icy
      
      integer :: i
      real(DP) :: order

      do i = 0,n-1
        order = fnu + real(i,DP)
        cy(i+1) = arb_K(order, z, kode)
      end do

      ierr = 0
      !!rcy = real(cy)
      !!icy = aimag(cy)
      !!
      !!if (any(abs(rcy) > huge(1.0_DP)) .or. any(abs(icy) > huge(1.0_DP))) then
      !!  ierr = 1
      !!end if
      !!
      !!if (any(rcy /= rcy) .or. any(icy /= icy)) then
      !!  ierr = 2
      !!end if
      nz = 0
      
    end subroutine cbesk
  
    subroutine cbesi(z, fnu, kode, n, cy, nz, ierr)
      COMPLEX (dp), INTENT(IN)   :: z
      REAL (dp), INTENT(IN)      :: fnu
      INTEGER, INTENT(IN)        :: kode
      INTEGER, INTENT(IN)        :: n
      COMPLEX (dp), INTENT(OUT)  :: cy(n)
      INTEGER, INTENT(OUT)       :: nz
      INTEGER, INTENT(OUT)       :: ierr
      
      !!real(DP), dimension(n) :: rcy, icy
      
      integer :: i
      real(DP) :: order

      do i=0,n-1
        order = fnu + real(i,DP)
        cy(i+1) = arb_I(order, z, kode)
      end do

      ierr = 0
      !!rcy = real(cy)
      !!icy = aimag(cy)
      !!
      !!if (any(abs(rcy) > huge(1.0_DP)) .or. any(abs(icy) > huge(1.0_DP))) then
      !!  ierr = 1
      !!end if
      !!
      !!if (any(rcy /= rcy) .or. any(icy /= icy)) then
      !!  ierr = 2
      !!end if
      nz = 0
      
    end subroutine cbesi  
  end Module Complex_Bessel
