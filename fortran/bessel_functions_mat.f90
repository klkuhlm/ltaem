module bessel_functions_deriv_mat

  ! this module is a wrapper for the complex Bessel funcitons
  ! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995

  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! this module returns SCALED bessel functions, which
  !! must be un-scaled by the calling program, according to
  !! how each routine works
  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  ! written by K Kuhlman, July 2006
  ! $Id: bessel_functions_mat.f90,v 1.1 2007/07/30 05:56:56 kris Exp kris $

  implicit none

  interface BesselI
     module procedure BesselI_val, BesselI_val_and_deriv
  end interface

  interface BesselK
     module procedure BesselK_val, BesselK_val_and_deriv
  end interface

  private
  public :: BesselI, BesselK

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains 

  subroutine BesselI_val(z,n,I)
    use constants, only : DP, SMALL
    use complex_bessel, only : cbesi

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I
    integer :: nz, ierr

    call cbesi(z, 0.0_DP, 2, n, I(0:n-1), nz, ierr)
    
    if (ierr /= 0) then
       select case(ierr)
       case(3)
          write(*,*) "CBESI: loss of precision, z=",z
       case(1)
          write(*,*) "CBESI: input error, z=",z," n=",n
          stop "CBESI: input error"
       case(2)
          write(*,*) "CBESI: overflow, z or order too large for unscaled output, z=",z," n=",n
          stop "CBESI: overflow, z or order too large for unscaled output"
       case(4)
          write(*,*) "CBESI: overflow, z or order too large, z=",z," n=",n
          stop "CBESI: overflow, z or order too large"
       case(5)
          stop "CBESI: algorithm termination not met"
       end select 
    end if

  end subroutine BesselI_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselI_val_and_deriv(z,n,I,ID)
    use constants, only : DP

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID

    call BesselI_val(z,n,I(0:n-1))

    ! middle
    ID(1:n-2) = 0.5_DP*(I(0:n-3) + I(2:n-1))

    ! low end
    ID(0) = I(1)

    ! high end
    ID(n-1) = I(n-2) - real(n-1,DP)/z*I(n-1)

  end subroutine BesselI_val_and_deriv

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val(z,n,K)
    use constants, only : DP, LARGE
    use complex_bessel, only : cbesk

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K
    integer :: nz, ierr

    call cbesk(z, 0.0_DP, 2, n, K(0:n-1), nz, ierr)
    
    if (ierr /= 0) then
       select case(ierr)
       case(3)
          write(*,*) "CBESK: loss of precision, z=",z
       case(1)
          write(*,*) "CBESK: input error, z=",z," n=",n
          stop "CBESK: input error"
       case(2)
          write(*,*) "CBESK: overflow, z too small or order too large for unscaled output, z=",z," n=",n
          stop "CBESK: overflow, z too small or order too large for unscaled output"
       case(4)
          write(*,*) "CBESK: overflow, z too small or order too large, z=",z," n=",n
          stop "CBESK: overflow, z too small or order too large"
       case(5)
          stop "CBESK: algorithm termination not met"
       end select
    end if

  end subroutine BesselK_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val_and_deriv(z,n,K,KD)
    use constants, only : DP

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD

    call BesselK_val(z,n,K(0:n-1))

    ! middle
    KD(1:n-2) = -0.5_DP*(K(0:n-3) + K(2:n-1))

    ! low end
    KD(0) = - K(1)

    ! high end
    KD(n-1) = -(K(n-2) + real(n-1,DP)/z*K(n-1))

  end subroutine BesselK_val_and_deriv

end module bessel_functions_deriv_mat

