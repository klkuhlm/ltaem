module calc_shared_data
  use constants, only : DP
  implicit none
  
  public

  ! the real/imaginary part of exp(-s)
  complex(DP), save :: intp
  logical, save :: intreal

  logical, save :: flux, match

  ! parameters of the line source
  real(DP), save :: LINq, LINalpha, LINSs, LINk, LINsf, LINeta0, LINgamma, LINEk

  logical, save :: calc_head, calc_velx
  
  real(DP), save :: calcX, calcY
  complex(DP), save, allocatable :: calcCoeff(:,:,:)
!!$  complex(DP), save, allocatable :: calcVn(:,:)
  

end module calc_shared_data
