module constants
  implicit none

  public
  
  ! real with range 300 orders of mag, 15 sig figs 
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  ! length of filenames 
  integer, parameter :: lenFN = 128

  ! useful? constants related to pi and ln
  real(DP), parameter :: PI = 4.0_DP*atan(1.0_DP) ! 3.1415926535897931...
  real(DP), parameter :: PIOV2 = 2.0_DP*atan(1.0_DP) ! 1.5707963267948966...

  ! these are both used in YNOT approximation to elliptical circumference
  real(DP), parameter :: LN2 = log(2.0_DP) ! 0.69314718055994529...
  real(DP), parameter :: LNPIOV2 = log(2.0_DP*atan(1.0_DP))  ! 0.45158270528945482...

  ! sqrt(-1)
  complex(DP), parameter :: EYE = cmplx(0.0_DP,1.0_DP,DP)

end module constants

