module constants
  implicit none
  
  public
  
  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind (p=15,r=300)

  ! length of filenames 
  integer, parameter :: lenFN = 128

  ! useful? constants related to pi and ln
  real(kind=DP), parameter :: PI =      4.0_DP*atan(1.0_DP) ! 3.14159...
  real(kind=DP), parameter :: PIOV2 =   2.0_DP*atan(1.0_DP) ! 1.57079...

  ! these are both used in YNOT approximation to elliptical circumference
  real(kind=DP), parameter :: LN2 =     log(2.0_DP) ! 0.693147...
  real(kind=DP), parameter :: LNPIOV2 = log(PIOV2)  ! 0.451582...

  ! sqrt(-1)
  complex(kind=DP), parameter :: EYE = cmplx(0.0,1.0,DP)

end module constants

