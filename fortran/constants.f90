module constants
  implicit none
  
  public
  
  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind (p=15,r=300)

  ! length of filenames 
  integer, parameter :: lenFN = 128

  ! useful? constants related to pi and ln
  real(kind=DP), parameter :: PI =      4.0_DP*atan(1.0_DP)
  real(kind=DP), parameter :: PIOV2 =   PI/2.0_DP

  ! Euler calculated to precision=34 using mathematica (could be quad precision)
  real(kind=DP), parameter :: EULER =   0.5772156649015328606065120900824025_DP

  ! these are both used in YNOT approximation to elliptical circumference
  real(kind=DP), parameter :: LN2 =     log(2.0_DP)
  real(kind=DP), parameter :: LNPIOV2 = log(PIOV2)

  ! sqrt(-1)
  complex(kind=DP), parameter :: EYE = cmplx(0.0,1.0,DP)

end module constants

