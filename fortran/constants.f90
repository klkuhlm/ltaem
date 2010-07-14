module constants
  implicit none
  
  public
  
  ! integers with range of 9 and 18 orders of magnitude
  integer, parameter :: I4B = selected_int_kind(r=9)   ! std 32-bit integer
  integer, parameter :: I8B = selected_int_kind(r=18)  ! std 64-bit integer

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind (p=15,r=300)

  ! real with 4800 orders of mag, and slightly higher precision (QR=10 on g95)
  ! on ifort (QR=16) this is quad precision, (p=33)
  integer, parameter :: QR = selected_real_kind (p=18,r=4800)

  ! length of filenames 
  integer, parameter :: lenFN = 128

  ! useful? constants related to pi and ln
  ! calculated to precision=34 using mathematica (could be quad precision)
  real(kind=DP), parameter :: PI =      4.0_DP*atan(1.0_DP)
  real(kind=DP), parameter :: EULER =   0.5772156649015328606065120900824025_DP

  ! these are both used in YNOT approximation to elliptical circumference
  real(kind=DP), parameter :: LN2 =     log(2.0_DP)
  real(kind=DP), parameter :: LNPIOV2 = log(PI/2.0_DP)

  ! sqrt(-1)
  complex(kind=DP), parameter :: EYE = cmplx(0.0,1.0,DP)

end module constants

