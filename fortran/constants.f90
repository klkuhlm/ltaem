! $Id: constants.f90,v 1.1 2007/03/15 00:17:59 kris Exp kris $
module constants
  implicit none
  
  public   ! just to be explicit
  
  ! integers with range of 9 and 18 orders of magnitude
  integer, parameter :: I4B = selected_int_kind(r=9)   ! std 32-bit integer
  integer, parameter :: I8B = selected_int_kind(r=18)  ! std 64-bit integer

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind (p=15,r=300)

  ! real with 4800 orders of mag, and slightly higher precision (QR=10 on g95)
  ! on ifort (QR=16) this is quad precision, (p=33)
  integer, parameter :: QR = selected_real_kind (p=18,r=4800)

  ! useful? constants related to pi and ln
  ! calculated to precision=33 using ifort or mathematica
  real(kind=DP), parameter :: PI =        3.141592653589793238462643383279503_DP
  real(kind=DP), parameter :: TWOPI =     6.28318530717958647692528676655901_DP
  real(kind=DP), parameter :: PIOV2 =     1.57079632679489661923132169163975_DP
  real(kind=DP), parameter :: PIOV3 =     1.04719755119659774615421446109317_DP
  real(kind=DP), parameter :: PIOV4 =     0.785398163397448309615660845819876_DP
  real(kind=DP), parameter :: ONEOVPI =   0.318309886183790671537767526745029_DP
  real(kind=DP), parameter :: LN10  =     2.30258509299404568401799145468436_DP
  real(kind=DP), parameter :: LN2 =       0.693147180559945309417232121458177_DP
  real(kind=DP), parameter :: LNPIOV2 =   0.451582705289454864726195229894882_DP 
  real(kind=DP), parameter :: EULER =     0.5772156649015328606065120900824025_DP
  real(kind=DP), parameter :: PISQ =      9.86960440108935861883449099987615_DP
  real(kind=DP), parameter :: SQRTPI =    1.77245385090551602729816748334115_DP
  real(kind=DP), parameter :: SQRTPIOV2 = 1.25331413731550025120788264240552_DP
  real(kind=DP), parameter :: INVLOG10 =  0.434294481903251827651128918916605_DP

  ! a very small number (if < SMALL, it is effectively zero, compared to one)
  real(kind=DP), parameter :: SMALL = epsilon(1.0_DP) ! approx 2.22E-16

  ! the largest number which can be represented accurately
  real(kind=DP), parameter :: LARGE = huge(1.0_DP)  ! approx 1.79E+308

  ! commonly used constants
  complex(kind=DP), parameter :: CZERO = (0.0_DP, 0.0_DP)
  complex(kind=DP), parameter :: CONE  = (1.0_DP, 0.0_DP)
  complex(kind=DP), parameter :: CTWO =  (2.0_DP, 0.0_DP)
  complex(kind=DP), parameter :: EYE = (0.0_DP, 1.0_DP)
  real(kind=DP),  parameter :: RZERO = 0.0_DP
  real(kind=DP),  parameter :: RONE  = 1.0_DP
  real(kind=DP),  parameter :: RTWO =  2.0_DP
  real(kind=DP),  parameter :: HALF =  0.5_DP
  

  ! x and y cartesian unit vectors
  real(kind=DP), dimension(2), parameter :: UNITI = (/ RONE, RZERO /)
  real(kind=DP), dimension(2), parameter :: UNITJ = (/ RZERO, RONE /)

end module constants

