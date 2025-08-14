!
! Copyright (c) 2011-2025 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module constants
  implicit none

  public

  ! NB: there are a pile of constants in cbessel, not used elsewhere.
  
  ! real with range 300 orders of mag, 15 sig figs
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  ! length of filenames
  integer, parameter :: lenFN = 128

  ! useful? constants related to pi and ln
  real(DP), parameter :: PI = 4.0_DP*atan(1.0_DP) ! 3.1415926535897931...
  real(DP), parameter :: TWOPI = 8.0_DP*atan(1.0_DP) 
  real(DP), parameter :: PIOV2 = 2.0_DP*atan(1.0_DP) ! 1.5707963267948966...
  real(DP), parameter :: PISQ = (4.0_DP*atan(1.0_DP))**2

  ! these are both used in YNOT approximation to elliptical circumference
  real(DP), parameter :: LN2 = log(2.0_DP) ! 0.69314718055994529...
  real(DP), parameter :: LNPIOV2 = log(2.0_DP*atan(1.0_DP))  ! 0.45158270528945482...

  ! sqrt(-1)
  complex(DP), parameter :: EYE = cmplx(0.0_DP,1.0_DP,DP)

end module constants

