!
! Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module geomConv
  implicit none

  private
  public :: xy2cR,xy2cA, xy2eR,xy2eA, c2xyR,c2xyA, e2xyR,e2xyA

contains

  ! ****************************************
  elemental function xy2cR(z) result(zeta)
    ! polar coordinate from local complex Cartesian with respect to circle
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: zeta
    zeta = cmplx(abs(z), atan2(aimag(z),real(z)),DP)
  end function xy2cR
  elemental function xy2cA(z,c) result(zeta)
    ! polar coordinate from global complex Cartesian
    use constants, only : DP
    use type_definitions, only : circle
    complex(DP), intent(in) :: z
    type(circle), intent(in) :: c
    complex(DP) :: zeta
    zeta = xy2cR(z - c%z)
  end function xy2cA

  elemental function xy2eR(z,e) result(zeta)
    ! elliptical coordinate from unrotated local complex Cartesian with respect to ellipse
    use constants, only : DP, EYE
    use type_definitions, only : ellipse
    use utility, only : cacosh
    complex(DP), intent(in) :: z
    type(ellipse), intent(in) :: e
    complex(DP) :: zeta
    zeta = cacosh(z*exp(-EYE*e%theta)/e%f)
  end function xy2eR
  elemental function xy2eA(z,e) result(zeta)
    ! elliptical coordinate from global complex Cartesian
    use constants, only : DP
    use type_definitions, only : ellipse
    complex(DP), intent(in) :: z
    type(ellipse), intent(in) :: e
    complex(DP) :: zeta
    zeta = xy2eR(z - e%z,e)
  end function xy2eA

  elemental function c2xyR(zeta) result(z)
    ! local Cartesian from elemental polar
    use constants, only : DP, EYE
    complex(DP), intent(in) :: zeta
    complex(DP) :: z
    z = real(zeta)*exp(EYE*aimag(zeta))
  end function c2xyR
  elemental function c2xyA(zeta,c) result(z)
    ! global Cartesian from elemental polar
    use constants, only : DP
    use type_definitions, only : circle
    complex(DP), intent(in) :: zeta
    type(circle), intent(in) :: c
    complex(DP) :: z
    z = c%z + c2xyR(zeta)
  end function c2xyA

  elemental function e2xyR(zeta,e) result(z)
    ! local unrotated Cartesian from elemental elliptical
    use constants, only : DP, EYE
    use type_definitions, only : ellipse
    use utility, only : ccosh
    complex(DP), intent(in) :: zeta
    type(ellipse), intent(in) :: e
    complex(DP) :: z
    z = e%f*ccosh(zeta)*exp(EYE*e%theta)
  end function e2xyR
  elemental function e2xyA(zeta,e) result(z)
    ! global Cartesian from elemental elliptical
    use constants, only : DP
    use type_definitions, only : ellipse
    complex(DP), intent(in) :: zeta
    type(ellipse), intent(in) :: e
    complex(DP) :: z
    z = e%z + e2xyR(zeta,e)
  end function e2xyA

end module geomConv
