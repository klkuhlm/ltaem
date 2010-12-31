module geomConv
  
  private 
  public :: xy2cR,xy2cA, xy2eR,xy2eA, c2xyR,c2xyA, e2xyR,e2xyA 
  
contains

  ! ****************************************
  elemental function xy2cR(z) result(zeta) 
    ! polar coord from local complex cartesian wrt circle
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: zeta
    zeta = cmplx(abs(z), atan2(aimag(z),real(z)),DP)
  end function xy2cR
  elemental function xy2cA(z,c) result(zeta)
    ! polar coord from global complex cartesian
    use constants, only : DP
    use type_definitions, only : circle
    complex(DP), intent(in) :: z
    type(circle), intent(in) :: c
    complex(DP) :: zeta
    zeta = xy2cR(z - c%z)
  end function xy2cA
  
  elemental function xy2eR(z,e) result(zeta)
    ! elliptical coord from unrotated local complex cartesian wrt ellipse
    use constants, only : DP, EYE
    use type_definitions, only : ellipse
    use utility, only : cacosh
    complex(DP), intent(in) :: z
    type(ellipse), intent(in) :: e
    complex(DP) :: zeta
    zeta = cacosh(z*exp(-EYE*e%theta)/e%f)
  end function xy2eR
  elemental function xy2eA(z,e) result(zeta)
    ! elliptical coord from global complex cartesian
    use constants, only : DP
    use type_definitions, only : ellipse
    complex(DP), intent(in) :: z
    type(ellipse), intent(in) :: e
    complex(DP) :: zeta
    zeta = xy2eR(z - e%z,e)
  end function xy2eA

  elemental function c2xyR(zeta) result(z)
    ! local cartesian from elemental polar 
    use constants, only : DP, EYE
    complex(DP), intent(in) :: zeta
    complex(DP) :: z
    z = real(zeta)*exp(EYE*aimag(zeta))
  end function c2xyR
  elemental function c2xyA(zeta,c) result(z)
    ! global cartesian from elemental polar
    use constants, only : DP
    use type_definitions, only : circle
    complex(DP), intent(in) :: zeta
    type(circle), intent(in) :: c
    complex(DP) :: z
    z = c%z + c2xyR(zeta)
  end function c2xyA
  
  elemental function e2xyR(zeta,e) result(z)
    ! local unrotated cartesian from elemental elliptical
    use constants, only : DP, EYE
    use type_definitions, only : ellipse
    use utility, only : ccosh
    complex(DP), intent(in) :: zeta
    type(ellipse), intent(in) :: e
    complex(DP) :: z
    z = e%f*ccosh(zeta)*exp(EYE*e%theta)
  end function e2xyR
  elemental function e2xyA(zeta,e) result(z)
    ! global cartesian from elemental elliptical
    use constants, only : DP
    use type_definitions, only : ellipse
    complex(DP), intent(in) :: zeta
    type(ellipse), intent(in) :: e
    complex(DP) :: z
    z = e%z + e2xyR(zeta,e)
  end function e2xyA

end module geomConv
