program test

complex(8) z
real(8) x,y

x=0.0
y=100.0

z=cmplx(x,y)

print*, z,acosh(z)




contains
  elemental function acosh(z) result(f)
    use constants, only : DP, RONE

    complex(DP), intent(in) :: z
    complex(DP) :: f

    f = log(z + sqrt((z-RONE)*(z+RONE)))

!!$    if(real(z) > 0.0) then
!!$       f = log(z + sqrt(z**2 - RONE))
!!$    else
!!$       f = log(z - sqrt(z**2 - RONE))
!!$    end if

  end function acosh
end program test
