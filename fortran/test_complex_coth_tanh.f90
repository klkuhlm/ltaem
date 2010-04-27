program test

  implicit none

  !! routines compare well with Mathematica for positive x and y quadrant

  integer, parameter :: DP = selected_real_kind(15,300), NUM=100

  complex(DP), dimension(NUM) :: z
  real(DP), dimension(NUM) :: x,y
  integer :: i

  x = [(real(i*1.0D-5,DP),i=0,NUM-1)] 
  y = [(real(i*5.0D-1,DP),i=0,NUM-1)] 

  z = cmplx(x,y,DP)
  
  open(unit=44,file='test_tanh_coth.out')

  do i=1,NUM
     write(44,'(6(ES15.7E3,1X))') z(i),tanh_double_complex(z(i)),coth_double_complex(z(i))
  end do
  
contains  

  !########################################
  elemental function tanh_double_complex(z) result(fz)
    !! definition in terms of real & imag components from
    !! Baker, "Less Complex Elementary Functions"
    complex(DP), intent(in) :: z
    complex(DP) :: fz
    real(DP) :: x, y

    x = real(z)
    y = aimag(z)

    fz = cmplx(tanh(2.0_DP*x)/(1.0_DP + cos(2.0_DP*y)/cosh(2.0_DP*x)), &
         & sin(2.0_DP*y)/(cosh(2.0_DP*x) + cos(2.0_DP*y)))

  end function tanh_double_complex
  
  elemental function coth_double_complex(z) result(fz)
    complex(DP), intent(in) :: z
    complex(DP) :: fz
    
    fz = 1.0_DP/tanh_double_complex(z)

  end function coth_double_complex
end program test
