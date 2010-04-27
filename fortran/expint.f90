module expint_approx
implicit none
private
public expint

contains

!----------EXPONENTIAL INTEGRAL-----------------------------
! checked against expint() function in Matlab over range
! 10^-5 < x < 1000.0 for an error < 10-7
function expint(x) result(e1)
  use constants, only : DP, RZERO, RONE ! for double precision constant 
  implicit none

  real(DP), intent(in) :: x
  real(DP)             :: e1
  real(DP), dimension(0:5), parameter :: a = (/ &
       & -0.57721566_DP,  0.99999193_DP, -0.24991055_DP, &
       &  0.05519968_DP, -0.00976004_DP,  0.00107857_DP /)
  real(DP), dimension(4), parameter :: c = (/ &
       & 8.5733287401_DP, 18.0590169730_DP, &
       & 8.6347608925_DP,  0.2677737343_DP /)
  real(DP), dimension(4), parameter ::  b = (/ &
       &  9.5733223454_DP, 25.6329561486_DP, &
       & 21.0996530827_DP,  3.9584969228_DP /) 

  if (x .ge. RZERO) then
     if (x .le. RONE) then
        ! equation 5.1.53 in Abramowitz & Stegun
        e1 = a(0) + x*(a(1) + x*(a(2) + x*(a(3) + &
             & x*(a(4) + x*a(5))))) - log(x)
     else 
        ! equation 5.1.56
        e1 = ((((x + c(1))*x + c(2))*x + c(3))*x + c(4)) / &
           &(((((x + b(1))*x + b(2))*x + b(3))*x + b(4))*x*exp(x))

     end if
  else
     write (*,*) "EXPINT: error. x must be greater than or equal to zero"
     stop
  end if
end function expint
end module expint_approx
