program test_ellipse_circumference
  use constants, only : DP
  implicit none
  
  real(DP), dimension(10) :: eta, f
  integer :: i,j

  eta = [(real(i,DP)*0.125_DP, i=0,9)]
  f = [(real(i,DP)*1.0_DP, i=0,9)]
  
  open(unit=22,file='ellipse-circumference.dat',status='replace',action='write')
  write(22,*) '# checking YNOT formulat for circumference of an ellipse'
  write(22,*) '# elliptical radisu (eta), semi-focal, semi-major, semi-minor, circumference'

  do i=1,10
     do j=1,10
        write(22,888) eta(i),f(j),circ(eta(i),f(j))
     end do
     write(22,'(//)') 
  end do
  
888 format(3(1X,ES22.14E3))

contains

  ! return the approximate circumference of an ellipse
  ! given the radius in elliptical coordinates, and 
  ! the semi-focal distance (YNOT formula)
  function circ(eta,f) result(P)
    use constants, only : DP, LN2, LNPIOV2, PI
    real(DP), intent(in) :: eta, f
    real(DP) :: P, y

    ! semi-major/minor length a:=f*cosh/sinh(eta)

    if(f < tiny(1.0D0)) then ! a circle
       P = PI*eta**2
    else
       y = log(2.0_DP)/log(PI/2.0_DP) ! YNOT constant
       P = 4.0_DP*((f*cosh(eta))**y + (f*sinh(eta))**y)**(1.0_DP/y)
    end if
    

  end function circ

end program test_ellipse_circumference
