program test_orthogonality

  use constants, only : DP, PI
  use mathieu_functions, only : mmatce, mmatse
  use shared_mathieu, only : A,B,q
  use mcn_matrix_method, only : mcn_eigenvalues
  use gauss_points, only : gauss_quad_setup

  integer, parameter :: ORD = 40, MAXI = 10
  real(DP), dimension(ORD) :: x,w,y
  complex(DP), dimension(ORD) :: ce,cec
  integer :: i, ms

  call gauss_quad_setup(ORD,x,w)

  q = cmplx(1.0,0.5,DP)
  ms = 20

  allocate(A(1:ms,0:ms,2),B(1:ms,0:ms,2))

  call mcn_eigenvalues(q,A,B,1)

  print *, 'matheiu parameter:', -q
  print *, 'GQ abscissa:',x
  y = (0.5_DP + x/2.0_DP)*2.0_DP*PI - PI
  print *, 'natural coordinates of abscissa:', y
  
  do i=0,MAXI
     ce = mmatce(i,y)
     cec = conjg(mmatce(i,y))
     print *, 'fcn values',i,ce
     print *, 'conjg fcn values',i,cec
     print *, 'z:',i,ce*cec
     print *, 'norm',i,sum(ce*cec*w)
     
     !! this is the norm/pi, so it should be 1 if correct

  end do
  
end program test_orthogonality

