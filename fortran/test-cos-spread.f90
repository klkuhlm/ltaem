program test
  implicit none

  integer, parameter :: DP = 8, N = 25
!  real(DP), parameter :: PI = atan(1.0)*4.0
  integer, dimension(N) :: vi
  real(DP), dimension(N,N) :: c
  real(DP), dimension(N) :: x
  integer :: i
  
  forall (i=0:N-1) vi(i+1) = i
  forall (i=0:N-1) x(i+1) = 0.2*i

  c = cos(spread(x,2,N)*spread(vi,1,N))

  print *, 'matrix ops'
  write(*,'(3X,12(I7,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F7.4,1X))') i,c(i,1:12)
  end do
  
  c = -999.9
  
  c(:,1) = 1.0
  c(:,2) = cos(x)
  forall(i=2:N)
     c(1:ceiling(N/real(i)),i+1) = c(1:N:i,2)
  end forall
  
  print *, 'manual recurrance'
  write(*,'(3X,12(I7,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F7.4,1X))') i,c(i,1:12)
  end do
  
  c = c + 10.0

  where(c < -2.0)
     c = cos(spread(x,2,N)*spread(vi,1,N))
  end where

  print *, 'fill-in'
  write(*,'(3X,12(I7,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F7.4,1X))') i,c(i,1:12)
  end do

  
end program test
