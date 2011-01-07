program test
  implicit none

  integer, parameter :: DP = 8, N = 25
  integer, dimension(N) :: vi
  real(DP), dimension(N,N) :: c
  real(DP), dimension(N) :: x,s
  complex(DP), dimension(N) :: z
  integer :: i
  
  forall (i=0:N-1) vi(i+1) = i
  forall (i=0:N-1) x(i+1) = 0.2*i

  ! outer product form	(simple, but no recurrence)
  c = cos(spread(x,2,N)*spread(vi,1,N))

  print *, 'matrix ops'
  write(*,'(3X,12(I21,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F21.18,1X))') i,c(i,1:12)
  end do
  
  ! de Moivre's recurrence (all orders computed from first one)
  ! (requires both sine & cosine, uses complex powers)
  ! http://www.trans4mind.com/personal_development/mathematics/trigonometry/deMoivre.htm
  c(:,1) = 1.0
  c(:,2) = cos(x)
  s(:) = sin(x)
  z = cmplx(c(:,2),s,DP)
  forall (i=2:N-1) c(:,i+1) = real(z**i)
  
  print *, "de Moivre's recurrance"
  write(*,'(3X,12(I21,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F21.18,1X))') i,c(i,1:12)
  end do
  
  ! Chebyshev 2-level recurrence  (no powers, but must be done in order)
  c = cos_recurrence(x,N)

  print *, "Chebyshev's recurrance"
  write(*,'(3X,12(I21,1X))') vi(1:12)
  do i=1,N
     write(*,'(I2,1X,12(F21.18,1X))') i,c(i,1:12)
  end do
  
  contains

  function cos_recurrence(x,n) result(c)
    real(DP), intent(in), dimension(:) :: x
    integer, intent(in) :: n
    real(DP), dimension(size(x),0:n-1) :: c
    integer :: i

    c(:,0) = 1.0
    c(:,1) = cos(x)
    do i = 2,N-2
       c(:,i) = 2.0*c(:,1)*c(:,i-1) - c(:,i-2)
    end do
  end function cos_recurrence

!!$  function sin_recurrence(x,n) result(s)
!!$    real(DP), intent(in), dimension(:) :: x
!!$    integer, intent(in) :: n
!!$    real(DP), dimension(size(x),1:n-1) :: s
!!$    real(DP), dimension(size(x)) :: c
!!$    integer :: i
!!$
!!$    s(:,1) = sin(x)
!!$    s(:,2) = sin(2.0*x)
!!$    c = cos(x)
!!$    do i = 3,N-2
!!$       s(:,i) = 2.0*c(:)*s(:,i-1) - s(:,i-2)
!!$    end do
!!$  end function sin_recurrence

end program test
