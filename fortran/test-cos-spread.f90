program test
  implicit none

  integer, parameter :: DP = 8, N = 200, M=N*10
  integer, dimension(N) :: vi
  real(DP), dimension(M,N) :: c
  real(DP), dimension(M) :: x !!,s
!!$  complex(DP), dimension(M) :: z
  integer :: i
  real(DP) :: t0,t1
  
  intrinsic :: cpu_time

!!$  forall (i=0:N-1) vi(i+1) = i
  forall (i=0:M-1) x(i+1) = 0.05*i

  call cpu_time(t0)

!!$
!!$  ! outer product form	(simple, but no recurrence)
!!$  c = cos(spread(x,2,N)*spread(vi,1,M))
!!$
  
!!$  ! de Moivre's recurrence (all orders computed from first one)
!!$  ! (requires both sine & cosine, uses complex powers)
!!$  ! http://www.trans4mind.com/personal_development/mathematics/trigonometry/deMoivre.htm
!!$  c(:,1) = 1.0
!!$  c(:,2) = cos(x)
!!$  s(:) = sin(x)
!!$  z = cmplx(c(:,2),s,DP)
!!$  forall (i=2:N-1) c(:,i+1) = real(z**i)
!!$  

  ! Chebyshev 2-level recurrence  (no powers, but must be done in order)
  c = cos_recurrence(x,N)

  call cpu_time(t1)

  print *, 'time:',t1-t0

  write(*,'(3X,200(I23,1X))') vi(1:N)
  do i=1,M
     write(*,'(I0,1X,200(F23.20,1X))') i,c(i,1:N)
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
