program test_bessel
  use constants, only : DP
  use bessel_functions, only : bK,dbk
  implicit none

  integer, parameter :: M = 6, N = 3
  complex(DP), allocatable :: ck(:,:),cdk(:,:)
  complex(DP), dimension(0:N-1) :: ck0
  complex(DP), dimension(M) :: z
  integer :: j

  ! vector of values parallel to imaginary axis
  z = [(cmplx(12.5,j*2.5,DP),j=1,M)]
  print *, 'z: ',z
  allocate(ck(M,0:N-1),cdk(M,0:N-1))
  ck(1:M,0:N-1) = bK(z,N)
!!$  print *, ' ' 
!!$
!!$  do j=1,M
!!$     print *, 'Re ck',j,real(ck(j,0:N-1))
!!$  end do
!!$  print *, ' '
!!$  do j=1,M
!!$     print *, 'Im ck',j,aimag(ck(j,0:N-1))
!!$  end do
!!$  print *, ' '
!!$  print *, 'Re ck0',real(ck0(0:N-1))
!!$  print *, ' '
!!$  print *, 'Im ck0',aimag(ck0(0:N-1))
!!$  print *, ' '
!!$  print *, ' '

  call dbk(z,N,ck(1:M,0:N-1),cdk(1:M,0:N-1))
  ck0(0:N-1) = bK(z(1),N)

!!$  do j=1,M
!!$     print *, 'Re ck',j,real(ck(j,0:N-1))
!!$  end do
!!$  print *, ' '
!!$  do j=1,M
!!$     print *, 'Im ck',j,aimag(ck(j,0:N-1))
!!$  end do
!!$  print *, ' '
!!$  do j=1,M
!!$     print *, 'Re Dck',j,real(cdk(j,0:N-1))
!!$  end do
!!$  print *, ' '
!!$  do j=1,M
!!$     print *, 'Im Dck',j,aimag(cdk(j,0:N-1))
!!$  end do
  deallocate(ck,cdk)

  print *, '##############################'
  print *, 'stage 2'
  print *, '##############################'
!  z = [(cmplx(12.5,j*2.5,DP),j=1,M)]
  print *, 'z: ',z
  allocate(ck(M,0:N-1),cdk(M,0:N-1))
  ck(1:M,0:N-1) = bK(z,N)
  print *, ' ' 

  do j=1,M
     print *, 'Re ck',j,real(ck(j,0:N-1))
  end do
  print *, ' '
  do j=1,M
     print *, 'Im ck',j,aimag(ck(j,0:N-1))
  end do

  call dbk(z,N,ck(1:M,0:N-1),cdk(1:M,0:N-1))
  ck0(0:N-1) = bK(z(1),N)
  call dbk(z,N,ck(1:M,0:N-1),cdk(1:M,0:N-1))
  ck0(0:N-1) = bK(z(1),N)
  call dbk(z,N,ck(1:M,0:N-1),cdk(1:M,0:N-1))
  ck0(0:N-1) = bK(z(1),N)
  print *, ' '
  print *, 'Re ck0',real(ck0(0:N-1))
  print *, ' '
  print *, 'Im ck0',aimag(ck0(0:N-1))
  print *, ' '
  print *, ' '

  do j=1,M
     print *, 'Re ck',j,real(ck(j,0:N-1))
  end do
  print *, ' '
  do j=1,M
     print *, 'Im ck',j,aimag(ck(j,0:N-1))
  end do
  print *, ' '
  do j=1,M
     print *, 'Re Dck',j,real(cdk(j,0:N-1))
  end do
  print *, ' '
  do j=1,M
     print *, 'Im Dck',j,aimag(cdk(j,0:N-1))
  end do


end program test_bessel
