program test
use constants, only : DP
use bessel_functions, only : bK,dbk
implicit none

integer, parameter :: nk = 2, nz = 6
complex(DP), dimension(nz,nk) :: ck,cdk
complex(DP), dimension(nz) :: z
integer :: j

z = [(cmplx(real(2*j,DP),real(j,DP),DP),j=1,nz)]
ck(1:nz,1:nk) = bK(z,nk)
print *, 'z',z
do j=1,nz
   print *, 'Re ck',j,real(ck(j,:))
end do
print *, ' '
do j=1,nz
   print *, 'Im ck',j,aimag(ck(j,:))
end do

print *, ' '
z = [(cmplx(real(2*j,DP),real(j,DP),DP),j=1,nz)]
!print *, 'z',z
call dbK(z,nk,ck(1:nz,1:nk),cdk(1:nz,1:nk))
!print *, 'z',z
do j=1,nz
   print *, 'Re ck',j,real(ck(j,:))
end do
print *, ' '
do j=1,nz
   print *, 'Im ck',j,aimag(ck(j,:))
end do
print *, '  '
print *, '  '
do j=1,nz
   print *, 'Re dck',j,real(cdk(j,:))
end do
print *, '  '
do j=1,nz
   print *, 'Im dck',j,aimag(cdk(j,:))
end do

end program test
