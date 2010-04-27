program test_mf
use mathieu_functions
implicit none

type(mathieu) :: m
integer :: i,n
real(8) :: x,y

write(*,'(A)') 'q (re,im), M'
read(*,*) x,y,n

m = mathieu_init(q=cmplx(x,y,8),MM=n)

open(22,file='Aev.out')
write(22,*) '#',m%q
do i=1,m%M
   write(22,*) real(m%A(i,:,0))
end do
write(22,'()') 
do i=1,m%M
   write(22,*) aimag(m%A(i,:,0))
end do
close(22)

open(22,file='Aod.out')
write(22,*) '#',m%q
do i=1,m%M
   write(22,*) real(m%A(i,:,1))
end do
write(22,'()')
do i=1,m%M
   write(22,*) aimag(m%A(i,:,1))
end do
close(22)

open(22,file='Bev.out')
write(22,*) '#',m%q
do i=1,m%M
   write(22,*) real(m%B(i,:,0))
end do
write(22,'()')
do i=1,m%M
   write(22,*) aimag(m%B(i,:,0))
end do
close(22)

open(22,file='Bod.out')
write(22,*) '#',m%q
do i=1,m%M
   write(22,*) real(m%B(i,:,1))
end do
write(22,'()') 
do i=1,m%M
   write(22,*) aimag(m%B(i,:,1))
end do
close(22)

open(22,file='mcn.out')
write(22,*) '#',m%q
do i=1,4
   write(22,*) real(m%mcn(m%M*(i-1)+1:m%M*i))
end do
write(22,'()')
do i=1,4
   write(22,*) aimag(m%mcn(m%M*(i-1)+1:(m%M*i)))
end do


end program test_mf


