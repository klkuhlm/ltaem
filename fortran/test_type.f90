program test_type

type test
   integer :: i,j,k
   real, allocatable :: x(:),y(:)
end type test

integer :: a,b,c
type(test), allocatable :: Z(:)

allocate(Z(5))

do a=1,5
   allocate(Z(a)%x(10),Z(a)%y(5))
   Z(a)%x = real(a)
   Z(a)%y = real(a)*10.0
end do

Z(1:5)%i = 100
Z(1:5)%j = 0
Z(1:5)%k = 22

do b=1,5
   print *, b,':',Z(b)%x,Z(b)%y,Z(b)%i,Z(b)%j,Z(b)%k
end do

end program test_type
