program test 
implicit none

integer :: i,j
real, dimension(5,-2:3) :: x
real, dimension(15) :: y
integer, dimension(2) :: ix

forall(i=1:5,j=-2:3)
   x(i,j) = real(i+10*j)
end forall

ix = shape(x)
print *, ix
y = reshape(x,[ix(1)*ix(2)])

write(*,*) 'x:'
do i=1,5
   write(*,*) ('(',i,',',j,')',x(i,j),j=-2,3)
end do

write(*,*) 'y:'
write(*,*) ('(',i,')',y(i),i=1,15)


end program test
