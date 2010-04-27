program test
implicit none
integer, dimension(5) :: w
integer, dimension(4) :: v
real, dimension(5) :: x
integer :: j

x = [1.0, 2.0, 3.0, 4.0, 5.0]
w = [1, 1, 2, 3, 4]
v = [3, 2, 1, 2]

print *, 'x:',x
print *, 'w:',w
print *, 'v:',v
print *, 'x(w):',x(w)
print *, 'x(w(v)):',x(w(v))

end program test
