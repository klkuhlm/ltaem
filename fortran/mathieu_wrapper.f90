program mathieu_wrapper

! this program is a thin wrapper around the mathieu_function library.  It reads in 
! an input file, computes the requested mathieu functions and writes a file of output

! this program reads the input format used by the Alhhargan driver program, and produces 
! comparable output to simplify testing and therfore is only setup to handle real q<0

use constants, only : DP, TWOPI
use utility, only : linspace
use mathieu_functions
implicit none

integer :: j,k,n,n0,kind
integer, parameter :: nsteps = 200
real(DP) :: radRange,angRange,h,q
real(DP), dimension(0:nsteps) :: arg
complex(DP), allocatable :: fcn(:,:)
character(1) :: typ,r
character(19) :: fmt

type(mathieu) :: mat

write(*,'(A)',advance='no') '(e)ven or (o)dd: '
read(*,'(A)')  typ
write(*,'(A)',advance='no') 'max order n: '
read(*,*)  n
write(*,'(A)',advance='no') '(c)icumferential or (r)adial: '
read(*,'(A)')  r
write(*,'(A)',advance='no') 'kind: First(1) or Second(2): '
read(*,*)  kind
write(*,'(A)') 'real part of Mathieu parameter, h: '
read(*,*) h

q = h**2/4.0_DP

radRange = 4.0_DP
angRange = TWOPI

if (typ == 'o') then
   n0 = 1
else
   n0 = 0
end if

allocate(fcn(n0:n, 0:nsteps))
mat = mathieu_init(cmplx(q,0.0,DP),MM=60)

if (r == 'r') then
   ! don't want the arg=0 step, since radial functions blow up there as q->0
   arg(0:nsteps) = linspace(0.0_DP,radRange,nsteps+1)
   
   if (kind == 1) then
      if (typ == 'e') then
         fcn(n0:n,:) = Ie(mat,[(j,j=n0,n)],arg(0:nsteps))
      else
         fcn(n0:n,:) = Io(mat,[(j,j=n0,n)],arg(0:nsteps))
      end if
   else
      if (typ == 'e') then
         fcn(n0:n,:) = Ke(mat,[(j,j=n0,n)],arg(0:nsteps))
      else
         fcn(n0:n,:) = Ko(mat,[(j,j=n0,n)],arg(0:nsteps))
      end if
   end if
else
   arg(0:nsteps) = linspace(0.0_DP,angRange,nsteps+1)

   ! second-kind angular functions not needed or handled
   if (typ == 'e') then
      fcn(n0:n,:) = ce(mat,[(j,j=n0,n)],arg(0:nsteps))
   else
      fcn(n0:n,:) = se(mat,[(j,j=n0,n)],arg(0:nsteps))
   end if
end if

fmt = '(XXX(ES22.14E3,1X))'
write(fmt(2:4),'(I3.3)') n-n0+2

! output results to file
do j=0,nsteps
   write(*,fmt) arg(j),(real(fcn(k,j)),k=n0,n)
end do

end program mathieu_wrapper
