program test

integer, parameter :: DP = selected_real_kind(14,300), NUM=5
real(DP), parameter :: PI = 4.0_DP*atan(1.0_DP)

real(DP), dimension(2,NUM) :: x

print *, 'linear:',linspace(-PI/2.0,PI/2.0,NUM)
print *, 'linear no ends:',linspace_no_ends(-PI/2.0,PI/2.0,NUM)
print *, '  '
x = chebspace(-PI/2.0,PI/2.0,NUM)
print *, 'chebyshev: int1',x(1,:)
print *, 'chebyshev: int2',x(2,:)
print *, '  '
x = chebspace(PI/2.0,-PI/2.0,NUM)
print *, 'chebyshev: int3',x(1,:)
print *, 'chebyshev: int4',x(2,:)


contains

  function linspace(lo,hi,num) result(v)
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v

    integer :: i
    real(DP) :: rnum, range

    if(lo >= hi) stop "LINSPACE: lower bound must be less than upper bound."

    rnum = real(num - 1,DP)
    range = hi - lo

    v = (/ ( lo + real(i,DP)*range/rnum, i=0,num-1) /)

  end function linspace
  function linspace_no_ends(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v

    integer :: i
    real(DP) :: rnum, range

    if(lo >= hi) stop "LINSPACE: lower bound must be less than upper bound."

    rnum = real(num + 1,DP)
    range = hi - lo

    v = (/ ( lo + real(1+i,DP)*range/rnum, i=0,num-1) /)

  end function linspace_no_ends

  !! given a low and high theta for the first interval, give
  !! two vectors of points distributed in a chebyshev-style 
  !! spacing over each interval (but covering all of the -pi to +pi 
  !! interval (taking wrapping into account)
  function chebspace(th1,th2,num) result(v)
    use constants, only : DP, PI, TWOPI
    real(DP), intent(in) ::  th1, th2
    integer, intent(in) :: num
    real(DP), dimension(2,num) :: v

    real(DP) :: range1, range2
    real(DP), dimension(num) :: cheb

    cheb(1:num) = 0.5_DP*(1.0_DP - cos(linspace(0.0_DP,PI,num)))

    if(th1 < th2) then
       !! +/-pi cut in second interval
       range1 = th2-th1
       v(1,:) = th1 + range1*cheb(:)
       
       range2 = TWOPI - range1
       !! th2 is now th1
       v(2,:) = th2 + range2*cheb(:)

       !! wrap interval around
       where (v(2,:) > PI)
          v(2,:) = v(2,:) - TWOPI
       end where
    else
       !! +/-pi cut in first interval
       range2 = th1-th2
       v(2,:) = th2 + range2*cheb(:)
       
       range1 = TWOPI - range2
       !! th2 is now th1
       v(1,:) = th1 + range1*cheb(:)

       !! wrap interval around
       where (v(1,:) > PI)
          v(1,:) = v(1,:) - TWOPI
       end where
    end if
  end function chebspace


 end program test
