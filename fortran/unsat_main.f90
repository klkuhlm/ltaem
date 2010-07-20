! this is the driver routine for steady-state unsaturated flow
! this can compute a steady-state solution either on a grid of locations.
! The solution depends on the given properties of circular and 
! elliptical elements.  Points are treated as special case circles, and lines are
! treated as special case ellipses.

program steady_unsat_main
  use constants, only : DP
  use unsat_type_definitions, only : domain, element, circle, ellipse, solution
  use unsat_file_ops, only : readinput, writeresults
  use unsat_solution_mod, only : matrix_solution
  use unsat_calc_routines, only : headCalc, velCalc
  use unsat_geometry, only : distanceAngleCalcs
  use mathieu_functions, only : mathieu_init
  implicit none

  ! types or "structs" that organize variables
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  integer :: i, j, ierr
  integer :: nc, ne    ! #-circles, #-ellipses
  complex(DP) :: calcZ ! calc-point-complex-coordinates

  intrinsic :: get_command_argument

  call get_command_argument(1,sol%inFName)
  if (len_trim(sol%infname) == 0) then
     write(*,'(A)') 'no command-line filename supplied, using default input file: input.in'
     sol%infname = 'input.in'
  end if

  ! read in data, initialize variables, allocate major structs
  call readInput(sol,dom,bg,c,e)
  nc = size(c,dim=1)
  ne = size(e,dim=1)
  
  ! compute element geometry from input
  call DistanceAngleCalcs(c,e,bg,dom,sol)
   
  ! initialize Mathieu function matrices
  if (ne > 0) then
     write(*,'(A)') 'Computing Mathieu coefficients ...'
     
     ! compute background
     bg%mat = mathieu_init(cmplx(bg%alpha/2.0,0.0,DP)**2, &
          & MM=bg%ms,CUTOFF=bg%cutoff)
     
     ! allocate/initialize each element for each value of p
     do j = 1, size(e,1)
        if (e(j)%ibnd == 0) then
           e(j)%mat = mathieu_init(cmplx(e(j)%alpha/2.0,0.0,DP)**2, &
                & MM=e(j)%MS,CUTOFF=e(j)%cutoff)
        end if
     end do
     
  end if

  call matrix_solution(c,e,dom)
 
  print *, 'coefficients computed, calculating solution'

  allocate(sol%h(sol%nx,sol%ny), sol%v(sol%nx,sol%ny,2), &
       & sol%smphi(sol%nx,sol%ny), sol%lgPHI(sol%nx,sol%ny), stat=ierr)
  if (ierr /= 0) stop 'unsat_main.f90 error allocating contour: sol%{h,v,smphi,lgphi}'

#ifdef DEBUG
  open(unit=303,file='unsat_calcloc.debug',status='replace',action='write')
  open(unit=404,file='unsat_calcloc.vdebug',status='replace',action='write')
#endif

  do j = 1,sol%nx
     print *, 'x:',sol%x(j)
     do i = 1,sol%ny
        calcZ = cmplx(sol%x(j),sol%y(i),DP)
        
        !! compute f(p) for all values of p
        !! not just one log cycle of time -- at this location 
        sol%h(j,i) =    real(headCalc(calcZ,dom,c,e,bg))
        sol%v(j,i,1:2) = real(velCalc(calcZ,dom,c,e,bg))
        
     end do
  end do

#ifdef DEBUG
  close(303)
  close(404)
#endif

  ! downward Z is -y
  sol%smphi = sol%h

  ! PHI = phi - Z
  sol%lgPHI = sol%h - spread(exp(-sol%y),1,sol%nx)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file
  call writeResults(sol)

  deallocate(c,e, stat=ierr)
  if (ierr /= 0) print *, 'WARNING: problem deallocating: c,e,sol%{h,v} at end of program'

end program steady_unsat_main

