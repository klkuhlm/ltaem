! this is the driver routine for steady-state unsaturated flow
! this can compute a steady-state solution either on a grid of locations and times,
! at a single point through time (i.e., hydrograph), or at a particle location through
! space and time.  The solution depends on the given properties of circular and 
! elliptical elements.  Wells or points are treated as special case circles, and lines are
! treated as special case ellipses.

program steady_unsat_main
  use constants, only : DP
  use unsat_type_definitions, only : domain, element, circle, ellipse, solution
  use unsat_file_ops, only : readinput, writeresults
  use unsat_solution_mod, only : matrix_solution
  use unsat_calc_routines, only : headCalc, velCalc
  use unsat_geometry, only : distanceAngleCalcs
  implicit none

  ! types or "structs" that organize variables
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  integer :: i, j, ierr
  integer :: nc, ne     ! #-circles, #-ellipses
  complex(DP) :: calcZ      ! calc-point-complex-coordinates

  intrinsic :: get_command_argument

  call get_command_argument(1,sol%inFName)
  if (len_trim(sol%infname) == 0) then
     write(*,'(A)') 'no command-line filename supplied, using default input file: input.in'
     sol%infname = 'input.in'
  end if

  ! read in data, initialize variables, allocate major structs
  call readInput(sol,dom,bg,c,e,part)
  nc = size(c,dim=1)
  ne = size(e,dim=1)
  
  ! compute element geometry from input
  call DistanceAngleCalcs(c,e,bg,dom,sol)
   
  ! initialize Mathieu function matrices
  if (ne > 0) then
     write(*,'(A)') 'Computing Mathieu coefficients ...'
     
     ! compute background
     bg%mat = mathieu_init((bg%alpha/2.0)**2,MM=bg%ms)
     
     ! allocate/initialize each element for each value of p
     do j = 1, size(e,1)
        if (e(j)%ibnd == 0 .or. e(j)%calcin) then
           e(j)%mat = mathieu_init((e(j)%alpha/2.0)**2,MM=e(j)%MS)
        end if
     end do
     
  end if

  call unsat_matrix_solution(c,e,dom)
 
  allocate(sol%h(sol%nx,sol%ny), sol%v(sol%nx,sol%ny,2), stat=ierr)
  if (ierr /= 0) stop 'ltaem_main.f90 error allocating contour: sol%h,sol%hp,sol%v,sol%vp'

#ifdef DEBUG
  open(unit=303,file='calcloc.debug',status='replace',action='write')
  open(unit=404,file='calcloc.vdebug',status='replace',action='write')
#endif

  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(sol,dom,c,e,bg)
  do j = 1,sol%nx
     do i = 1,sol%ny
        calcZ = cmplx(sol%x(j),sol%y(i),DP)
        
        !! compute f(p) for all values of p
        !! not just one log cycle of time -- at this location 
        sol%h(j,i) =    headCalc(calcZ,dom,c,e,bg)
        sol%v(j,i,1:2) = velCalc(calcZ,dom,c,e,bg)
        
     end do
  end do
  !$OMP END PARALLEL DO

#ifdef DEBUG
  close(303)
  close(404)
#endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file
  call writeResults(sol)

  deallocate(c,e,sol%h,sol%v, stat=ierr)
  if (ierr /= 0) print *, 'WARNING: problem deallocating: c,e,sol%{h,v} at end of program'

end program steady_unsat_main

