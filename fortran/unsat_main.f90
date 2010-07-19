! this is the driver routine for steady-state unsaturated flow
! this can compute a steady-state solution either on a grid of locations and times,
! at a single point through time (i.e., hydrograph), or at a particle location through
! space and time.  The solution depends on the given properties of circular and 
! elliptical elements.  Wells or points are treated as special case circles, and lines are
! treated as special case ellipses.

program steady_unsat_main
  use constants, only : DP
  use unsat_type_definitions, only : domain, element, circle, ellipse, solution
  use file_ops, only : readinput, writeresults
  use solution_mod, only : matrix_solution
  use calc_routines, only : headCalc, velCalc
  use geometry, only : distanceAngleCalcs
  implicit none

  ! types or "structs" that organize variables
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  integer :: i, j, ierr
  integer :: nc, ne     ! #-circles, #-ellipses
  integer :: crow, ccol ! coeff-#-row, coeff-#-col
  complex(DP) :: calcZ      ! calc-point-complex-coordinates
  character(6) :: elType    ! element-type {CIRCLE,ELLIPS}

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
  ! read in element heirarchy data from file
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

  call unsat_matrix_solution(c,e,dom,sol)
 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(sol%contour) then ! contour output (x,y locations outer product of x,y vectors)

     write(*,'(A)') 'compute solution for plotting contours'
     allocate(sol%h(sol%nx,sol%ny,sol%nt),   sol%hp(tnP), &
          &   sol%v(sol%nx,sol%ny,sol%nt,2), sol%vp(tnP,2), stat=ierr)
     if (ierr /= 0) stop 'ltaem_main.f90 error allocating contour: sol%h,sol%hp,sol%v,sol%vp'

#ifdef DEBUG
     open(unit=303,file='calcloc.debug',status='replace',action='write')
     open(unit=404,file='calcloc.vdebug',status='replace',action='write')
#endif

     do j = 1,sol%nx
        write (*,'(A,ES12.4E1)') 'x: ',sol%x(j)

        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(sol)
        do i = 1,sol%ny

           calcZ = cmplx(sol%x(j),sol%y(i),DP)

           !! compute f(p) for all values of p
           !! not just one log cycle of time -- at this location 
           sol%hp(1:tnp) =    headCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
           sol%vp(1:tnp,1:2) = velCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)

           !! invert solutions one log-cycle of t at a time
           do ilogt = iminlogt,imaxlogt-1
              
              !! group of times in current log cycle
              lot = 1 + sum(nt(iminlogt:ilogt-1))
              hit = sum(nt(iminlogt:ilogt))
              
              !! group of Laplace parameters corresponding to this logcycle
              lop = (ilogt - iminlogt)*size(s,dim=1) + 1
              hip = lop + size(s,dim=1) - 1

              sol%h(j,i,lot:hit) =     invlap(sol%t(lot:hit), tee(ilogt), &
                   & sol%hp(lop:hip), sol%INVLT)
              sol%v(j,i,lot:hit,1:2) = invlap(sol%t(lot:hit), tee(ilogt), &
                   & sol%vp(lop:hip,1:2), sol%INVLT)
           end do
        end do
     !$OMP END PARALLEL DO
     end do

#ifdef DEBUG
     close(303)
     close(404)
#endif

  else ! hydrograph output (x,y locations are in pairs; e.g. inner product)

     write(*,'(A)') 'compute solution for plotting hydrograph'

     allocate(sol%h(sol%nx,1,sol%nt), sol%hp(tnp), &
          &   sol%v(sol%nx,1,sol%nt,2), sol%vp(tnp,2), stat=ierr)
     if (ierr /= 0) stop 'ltaem_main.f90 error allocating hydrograph: sol%h,sol%hp,sol%v,sol%vp'
     
     do i = 1,sol%nx
        write(*,'(A,2(1X,ES14.7E1))') 'location:',sol%x(i),sol%y(i)

        calcZ = cmplx(sol%x(j),sol%y(i),DP)
        sol%hp(1:tnp) =  headCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
        sol%vp(1:tnp,:) = velCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
        
        do ilogt = iminlogt,imaxlogt-1
           lot = 1 + sum(nt(iminlogt:ilogt-1))
           hit = sum(nt(iminlogt:ilogt))

           lop = (ilogt - iminlogt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           ! don't need second dimension of results matricies
           sol%h(i,1,lot:hit) =   invlap(sol%t(lot:hit),tee(ilogt),sol%hp(lop:hip),sol%INVLT)
           sol%v(i,1,lot:hit,:) = invlap(sol%t(lot:hit),tee(ilogt),sol%vp(lop:hip,:),sol%INVLT)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file
  call writeResults(sol,part)

  deallocate(sol%h,sol%hp,sol%v,sol%vp, stat=ierr)
  if (ierr /= 0) print *, 'WARNING: problem deallocating: sol%{h,p,v,vp} at end of program'

  deallocate(c,e,part,s,tee,nt, stat=ierr)
  if (ierr /= 0) print *, 'WARNING: problem deallocating: c,e,part,s,tee,logt,nt at end of program'

end program steady_unsat_main

