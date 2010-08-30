! this is the driver routine for LT-AEM programs
! this can compute a time-domain solution either on a grid of locations and times,
! at a single point through time (i.e., hydrograph), or at a particle location through
! space and time.  The solution depends on the given properties of circular and 
! elliptical elements.  Wells or points are treated as special case circles, and lines are
! treated as special case ellipses.

program ltaem_main
  use constants, only : DP
  use type_definitions, only : domain, element, circle, ellipse, solution, INVLT, particle
  use file_ops, only : readinput, writeresults
  use inverse_Laplace_Transform, only : invlap => deHoog_invlap, pvalues => deHoog_pvalues
  use particle_integrate, only : rungeKuttaMerson, rungeKutta, fwdEuler
  use solution_mod, only : matrix_solution
  use calc_routines, only : headCalc, velCalc
  use geometry, only : distanceAngleCalcs
  use ellipse_mathieu_init, only : ellipse_init

#ifdef OMP
  use omp_lib, only : omp_get_thread_num
#endif

  implicit none

  ! types or "structs" that organize variables
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  type(particle), allocatable :: part(:)
  integer :: i, j, ierr
  integer :: tnp  ! total-number-p
  integer :: nc, ne     ! #-circles, #-ellipses
  integer :: crow, ccol ! coeff-#-row, coeff-#-col
  integer :: ilogt, iminlogt, imaxlogt     ! indexes related to log10 time
  integer :: lot, hit, lop, hip, lo        ! local hi and lo indices for each log cycle
  integer, allocatable :: nt(:), run(:)    ! #-times-each-log-cycle, 0,1,2,3...
  integer, allocatable :: parnumdt(:)      ! number of dt expected for each particle
  integer, allocatable :: idxmat(:,:)  ! matrix of indices for parallelization
  real(DP), allocatable :: logt(:), tee(:) ! log10-time(numt), Tmax-for-deHoog(num-log-cycles)
  complex(DP), allocatable :: s(:,:)       ! laplace-parameter(2*M-1,num-log-cycles)
  complex(DP) :: calcZ      ! calc-point-complex-coordinates
  character(6) :: elType    ! element-type {CIRCLE,ELLIPS}
  character(20) :: tmpfname
  complex(DP), allocatable :: hp(:),vp(:,:)

  real(DP), parameter :: TMAX_MULT = 2.0_DP  ! traditionally 2.0, but 4.0 could work...

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

  ! is the last value right on a log-cycle boundary?
  if (abs(sol%t(sol%nt) - 10.0**nint(log10(sol%t(sol%nt)))) < 0.001) then
     sol%t(sol%nt) = sol%t(sol%nt) - epsilon(sol%t(sol%nt))
  end if
  
  allocate(run(1:2*sol%m+1), stat=ierr)
  if (ierr /= 0) stop 'ltaem_main.f90 error allocating: run'
  forall(i=0:2*sol%m) run(i+1)=i

  ! compute element geometry from input
  ! read in element heirarchy data from file
  call DistanceAngleCalcs(c,e,bg,dom,sol)

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (sol%calc) then
     if (sol%particle) then   ! particle tracking
        
        ! TODO need to do something better here about particles starting at t=0?
        do i=1, sol%nPart
           if (part(i)%ti < 1.0D-5) then
              part(i)%ti = 1.0D-5
              write(*,'(A,I0,A)') '^^^^ start time for particle ',i,' reset to 1.0E-5 ^^^^'
           end if
        end do
        
        iminlogt = floor(  minval(log10(part(:)%ti)))
        imaxlogt = ceiling(maxval(log10(part(:)%tf)))

        allocate(s(2*sol%m+1,iminlogt:imaxlogt), nt(iminlogt:imaxlogt), &
             & tee(iminlogt:imaxlogt))

        nt(iminlogt:imaxlogt) = 1

        do ilogt = iminlogt, imaxlogt
           ! 0.999 is "most" of a log-cycle
           tee(ilogt) = min(10.0**(ilogt + 0.999), maxval(part(:)%tf))*TMAX_MULT
           s(:,ilogt) = pvalues(tee(ilogt),sol%INVLT)
        end do

        ! to make it possible for particles / contours to share code...
        imaxlogt = imaxlogt + 1
        deallocate(run)
        
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     else   ! contour maps / hydrographs
        allocate(logt(1:sol%nt), stat=ierr)
        if (ierr /= 0) stop 'ltaem_main.f90 error allocating: logt'
        
        lo = 1
        logt(:) = log10(sol%t(:))
        iminlogt =   floor(minval(logt(:)))
        imaxlogt = ceiling(maxval(logt(:)) + epsilon(1.0))

        allocate(s(2*sol%m+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1), stat=ierr)
        if (ierr /= 0) stop 'ltaem_main.f90 error allocating: s,nt,tee'

        do ilogt = iminlogt, imaxlogt-1
           ! number of times falling in this logcycle
           nt(ilogt) = count(logt >= real(ilogt,DP) .and. logt < real(ilogt+1,DP))
           ! T=2*max time in this logcycle
           tee(ilogt) = maxval(sol%t(lo:lo+nt(ilogt)-1))*TMAX_MULT
           s(:,ilogt) = pvalues(tee(ilogt),sol%INVLT)
           lo = lo + nt(ilogt)
        end do
        deallocate(logt,run, stat=ierr)
        if (ierr /= 0) stop 'ltaem_main.f90 error deallocating: logt,run'
     end if

     sol%totalnP = product(shape(s)) ! total number of Laplace parameters across all times
     tnp = sol%totalnP

     print *, 'shape(s)',shape(s),' tnp',tnp

     ! calculate coefficients for each value of Laplace parameter
     ! ** common between particle tracking and contours/hydrographs **
     
     ! initialize Mathieu function matrices
     if (ne > 0) then
        write(*,'(A)') 'Computing Mathieu coefficients ...'
        call ellipse_init(e,bg,s)
     end if

     allocate(idxmat(size(s,dim=1),lbound(s,dim=2):ubound(s,dim=2)))
     idxmat = reshape([(j,j=1,tnp)],shape(s))

     do ilogt = iminlogt,imaxlogt-1
        write(*,'(A,I0,A)') 'log t= 10^(',ilogt,')' 
        if (nt(ilogt) > 0) then
           do j = 1,2*sol%m+1
              write(*,'(I0,1X,2(A,ES10.3),A)') j, '(',real(s(j,ilogt)),',',aimag(s(j,ilogt)),')'
              call matrix_solution(c,e,dom,sol,s(j,ilogt),idxmat(j,ilogt))
           end do
        end if
     end do

     ! only write output if there is at least one matching element
     if(any(e(:)%match) .or. any(c(:)%match)) then

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! save coefficient matrices to file

        if (.not. sol%skipdump) then
           open(unit=77, file=sol%coefffname, status='replace', action='write', iostat=ierr)
           if (ierr /= 0) then
              print *, 'error writing intermediate coefficient matrix to file',sol%coefffname
              stop 000
           end if
           write(77,*) trim(sol%infname)//'.echo' ! file with all the input parameters
           write(77,*) iminlogt,imaxlogt,nc,ne
           write(77,*) nt(:)
           write(77,*) s(:,:)
           write(77,*) tee(:)
           do i = 1,nc
              write(77,*) 'CIRCLE',i,shape(c(i)%coeff)
              write(77,*) c(i)%coeff(:,:)
           end do
           do i = 1,ne
              write(77,*) 'ELLIPS',i,shape(e(i)%coeff)
              write(77,*) e(i)%coeff(:,:)
           end do
           close(77)
           write(*,'(A)') '  <matching finished>  '
        end if
     end if
  else ! do not re-calculate coefficients (this only makes sense if num matching elements > 0)

     if(any(e(:)%match) .or. any(c(:)%match)) then
        open(unit=77, file=sol%coefffname, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
           ! go back and recalc if no restart file
           write(*,'(A)') 'error opening restart file, recalculating...'
           sol%calc = .true.
           goto 111
        end if

        read(77,*) !! TODO not doing anything with input file, but I should
        read(77,*) iminlogt,imaxlogt,nc,ne ! scalars
        allocate(s(2*sol%m+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1),stat=ierr)
        if (ierr /= 0) stop 'ltaem_main.f90 error allocating: s,nt,tee'
        read(77,*) nt(:)
        read(77,*) s(:,:)
        sol%totalnP = product(shape(s))
        tnp = sol%totalnP
        read(77,*) tee(:)
        do i = 1,nc
           read(77,*) elType,j,crow,ccol
           if (elType == 'CIRCLE' .and. i == j) then
              allocate(c(i)%coeff(crow,ccol), stat=ierr)
              if (ierr /= 0) stop 'ltaem_main.f90 error allocating c(i)%coeff'
           else
              stop 'error reading in CIRCLE matching results'
           end if
           read(77,*) c(i)%coeff(:,:)
        end do
        do i = 1,ne
           read(77,*) elType,j,crow,ccol
           if (elType == 'ELLIPS' .and. i == j) then
              allocate(e(i)%coeff(crow,ccol), stat=ierr)
              if (ierr /= 0) stop 'ltaem_main.f90 error allocating e(i)%coeff'
           else
              stop 'error reading in ELLIPS matching results'
           end if
           read(77,*) e(i)%coeff(:,:)
        end do

        ! re-initialize Mathieu function matrices
        if (ne > 0) then
           write(*,'(A)') 'computing Mathieu coefficients ...'
           call ellipse_init(e,bg,s)
        end if

        write(*,'(A)') 'matching results successfully read from file'
     end if

  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(sol%particle) then ! integrate along particle path

     write(*,'(A)') 'compute solution for tracking particles'
     allocate(parnumdt(sol%nPart))

     parnumdt(:) = ceiling((part(:)%tf - part(:)%ti)/part(:)%dt)

     do i = 1,sol%nPart
        ! initialize particles
        allocate(part(i)%r(0:parnumdt(i),5))
        part(i)%r(:,:) = 0.0
        part(i)%id = i
     end do

     ! cycle over particles, integrating each
     ! re-shape matrix s into a vector
     tmpfname = 'calc_part_    .debug'

     !$OMP PARALLEL DO PRIVATE(j) SHARED(part)
     do j = 1, sol%nPart

#ifdef DEBUG
        write(tmpfname(11:14),'(I4.4)') j
        open(unit=77,file=tmpfname,action='write',status='replace')
#endif

        select case (part(j)%int)
        case (1)
           call rungekuttamerson(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case (2)
                 call rungekutta(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case (4)
                   call fwdEuler(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case default
           write(*,'(A,I0,1X,I0)') 'invalid integration code', j, part(j)%int
        end select

#ifdef DEBUG
        close(77)
#endif

     end do
     !$OMP END PARALLEL DO
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(sol%contour) then ! contour output (x,y locations outer product of x,y vectors)

     write(*,'(A)') 'compute solution for plotting contours'
     allocate(sol%h(sol%nx,sol%ny,sol%nt),   hp(tnP), &
          &   sol%v(sol%nx,sol%ny,sol%nt,2), vp(tnP,2), stat=ierr)
     if (ierr /= 0) stop 'ltaem_main.f90 error allocating contour: sol%h,hp,sol%v,vp'

#ifdef DEBUG
     open(unit=303,file='calcloc.debug',status='replace',action='write')
     open(unit=404,file='calcloc.vdebug',status='replace',action='write')
#endif

     !$OMP PARALLEL DO PRIVATE(i,j,calcZ,hp,vp,ilogt,lot,hit,lop,hip) SHARED(sol)
     do j = 1,sol%nx
#ifdef OMP
        write (*,'(A,ES13.5,A,I0)') 'x: ',sol%x(j),'  thr=',OMP_get_thread_num()
#else
        write (*,'(A,ES13.5)') 'x: ',sol%x(j) 
#endif
        do i = 1,sol%ny

           calcZ = cmplx(sol%x(j),sol%y(i),DP)

           !! compute f(p) for all values of p
           !! not just one log cycle of time -- at this location 
           hp(1:tnp) =    headCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
           vp(1:tnp,1:2) = velCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)

           !! invert solutions one log-cycle of t at a time
           do ilogt = iminlogt,imaxlogt-1
              
              !! group of times in current log cycle
              lot = 1 + sum(nt(iminlogt:ilogt-1))
              hit = sum(nt(iminlogt:ilogt))
              
              !! group of Laplace parameters corresponding to this logcycle
              lop = (ilogt - iminlogt)*size(s,dim=1) + 1
              hip = lop + size(s,dim=1) - 1

              sol%h(j,i,lot:hit) =     invlap(sol%t(lot:hit), tee(ilogt), hp(lop:hip), sol%INVLT)
              sol%v(j,i,lot:hit,1:2) = invlap(sol%t(lot:hit), tee(ilogt), vp(lop:hip,1:2), sol%INVLT)
           end do
        end do
     end do
     !$OMP END PARALLEL DO

#ifdef DEBUG
     close(303)
     close(404)
#endif

  else ! hydrograph output (x,y locations are in pairs; e.g. inner product)

     write(*,'(A)') 'compute solution for plotting hydrograph'
     allocate(sol%h(sol%nx,1,sol%nt), hp(tnp), &
          &   sol%v(sol%nx,1,sol%nt,2), vp(tnp,2), stat=ierr)
     if (ierr /= 0) stop 'ltaem_main.f90 error allocating hydrograph: sol%h,hp,sol%v,vp'
     
     do i = 1,sol%nx
        write(*,'(A,2(1X,ES14.7E1))') 'location:',sol%x(i),sol%y(i)

        calcZ = cmplx(sol%x(i),sol%y(i),DP)
        hp(1:tnp) =    headCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
        vp(1:tnp,1:2) = velCalc(calcZ,reshape(s,[tnp]),1,tnp,dom,c,e,bg)
        
        do ilogt = iminlogt,imaxlogt-1
           lot = 1 + sum(nt(iminlogt:ilogt-1))
           hit = sum(nt(iminlogt:ilogt))

           lop = (ilogt - iminlogt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           ! don't need second dimension of results matricies
           sol%h(i,1,lot:hit) =     invlap(sol%t(lot:hit),tee(ilogt),hp(lop:hip),sol%INVLT)
           sol%v(i,1,lot:hit,1:2) = invlap(sol%t(lot:hit),tee(ilogt),vp(lop:hip,1:2),sol%INVLT)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file
  call writeResults(sol,part)

end program ltaem_main

