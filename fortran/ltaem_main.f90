!
! Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

! this is the driver routine for LT-AEM programs, which
! can compute a time-domain solution either on a grid of locations and times,
! at a single point through time (i.e., time series), or at a particle location through
! space and time.  The solution depends on the given properties of circular and
! elliptical elements.  Wells or points are treated as special case circles, and lines are
! treated as special case ellipses.

program ltaem_main
  use constants, only : DP
  use type_definitions, only : domain, element, circle, ellipse, solution, INVLT, particle
  use file_ops, only : readinput, writeresults
  use inverse_Laplace_Transform, only : L => deHoog_invlap, pvalues => deHoog_pvalues
  use particle_integrate, only : rungeKuttaMerson, rungeKutta, fwdEuler, analytic
  use solution_mod, only : matrix_solution
  use calc_routines, only : headCalc, velCalc
  use geometry, only : distanceAngleCalcs
  use ellipse_mathieu_init, only : ellipse_init

  !!  for parallel execution using OpenMP
  !$  use omp_lib, only : omp_get_thread_num, omp_get_num_procs, omp_get_max_threads

  implicit none
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  type(particle), allocatable :: part(:)
  integer :: i, j, ierr
  integer :: tnp                              ! total-number-p (product of dimensions)
  integer :: nc, ne                           ! #-circles, #-ellipses
  integer :: crow, ccol                       ! coeff-#-row, coeff-#-col
  integer :: lt, minlt, maxlt                 ! indexes related to log10 time
  integer :: lot, hit, lop, hip, lo           ! local hi and lo indices for each log cycle
  integer, allocatable :: nt(:)               ! #-times-each-log-cycle
  integer, allocatable :: parnumdt(:)         ! number of dt expected for each particle
  integer, allocatable :: idxmat(:,:)         ! matrix of indices for parallelization
  real(DP), allocatable :: logt(:), tee(:)    ! log10-time(numt), Tmax-for-deHoog(num-log-cycles)
  complex(DP), allocatable :: s(:,:), stmp(:) ! Laplace-parameter(2*M-1,num-log-cycles)
  complex(DP) :: calcZ                        ! calc-point-complex-coordinates
  character(6) :: elType                      ! element-type {CIRCLE,ELLIPS}
  complex(DP), allocatable :: hp(:), vp(:,:)  ! Laplace-space head and velocity vectors

  ! some constants that shouldn't really need to be adjusted too often
  real(DP), parameter :: EARLIEST_PARTICLE = 1.0E-5, MOST_LOGT = 0.999
  real(DP), parameter :: TMAX_MULT = 2.0_DP  ! traditionally 2.0, but 4.0 could work (Mark Bakker)...

  intrinsic :: get_command_argument
  call get_command_argument(1,sol%inFName)
  if (len_trim(sol%infname) == 0) then
     write(*,'(A)') 'no command-line filename supplied, using default input file: input.in'
     sol%infname = 'input.in'
  end if

  !! some parallel-related statistics
  !$ write(*,'(2(A,I0))') 'OpenMP num procs:',omp_get_num_procs(), ' OpenMP max threads:',omp_get_max_threads()

  ! read in data, initialize variables, allocate major structs
  call readInput(sol,dom,bg,c,e,part)
  nc = size(c,dim=1)
  ne = size(e,dim=1)

  ! compute element geometry from input
  ! read in element hierarchy data from file
  call DistanceAngleCalcs(c,e,bg,dom,sol)

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (sol%calc) then
     if (sol%particle) then   ! particle tracking

        where (part(:)%ti < EARLIEST_PARTICLE)
           part(:)%ti = EARLIEST_PARTICLE
        end where

        minlt = floor(  minval(log10(part(:)%ti)))
        maxlt = ceiling(maxval(log10(part(:)%tf)))

        allocate(s(2*sol%m+1,minlt:maxlt), nt(minlt:maxlt), tee(minlt:maxlt))

        nt(minlt:maxlt) = 1

        forall(lt = minlt:maxlt)
           tee(lt) = min(10.0**(lt + MOST_LOGT), maxval(part(:)%tf))*TMAX_MULT
           s(:,lt) = pvalues(tee(lt),sol%INVLT)
        end forall

        ! to make it possible for particles / contours to share code...
        maxlt = maxlt + 1

     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     else   ! contour maps / time-series plots
        allocate(logt(1:sol%nt))

        lo = 1
        logt(:) = log10(sol%t(:))
        minlt =   floor(minval(logt(:)))
        ! add epsilon to ensure is bumped up to next log cycle if on fence
        maxlt = ceiling(maxval(logt(:)) + epsilon(1.0))

        allocate(s(2*sol%m+1,minlt:maxlt-1), nt(minlt:maxlt-1), tee(minlt:maxlt-1))

        do lt = minlt, maxlt-1
           ! number of times falling in this logcycle
           nt(lt) = count(logt >= real(lt,DP) .and. logt < real(lt+1,DP))
           ! T=2*max time in this logcycle
           tee(lt) = maxval(sol%t(lo:lo+nt(lt)-1))*TMAX_MULT
           s(:,lt) = pvalues(tee(lt),sol%INVLT)
           lo = lo + nt(lt)
        end do
        deallocate(logt)
     end if

     sol%totalnP = size(s) ! total number of Laplace parameters across all times
     tnp = sol%totalnP

     ! calculate coefficients for each value of Laplace parameter
     ! ** common between particle tracking and contours/time-series plots **

     ! initialize Mathieu function matrices
     if (ne > 0) then
        write(*,'(A)') 'Computing Mathieu coefficients ...'
        call ellipse_init(e,bg,s)
     end if

     ! create an index the same shape as s (which has non-standard lower bound on dim 2)
     allocate(idxmat(size(s,dim=1),lbound(s,dim=2):ubound(s,dim=2)))
     idxmat = reshape([(j,j=1,tnp)],shape(s))

     ! when in parallel, call once to initialize the results matrices in solution.f90
     !$ call matrix_solution(c,e,dom,sol,s(1,minlt),idxmat(1,minlt))
     !$ write(*,'(A/)') 'solution initialization...'

     !$ write(*,'(A)',advance="no") 'thr'
     write(*,'(A)') ' logt  idx            p'

     !$OMP PARALLEL DO PRIVATE(lt,j)
     do lt = minlt,maxlt-1
        if (nt(lt) > 0) then
           do j = 1,2*sol%m+1
              !$ write (*,'(I2,1X)',advance="no") OMP_get_thread_num()
              write(*,'(I4,1X,I4,1X,2(A,ES10.3),A)') lt,j, '(',&
                   & real(s(j,lt)),',',aimag(s(j,lt)),')'
              call matrix_solution(c,e,dom,sol,s(j,lt),idxmat(j,lt))
           end do
        end if
     end do
     !$OMP END PARALLEL DO

     ! only write output if there is at least one matching element
     if(any(e(:)%match) .or. any(c(:)%match)) then

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! save coefficient matrices to file (if sol%output < 100)

        if (.not. sol%skipdump) then
           open(unit=77, file=sol%coefffname, status='replace', action='write', iostat=ierr)
           if (ierr /= 0) then
              write(*,'(2A)') 'WARNING: error opening intermediate save file ',trim(sol%coefffname), &
                   & ' continuing without saving results'
           else
              write(77,'(A)') trim(sol%infname)//'.echo' ! file with all the input parameters
              write(77,'(4(I0,1X))') minlt,maxlt,nc,ne
              write(77,*) nt(:)
              write(77,*) s(:,:)
              write(77,*) tee(:)
              do i = 1,nc
                 write(77,'(A,3(I0,1X))') 'CIRCLE ',i,shape(c(i)%coeff)
                 write(77,*) c(i)%coeff(:,:)
              end do
              do i = 1,ne
                 write(77,'(A,3(I0,1X))') 'ELLIPS ',i,shape(e(i)%coeff)
                 write(77,*) e(i)%coeff(:,:)
              end do
              close(77)
           end if
           write(*,'(A)') '  <matching finished>  '
        end if
     end if
  else ! do not re-calculate coefficients (this only makes sense if num matching elements > 0)

     if(any(e(:)%match) .or. any(c(:)%match)) then
        open(unit=77, file=sol%coefffname, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
           ! go back and recalculate if no restart file
           write(*,'(A)') 'ERROR: cannot opening restart file, recalculating...'
           sol%calc = .true.
           goto 111
        end if

        read(77,*) !! TODO not doing anything with input file, but should I check inputs are same?
        read(77,*) minlt,maxlt,nc,ne ! scalars
        allocate(s(2*sol%m+1,minlt:maxlt-1), nt(minlt:maxlt-1), tee(minlt:maxlt-1))
        read(77,*) nt(:)
        read(77,*) s(:,:)
        sol%totalnP = product(shape(s))
        tnp = sol%totalnP
        read(77,*) tee(:)
        do i = 1,nc
           read(77,*) elType,j,crow,ccol
           if (elType == 'CIRCLE' .and. i == j) then
              allocate(c(i)%coeff(crow,ccol))
           else
              write(*,'(A)') 'ERROR reading in CIRCLE matching results, recalculating...'
              sol%calc = .true.
              goto 111
           end if
           read(77,*) c(i)%coeff(:,:)
        end do
        do i = 1,ne
           read(77,*) elType,j,crow,ccol
           if (elType == 'ELLIPS' .and. i == j) then
              allocate(e(i)%coeff(crow,ccol))
           else
              write(*,'(A)') 'ERROR reading in ELLIPS matching results, recalculating...'
              sol%calc = .true.
              goto 111
           end if
           read(77,*) e(i)%coeff(:,:)
        end do

        ! re-initialize Mathieu function matrices
        if (ne > 0) then
           write(*,'(A)') 're-computing Mathieu coefficients ...'
           call ellipse_init(e,bg,s)
        end if

        write(*,'(A)') ' <matching results successfully re-read from file> '
     end if

  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(sol%particle) then ! integrate along particle path

     write(*,'(A)') 'compute solution for tracking particles'
     allocate(parnumdt(sol%nPart))

     ! assuming constant time steps, determine max # steps required
     parnumdt(:) = ceiling((part(:)%tf - part(:)%ti)/part(:)%dt)

     do i = 1,sol%nPart
        ! initialize particles
        allocate(part(i)%r(0:parnumdt(i),5))
        part(i)%r(:,:) = 0.0
        part(i)%id = i
     end do

     !$OMP PARALLEL DO SHARED(part,sol)
     do j = 1, sol%nPart

        select case (part(j)%int)
        case (1)
           call rungeKuttaMerson(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case (2)
           call rungeKutta(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case (3)
           call analytic(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        case (4)
           call fwdEuler(s,tee,c,e,bg,sol,dom,part(j),lbound(tee,dim=1))
        end select
        ! invalid integration code checked in ltaem-io routine

     end do
     !$OMP END PARALLEL DO

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(sol%contour) then ! contour output (x,y locations outer product of x,y vectors)

     write(*,'(A)') 'compute solution for plotting contours'
     allocate(sol%h(sol%nx,sol%ny,sol%nt),   hp(tnP), &
          &   sol%v(sol%nx,sol%ny,sol%nt,2), vp(tnP,2), stmp(tnp))
     if (sol%deriv) then
        allocate(sol%dh(sol%nx,sol%ny,sol%nt))
     end if

     !$OMP PARALLEL DO PRIVATE(calcZ,hp,vp,lot,hit,lop,hip,stmp) SHARED(sol)
     do j = 1,sol%nx
        !$ write (*,'(I0,1X)',advance="no") OMP_get_thread_num()
        write (*,'(A,ES13.5)') 'x: ',sol%x(j)
        do i = 1,sol%ny

           calcZ = cmplx(sol%x(j),sol%y(i),DP)
           stmp(1:tnp) = reshape(s,[tnp])

           !! compute f(p) for all values of p at this location
           hp(1:tnp) =    headCalc(calcZ,stmp,1,tnp,dom,c,e,bg)
           vp(1:tnp,1:2) = velCalc(calcZ,stmp,1,tnp,dom,c,e,bg)

           !! invert solutions one log-cycle of t at a time
           do lt = minlt,maxlt-1

              !! group of times in current log cycle
              lot = 1 + sum(nt(minlt:lt-1))
              hit = sum(nt(minlt:lt))

              !! group of Laplace parameters corresponding to this logcycle
              lop = (lt - minlt)*size(s,dim=1) + 1
              hip = lop + size(s,dim=1) - 1

              sol%h(j,i,lot:hit) =     L(sol%t(lot:hit), tee(lt), hp(lop:hip), sol%INVLT)
              if (sol%deriv) then
                 sol%dh(j,i,lot:hit) = L(sol%t(lot:hit), tee(lt), hp(lop:hip)*stmp(lop:hip), sol%INVLT)*sol%t(lot:hit)
              end if
              sol%v(j,i,lot:hit,1:2) = L(sol%t(lot:hit), tee(lt), vp(lop:hip,1:2), sol%INVLT)
           end do
        end do
     end do
     !$OMP END PARALLEL DO

  else ! time-series output (x,y locations are in pairs; e.g. inner product)

     write(*,'(A)') 'compute solution for plotting time series'
     allocate(sol%h(sol%nx,1,sol%nt), hp(tnp), sol%v(sol%nx,1,sol%nt,2), vp(tnp,2), stmp(tnp))
     if (sol%deriv) then
        allocate(sol%dh(sol%nx,1,sol%nt))
     end if

     !$OMP PARALLEL DO PRIVATE(calcZ,hp,vp,lot,hit,lop,hip,stmp) SHARED(sol)
     do i = 1,sol%nx
        write(*,'(A,2(3X,ES14.7E1))') sol%obsname(i),sol%xshift+sol%x(i),sol%yshift+sol%y(i)

        stmp(1:tnp) = reshape(s,[tnp])

        calcZ = cmplx(sol%x(i),sol%y(i),DP)
        hp(1:tnp) =    headCalc(calcZ,stmp,1,tnp,dom,c,e,bg)
        vp(1:tnp,1:2) = velCalc(calcZ,stmp,1,tnp,dom,c,e,bg)

        do lt = minlt,maxlt-1
           lot = 1 + sum(nt(minlt:lt-1))
           hit = sum(nt(minlt:lt))

           lop = (lt - minlt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           ! don't need second dimension of results matrices
           sol%h(i,1,lot:hit) =     L(sol%t(lot:hit),tee(lt),hp(lop:hip),sol%INVLT)
           if (sol%deriv) then
              sol%dh(i,1,lot:hit) = L(sol%t(lot:hit),tee(lt),hp(lop:hip)*stmp(lop:hip),sol%INVLT)*sol%t(lot:hit)
           end if
           sol%v(i,1,lot:hit,1:2) = L(sol%t(lot:hit),tee(lt),vp(lop:hip,1:2),sol%INVLT)
        end do
     end do
     !$OMP END PARALLEL DO

  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file
  call writeResults(sol,part)

end program ltaem_main

