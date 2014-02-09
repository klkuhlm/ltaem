!
! Copyright (c) 2011-2014 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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
  use type_definitions, only : domain, element, matching, circle, ellipse, solution, INVLT, particle
  use file_ops, only : readinput, writeresults, read_coeff, dump_coeff
  use inverse_Laplace_Transform, only : L => deHoog_invlap, pvalues => deHoog_pvalues
  use particle_integrate, only : rungeKuttaMerson, rungeKutta, fwdEuler, analytic
  use solution_mod, only : matrix_solution
  use calc_routines, only : headCalc, velCalc, elementFlowrate
  use geometry, only : distanceAngleCalcs
  use ellipse_mathieu_init, only : ellipse_init

  use, intrinsic :: iso_fortran_env, only : stdout => output_unit

  !!  for parallel execution using OpenMP
  !$  use omp_lib, only : omp_get_thread_num, omp_get_num_procs, omp_get_max_threads

  implicit none
  type(domain)  :: dom
  type(element) :: bg
  type(circle),  allocatable :: c(:)
  type(ellipse), allocatable :: e(:)
  type(solution) :: sol
  type(particle), allocatable :: part(:)

  integer :: i, j, k
  integer :: tnp                              ! total # p (product of dimensions)
  integer :: nc, ne                           ! #circles, #ellipses
  integer :: lt, minlt, maxlt                 ! indexes related to log10 time
  integer :: lot, hit, lop, hip, lo           ! local hi and lo indices for log cycles
  integer, allocatable :: nt(:)               ! # times per log cycle
  integer, allocatable :: parnumdt(:)         ! number of dt expected for each particle
  integer, allocatable :: idxmat(:,:)         ! matrix of indices for parallelization
  real(DP), allocatable :: logt(:), tee(:)    ! log10 time, deHoog T
  complex(DP), allocatable :: s(:,:), stmp(:) ! Laplace parameter (different shapes)
  complex(DP) :: calcZ                        ! calculation point
  complex(DP), allocatable :: hp(:), vp(:,:)  ! Laplace-space head and velocity vectors
  complex(DP), allocatable :: qp(:,:) ! Laplace-space flowrate into/out of element
  logical :: fail

  ! constants that shouldn't be adjusted too often
  real(DP), parameter :: MOST_LOGT = 0.999, TMAX_MULT = 2.0_DP  

  intrinsic :: get_command_argument
  call get_command_argument(1,sol%inFName)
  if (len_trim(sol%infname) == 0) then
     write(stdout,'(A)') 'no command-line filename supplied, '&
          &//'using default input file: input.in'
     sol%infname = 'input.in'
  end if

  ! some parallel-related statistics ("!$" is special OMP directive)
  !$ write(stdout,'(2(A,I0))') 'OpenMP num procs:',omp_get_num_procs(), ' OpenMP max threads:', omp_get_max_threads()

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

        logcycles: do lt = minlt, maxlt-1
           ! number of times falling in this logcycle
           nt(lt) = count(logt >= real(lt,DP) .and. logt < real(lt+1,DP))
           if (nt(lt) == 0) then
              cycle logcycles
           end if
           
           ! T=2*max time in this logcycle
           !!print *, 'DBG:',lt,lo,nt(lt),lo+nt(lt)-1,maxval(sol%t(lo:lo+nt(lt)-1))
           tee(lt) = maxval(sol%t(lo:lo+nt(lt)-1))*TMAX_MULT
           s(:,lt) = pvalues(tee(lt),sol%INVLT)
           lo = lo + nt(lt)
        end do logcycles
        deallocate(logt)
     end if

     sol%totalnP = size(s) ! total number of Laplace parameters across all times
     tnp = sol%totalnP

     ! calculate coefficients for each value of Laplace parameter
     ! ** common between particle tracking and contours/time-series plots **

     ! initialize Mathieu function matrices
     if (ne > 0) then
        write(stdout,'(A)') 'Computing Mathieu coefficients ...'
        call ellipse_init(e,bg,s)
     end if

     ! create an index the same shape as s (which has non-standard lower bound on dim 2)
     allocate(idxmat(size(s,dim=1),lbound(s,dim=2):ubound(s,dim=2)))
     idxmat = reshape([(j,j=1,tnp)],shape(s))

     ! when in parallel, call once to initialize the results matrices in solution.f90
     !$ call matrix_solution(c,e,dom,sol,s(1,minlt),idxmat(1,minlt))
     !$ write(stdout,'(A/)') 'parallel solution initialization...'

     !$ write(stdout,'(A)',advance="no") 'thr'
     write(stdout,'(A)') ' logt  idx            p'

     !$OMP PARALLEL DO PRIVATE(lt,j)
     do lt = minlt,maxlt-1
        if (nt(lt) > 0) then
           do j = 1,2*sol%m+1
              !$ write (*,'(I2,1X)',advance="no") OMP_get_thread_num()
              write(stdout,"(I4,1X,I4,1X,'(',ES10.3,',',ES10.3,')')") lt,j,s(j,lt)
              call matrix_solution(c,e,dom,sol,s(j,lt),idxmat(j,lt))
           end do
        end if
     end do
     !$OMP END PARALLEL DO
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! save coefficient matrices to file (if sol%output < 100)
     
     if (.not. sol%skipdump) then
        call dump_coeff(sol,c,e,nt,s,tee,minlt,maxlt)
        write(stdout,'(A)') '  <coefficents dumped, matching finished>  '
     end if

  else   
     ! do not re-calculate coefficients  
     call read_coeff(sol,bg,c,e,nt,s,tee,minlt,maxlt,fail)
     tnp = sol%totalnP
     if (fail) goto 111
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(sol%particle) then ! integrate along particle path

     write(stdout,'(A)') 'compute solution for tracking particles'
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

     end do
     !$OMP END PARALLEL DO

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(sol%contour) then ! x,y locations outer product of x,y vectors

     write(stdout,'(A)') 'compute solution for plotting contours'
     allocate(sol%h(sol%nt,sol%nx,sol%ny),   hp(tnP), &
          &   sol%v(sol%nt,2,sol%nx,sol%ny), vp(tnP,2), stmp(tnP))
     if (sol%deriv) then
        allocate(sol%dh(sol%nt,sol%nx,sol%ny))
     end if

     stmp(1:tnP) = reshape(s,[tnP])

     ! TODO: this block of code is almost repeated below. Refactor into function?
     if (sol%Qcalc) then
        write(stdout,'(A)',advance='no') 'computing element boundary flowrates: '
        allocate(sol%Q(sol%nt,nc+ne),qp(tnP,nc+ne))
        if (sol%deriv) then
           allocate(sol%dQ(sol%nt,nc+ne))
        end if
        do j = 1,nc
           write(stdout,'(A,I0,1X)',advance='no') 'C',j
           qp(:,j) = elementFlowrate(c(j)%matching,stmp,1,tnP,dom,c,e,bg)
        end do
        do j = 1,ne
           write(stdout,'(A,I0,1X)',advance='no') 'E',j
           qp(:,nc+j) = elementFlowrate(e(j)%matching,stmp,1,tnP,dom,c,e,bg)
        end do
        write(stdout,'(/)')

        !$OMP PARALLEL DO PRIVATE(lt,lot,hit,lop,hip,k) shared(sol,qp)
        do lt = minlt,maxlt-1
           lot = 1 + sum(nt(minlt:lt-1))
           hit = sum(nt(minlt:lt))

           lop = (lt - minlt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           do k = 1,nc+ne
              sol%Q(lot:hit,k) = L(sol%t(lot:hit), tee(lt), qp(lop:hip,k), sol%INVLT)
              if (sol%deriv) then
                 sol%dQ(lot:hit,k) = L(sol%t(lot:hit), tee(lt), &
                      & qp(lop:hip,k)*stmp(lop:hip), sol%INVLT)*sol%t(lot:hit)
              end if
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     
     !$OMP PARALLEL DO PRIVATE(calcZ,hp,vp,lot,hit,lop,hip) SHARED(sol)
     do j = 1,sol%nx
        !$ write (*,'(I0,1X)',advance="no") OMP_get_thread_num()
        write (*,'(A,ES13.5)') 'x: ',sol%x(j)
        do i = 1,sol%ny

           calcZ = cmplx(sol%x(j),sol%y(i),DP)

           !! compute f(p) for all values of p at this location
           hp(1:tnP) =    headCalc(calcZ,stmp,1,tnP,dom,c,e,bg)
           vp(1:tnP,1:2) = velCalc(calcZ,stmp,1,tnP,dom,c,e,bg)

           !! invert solutions one log-cycle of t at a time
           do lt = minlt,maxlt-1
              
              if (nt(lt) == 0) then
                 ! empty logcycle
                 cycle
              end if

              !! group of times in current log cycle
              lot = 1 + sum(nt(minlt:lt-1))
              hit = sum(nt(minlt:lt))

              !! group of Laplace parameters corresponding to this logcycle
              lop = (lt - minlt)*size(s,dim=1) + 1
              hip = lop + size(s,dim=1) - 1

              sol%h(lot:hit,j,i) =     L(sol%t(lot:hit), tee(lt), hp(lop:hip), sol%INVLT)
              sol%v(lot:hit,1:2,j,i) = L(sol%t(lot:hit), tee(lt), vp(lop:hip,1:2), sol%INVLT)

              if (sol%deriv) then
                 ! multiply solution by p for derivative
                 hp(lop:hip) = hp(lop:hip)*stmp(lop:hip)
                 sol%dh(lot:hit,j,i) = L(sol%t(lot:hit), tee(lt), hp(lop:hip), &
                      & sol%INVLT)*sol%t(lot:hit)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(sol%timeseries) then ! x,y locations are in pairs; e.g. inner product

     write(stdout,'(A)') 'compute solution for plotting time series'
     allocate(sol%h(sol%nt,sol%nx,1), hp(tnP), sol%v(sol%nt,2,sol%nx,1), &
          & vp(tnP,2), stmp(tnP))
     if (sol%deriv) then
        allocate(sol%dh(sol%nt,sol%nx,1))
     end if

     stmp(1:tnp) = reshape(s,[tnP])

     ! TODO: see repeated code above
     if (sol%Qcalc) then
        write(stdout,'(A)',advance='no') 'computing element boundary flowrates: '
        allocate(sol%Q(sol%nt,nc+ne),qp(tnP,nc+ne))
        if (sol%deriv) then
           allocate(sol%dQ(sol%nt,nc+ne))
        end if
        do j = 1,nc
           write(stdout,'(A,I0,1X)',advance='no') 'C',j
           qp(:,j) = elementFlowrate(c(j)%matching,stmp,1,tnP,dom,c,e,bg)
        end do
        do j = 1,ne
           write(stdout,'(A,I0,1X)',advance='no') 'E',j
           qp(:,nc+j) = elementFlowrate(e(j)%matching,stmp,1,tnP,dom,c,e,bg)
        end do
        write(stdout,'(/)')

        !$OMP PARALLEL DO PRIVATE(lt,lot,hit,lop,hip,k) shared(sol,qp)
        do lt = minlt,maxlt-1
           lot = 1 + sum(nt(minlt:lt-1))
           hit = sum(nt(minlt:lt))

           lop = (lt - minlt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           do k = 1, nc+ne
              sol%Q(lot:hit,k) = L(sol%t(lot:hit), tee(lt), qp(lop:hip,k), sol%INVLT)
              if (sol%deriv) then
                 sol%dQ(lot:hit,k) = L(sol%t(lot:hit), tee(lt), &
                      & qp(lop:hip,k)*stmp(lop:hip), sol%INVLT)*sol%t(lot:hit)
              end if
           end do
        end do
        !$OMP END PARALLEL DO
     end if

     !$OMP PARALLEL DO PRIVATE(calcZ,hp,vp,lot,hit,lop,hip) SHARED(sol)
     do i = 1,sol%nx
        write(stdout,'(A,2(3X,ES14.7E1))') trim(sol%obsname(i)),sol%xshift+sol%x(i),&
             & sol%yshift+sol%y(i)

        calcZ = cmplx(sol%x(i),sol%y(i),DP)
        hp(1:tnP) =    headCalc(calcZ,stmp,1,tnP,dom,c,e,bg)
        vp(1:tnP,1:2) = velCalc(calcZ,stmp,1,tnP,dom,c,e,bg)

        do lt = minlt,maxlt-1
           lot = 1 + sum(nt(minlt:lt-1))
           hit = sum(nt(minlt:lt))

           lop = (lt - minlt)*size(s,1) + 1
           hip = lop + size(s,1) - 1

           ! don't need second dimension of results matrices
           sol%h(lot:hit,i,1) =     L(sol%t(lot:hit),tee(lt),hp(lop:hip),sol%INVLT)
           sol%v(lot:hit,1:2,i,1) = L(sol%t(lot:hit),tee(lt),vp(lop:hip,1:2),sol%INVLT)

           if (sol%deriv) then
              ! multiply solution by p for derivative
              hp(lop:hip) = hp(lop:hip)*stmp(lop:hip)
              sol%dh(lot:hit,i,1) = L(sol%t(lot:hit),tee(lt),hp(lop:hip)&
                   &,sol%INVLT)*sol%t(lot:hit)
           end if
        end do
     end do
     !$OMP END PARALLEL DO

  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! write output to file
  call writeResults(sol,part)

end program ltaem_main

