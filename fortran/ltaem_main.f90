! this is the driver routine for LT-AEM programs
! this can compute a time-domain solution either on a grid of locations and times,
! at a single point through time (i.e., hydrograph), or at a particle location through
! space and time.  The solution depends on the given properties of circular and 
! elliptical elements.  Wells or points are treated as special case circles, and lines are
! treated as special case ellipses.

program ltaem_main
  use constants, only : DP, PI
  use type_definitions, only : domain, element, circle, ellipse, solution, INVLT, particle
  use file_ops, only : readinput, writeresults
  use inverse_Laplace_Transform, only : invlap => deHoog_invlap, pvalues => deHoog_pvalues
  use solution, only : matrix_solution
  use calc_routines
  use element_geometry, only : distanceAngleCalcs
  use ellipse_mathieu_init, only : ellipse_init

  implicit none

  ! structs that organize variables
  type(domain) :: dom
  type(element) :: bg
  type(circle), allocatable :: c(:)
  type(ellipse), allocatable :: e(:)

  type(solution) :: sol
  type(INVLT) :: lap
  type(particle), allocatable :: p(:)
  
  integer :: i, j, part, tnp
  integer, dimension(4) :: ic
  integer, allocatable :: parnumdt(:)

  real(DP), allocatable :: logt(:), tee(:)
  integer, allocatable :: nt(:), run(:)
  complex(DP), allocatable :: s(:,:)
  integer :: ilogt, iminlogt, imaxlogt, lot, hit, lop, hip, ierr, lo
  character(20) :: tmpfname

  real(DP), parameter :: MOST = 0.99_DP  ! 'most' of a log-cycle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! either specify here or ask for at prompt
  sol%infname = 'input'

  ! read in data, initialize variables, allocate major structs
  call readInput(sol,lap,dom,bg,c,e,p)

  ! nudge times on 'edge' of logcycle down a tiny bit to increase accuracy
  where (abs(nint(sol%t(:)) - sol%t(:)) < epsilon(sol%t(:)))
     sol%t(:) = sol%t(:) - epsilon(sol%t(:))
  end where

  allocate(run(1:2*lap%m+1))
  forall(i=0:2*lap%m) run(i+1)=i

  ! independent of most choices
  call DistanceAngleCalcs()

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (sol%calc) then
     if (sol%particle) then   ! particle tracking
        
        ! need to do something better here about particles starting at t=0?
        do i=1, sol%nPart
           if (p(i)%ti < 1.0D-5) then
              p(i)%ti = 1.0D-5
              write(*,'(A,I0,A)') '^^^^ start time for particle ',i,' reset to 1.0E-5 ^^^^'
           end if
        end do
        
        iminlogt = floor(  minval(log10(p(:)%ti)))
        imaxlogt = ceiling(maxval(log10(p(:)%tf)))

        allocate(s(2*INVm+1,iminlogt:imaxlogt), nt(iminlogt:imaxlogt), &
             & tee(iminlogt:imaxlogt))

        nt(iminlogt:imaxlogt) = 1

        do ilogt = iminlogt, imaxlogt
           tee(ilogt) = min(10.0_DP**(ilogt + MOST), maxval(p(:)%tf))*2.0
           s(:,ilogt) = pvalues(tee(ilogt),lap)
        end do

        ! to make it possible for particles / contours to share code...
        imaxlogt = imaxlogt + 1

        deallocate(run)
        
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     else   ! contour maps / hydrographs
        allocate(logt(1:sol%numt))

        lo = 1
        logt(:) = log10(sol%t(:))
        iminlogt =   floor(minval(logt(:)))
        imaxlogt = ceiling(maxval(logt(:)))

        allocate(s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1))

        do ilogt = iminlogt, imaxlogt-1
           nt(ilogt) = count(logt >= real(ilogt,DP) .and. logt < real(ilogt+1,DP))
           tee(ilogt) = maxval(sol%t(lo:lo+nt(ilogt)-1))*2.0
           s(:,ilogt) = pvalues(t(ilogt),lap)
           lo = lo + nt(ilogt)
        end do
        deallocate(logt,run)
     end if

     ! only do matching if there is at least one matching element
     if(count(e(:)%match) + count(c(:)%match) > 0) then


!!$        allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
!!$             & Gm(1:4*CIn+2, 1:2*CIm, 1:CInum))
        
        coeff = CZERO
     
        ! calculate coefficients for each value of Laplace parameter
        ! ** common between particle tracking and contours/hydrographs **

        do ilogt = iminlogt,imaxlogt-1
           write(*,'(A,I3)') 'log t=',ilogt

           do j = 1,2*lap%m+1
              if (nt(ilogt) > 0) then
                 write(*,'(I4,1X,2(A,ES10.3),A)') j, '(',real(s(j,ilogt)),',',aimag(s(j,ilogt)),')'
                 ! compute generalized inverse of A matrix
                 call circInverse_old(s(j,ilogt),Gm)  
!!$               ! use direct approach to solve once
!!$              call circInverse_matrix(s(j,ilogt),coeff(j,ilogt,0:4*CIn+1,1:CInum))  

                 if (j > 1) then  ! copy previous as initial guess
                    coeff(j,ilogt,0:4*CIn+1,1:CInum) = coeff(j-1,ilogt,0:4*CIn+1,1:CInum)
                 end if
                 call matchIter(Gm, coeff(j,ilogt,0:4*CIn+1,1:CInum), s(j,ilogt),CImatchTol)
              end if
           end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! save coefficient matrix to file
        open(UNIT=77, FILE=sol%coeffFName, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
        if (ierr /= 0) then
           print *, 'error writing intermediate coefficient matrix to file',sol%coefffname
           stop 000
        end if
        write(77,*) iminlogt,imaxlogt,nt,s,tee
        write(77,*) coeff
        close(77)
        write(*,'(A)') '<matching finished>'
     end if

  else ! do not re-calculate coefficients (this only makes sense if num matching elements > 0)

     if(count(e(:)%match) + count(c(:)%match) > 0) then
        open(UNIT=77, FILE=BGCoeffFName, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
        if (ierr /= 0) then
           ! go back and recalc if no restart file
           write(*,'(A)') 'error opening restart file, recalculating...'
           sol%calc = .true.
           goto 111
        end if

        read(77,*) iminlogt,imaxlogt,nt,s,tee
        allocate(coeff(2*lap%m+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
             & s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1))

        ! read needed data calculated in previous run from file
        rewind(77)
        read(77)  iminlogt,imaxlogt,coeff,nt,s,tee,WLstorCoeff
        close(77)
        write(*,'(A)') 'matching results read from file'
     end if

  end if ! re-calculate coefficients

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(sol%contour) then ! contour output (x,y locations outerproduct of x,y vectors)

     write(*,'(A)') 'compute solution for plotting contours'
     allocate(sol%h(sol%nx,sol%ny,sol%nt), sol%hp(size(s)), &
          &   sol%vx(sol%nx,sol%ny,sol%nt), sol%vxp(size(s)), &
          &   sol%vy(sol%nx,sol%ny,sol%nt), sol%vyp(size(s)))

     ic = shape(coeff)  !! integer vector of matrix size in each dimension
     tnp = ic(1)*ic(2)  !! total number of Laplace parameters, over all times

     do j = 1,sol%nx
        write (*,'(A,ES14.6E1)') 'x: ',sol%x(j)
        do i = 1,sol%ny

           !! compute f(p) for all values of p -- not just one log cycle of time -- at this location 
           sol%hp(1:tnp) = headCalc(reshape(s,[tnp]),sol%x(j),sol%y(i),reshape(coeff,[tnp,ic(3:4)]))
           sol%vxp(1:tnp) = velxCalc(reshape(s,[tnp]),sol%x(j),sol%y(i),reshape(coeff,[tnp,ic(3:4)]))
           ! velx and vely share most of the calculations (velx _must_ be called first)
           sol%vyp(1:tnp) = velyCalc(reshape(s,[tnp]),sol%x(j),sol%y(i))

           !! invert solutions one log-cycle of t at a time
           do ilogt = iminlogt,imaxlogt-1
              
              !! group of times in current log cycle
              lot = 1 + sum(nt(iminlogt:ilogt-1))
              hit = sum(nt(iminlogt:ilogt))
              
              !! group of Laplace parameters corresponding to this logcycle
              lop = (ilogt - iminlogt)*ic(1) + 1
              hip = lop + ic(1) - 1

              sol%h(j,i,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%hp(lop:hip),lap)
              sol%vx(j,i,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%vxp(lop:hip),lap)
              sol%vy(j,i,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%vyp(lop:hip),lap)
           end do
        end do
     end do
     
  else ! hydrograph output (x,y locations are in pairs; e.g. inner product)

     write(*,'(A)') 'compute solution for plotting hydrograph'

     allocate(sol%h(sol%nx,1,sol%nt), sol%hp(size(s)), &
          &   sol%vx(sol%nx,1,sol%nt), sol%vxp(size(s)), &
          &   sol%vy(sol%nx,1,sol%nt), sol%vyp(size(s)))
     
     do i = 1,sol%nx
        write(*,'(A,2(1X,ES14.7E1))') 'location:',sol%x(i),sol%y(i)

        if (allocated(coeff)) then
           ic = shape(coeff)
           tnp = ic(1)*ic(2)   
        else
           tnp = size(s)
           ic(1:4) = [0,0,0,0]
        end if

        sol%hp(1:tnp) = headCalc(reshape(s,[tnp]),sol%x(i),sol%y(i),reshape(coeff,[tnp,ic(3:4)]))
        sol%vxp(1:tnp) = velxCalc(reshape(s,[tnp]),sol%x(i),sol%y(i),reshape(coeff,[tnp,ic(3:4)]))
        sol%vyp(1:tnp) = velyCalc(reshape(s,[tnp]),sol%x(i),sol%y(i))
        
        do ilogt = iminlogt,imaxlogt-1
           lot = 1 + sum(nt(iminlogt:ilogt-1))
           hit = sum(nt(iminlogt:ilogt))

           lop = (ilogt - iminlogt)*ic(1) + 1
           hip = lop + ic(1) - 1

           ! don't need second dimension of results matricies
           sol%h(i,1,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%hp(lop:hip),lap)
           sol%vx(i,1,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%vxp(lop:hip),lap)
           sol%vy(i,1,lot:hit) = invlap(sol%t(lot:hit),tee(ilogt),sol%vyp(lop:hip),lap)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file

  if (sol%particle) then 
     ! make sure correct output if doing particle tracking
     if (sol%output /= 4 .and. sol%output /= 5) then
        sol%output = 4
     end if
  else
     ! reset back to gnuplot output if not particle tracking
     if (sol%output == 4) sol%output = 1
  end if

  ! call subroutine to write output to file
  call writeResults(head,velx,vely,BGx,BGy,BGt,BGoutput,BGoutfname,Presult)

end program ltaem_main

