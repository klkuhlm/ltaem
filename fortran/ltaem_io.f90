! this module contains most of the I/O associated with LT-AEM

module file_ops
  implicit none  
  
  private
  public :: readInput, writeResults, writeGeometry

contains

  !##################################################
  ! this routine read the main input file, and allocates the main 
  ! data structures used to store data.
  subroutine readInput(sol,dom,bg,c,e,p)
    use type_definitions, only : solution, particle, domain, element, circle, ellipse
    use constants, only : DP, lenFN, PI

    type(solution), intent(inout) :: sol
    type(particle), intent(out), allocatable :: p(:)
    type(domain),   intent(out) :: dom
    type(element),  intent(out) :: bg
    type(circle),   intent(out), allocatable :: c(:)
    type(ellipse),  intent(out), allocatable :: e(:)

    character(4) :: chint
    character(20), dimension(3) :: fmt
    character(46) :: lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
    character(lenFN+5) :: echofname, circleFname, ellipseFname, particleFname
    integer :: ierr,j,ntot,nC,nE  ! #elements, #circles, #ellipses

    open(unit=15, file=trim(sol%infname), status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening input file ',trim(sol%infname)
       stop 100
    endif

    echofname = trim(sol%infname)//'.echo'
    open(unit=16, file=echofname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening echo file ',echofname
       stop 101
    else
       ! add a file variable to set Emacs to auto-revert mode
       write(16,'(A)') '-*-auto-revert-*-'
    endif   


    ! solution-specific and background aquifer parameters
    read(15,*,iostat=ierr) sol%calc, sol%particle, sol%contour, sol%output, &
         & sol%outFname, sol%coeffFName, sol%elemHfName, sol%geomfName
    if (ierr /= 0) stop 'error on line 1 of input file'
    ! types:: 3*logical, integer, 4*string
    ! some simple debugging of problem-type / output-type combinations
    if (sol%output < 1 .or. (sol%output > 5 .and. sol%output /= 10 &
         & .and. sol%output /= 11)) then
       print *, 'input file (line 1) sol%output must be in {1,2,3,4,5,10,11} ',sol%output
       stop
    end if
    if ((sol%output < 4 .or. sol%output > 5) .and. sol%particle) then
       print *, 'input file (line 1) if sol%particle==.True., sol%output should &
            &be in {4,5}',sol%output,sol%particle
       stop
    elseif (.not. sol%particle .and. (sol%output == 4 .or. sol%output == 5)) then
       print *, 'input file (line 1) sol%output should only be in {4,5} if &
            &sol%particle==.True.',sol%output,sol%particle
       stop
    end if
    if ((sol%output == 3 .or. sol%output == 11) .and. sol%contour) then
       print *, 'input file (line 1) sol%output should not be in {3,11} when contour &
            &output is selected',sol%output,sol%contour
       stop
    elseif ((sol%output /= 3 .and. sol%output /= 11) .and. .not. sol%contour) then
       print *, 'input file (line 1) sol%output should be in {3,11} when hydrograph &
            &(non-contour) output is selected',sol%output,sol%contour
       stop
    end if


    read(15,*,iostat=ierr) bg%por, bg%k, bg%ss, bg%leakFlag, &
         & bg%aquitardK, bg%aquitardSs, bg%aquitardb, bg%ms
    if (ierr /= 0) stop 'error on line 2 of input file'
    ! reals checked here, bg%ms checked in ellipse section, leakflag checked elsewhere
    if (any([bg%por,bg%k,bg%ss] < epsilon(0.0))) then
       print *, 'input file (line 2) bg%por, bg%k, bg%ss &
            &must all be > 0.0 :',bg%por,bg%k,bg%ss
       stop
    end if
    if (any([bg%aquitardK,bg%aquitardSs,bg%aquitardb] < epsilon(0.0))) then
       print *, 'input file (line 2) bg%aquitardK, bg%aquitardSs, bg%aquitardb &
            &must all be > 0.0 :',bg%aquitardK,bg%aquitardSs,bg%aquitardb 
       stop
    end if


    read(15,*,iostat=ierr) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b
    if (ierr /= 0) stop 'error on line 3 of input file'
    if (any([bg%Sy, bg%kz,bg%b] < epsilon(0.0))) then
       print *, 'input file (line 3) bg%Sy, bg%kz, bg%b must all be > 0.0 :', &
            & bg%Sy,bg%kz,bg%b
       stop
    end if

    ! echo input from first 3 lines to file
    write(16,'(3(1L,1X),I0,5(1X,A))') sol%calc, sol%particle, sol%contour, sol%output, &
         & trim(sol%outFname), trim(sol%coeffFName), trim(sol%elemHfName), &
         & trim(sol%geomFname),'  ||    re-calculate coefficients?, particle?, &
         & contour?, output, out/coeff/hierarchy/geometry file names'
    write(16,'(3(ES11.5,1X),1L,3(1X,ES11.5),1X,I0,A)') bg%por, bg%k, bg%ss, &
         & bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb, bg%ms, & 
         & '  ||    por, k, Ss, leaky flag, K2, Ss2, b2, ellipse MS'
    write(16,'(2(ES11.5,1X),1L,1X,ES11.5,A)') bg%Sy, bg%kz, bg%unconfinedFlag, &
         & bg%b, '  || Sy, Kz, unconfined?, BGb'
    

    ! desired solution points/times
    read(15,*,iostat=ierr) sol%nx, sol%ny, sol%nt
    if (ierr /= 0) stop 'error on line 4 of input file'
    if (any([sol%nx,sol%ny,sol%nt] < 1)) then
       print *, 'input file (line 4) sol%nx, sol%ny, sol%nt must be > 0:',&
             & sol%nx, sol%ny, sol%nt
       stop
    end if
    allocate(sol%x(sol%nx), sol%y(sol%ny), sol%t(sol%nt), stat=ierr)
    if (ierr /= 0) stop 'ltaem_io.f90 error allocating: sol%x,sol%y,sol%t'
    read(15,*,iostat=ierr) sol%x(:)
    if (ierr /= 0) stop 'error on line 5 (sol%x(:)) of input file'
    read(15,*,iostat=ierr) sol%y(:)
    if (ierr /= 0) stop 'error on line 6 (sol%y(:)) of input file'
    do j=1,sol%nt  !! modified to accommidate pest (make time a column)
       read(15,*,iostat=ierr) sol%t(j)
       if (ierr /= 0)  stop 'error reading time from input (> line 6)'
    end do
    if (.not. sol%particle) then
       fmt(1) = '(    (ES12.5,1X),A) '
       write(16,'(3(I0,1X),A)') sol%nx, sol%ny, sol%nt, '  ||    numX, numY, numt'
       write(fmt(1)(2:5),'(I4.4)') sol%nx
       write(16,fmt(1)) sol%x(:), '  ||    x Vector'
       write(fmt(1)(2:5),'(I4.4)') sol%ny
       write(16,fmt(1)) sol%y(:), '  ||    y Vector'
       write(fmt(1)(2:5),'(I4.4)') sol%nt
       write(16,fmt(1)) sol%t(:), '  ||    t Vector'
    endif

    ! inverse Laplace transform parameters
    read(15,*,iostat=ierr) sol%alpha, sol%tol, sol%m
    if (ierr /= 0) stop 'error on line 7 of input file'
    if (sol%M < 1) then
       print *, 'input file (line 7) sol%M must be > 0: ',sol%M
       stop
    end if
    if (sol%tol < epsilon(sol%tol)) then
       sol%tol = epsilon(sol%tol)
       write(*,'(A,ES11.5)') 'WARNING: increased INVLAP solution tolerance to ',sol%tol 
    end if
    write(16,'(2(ES11.5,1X),I0,A)') sol%alpha, sol%tol, sol%m,'  ||    alpha, tol, M'

    ! circular (includes wells)
    read(15,*) dom%num(1),circleFname
    nc = dom%num(1)
    if (nc > 0) then
       open(unit=22, file=trim(circleFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening circular data file for reading ',trim(circleFname)
          stop
       else
          write(16,'(A)') trim(circleFname)//' opened for reading circular data'
       end if

       allocate(c(nc), stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating c()'
       read(22,*) c(:)%n
       if (any(c%n < 1)) then
          print *, 'c%N must not be < 1',c%N
          stop
       end if

       read(22,*) c(:)%m
       if (any(c%M < 1)) then
          print *, 'c%M must not be < 1',c%M
          stop
       end if

       read(22,*) c(:)%ibnd
       if (any(c%ibnd < -1 .or. c%ibnd > 2)) then
          print *, 'c%ibnd must be in {-1,0,1,2}',c%ibnd
          stop
       end if
          
       read(22,*) c(:)%CalcIn ! calculate inside this element?
       read(22,*) c(:)%StorIn ! account for free-water storage effects inside element?

       read(22,*) c(:)%r
       if (any(c%r < epsilon(0.0))) then
          print *, 'c%r must be > 0.0',c%r
          stop
       end if

       read(22,*) c(:)%x ! any real number is ok
       read(22,*) c(:)%y

       read(22,*) c(:)%k
       if (any(c%k < epsilon(0.0))) then
          print *, 'c%K must be > 0.0',c%k
          stop
       end if      

       read(22,*) c(:)%Ss
       if (any(c%ss < epsilon(0.0))) then
          print *, 'c%Ss must be > 0.0',c%Ss
          stop
       end if      

       read(22,*) c(:)%por
       if (any(c%por < epsilon(0.0))) then
          print *, 'c%por must be > 0.0',c%por
          stop
       end if      

       read(22,*) c(:)%areaQ ! area source strength (flux)
       read(22,*) c(:)%bdryQ ! streng of specified value on bdry (head or flux) 
       do j=1,size(c,dim=1)
          read(22,'(I3)', advance='no') c(j)%AreaTime 
          if (c(j)%AreaTime > -1) then
             ! functional time behavior
             allocate(c(j)%ATPar(2),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating c%(j)%ATPAR(2)'
             read(22,*) c(j)%ATPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%AreaTime,c(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for circle ',j
          else
             ! piecewise-constant time behavior 
             allocate(c(j)%ATPar(-2*c(j)%AreaTime+1),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating c%(j)%ATPAR(:)'
             read(22,*) c(j)%ATPar(:)
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
             write(lfmt(8:11),'(I4.4)')  size(c(j)%ATPar(:-c(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%ATPar(-c(j)%AreaTime+2:),1)
             write(16,lfmt) c(j)%AreaTime,c(j)%ATPar(:-c(j)%AreaTime+1),'| ',&
                  & c(j)%ATPar(-c(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for circle ',j
          end if
       end do
       do j=1,size(c,dim=1)
          read(22,'(I3)', advance='no') c(j)%BdryTime
          if (c(j)%BdryTime > -1) then
             allocate(c(j)%BTPar(2),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating c%(j)%BTPAR(2)'
             read(22,*) c(j)%BTPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%BdryTime,c(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for circle ',j
          else
             allocate(c(j)%BTPar(-2*c(j)%BdryTime+1),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating c%(j)%BTPAR(:)'
             read(22,*) c(j)%BTPar(:)
             write(lfmt(8:11),'(I4.4)')  size(c(j)%BTPar(:-c(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%BTPar(-c(j)%BdryTime+2:),1)
             write(16,lfmt) c(j)%BdryTime,c(j)%BTPar(:-c(j)%BdryTime+1),' | ',&
                  & c(j)%BTPar(-c(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for circle ',j
          end if
       end do
       read(22,*) c(:)%leakFlag  ! checking handled elsewhere
       read(22,*) c(:)%aquitardK
       if (any(c%aquitardK < epsilon(0.0))) then
          print *, 'c%aquitardK must be > 0.0',c%aquitardk
          stop
       end if      

       read(22,*) c(:)%aquitardSs
       if (any(c%aquitardSS < epsilon(0.0))) then
          print *, 'c%aquitardSs must be > 0.0',c%aquitardSs
          stop
       end if      

       read(22,*) c(:)%aquitardb  !aquitard thickness
       if (any(c%aquitardB < epsilon(0.0))) then
          print *, 'c%aquitardB must be > 0.0',c%aquitardB
          stop
       end if      

       read(22,*) c(:)%unconfinedFlag ! checking handled elsewhere
       read(22,*) c(:)%Sy
       if (any(c%sy < epsilon(0.0))) then
          print *, 'c%Sy must be > 0.0',c%sy
          stop
       end if      

       read(22,*) c(:)%Kz
       if (any(c%kz < epsilon(0.0))) then
          print *, 'c%Kz must be > 0.0',c%kz
          stop
       end if      

       read(22,*) c(:)%b  ! aquifer thickness
       if (any(c%b < epsilon(0.0))) then
          print *, 'c%B must be > 0.0',c%b
          stop
       end if      

       read(22,*) c(:)%dskin ! dimensionless skin
       if (any(c%dskin < epsilon(0.0))) then
          print *, 'c%Dskin must be > 0.0',c%dskin
          stop
       end if            

       close(22)
   
       where (c(:)%ibnd == -1 .or. c(:)%ibnd == 0 .or. c(:)%ibnd == +1)
          c(:)%match = .true.
       elsewhere !( currently just 2)
          c(:)%match = .false.       
       end where
       
       if(any(c%match .and. c%storin)) then
          print *, 'WARNING: wellbore (free-water) storage only handled yet for ibnd==2'
          where(c%match .and. c%storin)
             c%storin = .false.
          end where
       end if

       write(chint,'(I4.4)') dom%num(1)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical 
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,'(I0,A)') dom%num(1), '  ||   number of circular elements (including wells)'
       write(16,fmt(1)) c(:)%n,'  ||   number of circular free parameter (Fourier coeffs)'
       write(16,fmt(1)) c(:)%m,'  ||   number of circular matching locations'
       write(16,fmt(1)) c(:)%ibnd, '  ||    circle ibnd array'
       write(16,fmt(2)) c(:)%match, '  ||    circle matching array'
       write(16,fmt(2)) c(:)%calcin, '  ||    calculate inside this circle?'
       write(16,fmt(2)) c(:)%storin, '  ||    calculate free-water storage effects of circle?'
       write(16,fmt(3)) c(:)%r, '  ||    circle radius'
       write(16,fmt(3)) c(:)%x, '  ||    circle center x'
       write(16,fmt(3)) c(:)%y, '  ||    circle center y'
       write(16,fmt(3)) c(:)%k, '  ||    circle aquifer k'
       write(16,fmt(3)) c(:)%ss, '  ||    circle aquifer Ss'
       write(16,fmt(3)) c(:)%por, '  ||    circle aquifer porosity'
       write(16,fmt(3)) c(:)%areaQ, '  ||    circle area rch rate'
       write(16,fmt(3)) c(:)%bdryQ, '  ||    circle boundry rch rate or head'
       write(16,fmt(2)) c(:)%leakFlag, '  ||     circle leaky type'
       write(16,fmt(3)) c(:)%aquitardK, '  ||     circle leaky aquitard K'
       write(16,fmt(3)) c(:)%aquitardSs, '  ||     circle leaky aquitard Ss'
       write(16,fmt(3)) c(:)%aquitardb, '  ||     circle leaky aquitard thickness'
       write(16,fmt(2)) c(:)%unconfinedFlag, '  ||     circle unconfined flag'
       write(16,fmt(3)) c(:)%Sy, '  ||     circle aquifer specific yield'
       write(16,fmt(3)) c(:)%Kz, '  ||     circle aquifer vertical K'
       write(16,fmt(3)) c(:)%b,'  ||    circle aquifer thickness'
       write(16,fmt(3)) c(:)%dskin,'  ||    circle boundary dimensionless skin factor'


       ! minor checking / correcting
       do j = 1,size(c,dim=1)
          if (c(j)%ibnd == 2 .and. c(j)%N /= 1) then
             write(*,'(A,I0,A)') '** wells (ibnd==2) must have N=1, fixing circle #',j,' to N=1'
             c(j)%N = 1
          end if
       end do
       
    else
       ! no circular elements
       allocate(c(0),stat=ierr)
       if (ierr /= 0) stop 'error allocating c(0)'
    end if

    ! elliptical (includes line sources/sinks)
    read(15,*) dom%num(2), ellipseFname
    ne = dom%num(2)
    if (ne > 0) then
       open(unit=33, file=trim(ellipseFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening elliptical data file for reading ',trim(ellipseFname)
          stop
       else
          write(16,'(A)') trim(ellipseFname)//' opened for reading elliptical data'
       end if

       allocate(e(ne),stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating e()'
       read(33,*) e(:)%n
       if (any(e%n < 1)) then
          print *, 'e%N must not be < 1',e%N
          stop
       end if

       read(33,*) e(:)%m
       if (any(e%m < 1)) then
          print *, 'e%M must not be < 1',e%M
          stop
       end if

       read(33,*) e(:)%ms ! bg%ms read above
       ! checked more carefully in Mathieu function library
       if (any(e%ms < e%n) .or. any(bg%ms < e%n)) then  
          print *, 'e%ms must not be < e%n + buffer: e%N',e%N,'e%MS',e%ms
          print *, 'bg%ms must not be < e%n + buffer: bg%MS',bg%ms
          stop
       end if

       read(33,*) e(:)%ibnd
       if (any(e%ibnd < -1 .or. e%ibnd > 2)) then
          print *, 'e%ibnd must be in {-1,0,1,2}',e%ibnd
          stop
       end if

       read(33,*) e(:)%CalcIn
       read(33,*) e(:)%StorIn
       
       read(33,*) e(:)%r   ! eta
       if (any(e%r < 0.0)) then
          print *, 'e%r must be >= 0.0',e%r
          stop
       end if

       read(33,*) e(:)%x
       read(33,*) e(:)%y

       read(33,*) e(:)%f
       if (any(e%f < epsilon(0.0))) then
          print *, 'e%f must be > 0.0',e%f
          stop
       end if

       read(33,*) e(:)%theta      
       if (any(e%theta < -PI) .or. any(e%theta > PI)) then
          print *, 'e%theta must be -pi <= theta <= PI',e%theta
          stop
       end if

       read(33,*) e(:)%k
       if (any(e%k < epsilon(0.0))) then
          print *, 'e%k must be > 0.0',e%k
          stop
       end if

       read(33,*) e(:)%Ss
       if (any(e%Ss < epsilon(0.0))) then
          print *, 'e%Ss must be > 0.0',e%Ss
          stop
       end if

       read(33,*) e(:)%por
       if (any(e%por < epsilon(0.0))) then
          print *, 'e%por must be > 0.0',e%por
          stop
       end if

       read(33,*) e(:)%areaQ
       read(33,*) e(:)%bdryQ

       do j=1,size(e,dim=1)
          read(33,'(I3)', advance='no') e(j)%AreaTime 
          if (e(j)%AreaTime > -1) then
             allocate(e(j)%ATPar(2),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating e%(j)%ATPAR(2)'
             read(33,*) e(j)%ATPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%AreaTime,e(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%ATPar(-2*e(j)%AreaTime+1),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating e%(j)%ATPAR(:)'
             read(33,*) e(j)%ATPar(:)
             write(lfmt(8:11),'(I3.3)')  size(e(j)%ATPar(:-e(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%ATPar(-e(j)%AreaTime+2:),1)
             write(16,lfmt) e(j)%AreaTime,e(j)%ATPar(:-e(j)%AreaTime+1),' | ',&
                  & e(j)%ATPar(-e(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for ellipse ',j
          end if
       end do
       do j=1,size(e,dim=1)
          read(33,'(I3)', advance='no') e(j)%BdryTime
          if (e(j)%BdryTime > -1) then
             allocate(e(j)%BTPar(2),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating e%(j)%BTPAR(2)'
             read(33,*) e(j)%BTPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%BdryTime,e(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%BTPar(-2*e(j)%BdryTime+1),stat=ierr)
             if (ierr /= 0) stop 'ltaem_io.f90 error allocating e%(j)%BTPAR(:)'
             read(33,*) e(j)%BTPar(:)
             write(lfmt(8:11),'(I3.3)')  size(e(j)%BTPar(:-e(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%BTPar(-e(j)%BdryTime+2:),1)
             write(16,lfmt) e(j)%BdryTime,e(j)%BTPar(:-e(j)%BdryTime+1),' | ',&
                  & e(j)%BTPar(-e(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for ellipse ',j
          end if
       end do

       read(33,*) e(:)%leakFlag  ! checking handled elsewhere

       read(33,*) e(:)%aquitardK
       if (any(e%aquitardK < epsilon(0.0))) then
          print *, 'e%aquitardK must be > 0.0',e%aquitardK
          stop
       end if

       read(33,*) e(:)%aquitardSs
       if (any(e%aquitardSs < epsilon(0.0))) then
          print *, 'e%aquitardSs must be > 0.0',e%aquitardSs
          stop
       end if

       read(33,*) e(:)%aquitardb  
       if (any(e%aquitardb < epsilon(0.0))) then
          print *, 'e%aquitardb must be > 0.0',e%aquitardb
          stop
       end if

       read(33,*) e(:)%unconfinedFlag

       read(33,*) e(:)%Sy
       if (any(e%Sy < epsilon(0.0))) then
          print *, 'e%Sy must be > 0.0',e%Sy
          stop
       end if

       read(33,*) e(:)%Kz
       if (any(e%Kz < epsilon(0.0))) then
          print *, 'e%Kz must be > 0.0',e%Kz
          stop
       end if

       read(33,*) e(:)%b
       if (any(e%b < epsilon(0.0))) then
          print *, 'e%b must be > 0.0',e%b
          stop
       end if

       read(33,*) e(:)%dskin
       if (any(e%dskin < epsilon(0.0))) then
          print *, 'e%dskin must be > 0.0',e%dskin
          stop
       end if
       
       close(33)
       
       where (e(:)%ibnd == -1 .or. e(:)%ibnd == 0 .or. e(:)%ibnd == +1)
          e(:)%match = .true.
       elsewhere
          e(:)%match = .false.       
       end where

       if(any(e%storin)) then
          print *, 'WARNING: wellbore (free-water) storage not handled for ellipses yet'
          e%storin = .false.
       end if

       write(chint,'(I4.4)') dom%num(2)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical 
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,'(I0,A)') dom%num(2), '  ||   number of elliptical elements (including lines)'
       write(16,fmt(1)) e(:)%n,'  ||   number of elliptical free parameter (Fourier coeffs)'
       write(16,fmt(1)) e(:)%m,'  ||   number of ellipse matching locations'
       write(16,fmt(1)) e(:)%ms,'  ||   size of "infinite" Mathieu matrices'
       write(16,fmt(1)) e(:)%ibnd, '  ||    ellipse ibnd array'
       write(16,fmt(2)) e(:)%match, '  ||    ellipse matching array'
       write(16,fmt(2)) e(:)%calcin, '  ||    calculate inside this ellipse?'
       ! free-water storage for ellipses not handled yet
       write(16,fmt(2)) e(:)%storin, '  ||    calculate free-water storage for this ellipse?' 
       write(16,fmt(3)) e(:)%r, '  ||    ellipse radius (eta)'
       write(16,fmt(3)) e(:)%x, '  ||    ellipse center x'
       write(16,fmt(3)) e(:)%y, '  ||    ellipse center y'
       write(16,fmt(3)) e(:)%f, '  ||    ellipse semi-focal length'
       write(16,fmt(3)) e(:)%theta, '  ||    ellipse angle rotation with +x axis'
       write(16,fmt(3)) e(:)%k, '  ||    ellipse aquifer k'
       write(16,fmt(3)) e(:)%ss, '  ||    ellipse aquifer Ss'
       write(16,fmt(3)) e(:)%por, '  ||    ellipse aquifer porosity'
       write(16,fmt(3)) e(:)%areaQ, '  ||     ellipse area rch rate'
       write(16,fmt(3)) e(:)%bdryQ, '  ||     ellipse boundary rch rate or head'
       write(16,fmt(2)) e(:)%leakFlag, '  ||     ellipse leaky type'
       write(16,fmt(3)) e(:)%aquitardK, '  ||     ellipse leaky aquitard K'
       write(16,fmt(3)) e(:)%aquitardSs, '  ||     ellipse leaky aquitard Ss'
       write(16,fmt(3)) e(:)%aquitardb, '  ||     ellipse leaky aquitard thickness'
       write(16,fmt(2)) e(:)%unconfinedFlag, '  ||     ellipse unconfined flag'
       write(16,fmt(3)) e(:)%Sy, '  ||     ellipse aquifer specific yield'
       write(16,fmt(3)) e(:)%Kz, '  ||     ellipse aquifer vertical K'
       write(16,fmt(3)) e(:)%b,'  ||    ellipse aquifer thickness'
       write(16,fmt(3)) e(:)%dskin,'  ||    ellipse boundary dimensionless skin factor'
    else
       allocate(e(0),stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating e(0)'
    end if
    
    ntot = sum(dom%num) ! total number of circular and elliptical elements
    if (ntot < 1) stop 'ltaem_io.f90 Need at least one circular (including well) or&
            & elliptical (including line) element.'

    ! compute secondary parameters
    bg%alpha = bg%K/bg%Ss
    c(:)%alpha = c(:)%K/c(:)%Ss
    e(:)%alpha = e(:)%K/e(:)%Ss
    bg%T = bg%K*bg%b
    c(:)%T = c(:)%K*c(:)%b
    e(:)%T = e(:)%K*e(:)%b

    write(chint,'(I4.4)') dom%num(1)
    fmt(2) = '('//chint//'(ES11.5,1X),A) ' ! circles
    write(chint,'(I4.4)') dom%num(2)
    fmt(3) = '('//chint//'(ES11.5,1X),A) ' ! ellipses

    if (dom%num(1) > 0) then
       write(16,fmt(2)) c(:)%alpha,'  ||    circle hydraulic diffusivity'
       write(16,fmt(2)) c(:)%T,'  ||    circle transmissivity'
    end if
    if (dom%num(2) > 0) then
       write(16,fmt(3)) e(:)%alpha,'  ||    ellipse hydraulic diffusivity'
       write(16,fmt(3)) e(:)%T,'  ||    ellipse transmissivity'
    end if
    write(16,'(ES11.5,A)') bg%alpha,'  ||    background hydraulic diffusivity'
    write(16,'(ES11.5,A)') bg%T,'  ||    background hydraulic diffusivity'

    ! particles
    if (sol%particle) then
       read(15,*) particleFname
       open(unit=44, file=particleFname, status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening particle data file for reading ',particleFname
          stop
       else
          write(16,'(A)') trim(particleFname)//' opened for reading particle data'
       end if

       read(44,*) sol%nPart,  sol%streakSkip
       allocate(p(sol%nPart),stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating p()'
       ! TODO: error checking for particle data
       read(44,*) p(:)%tol 
       read(44,*) p(:)%dt 
       read(44,*) p(:)%maxStep
       read(44,*) p(:)%x
       read(44,*) p(:)%y
       read(44,*) p(:)%ti
       read(44,*) p(:)%tf
       read(44,*) p(:)%int
       read(44,*) p(:)%InclIn
       close(44)

       write(chint,'(I4.4)') sol%nPart
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,fmt(3)) p(:)%tol,'  ||    particle solution tolerances'
       write(16,fmt(3)) p(:)%dt,'  ||    particle dt'
       write(16,fmt(3)) p(:)%maxStep,'  ||   particle max flux'
       write(16,fmt(3)) p(:)%x, '  ||    particle initial x'
       write(16,fmt(3)) p(:)%y, '  ||    particle initial y'
       write(16,fmt(3)) p(:)%ti, '  ||    particle initial t'
       write(16,fmt(3)) p(:)%tf, '  ||    particle maximum t'
       write(16,fmt(1)) p(:)%int, '  ||    particle integration method'
       write(16,fmt(2)) p(:)%InclIn, '  ||    particle begins inside CH/CF incl?'
    else
       ! no particles
       allocate(p(0),stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating: p(0)'
    endif
    close(15)
    close(16)
  end subroutine readInput

  !******************************************************
  subroutine writeResults(s,p)
    use constants, only : DP
    use type_definitions, only : solution, particle

    type(solution), intent(in) :: s
    type(particle), dimension(:), intent(in) :: p

    character(4), dimension(2) :: chint
    ! adjust the formats of time, location, and results here
    character(6) :: tfmt = 'ES13.5', xfmt = 'ES12.4'
    character(9) :: hfmt = 'ES22.14e3'
    integer :: i, j, k, nt

    select case (s%output)
    case (1) 
       ! ** gnuplot contour map friendly output **
       ! print results as x,y,z triplets with the given times separated by double blank lines

       open(unit=20, file=s%outfname, status='replace', action='write')
       write(20,*) '# ltaem contour map output'
       write(20,'(A,I0)') ' # t: ', s%nt
       write(20,'(A,I0)') ' # x: ', s%nx
       write(20,'(A,I0)') ' # y: ', s%ny
       write(20,'(A,I0)') ' # xy:', s%nx*s%ny     

       do i = 1, s%nt
          write(20,'(A,'//tfmt//')') ' # t= ',s%t(i)
          write(20,'(A)')   &
          & '#      X           Y               head&
          &                 velx                  vely'
          do j = 1, s%ny
             do k = 1, s%nx
                write(20,'(2('//xfmt//',1X),3('//hfmt//',1X))') &
                     & s%x(k), s%y(j), s%h(k,j,i), s%v(k,j,i,1:2)
             end do
          end do
          write(20,'(/)')
       end do       
       write(20,'(A)') '# EOF'
       close(20)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') ' gnuplot contour map style output written to ', trim(s%outfname)
       write(*,'(A)')  '***********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (2)
       
       ! ** matlab-friendly output **
       ! print results as matricies, with each variable (at each time) going to a separate file
       ! (x and y matricies - similar to results from matlab function meshgrid)
       
       ! x-matrix has same row repeated numy times
       open(unit=20, file=trim(s%outfname)//'_x.dat', status='replace', &
            & action='write')
       write(chint(1),'(i4.4)') s%nx
       do i = 1, s%ny
          write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%x(j), j=1,s%nx)
       end do
       close(20)

       ! y-matrix has same column repeated numx times
       open(unit=20, file=trim(s%outfname)//'_y.dat', status='replace', &
            & action='write')
       do i = 1, s%ny
          write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%y(i), j=1,s%nx)
       end do
       close(20)
       
       do k = 1, s%nt
          write(chint(2),'(i4.4)') k

          ! head-matrix
          open(unit=20, file=trim(s%outfname)//'_head_'//chint(2)//'.dat', &
               & status='replace', action='write')
          do i = 1, s%ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%h(j,i,k), j=1,s%nx)
          end do
          close(20)

          ! velx-matrix
          open(unit=20, file=trim(s%outfname)//'_velx_'//chint(2)//'.dat', &
               & status='replace', action='write')
          do i = 1, s%ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(j,i,k,1), j=1,s%nx)
          end do
          close(20)

          ! vely-matrix
          open(unit=20, file=trim(s%outfname)//'_vely_'//chint(2)//'.dat', &
               & status='replace', action='write')
          do i = 1, s%ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(j,i,k,2), j=1,s%nx)
          end do
          close(20)
       end do

       ! column of calculation times
       open(unit=20, file=trim(s%outfname)//'_t.dat', status='replace', &
            & action='write')
       write (20,'('//tfmt//')') (s%t(j), j=1,s%nt)
       close(20)

       write(*,'(/A)') '*********************************************************************'
       write(*,'(3A)') 'matlab output written to ', trim(s%outfname), &
              & '{x,y,t,head{1-n},velx{1-n},vely{1-n}}.dat'
       write(*,'(A)') '*********************************************************************'
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (3)

       ! ** gnuplot hydrograph-friendly output **
       ! column of time values at a location through time
       ! locations separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# ltaem hydrograph output'      
       do i = 1, s%nx
          write (20,'(2(A,'//xfmt//'))') ' # location: x=',s%x(i),' y=',s%y(i)
          write (20,'(A)')   '#     time              head&
               &                  velx                vely'
          do k = 1, s%nt
             write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                  & s%t(k),s%h(i,1,k),s%v(i,1,k,1:2)
          end do
          write (20,'(/)')
       end do       
       write(20,*) '# EOF'
       close(20)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') 'gnuplot style output written to ', trim(s%outfname)
       write(*,'(A)') '***********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (11)

       ! ** gnuplot hydrograph-friendly output **
       ! column of time values at a location through time
       ! locations separated by blank lines (grid no velocity)
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# ltaem hydrograph output'
       do i = 1, s%ny
          do j = 1, s%nx
             write (20,'(2(A,'//xfmt//'))') ' # location: x=',s%x(j),' y=',s%y(i)
             write (20,'(A)')   '#     time              head&
                  &                  velx                vely'
             do k = 1, s%nt
                write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                     & s%t(k),s%h(j,i,k),s%v(j,i,k,1:2)
             end do
             write (20,'(/)')
          end do
       end do       
       write(20,*) '# EOF'
       close(20)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') 'gnuplot style output written to ', trim(s%outfname)
       write(*,'(A)') '***********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (10)

       ! ** output for inverse in Matlab
       ! column of time values at a location through time
       ! locations separated by blank lines

       do i = 1, s%nx
          write(chint(1),'(I4.4)') i
          open(unit=20, file=trim(s%outfname)//'_'//chint(1), status='replace', &
               &action='write')
          do k = 1, s%nt
             write (20,'('//tfmt//',1X,'//hfmt//')') s%t(k),s%h(i,1,k)
          end do
          write(20,'(/)')
          close(20)
       end do       

       write(*,'(/A)') '***********************************************************'
       write(*,'(4A)') 'inverse output written to ',trim(s%outfname),'0000-',chint(1)
       write(*,'(A)') '***********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (4)

       ! ** pathline gnuplot-style output **
       ! columns of time values for starting locations
       ! particles separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# ltaem particle tracking output'
       do i = 1, size(p,dim=1)
          write (20,'(A,I0)') '# particle: ',i
          write (20,'(A)')   &
          & '#     time          x            y             velx          vely '
          do k=1,p(i)%numt
             write (20,'('//tfmt//',2'//xfmt//',2'//hfmt//')') &
                  & p(i)%result(1:5,k)
          end do
          write (20,'(/)')
       end do       
       write(20,'(A)') '# EOF'
       close(20)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') 'particle tracking output written to ', trim(s%outfname)
       write(*,'(A)') '***********************************************************'
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (5)

       ! ** streakline gnuplot-style output **
       ! this requires constant time steps (can't use adaptive integration)
       ! each block is a requested time, each row a particle

       open(unit=90, file=s%outfname, status='replace', action='write')
       write (90,'(A)') '# ltaem particle tracking streakfile output'

       ! max number of times for all particles
       nt = maxval(p(:)%numt,dim=1)

       do i = 1, nt, s%streakSkip
          ! use maxval to ensure a non-zero time is reported
          write (90,'(A'//tfmt//')') '# time:', maxval(p(i)%result(1,:)) 
          write (90,'(A)') '#  particle       x            y&
               &           velx         vely'
          do j = 1, size(p,dim=1)
             if (p(j)%result(1,i) > 0.0) then
                ! only write particle if it has non-zero data
                write (90,'(I0,4(1X'//hfmt//'))')  j,p(j)%result(2:5,i)
             end if
          end do
          write (90,'(/)')
       end do       
       write(90,'(A)') '# EOF'
       close(90)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') 'particle tracking streakfile written to ', trim(s%outfname)
       write(*,'(A)') '***********************************************************'

    case default
       write(*,'(A,I0)')  'invalid output code ',s%output
       stop 
    end select
  end subroutine writeResults  

  !##################################################
  subroutine writeGeometry(c,e,s)
    use constants, only : DP
    use type_definitions, only : circle, ellipse, solution
    
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: s
    integer :: nc, ne, i, j, ierr
    character(14) :: fmt

    nc = size(c,dim=1)
    ne = size(e,dim=1)
    fmt = '(2(ES13.5,1X))'

    ! write matching points to file
    open(unit=40, file=s%geomfname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       ! non-fatal error
       write(*,'(2A)') 'WARNING: ltaem_io.f90 error opening output file for writing &
            &element matching locations ',s%geomFname
    else    
       write(40,'(A)') '# points along circumference of circular and elliptical elements'
       do i = 1,nc
          write(40,'(A,I0)') '# circular element ',i
          write(40,fmt)  (c(i)%Zom(j), j=1,c(i)%M)
          if (c(i)%M > 1) then
             ! joins circle back up with beginning for plotting
             write(40,fmt)  c(i)%Zom(1)
          end if
          write(40,'(/)')    
       end do

       do i = 1,ne
          write(40,'(A,I0)') '# elliptical element ',i
          write(40,fmt)  (e(i)%Zom(j), j=1,e(i)%M)
          if (e(i)%M > 1) then
             ! joins ellipse back up with beginning for plotting
             write(40,fmt)  e(i)%Zom(1)
          end if
          write(40,'(/)')    
       end do
       write(40,'(A)') '# EOF'
       close(40)
    end if

  end subroutine writeGeometry
end module file_ops
