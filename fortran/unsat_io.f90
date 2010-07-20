! this module contains most of the I/O associated with unsaturated AEM

module unsat_file_ops
  implicit none  
  
  private
  public :: readInput, writeResults, writeGeometry

contains

  !##################################################
  ! this routine read the main input file, and allocates the main 
  ! data structures used to store data.
  subroutine readInput(sol,dom,bg,c,e)
    use constants, only : DP, lenFN, PI
    use unsat_type_definitions, only : solution, domain, element, circle, ellipse

    type(solution), intent(inout) :: sol
    type(domain),   intent(out) :: dom
    type(element),  intent(out) :: bg
    type(circle),   intent(out), allocatable :: c(:)
    type(ellipse),  intent(out), allocatable :: e(:)

    character(4) :: chint
    character(20), dimension(3) :: fmt
    character(lenFN+5) :: echofname, circleFname, ellipseFname
    integer :: ierr, ntot, nC, nE  ! #elements, #circles, #ellipses

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
    read(15,*,iostat=ierr) sol%output, sol%outFname, sol%geomfName
    if (ierr /= 0) stop 'error on line 1 of input file'

    read(15,*,iostat=ierr) bg%alpha, bg%qz0, bg%ms
    if (ierr /= 0) stop 'error on line 2 of input file'

    ! echo input from first 2 lines to file
    write(16,*) sol%output, trim(sol%outFname), &
         & trim(sol%geomFname),'  ||    output, out/geometry file names'

    write(16,*) bg%alpha, bg%qz0, bg%ms, '  ||    alpha, qz0, ellipse MS'

    ! desired solution points/times
    read(15,*,iostat=ierr) sol%nx, sol%ny
    if (ierr /= 0) stop 'error on line 3 of input file'
    if (any([sol%nx,sol%ny] < 1)) then
       print *, 'input file (line 3) sol%nx, sol%ny must be > 0:',&
             & sol%nx, sol%ny
       stop
    end if

    allocate(sol%x(sol%nx), sol%y(sol%ny), stat=ierr)
    if (ierr /= 0) stop 'ltaem_io.f90 error allocating: sol%x,sol%y'

    read(15,*,iostat=ierr) sol%x(:)
    if (ierr /= 0) stop 'error on line 4 (sol%x(:)) of input file'
    read(15,*,iostat=ierr) sol%y(:)
    if (ierr /= 0) stop 'error on line 4 (sol%y(:)) of input file'

    fmt(1) = '(    (ES12.5,1X),A) '
    write(16,'(2(I0,1X),A)') sol%nx, sol%ny, '  ||    numX, numY'
    write(fmt(1)(2:5),'(I4.4)') sol%nx
    write(16,fmt(1)) sol%x(:), '  ||    x Vector'
    write(fmt(1)(2:5),'(I4.4)') sol%ny
    write(16,fmt(1)) sol%y(:), '  ||    y Vector'

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

       read(22,*) c(:)%r
       if (any(c%r < epsilon(0.0))) then
          print *, 'c%r must be > 0.0',c%r
          stop
       end if

       read(22,*) c(:)%x ! any real number is ok
       read(22,*) c(:)%y

       read(22,*) c(:)%alpha
       if (any(c%alpha < epsilon(0.0))) then
          print *, 'c%alpha must be > 0.0',c%alpha
          stop
       end if      

       read(22,*) c(:)%bdryQ ! streng of specified value on bdry (head or flux) 
       close(22)
   
       where (c(:)%ibnd == -1 .or. c(:)%ibnd == 0 .or. c(:)%ibnd == +1)
          c(:)%match = .true.
       elsewhere !( currently just 2)
          c(:)%match = .false.       
       end where
       
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
       write(16,fmt(3)) c(:)%r, '  ||    circle radius'
       write(16,fmt(3)) c(:)%x, '  ||    circle center x'
       write(16,fmt(3)) c(:)%y, '  ||    circle center y'
       write(16,fmt(3)) c(:)%alpha, '  ||    circle aquifer alpha'
       write(16,fmt(3)) c(:)%bdryQ, '  ||    circle boundry rch rate or head'
       
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

       read(33,*) e(:)%alpha
       if (any(e%alpha < epsilon(0.0))) then
          print *, 'e%alpha must be > 0.0',e%alpha
          stop
       end if

       read(33,*) e(:)%bdryQ
       close(33)
       
       where (e(:)%ibnd == -1 .or. e(:)%ibnd == 0 .or. e(:)%ibnd == +1)
          e(:)%match = .true.
       elsewhere
          e(:)%match = .false.       
       end where

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
       write(16,fmt(3)) e(:)%r, '  ||    ellipse radius (eta)'
       write(16,fmt(3)) e(:)%x, '  ||    ellipse center x'
       write(16,fmt(3)) e(:)%y, '  ||    ellipse center y'
       write(16,fmt(3)) e(:)%f, '  ||    ellipse semi-focal length'
       write(16,fmt(3)) e(:)%theta, '  ||    ellipse angle rotation with +x axis'
       write(16,fmt(3)) e(:)%alpha, '  ||    ellipse aquifer alpha'
       write(16,fmt(3)) e(:)%bdryQ, '  ||     ellipse boundary rch rate or head'
    else
       allocate(e(0),stat=ierr)
       if (ierr /= 0) stop 'ltaem_io.f90 error allocating e(0)'
    end if
    
    ntot = sum(dom%num) ! total number of circular and elliptical elements
    if (ntot < 1) stop 'ltaem_io.f90 Need at least one circular (including well) or&
            & elliptical (including line) element.'
    close(15)
    close(16)
  end subroutine readInput

  !******************************************************
  subroutine writeResults(s)
    use constants, only : DP
    use unsat_type_definitions, only : solution

    type(solution), intent(in) :: s

    character(4) :: chint
    ! adjust the formats of location and results here
    character(6) :: xfmt = 'ES12.4'
    character(9) :: hfmt = 'ES22.14e3'
    integer :: i, j, k

    select case (s%output)
    case (1) 
       ! ** gnuplot contour map friendly output **
       ! print results as x,y,z triplets with the given times separated by double blank lines

       open(unit=20, file=s%outfname, status='replace', action='write')
       write(20,*) '# ltaem contour map output'
       write(20,'(A,I0)') ' # x: ', s%nx
       write(20,'(A,I0)') ' # y: ', s%ny
       write(20,'(A,I0)') ' # xy:', s%nx*s%ny     

       write(20,'(A)')   &
            & '#      X           Y               head&
            &                   velx                    vely&
            &                   phi                     PHI'
       do j = 1, s%ny
          do k = 1, s%nx
             write(20,'(2('//xfmt//',1X),5('//hfmt//',1X))') &
                  & s%x(k), s%y(j), s%h(k,j), s%v(k,j,1:2), s%smphi(k,j), s%lgPHI(k,j)
          end do
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
       write(chint,'(i4.4)') s%nx
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%x(j), j=1,s%nx)
       end do
       close(20)

       ! y-matrix has same column repeated numx times
       open(unit=20, file=trim(s%outfname)//'_y.dat', status='replace', &
            & action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%y(i), j=1,s%nx)
       end do
       close(20)
       
       ! head-matrix
       open(unit=20, file=trim(s%outfname)//'_head.dat', &
            & status='replace', action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%h(j,i), j=1,s%nx)
       end do
       close(20)

       ! velx-matrix
       open(unit=20, file=trim(s%outfname)//'_velx.dat', &
            & status='replace', action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%v(j,i,1), j=1,s%nx)
       end do
       close(20)

       ! vely-matrix
       open(unit=20, file=trim(s%outfname)//'_vely.dat', &
            & status='replace', action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%v(j,i,2), j=1,s%nx)
       end do
       close(20)

       ! phi-matrix
       open(unit=20, file=trim(s%outfname)//'_smphi.dat', &
            & status='replace', action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%smphi(j,i), j=1,s%nx)
       end do
       close(20)

       ! PHI-matrix
       open(unit=20, file=trim(s%outfname)//'_lgPHI.dat', &
            & status='replace', action='write')
       do i = 1, s%ny
          write (20,'('//chint//'(1x,'//hfmt//'))') (s%lgPHI(j,i), j=1,s%nx)
       end do
       close(20)

       write(*,'(/A)') '*********************************************************************'
       write(*,'(3A)') 'matlab output written to ', trim(s%outfname), &
              & '{x,y,head,velx,vely,smphi,lgPHI}.dat'
       write(*,'(A)') '*********************************************************************'
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case default
       write(*,'(A,I0)')  'invalid output code ',s%output
       stop 
    end select
  end subroutine writeResults  

  !##################################################
  subroutine writeGeometry(c,e,s)
    use constants, only : DP
    use unsat_type_definitions, only : circle, ellipse, solution
    
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
end module unsat_file_ops
