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
    use constants, only : DP, lenFN
    implicit none

    type(solution), intent(inout) :: sol
    type(particle), intent(out), allocatable :: p(:)
    type(domain), intent(out) :: dom
    type(element), intent(out) :: bg
    type(circle), intent(out), allocatable :: c(:)
    type(ellipse), intent(out), allocatable :: e(:)

    character(4) :: chint
    character(20), dimension(3) :: fmt
    character(46) :: lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
    character(lenFN+5) :: echofname
    integer :: ierr,j,ntot,nC,nE  ! #elements, #circles, #ellipses

    open(UNIT=15, FILE=sol%infname, STATUS='OLD', ACTION='READ')
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening input file ',sol%infname
       stop 100
    endif

    echofname = trim(sol%infname)//'.echo'
    open(UNIT=16, FILE=echofname, STATUS='REPLACE', ACTION='WRITE')
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening echo file ',echofname
       stop 101
    endif   
    
    ! solution-specific and background aquifer parameters
    read(15,*) sol%particle, sol%contour, sol%output, &
         & sol%outFname, sol%coeffFName, sol%elemHfName, sol%geomfName
    read(15,*) bg%por, bg%k, bg%ss, bg%leakFlag, &
         & bg%aquitardK, bg%aquitardSs, bg%aquitardb
    read(15,*) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b

    write(16,'(2(1L,1X),I0,5(1X,A))') sol%particle, sol%contour, sol%output, &
         & trim(sol%outFname), trim(sol%coeffFName), trim(sol%elemHfName), &
         & trim(sol%geomFname),'  ||    particle?, contour?, output,'// &
         & 'out/coeff/hierarchy/geometry file names'
    write(16,'(3(ES11.5,1X),1L,3(1X,ES11.5),A)') bg%por, bg%k, bg%ss, &
         & bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb, & 
         & '  ||    por, k, Ss, leaky flag, K2, Ss2, b2'
    write(16,'(2(ES11.5,1X),1L,1X,ES11.5,A)') bg%Sy, bg%kz, bg%unconfinedFlag, &
         & bg%b, '  || Sy, Kz, unconfined?, BGb'
    
    ! desired solution points/times
    read(15,*) sol%nx, sol%ny, sol%nt
    allocate(sol%x(sol%nx), sol%y(sol%ny), sol%t(sol%nt))
    read(15,*) sol%x(:)
    read(15,*) sol%y(:)
    do j=1,sol%nt  !! modified to accommidate pest (make time a column)
       read(15,*,iostat=ierr) sol%t(j)
       if(ierr /= 0) write(*,'(A,I0)') 'ERROR reading time from input',j
    end do
    if (.not. sol%particle) then
       write(16,'(3(I0,1X),A)') sol%nx, sol%ny, sol%nt, '  ||    numX, numY, numt'
       write(chint,'(I4.4)') sol%nx
       write(16,'('//chint//'(ES12.5,1X),A)') sol%x(:), '  ||    x Vector'
       write(chint,'(I4.4)') sol%ny
       write(16,'('//chint//'(ES12.5,1X),A)') sol%y(:), '  ||    y Vector'
       write(chint,'(I4.4)') sol%nt
       write(16,'('//chint//'(ES11.5,1X),A)') sol%t(:), '  ||    t Vector'
    endif

    ! inverse Laplace transform parameters
    read(15,*) sol%alpha, sol%tol, sol%m
    if (sol%tol < epsilon(sol%tol)) then
       sol%tol = epsilon(sol%tol)
       write(*,'(A,ES11.5)') 'WARNING: increased INVLAP solution tolerance to ',sol%tol 
    end if
    write(16,'(2(ES11.5,1X),I0,A)') sol%alpha, sol%tol, sol%m,'  ||    alpha, tol, M'

    ! circular (includes wells)
    read(15,*) dom%num(1)
    nc = dom%num(1)
    if (nc > 0) then
       allocate(c(nc))
       read(15,*) c(:)%n
       read(15,*) c(:)%m
       read(15,*) c(:)%ibnd
       read(15,*) c(:)%CalcIn
       read(15,*) c(:)%r
       read(15,*) c(:)%x
       read(15,*) c(:)%y
       read(15,*) c(:)%k
       read(15,*) c(:)%Ss
       read(15,*) c(:)%por
       read(15,*) c(:)%areaQ ! area source strength (flux)
       read(15,*) c(:)%bdryQ ! streng of specified value on bdry (head or flux) 
       do j=1,size(c,dim=1)
          read(15,'(I3)', advance='no') c(j)%AreaTime 
          if (c(j)%AreaTime > -1) then
             ! functional time behavior
             allocate(c(j)%ATPar(2))
             read(15,*) c(j)%ATPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%AreaTime,c(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for circle ',j
          else
             ! piecewise-constant time behavior 
             allocate(c(j)%ATPar(-2*c(j)%AreaTime+1))
             read(15,*) c(j)%ATPar(:)
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
             write(lfmt(8:11),'(I4.4)')  size(c(j)%ATPar(:-c(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%ATPar(-c(j)%AreaTime+2:),1)
             write(16,lfmt) c(j)%AreaTime,c(j)%ATPar(:-c(j)%AreaTime+1),'| ',&
                  & c(j)%ATPar(-c(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for circle ',j
          end if
       end do
       do j=1,size(c,dim=1)
          read(15,'(I3)', advance='no') c(j)%BdryTime
          if (c(j)%BdryTime > -1) then
             allocate(c(j)%BTPar(2))
             read(15,*) c(j)%BTPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%BdryTime,c(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for circle ',j
          else
             allocate(c(j)%BTPar(-2*c(j)%BdryTime+1))
             read(15,*) c(j)%BTPar(:)
             write(lfmt(8:11),'(I4.4)')  size(c(j)%BTPar(:-c(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%BTPar(-c(j)%BdryTime+2:),1)
             write(16,lfmt) c(j)%BdryTime,c(j)%BTPar(:-c(j)%BdryTime+1),' | ',&
                  & c(j)%BTPar(-c(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for circle ',j
          end if
       end do
       read(15,*) c(:)%leakFlag
       read(15,*) c(:)%aquitardK
       read(15,*) c(:)%aquitardSs
       read(15,*) c(:)%aquitardb  !aquitard thickness
       read(15,*) c(:)%unconfinedFlag
       read(15,*) c(:)%Sy
       read(15,*) c(:)%Kz
       read(15,*) c(:)%b  ! aquifer thickness
       read(15,*) c(:)%dskin ! dimensionless skin
   
       where (c(:)%ibnd == -1 .or. c(:)%ibnd == 0 .or. c(:)%ibnd == +1)
          c(:)%match = .true.
       elsewhere
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
    else
       allocate(c(0))
    end if

    ! elliptical (includes line sources/sinks)
    read(15,*) dom%num(2)
    ne = dom%num(2)
    if (ne > 0) then
       allocate(e(ne))
       read(15,*) e(:)%n
       read(15,*) e(:)%m
       read(15,*) e(:)%ms
       read(15,*) e(:)%ibnd
       read(15,*) e(:)%CalcIn
       read(15,*) e(:)%r   ! eta
       read(15,*) e(:)%x
       read(15,*) e(:)%y
       read(15,*) e(:)%f
       read(15,*) e(:)%theta      
       read(15,*) e(:)%k
       read(15,*) e(:)%Ss
       read(15,*) e(:)%por
       read(15,*) e(:)%areaQ
       read(15,*) e(:)%bdryQ
       do j=1,size(e,dim=1)
          read(15,'(I3)', advance='no') e(j)%AreaTime 
          if (e(j)%AreaTime > -1) then
             allocate(e(j)%ATPar(2))
             read(15,*) e(j)%ATPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%AreaTime,e(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%ATPar(-2*e(j)%AreaTime+1))
             read(15,*) e(j)%ATPar(:)
             write(lfmt(8:11),'(I3.3)')  size(e(j)%ATPar(:-e(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%ATPar(-e(j)%AreaTime+2:),1)
             write(16,lfmt) e(j)%AreaTime,e(j)%ATPar(:-e(j)%AreaTime+1),' | ',&
                  & e(j)%ATPar(-e(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for ellipse ',j
          end if
       end do
       do j=1,size(e,dim=1)
          read(15,'(I3)', advance='no') e(j)%BdryTime
          if (e(j)%BdryTime > -1) then
             allocate(e(j)%BTPar(2))
             read(15,*) e(j)%BTPar(:)
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%BdryTime,e(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%BTPar(-2*e(j)%BdryTime+1))
             read(15,*) e(j)%BTPar(:)
             write(lfmt(8:11),'(I3.3)')  size(e(j)%BTPar(:-e(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%BTPar(-e(j)%BdryTime+2:),1)
             write(16,lfmt) e(j)%BdryTime,e(j)%BTPar(:-e(j)%BdryTime+1),' | ',&
                  & e(j)%BTPar(-e(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for ellipse ',j
          end if
       end do
       read(15,*) e(:)%leakFlag
       read(15,*) e(:)%aquitardK
       read(15,*) e(:)%aquitardSs
       read(15,*) e(:)%aquitardb  !aquitard thickness
       read(15,*) e(:)%unconfinedFlag
       read(15,*) e(:)%Sy
       read(15,*) e(:)%Kz
       read(15,*) e(:)%b  ! aquifer thickness
       read(15,*) e(:)%dskin ! dimensionless skin

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
       allocate(e(0))
    end if
    
    ntot = sum(dom%num) ! total number of circular and elliptical elements
    if (ntot < 1) then
       write(*,'(A)') 'need at least one circular (including well) or'//&
            & 'elliptical (including line) element.'
       stop 102
    end if

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

    ! re-calculation parameter
    read(15,*) sol%calc
    write(16,'(1L,A)') sol%calc, '  || re-calculate coefficients?'

    ! particles
    if (sol%particle) then
       read(15,*) sol%nPart,  sol%streakSkip
       allocate(p(sol%nPart))
       read(15,*) p(:)%tol 
       read(15,*) p(:)%dt 
       read(15,*) p(:)%maxStep
       read(15,*) p(:)%x
       read(15,*) p(:)%y
       read(15,*) p(:)%ti
       read(15,*) p(:)%tf
       read(15,*) p(:)%int
       read(15,*) p(:)%InclIn

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
       allocate(p(0))
    endif
    close(15)
    close(16)
  end subroutine readInput

  !******************************************************
  subroutine writeResults(s,p)
    use constants, only : DP
    use type_definitions, only : solution, particle
    implicit none

    type(solution), intent(in) :: s
    type(particle), dimension(:), intent(in) :: p

    character(4), dimension(2) :: chint
    ! adjust the formats of time, location, and results here
    character(6) :: tfmt = 'ES13.5', xfmt = 'ES12.4'
    character(9) :: hfmt = 'ES22.14e3'
    
    integer :: i, j, k, nt, nx, ny

    select case (s%output)
    case (1) 
       ! ** gnuplot contour map friendly output **
       ! print results as x,y,z triplets with the given times separated by double blank lines

       open(UNIT=20, FILE=s%outFname, STATUS='REPLACE', ACTION='WRITE')
       write(20,*) '# ltaem contour map output'
       write(20,'(A,I0)') ' # t: ', nt
       write(20,'(A,I0)') ' # x: ', nx
       write(20,'(A,I0)') ' # y: ', ny
       write(20,'(A,I0)') ' #xy: ', nx*ny     

       do i = 1, nt
          write(20,'(A,'//tfmt//')') ' # t= ',s%t(i)
          write(20,'(A)')   &
          & '#      X          Y              head'//&
          & '               velx                vely'
          write(20,'(2('//xfmt//',1X),3('//hfmt//',1X))') &
               & ((s%x(k), s%y(j), s%h(k,j,i), s%v(k,j,i,1:2), k=1,nx), j=1,ny)
          write(20,'(//)')
       end do       
       write(20,'(A)') '# EOF'
       close(20)

       write(*,'(/A)') '***********************************************************'
       write(*,'(2A)') ' gnuplot contour map style output written to ', trim(s%outfname)
       write(*,'(A)') '***********************************************************'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (2)
       
       ! ** matlab-friendly output **
       ! print results as matricies, with each variable (at each time) going to a separate file
       ! (x and y matricies - similar to results from matlab function meshgrid)
       
       ! x-matrix has same row repeated numy times
       open(UNIT=20, FILE=trim(s%outFname)//'_x.dat', STATUS='REPLACE', &
            & ACTION='WRITE')
       write(chint(1),'(i4.4)') nx
       do i = 1, ny
          write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%x(j), j=1,nx)
       end do
       close(20)

       ! y-matrix has same column repeated numx times
       open(UNIT=20, FILE=trim(s%outFname)//'_y.dat', STATUS='REPLACE', &
            & ACTION='WRITE')
       do i = 1, ny
          write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%y(i), j=1,nx)
       end do
       close(20)
       
       do k = 1, nt
          write(chint(2),'(i4.4)') k

          ! head-matrix
          open(UNIT=20, FILE=trim(s%outFname)//'_head_'//chint(2)//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE')
          do i = 1, ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%h(j,i,k), j=1,nx)
          end do
          close(20)

          ! velx-matrix
          open(UNIT=20, FILE=trim(s%outfname)//'_velx_'//chint(2)//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE')
          do i = 1, ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(j,i,k,1), j=1,nx)
          end do
          close(20)

          ! vely-matrix
          open(UNIT=20, FILE=trim(s%outfname)//'_vely_'//chint(2)//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE')
          do i = 1, ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(j,i,k,2), j=1,nx)
          end do
          close(20)
       end do

       ! column of calculation times
       open(UNIT=20, FILE=trim(s%outfname)//'_t.dat', STATUS='REPLACE', &
            & ACTION='WRITE')
       write (20,'('//tfmt//')') (s%t(j), j=1,nt)
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
       open(UNIT=20, FILE=s%outfname, STATUS='REPLACE', ACTION='WRITE')
       write (20,'(A)') '# ltaem hydrograph output'      
       do i = 1, nx
          write (20,'(2(A,'//xfmt//'))') ' # location: x=',s%x(i),' y=',s%y(i)
          write (20,'(A)')   '#     time              head&
               &                  velx                vely'
          do k = 1, nt
             write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                  & s%t(k),s%h(i,1,k),s%v(i,1,k,1:2)
          end do
          write (20,'(//)')
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
       open(UNIT=20, FILE=s%outfname, STATUS='REPLACE', ACTION='WRITE')
       write (20,'(A)') '# ltaem hydrograph output'
       do i = 1, ny
          do j = 1, nx
             write (20,'(2(A,'//xfmt//'))') ' # location: x=',s%x(j),' y=',s%y(i)
             write (20,'(A)')   '#     time              head&
                  &                  velx                vely'
             do k = 1, nt
                write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                     & s%t(k),s%h(j,i,k),s%v(j,i,k,1:2)
             end do
             write (20,'(//)')
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

       do i = 1, nx
          write(chint(1),'(I4.4)') i
          open(UNIT=20, FILE=trim(s%outfname)//'_'//chint(1), STATUS='REPLACE', &
               &ACTION='WRITE')
          do k = 1, nt
             write (20,'('//tfmt//',1X,'//hfmt//')') s%t(k),s%h(i,1,k)
          end do
          write(20,'(//)')
          close(20)
       end do       

       write(*,'(/A)') '***********************************************************'
       write(*,'(4A)') 'inverse output written to ', trim(s%outfname) , '0000-',chint(1)
       write(*,'(A)') '***********************************************************'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (4)

       ! ** pathline gnuplot-style output **
       ! columns of time values for starting locations
       ! particles separated by blank lines
       open(UNIT=20, FILE=s%outFname, STATUS='REPLACE', ACTION='WRITE')
       write (20,'(A)') '# ltaem particle tracking output'
       do i = 1, size(p,dim=1)
          write (20,'(A,I0)') '# particle: ',i
          write (20,'(A)')   &
          & '#     time          x            y             velx          vely '
          do k=1,p(i)%numt
             write (20,'('//tfmt//',2'//xfmt//',2'//hfmt//')') &
                  & p(i)%result(1:5,k)
          end do
          write (20,'(//)')
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

       open(UNIT=90, FILE=s%outfname, STATUS='REPLACE', ACTION='WRITE')
       write (90,'(A)') '# ltaem particle tracking streakfile output'

       ! max number of times for all particles
       nt = maxval(p(:)%numt,dim=1)

       do i = 1, nt, s%streakSkip
          ! use maxval to ensure a non-zero time is reported
          write (90,'(A'//tfmt//')') '# time:', maxval(p(i)%result(1,:)) 
          write (90,'(A)') '#  particle       x            y'//&
               &'           velx         vely'
          do j = 1, size(p,dim=1)
             if (p(j)%result(1,i) > 0.0) then
                ! only write particle if it has non-zero data
                write (90,'(I0,4(1X'//hfmt//'))')  j,p(j)%result(2:5,i)
             end if
          end do
          write (90,'(//)')
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
    implicit none
    
    type(circle), dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: s
    integer :: nc, ne, ntot, i, j, ierr

    nc = size(c,1); ne = size(e,1)
    ntot = nc + ne

    ! write matching points to file
    open(UNIT=40, FILE=s%geomFname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) then
       write(*,'(2A)') 'error opening output file for writing &
            &element matching locations ',s%geomFname
       stop
    end if
    
    write(40,'(A)') '# points along circumference of circular and elliptical elements'
    do i = 1,nc
       write(40,'(A,I0)') '# circular element ',i
       do j = 1,c(i)%M
          write(40,'(2(ES13.5,1X))')  c(i)%Zom(j)
       end do
       if (c(i)%M > 1) then
          ! joins circle back up with beginning for plotting
          write(40,'(2(ES13.5,1X))')  c(i)%Zom(j)
       end if
       write(40,'(//)')    
    end do

    do i = 1,ne
       write(40,'(A,I0)') '# elliptical element ',i
       do j = 1,e(i)%M
          write(40,'(2(ES13.5,1X))')  e(i)%Zom(j)
       end do
       if (e(i)%M > 1) then
          ! joins ellipse back up with beginning for plotting
          write(40,'(2(ES13.5,1X))')  e(i)%Zom(j)
       end if
       write(40,'(//)')    
    end do
    write(40,'(A)') '# EOF'
    close(40)
  end subroutine writeGeometry
end module file_ops
