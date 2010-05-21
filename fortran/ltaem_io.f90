! this module contains most of the I/O associated with LT-AEM

module file_ops
  implicit none  
  
  private
  public :: readInput, writeResults, writeGeometry

contains

  !##################################################
  ! this routine read the main input file, and allocates the main 
  ! data structures used to store data.
  subroutine readInput(sol,lap,dom,bg,c,e,p)
    use element_specs
    use constants, only : DP
    implicit none

    type(solution), intent(inout) :: sol
    type(INVLT), intent(out) :: lap
    type(particle), intent(out), allocatable :: p(:)

    type(domain), intent(out) :: dom
    type(element), intent(out) :: bg
    type(circle), intent(out), allocatable :: c(:)
    type(ellipse), intent(out), allocatable :: e(:)

    character(128) :: subname = 'ReadInput', echofname
    integer :: ierr, j,nEl,nC,nE  ! #elements, #circles, #ellipses

    echofname = trim(sol%infname) + '.echo'
    open(UNIT=15, FILE=sol%infname, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    if (ierr /= 0) then
       print *, 'READINPUT: error opening input file ',sol%infname
       stop 100
    endif
    open(UNIT=16, FILE=echofname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) then
       print *, 'READINPUT: error opening echo file ',echofname
       stop 101
    endif   
    
    ! solution-specific and background aquifer parameters
    read(15,*) sol%particle, sol%contour, sol%output, &
         & sol%outFname, sol%coeffFName, sol%elemHfName
    read(15,*) bg%por, bg%k, bg%ss, bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb
    read(15,*) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b

    write(16,*) sol%particle, sol%contour, sol%output, &
         & trim(sol%outFname), trim(sol%coeffFName), trim(sol%elemHfName)&
         & '  ||    Lparticle, Lcontour, Ioutput, out/coeff/hierarchy file names'
    write(16,*) bg%por, bg%k, bg%ss, bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb
         & '  ||    por, k, Ss, leaky_type, k2, ss2, b2'
    write(16,*) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b, '  || BGSy, BGKz, unconfined flag, BGb'
    
    ! desired solution points/times
    read(15,*) sol%nx, sol%ny, sol%nt
    allocate(sol%x(sol%nx), sol%y(sol%ny), sol%t(sol%nt))
    read(15,*) sol%x(:)
    read(15,*) sol%y(:)
    do j=1,sol%numt  !! modified to accommidate pest (make time a column)
       read(15,*,iostat=ierr) sol%t(j)
       if(ierr /= 0) write(*,'(A,I0)') 'ERROR reading time ',j
    end do
    if (.not. sol%particle) then
       write(16,*) sol%nx, sol%ny, sol%nt, '  ||    numX, numY, numt'
       write(16,*) sol%x(:), '  ||    xVector'
       write(16,*) sol%y(:), '  ||    yVector'
       write(16,*) sol%t(:), '  ||    tVector'
    endif

    ! inverse Laplace transform parameters
    read(15,*) lap%alpha, lap%tol, lap%m
    if (lap%tol < epsilon(lap%tol)) lap%tol = epsilon(lap%tol)
    write(16,*) lap%alpha, lap%tol, lap%m,'  ||    alpha, tol, M'

    ! circular (includes wells)
    read(15,*) dom%num(1)
    nc = dom%num(1)
    if (nc > 0) then
       allocate(c(nc))
       read(15,*) c(:)%n
       read(15,*) c(:)%m
       read(15,*) c(:)%ibnd
       read(15,*) c(:)%spec
       read(15,*) c(:)%CalcIn
       read(15,*) c(:)%r
       read(15,*) c(:)%x
       read(15,*) c(:)%y
       read(15,*) c(:)%k
       read(15,*) c(:)%Ss
       read(15,*) c(:)%por
       read(15,*) c(:)%area
       do j=1,size(c,dim=1)
          read(15,'(I)', advance='no') c(j)%AreaTime 
          if (c(j)%AreaTime > -1) then
             allocate(c(j)%ATPar(2))
             read(15,*) c(j)%ATPar(:)
             write(15,*) c(j)%AreaTime,c(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for circle ',j
          else
             allocate(c(j)%ATPar(-2*c(j)%AreaTime+1))
             read(15,*) c(j)%ATPar(:)
             write(16,*) c(j)%AreaTime,c(j)%ATPar(:-c(j)%AreaTime+1),' | ',&
                  & c(j)%ATPar(-c(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for circle ',j
          end if
       end do
       do j=1,size(c,dim=1)
          read(15,'(I)', advance='no') c(j)%BdryTime
          if (c(j)%BdryTime > -1) then
             allocate(c(j)%BTPar(2))
             read(15,*) c(j)%BTPar(:)
             write(15,*) c(j)%BdryTime,c(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for circle ',j
          else
             allocate(c(j)%BTPar(-2*c(j)%BdryTime+1))
             read(15,*) c(j)%BTPar(:)
             write(16,*) c(j)%BdryTime,c(j)%BTPar(:-c(j)%BdryTime+1),' | ',&
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
          
       write(16,*) dom%num(1), '  ||   number of circular elements (including wells)'
       write(16,*) c(:)%n,'  ||   number of circular free parameter (Fourier coeffs)'
       write(16,*) c(:)%m,'  ||   number of circular matching locations'
       write(16,*) c(:)%ibnd, '  ||    circle ibnd array'
       write(16,*) c(:)%match, '  ||    circle matching array'
       write(16,*) c(:)%spec, '  ||    circle boundary specified head/flux strength'
       write(16,*) c(:)%calcin, '  ||    calculate inside this circle?'
       write(16,*) c(:)%r, '  ||    circle radius'
       write(16,*) c(:)%x, '  ||    circle center x'
       write(16,*) c(:)%y, '  ||    circle center y'
       write(16,*) c(:)%k, '  ||    circle aquifer k'
       write(16,*) c(:)%ss, '  ||    circle aquifer Ss'
       write(16,*) c(:)%por, '  ||    circle aquifer porosity'
       write(16,*) c(:)%area, '  ||    circle area rch rate'
       write(16,*) c(:)%leakFlag, '  ||     circle leaky type'
       write(16,*) c(:)%aquitardK, '  ||     circle leaky aquitard K'
       write(16,*) c(:)%aquitardSs, '  ||     circle leaky aquitard Ss'
       write(16,*) c(:)%aquitardb, '  ||     circle leaky aquitard thickness'
       write(16,*) c(:)%unconfinedFlag, '  ||     circle unconfined flag'
       write(16,*) c(:)%Sy, '  ||     circle aquifer specific yield'
       write(16,*) c(:)%Kz, '  ||     circle aquifer vertical K'
       write(16,*) c(:)%b,'  ||    circle aquifer thickness'
       write(16,*) c(:)%dskin,'  ||    circle boundary dimensionless skin factor'
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
       read(15,*) e(:)%spec
       read(15,*) e(:)%CalcIn
       read(15,*) e(:)%r
       read(15,*) e(:)%x
       read(15,*) e(:)%y
       read(15,*) e(:)%f
       read(15,*) e(:)%theta      
       read(15,*) e(:)%k
       read(15,*) e(:)%Ss
       read(15,*) e(:)%por
       read(15,*) e(:)%area
       do j=1,size(c,dim=1)
          read(15,'(I)', advance='no') e(j)%AreaTime 
          if (e(j)%AreaTime > -1) then
             allocate(e(j)%ATPar(2))
             read(15,*) e(j)%ATPar(:)
             write(15,*) e(j)%AreaTime,e(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%ATPar(-2*e(j)%AreaTime+1))
             read(15,*) e(j)%ATPar(:)
             write(16,*) e(j)%AreaTime,e(j)%ATPar(:-e(j)%AreaTime+1),' | ',&
                  & e(j)%ATPar(-e(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for ellipse ',j
          end if
       end do
       do j=1,size(c,dim=1)
          read(15,'(I)', advance='no') e(j)%BdryTime
          if (e(j)%BdryTime > -1) then
             allocate(e(j)%BTPar(2))
             read(15,*) e(j)%BTPar(:)
             write(15,*) e(j)%BdryTime,e(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for circle ',j
          else
             allocate(e(j)%BTPar(-2*e(j)%BdryTime+1))
             read(15,*) e(j)%BTPar(:)
             write(16,*) e(j)%BdryTime,e(j)%BTPar(:-e(j)%BdryTime+1),' | ',&
                  & e(j)%BTPar(-e(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for circle ',j
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

       write(16,*) dom%num(2), '  ||   number of elliptical elements (including lines)'
       write(16,*) e(:)%n,'  ||   number of elliptical free parameter (Fourier coeffs)'
       write(16,*) e(:)%m,'  ||   number of ellipse matching locations'
       write(16,*) e(:)%ms,'  ||   size of "infinite" Mathieu matrices'
       write(16,*) e(:)%ibnd, '  ||    ellipse ibnd array'
       write(16,*) e(:)%match, '  ||    ellipse matching array'
       write(16,*) e(:)%spec, '  ||    ellipse boundary specified head/flux strength'
       write(16,*) e(:)%calcin, '  ||    calculate inside this ellipse?'
       write(16,*) e(:)%r, '  ||    ellipse radius (eta)'
       write(16,*) e(:)%x, '  ||    ellipse center x'
       write(16,*) e(:)%y, '  ||    ellipse center y'
       write(16,*) e(:)%x, '  ||    ellipse semi-focal length'
       write(16,*) e(:)%y, '  ||    ellipse angle rotation with +x axis'
       write(16,*) e(:)%k, '  ||    ellipse aquifer k'
       write(16,*) e(:)%ss, '  ||    ellipse aquifer Ss'
       write(16,*) e(:)%por, '  ||    ellipse aquifer porosity'
       write(16,*) e(:)%area, '  ||     ellipse area rch rate'
       write(16,*) e(:)%leakFlag, '  ||     ellipse leaky type'
       write(16,*) e(:)%aquitardK, '  ||     ellipse leaky aquitard K'
       write(16,*) e(:)%aquitardSs, '  ||     ellipse leaky aquitard Ss'
       write(16,*) e(:)%aquitardb, '  ||     ellipse leaky aquitard thickness'
       write(16,*) e(:)%unconfinedFlag, '  ||     ellipse unconfined flag'
       write(16,*) e(:)%Sy, '  ||     ellipse aquifer specific yield'
       write(16,*) e(:)%Kz, '  ||     ellipse aquifer vertical K'
       write(16,*) e(:)%b,'  ||    ellipse aquifer thickness'
       write(16,*) e(:)%dskin,'  ||    ellipse boundary dimensionless skin factor'
    else
       allocate(e(0))
    end if
    
    nEl = sum(dom%num) ! total number of circular and elliptical elements
    if (nEl < 1) then
       print *, 'need at least one circular (including well) or'\\&
            & 'elliptical (including line) element.'
       stop 102
    end if

    ! compute secondary parameters
    c(:)%alpha = c(:)%K/c(:)%Ss
    e(:)%alpha = e(:)%K/e(:)%Ss

    write(16,*) c(:)%alpha,'  ||    circle hydraulic diffusivity'
    write(16,*) e(:)%alpha,'  ||    ellipse hydraulic diffusivity'

    ! re-calculation parameter
    read(15,*) sol%calc
    write(16,*) sol%calc, '  || re-calculate coefficients?'

    ! particles
    if (sol%particle) then
       read(15,*) sol%nPart
       allocate(p(sol%nPart))
       read(15,*) p(:)%tol 
       read(15,*) p(:)%dt 
       read(15,*) p(:)%maxStep
       read(15,*) p(:)%streakSkip
       read(15,*) p(:)%x
       read(15,*) p(:)%y
       read(15,*) p(:)%ti
       read(15,*) p(:)%tf
       read(15,*) p(:)%int
       read(15,*) p(:)%InclIn

       write(16,*) p(:)%tol,'  ||    particle solution tolerances'
       write(16,*) p(:)%dt,'  ||    particle dt'
       write(16,*) p(:)%maxStep,'  ||   particle max flux'
       write(16,*) p(:)%x, '  ||    particle initial x'
       write(16,*) p(:)%y, '  ||    particle initial y'
       write(16,*) p(:)%ti, '  ||    particle initial t'
       write(16,*) p(:)%tf, '  ||    particle maximum t'
       write(16,*) p(:)%int, '  ||    particle integration method'
       write(16,*) p(:)%InclIn, '  ||    particle begins inside CH/CF incl?'
    else
       allocate(p(0))
    endif
    close(15)
    close(16)
  end subroutine readInput

  !******************************************************
  subroutine writeResults(head,velx,vely,x,y,t,flag,filename,pout)
    use constants, only : DP
    use element_specs, only : parstreakskip
    use error_handler, only : fileerror
    implicit none
    
    character(2) :: cht
    character(4) :: chx
    character(25) :: fmtstr
    character(128) :: subname = 'WriteResults'
    real(DP), intent(in), dimension(:,:,:) :: head,velx,vely, pout
    real(DP), intent(in), dimension(:) :: x, y, t
    integer, intent(in) :: flag
    character(128), intent(in) :: filename
    integer :: ierr, i, j, k, numt, numx, numy

    numx = size(x,1); numy = size(y,1); numt = size(t,1)

    select case (flag)
    case (1) 
       ! ** gnuplot contour map friendly output **
       ! print results as x,y,z triplets with the given times separated by double blank lines

       open(UNIT=20, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(filename,ierr,subname,0)
       
       write (20,*) '# ltaem contour map output'
       write (20,234) ' # t:', numt
       write (20,234) ' # x:', numx
       write (20,234) ' # y:', numy
       write (20,234) ' #xy:', numx*numy     

       do i = 1, numt
          write (20,994) ' # t=',t(i)
          write (20,*)   &
          & '#      X          Y              head                velx                vely'
          write (20,995) ((x(k),y(j),head(k,j,i),velx(k,j,i),vely(k,j,i),k=1,numx),j=1,numy)
          write (20,*)  '  '
          write (20,*)  '  ' 
       end do       

       write(20,*) '# EOF'

       close(20)

       print *, '  '
       print *, '***********************************************************'
       print *, ' gnuplot contour map style output written to ', trim(filename)
       print *, '***********************************************************'

234    format (A5,i5) 
994    format (A5,ES11.5)
995    format (1X,ES13.5,1X,ES13.5,1X,ES22.14e3,1X,ES22.14e3,1X,ES22.14e3)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case (2)
       
       ! ** matlab-friendly output **
       ! print results as matricies, with each variable (at each time) going to a separate file
       ! (x and y matricies - similar to results from matlab function meshgrid)
       
       ! x-matrix has same row repeated numy times
       open(UNIT=20, FILE=trim(filename)//'x.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(trim(filename)//'x.dat',ierr,subname,0)
       
       write(chx,'(i4)') numx
       fmtstr = '('//chx//'(1x,ES22.14e3))'

       do i = 1, numy
          write (20,fmtstr) (x(j), j=1,numx)
       end do
       close(20)

       ! y-matrix has same column repeated numx times
       open(UNIT=20, FILE=trim(filename)//'y.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(trim(filename)//'y.dat',ierr,subname,0)
       
       do i = 1, numy
          write (20,fmtstr) (y(i), j=1,numx)
       end do
       close(20)
       
       do k = 1, numt
          write(cht,'(i2.2)') k

          ! head-matrix
          open(UNIT=20, FILE=trim(filename)//'head'//cht//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
          if (ierr /= 0) call fileError(trim(filename)//'head'//cht//'.dat',ierr,subname,0)

          do i = 1, numy
             write (20,fmtstr) (head(j,i,k), j=1,numx)
          end do
          close(20)

          ! velx-matrix
          open(UNIT=20, FILE=trim(filename)//'velx'//cht//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
          if (ierr /= 0) call fileError(trim(filename)//'velx'//cht//'.dat',ierr,subname,0)

          do i = 1, numy
             write (20,fmtstr) (velx(j,i,k), j=1,numx)
          end do
          close(20)

          ! vely-matrix
          open(UNIT=20, FILE=trim(filename)//'vely'//cht//'.dat', &
               & STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
          if (ierr /= 0) call fileError(trim(filename)//'vely'//cht//'.dat',ierr,subname,0)

          do i = 1, numy
             write (20,fmtstr) (vely(j,i,k), j=1,numx)
          end do
          close(20)
       end do

       ! column of calculation times
       open(UNIT=20, FILE=trim(filename)//'t.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(trim(filename)//'t.dat',ierr,subname,0)
       
       write (20,*) (t(j), j=1,numt)
       close(20)

       print *, ' '
       print *, '*********************************************************************'
       print *, 'matlab output written to ', trim(filename), &
              & '{x,y,t,head{1-n},velx{1-n},vely{1-n}}.dat'
       print *, '*********************************************************************'
       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (3)

       ! ** gnuplot hydrograph-friendly output **
       ! column of time values at a location through time
       ! locations separated by blank lines
       open(UNIT=20, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(filename,ierr,subname,0)
       
       write (20,*) '# ltaem hydrograph output'
      
       do i = 1, numx
          write (20,555) ' # location: x=',x(i),' y=',y(i)
          write (20,*)   '#     time              head                  velx                vely'
          do k = 1, numt
             write (20,556) t(k),head(i,1,k),velx(i,1,k),vely(i,1,k)
          end do
          write (20,*)  '  '
          write (20,*)  '  '
       end do       

       write(20,*) '# EOF'

       close(20)
       print *, ' '
       print *, '***********************************************************'
       print *, 'gnuplot style output written to ', trim(filename)
       print *, '***********************************************************'

555    format (2(A,ES12.5))
556    format (1X,ES13.5,3(1X,ES22.14e3))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (11)

       ! ** gnuplot hydrograph-friendly output **
       ! column of time values at a location through time
       ! locations separated by blank lines (grid no velocity)
       open(UNIT=20, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(filename,ierr,subname,0)
       
       write (20,*) '# ltaem hydrograph output'
      
       do i = 1, numy
          do j = 1, numx
             write (20,155) ' # location: x=',x(j),' y=',y(i)
             write (20,*)   '#     time              head                  velx                vely'
             do k = 1, numt
                write (20,156) t(k),head(j,i,k),velx(j,i,k),vely(j,i,k)
             end do
             write (20,*)  '  '
             write (20,*)  '  '
          end do
       end do
       
       write(20,*) '# EOF'

       close(20)
       print *, ' '
       print *, '***********************************************************'
       print *, 'gnuplot style output written to ', trim(filename)
       print *, '***********************************************************'

155    format (2(A,ES12.5))
156    format (1X,ES13.5,3(1X,ES22.14e3))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (10)

       ! ** output for inverse in Matlab
       ! column of time values at a location through time
       ! locations separated by blank lines

       do i = 1, numx
          write(chx,'(I4.4)') i
          open(UNIT=20, FILE=trim(filename)//chx, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
          if (ierr /= 0) call fileError(filename,ierr,subname,0)
          do k = 1, numt
             write (20,'(2(ES18.10,1X))') t(k),head(i,1,k)
          end do
          write(20,*)  '  '
          write(20,*)  '  '
       close(20)
       end do       

       print *, ' '
       print *, '***********************************************************'
       print *, 'inverse output written to ', trim(filename) , '0000-',chx
       print *, '***********************************************************'


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (4)

       ! ** pathline gnuplot-style output **
       ! columns of time values for starting locations
       ! particles separated by blank lines
       open(UNIT=20, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(filename,ierr,subname,0)
       
       write (20,*) '# ltaem particle tracking output'


       do i = 1, size(pout,dim=3)
          numt = count(pout(:,1,i) /= 0)
          write (20,*) '# particle:',i
          write (20,*)   &
          & '#     time          x            y             velx          vely '
          write (20,334) ((pout(k,j,i),j=1,5),k=1,numt)
          write (20,*)  '  '
          write (20,*)  '  '
       end do       

       write(20,*) '# EOF'

       close(20)
       print *, ' '
       print *, '***********************************************************'
       print *, 'particle tracking output written to ', trim(filename)
       print *, '***********************************************************'

334    format (5(1X,ES13.5))
       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (5)

       ! ** streakline gnuplot-style output **
       ! this requires constant time steps (can't use adaptive integration)
       ! each block is a requested time, each row a particle

       open(UNIT=90, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
       if (ierr /= 0) call fileError(filename,ierr,subname,0)
       
       write (90,*) '# ltaem particle tracking streakfile output'

       ! max number of times for all particles
       numt = maxval(count(pout(:,1,:) /= 0,dim=1))

       do i = 1, numt, PARstreakSkip
          write (90,*) '# time:',maxval(pout(i,1,:)) ! to ensure a non-zero time is reported
          write (90,*) '#  particle       x            y  '
          do j = 1, size(pout,dim=3)
             if (pout(i,1,j) > 0) then
                ! only write particle if it has non-zero data
                write (90,442)  j,pout(i,2:3,j)
             end if
          end do
          write (90,*)  '  '
          write (90,*)  '  '
       end do       

       write(90,*) '# EOF'

       close(90)
       print *, ' '
       print *, '***********************************************************'
       print *, 'particle tracking streakfile written to ', trim(filename)
       print *, '***********************************************************'

442    format (I3,1X,2(1X,ES13.5))

    case default
       print *, 'invalid output code', flag
    end select
    
  end subroutine writeResults  

  !##################################################
  subroutine writeGeometry(EXm,EYm,Wx,Wy,Wr,Wq)
    use constants, only : DP
    use error_handler, only : fileError
    implicit none
    
    real(DP), dimension(:,:), intent(in) :: EXm, EYm
    real(DP), dimension(:), intent(in) :: Wx, Wy, Wr, Wq
    character(128) :: subname = 'writeGeometry', inclfname = 'circles.dat', wellfname = 'well.dat'
    integer :: ierr, m, nm, ni, nw, incl, well
    
    ni = size(EXm,2)
    nm = size(EXm,1)
    nw = size(Wx,1)

    ! write inclusions to file
    open(UNIT=40, FILE=inclfname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) call fileError(inclfname,ierr,subname,0)
    
    write(40,*) '# ',nm,'points along circumference of',ni,'circular elements'
    do incl = 1,ni
       write(40,*) '# element ',incl
       do m = 1,nm
          write(40,666)  EXm(m,incl),EYm(m,incl)
       end do
       write(40,666)  EXm(1,incl),EYm(1,incl) ! joins circle back up with beginning
       write(40,*)    
       write(40,*)    
    end do

    write(40,*) '# EOF'

666 format (2(1x,ES13.6))
    close(40)
    
    ! write wells to file
    open(UNIT=40, FILE=wellfname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) call fileError(wellfname,ierr,subname,0)

    write(40,*) '# ',nw,'wells'
    write(40,*) '# id       x            y          well_r          well_Q'
    do well = 1,nw
       write(40,777)  well,Wx(well),Wy(well),Wr(well),Wq(well)
    end do
    
    write(40,*) '# EOF'
    
777 format (1x,I3,4(1x,ES13.6))
    close(40)

  end subroutine writeGeometry
end module file_ops
