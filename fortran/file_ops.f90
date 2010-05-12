! this module contains most of the I/O associated with LT-AEM

module file_ops
  implicit none  

  private
  public :: readInput, writeResults, writeGeometry, readEllipseInput

  contains

  !##################################################
  subroutine readInput(sol,lap,dom,bg,c,e,p)
    use element_specs
    use error_handler, only : fileerror
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
    integer :: ierr, j

    echofname = trim(sol%infname) + '.echo'

    open(UNIT=15, FILE=sol%infname, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    if (ierr /= 0) call fileError(filename,ierr,subname,1)
    
    open(UNIT=16, FILE=echofname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) call fileError(echofname,ierr,subname,0)
    
    ! solution-specific and background aquifer parameters
    read(15,*) sol%particle, sol%contour, sol%output, sol%outFname, sol%coeffFName
    read(15,*) bg%por, bg%k, bg%ss, bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb
    read(15,*) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b

    write(16,*) sol%particle, sol%contour, sol%output, trim(sol%outFname), trim(sol%coeffFName), &
         & '  ||    Lparticle, Lcontour, Ioutput, out/coeff fnames'
    write(16,*) bg%por, bg%k, bg%ss, bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb
         & '  ||    por, k, Ss, leaky_type, k2, ss2, b2'
    write(16,*) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b, '  || BGSy, BGKz, unconfined flag, BGb'
    
    ! desired solution points/times
    read(15,*) sol%numx, sol%numy, sol%numt
    allocate(sol%x(sol%numx), sol%y(sol%numy), sol%t(sol%numt))
    read(15,*) sol%x(:)
    read(15,*) sol%y(:)
    do j=1,sol%numt  !! modified to accommidate pest (make time a column)
       read(15,*,iostat=ierr) sol%t(j)
       if(ierr /= 0) write(*,'(A,I0)') 'ERROR reading time ',j
    end do
    if (.not. sol%particle) then
       write(16,*) sol%numx, sol%numy, sol%numt, '  ||    numX, numY, numt'
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
    allocate(c(dom%num(1)))
    read(15,*) c(:)%n
    read(15,*) c(:)%m
    read(15,*) c(:)%matchTol
    read(15,*) c(:)%matchOmega
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
    read(15,*) c(:)%Areatime
    read(15,*) c(:)%Atpar(1)  
    read(15,*) c(:)%Atpar(2)
    read(15,*) c(:)%Bdrytime(1:c(:)%num)
    read(15,*) c(:)%Btpar(1:c(:)%num,1)  
    read(15,*) c(:)%Btpar(1:c(:)%num,2)
    read(15,*) c(:)%aquitardLeak(1:c(:)%num) !! new leaky circle stuff
    read(15,*) c(:)%aquitardK(1:c(:)%num)
    read(15,*) c(:)%aquitardSs(1:c(:)%num)
    read(15,*) c(:)%aquitardb(1:c(:)%num)
    read(15,*) c(:)%unconfined(1:c(:)%num)  !! unconfined stuff
    read(15,*) c(:)%Sy(1:c(:)%num)
    read(15,*) c(:)%Kz(1:c(:)%num)

    CImatch = .false.
    where (CIibnd == -1 .or. CIibnd == 0 .or. CIibnd == +1)
       CImatch = .true.
    end where
       
    write(16,*) CInum, CIn, CIm, CImatchTol, CImatchOmega, &
         & '  ||    number of inclusions, N, M; match_tol, match_omega'
    write(16,*) CIibnd(1:CInum), '  ||    ibnd array'
    write(16,*) CImatch(1:CInum), '  ||    matching array'
    write(16,*) CIspec(1:CInum), '  ||    specified quantity'
    write(16,*) CIcalcin(1:CInum), '  ||    calculate inside this element?'
    write(16,*) CIr(1:CInum), '  ||    incl radius'
    write(16,*) CIx(1:CInum), '  ||    incl ctr x'
    write(16,*) CIy(1:CInum), '  ||    incl ctr y'
    write(16,*) CIk(1:CInum), '  ||    incl k'
    write(16,*) CIss(1:CInum), '  ||    incl Ss'
    write(16,*) CIpor(1:CInum), '  ||    incl porosity'
    write(16,*) CIarea(1:CInum), '  ||    incl area rch rate'
    write(16,*) CIAreaTime(1:CInum), '  ||    incl Area flux time fcn index'
    write(16,*) CIAtpar(1:CInum,1), '  ||    incl Area flux time fcn parameter "a"'
    write(16,*) CIAtpar(1:CInum,2), '  ||    incl Area flux time fcn parameter "b"'
    write(16,*) CIBdryTime(1:CInum), '  ||    incl Boundary head/flux time fcn index'
    write(16,*) CIBtpar(1:CInum,1), '  ||    incl Boundary head/flux time fcn parameter "a"'
    write(16,*) CIBtpar(1:CInum,2), '  ||    incl Boundary head/flux time fcn parameter "b"'
    write(16,*) CIaquitardLeak(1:CInum), '  ||     leaky type'
    write(16,*) CIaquitardK(1:CInum), '  ||     leaky aquitard K'
    write(16,*) CIaquitardSs(1:CInum), '  ||     leaky aquitard Ss'
    write(16,*) CIaquitardb(1:CInum), '  ||     leaky aquitard thickness'
    write(16,*) CIunconfined(1:CInum), '  ||     unconfined flag'
    write(16,*) CISy(1:CInum), '  ||     specific yield'
    write(16,*) CIKz(1:CInum), '  ||     vertical K'

    ! build up k,s, por and alpha vectors
    allocate(kv(0:CInum),sv(0:CInum),av(0:CInum),porv(0:CInum),k2v(0:CInum),&
         & S2v(0:CInum),a2v(0:CInum),leakv(0:CInum),b2v(0:CInum),&
         & syv(0:CInum),Kzv(0:CInum),unconfv(0:CInum))

    kv(0) = BGk; sv(0) = BGss; av(0) = BGk/BGss; porv(0) = BGpor;
    k2v(0) = BGk2; s2v(0) = BGs2; a2v(0) = BGk2/BGS2; b2v(0) = BGb2
    syv(0) = BGsy; kzv(0) = BGkz
    leakv(0) = BGaquitardleak
    unconfv(0) = BGunconfined
    kv(1:CInum) = CIk
    sv(1:CInum) = CIss
    av(1:CInum) = CIk/CIss
    porv(1:CInum) = CIpor
    k2v(1:CInum) = CIaquitardK
    s2v(1:CInum) = CIaquitardSs
    a2v(1:CInum) = CIaquitardK/CIaquitardSs
    b2v(1:CInum) = CIaquitardb
    leakv(1:CInum) = CIaquitardleak
    unconfv(1:CInum) = CIunconfined
    syv(1:CInum) = CIsy
    kzv(1:CInum) = CIkz

    write(16,*) kv, '  ||    k vector'
    write(16,*) av, '  ||    alpha vector'
    write(16,*) sv, '  ||    Ss vector'
    write(16,*) porv, '  ||    porosity vector'
    write(16,*) k2v, '  ||    aquitard k vector'
    write(16,*) a2v, '  ||    aquitard alpha vector'
    write(16,*) s2v, '  ||    aquitard Ss vector'
    write(16,*) b2v, '  ||    aquitard thickness vector'
    write(16,*) leakv, '  ||    aquitard leakyness type'
    write(16,*) unconfv, '  ||    unconfined flag'
    write(16,*) syv,'  ||    Sy vector'
    write(16,*) kzv,'  ||    Kz vector'

    ! wells
    read(15,*) WLnum
    write(16,*) WLnum, '  ||    number of wells'
    allocate(WLx(WLnum), WLy(WLnum), WLr(WLnum), WLq(WLnum), &
           & WLtime(WLnum), WLstor(WLnum),WLdskin(WLnum))
    ! this is sort of out of place
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(WLstorCoeff(WLnum))
    WLstorCoeff = (0.0_DP, 0.0_DP)
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    read(15,*) WLx(1:WLnum)
    write(16,*) WLx(1:WLnum), '  ||    well x'
    read(15,*) WLy(1:WLnum)
    write(16,*) WLy(1:WLnum), '  ||    well y'
    read(15,*) WLr(1:WLnum)
    write(16,*) WLr(1:WLnum), '  ||    well radius'
    read(15,*) WLq(1:WLnum)
    write(16,*) WLq(1:WLnum), '  ||    well pumping rate'
    read(15,*) WLtime(1:WLnum)
    write(16,*) WLtime(1:WLnum), '  ||    well time function index'
    if(any(WLtime > 0)) then  !! all or nothing
       ! previous way of doing things
       allocate(WLtpar(WLnum,2))
       read(15,*) WLtpar(1:WLnum,1)
       write(16,*) WLtpar(1:WLnum,1), '  ||    well time function parameter "a"'
       read(15,*) WLtpar(1:WLnum,2)
       write(16,*) WLtpar(1:WLnum,2), '  ||    well time function parameter "b"'
    else
       ! this way is "transpose" of other way
       ! new more flexible way for n piecewise constant sections
       ! n beginning times + 1 final end time + n pumping rates = 2n+1
       allocate(WLtpar(WLnum,2*maxval(abs(WLtime))+1))
       ! (some pumping wells may have more steps than others)
       WLtpar = -huge(1.0)
       do j=1,WLnum
          read(15,*) WLtpar(j,1:2*abs(WLtime(j))+1)
          write(16,*) WLtpar(j,1:abs(WLtime(j))+1),' | ',&
               & WLtpar(j,abs(WLtime(j))+2:2*abs(WLtime(j))+1), '  ||    well ti, tf | Q for well ',j
       end do
    end if
    
       
    read(15,*) WLstor(1:WLnum)
    write(16,*) WLstor(1:WLnum), '  ||    well bore storage computed?'
    read(15,*) WLdskin(1:WLnum)
    write(16,*) WLdskin(1:WLnum), '  ||    dimensionless well bore skin?'

    read(15,*) BGcalc
    write(16,*) BGcalc, '  || re-calculate coefficients'

    ! particles
    if (BGparticle) then
       read(15,*) PARnum, PARtol, PARdt, PARmaxStep, PARstreakSkip
       allocate(PARx(PARnum), PARy(PARnum), PARti(PARnum), PARtf(PARnum), &
              & PARint(PARnum), PARInclIn(PARnum))
       read(15,*) PARx(1:PARnum)
       read(15,*) PARy(1:PARnum)
       read(15,*) PARti(1:PARnum)
       read(15,*) PARtf(1:PARnum)
       read(15,*) PARint(1:PARnum)
       read(15,*) PARInclIn(1:PARnum)

       write(16,*) PARnum, PARtol, PARdt, PARmaxStep, '  ||    number, tolerance, dt, max flux for particles '
       write(16,*) PARx(1:PARnum), '  ||    part initial x'
       write(16,*) PARy(1:PARnum), '  ||    part initial y'
       write(16,*) PARti(1:PARnum), '  ||    part initial t'
       write(16,*) PARtf(1:PARnum), '  ||    part maximum t'
       write(16,*) PARint(1:PARnum), '  ||    part integration method'
       write(16,*) PARInclIn(1:PARnum), '  ||    part begins inside CH/CF incl?'
    else
       write(16,*) '  || no particle data read'
    endif

    close(15)
    close(16)

!!$  contains
!!$    integer function error(ierr) 
!!$      integer, intent(in) :: ierr
!!$      
!!$      write(*,'(A,I2)') 'ERROR reading input file, line:',ierr
!!$      error = 1
!!$      stop
!!$      
!!$    end function error

  end subroutine readInput

  !##################################################
  subroutine readEllipseInput(filename)
    use element_specs
    use error_handler, only : fileerror
    implicit none

    character(128), intent(in) :: filename
    character(128) :: subname = 'ReadInput', echofname = 'echo_input'
    integer :: ierr

    open(UNIT=15, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    if (ierr /= 0) call fileError(filename,ierr,subname,1)
    
    open(UNIT=16, FILE=echofname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
    if (ierr /= 0) call fileError(echofname,ierr,subname,0)

    ! general options / choices
    read(15,*) BGcontour, BGoutput, BGoutFname, BGcoeffFName
    read(15,*) BGpor, BGk, BGss

    write(16,*)  BGcontour, BGoutput,  trim(BGoutFname),' ',trim(BGcoeffFName), &
         & '  ||   Lcontour, Ioutput, out/coeff fnames'
    write(16,*) BGpor, BGk, BGss, '  ||    por, k, Ss'

    read(15,*) BGnumx, BGnumy, BGnumt
    allocate(BGx(BGnumx), BGy(BGnumy), BGt(BGnumt))
    read(15,*) BGx(1:BGnumx)
    read(15,*) BGy(1:BGnumy)
    read(15,*) BGt(1:BGnumt)
    
    ! inverse Laplace transform parameters
    read(15,*) INValpha, INVtol, INVm, INVsmooth
    if (INVtol < SMALL) INVtol = SMALL 

    write(16,*) INValpha, INVtol, INVm, INVsmooth,'  ||    alpha, tol, M, smooth'
     
    ! elliptical elements
    read(15,*) EInum, EIn, EIm, EIms, EImatchTol
    allocate(EIibnd(EInum), EImatch(EInum), EIspec(EInum), EItheta(EInum), &
           & EIeta(EInum), EIx(EInum), EIy(EInum), EIf(EInum), EIk(EInum), &
           & EIss(EInum), EIpor(EInum), EIarea(EInum), &
           & EIAreaTime(EInum), EIAtpar(EInum,2), EIBdryTime(EInum), EIBtpar(EInum,2))
    read(15,*) EIibnd(1:EInum)
    read(15,*) EIspec(1:EInum)
    read(15,*) EIeta(1:EInum)
    read(15,*) EIf(1:EInum)
    read(15,*) EIx(1:EInum)
    read(15,*) EIy(1:EInum)
    read(15,*) EItheta(1:EInum)
    read(15,*) EIk(1:EInum)
    read(15,*) EIss(1:EInum)
    read(15,*) EIpor(1:EInum)
    read(15,*) EIarea(1:EInum)
    read(15,*) EIAreatime(1:EInum)
    read(15,*) EIAtpar(1:EInum,1)  
    read(15,*) EIAtpar(1:EInum,2)
    read(15,*) EIBdrytime(1:EInum)
    read(15,*) EIBtpar(1:EInum,1)  
    read(15,*) EIBtpar(1:EInum,2)

    EImatch = .false.
    where (EIibnd == -1 .or. EIibnd == 0 .or. EIibnd == +1)
       EImatch = .true.  !! active element
    end where
       
    write(16,*) EInum, EIn, EIm, EImatchTol, &
         & '  ||    number of inclusions, N, M; match_tol '
    write(16,*) EIibnd(1:EInum), '  ||    ibnd array'
    write(16,*) EImatch(1:EInum), '  ||    matching array'
    write(16,*) EIspec(1:EInum), '  ||    specified quantity'
    write(16,*) EIeta(1:EInum), '  ||    incl eta_0'
    write(16,*) EIf(1:EInum), '  ||    incl semi-focal length'
    write(16,*) EIx(1:EInum), '  ||    incl ctr x'
    write(16,*) EIy(1:EInum), '  ||    incl ctr y'
    write(16,*) EItheta(1:EInum), '  ||    incl angle w/ x+ axis'
    write(16,*) EIk(1:EInum), '  ||    incl k'
    write(16,*) EIss(1:EInum), '  ||    incl Ss'
    write(16,*) EIpor(1:EInum), '  ||    incl porosity'
    write(16,*) EIarea(1:EInum), '  ||    incl area rch rate'
    write(16,*) EIAreaTime(1:EInum), '  ||    incl Area flux time fcn index'
    write(16,*) EIAtpar(1:EInum,1), '  ||    incl Area flux time fcn parameter "a"'
    write(16,*) EIAtpar(1:EInum,2), '  ||    incl Area flux time fcn parameter "b"'
    write(16,*) EIBdryTime(1:EInum), '  ||    incl Boundary head/flux time fcn index'
    write(16,*) EIBtpar(1:EInum,1), '  ||    incl Boundary head/flux time fcn parameter "a"'
    write(16,*) EIBtpar(1:EInum,2), '  ||    incl Boundary head/flux time fcn parameter "b"'


    ! build up k,s, por and alpha vectors
    allocate(kv(0:EInum),sv(0:EInum),av(0:EInum),porv(0:EInum))
    kv(0) = BGk; sv(0) = BGss; av(0) = BGk/BGss; porv(0) = BGpor;
    kv(1:EInum) = EIk
    sv(1:EInum) = EIss
    av(1:EInum) = EIk/EIss
    porv(1:EInum) = EIpor

    write(16,*) kv, '  ||    k vector'
    write(16,*) av, '  ||    alpha vector'
    write(16,*) sv, '  ||    Ss vector'
    write(16,*) porv, '  ||    porosity vector'

    ! elliptical inclusions

    ! infinite linear disconinuities

    ! wells
    read(15,*) WLnum
    write(16,*) WLnum, '  ||    number of wells'
    allocate(WLx(WLnum), WLy(WLnum), WLr(WLnum), WLq(WLnum), &
           & WLtime(WLnum), WLtpar(WLnum,2))
    read(15,*) WLx(1:WLnum)
    write(16,*) WLx(1:WLnum), '  ||    well x'
    read(15,*) WLy(1:WLnum)
    write(16,*) WLy(1:WLnum), '  ||    well y'
    read(15,*) WLr(1:WLnum)
    write(16,*) WLr(1:WLnum), '  ||    well radius'
    read(15,*) WLq(1:WLnum)
    write(16,*) WLq(1:WLnum), '  ||    well pumping rate'
    read(15,*) WLtime(1:WLnum)
    write(16,*) WLtime(1:WLnum), '  ||    well time function index'
    read(15,*) WLtpar(1:WLnum,1)
    write(16,*) WLtpar(1:WLnum,1), '  ||    well time function parameter "a"'
    read(15,*) WLtpar(1:WLnum,2)
    write(16,*) WLtpar(1:WLnum,2), '  ||    well time function parameter "b"'

    read(15,*) BGcalc
    write(16,*) BGcalc, '  || re-calculate coefficients'

    close(15)
    close(16)

  end subroutine readEllipseInput

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

!!$  !##################################################
!!$  function mergeSortLocations(v,v1) result(order)
!!$    use constants, only : DP
!!$
!!$    real(DP), intent(in), dimension(:) :: v, v1
!!$    integer, dimension(size(v)+size(v1)) :: order
!!$    integer :: old, neu, cnt, less
!!$    real(DP) :: difo, difn
!!$
!!$    ! three things to keep track of while working through the locations
!!$    ! index in old output => old -
!!$    ! index in new output => neu +
!!$    ! index in order vector => cnt
!!$
!!$    order = 0
!!$    old = 1; neu = 1; cnt = 1
!!$    less = 0
!!$
!!$    ! find first element
!!$    if (v(neu) < v1(old)) then ! new v is farther left/down
!!$       order(cnt) = neu
!!$       neu = neu + 1
!!$    elseif (v1(old)< v(neu)) then ! v1 is farther left/down
!!$       order(cnt) = -old
!!$       old = old + 1
!!$    else ! they are equal (advance both counters)
!!$       order(cnt) = neu
!!$       old = old + 1
!!$       neu = neu + 1
!!$       less = less + 1
!!$    end if
!!$
!!$    do cnt = 2, size(v)+size(v1)
!!$       difo = v1(old) - order(cnt)
!!$       difn = v(neu) - order(cnt)
!!$       if (difn < difo) then
!!$          order(cnt) = neu
!!$          neu = neu + 1
!!$       elseif (difo < difn) then
!!$          order(cnt) = -old
!!$          old = old + 1
!!$       else
!!$          order(cnt) = neu
!!$          old = old + 1
!!$          neu = neu + 1
!!$          less = less + 1
!!$       end if
!!$    end do
!!$
!!$    ! zero out values not used due to duplication
!!$    order(size(order)-less:) = 0
!!$
!!$  end function mergeSortLocations

end module file_ops
