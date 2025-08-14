!
! Copyright (c) 2011-2025 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

! this module contains most of the I/O associated with LT-AEM

module file_ops
  implicit none

  private
  public :: readInput, writeResults, writeGeometry, read_coeff, dump_coeff

contains

  !##################################################
  ! this routine read the main input file, and allocates the main
  ! data structures used to store data.
  subroutine readInput(s,dom,bg,c,e,p)

    use constants, only : DP, lenFN, PI
    use type_definitions, only : solution, particle, domain, element, circle, ellipse, explain_type
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit, stderr => error_unit

    type(solution), intent(inout) :: s
    type(particle), intent(out), allocatable :: p(:)
    type(domain),   intent(out) :: dom
    type(element),  intent(out) :: bg
    type(circle),   intent(out), allocatable :: c(:)
    type(ellipse),  intent(out), allocatable :: e(:)

    type(explain_type) :: explain
    character(4) :: chint
    character(20), dimension(3) :: fmt
    character(1024) :: buf 
    character(lenFN+5) :: echofname
    character(lenFN) :: circleFname, ellipseFname, particleFname
    integer :: ierr, j, ntot, nC, nE, ln, sln,s1,s2,slen, idx
    real(DP) :: tmp
    character(55) :: explain_txt

    ! unit numbers for input/output files
    integer, parameter :: UINPUT = 15, UECHO = 16, UCIRC = 22, UELIP = 33, UPAR = 44 

    ln = 0
    
    open(unit=UINPUT, file=s%infname, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(stderr,'(2A)') 'ERROR READINPUT: error opening input file ',trim(s%infname)
       stop 100
    end if

    ! search for '.' in filename from the end of the string
    idx = index(s%infName,'.',back=.true.)  
    if (idx > 0) then
       s%echofname = s%infname(1:idx-1)//'.echo'
       s%qfname = s%infname(1:idx-1)//'.Q'
    else
       s%echofname = trim(s%infname)//'.echo'
       s%qfname = trim(s%infname)//'.Q'
    end if

    open(unit=UECHO, file=s%echofname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       write(stderr,'(2A)') 'ERROR READINPUT: error opening echo file ',echofname
       stop 101
    else
       ! add a file variable at top of file to set Emacs to auto-revert mode
       write(UECHO,'(A)') '-*-auto-revert-*-'
    end if

    ! solution-specific and background aquifer parameters
    read(UINPUT,*,iostat=ierr) s%calc, s%particle, s%contour, s%deriv, s%Qcalc, s%output; ln=ln+1
    if (ierr /= 0) then
      write(stderr,*) 'ERROR: line ',ln,' (flags: calc, particle, contour, deriv, Qcalc, output) of ',s%infname
      write(stderr,*) 'WARNING: may have specified a secondary (i.e., circle, ellipse, or particle) '//&
           &'input file as main input file on command line'
      stop 110
    end if

    backspace(UINPUT) ! optional input parameter (read first input line again)
    read(UINPUT,*,iostat=ierr) s%calc, s%particle, s%contour, s%deriv, s%Qcalc, s%output, s%debug
    if (ierr /= 0) then
      s%debug = .false.
    end if
    
    if ((.not. s%particle) .and. (.not. s%contour)) then
       s%timeseries = .true.
    else
       s%timeseries = .false.
    end if

    ! if s%output > 100, then don't dump matching results to
    ! file for restart (a minor speedup?)
    if (s%output > 100) then
       s%skipdump = .true.
       s%output = s%output - 100
    else
       s%skipdump = .false.
    end if

    ! some simple debugging of problem-type / output-type combinations
    if (.not.(s%output == 1 .or. s%output == 2 .or. &
         & s%output == 10 .or. s%output == 11 .or. s%output == 12 .or. &
         & s%output == 20 .or. s%output == 21)) then
       write(stderr,*) 'ERROR: values (line ',ln,') s%output must be in '&
            &//'{1,2,10,11,12,20,21}: ',s%output
       stop 200
    end if
    if ((s%output >= 10  .and. s%contour) .or. (s%output < 10 .and. .not. s%contour)) then
       write(stderr,*) 'ERROR: values (line ',ln,') s%output should be in {1,2} '&
            &//'when contour output is selected: ',s%output,s%contour
       stop 2010
    end if

    if (((s%output < 10 .or. s%output >= 20) .and. s%timeseries) .or. &
         & (s%output >= 10 .and. s%output < 20 .and. .not. s%timeseries)) then
       write(stderr,*) 'ERROR: values (line ',ln,') s%output should be in {10,11,12}'&
            &//' when timeseries output is selected: ',s%output,s%timeseries
       stop 2020
    end if

    if ((s%output < 20 .and. s%particle) .or. (s%output >= 20  .and. .not. s%particle)) then
       write(stderr,*) 'ERROR: values (line ',ln,') iff s%particle==.True., '&
            &//'s%output should be in {20,21}: ',s%output,s%particle
       stop 2030
    end if    

    if (s%Qcalc .and. s%particle) then
       write(stdout,*) 'WARNING: resetting Qcalc to false for particle tracking'
       s%Qcalc = .false.
    end if
    
    ! read output filename
    read(UINPUT,*,iostat=ierr) s%outFname; ln=ln+1
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line ',ln,' (output filename) of ',s%infname
       stop 2032
    end if

    ! other names just based on output filename
    idx = index(s%outFname,'.') ! is there a . in filename?
    if (idx > 0) then
       s%elemHfName = s%outFname(1:idx-1)//'.elem'
       s%geomfName  = s%outFname(1:idx-1)//'.geom'
    else       
       s%elemHfName = trim(s%outFname)//'.elem' 
       s%geomfName  = trim(s%outFname)//'.geom'
    end if
    
    read(UINPUT,*,iostat=ierr) bg%por, bg%k, bg%ss; ln=ln+1
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',ln,'(basic background properties: por, k, ss) of',s%infname
       stop 204
    end if
    backspace(UINPUT)
    read(UINPUT,*,iostat=ierr) bg%por, bg%k, bg%ss, bg%wave
    if (ierr /= 0) then
      ! optional
      bg%wave = .false.
    end if
    
    if (any([bg%por,bg%k,bg%ss] <= 0.0)) then
       write(stderr,*) 'ERROR: value (line ',ln,') bg%por, bg%k, bg%ss '&
            & // 'must all be > 0.0: ',[bg%por,bg%k,bg%ss]
       stop 205
    end if
    
    read(UINPUT,*,iostat=ierr) bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb; ln=ln+1
    if (ierr /= 0 .and. bg%leakFlag /= 0) then
       write(stderr,*) 'ERROR: line ',ln,' (leaky aquitard props:' &
            & // 'leakFlag, aquitardK, aquitardSs, aquitardb) of ',s%infname
       stop 2050
    end if
    if (bg%leakFlag < 0 .or. bg%leakFlag > 3) then
       write(stderr,*) 'ERROR: value (line ',ln,') leak flag (bg%leakFlag) must be in {0,1,2,3}'
       stop 2051
    end if
    if (any([bg%aquitardK,bg%aquitardSs,bg%aquitardb] <= 0.0) .and. bg%leakFlag > 0) then
       write(stderr,*) 'ERROR: value (line ',ln,') bg%aquitardK, bg%aquitardSs, bg%aquitardb' &
            & // 'must all be > 0.0: ',[bg%aquitardK,bg%aquitardSs,bg%aquitardb]
       stop 206
    end if

    read(UINPUT,*,iostat=ierr) bg%unconfinedFlag, bg%Sy, bg%kz, bg%b; ln=ln+1
    if (ierr /= 0 .and. bg%unconfinedFlag) then
       write(stderr,*) 'ERROR: line',ln,'(unconfined props: unconfinedFlag, Sy, Kz, b) of',s%infname
       stop 2060
    end if

    if (any([bg%kz,bg%b] <= 0.0) .and. bg%unconfinedFlag) then
       write(stderr,*) 'ERROR: value (line',ln,&
            & ') vertical K (bg%kz), thickness (bg%b) must be > 0.0: ', &
            & [bg%kz,bg%b]
       stop 207
    end if

    ! Sy can be zero
    if (bg%Sy < 0.0 .and. bg%unconfinedFlag) then
       write(stderr,*) 'ERROR: value (line',ln,&
            &') specific yield (bg%Sy) must be non-negative:',bg%Sy
       stop 2071
    end if

    read(UINPUT,*,iostat=ierr) bg%dualPorosityFlag, bg%matrixSs, &
         & bg%lambda, bg%multiporosityDiffusion, bg%kappa, bg%NDiffTerms; ln=ln+1
    if (ierr /= 0 .and. bg%dualPorosityFlag) then
       write(stderr,*) 'ERROR: line',ln,'(dual porosity: dualPorosityFlag, matrixSs, ' &
            & // 'lambda, multiporosityDiffusion, kappa, NDiffTerms) of',s%infname
       stop 2072
    end if
    
    if (any([bg%matrixSs,bg%lambda,bg%kappa] <= 0.0) .and. bg%dualPorosityFlag) then
       write(stderr,*) 'ERROR: value (line',ln,') matrix specific storage (bg%matrixSs), ' &
            & // 'matrix/fracture connection factor (bg%lambda), and matrix/fracture k ratio (bg%kappa) ' &
            & // ' must be > 0.0:',&
            & [bg%matrixSs,bg%lambda,bg%kappa]
       stop 20721
    end if

    if (bg%NDiffTerms < 0) then
       write(stderr,*) 'ERROR: value (line',ln,&
            &') number of terms for multiporosity diffusion must be >=0:',bg%nDiffTerms
       stop 20722
    end if

    if (bg%multiporosityDiffusion < 0 .or. bg%multiporosityDiffusion > 3) then
       write(stderr,*) 'ERROR: value (line',ln,') matrix multiporosity diffusion order '//&
            &'must be 0 (off) or {1,2,3}:' ,bg%multiporosityDiffusion
       stop 20722
    end if
    
    ! echo input to file
    write(UECHO,'(A)') '=============== OVERALL OPTIONS ==============='
    write(UECHO,'(L1,1X,A)') s%calc, '  ||  re-calculate coefficients?'
    write(UECHO,'(L1,1X,A)') s%particle, '  ||  particle tracking?'
    write(UECHO,'(L1,1X,A)') s%contour, '  ||  compute solution for contour map? (many locations, few times)'
    write(UECHO,'(L1,1X,A)') s%timeseries, '  ||  compute solution for timeseries? (few locations, many times)'
    write(UECHO,'(L1,1X,A)') s%deriv, '  ||  compute log(t) derivative of heads?'
    write(UECHO,'(L1,1X,A)') s%qcalc, '  ||  compute element boundary flux? (saved to .Q file)'
    write(UECHO,'(I0,1X,A)') s%output,trim(explain%output(s%OEMap(s%output)))//&
         & '  ||    output flag (1:2 gnuplot/matlab contours, '//&
         & '10:11 gnuplot time series w/wo vel, 12: matlab time series, '//&
         &'20:21 pathline/streakline gnuplot)'
    write(UECHO,'(2A)') trim(s%outFname),'  ||    output file name'
    write(UECHO,'(L1,1X,A)') s%debug, '  || debugging output'
    write(UECHO,'(A)') '=============== BACKGROUND PROPERTIES ==============='
    write(UECHO,'(3(ES12.5,1X),L1,A)') bg%por, bg%k, bg%ss, bg%wave, '  ||   AQUIFER properties : '//&
         &'porosity, hydraulic conductivity, specific storage, wave eqn?'
    write(UECHO,'(I0,1X,A,1X,3(ES12.5,1X),A)') bg%leakFlag, trim(explain%leakFlag(bg%leakFlag)), &
         & bg%aquitardK, bg%aquitardSs, bg%aquitardb, &
         & '  ||   adjacent AQUITARD properties : leaky flag, hydraulic condictivity, specific storage, thickness'

    if (bg%unconfinedFlag) then
       explain_txt = '(unconfined aquifer ON)'
    else
       explain_txt = '(unconfined aquifer OFF)'
    end if
    write(UECHO,'(L1,1X,A,1X,3(ES12.5,1X),A)') bg%unconfinedFlag, trim(explain_txt), bg%Sy, bg%kz, bg%b, &
         & '  ||   UNCONFINED aquifer properties : unconfined?, specific yield, vertical hydraulic conductivity, thickness'

    if (bg%dualPorosityFlag) then
       explain_txt = '(dual porosity aquifer ON)'
    else
       explain_txt = '(dual porosity aquifer OFF)'
    end if
    write(UECHO,'(L1,1X,A,1X,2(ES12.5,1X),I0,1X,ES12.5,1X,I0,A)') bg%dualPorosityFlag, trim(explain_txt), &
         & bg%matrixSs, bg%lambda, bg%multiporosityDiffusion, bg%kappa, bg%NDiffTerms, &
         & '  ||   DUAL POROSITY aquifer properties : dual porosity?, matrix specific storage, matrix/fracture lambda, '//&
         & 'multiporosity diffusion index, matrix/fracture K ratio, number terms in diffusion series'

    ! desired solution points/times
    read(UINPUT,*,iostat=ierr) s%nx, s%ny, s%nt; ln=ln+1
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',ln, '(number solution locations: nx,ny,nt) of',s%infname
       stop 2073
    end if
    
    if (any([s%nx,s%ny] < 0) .or. s%nt < 1) then
       write(stderr,*) 'ERROR: value (line',ln,') number x,y observation locations '&
            &//'(s%nx,s%ny) must be >0  and number times must be >1', [s%nx,s%ny,s%nt]
       stop 208
    end if
    if (any([s%nx,s%ny] < 1) .and. .not. s%Qcalc) then
       write(stderr,*) 'ERROR: no calculation locations, and Qcalc is false, '//&
            &'therefore no output would be generated, stopping calculation.'
       stop 2081
    end if
    if (s%timeseries .and. s%nx /= s%ny) then
       write(stdout,*) 'WARNING: for time series output nx==ny.  nx=',s%nx,' ny=',s%ny
       write(stdout,*) '**** resetting ny to nx ****'
       s%ny = s%nx
    end if
    allocate(s%x(s%nx), s%y(s%ny), s%t(s%nt))
    if (s%timeseries) then
       allocate(s%obsname(s%nx))
       read(UINPUT,'(512A)',iostat=ierr) buf; ln=ln+1
       if (ierr /= 0) then
          write(stderr,*) 'ERROR: line ',ln,' (location names) of',s%infname
          stop 2081
       end if
       
       s1 = 1
       slen = len_trim(buf) ! don't include trailing blanks
       do j = 1,s%nx
          ! location names are separated by "|" character
          s2 = index(buf(s1:),'|')
          if (s2 == 0) then
             ! no "|" separator character found
             if (j == 1 .and. slen > 1) then
                ! first name no separator 
                s%obsname(1) = trim(buf)
                s1 = slen
             elseif(j == s%nx .and. s1 > 1) then
                ! last name no separator (normal)
                s%obsname(s%nx) = trim(buf(s1:))
             else
                write(chint,'(I4.4)') j
                s%obsname(j) = 'LOC-'//chint ! generic name
             end if
          else
             s%obsname(j) = buf(s1:s1+s2-2)
          end if
          s1 = s1+s2
       end do
    else
       allocate(s%obsname(0))
       read(UINPUT,*); ln=ln+1 ! read a single placeholder name line anyway
    end if

    read(UINPUT,*,iostat=ierr) tmp; ln=ln+1
    if (ierr == 0)  then
       ! read one real successfully
       backspace(UINPUT)
       read(UINPUT,*,iostat=ierr) s%x(:)
       if (ierr /= 0) then
          ! failed reading vector
          write(stderr,*) 'ERROR: line',ln,'x calc locations (s%x) of',s%infname
       end if
    else
       ! compute vector if first chars are 'linvec' or 'logvec'
       s%x(:) = computeVector(UINPUT,s%nx,ln)
    end if
    
    read(UINPUT,*,iostat=ierr) tmp; ln=ln+1
    if (ierr == 0)  then
       backspace(UINPUT)
       read(UINPUT,*,iostat=ierr) s%y(:)
       if (ierr /= 0) then
          write(stderr,*) 'ERROR: line',ln,'y calc locations (s%y) of',s%infname
       end if
    else
       s%y(:) = computeVector(UINPUT,s%ny,ln)
    end if

    ! shift x & y values to origin (useful when x and y are UTM coordinates)
    ! NOTE: shift is computed based on range of calculation locations,
    ! but there could also be a circle/ellipse element that is very big or very far away.
    s%xshift = (maxval(s%x) + minval(s%x))/2.0
    s%yshift = (maxval(s%y) + minval(s%y))/2.0
    s%x(:) = s%x(:) - s%xshift
    s%y(:) = s%y(:) - s%yshift

    read(UINPUT,*,iostat=ierr) tmp; ln=ln+1
    if (ierr == 0)  then
       backspace(UINPUT)
       read(UINPUT,*,iostat=ierr) s%t(:)
       if (ierr /= 0) then
          write(stderr,*) 'ERROR: line',ln,'t calc times (s%t) of',s%infname
       end if
    else
       s%t(:) = computeVector(UINPUT,s%nt,ln)
    end if

    if (any(s%t <= 0.0)) then
       write(stderr,*) 'ERROR: value (line ',ln,') requires t>0',s%t
       stop 2085
    end if

    if (s%nt >= 2) then
      if (any((s%t(2:s%nt) - s%t(1:s%nt-1)) <= 0.0)) then
         write(stderr,*) 'ERROR: value (line ',ln,') requires monotonically increasing times',s%t
         stop 2086
      end if
    end if
    
    if (.not. s%particle) then
       fmt(1) = '(    (ES12.5,1X),A) '
       write(UECHO,'(A)') '=============== CALCULATION LOCATIONS/TIMES ==============='
       write(UECHO,'(3(I0,1X),A)') s%nx, s%ny, s%nt, '  ||    numX, numY, numt'
       if (s%nx > 0 .and. s%ny > 0) then
          write(fmt(1)(2:5),'(I4.4)') s%nx
          write(UECHO,'(ES21.14,A)') s%xshift,'  ||    xshift'
          write(UECHO,fmt(1)) s%x(:)+s%xshift,'  ||    original x Vector'
          write(UECHO,fmt(1)) s%x(:),         '  ||    shifted x Vector'
          write(fmt(1)(2:5),'(I4.4)') s%ny
          write(UECHO,'(ES21.14,A)') s%yshift,'  ||    yshift'
          write(UECHO,fmt(1)) s%y(:)+s%yshift, '  ||    original y Vector'
          write(UECHO,fmt(1)) s%y(:),          '  ||    shifted y Vector'
       else
          write(UECHO,'(A)') '**WARNING: zero-length x or y output vector specified, only element flowrates computed **'
       end if
       write(fmt(1)(2:5),'(I4.4)') s%nt
       write(UECHO,fmt(1)) s%t(:), '  ||    t Vector'
    end if

    ! deHoog et al. inverse Laplace transform parameters
    read(UINPUT,*,iostat=ierr) s%alpha, s%tol, s%m; ln=ln+1
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',ln,'(deHoog: alpha, tol, m) of',s%infname
       stop 209
    end if

    if (s%M < 1) then
       write(stderr,*) 'ERROR: value (line ',ln,') s%M > 0 (typically >= 10)',s%M
       stop 2090
    end if
    if (s%tol < epsilon(s%tol)) then ! epsilon(1.0D+0) ~ 1.0E-16
       s%tol = epsilon(s%tol)
       write(stdout,'(A,ES12.5)') 'WARNING: increased deHoog INVLAP solution tolerance to ',s%tol
    end if
    if (s%alpha <= 0.0) then
       write(stdout,'(A,ES12.5)') 'WARNING: deHoog alpha typically > 0 ',s%alpha
    end if
    write(UECHO,'(A)') '=============== INVERSE LAPLACE TRANSFORM PARAMETERS ==============='
    write(UECHO,'(2(ES12.5,1X),I0,A)') s%alpha, s%tol, s%m,'  ||    deHoog: alpha, tol, M'

    ! circular (includes wells)
    read(UINPUT,*,iostat=ierr) dom%num(1),circleFname; ln=ln+1
    if (ierr /= 0) then
      write(stderr,*) 'ERROR: line ',ln,' (circles: number, circle input file name) of',s%infname
      stop 2091
    end if
    
    nc = dom%num(1)
    if (nc > 0) then
       sln = 1
       open(unit=UCIRC, file=trim(circleFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(2A)') 'ERROR: cannot open circular input file for reading ',trim(circleFname)
          stop 210
       else
          write(UECHO,'(A)') trim(circleFname)//' opened for circular input data'
       end if

       allocate(c(nc))
       read(UCIRC,*,iostat=ierr) c(1:nc)%N; sln=sln+1
       if (ierr /= 0) c(1:nc)%N = read_int(UCIRC,sln,'circle  ')
       if (any(c%N < 1)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) # Fourier terms (c%N) '//&
               &'must not be < 1 ',c%N
          stop 211
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%M; sln=sln+1
       if (ierr /= 0) c(1:nc)%M = read_int(UCIRC,sln,'circle  ')
       if (any(c%M < 1)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) # matching locations (c%M) '//&
               &'must not be < 1 ',c%M
          stop 212
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%ibnd; sln=sln+1
       if (ierr /= 0) c(1:nc)%ibnd = read_int(UCIRC,sln,'circle  ')
       if (any(c%ibnd < -1 .or. c%ibnd > 2)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) boundary type (c%ibnd) '//&
               &'must be in {-1,0,1,2} ',c%ibnd
          stop 213
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%CalcIn; sln=sln+1
       if (ierr /= 0) c(1:nc)%CalcIn = read_logical(UCIRC,sln,'circle  ')
       
       read(UCIRC,*,iostat=ierr) c(1:nc)%StorIn; sln=sln+1 
       if (ierr /= 0) c(1:nc)%StorIn = read_logical(UCIRC,sln,'circle  ')

       read(UCIRC,*,iostat=ierr) c(1:nc)%wave; sln=sln+1
       if (ierr /= 0) then
           backspace(UCIRC)
           read(UCIRC,*,iostat=ierr) c(1)%wave
           if (ierr /= 0) then
             ! optional, assume false
             c(1:nc)%wave = .false.
             backspace(UCIRC)
           else
             c(1:nc)%wave = c(1)%wave
           end if
       end if
       
       read(UCIRC,*,iostat=ierr) c(1:nc)%r; sln=sln+1
       if (ierr /= 0) c(1:nc)%r = read_real(UCIRC,sln,'circle  ')
       if (any(c%r <= 0.0)) then
          write(stderr,*) 'ERROR: value (line',sln,' cirlce input) circle radius (c%r) '//&
               & 'must be > 0.0 ',c%r
          stop 214
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%x; sln=sln+1
       if (ierr /= 0) c(1:nc)%x = read_real(UCIRC,sln,'circle  ')

       read(UCIRC,*,iostat=ierr) c(1:nc)%y; sln=sln+1
       if (ierr /= 0) c(1:nc)%y = read_real(UCIRC,sln,'circle  ')

       ! NOTE: shift is computed from range of calc locations
       c(1:nc)%x = c%x - s%xshift
       c(1:nc)%y = c%y - s%yshift
       c(1:nc)%z = cmplx(c%x, c%y, DP)

       read(UCIRC,*,iostat=ierr) c(1:nc)%k; sln=sln+1
       if (ierr /= 0) c(1:nc)%k = read_real(UCIRC,sln,'circle  ')
       if (any(c%k <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' circle input) hydraulic conductivity '//&
               &'(c%K) must be > 0.0 ',c%k
          stop 215
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%Ss; sln=sln+1
       if (ierr /= 0) c(1:nc)%Ss = read_real(UCIRC,sln,'circle  ')
       if (any(c%ss <= 0.0)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) specific storage (c%Ss) '//&
               &'must be > 0.0 ',c%Ss
          stop 216
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%por; sln=sln+1
       if (ierr /= 0) c(1:nc)%por = read_real(UCIRC,sln,'circle  ')
       if (any(c%por <= 0.0)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) porosity (c%por) '//&
               &'must be > 0.0 ',c%por
          stop 217
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%leakFlag; sln=sln+1
       if (ierr /= 0) c(1:nc)%leakFlag = read_int(UCIRC,sln,'circle  ')
       if (any(c%leakFlag < 0) .or. any(c%leakFlag > 3)) then
          write(stderr,*) 'ERROR: value (line',sln,' leaky flag) (c%leakFlag) circle '//&
               &'input; input must be in {0,1,2,3}'
          stop 2170
       end if       

       read(UCIRC,*,iostat=ierr) c(1:nc)%aquitardK; sln=sln+1
       if (ierr /= 0) c(1:nc)%aquitardK = read_real(UCIRC,sln,'circle  ')
       if (any(c%aquitardK <= 0.0 .and. c%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) aquitard vertical K '//&
               &'(c%aquitardK) must be > 0.0 ',c%aquitardk
          stop 218
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%aquitardSs; sln=sln+1
       if (ierr /= 0) c(1:nc)%aquitardSs = read_real(UCIRC,sln,'circle  ')
       if (any(c%aquitardSS <= 0.0 .and. c%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) aquitard specific storage '//&
               &'(c%aquitardSs) must be > 0.0 ',c%aquitardSs
          stop 219
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%aquitardb; sln=sln+1
       if (ierr /= 0) c(1:nc)%aquitardb = read_real(UCIRC,sln,'circle  ')
       if (any(c%aquitardB <= 0.0 .and. c%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' circle input) aquitard thickness '//&
               &'(c%aquitardB) must be > 0.0 ',c%aquitardB
          stop 220
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%unconfinedFlag; sln=sln+1
       if (ierr /= 0) c(1:nc)%unconfinedFlag = read_logical(UCIRC,sln,'circle  ')

       read(UCIRC,*,iostat=ierr) c(1:nc)%Sy; sln=sln+1
       if (ierr /= 0) c(1:nc)%Sy = read_real(UCIRC,sln,'circle  ')
       if (any(c%sy <= 0.0 .and. c%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) specific yield (c%Sy) '//&
               &'must be > 0.0 ',c%sy
          stop 221
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%Kz; sln=sln+1
       if (ierr /= 0) c(1:nc)%Kz = read_real(UCIRC,sln,'circle  ')
       if (any(c%kz <= 0.0 .and. c%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) aquifer vertical K '//&
               &'(c%Kz) must be > 0.0 ',c%kz
          stop 222
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%b; sln=sln+1
       if (ierr /= 0) c(1:nc)%b = read_real(UCIRC,sln,'circle  ')
       if (any(c%b <= 0.0 .and. c%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) aquifer thickness '//&
               &'(c%B) must be > 0.0 ',c%b
          stop 223
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%dualPorosityFlag; sln=sln+1
       if (ierr /= 0) c(1:nc)%dualPorosityFlag = read_logical(UCIRC,sln,'circle  ')

       read(UCIRC,*,iostat=ierr) c(1:nc)%matrixSs; sln=sln+1
       if (ierr /= 0) c(1:nc)%matrixSs = read_real(UCIRC,sln,'circle  ')
       if (any(c%matrixSs <= 0.0 .and. c%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) matrix specific storage '//&
               &'(c%matrixSs) must be > 0.0', c%matrixSs
          stop 2230
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%lambda; sln=sln+1
       if (ierr /= 0) c(1:nc)%lambda = read_real(UCIRC,sln,'circle  ')
       if (any(c%lambda < 0.0 .and. c%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) '//&
               & 'matrix/fracture connection (c%lambda) must be non-negative', c%lambda
          stop 2231
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%multiporosityDiffusion; sln=sln+1
       if (ierr /= 0) c(1:nc)%multiporosityDiffusion = read_int(UCIRC,sln,'circle  ')
       if (any((c%multiporosityDiffusion < 0 .or. c%multiporosityDiffusion > 3) &
            & .and. c%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) '//&
               & 'multiporosity diffusion index must be {0,1,2,3}', c%multiporosityDiffusion
          stop 2232
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%kappa; sln=sln+1
       if (ierr /= 0) c(1:nc)%kappa = read_real(UCIRC,sln,'circle  ')
       if (any(c%kappa < 0.0 .and. c%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) '//&
               & 'matrix/fracture hydraulic conductivity ratio (c%kappa) must be non-negative', c%kappa
          stop 2233
       end if
       
       read(UCIRC,*,iostat=ierr) c(1:nc)%NDiffTerms; sln=sln+1
       if (ierr /= 0) c(1:nc)%NDiffTerms = read_int(UCIRC,sln,'circle  ')
       if (any(c%NDiffTerms < 0 .and. c%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) '//&
               & 'number terms in leaky diffusion series (c%nDiffTerms) must be >=0', c%NDiffTerms
          stop 2234
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%dskin; sln=sln+1
       if (ierr /= 0) c(1:nc)%dskin = read_real(UCIRC,sln,'circle  ')
       if (any(c%dskin <= 0.0)) then
          write(stderr,*) 'ERROR: value (line',sln,' circle input) dimensionless skin '//&
               &'parameter (c%Dskin) must be > 0.0 ',c%dskin
          stop 224
       end if

       read(UCIRC,*,iostat=ierr) c(1:nc)%areaQ; sln=sln+1 
       if (ierr /= 0) c(1:nc)%areaQ = read_real(UCIRC,sln,'circle  ')
       ! source term can be positive, negative or zero (no need for check)
       
       read(UCIRC,*,iostat=ierr) c(1:nc)%bdryQ; sln=sln+1 
       if (ierr /= 0) c(1:nc)%bdryQ = read_real(UCIRC,sln,'circle  ')
       ! source term can be positive, negative or zero (no need for check)

       where (c(:)%ibnd == -1 .or. c(:)%ibnd == 0 .or. c(:)%ibnd == +1)
          c(:)%match = .true.
       elsewhere !( currently just 2)
          c(:)%match = .false.
       end where

       write(chint,'(I4.4)') dom%num(1)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(UECHO,'(A)') '=============== CIRCULAR ELEMENT PROPERTIES ==============='
       write(UECHO,'(I0,A)') dom%num(1),        '  ||   # circular elements (including wells)'
       write(UECHO,fmt(1)) c(:)%n,              '  ||   # circular free parameter (Fourier coefficients)'
       write(UECHO,fmt(1)) c(:)%m,              '  ||   # circular matching locations'
       write(UECHO,fmt(1)) c(:)%ibnd,           '  ||   circle ibnd array '//&
            &'(-1: spec head, 0: matching, :1 spec total flux, 2: spec elem flux)'
       write(UECHO,fmt(2)) c(:)%match,          '  ||   circle matching array'
       write(UECHO,fmt(2)) c(:)%calcin,         '  ||   calculate inside this circle?'
       write(UECHO,fmt(2)) c(:)%storin,         '  ||   calculate circle free-water storage effects?'
       write(UECHO,fmt(2)) c(:)%wave,           '  ||   calculate wave equation inside element?'
       write(UECHO,fmt(3)) c(:)%r,              '  ||   circle radius'
       write(UECHO,fmt(3)) c(:)%x+s%xshift,     '  ||   original circle center x'
       write(UECHO,fmt(3)) c(:)%x,              '  ||   shifted circle center x'
       write(UECHO,fmt(3)) c(:)%y+s%yshift,     '  ||   original circle center y'
       write(UECHO,fmt(3)) c(:)%y,              '  ||   shifted circle center y'
       write(UECHO,fmt(3)) c(:)%k,              '  ||   circle aquifer k'
       write(UECHO,fmt(3)) c(:)%ss,             '  ||   circle aquifer Ss'
       write(UECHO,fmt(3)) c(:)%por,            '  ||   circle aquifer porosity'
       write(UECHO,fmt(1)) c(:)%leakFlag,       '  ||   circle leaky type'
       write(UECHO,fmt(3)) c(:)%aquitardK,      '  ||   circle leaky aquitard K'
       write(UECHO,fmt(3)) c(:)%aquitardSs,     '  ||   circle leaky aquitard Ss'
       write(UECHO,fmt(3)) c(:)%aquitardb,      '  ||   circle leaky aquitard thickness'
       write(UECHO,fmt(2)) c(:)%unconfinedFlag, '  ||   circle unconfined flag'
       write(UECHO,fmt(3)) c(:)%Sy,             '  ||   circle aquifer specific yield'
       write(UECHO,fmt(3)) c(:)%Kz,             '  ||   circle aquifer vertical K'
       write(UECHO,fmt(3)) c(:)%b,              '  ||   circle aquifer thickness'
       write(UECHO,fmt(2)) c(:)%dualPorosityFlag,'  ||   circle dual porosity flag'
       write(UECHO,fmt(3)) c(:)%matrixSs,        '  ||   circle matrix Ss'
       write(UECHO,fmt(3)) c(:)%lambda,          '  ||   circle matrix/fracture connection lambda'
       write(UECHO,fmt(1)) c(:)%multiporosityDiffusion,  '  ||   circle multiporosity diffusion index'
       write(UECHO,fmt(3)) c(:)%kappa,          '  ||   circle matrix/fracture k ratio'
       write(UECHO,fmt(1)) c(:)%NDiffTerms,          '  ||   circle matrix diffusion series number terms'
       write(UECHO,fmt(3)) c(:)%dskin,           '  ||   circle boundary dimensionless skin factor'
       write(UECHO,fmt(3)) c(:)%areaQ,          '  ||   circle area rch rate'
       write(UECHO,fmt(3)) c(:)%bdryQ,          '  ||   circle boundary rch rate or head'
       
       ! area source behavior
       do j = 1,size(c,dim=1)
          read(UCIRC,*,iostat=ierr) c(j)%AreaTime; sln=sln+1
          if (ierr /= 0) then
             write(stderr,*) 'ERROR: line ',sln,' area time behavior (c%AreaTime) '//&
                  &'circle ',j,' input'
             stop 2201
          end if
          call read_time_behaviors(UCIRC,UECHO,c(j)%time,j,sln,'circle ',.true.)
       end do

       ! boundary source behavior
       do j = 1,size(c,dim=1)
          read(UCIRC,*,iostat=ierr) c(j)%BdryTime; sln=sln+1
          if (ierr /= 0) then
             write(stderr,*) 'ERROR: line ',sln,' boundary time behavior '//&
                  &'(c%BdryTime) circle ',j,' input'
             stop 2204
          end if
          call read_time_behaviors(UCIRC,UECHO,c(j)%time,j,sln,'circle ',.false.)
       end do

       close(UCIRC) ! circle input file

       if(any(c%match .and. c%storin)) then
          write(stdout,*) 'WARNING: wellbore (free-water) storage only implemented for ibnd==2'
          write(stdout,*) '**** resetting to false ****'
          where(c%match .and. c%storin)
             c%storin = .false.
          end where
       end if

       ! minor checking / correcting
       do j = 1,size(c,dim=1)
          if (c(j)%ibnd == 2 .and. c(j)%N /= 1) then
             write(stdout,'(A,I0,A)') 'WARNING wells (ibnd==2) must have N=1'
             write(stdout,*) '**** fixing circle #', j,' to N=1 ****'
             c(j)%N = 1
          end if
       end do

    else
       ! no circular elements
       allocate(c(0))
    end if

    ! elliptical (includes line sources/sinks)
    read(UINPUT,*,iostat=ierr) dom%num(2), bg%ms, ellipseFname; ln=ln+1
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line ',ln,' (ellipses: number, matrix size, ellipse file name) of',s%infname
       stop 2092
    end if

    ne = dom%num(2)
    sln = 1
    if (ne > 0) then
       open(unit=UELIP, file=trim(ellipseFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(2A)') 'ERROR: cannot open elliptical input file for reading ',&
               & trim(ellipseFname)
          stop 225
       else
          write(UECHO,'(A)') trim(ellipseFname)//' opened for elliptical input data'
       end if

       allocate(e(ne))
       read(UELIP,*,iostat=ierr) e(1:ne)%N; sln=sln+1
       if (ierr /= 0) e(1:ne)%N = read_int(UELIP,sln,'ellipse ')
       if (any(e%n < 1)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) # Fourier terms (e%N) '//&
                     &'must not be < 1 ',e%N
          stop 226
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%M; sln=sln+1
       if (ierr /= 0) e(1:ne)%M = read_int(UELIP,sln,'ellipse ')
       if (any(e%M < 1)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) # matching locations '//&
                     &'(e%M) must not be < 1 ',e%M
          stop 227
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%MS; sln=sln+1 ! bg%ms read above
       ! these are computed now (this value is now used as an alternate minimum),
       ! so checking this here doesn't mean as much as it once did.
       if (ierr /= 0) e(1:ne)%MS = read_int(UELIP,sln,'ellipse ')      

       read(UELIP,*,iostat=ierr) e(1:ne)%ibnd; sln=sln+1
       if (ierr /= 0) e(1:ne)%ibnd = read_int(UELIP,sln,'ellipse ')      
       if (any(e%ibnd < -1 .or. e%ibnd > 2)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) boundary type (e%ibnd) '//&
                     &'must be in {-1,0,1,2} ',e%ibnd
          stop 229
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%calcIn; sln=sln+1
       if (ierr /= 0) e(1:ne)%calcIn = read_logical(UELIP,sln,'ellipse ')      

       read(UELIP,*,iostat=ierr) e(1:ne)%storIn; sln=sln+1
       if (ierr /= 0) e(1:ne)%storIn = read_logical(UELIP,sln,'ellipse ')

       read(UELIP,*,iostat=ierr) e(1:ne)%wave; sln=sln+1
       if (ierr /= 0) then
           backspace(UELIP)
           read(UELIP,*,iostat=ierr) e(1)%wave
           if (ierr /= 0) then
             ! optional, assume false
             e(1:ne)%wave = .false.
             backspace(UELIP)
           else
             e(1:ne)%wave = e(1)%wave
           end if
       end if


       read(UELIP,*,iostat=ierr) e(1:ne)%r; sln=sln+1   ! eta
       if (ierr /= 0) e(1:ne)%r = read_real(UELIP,sln,'ellipse ')      
       if (any(e%r < 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) elliptical radius '//&
                     &'(e%r) must be >= 0.0 ',e%r
          stop 230
       end if

       where (e%r <= 0.0 .and. e%calcIn)
          ! can't calculate inside something of zero radius
          e%calcIn = .false.
       end where

       read(UELIP,*,iostat=ierr) e(1:ne)%x; sln=sln+1
       if (ierr /= 0) e(1:ne)%x = read_real(UELIP,sln,'ellipse ')       

       read(UELIP,*,iostat=ierr) e(1:ne)%y; sln=sln+1
       if (ierr /= 0) e(1:ne)%x = read_real(UELIP,sln,'ellipse ')
       
       ! NOTE: shift is computed from range of calc locations       
       e(1:ne)%x = e%x - s%xshift
       e(1:ne)%y = e%y - s%yshift
       e(1:ne)%z = cmplx(e%x, e%y, DP)

       read(UELIP,*,iostat=ierr) e(1:ne)%f; sln=sln+1
       if (ierr /= 0) e(1:ne)%f = read_real(UELIP,sln,'ellipse ')
       if (any(e%f <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) elliptical semi-focal '//&
                     &'length (e%f) must be > 0.0 ',e%f
          stop 231
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%theta; sln=sln+1
       if (ierr /= 0) e(1:ne)%theta = read_real(UELIP,sln,'ellipse ')
       if (any(e%theta < -PI) .or. any(e%theta > PI)) then
          write(stderr,*) 'WARNING: line ',sln,' ellipse input; eliptical angle wrt '//&
                     &'global cartesian (e%theta) must be -pi <= theta <= +pi',e%theta
          write(stderr,*) '*** resetting to inside range -pi <= theta <= +pi ***'
          where (e%theta < -PI)
            e%theta = modulo(e%theta,-PI)
          elsewhere (e%theta > PI)
            e%theta = modulo(e%theta,PI)
          end where
          write(stderr,*) 'modified angles:',e%theta
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%k; sln=sln+1
       if (ierr /= 0) e(1:ne)%k = read_real(UELIP,sln,'ellipse ')
       if (any(e%k <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) hydraulic conductivity '//&
                     &'(e%k) must be > 0.0 ',e%k
          stop 233
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%Ss; sln=sln+1
       if (ierr /= 0) e(1:ne)%Ss = read_real(UELIP,sln,'ellipse ')
       if (any(e%Ss <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) specific storage '//&
                     &'(e%Ss) must be > 0.0 ',e%Ss
          stop 234
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%por; sln=sln+1
       if (ierr /= 0) e(1:ne)%por = read_real(UELIP,sln,'ellipse ')
       if (any(e%por <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) porosity (e%por) '//&
                     &'must be > 0.0 ',e%por
          stop 235
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%leakFlag; sln=sln+1 
       if (ierr /= 0) e(1:ne)%leakFlag = read_int(UELIP,sln,'ellipse ')
       if (any(e%leakFlag < 0) .or. any(e%leakFlag > 3)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) leak flag (e%leakFlag) '&
               &//'must be in {0,1,2,3}'
          stop 2350
       end if       

       read(UELIP,*,iostat=ierr) e(1:ne)%aquitardK; sln=sln+1
       if (ierr /= 0) e(1:ne)%aquitardK = read_real(UELIP,sln,'ellipse ')
       if (any(e%aquitardK <= 0.0 .and. e%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line',sln,' ellipse input) aquitard vertical K '//&
                     &'(e%aquitardK) must be > 0.0 ',e%aquitardK
          stop 236
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%aquitardSs; sln=sln+1
       if (ierr /= 0) e(1:ne)%aquitardSs = read_real(UELIP,sln,'ellipse ')
       if (any(e%aquitardSs <= 0.0 .and. e%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line',sln,' ellipse input) aquitard specific '//&
                     &'storage (e%aquitardSs) must be > 0.0 ',e%aquitardSs
          stop 237
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%aquitardb; sln=sln+1
       if (ierr /= 0) e(1:ne)%aquitardb = read_real(UELIP,sln,'ellipse ')
       if (any(e%aquitardb <= 0.0 .and. e%leakFlag > 0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) aquitard thickness '//&
                     &'(e%aquitardB) must be > 0.0 ',e%aquitardb
          stop 238
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%unconfinedFlag; sln=sln+1
       if (ierr /= 0) e(1:ne)%unconfinedFlag = read_logical(UELIP,sln,'ellipse ')       

       read(UELIP,*,iostat=ierr) e(1:ne)%Sy; sln=sln+1
       if (ierr /= 0) e(1:ne)%Sy = read_real(UELIP,sln,'ellipse ')
       if (any(e%Sy <= 0.0 .and. e%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) specific storage '//&
                     &'(e%Sy) must be > 0.0 ',e%Sy
          stop 239
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%Kz; sln=sln+1
       if (ierr /= 0) e(1:ne)%Kz = read_real(UELIP,sln,'ellipse ')
       if (any(e%Kz <= 0.0 .and. e%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) aquifer vertical K '//&
                     &'(e%Kz) must be > 0.0 ',e%Kz
          stop 240
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%b; sln=sln+1
       if (ierr /= 0) e(1:ne)%b = read_real(UELIP,sln,'ellipse ')
       if (any(e%b <= 0.0 .and. e%unconfinedFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) aquifer thickness '//&
                     &'(e%b) must be > 0.0 ',e%b
          stop 241
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%dualPorosityFlag; sln=sln+1
       if (ierr /= 0) e(1:ne)%dualPorosityFlag = read_logical(UELIP,sln,'ellipse ')

       read(UELIP,*,iostat=ierr) e(1:ne)%matrixSs; sln=sln+1
       if (ierr /= 0) e(1:ne)%matrixSs = read_real(UELIP,sln,'ellipse ')
       if (any(e%matrixSs <= 0.0 .and. e%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) matrix specific '//&
                     &'storage (e%matrixSs) must be > 0.0', e%matrixSs
          stop 2411
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%lambda; sln=sln+1
       if (ierr /= 0) e(1:ne)%lambda = read_real(UELIP,sln,'ellipse ')
       if (any(e%lambda < 0.0 .and. e%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) matrix/fracture '//&
                     &'connection (e%lambda)'//&
               &' must be non-negative', e%lambda
          stop 2412
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%multiporosityDiffusion; sln=sln+1
       if (ierr /= 0) e(1:ne)%multiporosityDiffusion = read_int(UELIP,sln,'ellipse ')
       if (any((e%multiporosityDiffusion < 0 .or. e%multiporosityDiffusion > 3) &
            & .and. e%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) '//&
               &'multiporosity diffusion index must be {0,1,2,3}', e%multiporosityDiffusion
          stop 2412
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%kappa; sln=sln+1
       if (ierr /= 0) e(1:ne)%kappa = read_real(UELIP,sln,'ellipse ')
       if (any(e%kappa < 0.0 .and. e%dualPorosityFlag)) then
         write(stderr,*) 'ERROR: value (line',sln,' ellipse input) '//&
              & 'matrix/fracture hydraulic conductivity ratio (e%kappa) must be non-negative', e%kappa
         stop 2413
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%NDiffTerms; sln=sln+1
       if (ierr /= 0) e(1:ne)%NDiffTerms = read_int(UELIP,sln,'ellipse ')       
       if (any(e%NDiffTerms < 0 .and. e%dualPorosityFlag)) then
          write(stderr,*) 'ERROR: value (line',sln,' ellipse input) '//&
               & 'number terms in leaky diffusion series (e%nDiffTerms) must be >=0', e%NDiffTerms
          stop 2234
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%dskin; sln=sln+1
       if (ierr /= 0) e(1:ne)%dskin = read_real(UELIP,sln,'ellipse ')
       if (any(e%dskin <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' ellipse input) dimensionless '//&
                     &'skin (e%dskin) must be > 0.0 ',e%dskin
          stop 242
       end if

       read(UELIP,*,iostat=ierr) e(1:ne)%areaQ; sln=sln+1
       if (ierr /= 0) e(1:ne)%areaQ = read_real(UELIP,sln,'ellipse ') 

       read(UELIP,*,iostat=ierr) e(1:ne)%bdryQ; sln=sln+1
       if (ierr /= 0) e(1:ne)%bdryQ = read_real(UELIP,sln,'ellipse ') 

       where (e(:)%ibnd == -1 .or. e(:)%ibnd == 0 .or. e(:)%ibnd == +1)
          e(:)%match = .true.
       elsewhere
          e(:)%match = .false.
       end where

       write(chint,'(I4.4)') dom%num(2)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(UECHO,'(A)') '=============== ELLIPTICAL ELEMENT PROPERTIES ==============='
       write(UECHO,'(I0,A)') dom%num(2),        '  ||   # elliptical elements (including lines)'
       write(UECHO,fmt(1)) e(:)%n,              '  ||   # elliptical free parameter (Fourier coeffs)'
       write(UECHO,fmt(1)) e(:)%m,              '  ||   # ellipse matching locations'
       write(UECHO,fmt(1)) e(:)%ms,             '  ||   size "infinite" Mathieu matrices'
       write(UECHO,fmt(1)) e(:)%ibnd,           '  ||   ellipse ibnd array '//&
            &'(-1: spec head, 0: matching, :1 spec total flux, 2: spec elem flux)'
       write(UECHO,fmt(2)) e(:)%match,          '  ||   ellipse matching array'
       write(UECHO,fmt(2)) e(:)%calcin,         '  ||   calculate inside this ellipse?'
       write(UECHO,fmt(2)) e(:)%storin,         '  ||   calculate free-water storage for this ellipse?'
       write(UECHO,fmt(2)) e(:)%wave,           '  ||   calculate wave equation inside ellipse?'
       write(UECHO,fmt(3)) e(:)%r,              '  ||   ellipse radius (eta)'
       write(UECHO,fmt(3)) e(:)%x+s%xshift,     '  ||   original ellipse center x'
       write(UECHO,fmt(3)) e(:)%x,              '  ||   shifted ellipse center x'
       write(UECHO,fmt(3)) e(:)%y+s%yshift,     '  ||   original ellipse center y'
       write(UECHO,fmt(3)) e(:)%y,              '  ||   shifted ellipse center y'
       write(UECHO,fmt(3)) e(:)%f,              '  ||   ellipse semi-focal length'
       write(UECHO,fmt(3)) e(:)%theta,          '  ||   ellipse angle rotation with +x axis'
       write(UECHO,fmt(3)) e(:)%k,              '  ||   ellipse aquifer k'
       write(UECHO,fmt(3)) e(:)%ss,             '  ||   ellipse aquifer Ss'
       write(UECHO,fmt(3)) e(:)%por,            '  ||   ellipse aquifer porosity'
       write(UECHO,fmt(1)) e(:)%leakFlag,       '  ||   ellipse leaky type'
       write(UECHO,fmt(3)) e(:)%aquitardK,      '  ||   ellipse leaky aquitard K'
       write(UECHO,fmt(3)) e(:)%aquitardSs,     '  ||   ellipse leaky aquitard Ss'
       write(UECHO,fmt(3)) e(:)%aquitardb,      '  ||   ellipse leaky aquitard thickness'
       write(UECHO,fmt(2)) e(:)%unconfinedFlag, '  ||   ellipse unconfined flag'
       write(UECHO,fmt(3)) e(:)%Sy,             '  ||   ellipse aquifer specific yield'
       write(UECHO,fmt(3)) e(:)%Kz,             '  ||   ellipse aquifer vertical K'
       write(UECHO,fmt(3)) e(:)%b,              '  ||   ellipse aquifer thickness'
       write(UECHO,fmt(2)) e(:)%dualPorosityFlag,'  ||     ellipse dual porosity flag'
       write(UECHO,fmt(3)) e(:)%matrixSs,        '  ||     ellipse matrix Ss'
       write(UECHO,fmt(3)) e(:)%lambda,          '  ||     ellipse matrix/fracture connection lambda'
       write(UECHO,fmt(1)) e(:)%multiporosityDiffusion,  '  ||   ellipse multiporosity diffusion index'
       write(UECHO,fmt(3)) e(:)%kappa,          '  ||   ellipse matrix/fracture k ratio'
       write(UECHO,fmt(1)) e(:)%NDiffTerms,          '  ||   ellipse matrix diffusion series number terms'
       write(UECHO,fmt(3)) e(:)%dskin,           '  ||   ellipse boundary dimensionless skin factor'
       write(UECHO,fmt(3)) e(:)%areaQ,          '  ||   ellipse area rch rate'
       write(UECHO,fmt(3)) e(:)%bdryQ,          '  ||   ellipse boundary rch rate or head'

       do j = 1,size(e,dim=1)
          read(UELIP,*,iostat=ierr) e(j)%AreaTime; sln=sln+1
          if (ierr /= 0) then
             write(stderr,*) 'ERROR reading line ',sln,' area time behavior (e%AreaTime) &
                  &ellipse ',j,' input'
             stop 2422
          end if
          call read_time_behaviors(UELIP,UECHO,e(j)%time,j,sln,'ellipse',.true.)
       end do
       do j = 1,size(e,dim=1)
          read(UELIP,*,iostat=ierr) e(j)%BdryTime; sln=sln+1
          if (ierr /= 0) then
             write(stderr,*) 'ERROR reading line ',sln,' boundary time behavior '//&
                  &'(e%BdryTime) ellipse ',j,' input'
             stop 2425
          end if
          call read_time_behaviors(UELIP,UECHO,e(j)%time,j,sln,'ellipse',.false.)
       end do

       close(UELIP) ! ellipse input file
       
       ! TODO: handle free-water storage for ellipses?
       if(any(e%storin)) then
          write(stdout,*) 'WARNING: wellbore (free-water) storage not handled for ellipses yet'
          write(stdout,*) '**** resetting to false ****'
          e%storin = .false.
       end if

    else
       ! no elliptical elements
       allocate(e(0))
    end if

    ntot = sum(dom%num) ! total number of circular and elliptical elements
    if (ntot < 1) then
      stop 'READINPUT: Need at least one circular (including well) or'//&
           &'elliptical (including line) element.'
    end if
    
    ! compute secondary parameters
    bg%alpha = bg%K/bg%Ss
    c(:)%alpha = c%K/c%Ss
    e(:)%alpha = e%K/e%Ss
    
    bg%T = bg%K*bg%b
    c(:)%T = c%K*c%b
    e(:)%T = e%K*e%b

    write(chint,'(I4.4)') dom%num(1)
    fmt(2) = '('//chint//'(ES11.5,1X),A) ' ! circles
    write(chint,'(I4.4)') dom%num(2)
    fmt(3) = '('//chint//'(ES11.5,1X),A) ' ! ellipses

    write(UECHO,'(A)') '=============== COMPUTED PROPERTIES ==============='
    if (dom%num(1) > 0) then
       write(UECHO,fmt(2)) c(:)%alpha,'  ||    circle hydraulic diffusivity'
       write(UECHO,fmt(2)) c(:)%T,    '  ||    circle transmissivity'
    end if
    if (dom%num(2) > 0) then
       write(UECHO,fmt(3)) e(:)%alpha,'  ||    ellipse hydraulic diffusivity'
       write(UECHO,fmt(3)) e(:)%T,    '  ||    ellipse transmissivity'
    end if
    write(UECHO,'(ES12.5,A)') bg%alpha,'  ||    background hydraulic diffusivity'
    write(UECHO,'(ES12.5,A)') bg%T,    '  ||    background transmissivity'

    ! particles
    if (s%particle) then
       sln = 1
       read(UINPUT,*,iostat=ierr) particleFname; ln=ln+1
       if (ierr /= 0) then
          write(stderr,*) 'ERROR: line ',ln,' (particles: particle file name) of',s%infname
          stop 2500
       end if
       open(unit=UPAR, file=particleFname, status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(2A)') 'ERROR: cannot open particle input file for reading ',&
               & trim(particleFname)
          stop 243
       else
          write(UECHO,'(A)') trim(particleFname)//' opened for particle input data'
       end if

       read(UPAR,*,iostat=ierr) s%nPart,  s%streakSkip
       if (ierr /= 0 .or. (s%nPart < 1) .or. (s%streakSkip < 0 .and. s%output == 21)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) s%nPart and &
               &s%streakSkip must be >0',[s%nPart,s%streakSkip]
          stop 2430
       end if
       allocate(p(s%nPart))

       read(UPAR,*,iostat=ierr) p(:)%forward; sln=sln+1 
       if (ierr /= 0) p(:)%forward = read_logical(UPAR,sln,'particle')

       read(UPAR,*,iostat=ierr) p(:)%int; sln=sln+1
       if (ierr /= 0) p(:)%int = read_int(UPAR,sln,'particle')
       if (any(p%int < 1 .or. p%int > 4)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) integration '//&
                     &'method (p%int) must be {1,2,3,4} ',p%int
          stop 252
       end if

       read(UPAR,*,iostat=ierr) p(:)%tol; sln=sln+1
       if (ierr /= 0) p(:)%tol = read_real(UPAR,sln,'particle')
       if (any(p%tol <= 0.0 .and. p%int == 1)) then
          write(stderr,*) 'ERROR: value (line',sln,'particle input) RKM error '//&
                     &'tolerance (p%tol) must be > 0.0 ',p%tol
          stop 244
       end if

       read(UPAR,*,iostat=ierr) p(:)%maxL; sln=sln+1  
       if (ierr /= 0) p(:)%maxL = read_real(UPAR,sln,'particle')
       if (any(p%maxL <= 0.0 .and. p%int == 1)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) max RKM '//&
                     &'step length (p%maxL) must be > 0.0 ',p%maxL
          stop 245
       end if

       read(UPAR,*,iostat=ierr) p(:)%mindt; sln=sln+1
       if (ierr /= 0) p(:)%mindt = read_real(UPAR,sln,'particle')
       if (any(p%mindt <= 0.0 .and. p%int == 1)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) min RKM '//&
                     &'time step (p%mindt) must be > 0.0 ',p%mindt
          stop 246
       end if

       read(UPAR,*,iostat=ierr) p(:)%dt; sln=sln+1
       if (ierr /= 0) p(:)%dt = read_real(UPAR,sln,'particle')
       if (any(p%dt <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) initial '//&
                     &'time step (p%dt) must be > 0.0 ',p%dt
          stop 247
       end if

       read(UPAR,*,iostat=ierr) p(:)%x; sln=sln+1
       if (ierr /= 0) p(:)%x = read_real(UPAR,sln,'particle')

       read(UPAR,*,iostat=ierr) p(:)%y; sln=sln+1
       if (ierr /= 0) p(:)%y = read_real(UPAR,sln,'particle')

       p(:)%x = p%x - s%xshift
       p(:)%y = p%y - s%yshift

       read(UPAR,*,iostat=ierr) p(:)%ti; sln=sln+1
       if (ierr /= 0) p(:)%ti = read_real(UPAR,sln,'particle')
       if (any(p%ti <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) start '//&
                     &'time (p%ti) must be > 0.0 ',p%ti
          stop 248
       end if

       read(UPAR,*,iostat=ierr) p(:)%tf; sln=sln+1
       if (ierr /= 0) p(:)%tf = read_real(UPAR,sln,'particle')
       if (any(p%tf <= 0.0)) then
          write(stderr,*) 'ERROR: value (line ',sln,' particle input) end '//&
                     &'time (p%tf) must be > 0.0 ',p%tf
          stop 249
       end if

       if (any(p%tf < p%ti .and. p%forward)) then
          write(stderr,*) 'ERROR: value final t must be > initial t for forward tracking ti:',p%ti,' tf:',p%tf
          stop 250
       end if
       if (any(p%tf > p%ti .and. .not.p%forward)) then
          write(stderr,*) 'ERROR: value final t must be < initial t for backward tracking ti:',p%ti,' tf:',p%tf
          stop 251
       end if

       read(UPAR,*,iostat=ierr) p(:)%inclIn; sln=sln+1
       if (ierr /= 0) p(:)%inclIn = read_logical(UPAR,sln,'particle')

       close(UPAR) ! particle input file

       write(chint,'(I4.4)') s%nPart
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(UECHO,fmt(3)) p(:)%tol,    '  ||    particle sution tolerances (RKM only)'
       write(UECHO,fmt(3)) p(:)%dt,     '  ||    particle time step (RKM initial step)'
       write(UECHO,fmt(3)) p(:)%maxL,   '  ||    particle max step length (RKM only)'
       write(UECHO,fmt(3)) p(:)%mindt,  '  ||    particle min time step size (RKM only)'
       write(UECHO,fmt(3)) p(:)%x+s%xshift,  '  ||    original particle initial x'
       write(UECHO,fmt(3)) p(:)%x,           '  ||    shifted particle initial x'
       write(UECHO,fmt(3)) p(:)%y+s%yshift,  '  ||    original particle initial y'
       write(UECHO,fmt(3)) p(:)%y,           '  ||    shifted particle initial y'
       write(UECHO,fmt(3)) p(:)%ti,     '  ||    particle initial t'
       write(UECHO,fmt(3)) p(:)%tf,     '  ||    particle maximum t'
       write(UECHO,fmt(1)) p(:)%int,    '  ||    particle integration method '//&
            & '(1=RKM, 2=RK, 3=root, 4=Euler)'
       write(UECHO,fmt(2)) p(:)%InclIn, '  ||    particle begins inside CH/CF incl?'
    else
       ! no particles
       allocate(p(0))
    end if
    close(UINPUT) ! main input file
    close(UECHO) ! input echo file

    if (s%debug) then
      ! solution struct not always passed into all routines
      bg%debug = .true.
      c%debug = .true.
      e%debug = .true.
      p%debug = .true.
    else
      bg%debug = .false.
      c%debug = .false.
      e%debug = .false.
      p%debug = .false.
    end if
    
  end subroutine readInput

  subroutine read_time_behaviors(UIN,UECH,el,j,ln,tp,area) 
    use type_definitions, only : time, explain_type
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    integer, intent(in) :: UIN, UECH, ln, j
    type(time), intent(inout) :: el
    character(7), intent(in) :: tp
    logical, intent(in) :: area

    type(explain_type) :: explain
    character(46) :: lfmt
    character(8) :: lincon
    integer :: ierr, tsize

    lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
    
    if             (area .and. el%AreaTime < -100 .or. &
         & ((.not. area) .and. el%BdryTime < -100)) then
       lincon = 'linear'
    else
       lincon = 'constant'
    end if

    backspace(UIN)
    if (area) then
        if (el%AreaTime > -1) then
           ! functional time behavior
           allocate(el%ATPar(2))
           read(UIN,*,iostat=ierr) el%AreaTime,el%ATPar(:)
           if (ierr /= 0) then
              write(stderr,*) 'ERROR: line ',ln,' area functional time behavior '//&
                   &'(ATPar) ',tp,j,' input'
              stop 2202
           end if
           write(UECH,'(I0,2(1X,ES12.5),A,I0)') el%AreaTime, el%ATPar(:),&
                &'  ||  Area time behavior '//trim(explain%time(el%areaTime))//&
                &', par1, par2 for '//tp,j
        else
           ! piecewise-constant/linear time behavior           
           tsize = mod(abs(el%AreaTime),100) ! handle both piecewise constant (-100:-1) and piecewise linear (-infinity:-101)
           allocate(el%ATPar(2*tsize+1))
           read(UIN,*,iostat=ierr) el%AreaTime,el%ATPar(:)
           if (ierr /= 0) then
              write(stderr,*) 'ERROR: line ',ln,' area piecewise-'//&
                   &trim(lincon)//' time behavior (ATpar) ',tp,j,'input'
              stop 2203
           end if
           write(lfmt(8:11), '(I4.4)') size(el%ATPar(:tsize+1),1)
           write(lfmt(26:29),'(I4.4)') size(el%ATPar(tsize+2:),1)
           write(UECH,lfmt) el%AreaTime,el%ATPar(:tsize+1),'| ',&
                & el%ATPar(tsize+2:), &
                &'  ||    Area ti, tf | piecewise-'//trim(lincon)//' strength for '//tp,j
        end if
     else
        if (el%BdryTime > -1) then
           allocate(el%BTPar(2))
           read(UIN,*,iostat=ierr) el%BdryTime,el%BTPar(:)
           if (ierr /= 0) then
              write(stderr,*) 'ERROR: line ',ln,' boundary functional time '//&
                   &'behavior (BTPar) ',tp,j,' input'
              stop 2205
           end if
           write(UECH,'(I0,2(1X,ES12.5),A,I0)') el%BdryTime, el%BTPar(:),&
                &'  ||  Bdry time behavior '//trim(explain%time(el%BdryTime))//&
                &', par1, par2 for '//tp,j
        else
           tsize = mod(abs(el%BdryTime),100)
           allocate(el%BTPar(2*tsize+1))
           read(UIN,*,iostat=ierr) el%BdryTime,el%BTPar(:)
           if (ierr /= 0) then
              write(stderr,*) 'ERROR: line ',ln,' boundary piecewise- '//&
                   &trim(lincon)//' time behavior (BTpar) ',tp,j,'input'
              stop 2203
           end if
           write(lfmt(8:11), '(I4.4)') size(el%BTPar(:tsize+1),1)
           write(lfmt(26:29),'(I4.4)') size(el%BTPar(tsize+2:),1)
           write(UECH,lfmt) el%BdryTime,el%BTPar(:tsize+1),' | ',&
                & el%BTPar(tsize+2:), &
                &'  ||    Bdry ti, tf | piecewise-'//trim(lincon)//' strength for '//tp,j
        end if
     end if

   end subroutine read_time_behaviors

  !******************************************************
  function computeVector(unit,n,line) result(v)
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
        
    use constants, only : DP
    integer, intent(in) :: unit, n, line
    real(DP), dimension(n) :: v
    integer :: ierr,i
    real(DP) :: minv,maxv,delta
    character(6) :: vec

    ! assume file is positioned at beginning of line
    ! if 'linvec','LINVEC','logvec', or 'LOGVEC' are the first
    ! entry on the line, read two more integers (min/max)
    ! and create an evenly spaced vector n elements long
    ! line indicates the input files line, for error reporting

    backspace(unit)
    read(unit,*,iostat=ierr) vec

    if ((vec(4:6) /= 'vec' .and. vec(4:6) /= 'VEC') &
         & .or. ierr /= 0) then
       write(stderr,*) 'ERROR: line ',line,' input '//&
            &'({LIN,LOG}VEC min max) input' 
       stop 7770
    else
       backspace(unit)
       read(unit,*,iostat=ierr) vec,minv,maxv
       if (ierr /= 0) then
          write(stderr,*) 'ERROR: line ',line,&
               &' calc vector ({LIN,LOG}VEC min max) input'
          stop 7771
       else
          if (n == 1 ) then
             v = (maxv+minv)/2.0 ! average of min/max
          else
             delta = (maxv-minv)/real(n-1,DP)
             v = minv + real([(i,i=0,n-1)],DP)*delta
          end if

          if (vec(1:3) == 'LIN' .or. vec(1:3) == 'lin') then
             return
          elseif (vec(1:3) == 'LOG' .or. vec(1:3) == 'log') then
             v = 10.0_DP**v
             return
          else
             write(stderr,*) 'ERROR: vector type line ',line,&
                  & '{LIN,LOG}VEC or {lin,log}vec'
             stop 7772
          end if
       end if
    end if

  end function computeVector
  
  !******************************************************
  subroutine writeResults(s,p)

    use type_definitions, only : solution, particle, explain_type
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit, stderr => error_unit

    type(solution), intent(inout) :: s
    type(particle), dimension(:), intent(inout) :: p

    type(explain_type) :: explain
    character(4), dimension(2) :: chint
    integer :: i, j, k, nt

    ! adjust the formats of time, location, and results here
    character(6), parameter :: tfmt = 'ES13.5', xfmt = 'ES12.4'
    character(9), parameter :: hfmt = 'ES22.14e3'

    ! remove shift applied at beginning
    s%x(:) = s%x(:) + s%xshift
    s%y(:) = s%y(:) + s%yshift
    if (s%particle) then
       do i = 1,s%nPart
          p(i)%r(:,2) = p(i)%r(:,2) + s%xshift
          p(i)%r(:,3) = p(i)%r(:,3) + s%yshift
       end do
    end if
    
    select case (s%output)
    case (1)
       ! ** Gnuplot contour map friendly output **
       ! print results as x,y,{h,v,dh} "triplets" with the
       ! times separated by double blank lines

       open(unit=20, file=s%outfname, status='replace', action='write')
       write(20,'(A)') '# 1 LT-AEM contour map output     -*-auto-revert-*-'
       write(20,'(A,I0)') '# t: ', s%nt
       write(20,'(A,I0)') '# x: ', s%nx
       write(20,'(A,I0)') '# y: ', s%ny
       write(20,'(A,I0)') '# locations: ', s%nx*s%ny

       do i = 1, s%nt
          write(20,'(A,'//tfmt//')') '# t= ',s%t(i)
          write(20,'(A)',advance='no')   &
          & '#      X           Y               head'//&
          & '                velx                  vely'
          if (s%deriv) then
             write(20,'(A)') '                d(head)/d(log(t))'
          else
             write(20,'(A)') ''
          end if
          if (s%deriv) then
             do j = 1, s%ny
                do k = 1, s%nx
                   write(20,'(2('//xfmt//',1X),4('//hfmt//',1X))') &
                        & s%x(k), s%y(j), s%h(i,k,j), s%v(i,1:2,k,j), s%dh(i,k,j)
                end do
             end do
          else
             do j = 1, s%ny
                do k = 1, s%nx
                   write(20,'(2('//xfmt//',1X),3('//hfmt//',1X))') &
                        & s%x(k), s%y(j), s%h(i,k,j), s%v(i,1:2,k,j)
                end do
             end do
          end if
          write(20,'(/)')
       end do
       write(20,'(A)') '# EOF'
       close(20)

       write(stdout,'(/A)') '**************************************************'
       write(stdout,'(2A)') ' Gnuplot contour output => ', trim(s%outfname)
       write(stdout,'(A)')  '**************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (2)

       ! ** matrix/table contour output **
       ! print results as matrices, with each variable (at each time) going
       ! to a separate file x and y matrices - similar to results from Matlab
       ! function meshgrid()

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
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%h(k,j,i), j=1,s%nx)
          end do
          close(20)

          ! log-t derivative matrix
          if (s%deriv) then
             open(unit=20, file=trim(s%outfname)//'_dhead_'//chint(2)//'.dat', &
                  & status='replace', action='write')
             do i = 1, s%ny
                write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%dh(k,j,i), j=1,s%nx)
             end do
             close(20)
          end if

          ! velx-matrix
          open(unit=20, file=trim(s%outfname)//'_velx_'//chint(2)//'.dat', &
               & status='replace', action='write')
          do i = 1, s%ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(k,1,j,i), j=1,s%nx)
          end do
          close(20)

          ! vely-matrix
          open(unit=20, file=trim(s%outfname)//'_vely_'//chint(2)//'.dat', &
               & status='replace', action='write')
          do i = 1, s%ny
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%v(k,2,j,i), j=1,s%nx)
          end do
          close(20)
       end do

       ! column of calculation times
       open(unit=20, file=trim(s%outfname)//'_t.dat', status='replace', &
            & action='write')
       write (20,'('//tfmt//')') (s%t(j), j=1,s%nt)
       close(20)

       write(stdout,'(/A)') '*********************************************************'
       write(stdout,'(3A)') 'table output => ', trim(s%outfname), &
              & '{x,y,t,{d,}head{1-n},velx{1-n},vely{1-n}}.dat'
       write(stdout,'(A)') '**********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (10)

       ! ** Gnuplot-friendly time series output **
       ! column of time values at a location through time
       ! locations separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# 10 LT-AEM time series output   -*-auto-revert-*-'
       write (20,'(A,2(I0,A))') '# ',s%nx,' locations ',s%nt,' times'
       do i = 1, s%nx
          write (20,'(2(A,'//xfmt//'),3X,A)') '# location: x=',s%x(i),' y=',&
               & s%y(i),trim(s%obsname(i))
          write (20,'(A)',advance='no')  '#     time              head'//&
               & '                  velx                vely'
          if (s%deriv) then
             write(20,'(A)') '                 deriv'
          else
             write(20,'(A)') ''
          end if

          if (s%deriv) then
             do k = 1, s%nt
                write (20,'(1X,'//tfmt//',4(1X,'//hfmt//'))') &
                     & s%t(k),s%h(k,i,1),s%v(k,1:2,i,1),s%dh(k,i,1)
             end do             
          else
             do k = 1, s%nt
                write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                     & s%t(k),s%h(k,i,1),s%v(k,1:2,i,1)
             end do
          end if
          write (20,'(/)')
       end do
       write(20,*) '# EOF'
       close(20)

       write(stdout,'(/A)') '*****************************************************'
       write(stdout,'(2A)') 'Gnuplot timeseries output => ', trim(s%outfname)
       write(stdout,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (11)

       ! ** Gnuplot-friendly time series output **
       ! column of time values at a location through time
       ! locations separated by blank lines (no velocity)
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# 11 LT-AEM time series output    -*-auto-revert-*-'
       write (20,'(A,2(I0,A))') '# ',s%nx,' locations ',s%nt,' times'
       do j = 1, s%nx
          write (20,'(2(A,'//xfmt//'),3X,A)') '# location: x=',s%x(j),' y=',&
               &s%y(j),trim(s%obsname(j))
          if (s%deriv) then
             write (20,'(A)')   '#     time              head             deriv'
          else
             write (20,'(A)')   '#     time              head'
          end if
          if (s%deriv) then
             do k = 1, s%nt
                write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                     & s%t(k),s%h(k,j,1),s%dh(k,j,1)
             end do
          else
             do k = 1, s%nt
                write (20,'(1X,'//tfmt//',2(1X,'//hfmt//'))') &
                     & s%t(k),s%h(k,j,1)
             end do
          end if
          write (20,'(/)')
       end do
       write(20,*) '# EOF'
       close(20)

       write(stdout,'(/A)') '*****************************************************'
       write(stdout,'(2A)') 'Gnuplot timeseries output (no v) => ', trim(s%outfname)
       write(stdout,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (12)

       ! ** output for inverse in Matlab
       ! column of time values at a location through time
       ! locations separated by blank lines

       do i = 1, s%nx
          write(chint(1),'(I4.4)') i
          open(unit=20, file=trim(s%outfname)//'_'//chint(1), status='replace', &
               &action='write')
          if (s%deriv) then
             do k = 1, s%nt
                write (20,'('//tfmt//',2(1X,'//hfmt//'))') s%t(k),s%h(k,i,1),s%dh(k,i,1)
             end do
          else
             do k = 1, s%nt
                write (20,'('//tfmt//',1X,'//hfmt//')') s%t(k),s%h(k,i,1)
             end do
          end if
          write(20,'(/)')
          close(20)
       end do

       write(stdout,'(/A)') '*****************************************************'
       write(stdout,'(4A)') 'inverse output => ',trim(s%outfname),'0000-',chint(1)
       write(stdout,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (20)

       ! ** pathline Gnuplot-style output **
       ! columns of time values for starting locations
       ! particles separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# 20 LT-AEM particle tracking output  -*-auto-revert-*-'
       write (20,'(A,I0,A)') '# ',size(p,dim=1),' particles'
       do i = 1, size(p,dim=1)
          write (20,'(A,I0,A)') '# particle ',i,' '//trim(explain%particle(p(i)%int))
          write (20,'(A)')   &
          & '#     time              x                    y                  '//&
          &'velx                 vely '
          do k = 1,size(p(i)%r,dim=1)
             write (20,'('//tfmt//',4(1X,'//hfmt//'))') &
                  & p(i)%r(k,1:5)
          end do
          write (20,'(/)')
       end do
       write(20,'(A)') '# EOF'
       close(20)

       write(stdout,'(/A)') '*****************************************************'
       write(stdout,'(2A)') 'particle tracking output => ', trim(s%outfname)
       write(stdout,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (21)

       ! ** streakline Gnuplot-style output **
       ! this requires constant time steps (can't use adaptive integration)
       ! each block is a requested time, each row a particle

       open(unit=90, file=s%outfname, status='replace', action='write')
       write (90,'(A)') '# 21 LT-AEM particle tracking streakfile output   -*-auto-revert-*-'
       write (20,'(A,I0,A)') '# ',size(p,dim=1),' particles'

       ! max number of times for all particles
       nt = maxval(p(:)%numt,dim=1)
       
       do i = 1, nt, s%streakSkip
          ! use maxval to ensure a non-zero time is reported
          write (90,'(A,'//tfmt//')') '# time:', maxval(p(i)%r(:,1))
          write (90,'(A)') '#  particle       x            y&
               &           velx         vely'
          do j = 1, size(p,dim=1)
             if (p(j)%r(i,1) > 0.0) then
                ! only write particle if it has non-zero data
                write (90,'(I0,4(1X,'//hfmt//'))')  j,p(j)%r(i,2:5)
             end if
          end do
          write (90,'(/)')
       end do
       write(90,'(A)') '# EOF'
       close(90)

       write(stdout,'(/A)') '*****************************************************'
       write(stdout,'(2A)') 'particle tracking streakfile => ', trim(s%outfname)
       write(stdout,'(A)')  '*****************************************************'

    case default
       write(stderr,'(A,I0)')  'invalid output code ',s%output
       stop 300
    end select

    ! write total flowrate through each element to separate file
    ! each element includes a timeseries of flowrate
    if (s%Qcalc) then
       open(unit=91, file=s%qfname, status='replace', action='write')
       write(91,'(A)') '# LT-AEM element total flowrates  -*-auto-revert-*-'
       
       do j = 1,size(s%Q,dim=2)
          write(91,'(A,I0)') '# element ',j
          if (s%deriv) then
             write(91,'(A)') '#     tD      '//'          Q_D          '&
                  & //'   d(Q_D)/d(log(t))'
          else
             write(91,'(A)') '#     tD      '//'          Q_D          '
          end if
          
          if (s%deriv) then
             do i = 1,s%nt
                write(91,'('//tfmt//',2(1X,'//hfmt//'))') s%t(i),s%Q(i,j),s%dQ(i,j)
             end do
          else
             do i = 1,s%nt
                write(91,'('//tfmt//',1X,'//hfmt//')') s%t(i),s%Q(i,j)
             end do
          end if
          write(91,'(/)')
       end do
       write(91,'(A)') '# EOF'
       close(91)
       write(stdout,'(A)') '%% wrote element flowrates: '//trim(s%qfname)
    end if

  end subroutine writeResults

  !##################################################
  subroutine writeGeometry(c,e,s)
    use constants, only : DP
    use type_definitions, only : circle, ellipse, solution
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit

    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    type(solution), intent(in) :: s
    integer :: nc, ne, i, j, ierr
    complex(DP) :: origin
    character(14) :: fmt

    ! remove shift originally applied to coordinates
    origin = cmplx(s%xshift,s%yshift,DP)

    nc = size(c,dim=1)
    ne = size(e,dim=1)
    fmt = '(3(ES13.5,1X))'

    ! write matching points to file
    open(unit=40, file=s%geomfName, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       ! non-fatal error
       write(stdout,'(2A)') 'WARNING: writeGeometry error opening output file for writing &
            &element matching locations ',s%geomfName
    else
       write(40,'(A)') '# points along circumference circular '//&
            &'and elliptical elements  -*-auto-revert-*-'
       do i = 1,nc
          write(40,'(2(A,I0),A)') '# circular element ',i,' = ',c(i)%M,' points'
          write(40,fmt)  (origin + c(i)%Zom(j), c(i)%Pcm(j), j=1,c(i)%M)
          if (c(i)%M > 1) then
             ! joins circle back up with beginning for plotting
             write(40,fmt)  origin + c(i)%Zom(1), c(i)%Pcm(1)
          end if
          write(40,'(/)')
       end do

       do i = 1,ne
          write(40,'(2(A,I0),A)') '# elliptical element ',i,' = ',e(i)%M,' points'
          write(40,fmt)  (origin + e(i)%Zom(j), e(i)%Pcm(j), j=1,e(i)%M)
          if (e(i)%M > 1) then
             ! joins ellipse back up with beginning for plotting
             write(40,fmt)  origin + e(i)%Zom(1), e(i)%Pcm(1)
          end if
          write(40,'(/)')
       end do
       write(40,'(A)') '# EOF'
       close(40)
    end if

  end subroutine writeGeometry
  
  function read_int(unit,num,str) result(ival)
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    integer, intent(in) :: unit, num
    integer :: ival, ierr
    character(8) :: str
    backspace(unit)
    read(unit,*,iostat=ierr) ival
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',num,' ',trim(str),' input'
       stop 7771
    end if
  end function read_int

  function read_real(unit,num,str) result(fval)
    use constants, only : DP
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    integer, intent(in) :: unit, num
    real(DP) :: fval
    integer :: ierr
    character(8) :: str
    backspace(unit)
    read(unit,*,iostat=ierr) fval
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',num,' ',trim(str),' input'
       stop 7772
    end if
  end function read_real

  function read_logical(unit,num,str) result(lval)
    use, intrinsic :: iso_fortran_env, only : stderr => error_unit
    integer, intent(in) :: unit, num
    logical :: lval
    integer :: ierr
    character(8) :: str
    backspace(unit)
    read(unit,*,iostat=ierr) lval
    if (ierr /= 0) then
       write(stderr,*) 'ERROR: line',num,' ',trim(str),' input'
       stop 7773
    end if
  end function read_logical
  
  subroutine dump_coeff(sol,c,e,nt,s,tee,minlt,maxlt)
    use type_definitions, only : solution, circle, ellipse
    use constants, only : DP
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit

    type(circle), intent(in), dimension(:) :: c
    type(ellipse), intent(in), dimension(:) :: e
    type(solution), intent(in) :: sol
    integer, intent(in), dimension(:) :: nt
    complex(DP), intent(in), dimension(:,:) :: s
    real(DP), intent(in), dimension(:) :: tee
    integer, intent(in) :: minlt,maxlt

    integer :: ierr,i,nc,ne
    nc = size(c)
    ne = size(e)

    open(unit=77, file=sol%coefffname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       write(stdout,*) 'WARNING: error opening intermediate save file ',&
            & trim(sol%coefffname), ' continuing without saving results'
    else
       write(77,'(A)') sol%echofName ! file with all the input parameters
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
  end subroutine dump_coeff

  subroutine read_coeff(sol,bg,c,e,nt,s,tee,minlt,maxlt,fail)
    use type_definitions, only : solution, circle, ellipse, element
    use constants, only : DP
    use ellipse_mathieu_init, only : ellipse_init
    use, intrinsic :: iso_fortran_env, only : stdout => output_unit, stderr => error_unit

    type(solution), intent(inout) :: sol
    type(element), intent(inout) :: bg
    type(circle), intent(inout), dimension(:) :: c
    type(ellipse), intent(inout), dimension(:) :: e
    integer, intent(out), allocatable :: nt(:)
    complex(DP), intent(out), allocatable :: s(:,:)
    real(DP), intent(out), allocatable :: tee(:)
    integer, intent(out) :: minlt,maxlt
    logical, intent(out) :: fail

    character(6) :: elType  ! element type {CIRCLE,ELLIPS}
    integer :: ierr,i,j,nc,ne,crow,ccol
    nc = size(c)
    ne = size(e)

    open(unit=77, file=sol%coefffname, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       ! go back and recalculate if no restart file
       write(stderr,'(A)') 'ERROR: cannot opening restart file, recalculating...'
       sol%calc = .true.
       fail = .true.
    else
       fail = .false.
    end if

    read(77,*) !! TODO check inputs are same?
    read(77,*) minlt,maxlt,nc,ne 
    allocate(s(2*sol%m+1,minlt:maxlt-1), nt(minlt:maxlt-1), tee(minlt:maxlt-1))
    read(77,*) nt(:)
    read(77,*) s(:,:)
    sol%totalnP = product(shape(s))
    read(77,*) tee(:)
    do i = 1,nc
        read(77,*) elType,j,crow,ccol
        if (elType == 'CIRCLE' .and. i == j) then
           allocate(c(i)%coeff(crow,ccol))
        else
           write(stderr,'(A)') 'ERROR reading in CIRCLE matching '//&
                &'results, recalculating...'
           sol%calc = .true.
           fail = .true.
        end if
        read(77,*) c(i)%coeff(:,:)
     end do
     do i = 1,ne
        read(77,*) elType,j,crow,ccol
        if (elType == 'ELLIPS' .and. i == j) then
           allocate(e(i)%coeff(crow,ccol))
        else
           write(stderr,'(A)') 'ERROR reading in ELLIPS matching '//&
                &'results, recalculating...'
           sol%calc = .true.
           fail = .true.
        end if
        read(77,*) e(i)%coeff(:,:)
     end do

     ! re-initialize Mathieu function matrices
     if (ne > 0) then
        write(stdout,'(A)') 're-computing Mathieu coefficients ...'
        call ellipse_init(e,bg,s)
     end if

     write(stdout,'(A)') 'matching results successfully re-read from file'
   end subroutine read_coeff
   

end module file_ops

