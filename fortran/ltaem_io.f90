!
! Copyright (c) 2011,2012,2013 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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
  public :: readInput, writeResults, writeGeometry

contains

  !##################################################
  ! this routine read the main input file, and allocates the main
  ! data structures used to store data.
  subroutine readInput(s,dom,bg,c,e,p)
    use constants, only : DP, lenFN, PI
    use type_definitions, only : solution, particle, domain, element, circle, ellipse

    type(solution), intent(inout) :: s
    type(particle), intent(out), allocatable :: p(:)
    type(domain),   intent(out) :: dom
    type(element),  intent(out) :: bg
    type(circle),   intent(out), allocatable :: c(:)
    type(ellipse),  intent(out), allocatable :: e(:)

    character(4) :: chint
    character(20), dimension(3) :: fmt
    character(512)  :: input ! a "big enough" buffer to contain all location names?
    character(46) :: lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)'
    character(lenFN+5) :: echofname
    character(lenFN) :: circleFname, ellipseFname, particleFname
    integer :: ierr, j, ntot, nC, nE, ln = 1, sln, s1,s2,slen


    open(unit=15, file=trim(s%infname), status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening input file ',trim(s%infname)
       stop 100
    endif

    echofname = trim(s%infname)//'.echo'
    open(unit=16, file=echofname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       write(*,'(2A)') 'READINPUT: error opening echo file ',echofname
       stop 101
    else
       ! add a file variable at top of file to set Emacs to auto-revert mode
       write(16,'(A)') '-*-auto-revert-*-'
    endif

    ! solution-specific and background aquifer parameters
    read(15,*,iostat=ierr) s%calc, s%particle, s%contour, s%deriv, s%output, &
         & s%outFname, s%coeffFName, s%elemHfName, s%geomfName
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (flags & filenames) of input'
       stop 110
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
       write(*,*) 'input file (line ',ln,') s%output must be in '&
            &//'{1,2,10,11,12,20,21} ',s%output
       stop 200
    end if
    if (s%output >= 10  .and. s%contour) then
       write(*,*) 'input file (line ',ln,') s%output should be in {1,2} '&
            &//'when contour output is selected ',s%output,s%contour
       stop 2010
    elseif (s%output < 10 .and. .not. s%contour) then
       write(*,*) 'input file (line ',ln,') s%output should be only in '&
            &//'{1,2} when contour output is selected ',s%output,s%contour
       stop 2011
    end if

    if ((s%output < 10 .or. s%output >= 20) .and. s%timeseries) then
       write(*,*) 'input file (line ',ln,') s%output should be in {10,11,12}'&
            &//' when timeseries output is selected ',s%output,s%timeseries
       stop 2020
    elseif (s%output >= 10 .and. s%output < 20 .and. .not. s%timeseries) then
       write(*,*) 'input file (line ',ln,') s%output should only be in {10,'&
            &//'11,12} when contour output is selected ',s%output,s%timeseries
       stop 2021
    end if

    if (s%output < 20 .and. s%particle) then
       write(*,*) 'input file (line ',ln,') if s%particle==.True., '&
            &//'s%output should be in {20,21} ',s%output,s%particle
       stop 2030
    elseif (s%output >= 20  .and. .not. s%particle) then
       write(*,*) 'input file (line ',ln,') s%output should only be in '&
            &//'{20,21} if s%particle==.True. ',s%output,s%particle
       stop 2031
    end if    

    read(15,*,iostat=ierr) bg%por, bg%k, bg%ss, bg%leakFlag, &
         & bg%aquitardK, bg%aquitardSs, bg%aquitardb, bg%ms, bg%cutoff; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error on line ',ln,' of input file (background props)'
       stop 204
    end if
    
    ! reals checked here, bg%ms checked in ellipse section
    if (any([bg%por,bg%k,bg%ss] <= 0.0)) then
       write(*,*) 'input file (line ',ln,') bg%por, bg%k, bg%ss &
            &must all be > 0.0 ',[bg%por,bg%k,bg%ss]
       stop 205
    end if
    if (bg%leakFlag < 0 .or. bg%leakFlag > 3) then
       write(*,*) 'error on line',ln,' of input; leak flag (bg%leakFlag) must be in {0,1,2,3}'
       stop 2050
    end if
    if (any([bg%aquitardK,bg%aquitardSs,bg%aquitardb] <= 0.0) .and. bg%leakFlag > 0) then
       write(*,*) 'input file (line ',ln,') bg%aquitardK, bg%aquitardSs, bg%aquitardb &
            &must all be > 0.0 ',[bg%aquitardK,bg%aquitardSs,bg%aquitardb]
       stop 206
    end if

    read(15,*,iostat=ierr) bg%Sy, bg%kz, bg%unconfinedFlag, bg%b; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error on line ',ln,' of input file (unconfined)'
       stop 2060
    end if
    

    if (any([bg%kz,bg%b] <= 0.0) .and. bg%unconfinedFlag) then
       write(*,*) 'input (line ',ln,') vertical K (bg%kz), thickness (bg%b) must be > 0.0 ', &
            & [bg%kz,bg%b]
       stop 207
    end if

    if (bg%Sy < 0.0 .and. bg%unconfinedFlag) then
       write(*,*) 'input (line ',ln,') specific yield (bg%Sy) must be non-negative',bg%Sy
       stop 2071
    end if


    read(15,*,iostat=ierr) bg%matrixSs, bg%lambda, bg%dualPorosityFlag; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (dual porosity) of input'
       stop 2072
    end if
    
    if (any([bg%matrixSs,bg%lambda] <= 0.0) .and. bg%dualPorosityFlag) then
       write(*,*) 'line ',ln,' of input: matrix specific storage (bg%matrixSs) and '//&
            & 'matrix/fracture connection factor (bg%lambda) must be > 0.0',&
            & [bg%matrixSs,bg%lambda]
       stop 2072
    end if

    ! echo input from first 4 lines to file
    write(16,'(5(L1,1X),I0,5(1X,A))') s%calc, s%particle, s%contour, s%timeseries, s%deriv, &
         & s%output, trim(s%outFname), trim(s%coeffFName), trim(s%elemHfName), &
         & trim(s%geomFname),'  ||    re-calculate coefficients?, particle?, '//&
         & 'contour?, timeseries?, deriv?, output, out/coeff/hierarchy/geometry file names'
    write(16,'(3(ES12.5,1X),I0,3(1X,ES12.5),1X,I0,ES11.4,A)') bg%por, bg%k, bg%ss, &
         & bg%leakFlag, bg%aquitardK, bg%aquitardSs, bg%aquitardb, bg%ms, bg%cutoff, &
         & '  ||   background props: por, k, Ss, leaky flag, K2, Ss2, b2, MS, cutoff'
    write(16,'(2(ES12.5,1X),L1,1X,ES12.5,A)') bg%Sy, bg%kz, bg%unconfinedFlag, &
         & bg%b, '  || background props: Sy, Kz, unconfined?, BGb'
    write(16,'(2(ES12.5,1X),L1,A)') bg%matrixSs, bg%lambda, bg%dualPorosityFlag, &
         & '  || background props: matrixSs, matrix/fracture lambda, dual porosity?'

    ! desired sution points/times
    read(15,*,iostat=ierr) s%nx, s%ny, s%nt; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (nx,ny,nt) of input'
       stop 2073
    end if
    
    if (any([s%nx,s%ny,s%nt] < 1)) then
       write(*,*) 'input file (line ',ln,') number of x,y & t observation locations '//&
            & '(s%nx,s%ny,s%nt) must be > 0 ', [s%nx,s%ny,s%nt]
       stop 208
    end if
    if (s%timeseries .and. s%nx /= s%ny) then
       write(*,*) 'for time series output nx==ny.  nx=',s%nx,' ny=',s%ny
       stop 2080
    end if
    allocate(s%x(s%nx), s%y(s%ny), s%t(s%nt))
    if (s%timeseries) then
       allocate(s%obsname(s%nx))
       read(15,'(512A)',iostat=ierr) input; ln=ln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',ln,' (location names) of input'
          stop 2081
       end if
       
       s1 = 1
       slen = len_trim(input) ! don't include trailing blanks
       do j = 1,s%nx
          ! location names are separated by "|" character
          s2 = index(input(s1:),'|')
          if (s2 == 0) then
             ! no "|" separator character found
             if (j == 1 .and. slen > 1) then
                ! first name no separator 
                s%obsname(1) = trim(input)
                s1 = slen
             elseif(j == s%nx .and. s1 > 1) then
                ! last name no separator (normal)
                s%obsname(s%nx) = trim(input(s1:))
             else
                write(chint,'(I4.4)') j
                s%obsname(j) = 'LOC-'//chint ! generic name
             end if
          else
             s%obsname(j) = input(s1:s1+s2-2)
          end if
          s1 = s1+s2
       end do
    else
       allocate(s%obsname(0))
       read(15,*); ln=ln+1 ! read placeholder name line anyway
    end if
    read(15,*,iostat=ierr) s%x(:); ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' x calc locations (s%x) of input'
       stop 2082
    end if
    
    read(15,*,iostat=ierr) s%y(:); ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line',ln,' y calc locations (s%y) of input'
       stop 2083
    end if

    ! shift xy values (usefull when x and y are something like UTM coordinates)
    s%xshift = (maxval(s%x) + minval(s%x))/2.0
    s%yshift = (maxval(s%y) + minval(s%y))/2.0
    s%x(:) = s%x(:) - s%xshift
    s%y(:) = s%y(:) - s%yshift

    read(15,*,iostat=ierr) s%t(:); ln=ln+1
    if (ierr /= 0)  then
       write(*,*) 'error reading line ',ln,' calc times (s%t) of input' 
       stop 2084
    end if

    if (any(s%t <= 0.0)) then
       write(*,*) 'all times (line ',ln,') must be positive values',s%t
       stop 2085
    end if

    if (.not. s%particle) then
       fmt(1) = '(    (ES12.5,1X),A) '
       write(16,'(3(I0,1X),A)') s%nx, s%ny, s%nt, '  ||    numX, numY, numt'
       write(fmt(1)(2:5),'(I4.4)') s%nx
       write(16,'(ES21.14,A)') s%xshift,'  ||    xshift'
       write(16,fmt(1)) s%x(:)+s%xshift,'  ||    original x Vector'
       write(16,fmt(1)) s%x(:),         '  ||    shifted x Vector'
       write(fmt(1)(2:5),'(I4.4)') s%ny
       write(16,'(ES21.14,A)') s%yshift,'  ||    yshift'
       write(16,fmt(1)) s%y(:)+s%yshift, '  ||    original y Vector'
       write(16,fmt(1)) s%y(:),          '  ||    shifted y Vector'
       write(fmt(1)(2:5),'(I4.4)') s%nt
       write(16,fmt(1)) s%t(:), '  ||    t Vector'
    endif

    ! deHoog et al. inverse Laplace transform parameters
    read(15,*,iostat=ierr) s%alpha, s%tol, s%m; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (deHoog invlap) of input'
       stop 209
    end if

    if (s%M < 1) then
       write(*,*) 's%M (line ',ln,') must be > 0 (typically >= 10) ',s%M
       stop 2090
    end if
    if (s%tol < epsilon(s%tol)) then ! epsilon(1) ~ 1.0E-8
       s%tol = epsilon(s%tol)
       write(*,'(A,ES12.5)') 'WARNING: increased deHoog INVLAP solution tolerance to ',s%tol
    end if
    if (s%alpha <= 0.0) then
       write(*,'(A,ES12.5)') 'WARNING: deHoog alpha is typically > 0 ',s%alpha
    end if
    write(16,'(2(ES12.5,1X),I0,A)') s%alpha, s%tol, s%m,'  ||    deHoog: alpha, tol, M'

    ! circular (includes wells)
    read(15,*,iostat=ierr) dom%num(1),circleFname; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (circles) of main input'
       stop 2091
    end if
    
    nc = dom%num(1)
    if (nc > 0) then
       sln = 1
       open(unit=22, file=trim(circleFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening circular data file for reading ',trim(circleFname)
          stop 210
       else
          write(16,'(A)') trim(circleFname)//' opened for circular input data'
       end if

       allocate(c(nc))
       read(22,*,iostat=ierr) c(:)%n
       if (ierr /= 0 .or. any(c%n < 1)) then
          write(*,*) 'error reading line ',sln,' of circle input; # Fourier terms (c%N) '//&
               &'must not be < 1 ',c%N
          stop 211
       end if

       read(22,*,iostat=ierr) c(:)%m; sln=sln+1
       if (ierr /= 0 .or. any(c%M < 1)) then
          write(*,*) 'error reading line ',sln,' of circle input; # matching locations (c%M) '//&
               &'must not be < 1 ',c%M
          stop 212
       end if

       read(22,*,iostat=ierr) c(:)%ibnd; sln=sln+1
       if (ierr /= 0 .or. any(c%ibnd < -1 .or. c%ibnd > 2)) then
          write(*,*) 'error reading line ',sln,' of circle input; boundary type (c%ibnd) '//&
               &'must be in {-1,0,1,2} ',c%ibnd
          stop 213
       end if

       read(22,*,iostat=ierr) c(:)%CalcIn; sln=sln+1 
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' calc inside flag (c%CalcIn) of circle input'
          stop 2130
       end if
       
       read(22,*,iostat=ierr) c(:)%StorIn; sln=sln+1 
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' storage inside flag (c%StorIn) of circle input'
          stop 2131
       end if       

       read(22,*,iostat=ierr) c(:)%r; sln=sln+1
       if (ierr /= 0 .or. any(c%r <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of cirlce input; circle radius (c%r) '//&
               & 'must be > 0.0 ',c%r
          stop 214
       end if

       read(22,*,iostat=ierr) c(:)%x; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' center x (c%x) of circle input'
          stop 2140
       end if       

       read(22,*,iostat=ierr) c(:)%y; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' center y (c%y) of circle input'
          stop 2141
       end if       

       c(:)%z = cmplx(c%x-s%xshift, c%y-s%yshift, DP)

       read(22,*,iostat=ierr) c(:)%k; sln=sln+1
       if (ierr /= 0 .or. any(c%k <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of circle input; hydraulic conductivity '//&
               &'(c%K) must be > 0.0 ',c%k
          stop 215
       end if

       read(22,*,iostat=ierr) c(:)%Ss; sln=sln+1
       if (ierr /= 0 .or. any(c%ss <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of circle input; specific storage (c%Ss) '//&
               &'must be > 0.0 ',c%Ss
          stop 216
       end if

       read(22,*,iostat=ierr) c(:)%por; sln=sln+1
       if (ierr /= 0 .or. any(c%por <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of circle input; porosity (c%por) '//&
               &'must be > 0.0 ',c%por
          stop 217
       end if

       read(22,*,iostat=ierr) c(:)%leakFlag; sln=sln+1
       if (ierr /= 0 .or. any(c%leakFlag < 0) .or. any(c%leakFlag > 3)) then
          write(*,*) 'error reading line ',sln,' leaky flag (c%leakFlag) of circle '//&
               &'input; input must be in {0,1,2,3}'
          stop 2170
       end if       

       read(22,*,iostat=ierr) c(:)%aquitardK; sln=sln+1
       if (ierr /= 0 .or. any(c%aquitardK <= 0.0 .and. c%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of circle input; aquitard vertical K '//&
               &'(c%aquitardK) must be > 0.0 ',c%aquitardk
          stop 218
       end if

       read(22,*,iostat=ierr) c(:)%aquitardSs; sln=sln+1
       if (ierr /= 0 .or. any(c%aquitardSS <= 0.0 .and. c%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of circle input; aquitard specific storage '//&
               &'(c%aquitardSs) must be > 0.0 ',c%aquitardSs
          stop 219
       end if

       read(22,*,iostat=ierr) c(:)%aquitardb; sln=sln+1  
       if (ierr /= 0 .or. any(c%aquitardB <= 0.0 .and. c%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of circle input; aquitard thickness '//&
               &'(c%aquitardB) must be > 0.0 ',c%aquitardB
          stop 220
       end if

       read(22,*,iostat=ierr) c(:)%unconfinedFlag; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' unconfined flag (c%unconfinedFlag) '//&
               &'of circle input'
          stop 2200
       end if       

       read(22,*,iostat=ierr) c(:)%Sy; sln=sln+1
       if (ierr /= 0 .or. any(c%sy <= 0.0 .and. c%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of circle input; specific yield (c%Sy) '//&
               &'must be > 0.0 ',c%sy
          stop 221
       end if

       read(22,*,iostat=ierr) c(:)%Kz; sln=sln+1
       if (ierr /= 0 .or. any(c%kz <= 0.0 .and. c%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of circle input; aquifer vertical K '//&
               &'(c%Kz) must be > 0.0 ',c%kz
          stop 222
       end if

       read(22,*,iostat=ierr) c(:)%b; sln=sln+1
       if (ierr /= 0 .or. any(c%b <= 0.0 .and. c%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of circle input; aquifer thickness '//&
               &'(c%B) must be > 0.0 ',c%b
          stop 223
       end if

       read(22,*,iostat=ierr) c(:)%dualPorosityFlag; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' dual porosity flag (c%dualPorosityFlag) '//&
               &'of circle input'
          stop 2200
       end if       

       read(22,*,iostat=ierr) c(:)%matrixSs; sln=sln+1
       if (ierr /= 0 .or. any(c%matrixSs <= 0.0 .and. c%dualPorosityFlag)) then
          write(*,*) 'error reading line ',sln,' of circle input; matrix specific storage '//&
               &'(c%matrixSs) must be > 0.0', c%matrixSs
          stop 2230
       end if

       read(22,*,iostat=ierr) c(:)%lambda; sln=sln+1
       if (ierr /= 0 .or. any(c%lambda < 0.0 .and. c%dualPorosityFlag)) then
          write(*,*) 'error reading line ',sln,' of circle input; '//&
               & 'matrix/fracture connection (c%lambda) must be non-negative', c%lambda
          stop 2231
       end if

       read(22,*,iostat=ierr) c(:)%dskin; sln=sln+1
       if (ierr /= 0 .or. any(c%dskin <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of circle input; dimensionless skin '//&
               &'parameter (c%Dskin) must be > 0.0 ',c%dskin
          stop 224
       end if

       ! area source strength (flux)
       read(22,*,iostat=ierr) c(:)%areaQ; sln=sln+1 
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' area source (c%areaQ) of circle input'
          stop 2200
       end if       

       ! strength of specified value on boundary (head or flux)
       read(22,*,iostat=ierr) c(:)%bdryQ; sln=sln+1 
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' boundary source (c%bdryQ) of circle input'
          stop 2200
       end if       

       write(chint,'(I4.4)') dom%num(1)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,'(I0,A)') dom%num(1),        '  ||   # of circular elements (including wells)'
       write(16,fmt(1)) c(:)%n,              '  ||   # of circular free parameter (Fourier coefficients)'
       write(16,fmt(1)) c(:)%m,              '  ||   # of circular matching locations'
       write(16,fmt(1)) c(:)%ibnd,           '  ||   circle ibnd array'
       write(16,fmt(2)) c(:)%match,          '  ||   circle matching array'
       write(16,fmt(2)) c(:)%calcin,         '  ||   calculate inside this circle?'
       write(16,fmt(2)) c(:)%storin,         '  ||   calculate free-water storage effects of circle?'
       write(16,fmt(3)) c(:)%r,              '  ||   circle radius'
       write(16,fmt(3)) c(:)%x+s%xshift,     '  ||   original circle center x'
       write(16,fmt(3)) c(:)%x,              '  ||   shifted circle center x'
       write(16,fmt(3)) c(:)%y+s%yshift,     '  ||   original circle center y'
       write(16,fmt(3)) c(:)%y,              '  ||   shifted circle center y'
       write(16,fmt(3)) c(:)%k,              '  ||   circle aquifer k'
       write(16,fmt(3)) c(:)%ss,             '  ||   circle aquifer Ss'
       write(16,fmt(3)) c(:)%por,            '  ||   circle aquifer porosity'
       write(16,fmt(1)) c(:)%leakFlag,       '  ||   circle leaky type'
       write(16,fmt(3)) c(:)%aquitardK,      '  ||   circle leaky aquitard K'
       write(16,fmt(3)) c(:)%aquitardSs,     '  ||   circle leaky aquitard Ss'
       write(16,fmt(3)) c(:)%aquitardb,      '  ||   circle leaky aquitard thickness'
       write(16,fmt(2)) c(:)%unconfinedFlag, '  ||   circle unconfined flag'
       write(16,fmt(3)) c(:)%Sy,             '  ||   circle aquifer specific yield'
       write(16,fmt(3)) c(:)%Kz,             '  ||   circle aquifer vertical K'
       write(16,fmt(3)) c(:)%b,              '  ||   circle aquifer thickness'
       write(16,fmt(2)) c(:)%dualPorosityFlag,'  ||   circle dual porosity flag'
       write(16,fmt(3)) c(:)%matrixSs,        '  ||   circle matrix Ss'
       write(16,fmt(3)) c(:)%lambda,          '  ||   circle matrix/fracture connection lambda'
       write(16,fmt(3)) c(:)%dskin,           '  ||   circle boundary dimensionless skin factor'
       write(16,fmt(3)) c(:)%areaQ,          '  ||   circle area rch rate'
       write(16,fmt(3)) c(:)%bdryQ,          '  ||   circle boundary rch rate or head'
       
       ! area source behavior
       do j = 1,size(c,dim=1)
          read(22,*,iostat=ierr) c(j)%AreaTime; sln=sln+1
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area time behavior (c%AreaTime) &
                     &of circle ',j,' input'
                stop 2201
             end if
          backspace(22)
          if (c(j)%AreaTime > -1) then
             ! functional time behavior
             allocate(c(j)%ATPar(2))
             read(22,*,iostat=ierr) c(j)%AreaTime,c(j)%ATPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area functional time behavior '//&
               &'(c%ATPar) of circle ',j,' input'
                stop 2202
             end if
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%AreaTime,c(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for circle ',j
          else
             ! piecewise-constant time behavior
             allocate(c(j)%ATPar(-2*c(j)%AreaTime+1))
             read(22,*,iostat=ierr) c(j)%AreaTime,c(j)%ATPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area piecewise-constant time '//&
               &'behavior (c%ATpar) of circle ',j,'input'
                stop 2203
             end if
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)' 
             write(lfmt(8:11),'(I4.4)')  size(c(j)%ATPar(:-c(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%ATPar(-c(j)%AreaTime+2:),1)
             write(16,lfmt) c(j)%AreaTime,c(j)%ATPar(:-c(j)%AreaTime+1),'| ',&
                  & c(j)%ATPar(-c(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for circle ',j
          end if
       end do

       ! boundary source behavior
       do j = 1,size(c,dim=1)
          read(22,*,iostat=ierr) c(j)%BdryTime; sln=sln+1
          if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary time behavior '//&
                     &'(c%BdryTime) of circle ',j,' input'
                stop 2204
             end if
          backspace(22)
          if (c(j)%BdryTime > -1) then
             allocate(c(j)%BTPar(2))
             read(22,*,iostat=ierr) c(j)%BdryTime,c(j)%BTPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary functional time '//&
                     &'behavior (c%BTPar) of circle ',j,' input'
                stop 2205
             end if
             write(16,'(I0,2(1X,ES12.5),A,I0)') c(j)%BdryTime,c(j)%BTPar(:),&
                  &'  ||  Bdry time behavior, par1, par2 for circle ',j
          else
             allocate(c(j)%BTPar(-2*c(j)%BdryTime+1))
             read(22,*,iostat=ierr) c(j)%BdryTime,c(j)%BTPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary piecewise-constant time '//&
                     &'behavior (c%BTpar) of circle ',j,'input'
                stop 2203
             end if
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)' 
             write(lfmt(8:11),'(I4.4)')  size(c(j)%BTPar(:-c(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I4.4)') size(c(j)%BTPar(-c(j)%BdryTime+2:),1)
             write(16,lfmt) c(j)%BdryTime,c(j)%BTPar(:-c(j)%BdryTime+1),' | ',&
                  & c(j)%BTPar(-c(j)%BdryTime+2:), &
                  &'  ||    Bdry ti, tf | strength for circle ',j
          end if
       end do

       close(22) ! circle input file

       where (c(:)%ibnd == -1 .or. c(:)%ibnd == 0 .or. c(:)%ibnd == +1)
          c(:)%match = .true.
       elsewhere !( currently just 2)
          c(:)%match = .false.
       end where

       if(any(c%match .and. c%storin)) then
          write(*,*) 'WARNING: wellbore (free-water) storage only implemented for ibnd==2'
          where(c%match .and. c%storin)
             c%storin = .false.
          end where
       end if

       ! minor checking / correcting
       do j = 1,size(c,dim=1)
          if (c(j)%ibnd == 2 .and. c(j)%N /= 1) then
             write(*,'(A,I0,A)') '** wells (ibnd==2) must have N=1, fixing circle #',&
                  & j,' to N=1'
             c(j)%N = 1
          end if
       end do

    else
       ! no circular elements
       allocate(c(0))
    end if

    ! elliptical (includes line sources/sinks)
    read(15,*,iostat=ierr) dom%num(2), ellipseFname; ln=ln+1
    if (ierr /= 0) then
       write(*,*) 'error reading line ',ln,' (ellipses) of main input'
       stop 2092
    end if

    ne = dom%num(2)
    sln = 1
    if (ne > 0) then
       open(unit=33, file=trim(ellipseFname), status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening elliptical data file for reading ',trim(ellipseFname)
          stop 225
       else
          write(16,'(A)') trim(ellipseFname)//' opened for elliptical input data'
       end if

       allocate(e(ne))
       read(33,*,iostat=ierr) e(:)%n; sln=sln+1
       if (ierr /= 0 .or. any(e%n < 1)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; # Fourier terms (e%N) '//&
                     &'must not be < 1 ',e%N
          stop 226
       end if

       read(33,*,iostat=ierr) e(:)%m; sln=sln+1
       if (ierr /= 0 .or. any(e%m < 1)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; # matching locations '//&
                     &'(e%M) must not be < 1 ',e%M
          stop 227
       end if

       read(33,*,iostat=ierr) e(:)%ms; sln=sln+1 ! bg%ms read above
       ! these are computed now (this value is now used as an alternate minimum),
       ! so checking this here doesn't mean as much as it once did.
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' elliptical eigenvector matrix size '//&
                     &'(e%MS) of ellipse input'
          stop 2270
       end if       

       read(33,*,iostat=ierr) e(:)%cutoff; sln=sln+1
       if (ierr /= 0 .or. any(e%cutoff <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; matrix off-diagonal '//&
                     &'cutoff (e%cutoff) must be > 0.0 ',e%cutoff
          stop 228
       end if

       read(33,*,iostat=ierr) e(:)%ibnd; sln=sln+1
       if (ierr /= 0 .or. any(e%ibnd < -1 .or. e%ibnd > 2)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; boundary type (e%ibnd) '//&
                     &'must be in {-1,0,1,2} ',e%ibnd
          stop 229
       end if

       read(33,*,iostat=ierr) e(:)%CalcIn; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' calc inside flag (e%CalcIn) of '//&
                     &'ellipse input'
          stop 2290
       end if       

       read(33,*,iostat=ierr) e(:)%StorIn; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' storage inside flag (e%StorIn) of '//&
                     &'ellipse input'
          stop 2291
       end if       

       read(33,*,iostat=ierr) e(:)%r; sln=sln+1   ! eta
       if (ierr /= 0 .or. any(e%r < 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; elliptical radius '//&
                     &'(e%r) must be >= 0.0 ',e%r
          stop 230
       end if

       read(33,*,iostat=ierr) e(:)%x; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' center x (e%x) of ellipse input'
          stop 2300
       end if       

       read(33,*,iostat=ierr) e(:)%y; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' center y (e%y) of ellipse input'
          stop 2301
       end if       

       e(:)%z = cmplx(e%x-s%xshift, e%y-s%yshift, DP)

       read(33,*,iostat=ierr) e(:)%f; sln=sln+1
       if (ierr /= 0 .or. any(e%f <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; elliptical semi-focal '//&
                     &'length (e%f) must be > 0.0 ',e%f
          stop 231
       end if

       read(33,*,iostat=ierr) e(:)%theta; sln=sln+1
       if (ierr /= 0 .or. any(e%theta < -PI) .or. any(e%theta > PI)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; eliptical angle wrt '//&
                     &'global cartesian (e%theta) must be -pi <= theta <= PI ',e%theta
          stop 232
       end if

       read(33,*,iostat=ierr) e(:)%k; sln=sln+1
       if (ierr /= 0 .or. any(e%k <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; hydraulic conductivity '//&
                     &'(e%k) must be > 0.0 ',e%k
          stop 233
       end if

       read(33,*,iostat=ierr) e(:)%Ss; sln=sln+1
       if (ierr /= 0 .or. any(e%Ss <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; specific storage '//&
                     &'(e%Ss) must be > 0.0 ',e%Ss
          stop 234
       end if

       read(33,*,iostat=ierr) e(:)%por; sln=sln+1
       if (ierr /= 0 .or. any(e%por <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; porosity (e%por) '//&
                     &'must be > 0.0 ',e%por
          stop 235
       end if

       read(33,*,iostat=ierr) e(:)%leakFlag; sln=sln+1  ! checking handled elsewhere
       if (ierr /= 0 .or. any(e%leakFlag < 0) .or. any(e%leakFlag > 3)) then
          write(*,*) 'error reading line ',sln,' leak flag (e%leakFlag) of ellipse '//&
                     &'input; input must be in {0,1,2,3}'
          stop 2350
       end if       

       read(33,*,iostat=ierr) e(:)%aquitardK; sln=sln+1
       if (ierr /= 0 .or. any(e%aquitardK <= 0.0 .and. e%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; aquitard vertical K '//&
                     &'(e%aquitardK) must be > 0.0 ',e%aquitardK
          stop 236
       end if

       read(33,*,iostat=ierr) e(:)%aquitardSs; sln=sln+1
       if (ierr /= 0 .or. any(e%aquitardSs <= 0.0 .and. e%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; aquitard specific '//&
                     &'storage (e%aquitardSs) must be > 0.0 ',e%aquitardSs
          stop 237
       end if

       read(33,*,iostat=ierr) e(:)%aquitardb; sln=sln+1
       if (ierr /= 0 .or. any(e%aquitardb <= 0.0 .and. e%leakFlag > 0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; aquitard thickness '//&
                     &'(e%aquitardB) must be > 0.0 ',e%aquitardb
          stop 238
       end if

       read(33,*,iostat=ierr) e(:)%unconfinedFlag; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' unconfined flag (e%unconfinedFlag) '//&
                     &'of ellipse input'
          stop 2380
       end if       

       read(33,*,iostat=ierr) e(:)%Sy; sln=sln+1
       if (ierr /= 0 .or. any(e%Sy <= 0.0 .and. e%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; specific storage '//&
                     &'(e%Sy) must be > 0.0 ',e%Sy
          stop 239
       end if

       read(33,*,iostat=ierr) e(:)%Kz; sln=sln+1
       if (ierr /= 0 .or. any(e%Kz <= 0.0 .and. e%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; aquifer vertical K '//&
                     &'(e%Kz) must be > 0.0 ',e%Kz
          stop 240
       end if

       read(33,*,iostat=ierr) e(:)%b; sln=sln+1
       if (ierr /= 0 .or. any(e%b <= 0.0 .and. e%unconfinedFlag)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; aquifer thickness '//&
                     &'(e%b) must be > 0.0 ',e%b
          stop 241
       end if

       read(33,*,iostat=ierr) e(:)%dualPorosityFlag; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' dual porosity flag '//&
                     &'(e%dualPorosityFlag) of ellipse input'
          stop 2410
       end if       

       read(33,*,iostat=ierr) e(:)%matrixSs; sln=sln+1
       if (ierr /= 0 .or. any(e%matrixSs <= 0.0 .and. e%dualPorosityFlag)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; matrix specific '//&
                     &'storage (e%matrixSs) must be > 0.0', e%matrixSs
          stop 2411
       end if

       read(33,*,iostat=ierr) e(:)%lambda; sln=sln+1
       if (ierr /= 0 .or. any(e%lambda < 0.0 .and. e%dualPorosityFlag)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; matrix/fracture '//&
                     &'connection (e%lambda)'//&
               &' must be non-negative', e%lambda
          stop 2412
       end if

       read(33,*,iostat=ierr) e(:)%dskin; sln=sln+1
       if (ierr /= 0 .or. any(e%dskin <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of ellipse input; dimensionless '//&
                     &'skin (e%dskin) must be > 0.0 ',e%dskin
          stop 242
       end if

       read(33,*,iostat=ierr) e(:)%areaQ; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' area source (e%areaQ) of ellipse input'
          stop 2420
       end if       

       read(33,*,iostat=ierr) e(:)%bdryQ; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' boundary source (e%bdryQ) of ellipse input'
          stop 2421
       end if       

       write(chint,'(I4.4)') dom%num(2)
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,'(I0,A)') dom%num(2),        '  ||   # of elliptical elements (including lines)'
       write(16,fmt(1)) e(:)%n,              '  ||   # of elliptical free parameter (Fourier coeffs)'
       write(16,fmt(1)) e(:)%m,              '  ||   # of ellipse matching locations'
       write(16,fmt(1)) e(:)%ms,             '  ||   size of "infinite" Mathieu matrices'
       write(16,fmt(3)) e(:)%cutoff,         '  ||   desired accuracy cutoff for Mathieu function matrices'
       write(16,fmt(1)) e(:)%ibnd,           '  ||   ellipse ibnd array'
       write(16,fmt(2)) e(:)%match,          '  ||   ellipse matching array'
       write(16,fmt(2)) e(:)%calcin,         '  ||   calculate inside this ellipse?'
       write(16,fmt(2)) e(:)%storin,         '  ||   calculate free-water storage for this ellipse?'
       write(16,fmt(3)) e(:)%r,              '  ||   ellipse radius (eta)'
       write(16,fmt(3)) e(:)%x+s%xshift,     '  ||   original ellipse center x'
       write(16,fmt(3)) e(:)%x,              '  ||   shifted ellipse center x'
       write(16,fmt(3)) e(:)%y+s%yshift,     '  ||   original ellipse center y'
       write(16,fmt(3)) e(:)%y,              '  ||   shifted ellipse center y'
       write(16,fmt(3)) e(:)%f,              '  ||   ellipse semi-focal length'
       write(16,fmt(3)) e(:)%theta,          '  ||   ellipse angle rotation with +x axis'
       write(16,fmt(3)) e(:)%k,              '  ||   ellipse aquifer k'
       write(16,fmt(3)) e(:)%ss,             '  ||   ellipse aquifer Ss'
       write(16,fmt(3)) e(:)%por,            '  ||   ellipse aquifer porosity'
       write(16,fmt(1)) e(:)%leakFlag,       '  ||   ellipse leaky type'
       write(16,fmt(3)) e(:)%aquitardK,      '  ||   ellipse leaky aquitard K'
       write(16,fmt(3)) e(:)%aquitardSs,     '  ||   ellipse leaky aquitard Ss'
       write(16,fmt(3)) e(:)%aquitardb,      '  ||   ellipse leaky aquitard thickness'
       write(16,fmt(2)) e(:)%unconfinedFlag, '  ||   ellipse unconfined flag'
       write(16,fmt(3)) e(:)%Sy,             '  ||   ellipse aquifer specific yield'
       write(16,fmt(3)) e(:)%Kz,             '  ||   ellipse aquifer vertical K'
       write(16,fmt(3)) e(:)%b,              '  ||   ellipse aquifer thickness'
       write(16,fmt(2)) e(:)%dualPorosityFlag,'  ||     ellipse dual porosity flag'
       write(16,fmt(3)) e(:)%matrixSs,        '  ||     ellipse matrix Ss'
       write(16,fmt(3)) e(:)%lambda,          '  ||     ellipse matrix/fracture connection lambda'
       write(16,fmt(3)) e(:)%dskin,           '  ||   ellipse boundary dimensionless skin factor'
       write(16,fmt(3)) e(:)%areaQ,          '  ||   ellipse area rch rate'
       write(16,fmt(3)) e(:)%bdryQ,          '  ||   ellipse boundary rch rate or head'

       do j = 1,size(e,dim=1)
          read(33,*,iostat=ierr) e(j)%AreaTime; sln=sln+1
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area time behavior (e%AreaTime) &
                     &of ellipse ',j,' input'
                stop 2422
             end if
          backspace(33)
          if (e(j)%AreaTime > -1) then
             allocate(e(j)%ATPar(2))
             read(33,*,iostat=ierr) e(j)%AreaTime,e(j)%ATPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area functional time behavior '//&
                     &'(e%ATPar) of ellipse ',j,' input'
                stop 2423
             end if
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%AreaTime,e(j)%ATPar(:),&
                  &'  ||  Area time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%ATPar(-2*e(j)%AreaTime+1))
             read(33,*,iostat=ierr) e(j)%AreaTime,e(j)%ATPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' area piecewise-constant time '//&
                     &'behavior (e%ATPar) of circle ',j,'input'
                stop 2424
             end if
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)' 
             write(lfmt(8:11),'(I3.3)')  size(e(j)%ATPar(:-e(j)%AreaTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%ATPar(-e(j)%AreaTime+2:),1)
             write(16,lfmt) e(j)%AreaTime,e(j)%ATPar(:-e(j)%AreaTime+1),' | ',&
                  & e(j)%ATPar(-e(j)%AreaTime+2:), &
                  &'  ||    Area ti, tf | strength for ellipse ',j
          end if
       end do
       do j = 1,size(e,dim=1)
          read(33,*,iostat=ierr) e(j)%BdryTime; sln=sln+1
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary time behavior '//&
                     &'(e%BdryTime) of ellipse ',j,' input'
                stop 2425
             end if
          backspace(33)
          if (e(j)%BdryTime > -1) then
             allocate(e(j)%BTPar(2))
             read(33,*,iostat=ierr) e(j)%BdryTime,e(j)%BTPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary functional time '//&
                     &'behavior (e%BTPar) of ellipse ',j,' input'
                stop 2426
             end if
             write(16,'(I0,2(1X,ES12.5),A,I0)') e(j)%BdryTime,e(j)%BTPar(:),&
                  &'  ||  Boundary time behavior, par1, par2 for ellipse ',j
          else
             allocate(e(j)%BTPar(-2*e(j)%BdryTime+1))
             read(33,*,iostat=ierr) e(j)%BdryTime,e(j)%BTPar(:)
             if (ierr /= 0) then
                write(*,*) 'error reading line ',sln,' boundary piecewise-constant '//&
                     &'time behavior (e%BTPar) of circle ',j,'input'
                stop 2427
             end if
             lfmt = '(I0,1X,    (ES12.5,1X),A,    (ES12.5,1X),A,I0)' 
             write(lfmt(8:11),'(I3.3)')  size(e(j)%BTPar(:-e(j)%BdryTime+1),1)
             write(lfmt(26:29),'(I3.3)') size(e(j)%BTPar(-e(j)%BdryTime+2:),1)
             write(16,lfmt) e(j)%BdryTime,e(j)%BTPar(:-e(j)%BdryTime+1),' | ',&
                  & e(j)%BTPar(-e(j)%BdryTime+2:), &
                  &'  ||    Boundary ti, tf | strength for ellipse ',j
          end if
       end do

       close(33) ! ellipse input file

       where (e(:)%ibnd == -1 .or. e(:)%ibnd == 0 .or. e(:)%ibnd == +1)
          e(:)%match = .true.
       elsewhere
          e(:)%match = .false.
       end where
       
       ! TODO: handle free-water storage for ellipses?
       if(any(e%storin)) then
          write(*,*) 'WARNING: wellbore (free-water) storage not handled for ellipses yet'
          e%storin = .false.
       end if

    else
       ! no elliptical elements
       allocate(e(0))
    end if

    ntot = sum(dom%num) ! total number of circular and elliptical elements
    if (ntot < 1) stop 'READINPUT: Need at least one circular (including well) or'//&
         &'elliptical (including line) element.'

    ! compute secondary parameters
    bg%alpha = bg%K/bg%Ss
    bg%T = bg%K*bg%b
    c(:)%alpha = c(:)%K/c(:)%Ss
    e(:)%alpha = e(:)%K/e(:)%Ss
    c(:)%T = c(:)%K*c(:)%b
    e(:)%T = e(:)%K*e(:)%b

    write(chint,'(I4.4)') dom%num(1)
    fmt(2) = '('//chint//'(ES11.5,1X),A) ' ! circles
    write(chint,'(I4.4)') dom%num(2)
    fmt(3) = '('//chint//'(ES11.5,1X),A) ' ! ellipses

    if (dom%num(1) > 0) then
       write(16,fmt(2)) c(:)%alpha,'  ||    circle hydraulic diffusivity'
       write(16,fmt(2)) c(:)%T,    '  ||    circle transmissivity'
    end if
    if (dom%num(2) > 0) then
       write(16,fmt(3)) e(:)%alpha,'  ||    ellipse hydraulic diffusivity'
       write(16,fmt(3)) e(:)%T,    '  ||    ellipse transmissivity'
    end if
    write(16,'(ES12.5,A)') bg%alpha,'  ||    background hydraulic diffusivity'
    write(16,'(ES12.5,A)') bg%T,    '  ||    background transmissivity'

    ! particles
    if (s%particle) then
       sln = 1
       read(15,*,iostat=ierr) particleFname; ln=ln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',ln,' (particles) of main input'
          stop 2500
       end if
       open(unit=44, file=particleFname, status='old', action='read',iostat=ierr)
       if (ierr /= 0) then
          write(*,'(2A)') 'READINPUT: error opening particle data file for reading ',&
               & particleFname
          stop 243
       else
          write(16,'(A)') trim(particleFname)//' opened for particle input data'
       end if

       read(44,*,iostat=ierr) s%nPart,  s%streakSkip
       if (ierr /= 0 .or. (s%nPart < 1) .or. (s%streakSkip < 0 .and. s%output == 21)) then
          write(*,*) 'error reading line ',sln,' of particle input; s%nPart and &
               &s%streakSkip must be >0',[s%nPart,s%streakSkip]
          stop 2430
       end if
       allocate(p(s%nPart))

       read(44,*,iostat=ierr) p(:)%forward; sln=sln+1 
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' forward/backward particle tracking '//&
                     &'flag (p%forward) of particle input'
          stop 2431
       end if

       read(44,*,iostat=ierr) p(:)%tol; sln=sln+1
       if (ierr /= 0 .or. any(p%tol < spacing(1.0D0))) then
          write(*,*) 'error reading line ',sln,' of particle input; RKM error '//&
                     &'tolerance (p%tol) must be > 0.0 ',p%tol
          stop 244
       end if

       read(44,*,iostat=ierr) p(:)%maxL; sln=sln+1  
       if (ierr /= 0 .or. any(p%maxL < spacing(1.0))) then
          write(*,*) 'error reading line ',sln,' of particle input; max RKM '//&
                     &'step length (p%maxL) must be > 0.0 ',p%maxL
          stop 245
       end if

       read(44,*,iostat=ierr) p(:)%mindt; sln=sln+1
       if (ierr /= 0 .or. any(p%mindt < spacing(1.0))) then
          write(*,*) 'error reading line ',sln,' of particle input; min RKM '//&
                     &'time step (p%mindt) must be > 0.0 ',p%mindt
          stop 246
       end if

       read(44,*,iostat=ierr) p(:)%dt; sln=sln+1    
       if (ierr /= 0 .or. any(p%dt < spacing(1.0))) then
          write(*,*) 'error reading line ',sln,' of particle input; initial '//&
                     &'time step (p%dt) must be > 0.0 ',p%dt
          stop 247
       end if

       read(44,*,iostat=ierr) p(:)%x; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' start x (p%x) of particle input'
          stop 2470
       end if

       read(44,*,iostat=ierr) p(:)%y; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' start y (p%y) of particle input'
          stop 2471
       end if

       p(:)%x = p(:)%x - s%xshift
       p(:)%y = p(:)%y - s%yshift

       read(44,*,iostat=ierr) p(:)%ti; sln=sln+1
       if (ierr /= 0 .or. any(p%ti <= 0.0)) then
          write(*,*) 'error reading line ',sln,' of particle input; start '//&
                     &'time (p%ti) must be > 0.0 ',p%ti
          stop 248
       end if

       read(44,*,iostat=ierr) p(:)%tf; sln=sln+1
       if (ierr /= 0 .or. any(p%tf < spacing(1.0))) then
          write(*,*) 'error reading line ',sln,' of particle input; end '//&
                     &'time (p%tf) must be > 0.0 ',p%tf
          stop 249
       end if

       if (any(p%tf < p%ti .and. p%forward)) then
          write(*,*) 'p%tf must be > p%ti  for forward tracking ti:',p%ti,' tf:',p%tf
          stop 250
       end if
       if (any(p%tf > p%ti .and. .not.p%forward)) then
          write(*,*) 'p%tf must be < p%ti  for backward tracking ti:',p%ti,' tf:',p%tf
          stop 251
       end if

       read(44,*,iostat=ierr) p(:)%int; sln=sln+1
       if (ierr /= 0 .or. any(p%int < 1 .or. p%int > 4)) then
          write(*,*) 'error reading line ',sln,' of particle input; integration '//&
                     &'method (p%int) must be {1,2,3,4} ',p%int
          stop 252
       end if

       read(44,*,iostat=ierr) p(:)%InclIn; sln=sln+1
       if (ierr /= 0) then
          write(*,*) 'error reading line ',sln,' begin inside element flag '//&
                     &'(p%InclIn) of particle input'
          stop 2520
       end if

       close(44) ! particle input file

       write(chint,'(I4.4)') s%nPart
       fmt(1) = '('//chint//'(I0,1X),A)     ' ! integer
       fmt(2) = '('//chint//'(L1,1X),A)     ' ! logical
       fmt(3) = '('//chint//'(ES13.5,1X),A) ' ! real

       write(16,fmt(3)) p(:)%tol,    '  ||    particle sution tolerances (RKM only)'
       write(16,fmt(3)) p(:)%dt,     '  ||    particle time step (RKM initial step)'
       write(16,fmt(3)) p(:)%maxL,   '  ||    particle max step length (RKM only)'
       write(16,fmt(3)) p(:)%mindt,  '  ||    particle min time step size (RKM only)'
       write(16,fmt(3)) p(:)%x+s%xshift,  '  ||    original particle initial x'
       write(16,fmt(3)) p(:)%x,           '  ||    shifted particle initial x'
       write(16,fmt(3)) p(:)%y+s%yshift,  '  ||    original particle initial y'
       write(16,fmt(3)) p(:)%y,           '  ||    shifted particle initial y'
       write(16,fmt(3)) p(:)%ti,     '  ||    particle initial t'
       write(16,fmt(3)) p(:)%tf,     '  ||    particle maximum t'
       write(16,fmt(1)) p(:)%int,    '  ||    particle integration method'
       write(16,fmt(2)) p(:)%InclIn, '  ||    particle begins inside CH/CF incl?'
    else
       ! no particles
       allocate(p(0))
    endif
    close(15) ! main input file
    close(16) ! input echo file
  end subroutine readInput

  !******************************************************
  subroutine writeResults(s,p)

    use type_definitions, only : solution, particle

    type(solution), intent(inout) :: s
    type(particle), dimension(:), intent(inout) :: p

    character(4), dimension(2) :: chint
    integer :: i, j, k, nt

    ! adjust the formats of time, location, and results here
    character(6) :: tfmt = 'ES13.5', xfmt = 'ES12.4'
    character(9) :: hfmt = 'ES22.14e3'

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
       write(20,*) '# LT-AEM contour map output     -*-auto-revert-*-'
       write(20,'(A,I0)') '# t: ', s%nt
       write(20,'(A,I0)') '# x: ', s%nx
       write(20,'(A,I0)') '# y: ', s%ny
       write(20,'(A,I0)') '# locations:', s%nx*s%ny

       do i = 1, s%nt
          write(20,'(A,'//tfmt//')') ' # t= ',s%t(i)
          write(20,'(A)',advance='no')   &
          & '#      X           Y               head'//&
          & '                velx                  vely'
          if (s%deriv) then
             write(20,'(A)') '                d(head)/d(log(t))'
          else
             write(20,'(A)') ''
          end if
          do j = 1, s%ny
             do k = 1, s%nx
                if (s%deriv) then
                   write(20,'(2('//xfmt//',1X),4('//hfmt//',1X))') &
                        & s%x(k), s%y(j), s%h(k,j,i), s%v(k,j,i,1:2), s%dh(k,j,i)
                else
                   write(20,'(2('//xfmt//',1X),3('//hfmt//',1X))') &
                        & s%x(k), s%y(j), s%h(k,j,i), s%v(k,j,i,1:2)
                end if
             end do
          end do
          write(20,'(/)')
       end do
       write(20,'(A)') '# EOF'
       close(20)

       write(*,'(/A)') '**************************************************'
       write(*,'(2A)') ' Gnuplot contour output => ', trim(s%outfname)
       write(*,'(A)')  '**************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (2)

       ! ** Matlab-friendly contour output **
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
             write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%h(j,i,k), j=1,s%nx)
          end do
          close(20)

          ! log-t derivative matrix
          if (s%deriv) then
             open(unit=20, file=trim(s%outfname)//'_dhead_'//chint(2)//'.dat', &
                  & status='replace', action='write')
             do i = 1, s%ny
                write (20,'('//chint(1)//'(1x,'//hfmt//'))') (s%dh(j,i,k), j=1,s%nx)
             end do
             close(20)
          end if

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

       write(*,'(/A)') '*********************************************************'
       write(*,'(3A)') 'Matlab output => ', trim(s%outfname), &
              & '{x,y,t,{d,}head{1-n},velx{1-n},vely{1-n}}.dat'
       write(*,'(A)') '**********************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (10)

       ! ** Gnuplot-friendly time series output **
       ! column of time values at a location through time
       ! locations separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# LT-AEM time series output   -*-auto-revert-*-'
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

          do k = 1, s%nt
             if (s%deriv) then
                write (20,'(1X,'//tfmt//',4(1X,'//hfmt//'))') &
                     & s%t(k),s%h(i,1,k),s%v(i,1,k,1:2),s%dh(i,1,k)
             else
                write (20,'(1X,'//tfmt//',3(1X,'//hfmt//'))') &
                     & s%t(k),s%h(i,1,k),s%v(i,1,k,1:2)
             end if
          end do
          write (20,'(/)')
       end do
       write(20,*) '# EOF'
       close(20)

       write(*,'(/A)') '*****************************************************'
       write(*,'(2A)') 'Gnuplot timeseries output => ', trim(s%outfname)
       write(*,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (11)

       ! ** Gnuplot-friendly time series output **
       ! column of time values at a location through time
       ! locations separated by blank lines (no velocity)
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# LT-AEM time series output    -*-auto-revert-*-'
       do j = 1, s%nx
          write (20,'(2(A,'//xfmt//'),3X,A)') '# location: x=',s%x(j),' y=',&
               &s%y(j),trim(s%obsname(j))
          if (s%deriv) then
             write (20,'(A)')   '#     time              head             deriv'
          else
             write (20,'(A)')   '#     time              head'
          end if
          do k = 1, s%nt
             if (s%deriv) then
                write (20,'(1X,'//tfmt//',1(1X,'//hfmt//'))') &
                     & s%t(k),s%h(j,1,k),s%dh(j,1,k)
             else
                write (20,'(1X,'//tfmt//',2(1X,'//hfmt//'))') &
                     & s%t(k),s%h(j,1,k)
             end if
          end do
          write (20,'(/)')
       end do
       write(20,*) '# EOF'
       close(20)

       write(*,'(/A)') '*****************************************************'
       write(*,'(2A)') 'Gnuplot timeseries output (no v) => ', trim(s%outfname)
       write(*,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (12)

       ! ** output for inverse in Matlab
       ! column of time values at a location through time
       ! locations separated by blank lines

       do i = 1, s%nx
          write(chint(1),'(I4.4)') i
          open(unit=20, file=trim(s%outfname)//'_'//chint(1), status='replace', &
               &action='write')
          do k = 1, s%nt
             if (s%deriv) then
                write (20,'('//tfmt//',2(1X,'//hfmt//'))') s%t(k),s%h(i,1,k),s%dh(i,1,k)
             else
                write (20,'('//tfmt//',1X,'//hfmt//')') s%t(k),s%h(i,1,k)
             end if
          end do
          write(20,'(/)')
          close(20)
       end do

       write(*,'(/A)') '*****************************************************'
       write(*,'(4A)') 'inverse output => ',trim(s%outfname),'0000-',chint(1)
       write(*,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (20)

       ! ** pathline Gnuplot-style output **
       ! columns of time values for starting locations
       ! particles separated by blank lines
       open(unit=20, file=s%outfname, status='replace', action='write')
       write (20,'(A)') '# LT-AEM particle tracking output  -*-auto-revert-*-'
       do i = 1, size(p,dim=1)
          write (20,'(A)')   &
          & '#     time              x                    y                  '//&
          &'velx                 vely '
          do k = 1,p(i)%numt
             write (20,'('//tfmt//',4(1X,'//hfmt//'))') &
                  & p(i)%r(k,1:5)
          end do
          write (20,'(/)')
       end do
       write(20,'(A)') '# EOF'
       close(20)

       write(*,'(/A)') '*****************************************************'
       write(*,'(2A)') 'particle tracking output => ', trim(s%outfname)
       write(*,'(A)')  '*****************************************************'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case (21)

       ! ** streakline Gnuplot-style output **
       ! this requires constant time steps (can't use adaptive integration)
       ! each block is a requested time, each row a particle

       open(unit=90, file=s%outfname, status='replace', action='write')
       write (90,'(A)') '# ltaem particle tracking streakfile output   -*-auto-revert-*-'

       ! max number of times for all particles
       nt = maxval(p(:)%numt,dim=1)

       do i = 1, nt, s%streakSkip
          ! use maxval to ensure a non-zero time is reported
          write (90,'(A'//tfmt//')') '# time:', maxval(p(i)%r(:,1))
          write (90,'(A)') '#  particle       x            y&
               &           velx         vely'
          do j = 1, size(p,dim=1)
             if (p(j)%r(i,1) > 0.0) then
                ! only write particle if it has non-zero data
                write (90,'(I0,4(1X'//hfmt//'))')  j,p(j)%r(i,2:5)
             end if
          end do
          write (90,'(/)')
       end do
       write(90,'(A)') '# EOF'
       close(90)

       write(*,'(/A)') '*****************************************************'
       write(*,'(2A)') 'particle tracking streakfile => ', trim(s%outfname)
       write(*,'(A)')  '*****************************************************'

    case default
       write(*,'(A,I0)')  'invalid output code ',s%output
       stop 300
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
    complex(DP) :: origin
    character(14) :: fmt

    ! remove shift originally applied to coordinates
    origin = cmplx(s%xshift,s%yshift,DP)

    nc = size(c,dim=1)
    ne = size(e,dim=1)
    fmt = '(2(ES13.5,1X))'

    ! write matching points to file
    open(unit=40, file=s%geomfname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       ! non-fatal error
       write(*,'(2A)') 'WARNING: writeGeometry error opening output file for writing &
            &element matching locations ',s%geomFname
    else
       write(40,'(A)') '# points along circumference of circular '//&
            &'and elliptical elements  -*-auto-revert-*-'
       do i = 1,nc
          write(40,'(2(A,I0),A)') '# circular element ',i,' = ',c(i)%M,' points'
          write(40,fmt)  (origin + c(i)%Zom(j), j=1,c(i)%M)
          if (c(i)%M > 1) then
             ! joins circle back up with beginning for plotting
             write(40,fmt)  origin + c(i)%Zom(1)
          end if
          write(40,'(/)')
       end do

       do i = 1,ne
          write(40,'(2(A,I0),A)') '# elliptical element ',i,' = ',e(i)%M,' points'
          write(40,fmt)  (origin + e(i)%Zom(j), j=1,e(i)%M)
          if (e(i)%M > 1) then
             ! joins ellipse back up with beginning for plotting
             write(40,fmt)  origin + e(i)%Zom(1)
          end if
          write(40,'(/)')
       end do
       write(40,'(A)') '# EOF'
       close(40)
    end if

  end subroutine writeGeometry

end module file_ops

