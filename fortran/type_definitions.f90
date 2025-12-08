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

module type_definitions
  use constants, only : DP, lenFN, ASCII
  implicit none

  private
  public :: domain, time, element, match_result, matching, circle, ellipse, &
       & INVLT, solution, particle, mathieu, explain_type, print_match_result

  type :: mathieu
     ! things required to compute mathieu functions
     integer :: M,  buffer
     complex(DP) :: q
     complex(DP), allocatable :: mcn(:), A(:,:,:), B(:,:,:)  ! 4*M; 1:M, 0:M-1, 0:1
  end type mathieu

  type :: domain
     ! number of each type of element
     ! 1=circles (wells as special case), 2=ellipses (lines as special case)
     integer, dimension(2) :: num

     ! index of parent element
     integer, allocatable :: InclUp(:)

     ! matrix indicating if an element is inside or in the background of
     ! a current element
     logical, allocatable :: InclIn(:,:), InclBg(:,:)
  end type domain

  type :: explain_type
    ! explanation text (only used in ltaem_io.f90)
    character(kind=ASCII, len=55), dimension(10) :: time = [&
         & '[step on @ tpar(1)]                                    ',&
         & '[finite width pulse; tpar(1:2) = on/off times]         ',&
         & '[instantaneous pulse @ tpar(1)]                        ',&
         & '[stairs; stair width=tpar(1), off @ tpar(2)]           ',&
         & '[+only square wave; tpar(1)=period/2, off @ tpar(2)]   ',&
         & '[cos(tpar(1)*t); on @ tpar(2)]                         ',&
         & '[+only triangular wave; tpar(1)=period/4, on @ tpar(2)]',&
         & '[+/- square wave; tpar(1)=period/4, on @ tpar(2)]      ',&
         & '[piecewise-constant rate]                              ',&
         & '[piecewise-linear rate]                                ']

    character(kind=ASCII, len=39), dimension(0:3) :: leakFlag = &
         & ['(no leakage)                           ',&
         &  '(no-drawdown condition beyond aquitard)',&
         &  '(no-flow condition beyond aquitard)    ',&
         &  '(infinitely thick aquitard)            ']

    character(kind=ASCII, len=24), dimension(-1:2) :: ibnd = [&
         & '(specified total head)  ',&
         & '(matching)              ',&
         & '(specified total flux)  ',&
         & '(specified element flux)']

    character(kind=ASCII, len=34), dimension(7) :: output = [&
         & '(gnuplot contour map)             ',&
         & '(matlab contour map)              ',&
         & '(gnuplot time series w/ velocity) ',&
         & '(inverse time series w/o velocity)',&
         & '(gnuplot time series w/o velocity)',&
         & '(gnuplot pathlines)               ',&
         & '(gnuplot streakline)              ']

    character(kind=ASCII, len=39), dimension(4) :: particle = [&
         & 'Runge-Kutta-Merson (4th order adaptive)',&
         & 'Runge-Kutta  (4th order)               ',&
         & 'Analytical   (root-finding)            ',&
         & 'Forward Euler  (1st order)             ']
  end type explain_type

  type :: time
     ! all elements inherit the time behavior from this type

     ! time behavior / parameters
     ! 1 = step on,              tpar(1)  = on time
     ! 2 = finite pulse,         tpar(1:2) = on/off time
     ! 3 = instan. pulse         tpar(1) = pulse time
     ! 4 = stairs,               tpar(1) = time step (increasing Q by integer multiples
     !                                     @ integer multiples tpar(1)); tpar(2) =off time.
     ! 5 = + only square wave,   tpar(1) = 1/2 period of wave; tpar(2) = start time
     ! 6 = cosine(tpar(1)*t),    tpar(1) = frequency multiplier; tpar(2) = start time
     ! 7 = + only tri wave,      tpar(1) = 1/4 period of wave; tpar(2) = start time
     ! 8 = +/- square wave,      tpar(1) = 1/2 period of wave; tpar(2) = start time
     ! -100<n<0 = arbitrary piecewise constant rate, comprised of n steps from tpar(1) to tfinal
     !                tpar(1:n) = starting times of each step
     !                tpar(n+1) = final time of last step
     !                tpar(n+2:2*n+1) = strength at each of n steps
     ! (is multiplied by constant strength too -- you probably want to set that to unity)
     ! n<-100 = arbitrary piecewise linear rate, comprised of n steps from tpar(1) to tfinal
     !                tpar(1:n) = starting times of each linear segment
     !                tpar(n+1) = final time of last linear segment
     !                tpar(n+2:2*n+1) = strength at each of n linear segment
     ! (is multiplied by constant strength too -- you probably want to set that to unity)

     ! type of time behavior for AREA/Boundary Head/Flux (see above)
     integer :: AreaTime, BdryTime

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), allocatable :: ATPar(:), BTPar(:)
  end type time

  type, extends(time) :: element

    logical :: debug
    logical :: wave   ! solve wave problem (2nd time deriv)?

     ! global id for the current element
     integer :: id

     ! porosity, constant area source term
     ! main aquifer hydraulic conductivity and Ss for element
     real(DP) :: por, k, Ss, b, alpha, T

     ! leaky-related (adjoining aquitard/aquifer parameters)
     ! 0= no leakage
     ! 1= case I, no-drawdown condition at top of aquitard
     ! 2= case II, no-flow condition at top of aquitard
     ! 3= aquitard thickness -> infinity (no bc)
     integer :: leakFlag
     real(DP) :: aquitardK, aquitardSs, aquitardb

     ! unconfined-related (flag, specific yield, and vertical K)
     logical ::  unconfFlag
     real(DP) :: Sy, Kz

     ! dual-porosity related (flag, matrix storativity, matrix/fracture exchange param)
     logical :: dualPFlag
     integer :: multiPDiff ! 1=slab, 2=cylinder, 3=spherical
     real(DP) :: matrixSs, lambda, kappa
     integer :: NDiffterms

     ! specified value across area of element (including background)
     real(DP) :: areaQ

     ! whether to calculate solution (Helmholtz eqn) inside element
     ! and whether to compute storage (using mass conservation ODE) inside element
     ! StorIn is only checked if CalcIn is false for an element.
     logical :: CalcIn, StorIn

     ! pointer to the parent element
     type(element), pointer :: parent => null()

     ! structure containing matrices of mathieu function parameters
     type(mathieu), allocatable :: mat(:)
     integer :: ms ! not used in circle

  end type element

  type :: geom
     ! number of matching points along other elements
     real(DP), allocatable :: Rgm(:), Pgm(:), metric(:)
  end type geom

  type :: match_result
     ! structure for storing intermediate results
     complex(DP), allocatable :: LHS(:,:), RHS(:)
  end type match_result

  type, extends(element) :: matching
     ! number of FS terms, number of matching points on circles/ellipses
     ! for lines/wells can be one (e.g., borehole storage) or zero (known Q)
     integer :: n, m

     ! type of element: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
     !                  -2=specified head ELEMENT, +2=specified flux ELEMENT
     ! -2 doesn't really make sense from a physical perspective : not implemented
     integer :: ibnd

     ! whether inclusion is a matching(T) or specified(F) inclusion
     logical :: match

     ! specified value along boundary area of element
     real(DP) :: bdryQ

     ! dimensionless skin at boundary of element
     real(DP) :: dskin

     ! location of center of element
     real(DP) :: x, y
     complex(DP) :: z

     ! radius of element (eta for ellipses)
     ! semi-focal distance (zero for circles)
     ! angle of rotation (zero for circles)
     real(DP) :: r, f, theta

     ! vector of matching location angles
     ! theta for circles, psi for ellipses
     real(DP), allocatable :: Pcm(:)

     ! computed geometry from this element to all others
     complex(DP), allocatable :: Zom(:)
     type(geom), allocatable :: G(:) ! number of elements

     ! coefficients determined through matching
     complex(DP), allocatable :: coeff(:,:)
  end type matching

  type, extends(matching) :: circle
     ! Circular Inclusion related parameters
     ! well is special case of circle
  end type circle

  type, extends(matching) :: ellipse
     ! Elliptical Inclusion related parameters
     ! line is special case of ellipse
  end type ellipse

  type :: INVLT
     ! Inverse Laplace Transform parameters

     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha, tol

     ! number of Fourier series terms
     integer :: M
  end type INVLT

  ! things relating to the numerical solution, independent from flow elements
  type, extends(INVLT) :: solution

     logical :: debug  ! debugging output?

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle
     integer :: nPart

     ! compute flowrates across each element for contour/time series?
     logical :: Qcalc

     ! dump matching results to file for restart?
     logical :: skipDump

     ! total number of laplace parameters
     integer :: totalnP

     ! number of particle timesteps to skip when plotting streaklines
     ! (only one value for all particles)
     integer :: streakSkip

     ! calculate contours (thru space) vs. calculate time series
     logical :: contour
     ! calculate values at few locations through time
     logical :: timeseries
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc

     ! input/output filenames
     character(kind=ASCII, len=lenFN) :: outfName, infName
     character(kind=ASCII, len=lenFN) :: elemHfName, geomfName
     character(kind=ASCII, len=lenFN) :: echofName, qfName
     character(kind=ASCII, len=13) :: coefffName = 'dump-vars.out'

     ! output index
     ! ------------------- <10 = contour map output --------------------
     !  1= Gnuplot map (x,y,z triplets; times separated by blank lines);
     !  2= Matlab map (matrix output separate files);
     ! --------------->=10 <20 = time series output --------------------
     ! 10= Gnuplot time series with velocity (column of times; locs sep. by blank lines);
     ! 11= Gnuplot time serires no velocity (same as 10 no vel);
     ! 12= Matlab for SCEM-UA inverse (column of times, locs sep. by blank lines);
     ! --------------->=20 = particle track output ---------------------
     ! 20= pathline Gnuplot (column of times, particles separated by blank lines);
     ! 21= streakline Gnuplot (each block a requested time, each row a particle);

     integer, dimension(21) :: OEMap = [1,2,0,0,0,0,0,0,0,&
                                      & 3,4,5,0,0,0,0,0,0,0,&
                                      & 6,7]

     ! aquitardLeak and unconfined
     integer :: output, aquitardLeak, unconfined

     ! x-spacing vector, y-spacing vector, time vector
     integer :: nx, ny, nt
     real(DP) :: xshift, yshift
     real(DP), allocatable :: x(:), y(:), t(:)
     character(kind=ASCII, len=32), allocatable :: obsname(:)

     ! compute log-derivative of solution?
     logical :: deriv

     ! containers for time-domain final results (x,y,t,[i:j])
     real(DP),    allocatable :: h(:,:,:), v(:,:,:,:), dh(:,:,:), Q(:,:), dQ(:,:)

  end type solution

  ! particle related parameters (one for each particle)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type :: particle
     logical :: debug

     ! is particle tracked backwards or forwards?
     logical :: forward

     ! starting x&y location, inital & final times
     real(DP) :: x, y, ti, tf

     ! which integration scheme to use (for each particle)
     ! 1 = Runge-Kutta-Merson (4th order adaptive)
     ! 2 = Runge-Kutta        (4th order)
     ! 3 = Analytical        (root-finding method)
     ! 4 = Fwd Euler          (1st order)
     integer :: int, id

     ! error tolerance, minimum stepsize, and max step length for rkm
     real(DP) :: tol, mindt, maxL

     ! step size for other integration schemes (initial stepsize for rkm)
     real(DP) :: dt

     ! particle starts inside a constant head or constant flux inclusion?
     logical :: InclIn

     ! results from particle tracking (for each particle)
     ! first dimension is long enough to hold all needed times
     ! second dimension is always 5 (t,x,y,velx,vely)
     real(DP), allocatable :: r(:,:)

     ! number of time steps actuall used (some algorithms are adaptive)
     integer :: numt
  end type particle

contains
  subroutine print_match_result(r)
    use constants, only : ASCII
    type(match_result), intent(in) :: r
    integer :: row,col,i,j
    character(kind=ASCII, len=40), dimension(2) :: fmt

    write(*,'(A)') '** matching result **'
    row = size(r%RHS,dim=1)
    col = size(r%LHS,dim=2)
    if (row > 0 .and. col > 0) then
       write(*,'(2(A,2(1X,I0)))') 'Shape: r%LHS:',&
            & shape(r%LHS),' r%RHS:',shape(r%RHS)
       fmt(1) = '(A,I3,   (A,ES10.2E3,A,ES10.2E3,A))     '
       write(fmt(1)(7:9),'(I3.3)') col
       fmt(2) = '(7X,   (I12,12X))                       '
       write(fmt(2)(5:7),'(I3.3)') col
       write(*,fmt(2)) [(i,i=1,col)] ! headers
       do i = 1,row
          write(*,fmt(1)) 'LHS:',i,('(',real(r%LHS(i,j)),',',&
               & aimag(r%LHS(i,j)),') ',j=1,col)
       end do
       do i = 1,row
          write(*,'(A,I3,2(A,ES10.2E3),A)') 'RHS:',i,'(',&
               & real(r%RHS(i)),',',aimag(r%RHS(i)),')'
       end do

    elseif(row > 0) then
       do i = 1,row
          write(*,'(A,I3,2(A,ES10.2E3),A)') 'RHS:',i,'(',&
               & real(r%RHS(i)),',',aimag(r%RHS(i)),')'
       end do
    else
       write(*,*) '* nothing to print * row:',row,'col:',col
    end if
  end subroutine print_match_result
end module type_definitions
