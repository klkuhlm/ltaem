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

module type_definitions
  use constants, only : DP, lenFN
  implicit none

  private
  public :: domain, time, element, match_result, matching, circle, ellipse, &
       & INVLT, solution, particle, mathieu, explain_type, print_match_result

  type :: mathieu
     ! things required to compute mathieu functions
     integer :: M = -9999,  buffer = -9999
     complex(DP) :: q = (-9.9D+20,-9.9D+20)
     complex(DP), allocatable :: mcn(:), A(:,:,:), B(:,:,:)  ! 4*M; 1:M, 0:M-1, 0:1
  end type mathieu

  type :: domain
     ! number of each type of element
     ! 1=circles (wells as special case), 2=ellipses (lines as special case)
     integer, dimension(2) :: num = [-9999,-9999]

     ! index of parent element
     integer, allocatable :: InclUp(:)

     ! matrix indicating if an element is inside or in the background of
     ! a current element
     logical, allocatable :: InclIn(:,:), InclBg(:,:)
  end type domain

  type :: explain_type
    character(55), dimension(10) :: time = [&
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
    
    character(39), dimension(0:3) :: leakFlag = &
         & ['(no leakage)                           ',&
         &  '(no-drawdown condition beyond aquitard)',&
         &  '(no-flow condition beyond aquitard)    ',&
         &  '(infinitely thick aquitard)            ']
    
    character(24), dimension(-1:2) :: ibnd = [&
         & '(specified total head)  ',&
         & '(matching)              ',&
         & '(specified total flux)  ',&
         & '(specified element flux)']
    
    character(34), dimension(7) :: output = [&
         & '(gnuplot contour map)             ',&
         & '(matlab contour map)              ',&
         & '(gnuplot time series w/ velocity) ',&
         & '(inverse time series w/o velocity)',&
         & '(gnuplot time series w/o velocity)',&
         & '(gnuplot pathlines)               ',&
         & '(gnuplot streakline)              ']
    
    character(39), dimension(4) :: particle = [&
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
     integer :: AreaTime = -9999, BdryTime = -9999

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), allocatable :: ATPar(:), BTPar(:)
  end type time

  type, extends(time) :: element

     ! global id for the current element
     integer :: id = -9999

     ! porosity, constant area source term
     ! main aquifer hydraulic conductivity and Ss for element
     real(DP) :: por = -9.9D+20, k = -9.9D+20, Ss = -9.9D+20, b = -9.9D+20, alpha = -9.9D+20, T = -9.9D+20

     ! leaky-related (adjoining aquitard/aquifer parameters)
     ! 0= no leakage
     ! 1= case I, no-drawdown condition at top of aquitard
     ! 2= case II, no-flow condition at top of aquitard
     ! 3= aquitard thickness -> infinity (no bc)
     integer :: leakFlag  = -9999
     real(DP) :: aquitardK = -9.9D+20, aquitardSs = -9.9D+20, aquitardb = -9.9D+20

     ! unconfined-related (flag, specific yield, and vertical K)
     logical ::  unconfinedFlag = .false.
     real(DP) :: Sy = -9.9D+20, Kz = -9.9D+20

     ! dual-porosity related (flag, matrix storativity, matrix/fracture exchange param)
     logical :: dualPorosityFlag = .false.
     integer :: multiporosityDiffusion = 0 ! 1=slab, 2=cylinder, 3=spherical
     real(DP) :: matrixSs = -9.9D+20, lambda = -9.9D+20, kappa = -9.9D+20
     integer :: NDiffterms = -9999

     ! specified value across area of element (including background)
     real(DP) :: areaQ = -9.9D+20

     ! whether to calculate solution (Helmholtz eqn) inside element
     ! and whether to compute storage (using mass conservation ODE) inside element
     ! StorIn is only checked if CalcIn is false for an element.
     logical :: CalcIn = .false., StorIn = .false.

     ! pointer to the parent element
     type(element), pointer :: parent => null()

     ! structure containing matrices of mathieu function parameters
     type(mathieu), allocatable :: mat(:)
     integer :: ms = -9999 ! not used in circle

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
     integer :: n = -9999, m = -9999

     ! type of element: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
     !                  -2=specified head ELEMENT, +2=specified flux ELEMENT
     ! -2 doesn't really make sense from a physical perspective : not implemented
     integer :: ibnd = -9999

     ! whether inclusion is a matching(T) or specified(F) inclusion
     logical :: match = .false.

     ! specified value along boundary area of element
     real(DP) :: bdryQ = -9.9D+20

     ! dimensionless skin at boundary of element
     real(DP) :: dskin = -9.9D+20

     ! location of center of element
     real(DP) :: x = -9.9D+20, y =-9.9D+20
     complex(DP) :: z = (-9.9D+20,-9.9D+20)

     ! "radius" of element (eta for ellipses)
     ! semi-focal distance (zero for circles)
     ! angle of rotation (zero for circles)
     real(DP) :: r = -9.9D+20, f = -9.9D+20, theta = -9.9D+20

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
     real(DP) :: alpha = -9.9D+20, tol = -9.9D+20

     ! number of Fourier series terms
     integer :: M = -9999
  end type INVLT
  
  ! things relating to the numerical solution, independent from flow elements
  type, extends(INVLT) :: solution

     logical :: debug = .false.  ! debugging output?

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle = .false.
     integer :: nPart = -9999

     ! compute flowrates across each element for contour/time series?
     logical :: Qcalc = .false.

     ! dump matching results to file for restart?
     logical :: skipDump = .false.

     ! total number of laplace parameters
     integer :: totalnP = -9999

     ! number of particle timesteps to skip when plotting streaklines
     ! (only one value for all particles)
     integer :: streakSkip = -9999

     ! calculate contours (thru space) vs. calculate time series
     logical :: contour = .false.
     ! calculate values at few locations through time
     logical :: timeseries = .false.
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc = .false.

     ! input/output filenames
     character(lenFN) :: outfName='ERROR-unset', infName='ERROR-unset'
     character(lenFN) :: elemHfName='ERROR-unset', geomfName='ERROR-unset'
     character(lenFN) :: echofName='ERROR-unset', qfName='ERROR-unset'
     character(13) :: coefffName = 'dump-vars.out'

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
     integer :: output = -9999, aquitardLeak = -9999, unconfined= -9999

     ! x-spacing vector, y-spacing vector, time vector
     integer :: nx = -9999, ny = -9999, nt = -9999
     real(DP) :: xshift = -9.9D+20, yshift = -9.9D+20
     real(DP), allocatable :: x(:), y(:), t(:)
     character(32), allocatable :: obsname(:)

     ! compute log-derivative of solution?
     logical :: deriv = .false.

     ! containers for time-domain final results (x,y,t,[i:j])
     real(DP),    allocatable :: h(:,:,:), v(:,:,:,:), dh(:,:,:), Q(:,:), dQ(:,:)

  end type solution
  
  ! particle related parameters (one for each particle)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type :: particle

     ! is particle tracked backwards or forwards?
     logical :: forward = .true.

     ! starting x&y location, inital & final times
     real(DP) :: x = -9.9D+20, y = -9.9D+20, ti = -9.9D+20, tf = -9.9D+20

     ! which integration scheme to use (for each particle)
     ! 1 = Runge-Kutta-Merson (4th order adaptive)
     ! 2 = Runge-Kutta        (4th order)
     ! 3 = Analytical        (root-finding method)
     ! 4 = Fwd Euler          (1st order)
     integer :: int = -9999, id = -9999
     
     ! error tolerance, minimum stepsize, and max step length for rkm
     real(DP) :: tol = -9.9D+20, mindt = -9.9D+20, maxL = -9.9D+20

     ! step size for other integration schemes (initial stepsize for rkm)
     real(DP) :: dt = -9.9D+20

     ! particle starts inside a constant head or constant flux inclusion?
     logical :: InclIn = .false.

     ! results from particle tracking (for each particle)
     ! first dimension is long enough to hold all needed times
     ! second dimension is always 5 (t,x,y,velx,vely)
     real(DP), allocatable :: r(:,:)

     ! number of time steps actuall used (some algorithms are adaptive)
     integer :: numt = -9999
  end type particle

contains
  subroutine print_match_result(r)
    type(match_result), intent(in) :: r
    integer :: row,col,i,j
    character(40), dimension(2) :: fmt

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

