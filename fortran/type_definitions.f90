module type_definitions
  use constants, only : DP, lenFN
  use mathieu_functions, only : mathieu
  implicit none

  private 
  public :: domain, time, element, match_result, matching, circle, ellipse, &
       & INVLT, solution, particle, print_match_result
  
  type :: domain
     ! number of each type of element
     ! 1=circles (wells as special case), 2=ellipses (lines as special case)
     integer, dimension(2) :: num = [-999,-999]

     ! index of parent element
     integer, allocatable :: InclUp(:)

     ! matrix indicating if an element is inside or in the background of
     ! a current element
     logical, allocatable :: InclIn(:,:), InclBg(:,:)
  end type domain
  
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
     ! n<0 = arbitrary piecewise constant rate, comprised of n steps from tpar(1) to tfinal
     !                tpar(1:n) = starting times of each step
     !                tpar(n+1) = final time of last step
     !                tpar(n+2:2*n+1) = strength at each of n steps 
     ! (is multiplied by constant strength too -- you probably want to set that to unity)

     ! type of time behavior for AREA/Boundary Head/Flux (see above)
     integer :: AreaTime = -999, BdryTime = -999

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), allocatable :: AtPar(:), BtPar(:)
  end type time
 
  type, extends(time) :: element

     ! global id for the current element
     integer :: id = -999

     ! porosity, constant area source term
     ! main aquifer hydraulic conductivity and Ss for element
     real(DP) :: por = -999., k = -999., Ss = -999., b = -999., alpha = -999., T = -999.

     ! leaky-related (adjoining aquitard/aquifer parameters)
     ! 0= no leakage
     ! 1= case I, no-drawdown condition at top of aquifer
     ! 2= case II, no-flow condition at top of aquifer
     ! 3= aquitard thickness -> infinity (no bc)
     integer :: leakFlag  = -1
     real(DP) :: aquitardK = -999., aquitardSs = -999., aquitardb = -999.

     ! unconfined-related (flag, specific yield, and vertical K)
     logical ::  unconfinedFlag = .false.
     real(DP) :: Sy = -999., Kz = -999.

     ! specified value across area of element (including background)
     real(DP) :: areaQ = -999.   

     ! whether to calculate solution (Helmholtz eqn) inside element
     ! and whether to compute storage (using mass conservation ODE) inside element
     ! StorIn is only checked if CalcIn is false for an element.
     logical :: CalcIn = .false., StorIn = .false.

     ! the parent element
     type(element), pointer :: parent => null()

     ! structure containing matrices of mathieu function parameters
     type(mathieu), allocatable :: mat(:)
     integer :: ms = -999 ! not used in circle

  end type element
    
  type :: geom
     ! number of matching points along other elements
     complex(DP), allocatable :: Zgm(:)
     real(DP), allocatable :: Rgm(:), Pgm(:), metric(:)
  end type geom
  
  type :: match_result
     ! structure for storing intermediate results
     complex(DP), allocatable :: LHS(:,:), RHS(:)
  end type match_result

  type, extends(element) :: matching
     ! number of FS terms, number of matching points on circles/ellipses
     ! for lines/wells can be one (e.g., borehole storage) or zero (known Q)
     integer :: n = -999, m = -999
  
     ! type of element: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
     !                  -2=specified head ELEMENT, +2=specified flux ELEMENT
     ! -2 doesn't really make sense from a physical perspective : not implemented
     integer :: ibnd = -999
     
     ! whether inclusion is a matching(T) or specified(F) inclusion
     logical :: match = .false.

     ! specified value along boundary area of element
     real(DP) :: bdryQ = -999.

     ! dimensionless skin at boundary of element 
     real(DP) :: dskin = -999.

     ! location of center of element
     real(DP) :: x = -999., y =-999.

     ! "radius" of element (eta for ellipses)
     ! semi-focal distance (zero for circles)
     ! angle of rotation (zero for circles)
     real(DP) :: r = -999., f = -999., theta = -999.

     ! vector of matching location angles
     ! theta for circles, psi for ellipses
     real(DP), allocatable :: Pcm(:)

     ! computed geometry from this element to all others
     complex(DP), allocatable :: Zcm(:), Zom(:)
     type(geom), allocatable :: G(:) ! number of elements

     ! coefficients determined through matching
     complex(DP), allocatable :: coeff(:,:)
  end type matching

  type, extends(matching) :: circle
     ! Circular Inclusion related parameters
     ! well is special case of circle

     ! no special circle-only parameters
  end type circle

  type, extends(matching) :: ellipse
     ! Elliptical Inclusion related parameters
     ! line is special case of ellipse

  end type ellipse

  type :: INVLT
     ! Inverse Laplace Transform parameters

     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha = -999., tol = -999.

     ! number of Fourier series terms
     integer :: M = -999
  end type INVLT

  ! things relating to the numerical solution, independent from flow elements
  type, extends(INVLT) :: solution

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle = .false.
     integer :: nPart = -999

     ! total number of laplace parameters
     integer :: totalnP = -999

     ! number of particle timesteps to skip when plotting streaklines
     ! (only one value for all particles)
     integer :: streakSkip = -999

     ! calculate contours (thru space) vs. calculate hydrographs (thru time)
     logical :: contour = .false.
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc = .false.
     
     ! input/output filenames
     character(lenFN) :: outfname='unset', infname='unset', &
          & coefffname='unset', elemHfName='unset', geomFname='unset'
     
     ! output index (1= Gnuplot map (x,y,z triplets; times separated by blank lines);
     !               2= Matlab map (matrix output separate files);
     !               3= Gnuplot hydrograph with velocity (column of times; locs sep. by blank lines);
     !               4= pathline Gnuplot (column of times, particles separated by blank lines);
     !               5= streakline Gnuplot (each block a requested time, each row a particle);
     !               10= Matlab for SCEM-UA inverse (column of times, locs sep. by blank lines);
     !               11= Gnuplot hydrograph no velocity (same as 3 no vel);)
     ! aquitardLeak and unconfined
     integer :: output = -999, aquitardLeak = -999, unconfined= -999
     
     ! x-spacing vector, y-spacing vector, time vector
     integer :: nx = -999, ny = -999, nt = -999
     real(DP), allocatable :: x(:), y(:), t(:)

     ! containers for time-domain final results (x,y,t,[i:j])
     real(DP),    allocatable :: h(:,:,:), v(:,:,:,:)

  end type solution

  ! particle related parameters (one for each particle)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type :: particle
  
     ! starting x&y location, inital & final times
     real(DP) :: x = -999., y = -999., ti = -999., tf = -999.
  
     ! which integration scheme to use (for each particle??)
     ! 1 = Runge-Kutta-Merson (4th order adaptive)
     ! 2 = Runge-Kutta        (4th order)
     ! 3 = Richardson Extrapolation (2nd order)
     ! 4 = Fwd Euler          (1st order)
     integer :: int = -999
     
     ! error tolerance and minimum stepsize for rkm
     real(DP) :: tol = -999., min = -999., maxStep = -999.
     
     ! step size for other integration schemes (initial stepsize for rkm)
     real(DP) :: dt = -999.
     
     ! particle starts inside a constant head or constant flux inclusion?
     logical :: InclIn = .false.
     
     ! results from particle tracking (for each particle)
     ! first dimension is always 5 (t,x,y,velx,vely)
     ! second dimension is long enough to hold all needed times
     real(DP), allocatable :: result(:,:)

     ! number of time steps actuall used (some algorithms are adaptive)
     integer :: numt = -999
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
       do i=1,row
          write(*,fmt(1)) 'LHS:',i,('(',real(r%LHS(i,j)),',',&
               & aimag(r%LHS(i,j)),') ',j=1,col)
       end do
       do i=1,row
          write(*,'(A,I3,2(A,ES10.2E3),A)') 'RHS:',i,'(',&
               & real(r%RHS(i)),',',aimag(r%RHS(i)),')'
       end do

    elseif(row > 0) then
       do i=1,row
          write(*,'(A,I3,2(A,ES10.2E3),A)') 'RHS:',i,'(',&
               & real(r%RHS(i)),',',aimag(r%RHS(i)),')'
       end do
    else
       write(*,*) '* nothing to print * row:',row,'col:',col
    end if

  end subroutine print_match_result
  
end module type_definitions
