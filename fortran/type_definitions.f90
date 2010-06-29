module type_definitions
  use constants, only : DP, lenFN
  use mathieu_functions, only : mathieu
  implicit none

  public
  
  type, public :: domain
     ! number of each type of element
     ! 1=circles (wells as special case), 2=ellipses (lines as special case)
     integer, dimension(2) :: num

     ! index of parent element
     integer, allocatable :: InclUp(:)

     ! matrix indicating if an element is inside or in the background of
     ! a current element
     logical, allocatable :: InclIn(:,:), InclBg(:,:)
  end type domain
  
  type, public :: time
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
     integer :: AreaTime, BdryTime

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), allocatable :: AtPar(:), BtPar(:)
  end type time
 
  type, extends(time) :: element

     ! global id for the current element
     integer :: id

     ! porosity, constant area source term
     ! main aquifer hydraulic conductivity and Ss for element
     real(DP) :: por, k, Ss, b, alpha, T

     ! leaky-related (adjoining aquitard/aquifer parameters)
     ! 0= no leakage
     ! 1= case I, no-drawdown condition at top of aquifer
     ! 2= case II, no-flow condition at top of aquifer
     ! 3= aquitard thickness -> infinity (no bc)
     integer :: leakFlag  
     real(DP) :: aquitardK, aquitardSs, aquitardb

     ! unconfined-related (flag, specific yield, and vertical K)
     logical ::  unconfinedFlag
     real(DP) :: Sy, Kz

     ! specified value across area of element (including background)
     real(DP) :: areaQ    

     ! whether to calculate solution (Helmholtz eqn) inside element
     ! and whether to compute storage (using mass conservation ODE) inside element
     ! StorIn is only checked if CalcIn is false for an element.
     logical :: CalcIn, StorIn

     ! the parent element
     type(element), pointer :: parent => null()

     ! structure containing matrices of mathieu function parameters
     type(mathieu), allocatable :: mat(:)
     integer :: ms
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

     ! "radius" of element (eta for ellipses)
     ! semi-focal distance (zero for circles)
     ! angle of rotation (zero for circles)
     real(DP) :: r, f, theta

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
     real(DP) :: alpha, tol

     ! number of Fourier series terms
     integer :: M
  end type INVLT

  ! things relating to the numerical solution, independent from flow elements
  type, extends(INVLT) :: solution

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle
     integer :: nPart

     integer :: totalnP ! total number of laplace parameters

     ! number of particle timesteps to skip when plotting streaklines
     ! (only one value for all particles)
     integer :: streakSkip

     ! calculate contours (thru space) vs. calculate hydrographs (thru time)
     logical :: contour
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc
     
     ! input/output filenames
     character(lenFN) :: outfname, infname, coefffname, elemHfName, geomFname
     
     ! output index (1= Gnuplot map (x,y,z triplets; times separated by blank lines);
     !               2= Matlab map (matrix output separate files);
     !               3= Gnuplot hydrograph with velocity (column of times; locs sep. by blank lines);
     !               4= pathline Gnuplot (column of times, particles separated by blank lines);
     !               5= streakline Gnuplot (each block a requested time, each row a particle);
     !               10= Matlab for SCEM-UA inverse (column of times, locs sep. by blank lines);
     !               11= Gnuplot hydrograph no velocity (same as 3 no vel);)
     ! aquitardLeak and unconfined
     integer :: output, aquitardLeak, unconfined
     
     ! x-spacing vector, y-spacing vector, time vector
     integer :: nx, ny, nt
     real(DP), allocatable :: x(:), y(:), t(:)

     ! containers for time-domain final results (x,y,t,[i:j])
     real(DP),    allocatable :: h(:,:,:), v(:,:,:,:)

     ! container for Laplace-space intermediate results
     complex(DP), allocatable :: hp(:), vp(:,:)
  end type solution

  ! particle related parameters (one for each particle)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type :: particle
  
     ! starting x&y location, inital & final times
     real(DP) :: x, y, ti, tf
  
     ! which integration scheme to use (for each particle??)
     ! 1 = Runge-Kutta-Merson (4th order adaptive)
     ! 2 = Runge-Kutta        (4th order)
     ! 3 = Richardson Extrapolation (2nd order)
     ! 4 = Fwd Euler          (1st order)
     integer :: int
     
     ! error tolerance and minimum stepsize for rkm
     real(DP) :: tol, min, maxStep
     
     ! step size for other integration schemes (initial stepsize for rkm)
     real(DP) :: dt
     
     ! particle starts inside a constant head or constant flux inclusion?
     logical :: InclIn
     
     ! results from particle tracking (for each particle)
     ! first dimension is always 5 (t,x,y,velx,vely)
     ! second dimension is long enough to hold all needed times
     real(DP), allocatable :: result(:,:)

     ! number of time steps actuall used (some algorithms are adaptive)
     integer :: numt
  end type particle

contains
  subroutine print_match_result(r)
    type(match_result), intent(in) :: r
    integer :: row,col,i,j
    character(40),dimension(2) :: fmt

    write(*,'(A,/,2(A,1L))') '** matching result **','Allocated?  r%LHS:',allocated(r%LHS),' r%RHS:',allocated(r%RHS)
    if (allocated(r%LHS) .and. allocated(r%RHS)) then
       row = size(r%LHS,dim=1)
       col = size(r%LHS,dim=2)
       write(*,'(2(A,2(1X,I0)))') 'Shape: r%LHS:',shape(r%LHS),' r%RHS:',shape(r%RHS)
       fmt(1) = '(A,I3,   (A,ES13.5E3,A,ES13.5E3,A))     '
       write(fmt(1)(7:9),'(I3.3)') col
       fmt(2) = '(7X,   (I15,15X))                       '
       write(fmt(2)(5:7),'(I3.3)') col
       write(*,fmt(2)) [(i,i=1,col)] ! headers
       do i=1,row
          write(*,fmt(1)) 'LHS:',i,('(',real(r%LHS(i,j)),',',aimag(r%LHS(i,j)),') ',j=1,col)
       end do
       do i=1,row
          write(*,'(A,I3,(A,ES13.5E3,A,ES13.5E3,A))') 'RHS:',i,'(',real(r%RHS(i)),',',aimag(r%RHS(i)),')'
       end do
    end if
    
  end subroutine print_match_result
  
end module type_definitions
