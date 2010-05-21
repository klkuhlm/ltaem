module element_specs
  use constants, only : DP
  implicit none
  
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
     ! all element inherit this time behavior

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

     ! type of time behavior for AREA/Boundary Head/Flux 1=step, 2=pulse, 3=stair, ...
     integer :: AreaTime, BdryTime

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), allocatable :: AtPar(:), BtPar(:)
  end type time
 
  type, public, extends(time) :: element

     ! global id for the current element
     integer :: id

     ! porosity, constant area source term
     ! main aquifer hydraulic conductivity and Ss for element
     real(DP) :: por, area, k, Ss, b, alpha

     ! leaky-related (adjoining aquitard/aquifer parameters)
     integer :: leakFlag
     real(DP) :: aquitardK, aquitardSs, aquitardb

     ! unconfined-related (flag, specific yield, and vertical K)
     integer ::  unconfinedFlag
     real(DP) :: Sy, Kz

     ! whether to calculate solution (Helmholtz eqn) inside element
     ! and whether to compute storage (using mass conservation ODE) inside element
     ! StorIn is only checked if CalcIn is false for an element.
     logical :: CalcIn, StorIn

     ! the parent element
     type(element), pointer :: parent => null()

  end type element
    
  type, public, extends(element) :: matching
     ! number of FS terms, number of matching points on circles/ellipses
     ! for lines/wells can be one (e.g., borehole storage) or zero (known Q)
     integer :: n, m
  
     ! type of element: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
     !                  -2=specified head ELEMENT, +2=specified flux ELEMENT
     integer :: ibnd
     
     ! whether inclusion is a matching(T) or specified(F) inclusion
     logical :: match

     ! specified value on bdry of inclusion
     real(DP) :: spec    

     ! dimensionless skin at boundary of element 
     real(DP) :: dskin

     ! location of center of element
     real(DP) :: x, y 

  end type matching

  type, extends(matching) :: circle
     ! Circular Inclusion related parameters
     ! well is special case of circle

     ! inclusion radius
     real(DP) :: r 

     ! vector of matching locations 
     real(DP), allocatable :: Pcm(:)

  end type circle

  type, public, extends(matching) :: ellipse
     ! Elliptical Inclusion related parameters
     ! line is special case of ellipse

     ! size of MF infinite matrix
     integer, save :: ms

     ! inclusion elliptical 'radius', semi-focal dist, 
     ! angle major axis makes with Cartesian x axis
     real(DP) :: eta, f, theta

     ! vector of matching locations 
     real(DP), allocatable :: Pcm(:)

  end type ellipse

  type, public :: INVLT
     ! INVerse Laplace Transform parameters

     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha, tol
     integer :: m
  end type INVLT

  ! things relating to the numerical solution, independent from flow elements
  type, public, extends(laplace) :: solution

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle
     integer :: nPart

     ! calculate contours (thru space) vs. calculate hydrographs (thru time)
     logical :: contour
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc
     
     ! input/output filenames
     character(128) :: outfname, infname, coefffname, elemHfName
     
     ! output index (1= gnuplot; 2= matlab)
     integer :: output, aquitardLeak, unconfined
     
     ! x-spacing vector, y-spacing vector, time vector
     integer :: nx, ny, nt
     real(DP) allocatable :: x(:), y(:), t(:)

     ! containers for time-domain final results
     real(DP),    allocatable :: h(:,:,:), vx(:,:,:), vy(:,:,:)

     ! container for Laplace-space intermediate results
     complex(DP), allocatable :: hp(:), vxp(:), vyp(:)
     complex(DP), allocatable :: coeff(:,:,:,:) 
     complex(DP), allocatable :: Gm(:,:,:)
  end type solution

  ! PARticle related parameters (one for each particle)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, public :: particle
  
     ! number of particles
     integer :: streakSkip

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
     real(DP), allocatable :: result(:,:)

  end type particle

end module element_specs
