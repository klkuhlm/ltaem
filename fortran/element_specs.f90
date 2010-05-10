!##################################################
module element_specs
  use constants, only : DP
  implicit none
  
  type, private :: time
     ! time behavior / parameters (WLtpar())
     ! 1=step on,              WLtpar(1)  =on time
     ! 2=finite pulse,         WLtpar(1:2)=on/off time
     ! 3=instan. pulse,  -  no parameters
     ! 4=stairs,               WLtpar(1:2)=time step, off time
     ! 5=1/2 square wave,      WLtpar(1)  =1/2 period of wave
     ! 6=sine,                 WLtpar(1)=1/amplitude & wave number
     ! 7=1/2 triangular wave,  WLtpar(1)=1/4 period of wave

     ! type of time behavior for AREA FLUX 1=step, 2=pulse, 3=stair, ...
     integer :: AreaTime

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), dimension(2) :: Atpar
     
     ! type of time behavior for BOUNDARY HEAD/FLUX 1=step, 2=pulse, 3=stair, ...
     integer :: BdryTime

     ! parameters related to different time behaviors (on, off, etc)
     real(DP), dimension(2) :: Btpar
  end type time

  ! parameters related to the entire domain (background)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  type, public :: domain

     ! vectors w/ bg & circular (elliptical too, later) properties
     real(DP),  allocatable :: kv(:), sv(:), av(:), porv(:), &
          & k2v(:), s2v(:),a2v(:), b2v(:), syv(:), kzv(:)
     integer, allocatable :: leakv(:), unconfv(:)
     real :: trash

     ! given hierarchy arrays (eventually to be calculated)
     integer, allocatable :: InclUp(:), WellUp(:)

     logical, allocatable :: InclIn(:,:), WellIn(:,:), WellBg(:,:), &
          & InclBg(:,:), CalcIn(:)

     ! number of each type of element
     ! 1=wells, 2=circles, 3=ellipses
     integer, dimension(3) :: num

     ! porosity, k, Ss, etc of background material
     real(DP), save :: por, k, Ss, k2, S2, b2, Sy, Kz, b
  end type domain
  
  type, private, extends(time) :: matching
     ! number of FS terms, number of matching points on circles
     ! for wells can be one (e.g., borehole storage) or zero (known Q)
     integer :: n,m

     ! tolerance to use in iterative solution for coefficients
     real(DP) :: matchTol

     ! SOR parameter for iterative solution for coefficients
     real(DP) :: matchOmega
  
     ! type of inclusion: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
     !                    -2=specified head ELEMENT, +2=specified flux ELEMENT
     integer :: ibnd

     !! leaky-related 
     real(DP) :: aquitardK, aquitardSs, aquitardb, Sy, Kz
     integer :: aquitardLeak , unconfined 

     ! whether inclusion is a matching(T) or specified(F) inclusion
     logical :: match

     ! specified value on bdry of inclusion
     real(DP) :: spec    

     ! porosity, starting time (step), and constant area flux for inclusion
     ! hydraulic conductivity and Ss inside matching element
     real(DP) :: por, area, k, Ss

     ! location of center of element
     real(DP) ::x, y 
  end type matching
    

  ! Circular Inclusion related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, extends(matching) :: circle
     ! inclusion radius
     real(kind=DP) :: r 
  end type circle
  
  ! Elliptical Inclusion related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, public, extends(matching) :: ellipse

     ! size of MF infinite matrix
     integer, save :: ms

     ! inclusion elliptical 'radius'
     ! semi-focal dist, angle major axis makes with Cartesian x axis,
     real(DP) :: eta, f, theta
  end type ellipse

  ! WeLl related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, public, extends(matching) :: well

     ! wellbore radius, pumping rate, dimensionless skin
     real(DP) :: r, q, dskin
     complex(DP) :: storCoeff
  end type well

  ! INVerse Laplace Transform parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, private :: laplace
     
     ! abcissa of convergence, LT tolerance
     real(kind=DP), save :: alpha, tol
     integer, save :: m
     logical, save :: smooth
  end type laplace

  ! things relating to the numerical solution, independent from flow elements
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type, public, extends(laplace) :: solution

     ! integrate particle lines vs. calculate at set locations/times
     logical :: particle
     integer :: numParticles

     ! calculate contours (thru space) vs. calculate hydrographs (thru time)
     logical :: contour
     ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
     logical :: calc
     
     ! input/output filenames
     character(128) :: outfname, infname, coefffname
     
     ! output index (1= gnuplot; 2= matlab)
     integer :: output, aquitardLeak, unconfined
     
     ! x-spacing vector, y-spacing vector, time vector
     integer :: numx, numy, numt
     real(DP) allocatable :: x(:), y(:), t(:)

     ! containers for time-domain results
     real(DP),    allocatable :: head(:,:,:), velx(:,:,:), vely(:,:,:)

     ! container for Laplace-space results
     complex(DP), allocatable :: headp(:), velxp(:), velyp(:)
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
  end type particle

end module element_specs
