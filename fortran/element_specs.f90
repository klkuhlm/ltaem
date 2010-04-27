! $Id: element_specs.f90,v 1.5 2007/09/21 06:47:09 kris Exp kris $

!##################################################
module element_specs
  use constants
  implicit none
  
  public
  
  ! vectors w/ bg & circular (elliptical too, later) properties
  real(kind=DP), save, allocatable :: kv(:), sv(:), av(:), porv(:), &
       & k2v(:), s2v(:),a2v(:), b2v(:), syv(:), kzv(:)
  integer, save, allocatable :: leakv(:), unconfv(:)

  real :: trash

  ! Circular Inclusion related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! number of circular inclusion
  integer, save :: CInum

  ! number of FS terms, number of matching points on circles
  integer, save :: CIn, CIm

  ! tolerance to use in iterative solution for coefficients
  real(kind=DP), save :: CImatchTol

  ! SOR parameter for iterative solution for coefficients
  real(kind=DP), save :: CImatchOmega
  
  ! type of inclusion: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
  !                    -2=specified head ELEMENT, +2=specified flux ELEMENT
  integer, save, allocatable :: CIibnd(:)

  !! leaky-related 
  real(kind=DP), save, allocatable :: CIaquitardK(:), CIaquitardSs(:),&
       & CIaquitardb(:), CISy(:), CIKz(:)
  integer, save, allocatable :: CIaquitardLeak(:) , CIunconfined(:) 

  ! whether inclusion is a matching(T) or specified(F) inclusion
  logical, save, allocatable :: CImatch(:)

  ! specified value on bdry of inclusion
  real(kind=DP), save, allocatable :: CIspec(:)    

  ! inclusion radius, x&y location of center, k and Ss of inclusion
  real(kind=DP), save, allocatable ::CIr(:), CIx(:), CIy(:), CIk(:), CIss(:)

  ! porosity, starting time (step), and constant area flux for inclusion
  real(kind=DP), save, allocatable :: CIpor(:), CIarea(:)

  ! type of time behavior for AREA FLUX 1=step, 2=pulse, 3=stair, ...
  integer, save, allocatable :: CIAreaTime(:)
  ! parameters related to different time behaviors (on, off, etc)
  real(kind=DP), save, allocatable :: CIAtpar(:,:)

  ! type of time behavior for BOUNDARY HEAD/FLUX 1=step, 2=pulse, 3=stair, ...
  integer, save, allocatable :: CIBdryTime(:)
  ! parameters related to different time behaviors (on, off, etc)
  real(kind=DP), save, allocatable :: CIBtpar(:,:)

  ! given hierarchy arrays (eventually to be calculated)
  integer, save, allocatable :: CIInclUp(:), CIWellUp(:)
  logical, save, allocatable :: CIInclIn(:,:), CIWellIn(:,:), CIWellBg(:,:), &
       & CIInclBg(:,:), CICalcIn(:)

  ! Elliptical Inclusion related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! number of elliptical elements
  integer, save :: EInum

  ! number of FS terms, number of matching points on circles
  ! size of MF infinite matrix
  integer, save :: EIn, EIm, EIms

  ! tolerance to use in iterative solution for coefficients
  real(kind=DP), save :: EImatchTol

  ! type of inclusion: -1=specified head TOTAL, 0=match, +1=specified flux TOTAL
  !                    -2=specified head ELEMENT, +2=specified flux ELEMENT
  integer, save, allocatable :: EIibnd(:)

  ! whether inclusion is a matching(T) or specified(F) inclusion
  logical, save, allocatable :: EImatch(:)

  ! specified value on bdry of inclusion
  real(kind=DP), save, allocatable :: EIspec(:)    

  ! inclusion radius, x&y location of center,
  ! semi-focal dist, angle major axis makes with Cartesian x axis,
  ! k and Ss of inclusion
  real(kind=DP), save, allocatable :: EIeta(:), EIf(:), EIx(:), EIy(:), EItheta(:)
  real(kind=DP), save, allocatable :: EIk(:), EIss(:)

  ! porosity, starting time (step), and constant area flux for inclusion
  real(kind=DP), save, allocatable :: EIpor(:), EIarea(:)

  ! type of time behavior for AREA FLUX 1=step, 2=pulse, 3=stair, ...
  integer, save, allocatable :: EIAreaTime(:)
  ! parameters related to different time behaviors (on, off, etc)
  real(kind=DP), save, allocatable :: EIAtpar(:,:)

  ! type of time behavior for BOUNDARY HEAD/FLUX 1=step, 2=pulse, 3=stair, ...
  integer, save, allocatable :: EIBdryTime(:)
  ! parameters related to different time behaviors (on, off, etc)
  real(kind=DP), save, allocatable :: EIBtpar(:,:)

  ! given hierarchy arrays (eventually to be calculated)
  integer, save, allocatable :: EIInclUp(:), EIWellUp(:)
  logical, save, allocatable :: EIInclIn(:,:), EIWellIn(:,:), EIWellBg(:,:), &
       & EIInclBg(:,:), EICalcIn(:)

  ! BackGround / GENeral parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! integrate particle lines vs. calculate at set locations/times
  logical, save :: BGparticle
  ! calculate contours (thru space) vs. calculate hydrographs (thru time)
  logical, save :: BGcontour
  ! re-calculate coefficient matrix, or load saved file (BGcoefffname)
  logical, save :: BGcalc

  ! input/output filenames
  character(128), save :: BGoutfname, BGinfname, BGcoefffname

  ! output index (1= gnuplot; 2= matlab)
  integer, save :: BGoutput, BGaquitardLeak, BGunconfined

  ! porosity, k and Ss of background material
  real(kind=DP), save :: BGpor, BGk, BGss, BGk2, BGS2, BGb2, BGSy, BGKz, bgb

  ! x-spacing vector, y-spacing vector, time vector
  integer, save :: BGnumx, BGnumy, BGnumt
  real(kind=DP), save, allocatable :: BGx(:), BGy(:), BGt(:)

  ! INVerse Laplace Transform parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! abcissa of convergence, LT tolerance
  real(kind=DP), save :: INValpha, INVtol
  integer, save :: INVm
  logical, save :: INVsmooth


  ! WeLl related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! number of wells
  integer, save :: WLnum

  ! x&y location, wellbore radius, pumping rate, dimensionless skin
  real(kind=DP), save, allocatable :: WLx(:), WLy(:), WLr(:), WLq(:), WLdskin(:)

  ! simulate wellbore storage in pumping well?
  logical, save, allocatable :: WLstor(:)
  complex(kind=DP), save, allocatable :: WLstorCoeff(:)

  ! time behavior / parameters (WLtpar())
  ! 1=step on,              WLtpar(1)  =on time
  ! 2=finite pulse,         WLtpar(1:2)=on/off time
  ! 3=instan. pulse,  -  no parameters
  ! 4=stairs,               WLtpar(1:2)=time step, off time
  ! 5=1/2 square wave,      WLtpar(1)  =1/2 period of wave
  ! 6=sine,                 WLtpar(1)=1/amplitude & wave number
  ! 7=1/2 triangular wave,  WLtpar(1)=1/4 period of wave
  integer, save, allocatable :: WLtime(:)
  real(kind=DP), save, allocatable :: WLtpar(:,:)

  
  ! PARticle related parameters
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! number of particles
  integer, save :: PARnum, PARstreakSkip

  ! starting x&y location, inital & final times
  real(kind=DP), save, allocatable :: PARx(:), PARy(:), PARti(:), PARtf(:)
  
  ! which integration scheme to use (for each particle??)
  ! 1 = Runge-Kutta-Merson (4th order adaptive)
  ! 2 = Runge-Kutta        (4th order)
  ! 3 = Richardson Extrapolation (2nd order)
  ! 4 = Fwd Euler          (1st order)
  integer, save, allocatable :: PARint(:)
  
  ! error tolerance and minimum stepsize for rkm
  real(kind=DP), save :: PARtol, PARmin, PARmaxStep

  ! step size for other integration schemes (initial stepsize for rkm)
  real(kind=DP), save :: PARdt

  ! particle starts inside a constant head or constant flux inclusion?
  logical, save, allocatable :: PARInclIn(:)

end module element_specs
