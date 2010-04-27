 ! $Id: mathieu_line_const.f90,v 1.3 2010/02/24 21:31:05 klkuhlm Exp klkuhlm $
module mathieu_line_const
  use constants, only : DP
  implicit none

  ! this module is for the special case of constant strength line
  ! sources (both head and flux) using simplifications which only
  ! apply to the constant case

  private
  public :: mathieulinehead, cacosh, pointhead

  type, public :: line_calc
     ! the real/imaginary part of exp(-s)
     complex(DP) :: intp
     logical :: intreal, flux, match, calc_head, calc_velx

     ! parameters of the line source
     real(DP) :: q, alpha, Ss, k, sf, eta0, gamma, Ek
     real(DP) :: calcX, calcY
     complex(DP), allocatable :: calcCoeff(:,:,:)

     ! results from particle tracking
     real(DP),  allocatable :: Presult(:,:,:)
  end type line_calc
 
contains

  function mathieulinehead(LC,p) result(F)
    use mathieu_functions, only : mmatce, mmatse, mmatKe, mmatKo, mmatDKe, &
         & mmatDKo, mathieu_init, mathieu
    use constants, only : DP, TWOPI, PI, CZERO, RONE, PIOV2
    use utility, only : cacosh, linspace, circ => ellipse_circumference

    complex(DP), dimension(:), intent(in) :: p
    type(line_calc), intent(in) :: LC
    complex(DP), dimension(size(p)) :: F

    ! elliptical coords, semi-focal length
    real(DP) :: eta, psi
    complex(DP) :: W

    ! mathieu derived type for storing each set of Mathieu parameters, etc
    type(mathieu), save, dimension(size(p)) :: mathieu_vect
    type(mathieu) :: Z
    integer :: nP, i, k, j

    ! m: number of matching points
    ! n: number of terms in generalized fourier series
    ! ms: size of "infinite" matrix 
    !   used to estimate eigen{values,vectors} (calculated below)
    integer :: MS, N, M
    real(DP), allocatable :: vi(:), vs(:)

    ! coefficient matrix used in least-squares process for finding x (strengths of gFc)
    complex(DP), allocatable :: Am(:,:)
    complex(DP), allocatable :: bm(:)        ! RHS vector
    real(DP), allocatable :: Diag(:)     ! norm vector

    ! vectors of a_n and b_n for each value of p
    complex(DP), allocatable, save :: x(:,:)

    ! is this the first time through this subroutine?
    logical, save :: first = .true. 

    ! Mathieu coefficients (eigenvectors) for each value of p
    real(DP) :: hsq ! scale factors for coordinate system
    real(DP), allocatable, save :: psiM(:) ! vector of matching locations

    nP = size(p)

    ! acosh(Z/f) for complex argument
    W = cacosh(cmplx(LC%calcX,LC%calcY,DP)/LC%sf)

    ! radial and angular elliptical components
    eta = real(W)  ! radial
    psi = aimag(W) ! azimuthal

    ! mathieu parameter (McLachlan definition)
    ! NB: negative removed from definition !!
    ! (mathieu function re-defined accordingly)
    q = LC%sf**2*p(1:np)/(2.0_DP*LC%alpha)

    N = max(floor(3.0*maxval(abs(qq))),12)
    M = 2*(2*N + 1)
    allocate(vi(0:N))
    forall (j=0:N) vi(j) = real(j,DP)
    vs = 1.0_DP
    where(mod(vi,2)==1) vs = -1.0_DP

    if(first .and. LC%flux) then
       print *, 'specified flux element; q=',LC%q
    elseif(first .and. .not. LC%flux) then
       print *, 'specified head element; h=',LC%q
    end if
    
    if(first) then
       print *, 'N:',N,' M: ',M
       allocate(Am(1:M,0:2*N),bm(1:M),Diag(0:2*N),psiM(1:M))
    end if
    
    ! matching points along circumference of ellipse
    if (first) psiM(1:M) = linspace(-PI,+PI,M)
    
    ! some heuristic to estimate an appropriate matrix size
    MS = floor(1.5*maxval(abs(qq))) +N+25
    if (first) print *, 'N',N,'M',M,'MS',MS

    if (first) then
       ! coefficients for all values of p needed for inversion
       allocate(x(0:2*N,1:nP)) 
    end if

    ! for checking size of q vs size of matrix
    if (first) then
       print *, 'largest part of q:',maxval(abs(qq))
    end if

    do i = 1,np
       if (first) then
          print *, 'p:',p(i), 'i:',i, ' calculating MF coefficients...'
          ! calcualte eigenvectors for a single value of p         
          mathieu_vect(i) = mathieu_init(qq(i),MS)
       end if
       Z = mathieu_vect(i)

       if (LC%calc_head) then
          if (.not. LC%flux) then

             ! head effects of a constant strength
             ! specified head line element
             
             ! a_2 -> a_n (effectively stepping by twos)
             F(i) = TWOPI*LC%q*LC%k*sum(vs(0:N)*A(1:N+1,0,1)*mmatce(Z,2*vi(0:N),psi)* &
                  & mmatKe(Z,2*vi(0:N),eta)/mmatKe(Z,2*vi(0:N),LC%eta0))

          else
             
             ! head effects of a consant strength
             ! specified flux line element
             
             F(i) = LC%q/(4.0_DP*LC%sf)*mmatce(Z,0,psi)*mmatKe(Z,0,eta)/ &
                  & (mmatDKe(Z,0,LC%eta0)*mmatce(Z,0,PIOV2))

          endif
       else
          
          print *, 'flux effects not implemented yet'
       end if
    end do

    ! apply geometry scaling factors to all np terms at once
    if (.not. LC%calc_head) then
       hsq = LC%sf**2/2.0_DP*(cosh(2.0_DP*eta) - cos(2.0_DP*psi))
       if (LC%calc_velx) then
          F(1:np) = F(1:np)*sinh(eta)*cos(psi)/hsq
       else
          F(1:np) = F(1:np)*cosh(eta)*sin(psi)/hsq
       end if
    end if

    first = .false.
  end function mathieulinehead  

  ! an approximation to the constant line source using a sum of point sources
  function pointhead(LC,p) result(PH)
    use constants, only : DP, CZERO, TWOPI
    use complex_bessel, only : cbesk
    implicit none

    complex(DP), intent(in), dimension(:) :: p
    type(line_calc), intent(in) :: LC
    complex(DP), dimension(size(p)) :: PH
    complex(DP), dimension(size(p)) :: q, k0
    integer, parameter :: pts = 20 ! number of points to approximate line
    real(DP) :: r, begin, delxi, xi
    integer :: i, j, np, nz, ierr
    
    np = size(p)
    q(1:np) = sqrt(p(:)/LC%alpha)

    begin = -LC%sf
    delxi = 2.0_DP*LC%sf/real(pts-1,DP)

    PH = CZERO

    xi = begin
    do i=1,pts
       r = sqrt((calcX-xi)**2 + calcY**2)
       do j=1, np
          call cbesk(r*q(j),0.0_DP,1,1,k0(j),nz,ierr)
          if (nz /=0 .or. ierr /= 0) print *, 'bessel function error', nz, ierr
       end do
       PH(1:np) = PH + LC%q/(TWOPI*real(pts,DP)*p(:))*k0(:)
       xi = xi + delxi
    end do
  end function pointhead

end module mathieu_line_const
