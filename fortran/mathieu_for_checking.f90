module mathieu_functions
  use constants, only : DP
  implicit none

  type, public :: mathieu
     integer :: M,  buffer
     real(DP) :: CUTOFF
     complex(DP) :: q
     complex(DP), allocatable :: mcn(:)  ! 4*M
     complex(DP), allocatable :: A(:,:,:), B(:,:,:)  ! 1:M, 0:M-1, 0:1
  end type mathieu

  ! $Id: mathieu_functions2.f90,v 1.18 2010/02/26 15:43:23 klkuhlm Exp klkuhlm $
  ! updated Feb 2010

  private  !! only interfaces and mathieu_init are publicly callable
  public :: ce, Dce, se, Dse, Ke, Ko, DKe, DKo, Ie, Io, DIe, DIo, mathieu_init

  interface ce   ! even first kind angular MF
     module procedure ce_scalar_nz, ce_scalar_n, ce_scalar_z, ce_vect_nz
  end interface

  interface Dce  ! derivative of even first kind angular MF
     module procedure Dce_scalar_nz, Dce_scalar_n, Dce_scalar_z, Dce_vect_nz
  end interface

  interface se  ! odd first kind angular MF
     module procedure se_scalar_nz, se_scalar_n, se_scalar_z, se_vect_nz
  end interface

  interface Dse  ! derivative of odd first kind angular MF
     module procedure Dse_scalar_nz, Dse_scalar_n, Dse_scalar_z, Dse_vect_nz
  end interface

  interface Ie  ! even first kind radial MF
     module procedure Ie_scalar_nz, Ie_scalar_n, Ie_scalar_z, Ie_vect_nz
  end interface

  interface Io  ! odd first kind radial MF
     module procedure Io_scalar_nz, Io_scalar_n, Io_scalar_z, Io_vect_nz
  end interface

  interface Ke  ! even second kind radial MF
     module procedure Ke_scalar_nz, Ke_scalar_n, Ke_scalar_z, Ke_vect_nz
  end interface

  interface Ko  ! odd second kind radial MF
     module procedure Ko_scalar_nz, Ko_scalar_n, Ko_scalar_z, Ko_vect_nz
  end interface

  interface DIe  ! derivative of even first kind radial MF
     module procedure DIe_scalar_nz, DIe_scalar_n, DIe_scalar_z, DIe_vect_nz
  end interface

  interface DIo  ! derivative of odd first kind radial MF
     module procedure DIo_scalar_nz, DIo_scalar_n, DIo_scalar_z, DIo_vect_nz
  end interface

  interface DKe  ! derivative of even second kind radial MF
     module procedure DKe_scalar_nz, DKe_scalar_n, DKe_scalar_z, DKe_vect_nz
  end interface

  interface DKo  ! derivative of odd second kind radial MF
     module procedure DKo_scalar_nz, DKo_scalar_n, DKo_scalar_z, DKo_vect_nz
  end interface
  
contains

  function mathieu_init(q,MM,CUTOFF) result(mat)
    ! this subroutine computes the eigenvalues given at least a value for the
    ! Mathieu parameter (q), the norm and matrix size (M) are optional

    use constants, only : DP, CZERO
    use utility, only : diagonal, signz

    interface ! LAPACK 3.2.1 eigenvalue/eigenfunction routine
       subroutine ZGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK, &
            & LWORK,RWORK,INFO)
         integer, intent(in) :: LDVR, LDVL, LDA, N, LWORK
         character(LEN=1), intent(in) :: JOBVL, JOBVR
         complex(KIND=8), intent(inout) :: A(N,N), WORK(LWORK)
         complex(KIND=8), intent(out) :: W(N), VL(LDVL,N), VR(LDVR,N)
         real(KIND=8), intent(inout) :: RWORK(2*N)
         integer, intent(out) :: INFO
       end subroutine ZGEEV
    end interface

    complex(DP), intent(in) :: q
    integer, optional, intent(in) :: MM
    real(DP), optional, intent(in) :: cutoff

    ! resulting structure defined at top of containing module
    type(mathieu) :: mat
    integer :: i, j, M

    real(DP) :: diag

    ! just one matrix of recursion coefficients (used 4 times)
    complex(DP), allocatable :: coeff(:,:)
    
    ! parameters for lapack routine
    !! lwork=570 Intel PentiumD 64-bit desktop
    !! lwork=330 Intel pro Mac 32-bit workstation
    !! lwork=330 Intel eagles 32-bit workstation
!!$    integer, parameter :: lwork = 330
    complex(DP), allocatable :: work(:)
    complex(DP), dimension(1) :: dc
    real(DP), allocatable :: rwork(:) ! 1:2*M
    integer :: info, di, lwork

    !! used in stratton normalization only
    complex(DP), allocatable :: w(:) ! 1:M

    mat%q = q

    ! hopefully a better heuristic for M will be figured out some day
    if (present(MM)) then
       mat%M = MM
    else
       ! useful? default
       mat%M = max(16,ceiling(2.0*abs(q)))
    end if
    M = mat%M

    if (present(CUTOFF)) then
       mat%CUTOFF = CUTOFF
    else
       mat%CUTOFF = 1.0D-12
    endif

    ! A/B 1st dimension: subscript in McLachlan notation,
    ! i.e., the position in the infinite sum
    
    ! A/B 2nd dimension: superscript in McLachlan notation,
    ! i.e., the n in the order 2n+1 of the Mathieu function which it is associated with
    
    ! A/B 3rd dimension: 0(even) or 1(odd) cases of the second dimension
    
    allocate(mat%mcn(4*M), mat%A(1:M,0:M-1,0:1), mat%B(1:M,0:M-1,0:1))
    allocate(coeff(M,M), rwork(2*M), w(M))

    di = 1 ! dummy integer for lapack
    dc(1) = CZERO ! dummy complex for lapack

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of even order
    ! De_{2n} in eq 3.12 of Stamnes & Spjelkavik, 
    ! Pure and Applied Optics, 4(3), 251-262, 1995
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO

    ! main diagonal/ r counting from 0:m-1 like McLachlan
    forall(i=1:M-1) 
       Coeff(i+1,i+1) = cmplx((2*i)**2, 0, DP)
    end forall
    
    ! off diagonals
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) 
       Coeff(i,j) = q
    end forall
        
    ! special case
    Coeff(2,1) = 2.0_DP*q

    ! compute optimium work-size array
    allocate(work(1000))
    call zgeev('N','V',M,Coeff,M,mat%mcn(1:m),dc,di,mat%A(1:m,0:m-1,0),M,&
         & work,-1,rwork,info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    call zgeev('N','V',M,Coeff,M,mat%mcn(1:m),dc,di,mat%A(1:m,0:m-1,0),M,&
         & work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of even order"

    ! Morse norm  ce_2n(psi=0)
    w = sum(mat%A(:,:,0),dim=1)
    
    ! normalize so |ce_2n(psi=0)| = +1
    mat%A(:,:,0) = mat%A(:,:,0)/spread(abs(w)*signz(w),1,m)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of odd order
    ! De_{2n+1} in eqn 3.14 of St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO
    Coeff(1,1) = 1.0_DP + q
    
    forall(i=1:m-1) 
       Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    end forall
    
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) 
       Coeff(i,j) = q    
    end forall
    
    call zgeev('N','V',M,Coeff,M,mat%mcn(m+1:2*m),dc,di,mat%A(1:m,0:m-1,1),M, &
         & work,lwork,rwork,info)
         
    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of odd order"

    ! se'_2n+1(psi=0)
    w = sum(spread(2*(\(i,i=0,m-1)/)+1,2,m)*mat%A(:,:,1),dim=1)
    
    ! normalize so |se'_2n+1(psi=0)|=+1
    mat%A(:,:,1) = mat%A(:,:,1)/spread(abs(w)*signz(w),1,m)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of even order
    ! Do_{2n+2} in eq 3.16 of St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO 
    
    ! this one not shifted by one, since 2n+2 -> 2n 
    !! (but starting from one, rather than zero)
    forall(i=1:m) 
       Coeff(i,i) = cmplx((2*i)**2, 0, DP)
    end forall
    
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) 
       Coeff(i,j) = q 
    end forall
    
    call zgeev('N','V',M,Coeff,M,mat%mcn(2*m+1:3*m),dc,di,mat%B(1:m,0:m-1,0),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of even order"

    ! ce_2n+2(psi=0)
    w = sum(spread(2*(\(i,i=1,m)/),dim=2,ncopies=m)*mat%B(:,:,0),dim=1)

    ! normalize so |ce_2n+2(psi=0)| = +1
    mat%B(:,:,0) = mat%B(:,:,0)/spread(abs(w)*signz(w),1,m)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of odd order
    ! Do_{2n+1} of eqn 3.18 in St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO
    Coeff(1,1) = 1.0_DP - q
    
    forall(i=1:m-1) 
       Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    end forall
    
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) 
       Coeff(i,j) = q 
    end forall
    
    call zgeev('N','V',M,Coeff,M,mat%mcn(3*m+1:4*m),dc,di,mat%B(1:m,0:m-1,1),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of odd order"

    ! ce_2n+1(psi=0)
    w = sum(mat%B(:,:,1),dim=1)
    
    mat%B(:,:,1) = mat%B(:,:,1)/spread(abs(w)*signz(w),1,m)

    deallocate(coeff,rwork,w)

    ! perform some heuristic checking of max allowable order
    diag = sum(abs(diagonal(mat%A(:,:,0),0)))/m
    mat%buffer = -999

    bufcheck: do j=1,m
       if    (sum(abs(diagonal(mat%A(:,:,0), j)))/(diag*(m-j)) < mat%CUTOFF .and. &
            & sum(abs(diagonal(mat%A(:,:,0),-j)))/(diag*(m-j)) < mat%CUTOFF) then
           mat%buffer = j
           exit bufcheck
        end if
        if (j == m) then
           write(*,*) 'MATHIEU_INIT: matrix not large enough ',m,&
                & ' to meet tolerance ', mat%CUTOFF,' for given q', mat%q
           stop 'increase M or decrease CUTOFF in mathieu_init'
        end if
     end do bufcheck
     
  end function mathieu_init

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !############################################################
  ! even angular modified mathieu function (q<0)
  ! for vector order and argument (retuns an outer-product type result)
  ! ce(q) is called Se(-q) by Blanch, or Qe(q) by Alhargan
  ! functions here use identities in 7.02 of Blanch's AMS#59 publication
  function ce_vect_nz(mf,n,z) result(ce)
    use constants, only : DP, PIOV2
    use utility, only : outerprod

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: ce

    ! internal variables
    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M) :: v, vi
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    integer :: nz, nje, njo

    nz = size(z)
    call angfcnsetup(mf,n,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD); j = floor(n/2.0)

    ce(EV,:) = sum(spread(mf%A(:,j(EV),0),3,nz)* &
         & spread(cos(outerprod(2*v,PIOV2-z)),2,nje),dim=1)/ &
         & spread(sum(spread(vi,2,nje)*mf%A(:,j(EV),0),dim=1),2,nz)

    ce(OD,:) = sum(spread(mf%B(:,j(OD),1),3,nz)* &
         & spread(sin(outerprod(2*v+1,PIOV2-z)),2,njo),dim=1)/ &
         & spread(sum(spread(vi,2,njo)*mf%B(:,j(OD),1),dim=1),2,nz)

  end function ce_vect_nz

  !############################################################
  ! odd angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  ! se is called So(-q) by Blanch and Qo(q) by Alhargan
  function se_vect_nz(mf,n,z) result(se)
    use constants, only : DP, PIOV2
    use utility, only : outerprod
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: se

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    integer :: nz, nje, njo

    nz = size(z)
    call angfcnsetup(mf,n,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)

    se(EV,:) = sum(spread(mf%B(:,j(EV),0),3,nz)* &
         & spread(sin(outerprod(2*v+2,PIOV2-z)),2,nje),dim=1)/ &
         & spread(sum(spread(vi*(2*v+2),2,nje)*mf%B(:,j(EV),0),dim=1),2,nz)

    se(OD,:) = sum(spread(mf%A(:,j(OD),1),3,nz)* &
         & spread(cos(outerprod(2*v+1,PIOV2-z)),2,njo),dim=1)/ &
         & spread(sum(spread(vi*(2*v+1),2,njo)*mf%A(:,j(OD),1),dim=1),2,nz)

    where (spread(n,2,nz) == 0) se = -huge(1.0)  ! se_0() is invalid
  end function se_vect_nz

  !############################################################
  ! derivative of even angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  function Dce_vect_nz(mf,n,z) result(Dce)
    use constants, only : DP, PIOV2
    use utility, only : outerprod
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Dce

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    integer :: nz, nje, njo

    nz = size(z)
    call angfcnsetup(mf,n,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor(n/2.0)

    Dce(EV,:) = sum(spread(spread(2*v,2,nje)*mf%A(:,j(EV),0),3,nz)* &
         & spread(sin(outerprod(2*v,PIOV2-z)),2,nje),dim=1)/ &
         & spread(sum(spread(vi,2,nje)*mf%A(:,j(EV),0),dim=1),2,nz)

    Dce(OD,:) = -sum(spread(spread(2*v+1,2,njo)*mf%B(:,j(OD),1),3,nz)* &
         & spread(cos(outerprod(2*v+1,PIOV2-z)),2,njo),dim=1)/ &
         & spread(sum(spread(vi,2,njo)*mf%B(:,j(OD),1),dim=1),2,nz)

  end function Dce_vect_nz

  !############################################################
  ! derivative of odd angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  function Dse_vect_nz(mf,n,z) result(Dse)
    use constants, only : DP, PIOV2
    use utility, only : outerprod
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Dse

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    integer :: nz, nje, njo

    nz = size(z)
    call angfcnsetup(mf,n,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)

    Dse(EV,:) = -sum(spread(spread(2*v+2,2,nje)*mf%B(:,j(EV),0),3,nz)* &
         & spread(cos(outerprod(2*v+2,PIOV2-z)),2,nje),dim=1)/ &
         & spread(sum(spread(vi*(2*v+2),2,nje)*mf%B(:,j(EV),0),dim=1),2,nz)

    Dse(OD,:) = sum(spread(spread(2*v+1,2,njo)*mf%A(:,j(OD),1),3,nz)* &
         & spread(sin(outerprod(2*v+1,PIOV2-z)),2,njo),dim=1)/ &
         & spread(sum(spread(vi*(2*v+1),2,njo)*mf%A(:,j(OD),1),dim=1),2,nz)

    where (spread(n,2,nz) == 0) Dse = -huge(1.0)  ! se_0() is undefined
  end function Dse_vect_nz

  !############################################################
  ! even radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ie_vect_nz(mf,n,z) result(Ie)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Ie

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+1,size(z)) :: I1, I2
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, M

    nz = size(z); M = mf%M
    call radfcnsetup(mf,n,z,sqrtq,v1,v2,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor(n/2.0)
    call BesselI_val(v1,M+2,I1(0:M+1,1:nz))
    call BesselI_val(v2,M+2,I2(0:M+1,1:nz))

    Ie(EV,:) = sum(spread(spread(vi,2,nje)*mf%A(:,j(EV),0),3,nz)* &
         & spread(I1(0:m-1,:)*I2(0:m-1,:),2,nje),dim=1)/spread(mf%A(1,j(EV),0),2,nz)

    Ie(OD,:) = sum(spread(spread(vi,2,njo)*mf%B(:,j(OD),1),3,nz)* &
         & spread(I1(0:m-1,:)*I2(1:m,:) + I1(1:m,:)*I2(0:m-1,:),2,njo),dim=1)/ &
         & spread(mf%B(1,j(OD),1),2,nz)

    Ie = Ie*spread(exp(abs(real(v1)) + abs(real(v2))),1,size(n))
  end function Ie_vect_nz

  !############################################################
  ! odd radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Io_vect_nz(mf,n,z) result(Io)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Io

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+2,size(z)) :: I1, I2
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, M

    nz = size(z); M = mf%M
    call radfcnsetup(mf,n,z,sqrtq,v1,v2,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)
    call BesselI_val(v1,m+3,I1(0:m+2,:))
    call BesselI_val(v2,m+3,I2(0:m+2,:))

    Io(EV,:) = sum(spread(spread(vi,2,nje)*mf%B(:,j(EV),0),3,nz)* &
         & spread(I1(0:m-1,:)*I2(2:m+1,:) - I1(2:m+1,:)*I2(0:m-1,:),2,nje),dim=1)/ &
         & spread(mf%B(1,j(EV),0),2,nz)

    Io(OD,:) = sum(spread(spread(vi,2,njo)*mf%A(:,j(OD),1),3,nz)* &
         & spread(I1(0:m-1,:)*I2(1:m,:) - I1(1:m,:)*I2(0:m-1,:),2,njo),dim=1)/ &
         & spread(mf%A(1,j(OD),1),2,nz)

    Io = Io*spread(exp(abs(real(v1)) + abs(real(v2))),1,size(n))
  end function Io_vect_nz

  !############################################################
  ! even radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ke_vect_nz(mf,n,z) result(Ke)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Ke

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+1,size(z)) :: I, K
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radfcnsetup(mf,n,z,sqrtq,v1,v2,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor(n/2.0)
    call BesselI_val(v1,m+2,I(0:m+1,:))
    call BesselK_val(v2,m+2,K(0:m+1,:))

    Ke(EV,:) = sum(spread(mf%A(:,j(EV),0),3,nz)* &
         & spread(I(0:m-1,:)*K(0:m-1,:),2,nje),dim=1)/spread(mf%A(1,j(EV),0),2,nz)

    Ke(OD,:) = sum(spread(mf%B(:,j(OD),1),3,nz)* &
         & spread(I(0:m-1,:)*K(1:m,:) - I(1:m,:)*K(0:m-1,:),2,njo),dim=1)/ &
         & spread(mf%B(1,j(OD),1),2,nz)

    Ke = Ke*spread(exp(abs(real(v1)) - v2),1,size(n))
  end function Ke_vect_nz

  !############################################################
  ! odd radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ko_vect_nz(mf,n,z) result(Ko)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: Ko

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+2,size(z)) :: I, K
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radfcnsetup(mf,n,z,sqrtq,v1,v2,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)
    call BesselI_val(v1,m+3,I(0:m+2,:))
    call BesselK_val(v2,m+3,K(0:m+2,:))

    Ko(EV,:) = sum(spread(mf%B(:,j(EV),0),3,nz)* &
         & spread(I(0:m-1,:)*K(2:m+1,:) - I(2:m+1,:)*K(0:m-1,:),2,nje),dim=1)/ &
         & spread(mf%B(1,j(EV),0),2,nz)

    Ko(OD,:) = sum(spread(mf%A(:,j(OD),1),3,nz)* &
         & spread(I(0:m-1,:)*K(1:m,:) - I(1:m,:)*K(0:m-1,:),2,njo),dim=1)/ &
         & spread(mf%A(1,j(OD),1),2,nz)

    Ko = Ko*spread(exp(abs(real(v1)) - v2),1,size(n))
    where (spread(n,2,nz) == 0) Ko = -huge(1.0)  ! Ko_0() is invalid
  end function Ko_vect_nz

  !############################################################
  ! derivative of even radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DIe_vect_nz(mf,n,z) result(DIe)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: DIe

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+1,size(z)) :: I1, I2, DI1, DI2
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(mf%M,size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radderivfcnsetup(mf,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor(n/2.0)
    call BesselI_val_and_deriv(v1,m+2,I1(0:m+1,:),DI1(0:m+1,:))
    call BesselI_val_and_deriv(v2,m+2,I2(0:m+1,:),DI2(0:m+1,:))

    DIe(EV,:) = spread(sqrtq/mf%A(1,j(EV),0),2,nz)* &
         & sum(spread(spread(vi,2,nje)*mf%A(:,j(EV),0),3,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(0:m-1,:) - enz*DI1(0:m-1,:)*I2(0:m-1,:),2,nje),dim=1)

    DIe(OD,:) = spread(sqrtq/mf%B(1,j(OD),1),2,nz)* &
         & sum(spread(spread(vi,2,njo)*mf%B(:,j(OD),1),3,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(1:m,:) - enz*DI1(0:m-1,:)*I2(1:m,:) + &
         &        epz*I1(1:m,:)*DI2(0:m-1,:) - enz*DI1(1:m,:)*I2(0:m-1,:),2,njo),dim=1)

    DIe = DIe*spread(exp(abs(real(v1)) + abs(real(v2))),1,size(n))
  end function DIe_vect_nz

  !############################################################
  ! derivative of odd radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DIo_vect_nz(mf,n,z) result(DIo)
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: DIo

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+2,size(z)) :: I1, I2, DI1, DI2
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(mf%M,size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radderivfcnsetup(mf,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)
    call BesselI_val_and_deriv(v1,m+3,I1(0:m+2,:),DI1(0:m+2,:))
    call BesselI_val_and_deriv(v2,m+3,I2(0:m+2,:),DI2(0:m+2,:))

    DIo(EV,:) = spread(sqrtq/mf%B(1,j(EV),0),2,nz)* &
         & sum(spread(spread(vi,2,nje)*mf%B(:,j(EV),0),3,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(2:m+1,:) - enz*DI1(0:m-1,:)*I2(2:m+1,:) - &
         &        epz*I1(2:m+1,:)*DI2(0:m-1,:) - enz*DI1(2:m+1,:)*I2(0:m-1,:),2,nje),dim=1)

    DIo(OD,:) = spread(sqrtq/mf%A(1,j(OD),1),2,nz)* &
         & sum(spread(spread(vi,2,njo)*mf%A(:,j(OD),1),3,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(1:m,:) - enz*DI1(0:m-1,:)*I2(1:m,:) - &
         &        epz*I1(1:m,:)*DI2(0:m-1,:) - enz*DI1(1:m,:)*I2(0:m-1,:),2,njo),dim=1)

    DIo = DIo*spread(exp(abs(real(v1)) + abs(real(v2))),1,size(n))
    where (spread(n,2,nz) == 0) DIo = -huge(1.0)  ! DIo_0() is invalid
  end function DIo_vect_nz

  !############################################################
  ! derivative of even radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DKe_vect_nz(mf,n,z) result(DKe)
    ! mathieu coefficients corresponding to q passed via module
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: DKe

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+1,size(z)) :: I, K, DI, DK
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(mf%M,size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radderivfcnsetup(mf,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor(n/2.0)
    call BesselI_val_and_deriv(v1,m+2,I(0:m+1,:),DI(0:m+1,:))
    call BesselK_val_and_deriv(v2,m+2,K(0:m+1,:),DK(0:m+1,:))

    DKe(EV,:) = spread(sqrtq/mf%A(1,j(EV),0),2,nz)*sum(spread(mf%A(:,j(EV),0),3,nz)* &
         & spread(epz*I(0:m-1,:)*DK(0:m-1,:) - enz*DI(0:m-1,:)*K(0:m-1,:),2,nje),dim=1)

    DKe(OD,:) = spread(sqrtq/mf%B(1,j(OD),1),2,nz)*sum(spread(mf%B(:,j(OD),1),3,nz)* &
         & spread(epz*I(0:m-1,:)*DK(1:m,:) - enz*DI(0:m-1,:)*K(1:m,:) - &
         &        epz*I(1:m,:)*DK(0:m-1,:) - epz*DI(1:m,:)*K(0:m-1,:),2,njo),dim=1)

    DKe = DKe*spread(exp(abs(real(v1)) - v2),1,size(n))
  end function DKe_vect_nz

  !############################################################
  ! derivative of odd radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DKo_vect_nz(mf,n,z) result(DKo)
    ! mathieu coefficients corresponding to q passed via module
    use constants, only : DP
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n),size(z)) :: DKo

    integer, dimension(size(n)) :: j
    real(DP), dimension(mf%M):: v, vi
    complex(DP), dimension(0:mf%M+2,size(z)) :: I, K, DI, DK
    integer, dimension(count(mod(n,2)==0)) :: EV
    integer, dimension(count(mod(n,2)==1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(mf%M,size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z); m = mf%M
    call radderivfcnsetup(mf,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV);  njo = size(OD);  j = floor((n-1)/2.0)
    call BesselI_val_and_deriv(v1,m+3,I(0:m+2,:),DI(0:m+2,:))
    call BesselK_val_and_deriv(v2,m+3,K(0:m+2,:),DK(0:m+2,:))

    DKo(EV,:) = spread(sqrtq/mf%B(1,j(EV),0),2,nz)*sum(spread(mf%B(:,j(EV),0),3,nz)* &
         & spread(epz*I(0:m-1,:)*DK(2:m+1,:) - enz*DI(0:m-1,:)*K(2:m+1,:) - &
         &        epz*I(2:m+1,:)*DK(0:m-1,:) - epz*DI(2:m+1,:)*K(0:m-1,:),2,nje),dim=1)

    DKo(OD,:) = spread(sqrtq/mf%A(1,j(OD),1),2,nz)*sum(spread(mf%A(:,j(OD),1),3,nz)* &
         & spread(epz*I(0:m-1,:)*DK(1:m,:) - enz*DI(0:m-1,:)*K(1:m,:) + &
         &        epz*I(1:m,:)*DK(0:m-1,:) - enz*DI(1:m,:)*K(0:m-1,:),2,njo),dim=1)

    DKo = DKo*spread(exp(abs(real(v1)) - v2),1,size(n))
    where (spread(n,2,nz) == 0) DKo = -huge(1.0)  ! DKo_0() is invalid
  end function DKo_vect_nz

  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! wrapper functions allowing vector functions and subroutines
  ! to be called for a scalar order or arguments or both

  function ce_scalar_n(mf,n,z) result(ce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: ce
    ce = sum(ce_vect_nz(mf,(\n/),z),dim=1)
  end function ce_scalar_n
  function ce_scalar_z(mf,n,z) result(ce)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: ce
    ce = sum(ce_vect_nz(mf,n,(\z/)),dim=2)
  end function ce_scalar_z
  function ce_scalar_nz(mf,n,z) result(ce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: ce
    ce = sum(ce_vect_nz(mf,(\n/),(\z/)))
  end function ce_scalar_nz

  function se_scalar_n(mf,n,z) result(se)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: se
    se = sum(se_vect_nz(mf,(\n/),z),dim=1)
  end function se_scalar_n
  function se_scalar_z(mf,n,z) result(se)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: se
    se = sum(se_vect_nz(mf,n,(\z/)),dim=2)
  end function se_scalar_z
  function se_scalar_nz(mf,n,z) result(se)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: se
    se = sum(se_vect_nz(mf,(\n/),(\z/)))
  end function se_scalar_nz

  function Dce_scalar_n(mf,n,z) result(Dce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Dce
    Dce = sum(Dce_vect_nz(mf,(\n/),z),dim=1)
  end function Dce_scalar_n
  function Dce_scalar_z(mf,n,z) result(Dce)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Dce
    Dce = sum(Dce_vect_nz(mf,n,(\z/)),dim=2)
  end function Dce_scalar_z
  function Dce_scalar_nz(mf,n,z) result(Dce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Dce
    Dce = sum(Dce_vect_nz(mf,(\n/),(\z/)))
  end function Dce_scalar_nz

  function Dse_scalar_n(mf,n,z) result(Dse)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Dse
    Dse = sum(Dse_vect_nz(mf,(\n/),z),dim=1)
  end function Dse_scalar_n
  function Dse_scalar_z(mf,n,z) result(Dse)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Dse
    Dse = sum(Dse_vect_nz(mf,n,(\z/)),dim=2)
  end function Dse_scalar_z
  function Dse_scalar_nz(mf,n,z) result(Dse)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Dse
    Dse = sum(Dse_vect_nz(mf,(\n/),(\z/)))
  end function Dse_scalar_nz

  function Ie_scalar_n(mf,n,z) result(Ie)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Ie
    Ie = sum(Ie_vect_nz(mf,(\n/),z),dim=1)
  end function Ie_scalar_n
  function Ie_scalar_z(mf,n,z) result(Ie)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Ie
    Ie = sum(Ie_vect_nz(mf,n,(\z/)),dim=2)
  end function Ie_scalar_z
  function Ie_scalar_nz(mf,n,z) result(Ie)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Ie
    Ie = sum(Ie_vect_nz(mf,(\n/),(\z/)))
  end function Ie_scalar_nz

  function Io_scalar_n(mf,n,z) result(Io)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Io
    Io = sum(Io_vect_nz(mf,(\n/),z),dim=1)
  end function Io_scalar_n
  function Io_scalar_z(mf,n,z) result(Io)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Io
    Io = sum(Io_vect_nz(mf,n,(\z/)),dim=2)
  end function Io_scalar_z
  function Io_scalar_nz(mf,n,z) result(Io)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Io
    Io = sum(Io_vect_nz(mf,(\n/),(\z/)))
  end function Io_scalar_nz

  function Ke_scalar_n(mf,n,z) result(Ke)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Ke
    Ke = sum(Ke_vect_nz(mf,(\n/),z),dim=1)
  end function Ke_scalar_n
  function Ke_scalar_z(mf,n,z) result(Ke)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Ke
    Ke = sum(Ke_vect_nz(mf,n,(\z/)),dim=2)
  end function Ke_scalar_z
  function Ke_scalar_nz(mf,n,z) result(Ke)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Ke
    Ke = sum(Ke_vect_nz(mf,(\n/),(\z/)))
  end function Ke_scalar_nz

  function Ko_scalar_n(mf,n,z) result(Ko)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: Ko
    Ko = sum(Ko_vect_nz(mf,(\n/),z),dim=1)
  end function Ko_scalar_n
  function Ko_scalar_z(mf,n,z) result(Ko)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: Ko
    Ko = sum(Ko_vect_nz(mf,n,(\z/)),dim=2)
  end function Ko_scalar_z
  function Ko_scalar_nz(mf,n,z) result(Ko)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: Ko
    Ko = sum(Ko_vect_nz(mf,(\n/),(\z/)))
  end function Ko_scalar_nz

  function DIe_scalar_n(mf,n,z) result(DIe)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: DIe
    DIe = sum(DIe_vect_nz(mf,(\n/),z),dim=1)
  end function DIe_scalar_n
  function DIe_scalar_z(mf,n,z) result(DIe)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: DIe
    DIe = sum(DIe_vect_nz(mf,n,(/z/)),dim=2)
  end function DIe_scalar_z
  function DIe_scalar_nz(mf,n,z) result(DIe)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: DIe
    DIe = sum(DIe_vect_nz(mf,(/n/),(/z/)))
  end function DIe_scalar_nz

  function DIo_scalar_n(mf,n,z) result(DIo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: DIo
    DIo = sum(DIo_vect_nz(mf,(/n/),z),dim=1)
  end function DIo_scalar_n
  function DIo_scalar_z(mf,n,z) result(DIo)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: DIo
    DIo = sum(DIo_vect_nz(mf,n,(/z/)),dim=2)
  end function DIo_scalar_z
  function DIo_scalar_nz(mf,n,z) result(DIo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: DIo
    DIo = sum(DIo_vect_nz(mf,(/n/),(/z/)))
  end function DIo_scalar_nz

  function DKe_scalar_n(mf,n,z) result(DKe)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: DKe
    DKe = sum(DKe_vect_nz(mf,(/n/),z),dim=1)
  end function DKe_scalar_n
  function DKe_scalar_z(mf,n,z) result(DKe)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: DKe
    DKe = sum(DKe_vect_nz(mf,n,(/z/)),dim=2)
  end function DKe_scalar_z
  function DKe_scalar_nz(mf,n,z) result(DKe)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: DKe
    DKe = sum(DKe_vect_nz(mf,(/n/),(/z/)))
  end function DKe_scalar_nz

  function DKo_scalar_n(mf,n,z) result(DKo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(z)) :: DKo
    DKo = sum(DKo_vect_nz(mf,(/n/),z),dim=1)
  end function DKo_scalar_n
  function DKo_scalar_z(mf,n,z) result(DKo)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP), dimension(size(n)) :: DKo
    DKo = sum(DKo_vect_nz(mf,n,(/z/)),dim=2)
  end function DKo_scalar_z
  function DKo_scalar_nz(mf,n,z) result(DKo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    type(mathieu), intent(in) :: mf
    complex(DP) :: DKo
    DKo = sum(DKo_vect_nz(mf,(/n/),(/z/)))
  end function DKo_scalar_nz

  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! some utility routines

  function check_cutoff(mf,n) result(err)
    type(mathieu), intent(in) :: mf
    integer, dimension(:), intent(in) :: n
    integer :: err

    if (maxval(n,1) > mf%m - mf%buffer) then
       write(*,*) 'CHECK_CUTOFF: max order requested ',maxval(n,1), &
            & ' too large for q ',mf%q, ' given current CUTOFF ',mf%CUTOFF, &
            & ' since buffer =',mf%buffer, 'and M=',mf%M
       stop 'CHECK_CUTOFF: increase M or decrease CUTOFF'
    else
       err = 0
    end if
  end function check_cutoff

  subroutine angfcnsetup(mf,n,v,vi,EV,OD)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    type(mathieu), intent(in) :: mf
    real(DP), intent(out), dimension(mf%m) :: v,vi
    integer, intent(out), dimension(count(mod(n,2)==0)) :: EV
    integer, intent(out), dimension(count(mod(n,2)==1)) :: OD
    integer, dimension(mf%m) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    j = check_cutoff(mf,n)
    ! compute vector of integers counting up
    forall (j=0:mf%m-1) i(j+1) = j
    v = real(i,DP) 

    ! compute the "sign" vector
    vi = 1.0_DP
    where (mod(i,2)==1)
       vi = -1.0_DP
    end where

    ! indexing vectors based on even/odd-ness of n
    forall (j=1:size(n)) nn(j) = j
    EV = pack(nn,mod(n,2)==0)
    OD = pack(nn,mod(n,2)==1)
  end subroutine angfcnsetup

  subroutine radfcnsetup(mf,n,z,sqrtq,v1,v2,v,vi,EV,OD)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    real(DP), intent(out), dimension(mf%m) :: v,vi
    complex(DP), intent(out) :: sqrtq
    complex(DP), intent(out), dimension(size(z)) :: v1,v2
    integer, intent(out), dimension(count(mod(n,2)==0)) :: EV
    integer, intent(out), dimension(count(mod(n,2)==1)) :: OD
    integer, dimension(mf%m) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    j = check_cutoff(mf,n)
    ! compute vector of integers counting up
    forall (j=0:mf%m-1) i(j+1) = j
    v = real(i,DP)

    ! compute the "sign" vector
    vi = 1.0_DP
    where (mod(i,2)==1)
       vi = -1.0_DP
    end where

    ! indexing vectors based on even/odd-ness of n
    forall (j=1:size(n)) nn(j) = j
    EV = pack(nn,mod(n,2)==0)
    OD = pack(nn,mod(n,2)==1)
       
    sqrtq = sqrt(mf%q)
    v1 = sqrtq*exp(-z)
    v2 = sqrtq*exp(z)
  end subroutine radfcnsetup

  subroutine radderivfcnsetup(mf,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    use constants, only : DP
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in), dimension(:) :: z
    type(mathieu), intent(in) :: mf
    real(DP), intent(out), dimension(mf%m) :: v, vi
    complex(DP), intent(out) :: sqrtq
    complex(DP), intent(out), dimension(size(z)) :: v1, v2
    complex(DP), intent(out), dimension(mf%m,size(z)) :: enz, epz
    integer, intent(out), dimension(count(mod(n,2)==0)) :: EV
    integer, intent(out), dimension(count(mod(n,2)==1)) :: OD
    integer, dimension(mf%m) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    j = check_cutoff(mf,n)
    ! compute vector of integers counting up
    forall (j=0:mf%m-1) i(j+1) = j
    v = real(i,DP)

    ! compute the "sign" vector
    vi = 1.0_DP
    where (mod(i,2)==1)
       vi = -1.0_DP
    end where

    ! indexing vectors based on even/odd-ness of n
    forall (j=1:size(n)) nn(j) = j
    EV = pack(nn,mod(n,2)==0)
    OD = pack(nn,mod(n,2)==1)

    enz = spread(exp(-z),dim=1,ncopies=mf%m)
    epz = spread(exp( z),dim=1,ncopies=mf%m)
    sqrtq = sqrt(mf%q)
    v1 = enz(1,:)*sqrtq
    v2 = epz(1,:)*sqrtq
  end subroutine radderivfcnsetup

  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! these subroutnies are wrappers for the complex Bessel funcitons
  !! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995
  !!
  !! these are SCALED results, un-scaling is done in the MF routines

  subroutine BesselI_val(z,n,I)
    use constants, only : DP, SMALL
    use complex_bessel, only : cbesi
    complex(DP), intent(in), dimension(:) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(z)) :: I
    integer :: numzero, ierr, j

    ! scaling for I BF:: cy = I_fnu(z)*exp(-abs(x))
    ! where z = x + iy

    do j=1,size(z)
       call cbesi(z=z(j), fnu=0.0_DP, kode=2, n=n, cy=I(0:n-1,j), nz=numzero, ierr=ierr)

       if (ierr /= 0) then
          select case(ierr)
          case(1)
             write(*,*) "CBESI: input error, z=",z," n=",n
             stop "CBESI: input error"
          case(2)
             write(*,*) "CBESI: overflow, z or order too &
                  &large for unscaled output, z=",z," n=",n
             stop "CBESI: overflow, z or order too large &
                  &for unscaled output"
          case(3)
             write(*,*) "CBESI: loss of precision, z=",z,numzero
          case(4)
             write(*,*) "CBESI: overflow, z or order too &
                  &large, z=",z," n=",n
             stop "CBESI: overflow, z or order too large"
          case(5)
             stop "CBESI: algorithm termination not met"
          end select
       end if
    end do
  end subroutine BesselI_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  ! that is handled separately in the MF routines
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselI_val_and_deriv(arg,n,I,ID)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: I, ID

    call BesselI_val(arg,n,I(0:n-1,:))
    ID(1:n-2,:) = 0.5_DP*(I(0:n-3,:) + I(2:n-1,:)) ! middle
    ID(0,:) = I(1,:) ! low end
    ID(n-1,:) = I(n-2,:) - real(n-1,DP)/arg*I(n-1,:) ! high end
  end subroutine BesselI_val_and_deriv

  subroutine BesselK_val(arg,n,K)
    use constants, only : DP, LARGE
    use complex_bessel, only : cbesk
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: K
    integer :: numzero, ierr, j

    ! scaling for K BF :: cy = K_fnu(z)*exp(z)

    do j=1,size(arg)
       call cbesk(z=arg(j), fnu=0.0_DP, kode=2, n=n, cy=K(0:n-1,j), nz=numzero, ierr=ierr)

       if (ierr /= 0) then
          select case(ierr)
          case(1)
             write(*,*) "CBESK: input error, z=",arg," n=",n
             stop "CBESK: input error"
          case(2)
             write(*,*) "CBESK: overflow, z too small or order &
                  &too large for unscaled output, z=",arg," n=",n
             stop "CBESK: overflow, z too small or order too &
                  &large for unscaled output"
          case(3)
             write(*,*) "CBESK: loss of precision, z=",arg,numzero
          case(4)
             write(*,*) "CBESK: overflow, z too small or order &
                  &too large, z=",arg," n=",n
             stop "CBESK: overflow, z too small or order too large"
          case(5)
             stop "CBESK: algorithm termination not met"
          end select
       end if
    end do

  end subroutine BesselK_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val_and_deriv(arg,n,K,KD)
    use constants, only : DP
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: K, KD

    call BesselK_val(arg,n,K(0:n-1,:))
    KD(1:n-2,:) = -0.5_DP*(K(0:n-3,:) + K(2:n-1,:)) ! middle
    KD(0,:) = -K(1,:) ! low end
    KD(n-1,:) = -(K(n-2,:) + real(n-1,DP)/arg*K(n-1,:)) ! high end
  end subroutine BesselK_val_and_deriv

end module mathieu_functions
