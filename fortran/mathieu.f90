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

! this Fortran90 module is designed to not depend on the LT-AEM.  It
! is also not called from LT-AEM (it uses mathieu_functions2.f90).
! This module, the modified version of the Amos Bessel function
! library (cbessel.f90), and the fmathieu.pyf definition file should
! be used to create a python-callable complex-parameter modified
! Mathieu function library.
!
! f2py --fcompiler=gnu95 -llapack -lblas -c fmathieu.pyf cbessel.f90 mathieu.f90

module mf

  implicit none

  private  !! only interfaces and mathieu_init are publicly callable
  public :: ce, Dce, se, Dse, Ke, Ko, DKe, DKo, Ie, Io, DIe, DIo, mathieu_init

  integer, parameter, private :: DP = 8 !! selected_real_kind(p=15,r=300)
  real(DP), parameter, private :: PIOV2 = atan(1.0_DP)*2.0_DP

contains
  subroutine mathieu_init(q,M,mcn,A,B) 
    ! this subroutine computes the eigenvalues given at least a value for the
    ! Mathieu parameter (q), the norm and matrix size (M) are optional

    interface ! LAPACK eigenvalue/eigenfunction routine interface
       subroutine ZGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK, &
            & LWORK,RWORK,INFO)
         integer, intent(in) :: LDVR, LDVL, LDA, N, LWORK
         character(LEN=1), intent(in) :: JOBVL, JOBVR
         complex(KIND=8), intent(inout) :: A(N,N), WORK(LWORK)
         complex(KIND=8), intent(out) :: W(N), VL(LDVL,N), VR(LDVR,N)
         real(KIND=8), intent(inout) :: RWORK(33*N)
         integer, intent(out) :: INFO
       end subroutine ZGEEV
    end interface

    ! externally visible things
    complex(DP), intent(in) :: q
    integer, intent(in) :: M
    complex(DP), intent(out), dimension(4*M) :: mcn 
    complex(DP), intent(out), dimension(M,M,2) :: A, B

    ! just one matrix of recursion coefficients (used 4 times)
    complex(DP), allocatable :: coeff(:,:)

    ! parameters for lapack routine
    complex(DP), allocatable :: work(:) ! 1:lwork
    complex(DP), dimension(1) :: dc
    real(DP), allocatable :: rwork(:) ! 1:33*MS
    integer :: info, di, lwork
    real(DP), allocatable :: v(:), vi(:) ! 0:M

    !! used in stratton normalization only
    complex(DP), allocatable :: w(:) ! 1:M

    integer :: i, j

    allocate(v(0:M),vi(0:M))
    do concurrent (j = 0:M)
      v(j) = j
    end do
    where(mod([(j,j=0,M)],2) == 0)
       vi = 1.0_DP
    elsewhere
       vi = -1.0_DP
    end where

    ! A/B 1st dimension: subscript in McLachlan notation,
    ! i.e., the position in the infinite sum

    ! A/B 2nd dimension: superscript in McLachlan notation,
    ! i.e., the n in the order 2n+1 of the Mathieu function which it is associated with

    ! A/B 3rd dimension: 0(even) or 1(odd) cases of the second dimension

    lwork = 2*M+1
    allocate(coeff(M,M), rwork(33*M), work(lwork), w(M))

    di = 1 ! dummy integer for lapack
    dc(1) = cmplx(0,0,DP) ! dummy complex for lapack

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of even order
    ! De_{2n} in eq 3.12 of Stamnes & Spjelkavik,
    ! Pure and Applied Optics, 4(3), 251-262, 1995
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ! main diagonal/ r counting from 0:m-1 like McLachlan
    Coeff(:,:) = 0.0
    do concurrent(i = 1:M-1)
      Coeff(i+1,i+1) = cmplx((2*i)**2, 0, DP)
    end do

    ! off diagonals
    do concurrent(i = 1:m, j = 1:m, j == i+1 .or. j == i-1)
      Coeff(i,j) = q
    end do

    ! special case
    Coeff(2,1) = 2.0*q

    call ZGEEV(JOBVL='N', JOBVR='V',N=M, A=Coeff, LDA=M, W=mcn(1:m), &
         & VL=dc, LDVL=di, VR=A(1:m,1:m,1), LDVR=M, &
         & WORK=work, LWORK=lwork, RWORK=rwork, INFO=info)

    if (info /= 0) write (*,'(A,I0,A)') "ZGEEV ERROR ",info, &
         &" calculating even coefficients of even order"

    ! Morse norm  ce_2n(psi=0)
    w = sum(spread(vi(1:M),2,M)*A(:,:,1),dim=1)
    A(:,:,1) = A(:,:,1)/spread(w,1,M)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of odd order
    ! De_{2n+1} in eqn 3.14 of St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = 0.0
    Coeff(1,1) = 1.0_DP + q

    do concurrent (i = 1:m-1)
      Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    end do
    do concurrent (i = 1:m, j = 1:m, j == i+1 .or. j == i-1)
      Coeff(i,j) = q
    end do

    call ZGEEV(JOBVL='N', JOBVR='V', N=M, A=Coeff, LDA=M, W=mcn(m+1:2*m),&
         & VL=dc, LDVL=di, VR=A(1:m,1:m,2), LDVR=M, &
         & WORK=work, LWORK=lwork, RWORK=rwork, INFO=info)

    if (info /= 0) write(*,'(A,I0,A)') "ZGEEV ERROR ",info, &
         &" calculating even coefficients of odd order"

    ! se'_2n+1(psi=0)
    w = sum(spread((2*v(0:M-1)+1)*vi(0:M-1),2,M)*A(:,:,2),dim=1)
    A(:,:,2) = A(:,:,2)/spread(w,1,M)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of even order
    ! Do_{2n+2} in eq 3.16 of St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ! this one not shifted by one, since 2n+2 -> 2n
    !! (but starting from one, rather than zero)
    Coeff(:,:) = 0.0

    do concurrent (i = 1:m)
      Coeff(i,i) = cmplx((2*i)**2, 0, DP)
    end do
    
    do concurrent (i = 1:m, j = 1:m, j == i+1 .or. j == i-1)
      Coeff(i,j) = q
    end do

    call ZGEEV(JOBVL='N', JOBVR='V', N=M, A=Coeff, LDA=M, W=mcn(2*m+1:3*m),&
         & VL=dc, LDVL=di, VR=B(1:m,1:m,1),LDVR=M,&
         & WORK=work, LWORK=lwork, RWORK=rwork, INFO=info)

    if (info /= 0) write (*,'(A,I0,A)') "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of even order"

    ! ce_2n+2(psi=0)
    w = sum(spread((2*v(0:M-1)+2)*vi(0:M-1),2,M)*B(:,:,1),dim=1)
    B(:,:,1) = B(:,:,1)/spread(w,1,M)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of odd order
    ! Do_{2n+1} of eqn 3.18 in St&Sp
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = 0.0
    Coeff(1,1) = 1.0_DP - q

    do concurrent (i = 1:m-1)
      Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    end do
    do concurrent (i = 1:m, j = 1:m, j == i+1 .or. j == i-1)
      Coeff(i,j) = q
    end do
    
    call ZGEEV(JOBVL='N', JOBVR='V', N=M, A=Coeff, LDA=M, W=mcn(3*m+1:4*m),&
         & VL=dc, LDVL=di, VR=B(1:m,1:m,2), LDVR=M,&
         & WORK=work, LWORK=lwork, RWORK=rwork, INFO=info)

    if (info /= 0) write (*,'(A,I0,A)') "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of odd order"

    ! ce_2n+1(psi=0)
    w = sum(B(:,:,2)*spread(vi(0:M-1),2,M),dim=1)
    B(:,:,2) = B(:,:,2)/spread(w,1,m)

    deallocate(coeff,rwork,w,work)

  end subroutine mathieu_init

  !############################################################
  ! even angular modified mathieu function (q<0)
  ! for vector order and argument (retuns an outer-product type result)
  ! ce(q) is called Se(-q) by Blanch, or Qe(q) by Alhargan
  ! functions here use identities in 7.02 of Blanch's AMS#59 publication
  function ce(n,z,A,B) 

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: ce

    ! internal variables
    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)) :: v, vi
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    integer :: M, nz, nje, njo

    M = size(A,dim=1)
    nz = size(z)
    call angfcnsetup(M,n,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(n/2.0)+1

    ce(:,EV) = sum(spread(A(:,j(EV),1),2,nz)* &
         & spread(cos(outer(2*v,PIOV2-z)),3,nje),dim=1)

    ce(:,OD) = sum(spread(B(:,j(OD),2),2,nz)* &
         & spread(sin(outer(2*v+1,PIOV2-z)),3,njo),dim=1)

  end function ce

  !############################################################
  ! odd angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  ! se is called So(-q) by Blanch and Qo(q) by Alhargan
  function se(n,z,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: se

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    integer :: M, nz, nje, njo

    M = size(A,dim=1)
    nz = size(z)
    call angfcnsetup(M,n,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where (n == 0) j = 0
    j = j+1
    
    se(:,EV) = sum(spread(B(:,j(EV),1),2,nz)* &
         & spread(sin(outer(2*v+2,PIOV2-z)),3,nje),dim=1)

    se(:,OD) = sum(spread(A(:,j(OD),2),2,nz)* &
         & spread(cos(outer(2*v+1,PIOV2-z)),3,njo),dim=1)

    where (spread(n,1,nz) == 0) se = -huge(1.0)  ! se_0() is invalid
  end function se

  !############################################################
  ! derivative of even angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  function Dce(n,z,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Dce

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    integer :: M, nz, nje, njo

    M = size(A,dim=1)
    nz = size(z)
    call angfcnsetup(M,n,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(n/2.0)+1

    Dce(:,EV) = sum(spread(spread(2*v,2,nje)*A(:,j(EV),1),2,nz)* &
         & spread(sin(outer(2*v,PIOV2-z)),3,nje),dim=1)

    Dce(:,OD) = -sum(spread(spread(2*v+1,2,njo)*B(:,j(OD),2),2,nz)* &
         & spread(cos(outer(2*v+1,PIOV2-z)),3,njo),dim=1)

  end function Dce

  !############################################################
  ! derivative of odd angular modified mathieu function (q<0)
  ! for vector order and argument (returns an outer-product type result)
  function Dse(n,z,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Dse

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    integer :: M, nz, nje, njo

    M = size(A,dim=1)
    nz = size(z)
    call angfcnsetup(M,n,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where (n == 0) j = 0
    j = j+1

    Dse(:,EV) = -sum(spread(spread(2*v+2,2,nje)*B(:,j(EV),1),2,nz)* &
         & spread(cos(outer(2*v+2,PIOV2-z)),3,nje),dim=1)

    Dse(:,OD) = sum(spread(spread(2*v+1,2,njo)*A(:,j(OD),2),2,nz)* &
         & spread(sin(outer(2*v+1,PIOV2-z)),3,njo),dim=1)

    where (spread(n,1,nz) == 0) Dse = -huge(1.0)  ! se_0() is undefined
  end function Dse

  !############################################################
  ! even radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ie(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Ie

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1),size(z)) :: I1, I2
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    integer :: nje, njo, nz, M

    nz = size(z)
    M = size(A,dim=1)
    call radfcnsetup(M,q,n,z,v1,v2,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(n/2.0)+1
    call BesselI_val(v1,M+1,I1(0:M,1:nz))
    call BesselI_val(v2,M+1,I2(0:M,1:nz))

    Ie(:,EV) = sum(spread(spread(vi,2,nje)*A(:,j(EV),1),2,nz)* &
         & spread(I1(0:m-1,:)*I2(0:m-1,:),3,nje),dim=1)/spread(A(1,j(EV),1),1,nz)

    Ie(:,OD) = sum(spread(spread(vi,2,njo)*B(:,j(OD),2),2,nz)* &
         & spread(I1(0:m-1,:)*I2(1:m,:) + I1(1:m,:)*I2(0:m-1,:),3,njo),dim=1)/ &
         & spread(B(1,j(OD),2),1,nz)

    Ie = Ie*spread(exp(abs(real(v1)) + abs(real(v2))),2,size(n))
  end function Ie

  !############################################################
  ! odd radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Io(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Io

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1)+1,size(z)) :: I1, I2
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    integer :: nje, njo, nz, M

    nz = size(z)
    M = size(A,dim=1)
    call radfcnsetup(M,q,n,z,v1,v2,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where (n == 0) j = 0
    j = j+1
    call BesselI_val(v1,m+2,I1(0:m+1,:))
    call BesselI_val(v2,m+2,I2(0:m+1,:))

    Io(:,EV) = sum(spread(spread(vi,2,nje)*B(:,j(EV),1),2,nz)* &
         & spread(I1(0:m-1,:)*I2(2:m+1,:) - I1(2:m+1,:)*I2(0:m-1,:),3,nje),dim=1)/ &
         & spread(B(1,j(EV),1),1,nz)

    Io(:,OD) = sum(spread(spread(vi,2,njo)*A(:,j(OD),2),2,nz)* &
         & spread(I1(0:m-1,:)*I2(1:m,:) - I1(1:m,:)*I2(0:m-1,:),3,njo),dim=1)/ &
         & spread(A(1,j(OD),2),1,nz)

    Io = Io*spread(exp(abs(real(v1)) + abs(real(v2))),2,size(n))
  end function Io

  !############################################################
  ! even radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ke(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Ke

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1),size(z)) :: I, K
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radfcnsetup(M,q,n,z,v1,v2,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(real(n,DP)/2.0)+1
    call BesselI_val(v1,m+1,I(0:m,:))
    call BesselK_val(v2,m+1,K(0:m,:))

    Ke(:,EV) = sum(spread(A(:,j(EV),1),2,nz)* &
         & spread(I(0:m-1,:)*K(0:m-1,:),3,nje),dim=1)/spread(A(1,j(EV),1),1,nz)

    Ke(:,OD) = sum(spread(B(:,j(OD),2),2,nz)* &
         & spread(I(0:m-1,:)*K(1:m,:) - I(1:m,:)*K(0:m-1,:),3,njo),dim=1)/ &
         & spread(B(1,j(OD),2),1,nz)

    Ke = Ke*spread(exp(abs(real(v1)) - v2),2,size(n))
  end function Ke

  !############################################################
  ! odd radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function Ko(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: Ko

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1)+1,size(z)) :: I, K
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radfcnsetup(M,q,n,z,v1,v2,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where(n == 0) j = 0
    j = j+1
    call BesselI_val(v1,m+2,I(0:m+1,:))
    call BesselK_val(v2,m+2,K(0:m+1,:))

    Ko(:,EV) = sum(spread(B(:,j(EV),1),2,nz)* &
         & spread(I(0:m-1,:)*K(2:m+1,:) - I(2:m+1,:)*K(0:m-1,:),3,nje),dim=1)/ &
         & spread(B(1,j(EV),1),1,nz)

    Ko(:,OD) = sum(spread(A(:,j(OD),2),2,nz)* &
         & spread(I(0:m-1,:)*K(1:m,:) + I(1:m,:)*K(0:m-1,:),3,njo),dim=1)/ &
         & spread(A(1,j(OD),2),1,nz)

    Ko = Ko*spread(exp(abs(real(v1)) - v2),2,size(n))
    where (spread(n,1,nz) == 0) Ko = -huge(1.0)  ! Ko_0() is invalid
  end function Ko

  !############################################################
  ! derivative of even radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DIe(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: DIe

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1),size(z)) :: I1, I2, DI1, DI2
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(size(A,dim=1),size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radderivfcnsetup(M,q,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(n/2.0)+1
    call BesselI_val_and_deriv(v1,m+1,I1(0:m,:),DI1(0:m,:))
    call BesselI_val_and_deriv(v2,m+1,I2(0:m,:),DI2(0:m,:))

    DIe(:,EV) = sqrtq*sum(spread(spread(vi,2,nje)*A(:,j(EV),1),2,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(0:m-1,:) - enz*DI1(0:m-1,:)*I2(0:m-1,:),3,nje),dim=1)/ &
         & spread(A(1,j(EV),1),1,nz)

    DIe(:,OD) = sqrtq*sum(spread(spread(vi,2,njo)*B(:,j(OD),2),2,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(1:m,:) - enz*DI1(0:m-1,:)*I2(1:m,:) + &
         &        epz*I1(1:m,:)*DI2(0:m-1,:) - enz*DI1(1:m,:)*I2(0:m-1,:),3,njo),dim=1)/&
         & spread(B(1,j(OD),2),1,nz)

    DIe = DIe*spread(exp(abs(real(v1)) + abs(real(v2))),2,size(n))
  end function DIe

  !############################################################
  ! derivative of odd radial first kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DIo(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: DIo

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1)+1,size(z)) :: I1, I2, DI1, DI2
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(size(A,dim=1),size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radderivfcnsetup(M,q,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where (n == 0) j = 0
    j = j+1
    call BesselI_val_and_deriv(v1,m+2,I1(0:m+1,:),DI1(0:m+1,:))
    call BesselI_val_and_deriv(v2,m+2,I2(0:m+1,:),DI2(0:m+1,:))

    DIo(:,EV) = sqrtq*sum(spread(spread(vi,2,nje)*B(:,j(EV),1),2,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(2:m+1,:) - enz*DI1(0:m-1,:)*I2(2:m+1,:) - &
         &        epz*I1(2:m+1,:)*DI2(0:m-1,:) - enz*DI1(2:m+1,:)*I2(0:m-1,:),3,nje),dim=1)/&
         & spread(B(1,j(EV),1),1,nz)

    DIo(:,OD) = sqrtq*sum(spread(spread(vi,2,njo)*A(:,j(OD),2),2,nz)* &
         & spread(epz*I1(0:m-1,:)*DI2(1:m,:) - enz*DI1(0:m-1,:)*I2(1:m,:) - &
         &        epz*I1(1:m,:)*DI2(0:m-1,:) - enz*DI1(1:m,:)*I2(0:m-1,:),3,njo),dim=1)/&
         & spread(A(1,j(OD),2),1,nz)

    DIo = DIo*spread(exp(abs(real(v1)) + abs(real(v2))),2,size(n))
    where (spread(n,1,nz) == 0) DIo = -huge(1.0)  ! DIo_0() is invalid
  end function DIo

  !############################################################
  ! derivative of even radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DKe(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: DKe

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1),size(z)) :: I, K, DI, DK
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(size(A,dim=1),size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radderivfcnsetup(M,q,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor(n/2.0)+1
    call BesselI_val_and_deriv(v1,m+1,I(0:m,:),DI(0:m,:))
    call BesselK_val_and_deriv(v2,m+1,K(0:m,:),DK(0:m,:))

    DKe(:,EV) = sqrtq*sum(spread(A(:,j(EV),1),2,nz)* &
         & spread(epz*I(0:m-1,:)*DK(0:m-1,:) - enz*DI(0:m-1,:)*K(0:m-1,:),3,nje),dim=1)/&
         & spread(A(1,j(EV),1),1,nz)

    DKe(:,OD) = sqrtq*sum(spread(B(:,j(OD),2),2,nz)* &
         & spread(epz*I(0:m-1,:)*DK(1:m,:) - enz*DI(0:m-1,:)*K(1:m,:) - &
         &        epz*I(1:m,:)*DK(0:m-1,:) - epz*DI(1:m,:)*K(0:m-1,:),3,njo),dim=1)/ &
         & spread(B(1,j(OD),2),1,nz)

    DKe = DKe*spread(exp(abs(real(v1)) - v2),2,size(n))
  end function DKe

  !############################################################
  ! derivative of odd radial second kind modified mathieu functions (q<0)
  ! for vector order and argument
  function DKo(n,z,q,A,B) 
    integer, dimension(:), intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), intent(in) :: q
    complex(DP), dimension(:,:,:), intent(in) :: A,B
    complex(DP), dimension(size(z),size(n)) :: DKo

    integer, dimension(size(n)) :: j
    complex(DP), dimension(size(A,dim=1)):: v, vi
    complex(DP), dimension(0:size(A,dim=1)+1,size(z)) :: I, K, DI, DK
    integer, dimension(count(mod(n,2) == 0)) :: EV
    integer, dimension(count(mod(n,2) == 1)) :: OD
    complex(DP), dimension(size(z)) :: v1, v2
    complex(DP), dimension(size(A,dim=1),size(z)) :: enz, epz
    complex(DP) :: sqrtq
    integer :: nje, njo, nz, m

    nz = size(z)
    m = size(A,dim=1)
    call radderivfcnsetup(M,q,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    nje = size(EV)
    njo = size(OD)
    j = floor((n-1)/2.0)
    where (n == 0) j = 0
    j = j+1
    call BesselI_val_and_deriv(v1,m+2,I(0:m+1,:),DI(0:m+1,:))
    call BesselK_val_and_deriv(v2,m+2,K(0:m+1,:),DK(0:m+1,:))

    DKo(:,EV) = sqrtq*sum(spread(B(:,j(EV),1),2,nz)* &
         & spread(epz*I(0:m-1,:)*DK(2:m+1,:) - enz*DI(0:m-1,:)*K(2:m+1,:) - &
         &        epz*I(2:m+1,:)*DK(0:m-1,:) - epz*DI(2:m+1,:)*K(0:m-1,:),3,nje),dim=1)/&
         & spread(B(1,j(EV),1),1,nz)

    DKo(:,OD) = sqrtq*sum(spread(A(:,j(OD),2),2,nz)* &
         & spread(epz*I(0:m-1,:)*DK(1:m,:) - enz*DI(0:m-1,:)*K(1:m,:) + &
         &        epz*I(1:m,:)*DK(0:m-1,:) - enz*DI(1:m,:)*K(0:m-1,:),3,njo),dim=1)/&
         & spread(A(1,j(OD),2),1,nz)

    DKo = DKo*spread(exp(abs(real(v1)) - v2),2,size(n))
    where (spread(n,1,nz) == 0) DKo = -huge(1.0)  ! DKo_0() is invalid
  end function DKo

  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! some utility routines

  ! possibly disable this checking subroutine for extra speed??
  subroutine check_order(n)
    integer, dimension(:), intent(in) :: n
    integer :: nn

    nn = size(n)
    if (any(n < 0)) then
       write(*,*) 'CHECK_ORDER1: order of Mathieu functions must be >= 0',n
       stop 'CHECK_ORDER1: only non-negative Mathieu function orders'
    end if

    if (nn > 1) then
       ! since vectors of integer indexing vectors are used on the LHS of
       ! expressions, we must test for many-to-one conditions, since it is a no-no

       ! check integer orders are monotonically increasing
       if (any(n(2:nn) - n(1:nn-1) <= 0)) then
          write(*,*) 'CHECK_ORDER2: cannot have out-of order or repeated '//&
               & 'indices in vector order',n
          stop 'CHECK_ORDER2: eliminate out-of-order or repeated indices '//&
               &'to Mathieu functions'
       end if
    end if

  end subroutine check_order

  subroutine angfcnsetup(M,n,v,vi,EV,OD)
    integer, intent(in) :: M
    integer, intent(in), dimension(:) :: n
    complex(DP), intent(out), dimension(M) :: v,vi
    integer, intent(out), dimension(count(mod(n,2) == 0)) :: EV
    integer, intent(out), dimension(count(mod(n,2) == 1)) :: OD
    integer, dimension(M) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    call check_order(n)

    ! compute vector of integers counting up
    do concurrent (j = 0:M-1)
      i(j+1) = j
    end do
    v = cmplx(i,0,DP)

    ! compute the "sign" vector
    vi = cmplx(1,0,DP)
    where (mod(i,2) == 1) vi = cmplx(-1,0,DP)

    ! indexing vectors based on even/odd-ness of n
    do concurrent (j = 1:size(n))
      nn(j) = j
    end do

    EV = pack(nn,mod(n,2) == 0)
    OD = pack(nn,mod(n,2) == 1)

  end subroutine angfcnsetup

  subroutine radfcnsetup(M,q,n,z,v1,v2,v,vi,EV,OD)
    integer, intent(in) :: M
    complex(DP), intent(in) :: q
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in), dimension(:) :: z
    complex(DP), intent(out), dimension(M) :: v,vi
    complex(DP), intent(out), dimension(size(z)) :: v1,v2
    integer, intent(out), dimension(count(mod(n,2) == 0)) :: EV
    integer, intent(out), dimension(count(mod(n,2) == 1)) :: OD
    integer, dimension(M) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    call check_order(n)

    ! compute vector of integers counting up
    do concurrent (j = 0:M-1)
      i(j+1) = j
    end do
    v = cmplx(i,0,DP)

    ! compute the "sign" vector
    vi = cmplx(1,0,DP)
    where (mod(i,2) == 1) vi = cmplx(-1,0,DP)

    v1 = sqrt(q)*exp(-z)
    v2 = sqrt(q)*exp( z)

    ! indexing vectors based on even/odd-ness of n
    do concurrent (j = 1:size(n))
      nn(j) = j
    end do

    EV = pack(nn,mod(n,2) == 0)
    OD = pack(nn,mod(n,2) == 1)

  end subroutine radfcnsetup

  subroutine radderivfcnsetup(M,q,n,z,sqrtq,v1,v2,enz,epz,v,vi,EV,OD)
    integer, intent(in) :: M
    complex(DP), intent(in) :: q
    integer, intent(in), dimension(:) :: n
    real(DP), intent(in), dimension(:) :: z
    complex(DP), intent(out), dimension(M) :: v, vi
    complex(DP), intent(out) :: sqrtq
    complex(DP), intent(out), dimension(size(z)) :: v1, v2
    complex(DP), intent(out), dimension(M,size(z)) :: enz, epz
    integer, intent(out), dimension(count(mod(n,2) == 0)) :: EV
    integer, intent(out), dimension(count(mod(n,2) == 1)) :: OD
    integer, dimension(M) :: i
    integer, dimension(size(n)) :: nn
    integer :: j

    call check_order(n)

    ! compute vector of integers counting up
    do concurrent (j = 0:M-1)
      i(j+1) = j
    end do
    v = cmplx(i,0,DP)

    ! compute the "sign" vector
    vi = cmplx(1,0,DP)
    where (mod(i,2) == 1) vi = cmplx(-1,0,DP)

    enz = spread(exp(-z),dim=1,ncopies=M)
    epz = spread(exp( z),dim=1,ncopies=M)
    sqrtq = sqrt(q)
    v1 = enz(1,:)*sqrtq
    v2 = epz(1,:)*sqrtq

    ! indexing vectors based on even/odd-ness of n
    do concurrent (j = 1:size(n))
      nn(j) = j
    end do

    EV = pack(nn,mod(n,2) == 0)
    OD = pack(nn,mod(n,2) == 1)

  end subroutine radderivfcnsetup

  pure function outer(a,b) result(c)
    complex(DP), intent(in), dimension(:) :: a
    real(DP), intent(in), dimension(:) :: b
    complex(DP), dimension(size(a),size(b)) :: c
    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  end function outer

  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! these subroutnies are wrappers for the complex Bessel funcitons
  !! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995
  !!
  !! these are SCALED results, un-scaling is done in the MF routines

  subroutine BesselI_val(arg,n,I)
    use complex_bessel, only : cbesi
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: I
    integer :: numzero, ierr, j
    integer, parameter :: NPRINT = 5

!!$    character(33) :: fmt

    ! scaling for I BF:: cy = I_fnu(z)*exp(-abs(x))
    ! where z = x + iy
    ! only print up to the first 5 elements of the z vectors during error reporting

    do j = 1,size(arg)
       call cbesi(z=arg(j), fnu=0.0_DP, kode=2, n=n, cy=I(0:n-1,j), nz=numzero, ierr=ierr)

       if (ierr /= 0) then
          select case(ierr)
          case(1)
             write(*,*) "CBESI: input error, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
             stop "CBESI: input error"
          case(2)
             write(*,*) "CBESI: overflow, z or order too" //&
                  &"large for unscaled output, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
             stop "CBESI: overflow, z or order too large for unscaled output"
!!$          case(3)
!!$             fmt = '(A, (ES11.3E3,1X,ES11.3E3,3X),I0)'
!!$             write(fmt(4:4),'(I1)') min(ubound(arg,1),NPRINT)
!!$             write(*,fmt) "CBESI: loss of precision, z=",arg(1:min(ubound(arg,1),NPRINT)),numzero
          case(4)
             write(*,*) "CBESI: overflow, z or order too &
                  &large, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
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
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: I, ID

    ! the recurrence relationships assume at least three entries
    if (n < 3) stop 'mathieu_functions2.f90 : besseli_val_and_deriv, n must be > 3'
    call BesselI_val(arg,n,I(0:n-1,:))

    ID(1:n-2,:) = 0.5_DP*(I(0:n-3,:) + I(2:n-1,:)) ! middle
    ID(0,:) = I(1,:) ! low end
    ID(n-1,:) = I(n-2,:) - real(n-1,DP)/arg*I(n-1,:) ! high end

  end subroutine BesselI_val_and_deriv

  subroutine BesselK_val(arg,n,K)
    use complex_bessel, only : cbesk
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: K
    integer :: numzero, ierr, j
    integer, parameter :: NPRINT = 5

!!$    character(33) :: fmt

    ! scaling for K BF :: cy = K_fnu(z)*exp(z)
    do j = 1,size(arg)
       call cbesk(z=arg(j), fnu=0.0_DP, kode=2, n=n, cy=K(0:n-1,j), nz=numzero, ierr=ierr)

       if (ierr /= 0) then
          select case(ierr)
          case(1)
             write(*,*) "CBESK: input error, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
             stop "CBESK: input error"
          case(2)
             write(*,*) "CBESK: overflow, z too small or order " // &
                  &"too large for unscaled output, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
             stop "CBESK: overflow, z too small or order too &
                  &large for unscaled output"
!!$          case(3)
!!$             fmt = '(A, (ES11.3E3,1X,ES11.3E3,3X),I0)'
!!$             write(fmt(4:4),'(I1)') min(ubound(arg,1),NPRINT)
!!$             write(*,fmt) "CBESK: loss of precision, z=",arg(1:min(ubound(arg,1),NPRINT)),numzero
          case(4)
             write(*,*) "CBESK: overflow, z too small or order " //&
                  &"too large, z=",arg(1:min(ubound(arg,1),NPRINT))," n=",n
             stop "CBESK: overflow, z too small or order too large"
          case(5)
             stop "CBESK: algorithm termination not met"
          end select
       end if
    end do

  end subroutine BesselK_val

  ! use recurrence relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val_and_deriv(arg,n,K,KD)
    complex(DP), intent(in), dimension(:) :: arg
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1,size(arg)) :: K, KD

    ! the recurrence relationships assume at least three entries
    if (n < 3) stop 'mathieu_functions2.f90 : besselk_val_and_deriv, n must be > 3'
    call BesselK_val(arg,n,K(0:n-1,:))

    KD(1:n-2,:) = -0.5_DP*(K(0:n-3,:) + K(2:n-1,:)) ! middle
    KD(0,:) = -K(1,:) ! low end
    KD(n-1,:) = -(K(n-2,:) + real(n-1,DP)/arg*K(n-1,:)) ! high end

  end subroutine BesselK_val_and_deriv

end module mf
