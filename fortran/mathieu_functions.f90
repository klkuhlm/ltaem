module mathieu_functions
  implicit none

  ! $Id: mathieu_functions.f90,v 1.14 2007/09/12 20:22:28 kris Exp kris $

  !! only interfaces are publicly callable
  private
  public :: mmatce, mmatDce, mmatse, mmatDse, mmatKeKo, &
       & mmatDKeDKo, mmatIeIo, mmatDIeDIo

  ! even angular MF
  interface mmatce
     module procedure mmatce_scalar, mmatce_vect_z, mmatce_vect_n
  end interface

  ! derivative of even angular MF
  interface mmatDce
     module procedure mmatDce_scalar, mmatDce_vect_z, mmatDce_vect_n
  end interface

  ! odd angular MF
  interface mmatse
     module procedure mmatse_scalar, mmatse_vect_z, mmatse_vect_n
  end interface

  ! derivative of odd angular MF
  interface mmatDse
     module procedure mmatDse_scalar, mmatDse_vect_z, mmatDse_vect_n
  end interface

  ! even and odd second kind radial MF
  interface mmatKeKo
     module procedure mmatKeKo_scalar, mmatKeKo_vect_n
  end interface

  ! derivative of even and odd second kind radial MF
  interface mmatDKeDKo
     module procedure mmatDKeDKo_scalar, mmatDKeDko_vect_n
  end interface

  ! even and odd first kind radial MF
  interface mmatIeIo
     module procedure mmatIeIo_scalar, mmatIeIo_vect_n
  end interface

  ! derivative of even and odd first kind radial MF
  interface mmatDIeDIo
     module procedure mmatDIeDIo_scalar, mmatDIeDIo_vect_n
  end interface

  ! matrix size must be >= BUFFER + max(n)
  integer, parameter :: BUFFER = 5

contains

  !############################################################
  ! even angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatce_vect_n(n,z) result(ce)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)) :: ce

    ! internal variables
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, k, nm

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    do k=1,nm
       if (mod(n(k),2) == 0) then ! even order
          j = n(k)/2     ! a_{2m}(-q) = a_{2m}(q)
          vr(1:m) = real((/(2*i, i=0,m-1)/),DP)
          ce(k) = sgn(j)*sum(vi(1:m)*A(1:m,j,1)*cos(z*vr(1:m)))
       else                    ! odd order
          j = (n(k)-1)/2 ! a_{2m+1}(-q) = b_{2m+1}(q)
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          ce(k) = sgn(j)*sum(vi(1:m)*B(1:m,j,2)*cos(z*vr(1:m)))
       end if
    end do
 
  end function mmatce_vect_n

  !############################################################
  ! even angular modified mathieu function (q<0)
  ! for real vector of angular locations
  function mmatce_vect_z(n,z) result(ce)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP

    ! external arguments
    integer, intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(size(z)) :: ce

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, nz

    nz = size(z)
    m = size(A,dim=1)
    vi(1:m) = check_n(n,m)

    if (mod(n,2) == 0) then ! even order
       j = n/2
       vr(1:m) = real((/(2*i, i=0,m-1)/),DP)
       ce(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*A(1:m,j,1),dim=2,ncopies=nz)*&
            & cos(outerprod(vr(1:m),z(1:nz))),dim=1)
    else                    ! odd order
       j = (n-1)/2
       vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
       ce(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*B(1:m,j,2),dim=2,ncopies=nz)*&
            & cos(outerprod(vr(1:m),z(1:nz))),dim=1)
    end if
  end function mmatce_vect_z

  !############################################################
  ! derivative of even angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatDce_vect_n(n,z) result(Dce)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)) :: Dce

    ! internal variables
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, k ,nm

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    do k = 1,nm
       if (mod(n(k),2) == 0) then ! even order
          j = n(k)/2     ! a_{2m}(-q) = a_{2m}(q)
          vr(1:m) = real((/(2*i, i=0,m-1)/),DP)
          Dce(k) = -sgn(j)*sum(vi(1:m)*vr(1:m)*A(1:m,j,1)*sin(z*vr(1:m)))
       else                    ! odd order
          j = (n(k)-1)/2 ! a_{2m+1}(-q) = b_{2m+1}(q)
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          Dce(k) = -sgn(j)*sum(vi(1:m)*vr(1:m)*B(1:m,j,2)*sin(z*vr(1:m)))
       end if
    end do
 
  end function mmatDce_vect_n

  !############################################################
  ! derivative of even angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatDce_vect_z(n,z) result(Dce)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP

    ! external arguments
    integer, intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(size(z)) :: Dce

    ! internal variables
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, nz

    nz = size(z)
    m = size(A,dim=1)
    vi(1:m) = check_n(n,m)

    if (mod(n,2) == 0) then ! even order
       j = n/2     ! a_{2m}(-q) = a_{2m}(q)
       vr(1:m) = real((/(2*i, i=0,m-1)/),DP)
       Dce(1:nz) = -sgn(j)*sum(&
            & spread(vi(1:m)*vr(1:m)*A(1:m,j,1),dim=2,ncopies=nz)*&
            & sin(outerprod(vr(1:m),z(1:nz))),dim=1)
    else                    ! odd order
       j = (n-1)/2 ! a_{2m+1}(-q) = b_{2m+1}(q)
       vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
       Dce(1:nz) = -sgn(j)*sum(&
            & spread(vi(1:m)*vr(1:m)*B(1:m,j,2),dim=2,ncopies=nz)*&
            & sin(outerprod(vr(1:m),z(1:nz))),dim=1)
    end if
  end function mmatDce_vect_z

  !############################################################
  ! odd angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatse_vect_n(n,z) result(se)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP, CZERO

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)) :: se

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, k, nm

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    do k = 1,nm
       if (n(k) == 0) then
          se(k) = CZERO ! no se_0 defined, just return zero
       elseif (mod(n(k),2) == 0) then ! even order
          j = (n(k)-2)/2
          vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
          se(k) = sgn(j)*sum(vi(1:m)*B(1:m,j,1)*sin(z*vr(1:m)))
       else                        ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          se(k) = sgn(j)*sum(vi(1:m)*A(1:m,j,2)*sin(z*vr(1:m)))
       end if
    end do
 
  end function mmatse_vect_n

  !############################################################
  ! odd angular modified mathieu function (q<0)
  ! for real vector of angular locations
  function mmatse_vect_z(n,z) result(se)
    ! mathieu coefficients passed via module

    use shared_mathieu, only : A, B
    use constants, only : DP, CZERO

    ! external arguments
    integer, intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(size(z)) :: se

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, nz

    nz = size(z)
    m = size(A,dim=1)
    vi(1:m) = check_n(n,m)

    if (n == 0) then
       se(1:nz) = CZERO
    elseif (mod(n,2) == 0) then ! even order
       j = (n-2)/2
       vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
       se(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*B(1:m,j,1),dim=2,ncopies=nz)*&
            & sin(outerprod(vr(1:m),z(1:nz))),dim=1)
    else                    ! odd order
       j = (n-1)/2
       vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
       se(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*A(1:m,j,2),dim=2,ncopies=nz)*&
            & sin(outerprod(vr(1:m),z(1:nz))),dim=1)
    end if
  end function mmatse_vect_z

  !############################################################
  ! derivative of odd angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatDse_vect_n(n,z) result(Dse)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP, CZERO

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)) :: Dse

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, k, nm

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    do k = 1,nm
       if (n(k) == 0) then
          Dse(k) = CZERO ! no se_0 defined, just return zero
       elseif (mod(n(k),2) == 0) then ! even order
          j = (n(k)-2)/2
          vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
          Dse(k) = sgn(j)*sum(vi(1:m)*vr(1:m)*B(1:m,j,1)*cos(z*vr(1:m)))
       else                        ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          Dse(k) = sgn(j)*sum(vi(1:m)*vr(1:m)*A(1:m,j,2)*cos(z*vr(1:m)))
       end if
    end do
 
  end function mmatDse_vect_n

  !############################################################
  ! derivative of odd angular modified mathieu function (q<0)
  ! for scalar arguments
  function mmatDse_vect_z(n,z) result(Dse)
    ! mathieu coefficients passed via module
    use shared_mathieu, only : A, B
    use constants, only : DP, CZERO

    ! external arguments
    integer, intent(in) :: n
    real(DP), dimension(:), intent(in) :: z
    complex(DP), dimension(size(z)) :: Dse

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    integer :: i, j, m, nz

    nz = size(z)
    m = size(A,dim=1)
    vi(1:m) = check_n(n,m)

    if (n == 0) then
       Dse(1:nz) = CZERO ! no se_0 defined, just return zero
    elseif (mod(n,2) == 0) then ! even order
       j = (n-2)/2
       vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
       Dse(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*vr(1:m)*B(1:m,j,1),dim=2,ncopies=nz)*&
            & cos(outerprod(vr(1:m),z(1:nz))),dim=1)
    else                        ! odd order
       j = (n-1)/2
       vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
       Dse(1:nz) = sgn(j)*sum(&
            & spread(vi(1:m)*vr(1:m)*A(1:m,j,2),dim=2,ncopies=nz)*&
            & cos(outerprod(vr(1:m),z(1:nz))),dim=1)
    end if

  end function mmatDse_vect_z
  
  !############################################################
  ! even _and_ odd radial second kind modified mathieu functions (q<0)
  ! for a vector of orders (they both use same Bessel functions)
  subroutine mmatKeKo_vect_n(n,z,Ke,Ko)
    ! mathieu coefficients corresponding to q passed via module
    use shared_mathieu, only : A, B, q
    use constants, only : DP, PI, CZERO

    ! external arguments
    integer, dimension(:),intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)), intent(out) :: Ke, Ko
    
    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    complex(DP), dimension(0:size(A,dim=1)+1) :: besi, besk

    integer :: i, j, m, k, nm
    complex(DP) :: v1, v2, scale, p2n, p2n1, s2n1, s2n2

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    v1 = sqrt(q)*exp(-z)
    v2 = sqrt(q)*exp(+z)
    scale = exp(abs(real(v1)) - v2)
    call BesselI_val(v1,m+2,besi(0:m+1))
    call BesselK_val(v2,m+2,besk(0:m+1))

    do k = 1,nm
       if(mod(n(k),2) == 0) then ! even order
          j = n(k)/2
          ! pn2 := ce_2n(0,+q) * ce_2n(pi/2,+q)
          p2n = sum(A(1:m,j,1)) * sum(vi(1:m)*A(1:m,j,1))

          ! Ke even orders (Fek)
          Ke(k) = scale*sgn(j)*p2n/(PI*A(1,j,1)**2)* &
               & sum(A(:,j,1)*besi(0:m-1)*besk(0:m-1))

          if(n(k) == 0) then
             Ko(k) = CZERO
          else
             j = (n(k)-2)/2
             vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
             ! s2n2 := se'_2n+2(0,q) * se'_2n+2(pi/2,q)
             s2n2 = sum(B(1:m,j,1)*vr(1:m)) * sum(-vi(1:m)*B(1:m,j,1)*vr(1:m))

             ! Ko even orders (Gek)
             Ko(k) = scale*sgn(j+1)*s2n2/(PI*q*B(1,j,1)**2)*&
                  & sum(B(:,j,1)*(besi(0:m-1)*besk(2:m+1) - besi(2:m+1)*besk(0:m-1)))
          end if
       else ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          ! s2n1 := se'_2n+1(0,q) * se_2n+1(pi/2,q)
          s2n1 = sum(B(1:m,j,2)*vr(1:m)) * sum(vi(1:m)*B(1:m,j,2))

          ! Ke odd orders (Fek)
          Ke(k) = scale*sgn(j)*s2n1/(PI*sqrt(q)*B(1,j,2)**2)*&
               & sum(B(:,j,2)*(besi(0:m-1)*besk(1:m) - besi(1:m)*besk(0:m-1)))

          ! same definitions of j and vr here
          ! p2n1 := ce_2n+1(0,q) * ce'_2n+1(pi/2,q)
          p2n1 = sum(A(1:m,j,2)) * sum(-vi(1:m)*A(1:m,j,2)*vr(1:m))

          ! Ko odd orders (Gek)
          Ko(k) = scale*sgn(j+1)*p2n1/(PI*sqrt(q)*A(1,j,2)**2)* &
               & sum(A(:,j,2)*(besi(0:m-1)*besk(1:m) + besi(1:m)*besk(0:m-1)))
       end if
    end do
  end subroutine mmatKeKo_vect_n

  !############################################################
  ! derivative of even _and_ odd radial 
  ! second kind modified mathieu functions (q<0)
  ! for scalar arguments (they both use same Bessel functions)
  subroutine mmatDKeDKo_vect_n(n,z,DKe,DKo)
    ! mathieu coefficients corresponding to q passed via module
    use shared_mathieu, only : A, B, q
    use constants, only : DP, PI, CZERO

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)), intent(out) :: DKe, DKo
    
    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    complex(DP), dimension(0:size(A,dim=1)+1) :: besi, besk, Dbesi, Dbesk

    integer :: i, j, m, k, nm
    complex(DP) :: v1, v2, scale, p2n, p2n1, s2n1, s2n2, sqrtq, exppz, expnz
    
    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    sqrtq = sqrt(q)
    expnz = exp(-z)
    exppz = exp(+z)
    v1 = sqrtq*expnz
    v2 = sqrtq*exppz
    scale = exp(abs(real(v1)) - v2)
    call BesselI_val_and_deriv(v1,m+2,besi(0:m+1),Dbesi(0:m+1))
    call BesselK_val_and_deriv(v2,m+2,besk(0:m+1),Dbesk(0:m+1))

    do k = 1, nm
       if(mod(n(k),2) == 0) then ! even order
          j = n(k)/2
          ! pn2 := ce_2n(0,q) * ce_2n(pi/2,q)
          p2n = sum(A(1:m,j,1)) * sum(vi(1:m)*A(1:m,j,1))

          ! DKe even orders (Fek)
          DKe(k) = scale*sqrtq*sgn(j)*p2n/(PI*A(1,j,1)**2)* &
               & sum(A(:,j,1)*(exppz*besi(0:m-1)*Dbesk(0:m-1) - &
               &              expnz*Dbesi(0:m-1)*besk(0:m-1)))

          if(n(k) == 0) then
             DKo(k) = CZERO
          else
             j = (n(k)-2)/2
             vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
             ! s2n2 := se'_2n+2(0,q) * se'_2n+2(pi/2,q)
             s2n2 = sum(B(1:m,j,1)*vr(1:m)) * sum(-vi(1:m)*B(1:m,j,1)*vr(1:m))

             ! DKo even orders (Gek) - sqrt(q) cancels
             DKo(k) = scale*sgn(j+1)*s2n2/(PI*sqrtq*B(1,j,1)**2)* sum(B(1:m,j,1)*(&
                  & exppz*besi(0:m-1)*Dbesk(2:m+1) - expnz*Dbesi(0:m-1)*besk(2:m+1) - &
                  & exppz*besi(2:m+1)*Dbesk(0:m-1) + expnz*Dbesi(2:m+1)*besk(0:m-1)))
          end if
       else ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          ! s2n1 := se'_2n+1(0,q) * se_2n+1(pi/2,q)
          s2n1 = sum(B(1:m,j,2)*vr(1:m)) * sum(vi(1:m)*B(1:m,j,2))

          ! DKe odd orders (Fek) - sqrt(q) cancels
          DKe(k) = scale*sgn(j)*s2n1/(PI*B(1,j,2)**2)* sum(B(:,j,2)*(&
               & exppz*besi(0:m-1)*Dbesk(1:m) - expnz*Dbesi(0:m-1)*besk(1:m) - &
               & exppz*besi(1:m)*Dbesk(0:m-1) + expnz*Dbesi(1:m)*besk(0:m-1)))

          ! same definitions of j and vr here
          ! p2n1 := ce_2n+1(0,q) * ce'_2n+1(pi/2,q)
          p2n1 = sum(A(1:m,j,2)) * sum(-vi(1:m)*A(1:m,j,2)*vr(1:m))

          ! DKo odd orders (Gek) - sqrt(q) cancels
          DKo(k) = scale*sgn(j+1)*p2n1/(PI*A(1,j,2)**2)* sum(A(1:m,j,2)*(&
               & exppz*besi(0:m-1)*Dbesk(1:m) - expnz*Dbesi(0:m-1)*besk(1:m) + &
               & exppz*besi(1:m)*Dbesk(0:m-1) - expnz*Dbesi(1:m)*besk(0:m-1)))
       end if
    end do
  end subroutine mmatDKeDKo_vect_n

  !############################################################
   !! ---- based on Bessel function product series ----
  ! even _and_ odd radial first kind modified mathieu functions (q<0)
  ! for scalar arguments (they both use same Bessel functions)
  subroutine mmatIeIo_vect_n(n,z,Ie,Io)
    ! mathieu coefficients corresponding to q passed via module
    use shared_mathieu, only : A, B, q
    use constants, only : DP, CZERO

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)), intent(out) :: Ie, Io

    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    complex(DP), dimension(0:size(A,dim=1)+1) :: besiv1, besiv2

    integer :: i, j, m, k, nm
    complex(DP) :: v1, v2, scale, p2n, p2n1, s2n1, s2n2

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)

    v1 = sqrt(q)*exp(-z)
    v2 = sqrt(q)*exp(+z)
    scale = exp(abs(real(v1)) + abs(real(v2)))
    call BesselI_val(v1,m+2,besiv1(0:m+1))
    call BesselI_val(v2,m+2,besiv2(0:m+1))

    do k = 1, nm
       if(mod(n(k),2) == 0) then ! even order
          j = n(k)/2
          ! pn2 := ce_2n(0,q) * ce_2n(pi/2,q)
          p2n = sum(A(1:m,j,1))*sum(vi(1:m)*A(1:m,j,1))

          ! Ie even orders (Ce_2n)
          Ie(k) = sgn(j)*p2n/A(1,j,1)**2* &
               & sum(vi(1:m)*A(1:m,j,1)*besiv1(0:m-1)*besiv2(0:m-1))

          if(n(k) == 0) then
             Io(k) = CZERO
          else
             j = (n(k)-2)/2
             vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
             ! s2n2 := se'_2n+2(0,q) * se'_2n+2(pi/2,q)
             s2n2 = sum(B(1:m,j,1)*vr(1:m)) * sum(-vi(1:m)*B(1:m,j,1)*vr(1:m))

             ! Io even orders (Se_2n+2)
             Io(k) = sgn(j+1)*s2n2/(q*B(1,j,1)**2)*&
                  & sum(vi(1:m)*B(1:m,j,1)*(besiv1(0:m-1)*besiv2(2:m+1) - &
                  &                         besiv1(2:m+1)*besiv2(0:m-1)))
          end if
       else ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          ! s2n1 := se'_2n+1(0,q) * se_2n+1(pi/2,q)
          s2n1 = sum(B(1:m,j,2)*vr(1:m)) * sum(vi(1:m)*B(1:m,j,2)*vr(1:m))

          ! Ie odd orders (Ce_2n+1)
          Ie(k) = sgn(j)*s2n1/(sqrt(q)*B(1,j,2)**2)*&
               & sum(vi(1:m)*B(1:m,j,2)*(besiv1(0:m-1)*besiv2(1:m) + &
               &                         besiv1(1:m)*besiv2(0:m-1)))

          ! same definitions of j and vr here
          ! p2n1 := ce_2n+1(0,q) * ce'_2n+1(pi/2,q)
          p2n1 = sum(A(1:m,j,2)) * sum(-vi(1:m)*A(1:m,j,2)*vr(1:m))

          ! Io odd orders (Se_2n+1)
          Io(k) = sgn(j+1)*p2n1/(sqrt(q)*A(1,j,2)**2)* &
               & sum(vi(1:m)*A(1:m,j,2)*(besiv1(0:m-1)*besiv2(1:m) - &
               &                         besiv1(1:m)*besiv2(0:m-1)))
       end if
    end do
    
    Ie = scale*Ie
    Io = scale*Io

  end subroutine mmatIeIo_vect_n

  !############################################################
  !! ---- based on Bessel function product series ----
  ! derivative of even _and_ odd radial first kind modified mathieu functions (q<0)
  ! for scalar arguments (they both use same Bessel functions)
  subroutine mmatDIeDIo_vect_n(n,z,DIe,DIo)
    ! mathieu coefficients corresponding to q passed via module
    use shared_mathieu, only : A, B, q
    use constants, only : DP, CZERO

    ! external arguments
    integer, dimension(:), intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), dimension(size(n)), intent(out) :: DIe, DIo
    
    ! internal variable
    real(DP), dimension(size(A,dim=1)):: vr, vi
    complex(DP), dimension(0:size(A,dim=1)+1) :: besiv1, Dbesiv1, besiv2, Dbesiv2

    integer :: i, j, m, k, nm
    complex(DP) :: v1, v2, scale, p2n, p2n1, s2n1, s2n2, sqrtq, exppz, expnz

    nm = size(n)
    m = size(A,dim=1)
    vi(1:m) = check_n(maxval(n),m)
    sqrtq = sqrt(q)
    expnz = exp(-z)
    exppz = exp(+z)

    v1 = sqrtq*expnz
    v2 = sqrtq*exppz

    scale = exp(abs(real(v1)) + abs(real(v2)))
    call BesselI_val_and_deriv(v1,m+2,besiv1(0:m+1),Dbesiv1(0:m+1)) ! also -v1 on D
    call BesselI_val_and_deriv(v2,m+2,besiv2(0:m+1),Dbesiv2(0:m+1)) ! also +v2 on D

    do k = 1, nm
       if(mod(n(k),2) == 0) then ! even order
          j = n(k)/2
          ! pn2 := ce_2n(0,q) * ce_2n(pi/2,q)
          p2n = sum(A(1:m,j,1)) * sum(vi(1:m)*A(1:m,j,1))

          ! DIe even orders (Ce_2n)
          DIe(k) = sqrtq*sgn(j)*p2n/A(1,j,1)**2* &
               & sum(vi(1:m)*A(1:m,j,1)*(&
               & exppz*besiv1(0:m-1)*Dbesiv2(0:m-1) - &
               & expnz*Dbesiv1(0:m-1)*besiv2(0:m-1)))

          if(n(k) == 0) then
             DIo(k) = CZERO
          else
             j = (n(k)-2)/2
             vr(1:m) = real((/(2*i+2, i=0,m-1)/),DP)
             ! s2n2 := se'_2n+2(0,q) * se'_2n+2(pi/2,q)
             s2n2 = sum(B(1:m,j,1)*vr(1:m)) * sum(-vi(1:m)*B(1:m,j,1)*vr(1:m))

             ! DIo even orders (Se_2n+2)
             DIo(k) = sgn(j+1)*s2n2/(sqrtq*B(1,j,1)**2)*&
                  & sum(vi(1:m)*B(1:m,j,1)*(exppz*besiv1(0:m-1)*Dbesiv2(2:m+1) - &
                  & expnz*Dbesiv1(0:m-1)*besiv2(2:m+1) - &
                  & exppz*besiv1(2:m+1)*Dbesiv2(0:m-1) + &
                  & expnz*Dbesiv1(2:m+1)*besiv2(0:m-1)))
          end if
       else ! odd order
          j = (n(k)-1)/2
          vr(1:m) = real((/(2*i+1, i=0,m-1)/),DP)
          ! s2n1 := se'_2n+1(0,q) * se_2n+1(pi/2,q)
          s2n1 = sum(B(1:m,j,2)*vr(1:m)) * sum(vi(1:m)*B(1:m,j,2)*vr(1:m))

          ! DIe odd orders (Ce_2n+1)
          DIe(k) = sgn(j)*s2n1/B(1,j,2)**2 *&
               & sum(vi(1:m)*B(1:m,j,2)*(exppz*besiv1(0:m-1)*Dbesiv2(1:m) - &
               & expnz*Dbesiv1(0:m-1)*besiv2(1:m) + &
               & exppz*besiv1(1:m)*Dbesiv2(0:m-1) - &
               & expnz*Dbesiv1(1:m)*besiv2(0:m-1)))

          ! same definitions of j and vr here
          ! p2n1 := ce_2n+1(0,q) * ce'_2n+1(pi/2,q)
          p2n1 = sum(A(1:m,j,2)) * sum(-vi(1:m)*A(1:m,j,2)*vr(1:m))

          ! DIo odd orders (Se_2n+1)
          DIo(k) = sgn(j+1)*p2n1/A(1,j,2)**2 * &
               & sum(vi(1:m)*A(1:m,j,2)*(exppz*besiv1(0:m-1)*Dbesiv2(1:m) - &
               & expnz*Dbesiv1(0:m-1)*besiv2(1:m) - &
               & exppz*besiv1(1:m)*Dbesiv2(0:m-1) + &
               & expnz*Dbesiv1(1:m)*besiv2(0:m-1)))
       end if
    end do

    DIo = scale*DIo
    DIe = scale*DIe

  end subroutine mmatDIeDIo_vect_n


  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! wrapper functions allowing vector functions and subroutines 
  ! to be called for a scalar order 
  
  function mmatce_scalar(n,z) result(ce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP) :: ce
    ce = sum(mmatce_vect_n((/n/),z))
  end function mmatce_scalar

  function mmatDce_scalar(n,z) result(Dce)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP) :: Dce
    Dce = sum(mmatDce_vect_n((/n/),z))
  end function mmatDce_scalar

  function mmatse_scalar(n,z) result(se)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP) :: se
    se = sum(mmatse_vect_n((/n/),z))
  end function mmatse_scalar

  function mmatDse_scalar(n,z) result(Dse)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP) :: Dse
    Dse = sum(mmatDse_vect_n((/n/),z))
  end function mmatDse_scalar
  
  subroutine mmatKeKo_scalar(n,z,Ke,Ko)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), intent(out) :: Ke, Ko
    complex(DP), dimension(1) :: Ke_temp, Ko_temp
    call mmatKeKo_vect_n((/n/),z,Ke_temp(1),Ko_temp(1))
    Ke = Ke_temp(1); Ko = Ko_temp(1)
  end subroutine mmatKeKo_scalar

  subroutine mmatDKeDKo_scalar(n,z,DKe,DKo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), intent(out) :: DKe, DKo
    complex(DP), dimension(1) :: DKe_temp, DKo_temp
    call mmatDKeDKo_vect_n((/n/),z,DKe_temp(1),DKo_temp(1))
    DKe = DKe_temp(1); DKo = DKo_temp(1);
  end subroutine mmatDKeDKo_scalar

  subroutine mmatIeIo_scalar(n,z,Ie,Io)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), intent(out) :: Ie, Io
    complex(DP), dimension(1) :: Ie_temp, Io_temp    
    call mmatIeIo_vect_n((/n/),z,Ie_temp(1),Io_temp(1))
    Ie = Ie_temp(1); Io = Io_temp(1)
  end subroutine mmatIeIo_scalar

  subroutine mmatDIeDIo_scalar(n,z,DIe,DIo)
    use constants, only : DP
    integer, intent(in) :: n
    real(DP), intent(in) :: z
    complex(DP), intent(out) :: DIe, DIo
    complex(DP), dimension(1) :: DIe_temp, DIo_temp
    call mmatDIeDIo_vect_n((/n/),z,DIe_temp(1),DIo_temp(1))
    DIe = DIe_temp(1); DIo = DIo_temp(1)
  end subroutine mmatDIeDIo_scalar
 

  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! some utility routines
  
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! shorthand notation for real((-1)**i,8)
  elemental function sgn(i) result(x)
    use constants, only : DP, RONE
    integer, intent(in) :: i
    real(DP) :: x

    if (mod(i,2) == 0) then
       x = +RONE
    else
       x = -RONE
    end if
  end function sgn

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! shorthand notation for outer (tensor) product
  function outerprod(a,b) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: a,b
    real(DP), dimension(size(a),size(b)) :: c

    c = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
  end function outerprod

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! centralize error-checking, generate vi vector
  function check_n(n,ms) result(v) 
    use constants, only : DP
    integer, intent(in) :: n,ms
    real(DP), dimension(1:ms) :: v
    integer :: j

    if(ms < (n+BUFFER)) then
       write(*,'(3(A,I4))') "increase matrix &
            &size, n:",n,' m:',ms, ' buffer:',BUFFER
    end if

    v(1:ms) = sgn((/(j, j=0,ms-1)/))

  end function check_n

  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !! these subroutnies are wrappers for the complex Bessel funcitons
  !! implemented by Amos, Algorithm 644 TOMS, Vol 21, No 4, 1995
  !!
  !! these are SCALED results, un-scaling is done in the MF routines
  !!

  subroutine BesselI_val(z,n,I)
    use constants, only : DP, SMALL
    use complex_bessel, only : cbesi

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I
    integer :: nz, ierr

    ! scaling for I BF:: cy = I_fnu(z)*exp(-abs(x))
    ! where z = x + iy

    call cbesi(z=z, fnu=0.0_DP, kode=2, n=n, cy=I(0:n-1), nz=nz, ierr=ierr)
    
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
          write(*,*) "CBESI: loss of precision, z=",z
       case(4)
          write(*,*) "CBESI: overflow, z or order too &
               &large, z=",z," n=",n
          stop "CBESI: overflow, z or order too large"
       case(5)
          stop "CBESI: algorithm termination not met"
       end select 
    end if

  end subroutine BesselI_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  ! that is handled separately in the MF routines
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselI_val_and_deriv(z,n,I,ID)
    use constants, only : DP

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: I, ID

    call BesselI_val(z,n,I(0:n-1))

    ! middle
    ID(1:n-2) = 0.5_DP*(I(0:n-3) + I(2:n-1))

    ! low end
    ID(0) = I(1)

    ! high end
    ID(n-1) = I(n-2) - real(n-1,DP)/z*I(n-1)

  end subroutine BesselI_val_and_deriv

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val(z,n,K)
    use constants, only : DP, LARGE
    use complex_bessel, only : cbesk

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K
    integer :: nz, ierr

    ! scaling :: cy = K_fnu(z)*exp(z)
    call cbesk(z=z, fnu=0.0_DP, kode=2, n=n, cy=K(0:n-1), nz=nz, ierr=ierr)
    
    if (ierr /= 0) then
       select case(ierr)
       case(1)
          write(*,*) "CBESK: input error, z=",z," n=",n
          stop "CBESK: input error"
       case(2)
          write(*,*) "CBESK: overflow, z too small or order &
               &too large for unscaled output, z=",z," n=",n
          stop "CBESK: overflow, z too small or order too &
               &large for unscaled output"
       case(3)
          write(*,*) "CBESK: loss of precision, z=",z
       case(4)
          write(*,*) "CBESK: overflow, z too small or order &
               &too large, z=",z," n=",n
          stop "CBESK: overflow, z too small or order too large"
       case(5)
          stop "CBESK: algorithm termination not met"
       end select
    end if

  end subroutine BesselK_val

  ! use recurrance relationships for derivatives
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! this does not include the derivative of the argument
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine BesselK_val_and_deriv(z,n,K,KD)
    use constants, only : DP

    complex(DP), intent(in) :: z
    integer, intent(in) :: n
    complex(DP), intent(out), dimension(0:n-1) :: K, KD

    call BesselK_val(z,n,K(0:n-1))

    ! middle
    KD(1:n-2) = -0.5_DP*(K(0:n-3) + K(2:n-1))

    ! low end
    KD(0) = -K(1)

    ! high end
    KD(n-1) = -(K(n-2) + real(n-1,DP)/z*K(n-1))

  end subroutine BesselK_val_and_deriv

end module mathieu_functions
