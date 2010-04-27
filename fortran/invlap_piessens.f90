module invlap_piessens
! $Id: invlap_piessens.f90,v 1.9 2007/09/18 19:28:44 kris Exp kris $

private
public invlap, pvalues

contains

function pvalues(b,c,N) result(p)
  use constants, only : DP, PI
  implicit none

  real(DP), intent(in) :: b,c
  integer, intent(in) :: N
  real(DP), dimension(0:N-1) :: p,z

  real(DP), dimension(0:N-1) :: iv 
  integer :: j

  iv(0:N-1) = real((/(j, j=0,N-1)/),DP)
  z(0:N-1) = cos((iv(0:N-1) + 0.5_DP)*PI/real(N,DP))
  p(0:N-1) = 2.0_DP*b/(1.0_DP - z(0:N-1)) + c   

end function pvalues

function invlap(b,c,p,fp,t,N) result(ft)
  use constants, only : DP, PI
  implicit none

  real(DP), parameter :: S = 1.0_DP  !! not changing this now
  real(DP), intent(in) :: b,c
  integer, intent(in) :: N
  real(DP), intent(in), dimension(0:N-1) :: p,fp
  real(DP), intent(in), dimension(:) :: t
  real(DP), dimension(size(t)) :: ft
  real(DP), dimension(0:N-1) :: phi, z

  real(DP), dimension(0:N-1) :: a, vi
  integer :: i,nt

  intrinsic :: dgamma

  !! check parameters
  if(N <= 3) then
     print *, 'N=',N
     stop 'N must be > 3'
  end if
!!$  if(s < 1.0_DP) then
!!$     print *, 's=',s
!!$     stop 's must be >= 1'
!!$  end if
  if(b < c) then
     print *, 'b=',b,' c=',c
     stop 'b must be >= abcissa of convergence'
  end if
  
  vi(0:N-1) = real((/(i, i=0,N-1)/),DP)
  nt = size(t)
  
  !! compute Chebyshev coefficients
  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ! z is independent of i
  z(0:N-1) = cos((vi(0:N-1) + 0.5_DP)*PI/real(N,DP))

  ! zero term is halved and cos(0)=1
  a(0) = sum( fp(0:N-1)*(2.0_DP*b/(1.0_DP - z(0:N-1)))**s )/real(N,DP)

  do i=1,N-1
     a(i) = 2.0_DP/real(N,DP)*sum( fp(0:N-1)*&
          & (2.0_DP*b/(1.0_DP - z(0:N-1)))**s *&
          & cos((vi(0:N-1) + 0.5_DP)*real(i,DP)*PI/real(N,DP)) )
  end do
 
  !! compute solution
  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  do i=1,nt
     phi(0:N-1) = hypergeom22(b*t(i),N-1,s)
     ft(i) = exp(c*t(i))*t(i)**(s - 1.0_DP)/dgamma(s)*&
          & sum(a(0:N-1)*phi(0:N-1))
  end do

end function invlap


function hypergeom22(x,N,s) result(phi)  
  use constants, only : DP
  implicit none

  real(DP), intent(in) :: x,s
  integer, intent(in) :: N
  real(DP), dimension(0:N) :: phi

  real(DP), dimension(3:N) :: A,B,C,D,E, vi
  integer :: j

  vi(3:N) = real((/(j, j=3,N)/),DP)

  !! recurrence taken from Davies, 2002 (Chap 19, p 336), other versions of this
  !! (Davies, 1979 + Piessens, 1972) seem to have typos in them

  ! compute coefficients vector-wise
  A(3:N) = (3.0_DP*vi**2 - 9.0_DP*vi + s*vi - 3.0_DP*s + 6.0_DP)/&
       & ((vi - 2.0_DP)*(s + vi - 1.0_DP))
  B(3:N) = -4.0_DP/(s + vi - 1.0_DP)
  C(3:N) = -(3.0_DP*vi**2 - 9.0_DP*vi - s*vi + 6.0_DP)/&
       & ((vi - 2.0_DP)*(s + vi - 1.0_DP))
  D(3:N) = -4.0_DP*(vi - 1.0_DP)/((vi - 2.0_DP)*(s + vi - 1.0_DP))
  E(3:N) = -(vi - 1.0_DP)*(s - vi + 2.0_DP)/((vi - 2.0_DP)*(s + vi - 1.0_DP))

  ! seed recurrence
  phi(0) = 1.0_DP
  phi(1) = 1.0_DP - 2.0_DP*x/s
  phi(2) = 1.0_DP - 8.0_DP*x/s + 8.0_DP*x**2/(s*(s + 1.0_DP))

  ! compute recurrence
  do j=3,N
     phi(j) = (A(j) + B(j)*x)*phi(j-1) + &
          & (C(j) + D(j)*x)*phi(j-2) + E(j)*phi(j-3)
  end do
    
end function hypergeom22
end module invlap_piessens
