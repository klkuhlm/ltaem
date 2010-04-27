! $Id: shared_mathieu.f90,v 1.1 2007/07/27 16:27:26 kris Exp kris $
module shared_mathieu
  use constants, only : DP
  implicit none

  public :: A,B,q

  complex(DP) :: q
  complex(DP), allocatable :: mcn(:)
  complex(DP), allocatable :: A(:,:,:), B(:,:,:)

end module shared_mathieu
