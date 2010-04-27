program bes_deriv

use complex_bessel, only : cbesk, cbesi
implicit none

integer, parameter :: DP = selected_real_kind(15,300), NUM=10
complex(DP), dimension(0:NUM) :: ky, iy, Dky, Diy
complex(DP) :: arg
integer :: nz,ierr,i

arg = (14.0_DP,5.0_DP)

call cbesk(arg,0.0_DP,1,NUM+1,ky(0:NUM),nz,ierr)
call cbesi(arg,0.0_DP,1,NUM+1,iy(0:NUM),nz,ierr)
Dky(0:NUM) = besk_deriv(ky(0:NUM),arg,(1.0_DP,0.0_DP))
Diy(0:NUM) = besi_deriv(iy(0:NUM),arg,(1.0_DP,0.0_DP))


do i=0,NUM
   write(*,'(i2,8(1x,ES12.4E2))') i,real(ky(i)),aimag(ky(i)),real(iy(i)),aimag(iy(i)),&
        & real(Dky(i)),aimag(Dky(i)),real(Diy(i)),aimag(Diy(i))
end do


contains

  !! takes a vector of K bessel functions, (order 0 to nk)
  !! computes a vector of K' using recurrence
  function besk_deriv(bk,arg,Darg) result(Dbk)
    use constants, only : DP
    complex(DP), intent(in), dimension(0:) :: bk
    complex(DP), intent(in) :: arg, Darg
    complex(DP), dimension(0:size(bk)-1) :: Dbk

    integer :: nk

    nk = size(bk)-1

    Dbk(0) = -bk(1)*Darg
    Dbk(1:nk-1) = -0.5_DP*(bk(0:nk-2) + bk(2:nk))*Darg
    Dbk(nk) = -(bk(nk-1) + real(nk,DP)/arg*bk(nk))*Darg
    
  end function besk_deriv
 
   function besi_deriv(bi,arg,Darg) result(Dbi)
    use constants, only : DP
    complex(DP), intent(in), dimension(0:) :: bi
    complex(DP), intent(in) :: arg, Darg
    complex(DP), dimension(0:size(bi)-1) :: Dbi

    integer :: ni

    ni = size(bi)-1

    Dbi(0) = bi(1)*Darg
    Dbi(1:ni-1) = 0.5_DP*(bi(0:ni-2) + bi(2:ni))*Darg
    Dbi(ni) = (bi(ni-1) - real(ni,DP)/arg*bi(ni))*Darg
    
  end function besi_deriv
end program bes_deriv
