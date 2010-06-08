! this module contains the routine that calls the individual
! element sub-matrix routines, and constructs the larger
! matrix used to compute the solution via LAPACK.

module solution

contains
subroutine matrix_solution(c,e,dom,p,idx)
  use constants, only : DP
  use type_definitions, only : circle, ellipse, match_result
  use circular_elements
  use elliptical_elements
  
  type(circle), dimension(:), intent(inout) :: c
  type(ellipse), dimension(:), intent(inout) :: e
  type(domain), intent(in) :: dom
  complex(DP), intent(in) :: p
  
  complex(DP), allocatable :: A(:,:)
  type(match_result), allocatable(:,:,:) :: res 
  integer, allocatable :: row(:,:), col(:,:)
  integer :: nc, ne, ntot, i, j
  
  nc = size(c,1)
  ne = size(e,1)
  ntot = nc + ne
  
  allocate(res(ntot,ntot,2),row(ntot,0:2),col(ntot,0:2))

  ! accumulate results into matrices of structures
  do i = 1, nc
     ! circle on self
     res(i,i,1) = circle_head(c(i),p)
     res(i,i,2) = circle_flux(c(i),p)
     
     col(i,1) = size(res(i,i,1)%LHS,2)
     
     ! circle on other circle
     do j = 1, nc
        if (i /= j) then
           res(i,j,1) = circle_head(c(i),c(j)%matching,dom,p)
           res(i,j,2) = circle_flux(c(i),c(j)%matching,dom,p)
        end if
     end do

     ! circle on other ellipse
     do j = 1, ne
        res(i,j+nc,1) = circle_head(c(i),e(j)%matching,dom,p)
        res(i,j+nc,2) = circle_flux(c(i),e(j)%matching,dom,p)
     end do
  end do

  where (c(1:nc)%ibnd == 0)
     row(1:nc,1) = 2*c(1:nc)%M
  elsewhere
     row(1:nc,1) = c(1:nc)%M
  end where

  do i = 1, ne
     ! ellipse on self
     res(nc+i,nc+i,1) = ellipse_head(e(i),p)
     res(nc+i,nc+i,2) = ellipse_flux(e(i),p)
     
     col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)

     ! ellipse on other circle
     do j = 1, nc
        res(nc+i,j,1) = circle_head(e(i),c(j)%matching,dom,p)
        res(nc+i,j,2) = circle_flux(e(i),c(j)%matching,dom,p)
     end do

     ! ellipse on other ellipse
     do j = 1, ne
        if (i /= j) then
           res(nc+i,j+nc,1) = circle_head(e(i),e(j)%matching,dom,p)
           res(nc+i,j+nc,2) = circle_flux(e(i),e(j)%matching,dom,p)
        end if
     end do
  end do
  
  where (e(1:ne)%ibnd == 0)
     row(nc+1:ntot,1) = 2*e(1:ne)%M
  elsewhere
     row(nc+1:ntot,1) = e(1:ne)%M
  end where

  allocate(A(sum(row),sum(col)))

  forall (i=1:ntot)
     row(i,0) = 1 + sum(row(1:i-1))  ! lower bound
     row(i,2) = sum(row(1:i))        ! upper bound
     col(i,0) = 1 + sum(col(1:i-1))
     col(i,2) = sum(col(1:i))
  end forall

  ! convert structures into single matrix for solution via least squares
  do i=1,ntot
     do j=1,ntot
        
        A(row(i,0),row(j,2)) = res(i,j)
     end do
  end do
  
end subroutine matrix_solution

end module solution
