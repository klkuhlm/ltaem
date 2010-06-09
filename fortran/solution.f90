! this module contains the routine that calls the individual
! element sub-matrix routines, and constructs the larger
! matrix used to compute the solution via LAPACK.

module solution

contains
  subroutine matrix_solution(c,e,dom,p,idx)
    use constants, only : DP
    use type_definitions, only : circle, ellipse, domain, match_result
    use circular_elements
    use elliptical_elements

    !! solve over-determined system via least-squares
    interface
       subroutine ZGELS(TRANSA, M, N, NRHS, A, LDA, B, LDB, WORK, &
            & LDWORK, INFO)
         integer, intent(in) :: M, N, NRHS, LDA, LDB, LDWORK
         character(LEN=1), intent(in) :: TRANSA
         complex(8), intent(inout), dimension(LDWORK) :: WORK
         complex(8), intent(inout), dimension(LDA,N) :: A
         complex(8), intent(inout), dimension(LDB,NRHS) ::  B
         integer, intent(out) :: INFO
       end subroutine ZGELS
    end interface

    type(circle), dimension(:), intent(inout) :: c
    type(ellipse), dimension(:), intent(inout) :: e
    type(domain), intent(in) :: dom
    complex(DP), intent(in) :: p

    complex(DP), allocatable :: A(:,:), b(:)
    type(match_result), allocatable :: res(:,:) ! results
    integer, allocatable :: row(:,:), col(:,:) ! row/column trackers
    integer :: nc, ne, ntot, i, j, bigM, bigN

    ! things only needed for LAPACK routine
    integer, parameter :: iwork = 10250   ! <- optimal??? probably not
    complex(DP), dimension(iwork) :: WORK
    integer :: IERR

    nc = size(c,1);  ne = size(e,1)
    ntot = nc + ne
    allocate(res(ntot,ntot),row(ntot,0:2),col(ntot,0:2))

    ! accumulate results into matrices of structures
    do i=1,nc
       ! circle on self
       res(i,i) = circle_match(c(i),p)
       row(i,1) = size(res(i,i)%LHS,1)
       col(i,1) = size(res(i,i)%LHS,2)

       ! circle on other circle
       do j=1,nc
          if(i/=j) then
             res(i,j) = circle_match(c(i),c(j)%matching,dom,p)
          end if
       end do

       ! circle on other ellipse
       do j=1,ne
          res(i,j+nc) = circle_match(c(i),e(j)%matching,dom,p)
       end do
    end do

    do i = 1, ne
       ! ellipse on self
       res(nc+i,nc+i) = ellipse_match(e(i),p,idx)
       row(i+nc,1) = size(res(nc+i,nc+i)%LHS,1)
       col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)

       ! ellipse on other circle
       do j = 1, nc
          res(nc+i,j) = ellipse_match(e(i),c(j)%matching,dom,p,idx)
       end do

       ! ellipse on other ellipse
       do j = 1, ne
          if (i /= j) then
             res(nc+i,j+nc) = ellipse_match(e(i),e(j)%matching,dom,p,idx)
          end if
       end do
    end do

    bigM = sum(row(:,1))
    bigN = sum(col(:,1))
    allocate(A(bigM,bigN), b(bigM))

    forall (i=1:ntot)
       row(i,0) = 1 + sum(row(1:i-1,1))  ! lower bound
       row(i,2) = sum(row(1:i,1))        ! upper bound
       col(i,0) = 1 + sum(col(1:i-1,1))
       col(i,2) = sum(col(1:i,1))
    end forall

    ! convert structures into single matrix for solution via least squares
    do i=1,ntot
       do j=1,ntot
          A(row(i,0):row(i,2),col(j,0):col(j,2)) = res(i,j)%LHS
          b(row(i,0):row(i,2)) = res(i,j)%RHS
       end do
    end do
    deallocate(res)

    ! use LAPACK routine to solve least-squares via Q-R decomposition
    call ZGELS('N',bigM,bigN,1,A(:,:),bigM,b(:),bigM,WORK,IWORK,ierr)
    if (ierr /= 0) then
       write(*,'(A,I0,A,ES10.3,A,ES10.3)') 'ZGELS error: ',ierr,' p:',real(p),'+i',aimag(p)
       stop 
    end if
       
    ! put results into local coeff variables
    do i=1,nc
       ! circles
       allocate( c(i)%coeff(1:col(i,1)) )
       c(i)%coeff(:) = b(col(i,0):col(i,2))
    end do
    do i=1,ne
       ! ellipses
       allocate(e(i)%coeff(1:col(nc+i,1)))
       e(i)%coeff(:) = b(col(nc+i,0):col(nc+i,2))
    end do
    deallocate(A,b,row,col)

  end subroutine matrix_solution
end module solution
