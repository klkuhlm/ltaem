! this module contains the routine that calls the individual
! element sub-matrix routines, and constructs the larger
! matrix used to compute the solution via LAPACK.

module solution_mod
  implicit none

contains
  subroutine matrix_solution(c,e,dom,sol,p,idx)
    use constants, only : DP
    use type_definitions, only : circle, solution, ellipse, domain, match_result, print_match_result
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
    type(solution), intent(in) :: sol
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx

    complex(DP), allocatable :: A(:,:), b(:)
    type(match_result), allocatable :: res(:,:) ! results
    integer, allocatable :: row(:,:), col(:,:) ! row/column trackers
    integer :: nc, ne, ntot, i, j, bigM, bigN

    ! things only needed for LAPACK routine
    integer, parameter :: iwork = 10250   ! <- optimal??? probably not
    complex(DP), dimension(iwork) :: WORK
    integer :: IERR

#ifdef DEBUG
    print *, 'matrix_solution: c:',c%id,' e:',e%id,' dom:',dom%num,' p:',p,' idx:',idx
#endif

    nc = size(c,1);  ne = size(e,1)
    ntot = nc + ne
    allocate(res(ntot,ntot), row(ntot,0:2), col(ntot,0:2))

    ! accumulate result into matrices of structures
    do i=1,nc
       ! circle on self
       write(*,'(A,I0)') 'before c on self: ',i
       res(i,i) = circle_match(c(i),p)
       call print_match_result(res(i,i))
       row(i,1) = size(res(i,i)%LHS,1)
       col(i,1) = size(res(i,i)%LHS,2)
       write(*,'(2(A,I0))') 'row ',row(i,1),' col ',col(i,1)

       ! circle on other circle
       do j=1,nc
          if(i/=j) then
             write(*,'(A,2(1X,I0))') 'before c on c:',i,j
             res(i,j) = circle_match(c(i),c(j)%matching,dom,p)
             call print_match_result(res(i,j))
             write(*,'(A,2(1X,I0))') 'after c on c:',i,j
          end if
       end do

       ! circle on other ellipse
       do j=1,ne
          write(*,'(A,2(1X,I0))') 'before c on e:',i,j+nc
          res(i,j+nc) = circle_match(c(i),e(j)%matching,dom,p)
          call print_match_result(res(i,j+nc))
          write(*,'(A,2(1X,I0))') 'after c on e:',i,j+nc
       end do
    end do

    do i = 1, ne
       ! ellipse on self
       write(*,'(A,I0)') 'before e on self: ',nc+i
       res(nc+i,nc+i) = ellipse_match(e(i),p,idx)
       call print_match_result(res(nc+i,nc+i))
       row(i+nc,1) = size(res(nc+i,nc+i)%LHS,1)
       col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)
       write(*,'(2(A,I0))') 'row ',row(i+nc,1),' col ',col(i+nc,1)

       ! ellipse on other circle
       do j = 1, nc
          write(*,'(A,2(1X,I0))') 'before e on c:',nc+i,j
          res(nc+i,j) = ellipse_match(e(i),c(j)%matching,dom,p,idx)
          call print_match_result(res(nc+i,j))
          write(*,'(A,2(1X,I0))') 'after e on c:',nc+i,j
       end do

       ! ellipse on other ellipse
       do j = 1, ne
          if (i /= j) then
             write(*,'(A,2(1X,I0))') 'before e on e:',nc+i,nc+j
             res(nc+i,j+nc) = ellipse_match(e(i),e(j)%matching,dom,p,idx)
             call print_match_result(res(nc+i,nc+j))
             write(*,'(A,2(1X,I0))') 'after e on e:',nc+i,nc+j
          end if
       end do
    end do

    bigM = sum(row(:,1))
    bigN = sum(col(:,1))
     allocate(A(bigM,bigN), b(bigM))
     print *, 'N,M',bigN,bigM,shape(A),'::',shape(b)

    forall (i=1:ntot)
       row(i,0) = 1 + sum(row(1:i-1,1))  ! lower bound
       row(i,2) = sum(row(1:i,1))        ! upper bound
       col(i,0) = 1 + sum(col(1:i-1,1))
       col(i,2) = sum(col(1:i,1))
    end forall

    do i=1,ntot
       print '(A,I0,A,3(1X,I3))', 'row(',i,',0:2)',row(i,0:2)
       print '(A,I0,A,3(1X,I3))', 'col(',i,',0:2)',col(i,0:2)
    end do

    ! convert structures into single matrix for solution via least squares
    do i=1,ntot
       do j=1,ntot
          print '(A,2(1X,I0))', 'row i,j',i,j
          print '(2(A,2(1X,I0)))','row:',row(j,0),row(j,2),' col:',col(i,0),col(i,2)
          print '(A,2(1X,I0))', 'LHS:',shape(res(i,j)%LHS)
          A(row(j,0):row(j,2),col(i,0):col(i,2)) = res(i,j)%LHS
          print '(A,2(1X,I0))', 'RHS',shape(res(i,j)%RHS)
          b(row(j,0):row(j,2)) = res(i,j)%RHS
       end do
    end do
    deallocate(res,stat=ierr)
    if (ierr /= 0) stop 'solution.f90: error deallocating res'

    ! use LAPACK routine to solve least-squares via Q-R decomposition
    call ZGELS('N',bigM,bigN,1,A(:,:),bigM,b(:),bigM,WORK,IWORK,ierr)
    if (ierr /= 0) then
       write(*,'(A,I0,A,ES10.3,A,ES10.3)') 'ZGELS error: ',ierr,' p:',real(p),'+i',aimag(p)
       stop 
    end if
       
    ! put result into local coeff variables
    do i=1,nc
       ! circles
       if (.not. allocated(c(i)%coeff)) then
          allocate(c(i)%coeff(sol%totalnP,col(i,1)))
       end if
       c(i)%coeff(idx,:) = b(col(i,0):col(i,2))
    end do
    do i=1,ne
       ! ellipses
       if (.not. allocated(e(i)%coeff)) then
          allocate(e(i)%coeff(sol%totalnp,col(nc+i,1)))
       end if
       e(i)%coeff(idx,:) = b(col(nc+i,0):col(nc+i,2))
    end do
    deallocate(A,b,row,col,stat=ierr)
    if (ierr /= 0) stop 'solution.f90: error deallocating A,B,row,col'

  end subroutine matrix_solution
end module solution_mod
