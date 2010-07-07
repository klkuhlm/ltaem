! this module contains the routine that calls the individual
! element sub-matrix routines, and constructs the larger
! matrix used to compute the solution via LAPACK.

module solution_mod
  implicit none

  private
  public :: matrix_solution

contains
  subroutine matrix_solution(c,e,dom,sol,p,idx)
    use constants, only : DP
    use type_definitions, only : circle, solution, ellipse, domain, match_result, print_match_result
    use circular_elements, only : circle_match
    use elliptical_elements, only : ellipse_match

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
    type(match_result), allocatable :: res(:,:) ! results(target_id,source_id)
    integer, allocatable :: row(:,:), col(:,:) ! row/column trackers
    integer :: nc, ne, ntot, i, j, bigM, bigN, rr,cc

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
       print '(4(A,I0))', 'before c on self: ',i,' N:',c(i)%N,' M:',c(i)%M,' ibnd:',c(i)%ibnd
       res(i,i) = circle_match(c(i),p)
       call print_match_result(res(i,i))
       row(i,1) = size(res(i,i)%LHS,1)
       col(i,1) = size(res(i,i)%LHS,2)
       print '(2(A,I0))', 'row ',row(i,1),' col ',col(i,1)

       ! circle on other circle
       do j=1,nc
          if(i/=j) then
             print '(A,2(1X,I0),4(A,I0))', 'before c on c:',i,j,' N:',c(i)%N,' M:',c(j)%M, &
                  &' <-ibnd:',c(i)%ibnd,' ->ibnd:',c(j)%ibnd
             res(j,i) = circle_match(c(i),c(j)%matching,dom,p)
             call print_match_result(res(j,i))
             print '(A,2(1X,I0))', 'after c on c:',i,j
          end if
       end do

       ! circle on other ellipse
       do j=1,ne
          print '(A,2(1X,I0),4(A,I0))', 'before c on e:',i,j+nc,' N:',c(i)%N,' M:',e(j)%M, &
               &' <-ibnd:',c(i)%ibnd,' ->ibnd:',e(j)%ibnd
          res(j+nc,i) = circle_match(c(i),e(j)%matching,dom,p)
          call print_match_result(res(j+nc,i))
         print '(A,2(1X,I0))', 'after c on e:',i,j+nc
       end do
    end do

    do i = 1, ne
       ! ellipse on self
       print '(5(A,I0))', 'before e on self: ',nc+i,' MS:',e(i)%ms,&
            &' N:',e(i)%N,' M:',e(i)%M,' ibnd:',e(i)%ibnd
       res(nc+i,nc+i) = ellipse_match(e(i),p,idx)
       call print_match_result(res(nc+i,nc+i))
       row(i+nc,1) = size(res(nc+i,nc+i)%LHS,1)
       col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)
       print '(2(A,I0))', 'row ',row(i+nc,1),' col ',col(i+nc,1)

       ! ellipse on other circle
       do j = 1, nc
          print '(A,2(1X,I0))', 'before e on c:',nc+i,j
          res(j,nc+i) = ellipse_match(e(i),c(j)%matching,dom,p,idx)
          call print_match_result(res(j,nc+i))
          print '(A,2(1X,I0))', 'after e on c:',nc+i,j
       end do

       ! ellipse on other ellipse
       do j = 1, ne
          if (i /= j) then
             print '(A,2(1X,I0))', 'before e on e:',nc+i,nc+j
             res(nc+j,nc+i) = ellipse_match(e(i),e(j)%matching,dom,p,idx)
             call print_match_result(res(nc+j,nc+i))
             print '(A,2(1X,I0))', 'after e on e:',nc+i,nc+j
          end if
       end do
    end do

    bigM = sum(row(:,1))
    bigN = sum(col(:,1))
     allocate(A(bigM,bigN), b(bigM))
     print *, 'N,M',bigN,bigM,'::',shape(A),'::',shape(b)

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

    print '(A,I0,1X,I0)','shape(res): ',shape(res)

    ! convert structures into single matrix for solution via least squares
    do rr=1,ntot
       do cc=1,ntot
          print '(2(A,I0))', 'row ',rr,' col ',cc
          print '(2(A,2(1X,I0)))','row lo:hi',row(rr,0),row(rr,2),'  col lo:hi',col(cc,0),col(cc,2)
          print '(A,2(1X,I0))', 'LHS shape:',shape(res(rr,cc)%LHS)
          A(row(rr,0):row(rr,2),col(cc,0):col(cc,2)) = res(rr,cc)%LHS
          print '(A,2(1X,I0))', 'RHS shape:',shape(res(rr,cc)%RHS)
          b(row(rr,0):row(rr,2)) = res(rr,cc)%RHS
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
          print *, 'allocate c',i
          allocate(c(i)%coeff(sol%totalnP,col(i,1)))
       end if
       c(i)%coeff(idx,:) = b(col(i,0):col(i,2))
    end do
    do i=1,ne
       ! ellipses
       if (.not. allocated(e(i)%coeff)) then
          print *, 'allocate e',i
          allocate(e(i)%coeff(sol%totalnp,col(nc+i,1)))
       end if
       e(i)%coeff(idx,:) = b(col(nc+i,0):col(nc+i,2))
    end do
    deallocate(A,b,row,col,stat=ierr)
    if (ierr /= 0) stop 'solution.f90: error deallocating A,B,row,col'

  end subroutine matrix_solution
end module solution_mod
