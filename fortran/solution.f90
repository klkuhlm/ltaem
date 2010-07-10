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
    use circular_elements, only : circle_match, well
    use elliptical_elements, only : ellipse_match ! line() is repeated below to accommodate gfortran bug

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
       row(i,1) = size(res(i,i)%RHS,1)
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
       row(i+nc,1) = size(res(nc+i,nc+i)%RHS,1)
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

    bigM = sum(row(:,1)) ! total number rows/cols
    bigN = sum(col(:,1))
    allocate(A(bigM,bigN), b(bigM))
    b = 0.0
    print *, 'N,M',bigN,bigM,'::',shape(A),'::',shape(b)

    forall (i=1:ntot)
       row(i,0) = 1 + sum(row(1:i-1,1))  ! lower bound
       row(i,2) = sum(row(1:i,1))        ! upper bound
       col(i,0) = 1 + sum(col(1:i-1,1))
       col(i,2) = sum(col(1:i,1))
    end forall

    print '(A,I0,1X,I0)','shape(res): ',shape(res)

    ! convert structures into single matrix for solution via least squares
    do rr=1,ntot
       do cc=1,ntot
          print '(2(A,I0))', 'row ',rr,' col ',cc
          print '(2(A,2(1X,I0)))','row lo:hi',row(rr,0),row(rr,2),'  col lo:hi',col(cc,0),col(cc,2)
          print '(A,2(1X,I0))', 'LHS shape:',shape(res(rr,cc)%LHS)
          A(row(rr,0):row(rr,2),col(cc,0):col(cc,2)) = res(rr,cc)%LHS
          print '(A,2(1X,I0))', 'RHS shape:',shape(res(rr,cc)%RHS)
          b(row(rr,0):row(rr,2)) = b(row(rr,0):row(rr,2)) + res(rr,cc)%RHS
       end do
    end do
    deallocate(res,stat=ierr)
    if (ierr /= 0) stop 'solution.f90: error deallocating res'

    ! use LAPACK routine to solve least-squares via Q-R decomposition
    call ZGELS('N',bigM,bigN,1,A(:,:),bigM,b(:),bigM,WORK,IWORK,ierr)
    if (ierr /= 0) then
       write(*,'(A,I0,2(A,ES10.3))') 'ZGELS error: ',ierr,' p:',real(p),'+i',aimag(p)
       stop 
    end if
       
    ! put result into local coeff variables
    do i=1,nc
       ! circles
       if (.not. allocated(c(i)%coeff)) then
          allocate(c(i)%coeff(sol%totalnP,col(i,1)), stat=ierr)
          if (ierr /= 0) then
             print '(A,I0,A)', 'solution.f90: error allocating c(',i,')%coeff'
             stop
          end if
       end if
       if (.not. (c(i)%ibnd == 2 .and. (.not. c(i)%storin))) then
          ! coefficients come from least-squares solution above
          c(i)%coeff(idx,:) = b(col(i,0):col(i,2))
       else
          ! a specified-flux point source (known strength, and zero unknowns)
          if (size(c(i)%coeff,dim=2) == 0) then
             ! fix size of coefficient container
             deallocate(c(i)%coeff, stat=ierr)
             if (ierr /= 0) stop 'solution.f90: error deallocating c(i)%coeff'
             allocate(c(i)%coeff(sol%totalnP,1), stat=ierr)
             if (ierr /= 0) stop 'solution.f90: error re-allocating c(i)%coeff'
          end if
          ! get coefficients from well routine
          c(i)%coeff(idx,1) = well(c(i),p)
       end if
    end do
    do i=1,ne
       ! ellipses
       if (.not. allocated(e(i)%coeff)) then
          allocate(e(i)%coeff(sol%totalnp,col(nc+i,1)), stat=ierr)
          if (ierr /= 0) then
              print '(A,I0,A)', 'solution.f90: error allocating e(',i,')%coeff'
              stop
           end if
       end if
       if (.not. e(i)%ibnd == 2) then
          ! coefficients from least-squares solution above
          e(i)%coeff(idx,:) = b(col(nc+i,0):col(nc+i,2))
       else
          if (size(e(i)%coeff,dim=2) == 0) then
             ! fix size of coefficient container
             deallocate(e(i)%coeff, stat=ierr)
             if (ierr /= 0) stop 'solution.f90: error deallocating e(i)%coeff'
             allocate(e(i)%coeff(sol%totalnP,2*e(i)%N-1), stat=ierr)
             if (ierr /= 0) stop 'solution.f90: error re-allocating e(i)%coeff'
          end if
          ! get coefficients from line routine
          e(i)%coeff(idx,:) = 0.0 
          e(i)%coeff(idx,1:e(i)%N:2) = line(e(i),p,idx) ! a_(2n)
       end if
    end do
    deallocate(A,b,row,col,stat=ierr)
    if (ierr /= 0) stop 'solution.f90: error deallocating A,B,row,col'

  end subroutine matrix_solution

  ! this function should be in elliptical_elements module, but gfortran bug forces it to be here.
  function line(e,p,idx) result(a2n)
    ! this function returns the coefficients for a specified-flux line source
    use constants, only : DP, PI
    use time_mod, only : time
    use type_definitions, only : ellipse
    use mathieu_functions, only : Ke,dKe
    type(ellipse), intent(in) :: e
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx
    complex(DP), dimension(ceiling(e%N/2.0)) :: a2n ! only even coefficients of even order
    real(DP), dimension(0:e%ms-1) :: vs
    real(DP), dimension(1:e%ms,ceiling(e%N/2.0)) :: arg
    integer, dimension(0:e%ms-1) :: vi
    integer :: i, N, MS, nmax

    N = e%N 
    MS = e%ms
    nmax = ceiling(e%N/2.0)
    vi(0:MS-1) = [(i,i=0,MS-1)]  ! integer vector
    where (mod(vi,2)==0) ! sign vector
       vs = 1.0_DP
    elsewhere
       vs = -1.0_DP
    end where
    arg(1:MS,1:nmax) = spread(vs(0:MS-1)/real(1-(2*vi(0:MS-1))**2,DP),dim=2,ncopies=nmax)
    
    ! factor of 4 different from Kuhlman&Neuman paper
    ! include Radial/dRadial MF here to balance with those in general solution
    a2n(1:nmax) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI)* &
            & Ke(e%parent%mat(idx), vi(0:N-1:2), e%r) / dKe(e%parent%mat(idx), vi(0:N-1:2), e%r)* &
            & (-vs(0:N-1:2))*sum(arg(1:MS,1:nmax)*conjg(e%mat(idx)%A(1:MS,0:nmax-1,0)),dim=1)

  end function line


end module solution_mod
