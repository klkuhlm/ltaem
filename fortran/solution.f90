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
    use type_definitions, only : circle, solution, ellipse, domain, match_result
    use circular_elements, only : circle_match, well
    use elliptical_elements, only : ellipse_match !!, line is repeated below to accommodate gfortran bug
#ifdef DEBUG
    use type_definitions, only : print_match_result
#endif

    interface  !! solve over-determined system via least-squares in LAPACK
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

#ifdef DEBUG
    character(13) :: fmt
#endif

    type(circle),  dimension(:), intent(inout) :: c
    type(ellipse), dimension(:), intent(inout) :: e
    type(domain), intent(in) :: dom
    type(solution), intent(in) :: sol
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx

    complex(DP), allocatable :: A(:,:), b(:)
    type(match_result), allocatable :: res(:,:) ! results(target_id,source_id)
    integer, allocatable :: row(:,:), col(:,:) ! row/column trackers
    integer :: nc, ne, ntot, i, j, bigM, bigN, rr,cc, ierr

    ! only needed for LAPACK routine
    ! size(work) should be ~ 33xbigN? (32- & 64-bit linux)
    complex(DP), allocatable :: WORK(:)

#ifdef DEBUG
    print *, 'matrix_solution: c:',c%id,' e:',e%id,' #tot el:',sum(dom%num(1:2)),' p:',p,' idx:',idx
#endif

    nc = size(c,dim=1)
    ne = size(e,dim=1)
    ntot = nc + ne
    allocate(res(ntot,ntot), row(ntot,0:2), col(ntot,0:2))

    ! accumulate result into matrices of structures
    do i=1,nc
       ! circle on self
       res(i,i) = circle_match(c(i),p)
       row(i,1) = size(res(i,i)%RHS,1)
       col(i,1) = size(res(i,i)%LHS,2)

#ifdef DEBUG
       print *, 'circle on self',i,' row',row(i,1),' col',col(i,1)
!!$       call print_match_result(res(i,i))
#endif

       ! circle on other circle
       do j=1,nc
          if(i/=j) then
             res(j,i) = circle_match(c(i),c(j)%matching,dom,p)
#ifdef DEBUG
             print *, 'circle',i,'on other circle',j
!!$             call print_match_result(res(j,i))
#endif
          end if
       end do

       ! circle on other ellipse
       do j=1,ne
          res(j+nc,i) = circle_match(c(i),e(j)%matching,dom,p)
#ifdef DEBUG
          print *, 'circle',i,'on other ellipse',j+nc
!!$          call print_match_result(res(j+nc,i))
#endif
       end do
    end do

    do i = 1, ne
       ! ellipse on self
       res(nc+i,nc+i) = ellipse_match(e(i),p,idx)
       row(i+nc,1) = size(res(nc+i,nc+i)%RHS,1)
       col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)

#ifdef DEBUG
       print *, 'ellipse on self',nc+i,' row',row(i+nc,1),' col',col(i+nc,1)
!!$       call print_match_result(res(nc+i,nc+i))
#endif

       ! ellipse on other circle
       do j = 1, nc
          res(j,nc+i) = ellipse_match(e(i),c(j)%matching,dom,p,idx)
#ifdef DEBUG
          print *, 'ellipse',i+nc,'on other circle',j
!!$          call print_match_result(res(j,nc+i))
#endif
       end do

       ! ellipse on other ellipse
       do j = 1, ne
          if (i /= j) then
             res(nc+j,nc+i) = ellipse_match(e(i),e(j)%matching,dom,p,idx)
#ifdef DEBUG
             print *, 'ellipse',i+nc,'on other ellipse',j+nc
!!$             call print_match_result(res(nc+j,nc+i))
#endif
          end if
       end do
    end do

    bigM = sum(row(:,1)) ! total number rows/cols
    bigN = sum(col(:,1))

    allocate(A(bigM,bigN), b(bigM))
    b = 0.0

    if (any(c%match) .or. any(e%match)) then
       allocate(work(33*bigN))       
    end if

    !$OMP WORKSHARE
    forall (i=1:ntot)
       row(i,0) = 1 + sum(row(1:i-1,1))  ! lower bound
       row(i,2) = sum(row(1:i,1))        ! upper bound
       col(i,0) = 1 + sum(col(1:i-1,1))
       col(i,2) = sum(col(1:i,1))
    end forall
    !$OMP END WORKSHARE

#ifdef DEBUG
    print '(A,I0,1X,I0)','shape(res): ',shape(res)
#endif

    ! convert structures into single matrix for solution via least squares
    do rr=1,ntot
       do cc=1,ntot
#ifdef DEBUG
          print '(2(A,I0))', 'row ',rr,' col ',cc
          print '(2(A,2(1X,I0)))','row lo:hi',row(rr,0),row(rr,2),'  col lo:hi',col(cc,0),col(cc,2)
          print '(A,2(1X,I0))', 'LHS shape:',shape(res(rr,cc)%LHS)
          print '(A,2(1X,I0))', 'RHS shape:',shape(res(rr,cc)%RHS)
#endif
          
          A(row(rr,0):row(rr,2),col(cc,0):col(cc,2)) = res(rr,cc)%LHS
          b(row(rr,0):row(rr,2)) = b(row(rr,0):row(rr,2)) + res(rr,cc)%RHS
          
       end do
    end do

    deallocate(res)

    if (any(c%match) .or. any(e%match)) then
       ! use LAPACK routine to solve least-squares via Q-R decomposition
       ! this routine works for all three potential use cases
       ! M>N (overdetermined), M==N (even-determined), and M<N (underdetermined)
#ifdef DEBUG
       print *, 'bigM',bigM,' bigN',bigN,' shape(A)',shape(A),' shape(b)',shape(b)
#endif
       call ZGELS(TRANSA='N',M=bigM,N=bigN,NRHS=1,A=A(:,:),LDA=bigM,B=b(:),LDB=bigM,&
            & WORK=work,LDWORK=size(work,dim=1),INFO=ierr)
       if (ierr /= 0) then
          write(*,'(A,I0,2(A,ES10.3))') 'ZGELS error: ',ierr,' p:',real(p),'+i',aimag(p)
          stop 
       else
#ifdef DEBUG
          print *, 'ZEGLS successful'
#endif
       end if
    end if

    ! put result into local coeff variables
    do i=1,nc
#ifdef DEBUG       
       print *, 'circle',i
#endif
       ! Circles -- ensure container for results is allocated
       if (.not. allocated(c(i)%coeff)) then
          ! solution for each value of p, saved as a 2D matrix
          allocate(c(i)%coeff(sol%totalnP,col(i,1)))
       end if

       if (.not. (c(i)%ibnd == 2 .and. (.not. c(i)%storin))) then
          ! coefficients come from least-squares solution above
          c(i)%coeff(idx,:) = b(col(i,0):col(i,2))
       else
          ! a specified-flux point source (known strength, and zero unknowns)
          if (size(c(i)%coeff,dim=2) == 0) then
             ! fix size of coefficient container
             deallocate(c(i)%coeff)
             allocate(c(i)%coeff(sol%totalnP,1))
          end if
          ! get a0 coefficient from well routine
          c(i)%coeff(idx,1) = well(c(i),p)
       end if
    end do

    do i=1,ne
       ! ellipses
       if (.not. allocated(e(i)%coeff)) then
          allocate(e(i)%coeff(sol%totalnp,col(nc+i,1)))
       end if
       if (.not. e(i)%ibnd == 2) then
          ! coefficients from least-squares solution above
          e(i)%coeff(idx,:) = b(col(nc+i,0):col(nc+i,2))
       else
          if (size(e(i)%coeff,dim=2) == 0) then
             ! fix size of coefficient container
             deallocate(e(i)%coeff)
             ! allocate space for all the even (a_n) coefficients
             allocate(e(i)%coeff(sol%totalnP,2*e(i)%N-1))
          end if
          ! get coefficients from line routine (only even-order, even coeff used)
          e(i)%coeff(idx,:) = 0.0 
          e(i)%coeff(idx,1:e(i)%N:2) = line(e(i),p,idx) ! a_(2n)
#ifdef DEBUG
          print *, 'line source coefficients: N:',e(i)%N,' shape(coeff)',&
               & shape(e(i)%coeff),' shape(line)',shape(line(e(i),p,idx))
          fmt = '(  (I12,12X))'
          write(fmt(2:3),'(I2.2)') size(e(i)%coeff,dim=2)
          write(*,fmt) (j,j=1,size(e(i)%coeff,dim=2))
          do j=1,size(e(i)%coeff,dim=2)
             write(*,'(2(A,ES10.2E3),A)',advance='no'),'(',real(e(i)%coeff(idx,j)),',',aimag(e(i)%coeff(idx,j)),') '
          end do
          write(*,*)
#endif
       end if
    end do
    deallocate(A,b,row,col)

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
    vi(0:MS-1) = [(i,i=0,MS-1)]  ! integer index vector
    vs = -1.0 ! sign vector
    where (mod(vi,2)==0) vs = 1.0

    arg(1:MS,1:nmax) = spread(vs(0:MS-1)/real(1-(2*vi(0:MS-1))**2,DP),2,nmax)

    ! factor of 4 different from Kuhlman&Neuman paper
    ! include Radial/dRadial MF here to balance with those in general solution
    a2n(1:nmax) = time(p,e%time,.false.)*e%bdryQ/(2.0*PI)* &
            & Ke(e%parent%mat(idx), vi(0:N-1:2), e%r) / &
            & dKe(e%parent%mat(idx), vi(0:N-1:2), e%r)* &
            & (-vs(0:N-1:2))*sum( arg(1:MS,1:nmax)* &
            & conjg(e%parent%mat(idx)%A(1:MS,0:nmax-1,0)), dim=1)

  end function line
end module solution_mod
