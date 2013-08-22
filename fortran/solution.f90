!
! Copyright (c) 2011,2012,2013 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

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
    use elliptical_elements, only : ellipse_match, line

    interface  !! solve over-determined system via least-squares in LAPACK
       subroutine ZGELS(TRANSA, M, N, NRHS, A, LDA, B, LDB, WORK, &
            & LDWORK, INFO)
         integer, intent(in) :: M, N, NRHS, LDA, LDB, LDWORK
         character(LEN=1), intent(in) :: TRANSA
         complex(KIND=8), intent(inout), dimension(LDWORK) :: WORK
         complex(KIND=8), intent(inout), dimension(LDA,N) :: A
         complex(KIND=8), intent(inout), dimension(LDB,NRHS) ::  B
         integer, intent(out) :: INFO
       end subroutine ZGELS
    end interface

    type(circle),  dimension(:), intent(inout) :: c
    type(ellipse), dimension(:), intent(inout) :: e
    type(domain), intent(in) :: dom
    type(solution), intent(in) :: sol
    complex(DP), intent(in) :: p
    integer, intent(in) :: idx

    complex(DP), allocatable :: A(:,:), b(:)
    type(match_result), allocatable :: res(:,:) ! results(target_id,source_id)
    integer, allocatable :: row(:,:), col(:,:) ! row/column trackers
    integer :: nc, ne, ntot, i, j, k, bigM, bigN, rr,cc, ierr

    ! only needed for LAPACK routine; size(work) = 33*bigN
    complex(DP), allocatable :: WORK(:)
    
    integer :: nrow, ncol
    character(41) :: fmt
    fmt = "(I0,1X,   ('(',ES8.1,',',ES8.1,')',1X))"

    nc = size(c,dim=1)
    ne = size(e,dim=1)
    ntot = nc + ne

    ! row/col have three entries per row or column
    ! 0 : lower bound of sub-block in A matrix
    ! 1 : size of sub-block 
    ! 2 : upper bound of sub-block in A matrix

    allocate(res(ntot,ntot), row(ntot,0:2), col(ntot,0:2))

    ! accumulate result into matrices of structures
    do i = 1,nc
       ! circle on self
       res(i,i) = circle_match(c(i),p)
       row(i,1) = size(res(i,i)%RHS,1) 
       col(i,1) = size(res(i,i)%LHS,2)

       print *, 'SOL circ-self i:',i,'LHS shape:',shape(res(i,i)%LHS),'&
            &RHS shape:',shape(res(i,i)%RHS),'row:',row(i,1),'col:',col(i,1)

       print *, 'SOL circ-self i:',i,'LHS'
       nrow = size(res(i,i)%LHS,1)
       ncol = size(res(i,i)%LHS,2)
       if (ncol > 0) then 
          write(fmt(8:10),'(I3.3)') ncol
          do j = 1, nrow
             write(*,fmt) j,res(i,i)%LHS(j,:)
          end do
       else
          write(*,'(A)') 'no output ncol==0'
       end if
       
       print *, 'SOL circ-self i:',i,'RHS'
       nrow = size(res(i,i)%RHS)
       if (nrow > 0) then
          write(fmt(8:10),'(I3.3)') nrow
          write(*,fmt) 1,res(i,i)%RHS(:)
       else
          write(*,'(A)') 'no output nrow==0'
       end if
       
       ! circle on other circle
       do j = 1,nc
          if(i /= j) then
             res(j,i) = circle_match(c(i),c(j)%matching,dom,p)
             print *, 'SOL circ-circ i,j:',i,j,'LHS shape:',shape(res(j,i)%LHS),&
                  &'RHS shape:',shape(res(j,i)%RHS)

             print *, 'SOL circ-circ i,j:',i,j,'LHS'
             nrow = size(res(j,i)%LHS,1) 
             ncol = size(res(j,i)%LHS,2) 
             if (ncol > 0) then
                write(fmt(8:10),'(I3.3)') ncol
                do k = 1, nrow
                   write(*,fmt) k,res(j,i)%LHS(k,:)
                end do
             else
                write(*,'(A)') 'no output ncol==0'
             end if
             
             print *, 'SOL circ-circ i,j:',i,j,'RHS'
             nrow = size(res(j,i)%RHS(:)) 
             if (nrow > 0) then
                write(fmt(8:10),'(I3.3)') nrow
                write(*,fmt) 1,res(j,i)%RHS(:)
             else
                write(*,'(A)') 'no output nrow==0'
             end if
             
          end if
       end do

       ! circle on other ellipse
       do j = 1,ne
          res(j+nc,i) = circle_match(c(i),e(j)%matching,dom,p)
          print *, 'SOL circ-ellip i,j:',i,j,'LHS shape:',shape(res(j+nc,i)%LHS),&
               &'RHS shape:',shape(res(j+nc,i)%RHS)
       end do
    end do

    do i = 1,ne
       ! ellipse on self
       res(nc+i,nc+i) = ellipse_match(e(i),p,idx)
       row(i+nc,1) = size(res(nc+i,nc+i)%RHS,1)
       col(i+nc,1) = size(res(nc+i,nc+i)%LHS,2)

       print *, 'SOL ellip-self i:',i,'LHS shape:',shape(res(nc+i,nc+i)%LHS),&
            &'RHS shape:',shape(res(nc+i,nc+i)%RHS),&
            & 'row:',row(i+nc,1),'col:',col(i+nc,1)

       ! ellipse on other circle
       do j = 1,nc
          res(j,nc+i) = ellipse_match(e(i),c(j)%matching,dom,p,idx)
          print *, 'SOL ellip-circ i,j:',i,j,'LHS shape:',shape(res(j,nc+i)%LHS),&
               &'RHS shape:',shape(res(j,nc+i)%RHS)
       end do

       ! ellipse on other ellipse
       do j = 1,ne
          if (i /= j) then
             res(nc+j,nc+i) = ellipse_match(e(i),e(j)%matching,dom,p,idx)
             print *, 'SOL ellip-ellip i,j:',i,j,'LHS shape:',shape(res(j+nc,nc+i)%LHS),&
                  &'RHS shape:',shape(res(j+nc,nc+i)%RHS)
          end if
       end do
    end do

    bigM = sum(row(:,1)) ! total number rows/cols
    bigN = sum(col(:,1))

    print *, 'SOL bigM:',bigM,' bigN:',bigN

    allocate(A(bigM,bigN), b(bigM))
    b = cmplx(0,0,DP)

    if (any(c%match) .or. any(e%match)) then
       allocate(work(33*bigN))
    end if

    forall (i = 1:ntot)
       row(i,0) = 1 + sum(row(1:i-1,1))  ! lower bound
       row(i,2) = sum(row(1:i,1))        ! upper bound
       col(i,0) = 1 + sum(col(1:i-1,1))
       col(i,2) = sum(col(1:i,1))
    end forall

    print *, 'SOL row:'
    do i=1,size(row,1)
       print *, i,':',row(i,:)
    end do
    
    print *, 'SOL col:'
    do i=1,size(col,1)
       print *, i,':',col(i,:)
    end do
    
    ! convert structures into single matrix for solution via least squares
    do rr = 1,ntot
       do cc = 1,ntot
          print *, 'convert',rr,cc,' row range:',row(rr,0),row(rr,2),&
                                 & ' col range:',col(cc,0),col(cc,2)
          A(row(rr,0):row(rr,2),col(cc,0):col(cc,2)) = res(rr,cc)%LHS
          b(row(rr,0):row(rr,2)) = b(row(rr,0):row(rr,2)) + res(rr,cc)%RHS
       end do
    end do

    deallocate(res)

    if (any(c%match) .or. any(e%match)) then
       ! use LAPACK routine to solve least-squares via Q-R decomposition
       ! this routine works for all three potential use cases
       ! M>N (overdetermined), M==N (even-determined), and M<N (underdetermined)

       print *, 'ZGELS debug',bigM,bigN,size(work),':: bshape',shape(b),':: Ashape',shape(A)
       call ZGELS(TRANSA='N', M=bigM, N=bigN, NRHS=1, A=A(:,:), LDA=bigM, B=b(:), &
            & LDB=bigM, WORK=work, LDWORK=size(work), INFO=ierr)
       if (ierr /= 0) then
          write(*,'(A,I0,2(A,ES10.3))') 'ZGELS error: ',ierr,' p:',real(p),'+i',aimag(p)
          stop 999
       else
       end if
    end if

    ! put result into local coeff variables
    do i = 1,nc
       ! Circles -- ensure container for results is allocated
       ! allocation should only happen on first pass through this subroutine
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

    do i = 1,ne
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
       end if
    end do
    deallocate(A,b,row,col)

  end subroutine matrix_solution

end module solution_mod
