module mcn_matrix_method
  implicit none

  ! $Id: mcn_eigenvalues.f90,v 1.8 2008/03/28 18:12:39 kris Exp kris $
  ! written by K Kuhlman, July 2006

  private
  public :: mcn_eigenvalues

  contains
  subroutine mcn_eigenvalues(q,AA,BB,norm)
    use constants, only : DP, CZERO, RONE

!! things in pre-processor macros required for ode debugging program
#ifdef ODECHECK
    use shared_mathieu, only : mcn
#endif

    interface ! LAPACK 3.3.1 eigenvalue/eigenfunction routine
       subroutine ZGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,&
            & LWORK,RWORK,INFO)
         integer, intent(in) :: LDVR, LDVL, LDA, N, LWORK
         character(LEN=1), intent(in) :: JOBVL, JOBVR
         complex(KIND=8), intent(inout) :: A(N,N), WORK(LWORK)
         complex(KIND=8), intent(out) :: W(N), VL(LDVL,N), VR(LDVR,N)
         real(KIND=8), intent(inout) :: RWORK(2*N)
         integer, intent(out) :: INFO
       end subroutine ZGEEV
    end interface

#ifdef DEBUG
    character(30) :: fmt
#endif

    complex(DP), intent(in) :: q

      ! normalization as defined by 
      ! McLachlan (1), Stratton (2), or none (other)
    integer, intent(in) :: norm

      ! Coeff 1st dimension: subscript in McLachlan notation,
      ! i.e., the position in the infinite sum
    
      ! Coeff 2nd dimension: superscript in McLachlan notation,
      ! i.e., the order of the Mathieu function which it is associated with

      ! Coeff 3rd dimension: 1(even) or 2(odd) cases of the second dimension
      ! allowing matrix to be stored naturally without large number of zeros

    complex(DP), intent(inout) :: AA(1:,0:,1:), BB(1:,0:,1:)

#ifndef ODECHECK
      ! the mcns are not really utilized, so it is sort of a dummy variable
    complex(DP), dimension(size(AA,dim=1)*4) :: mcn
#endif

    integer :: M, di,i,j
    complex(DP), dimension(1) :: dc ! dummy arg

    ! just one matrix of recursion coefficients (used 4 times)
    complex(DP), dimension(size(AA,dim=1),size(AA,dim=1)) :: Coeff
    
    ! parameters for lapack routine (optimal?)
    !! lwork=410 for celeron laptop
    !! lwork=570 for PentiumD 64-bit desktop
    integer, parameter :: lwork = 410
    complex(DP), dimension(lwork) :: work
    real(DP), dimension(2*size(AA,dim=1)) :: rwork
    integer :: info

    !! used in stratton normalization only
    real(DP), dimension(size(AA,dim=1)) :: vc2, vc3

#ifdef ODECHECK    
    if(.not. allocated(mcn)) allocate(mcn(size(AA,dim=1)*4))
    mcn = cmplx(0.0,0.0,DP)
#endif

    ! check to provide helpful message when I forget to allocate arrays
    if(size(AA,dim=1) < 2 .or. size(BB,dim=1) < 2) then
       print *, 'A and B shared arrays must be allocated&
            & before calling mcn_eigenvalues'
       stop
    end if
    
    M = size(AA,dim=1)
    if(norm == 2) then
       vc2 = (/( real(2*i+1,DP), i=0,m-1)/)
       vc3 = (/( real(2*i+2,DP), i=0,m-1)/)
    end if
    
    di = 1 ! dummy integer
    dc(1) = CZERO ! dummy complex

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of even order
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO

    ! main diagonal/ r counting from 0:m-1 like McLachlan
    forall(i=1:M-1) Coeff(i+1,i+1) = cmplx((2*i)**2, 0, DP)
    
    ! off diagonals
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) Coeff(i,j) = q
    
    ! special cases
    Coeff(2,1) = sqrt(2.0_DP)*q
    coeff(1,2) = Coeff(2,1)

#ifdef DEBUG
    open(unit=99,file='recurrence_a_ev.dat')
    fmt = '(   (2(ES9.2,1X)))'
    write(fmt(2:4),'(I3.3)') M
    write(99,*) '# recurrence matrix for even a, q=',q
    do i=1,M
       write(99,fmt) coeff(i,:)
    end do
    close(99)
#endif
    
    call zgeev('N','V',M,Coeff,M,mcn(1:m),dc,di,AA(1:m,0:m-1,1),M,&
         & work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of even order"

    if (norm == 1) then 
       ! McLachlan, defined in terms of norm which has complex counterpart
       ! use 'dot product' to compute length^2, divide by length
       
       ! ZGEEV returns eigenvectors of unit euclidian length already
       ! only the first eigenvector must be re-scaled, then

!!$       AA(1:m,1:m-1,1) = AA(1:m,1:m-1,1)/&
!!$            spread(sqrt(sum(AA(1:m,1:m-1,1)*conjg(AA(1:m,1:m-1,1)),dim=1)),&
!!$            & dim=1,ncopies=m)
       !! handle ce_0 separately
       AA(1:m,0,1) = AA(1:m,0,1)/sqrt(2.0_DP*AA(1,0,1)*conjg(AA(1,0,1)) + &
            & sum(AA(2:m,0,1)*conjg(AA(2:m,0,1))))

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(AA(i,i-1,1)) < 0.0) then
             AA(:,i-1,1) = -AA(:,i-1,1)
          end if
       end do

    elseif (norm == 2) then 
       ! Stratton (sum of coefficients is unity)
       !! does this make sense for __complex__ coefficients??
       AA(1:m,0:m-1,1) = AA(:,:,1)/ &
            & spread(sum(AA(:,:,1),dim=1),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of odd order
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO
    Coeff(1,1) = RONE + q
    
    forall(i=1:m-1) Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) Coeff(i,j) = q    

#ifdef DEBUG
    open(unit=99,file='recurrence_a_od.dat')
    write(99,*) '# recurrence matrix for odd a, q=',q
    do i=1,M
       write(99,fmt) coeff(i,:)
    end do
    close(99)
#endif

    call zgeev('N','V',M,Coeff,M,mcn(m+1:2*m),dc,di,AA(1:m,0:m-1,2),M,&
         &work,lwork,rwork,info)
         
    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of odd order"

    if (norm == 1) then ! McLachlan 
       !! ZGEEV returns eigenvectors with unit euclidian norm
!!$       AA(1:m,0:m-1,2) = AA(1:m,0:m-1,2)/&
!!$            spread(sqrt(sum(AA(1:m,0:m-1,2)*conjg(AA(1:m,0:m-1,2)),dim=1)),&
!!$            & dim=1,ncopies=m)

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(AA(i,i-1,2)) < 0.0) then
             AA(:,i-1,2) = -AA(:,i-1,2)
          end if
       end do
    elseif (norm == 2) then ! Stratton
       AA(1:m,0:m-1,2) = AA(:,:,2)/ &
            & spread(sum(AA(:,:,2),dim=1),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of even order
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO 
    
    ! this one not shifted by one, since 2n+2 -> 2n 
    !! (but starting from one, rather than zero)
    forall(i=1:m) Coeff(i,i) = cmplx((2*i)**2, 0, DP)
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) Coeff(i,j) = q 

#ifdef DEBUG
    open(unit=99,file='recurrence_b_ev.dat')
    write(99,*) '# recurrence matrix for even b, q=',q
    do i=1,M
       write(99,fmt) coeff(i,:)
    end do
    close(99)
#endif
   
    call zgeev('N','V',M,Coeff,M,mcn(2*m+1:3*m),dc,di,BB(1:m,0:m-1,1),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of even order"

    if (norm == 1) then ! McLachlan
!!$       BB(1:m,0:m-1,1) = BB(1:m,0:m-1,1)/&
!!$            spread(sqrt(sum(BB(1:m,0:m-1,1)*conjg(BB(1:m,0:m-1,1)),dim=1)),&
!!$            & dim=1,ncopies=m)

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(BB(i,i-1,1)) < 0.0) then
             BB(:,i-1,1) = -BB(:,i-1,1)
          end if
       end do
    elseif (norm == 2) then ! Stratton
       BB(1:m,0:m-1,1) = BB(:,:,1)/ &
            & spread(sum(spread(vc3,dim=2,ncopies=m)*&
            & BB(:,:,1),dim=1),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of odd order
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO
    Coeff(1,1) = RONE - q
    
    forall(i=1:m-1) Coeff(i+1,i+1) = cmplx((2*i+1)**2, 0, DP)
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) Coeff(i,j) = q 

#ifdef DEBUG
    open(unit=99,file='recurrence_b_od.dat')
    write(99,*) '# recurrence matrix for odd b, q=',q
    do i=1,M
       write(99,fmt) coeff(i,:)
    end do
    close(99)
#endif

    call zgeev('N','V',M,Coeff,M,mcn(3*m+1:4*m),dc,di,BB(1:m,0:m-1,2),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of odd order"

    if (norm == 1) then ! McLachlan
!!$       BB(1:m,0:m-1,2) = BB(1:m,0:m-1,2)/&
!!$            spread(sqrt(sum(BB(1:m,0:m-1,2)*conjg(BB(1:m,0:m-1,2)),dim=1)),&
!!$            & dim=1,ncopies=m)

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(BB(i,i-1,2)) < 0.0) then
             BB(:,i-1,2) = -BB(:,i-1,2)
          end if
       end do

    elseif (norm == 2) then ! Stratton
       BB(1:m,0:m-1,2) = BB(:,:,2)/ &
            & spread(sum(spread(vc2,dim=2,ncopies=m)*&
            & BB(:,:,2),dim=1),dim=1,ncopies=m)
    end if

  end subroutine mcn_eigenvalues
  
end module mcn_matrix_method

