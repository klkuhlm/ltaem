module mcn_matrix_method
  implicit none

  ! $Id: mcn_eigenvalues2.f90,v 1.6 2010/02/19 23:17:40 klkuhlm Exp klkuhlm $
  ! written by K Kuhlman, July 2006
  ! updated to fix Stratton norm Feb 2010

  private
  public :: mcn_eigenvalues

  contains
  subroutine mcn_eigenvalues(q,AA,BB,norm)
    use constants, only : DP, CZERO, RONE
    use shared_mathieu, only : mcn

    interface ! LAPACK 3.2.1 eigenvalue/eigenfunction routine
       subroutine ZGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK, &
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

    complex(DP), intent(inout) :: AA(1:,0:,0:), BB(1:,0:,0:)

    integer :: M, di,i,j
    complex(DP), dimension(1) :: dc ! dummy arg

    ! just one matrix of recursion coefficients (used 4 times)
    complex(DP), dimension(size(AA,dim=1),size(AA,dim=1)) :: coeff
    
    ! parameters for lapack routine (optimal?)
    !! lwork=410 for celeron laptop
    !! lwork=570 for PentiumD 64-bit desktop
    !! lwork=330 for mac desktop
    integer, parameter :: lwork = 330
    complex(DP), dimension(lwork) :: work
    real(DP), dimension(2*size(AA,dim=1)) :: rwork
    integer :: info

    !! used in stratton normalization only
    complex(DP), dimension(size(AA,dim=1)) :: w

    if(.not. allocated(mcn)) allocate(mcn(size(AA,dim=1)*4))
    mcn = cmplx(0.0,0.0,DP)

    ! check to provide helpful message when I forget to allocate arrays
    if(size(AA,dim=1) < 2 .or. size(BB,dim=1) < 2) then
       print *, 'A and B shared arrays must be allocated&
            & before calling mcn_eigenvalues'
       stop
    end if
    
    M = size(AA,dim=1)
    
    di = 1 ! dummy integer
    dc(1) = CZERO ! dummy complex

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of even order
    ! De_{2n} in eq 3.12 of Stamnes & Spjelkavik, 
    ! Pure and Applied Optics, 4(3), 251-262, 1995
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Coeff(:,:) = CZERO

    ! main diagonal/ r counting from 0:m-1 like McLachlan
    forall(i=1:M-1) Coeff(i+1,i+1) = cmplx((2*i)**2, 0, DP)
    
    ! off diagonals
    forall(i=1:m, j=1:m, j==i+1 .or. j==i-1) Coeff(i,j) = q
    
    ! special case
    Coeff(2,1) = 2.0_DP*q

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
    
    call zgeev('N','V',M,Coeff,M,mcn(1:m),dc,di,AA(1:m,0:m-1,0),M,&
         & work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of even order"

    if (norm == 1) then 
       ! McLachlan, defined in terms of norm which has complex counterpart
       ! use 'dot product' to compute length^2, divide by length
       
       ! ZGEEV returns eigenvectors of unit euclidian length already
       ! only the first eigenvector must be re-scaled, then

       AA(1:m,0,1) = AA(1:m,0,1)/sqrt(2.0_DP*AA(1,0,0)*conjg(AA(1,0,0)) + &
            & sum(AA(2:m,0,0)*conjg(AA(2:m,0,0))))

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(AA(i,i-1,0)) < 0.0) then
             AA(:,i-1,0) = -AA(:,i-1,0)
          end if
       end do

    elseif (norm == 2) then 
       ! Morse norm
       ! ce_2n(psi=0)
       w = sum(AA(:,:,0),dim=1)
       
       ! normalize so |ce_2n(psi=0)| = +1
       AA(:,:,0) = AA(:,:,0)/ &
            &spread(abs(w)*signz(w),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! even coefficients (a) of odd order
    ! De_{2n+1} in eqn 3.14 of St&Sp
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

    call zgeev('N','V',M,Coeff,M,mcn(m+1:2*m),dc,di,AA(1:m,0:m-1,1),M,&
         &work,lwork,rwork,info)
         
    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating even coefficients of odd order"

    if (norm == 1) then ! McLachlan 

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(AA(i,i-1,1)) < 0.0) then
             AA(:,i-1,1) = -AA(:,i-1,1)
          end if
       end do
    elseif (norm == 2) then ! Stratton
       ! se'_2n+1(psi=0)
       w = sum(spread(2*[(i,i=0,m)]+1,dim=2,ncopies=m)*AA(:,:,1),dim=1)
       
       ! normalize so |se'_2n+1(psi=0)|=+1
       AA(:,:,1) = AA(:,:,1)/ &
            &spread(abs(w)*signz(w),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of even order
    ! Do_{2n+2} in eq 3.16 of St&Sp
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
   
    call zgeev('N','V',M,Coeff,M,mcn(2*m+1:3*m),dc,di,BB(1:m,0:m-1,0),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of even order"

    if (norm == 1) then ! McLachlan

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(BB(i,i-1,0)) < 0.0) then
             BB(:,i-1,0) = -BB(:,i-1,0)
          end if
       end do
    elseif (norm == 2) then ! Stratton
       ! ce_2n+2(psi=0)
       w = sum(spread(2*[(i,i=1,m)],dim=2,ncopies=m)*BB(:,:,0),dim=1)

       ! normalize so |ce_2n+2(psi=0)| = +1
       BB(:,:,0) = BB(:,:,0)/ &
            & spread(abs(w)*signz(w),dim=1,ncopies=m)
    end if

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! odd coefficients (b) of odd order
    ! Do_{2n+1} of eqn 3.18 in St&Sp
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

    call zgeev('N','V',M,Coeff,M,mcn(3*m+1:4*m),dc,di,BB(1:m,0:m-1,1),M,&
         &work,lwork,rwork,info)

    if (info /= 0) write (*,*) "ZGEEV ERROR ",info, &
         &" calculating odd coefficients of odd order"

    if (norm == 1) then ! McLachlan

       !! McLachlan has ambiguous sign, make real part positive
       do i = 1,m
          if (real(BB(i,i-1,1)) < 0.0) then
             BB(:,i-1,1) = -BB(:,i-1,1)
          end if
       end do

    elseif (norm == 2) then ! Stratton
       ! ce_2n+1(psi=0)
       w = sum(BB(:,:,1),dim=1)

       BB(:,:,1) = BB(:,:,1)/ &
            & spread(abs(w)*signz(w),dim=1,ncopies=m)
    end if

  end subroutine mcn_eigenvalues
  
  elemental function signz(z) result(s)
    ! sign built-in function only works for reals & ints
    use constants, only : DP
    complex(DP), intent(in) :: z
    complex(DP) :: s
    
    s = cmplx(sign(1.0_DP,real(z)),sign(1.0_DP,aimag(z)),DP)
  end function signz
  
  
end module mcn_matrix_method

