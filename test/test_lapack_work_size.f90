program test

integer, parameter :: DP = selected_real_kind(14,300)
integer :: info
!!integer, parameter :: N=10, LDA=10, LDVL=N, LDVR=N
integer:: M, N, LDA, LDVL, LDVR
complex(DP), dimension(1) :: dc
complex(DP), allocatable :: A(:,:),B(:,:),VL(:,:),VR(:,:)
complex(DP), allocatable :: W(:)
complex(DP), dimension(100) :: work

!##      SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
!##     $                  WORK, LWORK, RWORK, INFO )
!##*
!##*  -- LAPACK driver routine (version 3.2) --
!##*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!##*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!##*     November 2006
!##*
!##*     .. Scalar Arguments ..
!##      CHARACTER          JOBVL, JOBVR
!##      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!##*     ..
!##*     .. Array Arguments ..
!##      DOUBLE PRECISION   RWORK( * )
!##      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
!##     $                   W( * ), WORK( * )
!##*     ..
!##*
!##*  Purpose
!##*  =======
!##*
!##*  ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
!##*  eigenvalues and, optionally, the left and/or right eigenvectors.
!##*
!##*  The right eigenvector v(j) of A satisfies
!##*                   A * v(j) = lambda(j) * v(j)
!##*  where lambda(j) is its eigenvalue.
!##*  The left eigenvector u(j) of A satisfies
!##*                u(j)**H * A = lambda(j) * u(j)**H
!##*  where u(j)**H denotes the conjugate transpose of u(j).
!##*
!##*  The computed eigenvectors are normalized to have Euclidean norm
!##*  equal to 1 and largest component real.
!##*
!##*  Arguments
!##*  =========
!##*
!##*  JOBVL   (input) CHARACTER*1
!##*          = 'N': left eigenvectors of A are not computed;
!##*          = 'V': left eigenvectors of are computed.
!##*
!##*  JOBVR   (input) CHARACTER*1
!##*          = 'N': right eigenvectors of A are not computed;
!##*          = 'V': right eigenvectors of A are computed.
!##*
!##*  N       (input) INTEGER
!##*          The order of the matrix A. N >= 0.
!##*
!##*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!##*          On entry, the N-by-N matrix A.
!##*          On exit, A has been overwritten.
!##*
!##*  LDA     (input) INTEGER
!##*          The leading dimension of the array A.  LDA >= max(1,N).
!##*
!##*  W       (output) COMPLEX*16 array, dimension (N)
!##*          W contains the computed eigenvalues.
!##*
!##*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
!##*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!##*          after another in the columns of VL, in the same order
!##*          as their eigenvalues.
!##*          If JOBVL = 'N', VL is not referenced.
!##*          u(j) = VL(:,j), the j-th column of VL.
!##*
!##*  LDVL    (input) INTEGER
!##*          The leading dimension of the array VL.  LDVL >= 1; if
!##*          JOBVL = 'V', LDVL >= N.
!##*
!##*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
!##*          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!##*          after another in the columns of VR, in the same order
!##*          as their eigenvalues.
!##*          If JOBVR = 'N', VR is not referenced.
!##*          v(j) = VR(:,j), the j-th column of VR.
!##*
!##*  LDVR    (input) INTEGER
!##*          The leading dimension of the array VR.  LDVR >= 1; if
!##*          JOBVR = 'V', LDVR >= N.
!##*
!##*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
!##*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!##*
!##*  LWORK   (input) INTEGER
!##*          The dimension of the array WORK.  LWORK >= max(1,2*N).
!##*          For good performance, LWORK must generally be larger.
!##*
!##*          If LWORK = -1, then a workspace query is assumed; the routine
!##*          only calculates the optimal size of the WORK array, returns
!##*          this value as the first entry of the WORK array, and no error
!##*          message related to LWORK is issued by XERBLA.
!##*
!##*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)
!##*
!##*  INFO    (output) INTEGER
!##*          = 0:  successful exit
!##*          < 0:  if INFO = -i, the i-th argument had an illegal value.
!##*          > 0:  if INFO = i, the QR algorithm failed to compute all the
!##*                eigenvalues, and no eigenvectors have been computed;
!##*                elements and i+1:N of W contain eigenvalues which have
!##*                converged.
!##*


do n=30,600,30
   lda = n
   ldvl = n
   ldvr = n

   if (allocated(A)) deallocate(A,B,VL,VR,W)
   allocate(a(n,n),b(n,n),vl(n,n),vr(n,n),w(n))
   
   a = cmplx(2.0,3.0,DP); b = cmplx(3.0,-4.0,DP); 

   call zgeev('N','V',N,A,LDA,W,VL,LDVL,VR,LDVR,work,-1,rwork,info)
   print '(3(A,I0))', 'optimum ZGEEV work size: n=',n,' iwork=',nint(real(work(1))),' factor:',nint(real(work(1))/n)
end do


!!$!!      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
!!$!!     $                  INFO )
!!$!!*
!!$!!*  -- LAPACK driver routine (version 3.2) --
!!$!!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!$!!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!$!!*     November 2006
!!$!!*
!!$!!*     .. Scalar Arguments ..
!!$!!      CHARACTER          TRANS
!!$!!      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!!$!!*     ..
!!$!!*     .. Array Arguments ..
!!$!!      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!!$!!*     ..
!!$!!*
!!$!!*  Purpose
!!$!!*  =======
!!$!!*
!!$!!*  ZGELS solves overdetermined or underdetermined complex linear systems
!!$!!*  involving an M-by-N matrix A, or its conjugate-transpose, using a QR
!!$!!*  or LQ factorization of A.  It is assumed that A has full rank.
!!$!!*
!!$!!*  The following options are provided:
!!$!!*
!!$!!*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!!$!!*     an overdetermined system, i.e., solve the least squares problem
!!$!!*                  minimize || B - A*X ||.
!!$!!*
!!$!!*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!!$!!*     an underdetermined system A * X = B.
!!$!!*
!!$!!*  3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
!!$!!*     an undetermined system A**H * X = B.
!!$!!*
!!$!!*  4. If TRANS = 'C' and m < n:  find the least squares solution of
!!$!!*     an overdetermined system, i.e., solve the least squares problem
!!$!!*                  minimize || B - A**H * X ||.
!!$!!*
!!$!!*  Several right hand side vectors b and solution vectors x can be
!!$!!*  handled in a single call; they are stored as the columns of the
!!$!!*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!!$!!*  matrix X.
!!$!!*
!!$!!*  Arguments
!!$!!*  =========
!!$!!*
!!$!!*  TRANS   (input) CHARACTER*1
!!$!!*          = 'N': the linear system involves A;
!!$!!*          = 'C': the linear system involves A**H.
!!$!!*
!!$!!*  M       (input) INTEGER
!!$!!*          The number of rows of the matrix A.  M >= 0.
!!$!!*
!!$!!*  N       (input) INTEGER
!!$!!*          The number of columns of the matrix A.  N >= 0.
!!$!!*
!!$!!*  NRHS    (input) INTEGER
!!$!!*          The number of right hand sides, i.e., the number of
!!$!!*          columns of the matrices B and X. NRHS >= 0.
!!$!!*
!!$!!*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!!$!!*          On entry, the M-by-N matrix A.
!!$!!*            if M >= N, A is overwritten by details of its QR
!!$!!*                       factorization as returned by ZGEQRF;
!!$!!*            if M <  N, A is overwritten by details of its LQ
!!$!!*                       factorization as returned by ZGELQF.
!!$!!*
!!$!!*  LDA     (input) INTEGER
!!$!!*          The leading dimension of the array A.  LDA >= max(1,M).
!!$!!*
!!$!!*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!!$!!*          On entry, the matrix B of right hand side vectors, stored
!!$!!*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!!$!!*          if TRANS = 'C'.
!!$!!*          On exit, if INFO = 0, B is overwritten by the solution
!!$!!*          vectors, stored columnwise:
!!$!!*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!!$!!*          squares solution vectors; the residual sum of squares for the
!!$!!*          solution in each column is given by the sum of squares of the
!!$!!*          modulus of elements N+1 to M in that column;
!!$!!*          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!!$!!*          minimum norm solution vectors;
!!$!!*          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
!!$!!*          minimum norm solution vectors;
!!$!!*          if TRANS = 'C' and m < n, rows 1 to M of B contain the
!!$!!*          least squares solution vectors; the residual sum of squares
!!$!!*          for the solution in each column is given by the sum of
!!$!!*          squares of the modulus of elements M+1 to N in that column.
!!$!!*
!!$!!*  LDB     (input) INTEGER
!!$!!*          The leading dimension of the array B. LDB >= MAX(1,M,N).
!!$!!*
!!$!!*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
!!$!!*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!$!!*
!!$!!*  LWORK   (input) INTEGER
!!$!!*          The dimension of the array WORK.
!!$!!*          LWORK >= max( 1, MN + max( MN, NRHS ) ).
!!$!!*          For optimal performance,
!!$!!*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!!$!!*          where MN = min(M,N) and NB is the optimum block size.
!!$!!*
!!$!!*          If LWORK = -1, then a workspace query is assumed; the routine
!!$!!*          only calculates the optimal size of the WORK array, returns
!!$!!*          this value as the first entry of the WORK array, and no error
!!$!!*          message related to LWORK is issued by XERBLA.
!!$!!*
!!$!!*  INFO    (output) INTEGER
!!$!!*          = 0:  successful exit
!!$!!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$!!*          > 0:  if INFO =  i, the i-th diagonal element of the
!!$!!*                triangular factor of A is zero, so that A does not have
!!$!!*                full rank; the least squares solution could not be
!!$!!*                computed.
!!$!!*
!!$!!*  =====================================================================
!!$!!*
!!$!!

do n=30,600,30
   do m = 2*n,3*n,n
      if (allocated(A)) deallocate(A,B)
      allocate(a(m,n),b(m,1))
   
      a = cmplx(4.0,-119.9,DP); b = cmplx(1.0,1.0,DP); 

!!$!!      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
!!$!!     $                  INFO )

      call zgels('N',M,N,1,A,M,b,m,work,-1,info)
      print '(5(A,I0))', 'optimum ZGELS work size: n=',n,' m=',m,' iwork=',nint(real(work(1))),&
           & ' n-factor:',nint(real(work(1))/n),' m-factor:',nint(real(work(1))/m)
   end do
end do

end program test
