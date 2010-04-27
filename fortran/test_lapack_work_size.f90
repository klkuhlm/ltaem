program test

integer, parameter :: DP = selected_real_kind(14,300)
integer :: info
integer, parameter :: N=10, LDA=10, LDVL=N, LDVR=N
complex(DP),dimension(1) :: dc
complex(DP), dimension(LDA,N) :: A,B,VL,VR
complex(DP), dimension(N) :: W
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



call zgeev('N','V',N,A,LDA,W,VL,LDVL,VR,LDVR,work,-1,rwork,info)

print *, 'optimum work size:',real(work(1))

end program test
