	PROGRAM MCERROR
C
C       ============================================================
C       Purpose: This program computes the error function erf(z) 
C                for a complex argument using subroutine CERROR
C       Input :  x   --- Real part of z
C                y   --- Imaginary part of z  ( y ó 3.0 )
C       Output:  ERR --- Real part of erf(z)
C                ERI --- Imaginary part of erf(z)
C       Example:
C                   x       y       Re[erf(z)]      Im[erf(z)]
C                 ---------------------------------------------
C                  1.0     2.0      -.53664357     -5.04914370
C                  2.0     2.0      1.15131087       .12729163
C                  3.0     2.0       .99896328      -.00001155
C                  4.0     2.0      1.00000057      -.00000051
C                  5.0     2.0      1.00000000       .00000000
C       ============================================================
C
	IMPLICIT COMPLEX *16 (C,Z)  
	DOUBLE PRECISION X,Y
	WRITE(*,*)'X,Y=?'
	READ(*,*)X,Y
	WRITE(*,*)'   x      y      Re[erf(z)]      Im[erf(z)]'
	WRITE(*,*)' ---------------------------------------------'
	Z=CMPLX(X,Y)
	CALL CERROR(Z,CER)
	WRITE(*,10) Z,CER
10      FORMAT(1X,F5.1,2X,F5.1,1X,2E16.8)
	END


