program plot_q_surface

  use constants, only : DP, PI
  use shared_mathieu, only : A,B,mcn
  use mcn_matrix_method, only : mcn_eigenvalues

  implicit none

  character(2) :: chOrd
  character(4) :: chq
  character(30) :: fmt1,fmt2
  integer, parameter :: NUMR = 41, NUMI = 21, MAXORD = 8, MS = MAXORD+10
  real(DP), parameter :: DELR = 5.0D-1, DELI = 5.0D-1
  real(DP), parameter :: MINR = -1.0D1, MINI = -5.0D0

  complex(DP), dimension(NUMR,NUMI) :: qq
  complex(DP), dimension(0:2*MAXORD+1,NUMR,NUMI) :: mcna
  complex(DP), dimension(1:2*MAXORD+2,NUMR,NUMI) :: mcnb

  integer :: i,j,k,ii

  qq = cmplx(spread((/(MINR + i*DELR,i=0,NUMR-1)/),dim=2,ncopies=NUMI),&
           & spread((/(MINI + j*DELI,j=0,NUMI-1)/),dim=1,ncopies=NUMR),DP)
  

  qq(1,1) = cmplx(0.0_DP,-4.8_DP,DP)  ! examples from McLachlan for comparison
  qq(1,2) = cmplx(0.0_DP,-0.16_DP,DP)
  qq(1,3) = cmplx(0.0_DP,+4.8_DP,DP)
  qq(1,4) = cmplx(0.0_DP,+0.16_DP,DP)
  
  !! double points of Mathieu's equation for investigation of
  !! numerical behavior 

  qq(1, 5) = cmplx( 0.000000000,  1.468768610,DP)   ! a_{0,2}
  qq(1, 6) = cmplx( 1.931392507,  3.237638408,DP)   ! a_{1,3}
  qq(1, 7) = cmplx( 5.174241347,  5.104232136,DP)   ! a_{2,4}
  qq(1, 8) = cmplx( 9.687743644,  7.044513756,DP)   ! a_{3,5}
  qq(1, 9) = cmplx( 0.000000000, 16.471165890,DP)   ! a_{4,6}a
  qq(1,10) = cmplx(15.456855678,  9.042769248,DP)   ! a_{4,6}b
  qq(1,11) = cmplx( 4.851579879, 22.334379185,DP)   ! a_{5,7}a
  qq(1,12) = cmplx(22.474249041, 11.088535471,DP)   ! a_{5,7}b
  qq(1,13) = cmplx(11.084594179, 28.336502137,DP)   ! a_{6.8}a
  qq(1,14) = cmplx(30.735777246, 13.174441722,DP)   ! a_{6,8}b
  
  qq(2, 5) = cmplx( 0.000000000,  6.928954760,DP)   ! b_{2,4}
  qq(2, 6) = cmplx( 3.396646467, 10.746995009,DP)   ! b_{3,5}
  qq(2, 7) = cmplx( 8.151749880, 14.693291045,DP)   ! b_{4,6}
  qq(2, 8) = cmplx(14.220476068, 18.752577549,DP)   ! b_{5,7}
  qq(2, 9) = cmplx( 0.000000000, 30.096772840,DP)   ! b_{6,8}a
  qq(2,10) = cmplx(21.579035313, 22.910519866,DP)   ! b_{6,8}b
  qq(2,11) = cmplx( 6.303513622, 38.003703337,DP)   ! b_{7,9}a
  qq(2,12) = cmplx(30.213105680, 27.155319753,DP)   ! b_{7,9}b
  qq(2,13) = cmplx(14.002094660, 46.054886860,DP)   ! b_{8,10}a  
  qq(2,14) = cmplx(40.113332889, 31.477416947,DP)   ! b_{8,10}b

  !! not used here
  allocate(A(1:MS,0:MS-1,1:2),B(1:MS,0:MS-1,1:2))

  fmt2 = '(  (ES10.2E2,",",ES10.2E2,3X))'
  write(fmt2(2:3),'(I2.2)') MS
  
  do i=1,NUMR
     write(chq(1:2),'(I2.2)') i
     do j=1,NUMI
        write(chq(3:4),'(I2.2)') j

        call mcn_eigenvalues(qq(i,j),A(:,:,:),B(:,:,:),1)
        
        open(unit=77,file='eigv_A_'//chq//'.out',&
             & action='write',status='replace')
        write(77,*) '# A even; q=',qq(i,j)
        do ii = 1,MS
           write(77,fmt2) A(ii,:,1)
        end do
        write(77,'(///)')
        write(77,*) '# A odd; q=',qq(i,j)
        do ii = 1,MS
           write(77,fmt2) A(ii,:,2)
        end do
        close(77)

        open(unit=77,file='eigv_B_'//chq//'.out',&
             & action='write',status='replace')
        write(77,*) '# B even; q=',qq(i,j)
        do ii = 1,MS
           write(77,fmt2) B(ii,:,1)
        end do
        write(77,'(///)')
        write(77,*) '# B odd; q=',qq(i,j)
        do ii = 1,MS
           write(77,fmt2) B(ii,:,2)
        end do

        mcna(0:2*MAXORD:2,i,j) = mcn(1:MAXORD+1) !! even a; ce_2n
        mcna(1:2*MAXORD+1:2,i,j) = mcn(MS+1:MS+MAXORD+1) !! odd a; se_2n+1    

        mcnb(2:2*MAXORD+2:2,i,j) = mcn(2*MS+1:2*MS+MAXORD+1) !! even b; se_2n+2
        mcnb(1:2*MAXORD+1:2,i,j) = mcn(3*MS+1:3*MS+MAXORD+1) !! odd b; ce_2n+1    
     end do
  end do

  fmt1 = '(   (ES14.6,1X))'
  write(fmt1(2:4),'(I3.3)') NUMR
  
  !! matlab-style plotting (coordinates)
  open(unit=44,file='mathieu_re_q.dat')
  do j=1,NUMI
     write(44,fmt1) real(qq(:,j))
  end do
  close(44)

  open(unit=44,file='mathieu_im_q.dat')
  do j=1,NUMI
     write(44,fmt1) aimag(qq(:,j))
  end do
  close(44)


  !! write output (matlab-style flat matricies)
  do i=0,2*MAXORD
     write(chOrd,'(I2.2)') i

     open(unit=22,file='mathieu_q_r_a'//chOrd//'.dat',&
          &action='write',status='replace')
!!     write(22,'(A)') '# Re{a mcn}'
     do j=1,NUMI
        write(22,fmt1) real(mcna(i,:,j))
     end do
     close(22)

     open(unit=22,file='mathieu_q_i_a'//chOrd//'.dat',&
          &action='write',status='replace')
!!     write(22,'(A)') '# Im{a mcn}'
     do j=1,NUMI
        write(22,fmt1) aimag(mcna(i,:,j))
     end do
     close(22)

     write(chOrd,'(I2.2)') i+1
     open(unit=22,file='mathieu_q_r_b'//chOrd//'.dat',&
          &action='write',status='replace')
!!     write(22,'(A)') '# Re{b mcn}'
     do j=1,NUMI
        write(22,fmt1) real(mcnb(i+1,:,j))
     end do
     close(22)

     open(unit=22,file='mathieu_q_i_b'//chOrd//'.dat',&
          &action='write',status='replace')
!!     write(22,'(A)') '# Im{b mcn}'
     do j=1,NUMI
        write(22,*) aimag(mcnb(i+1,:,j))
     end do
     close(22)
  end do
  
  !! write output gnuplot-style (xyz)
  open(unit=33,file='mathieu_q_gnuplot.dat')
  do k=0,2*MAXORD
     write(33,'(A,I2)') '# re(q), aimag(q), re(a), aimag(a), &
          &abs(a), re(b), aimag(b), abs(b),', k
     do i=1,NUMR
        do j=1,NUMI
           write(33,333) qq(i,j),mcna(k,i,j),abs(mcna(k,i,j)),&
                &mcnb(k+1,i,j),abs(mcnb(k+1,i,j))
        end do
     end do
     write(33,'(//)')
  end do
  
333 format(2(F7.3,1X),6(ES16.8,1X))

end program plot_q_surface
