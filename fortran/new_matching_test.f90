program new_matching_test
  ! $Id: new_matching_test.f90,v 1.1 2006/08/08 20:29:32 kris Exp kris $
  use constants, only     : DP, PI, RTWO, RONE
  use calc_shared_data, only : calcX, calcY, LINk, LINalpha, LINSs
  use file_ops, only : writeResults
  use inverse_Laplace_Transform
  use new_matching, only : new_match

  implicit none 

  integer :: i, j, k, jj
  integer, parameter :: NUM = 40, NUMT = 1, NUMFT = 20
  real(DP), dimension(NUMT) :: t
  real(DP), dimension(NUM,NUM,NUMT) :: head, bh, velx, vely, Presult
  real(DP), dimension(NUM) :: x, y
  integer, dimension(5) :: is
  character(1) :: chslice
  logical :: lslice = .false.
  character(128) :: BGoutfname = 'line_mathieu.out', bout='line_mathieu_bench_'

  y = [(1.0D-1*real(i,DP) - 2.0_DP, i=0,NUM-1 )]
  x = y ! + 1.0_DP
  x=x+0.05
  y=y+0.05

!!$  LINq = RONE

  ! k = 1 if Ss=1.0E-3, alpha=1.0E+3
  LINK = 1.0D1
  LINSs = 1.0D-3
  LINalpha = LINk/LINSs

  velx = 0.0; vely = 0.0; Presult = 0.0;
  t(1) = 2.0_DP
  
  bh = 0.0_DP
!!$  LINsf = 0.5_DP ! semi-focal distance
!!$  LINpsi0 = 0.0_DP ! radial thickness of "line"

  do j = 1,NUM
     write(*,'(A2,1x,f7.4)') 'x:',x(j)
     do i = 1,NUM

        ! pass x, y and coefficients through module variables
        calcX = x(j)
        calcY = y(i)

        head(j,i,:) = invlap(new_match, t, 1.0D-9, 1.0D-12, NUMFT, .false.)
     end do
  end do

  is = [ 1, 4, 7, 10, 13 ]

  ! set to false to get contour maps
  lslice = .false.

  if (lslice) then
     ! write out results as x-slices for plotting in 2D
     do jj = 1,5
        write(chslice,'(i1)') jj
        open(unit=22, file=trim(bout)//'x'//chslice//'.out', status='replace', action='write')
        write (22,*) '# instantaneous line source benchmark'

        do i = 1, NUMT
           write (22,*) '# t=',t(i)
           write (22,*) '# x-slice',jj,'at x=',x(is(jj))
           write (22,*) '# Y, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
           write (22,995) (y(j),bh(is(jj),j,i), bh(is(jj),j,i)-head(is(jj),j,i),j=1,num)
           write (22,'(//)')
        end do

        write(22,*) '# EOF'
        close(22)
     end do
     
     ! write out results as y-slices for plotting in 2D
     do jj = 1,5
        write(chslice,'(i1)') jj
        open(unit=22, file=trim(bout)//'y'//chslice//'.out', status='replace', action='write')
        write (22,*) '# instantaneous line source benchmark'

        do i = 1, NUMT
           write (22,*) '# t=',t(i)
           write (22,*) '# y-slice',jj,'at y=',y(is(jj))
           write (22,*) '# X, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
           write (22,995) (x(k),bh(k,is(jj),i),bh(k,is(jj),i)-head(k,is(jj),i),k=1,num)
           write (22,'(//)')
        end do

        write(22,*) '# EOF'
        close(22)
     end do

  else
     ! write out all results in x-y-z style, for gnuplot surface plotting
     open(unit=22, file='line_bench.out', status='replace', action='write')
     write (22,*) '# instantaneous line source benchmark'
     
     do i = 1, NUMT
        write (22,*) '# t=',t(i)
        write (22,*) '# X, Y, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
        write (22,992) ((x(k),y(j),bh(k,j,i),bh(k,j,i)-head(k,j,i),k=1,num),j=1,num)
        write (22,'(//)')
     end do
     
     write(22,*) '# EOF'
     close(22)
  end if
992 format(1X,ES13.5,1X,ES13.5,1X,ES22.14e3,1X,ES22.14e3)
995 format(1X,ES13.5,1X,ES22.14e3,1X,ES22.14e3)
  
  ! call subroutine to write output to file
  call writeResults(head,velx,vely,x,y,t,1,BGoutfname,Presult)

end program new_matching_test
