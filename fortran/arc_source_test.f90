program arc_source_test
  ! $Id: arc_source_test.f90,v 1.2 2007/07/13 18:38:08 kris Exp kris $
  use constants, only     : DP, PI, RTWO, RONE
  use file_ops, only : writeResults
  use inverse_Laplace_Transform
  use arc_test, only : arc_head, pointhead
  implicit none 

  integer :: i, j, k, jj
  real(DP), parameter :: TOL = 1.0D-10, DHAL = 1.0D-7  !! dehoog tolerance and alpha
  integer, parameter :: NUMX = 30, NUMY = NUMX, NUMT = 1, NUMFT = 20
  real(DP), dimension(NUMT) :: t
  real(DP), dimension(NUMX,NUMY,NUMT) :: head,velx,vely, Presult, bh,bvx,bvy
  real(DP), dimension(NUMX) :: x
  real(DP), dimension(NUMY) :: y
  integer, dimension(5) :: is
  character(1) :: chslice
  logical ::  passive, lslice = .false.
  real(DP) :: xmin, xmax, ymin, ymax, delx, dely, tee
  character(128) :: BGoutfname = 'arc_source.out', bout='arc_bench_'

  ! parameters of the line source
  real(DP) :: wellQ,arcQ, alpha, Ss, kbg, arcth1, arcth2, arcx1, arcy1, arcr1

  complex(DP), dimension(2*NUMFT+1) :: p  !!,fp

  xmin = -1.0_DP
  xmax = +1.1_DP
  ymin = -0.5_DP
  ymax = +1.1_DP

  !! setup a grid of observation points NUMxNUM
  delx = (xmax - xmin)/real(NUMX - 1,DP)
  dely = (ymax - ymin)/real(NUMY - 1,DP)

  forall(i=1:NUMX-1)
     x(i) = xmin + real(i-1,DP)*delx
  end forall
  x(NUMX) = xmax
  
  forall(i=1:NUMY-1)
     y(i) = ymin + real(i-1,DP)*dely
  end forall
  y(NUMY) = ymax

  open(unit=23,file='arc.in',action='read',status='old')
  read(23,*) passive
  read(23,*) wellQ, arcQ
  read(23,*) Kbg, SS
  read(23,*) arcr1, arcth1, arcth2
  read(23,*) arcx1, arcy1
  read(23,*) t(1)
  close(23)

  alpha = Kbg/Ss
  Presult = 0.0

  call invlap_setup(p,tee,NUMFT,DHAL,TOL,t(1))

  do j = 1,NUMX
     write(*,'(A2,1x,f7.4,1X)',advance='no') 'x:',x(j)
     do i = 1,NUMY

!!$        fp = arc_head(p(:),x(j),y(i),alpha,arcr1,arcth1,arcth2,.true.,.false., &
!!$             & wellQ,arcQ,kbg,passive,arcX1,arcY1)
!!$        open(unit=77,file='debug_fp.out')
!!$        do jj=1,size(p)
!!$           write(77,*) jj,p(jj),fp(jj) 
!!$        end do
!!$        close(77)

        ! line source potential
        head(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & arc_head(p(:),x(j),y(i),alpha,arcr1,arcth1,arcth2,.true.,.false., &
             & wellQ,arcQ,kbg,passive,arcX1,arcY1),NUMFT)

        if(passive) then
           bh(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & pointhead(p(:),x(j),y(i),alpha,kbg,arcq,arcr1,arcX1,arcY1, &
             & arcth1,arcth2,.true.,.false.),NUMFT)
        end if
        
        ! line source x-flux
        velx(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & arc_head(p(:),x(j),y(i),alpha,arcr1,arcth1,arcth2,.false.,.true., &
             & wellQ,arcQ,kbg,passive,arcX1,arcY1),NUMFT)

        if(passive) then
           bvx(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & pointhead(p(:),x(j),y(i),alpha,kbg,arcq,arcr1,arcX1,arcY1, &
             & arcth1,arcth2,.false.,.true.),NUMFT)
        end if

        ! line source y-flux
        vely(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & arc_head(p(:),x(j),y(i),alpha,arcr1,arcth1,arcth2,.false.,.false., &
             & wellQ,arcQ,kbg,passive,arcX1,arcY1),NUMFT)

        if(passive) then
           bvy(j,i,1:NUMT) = deHoog_invlap(DHAL,TOL,t(1:NUMT),tee,&
             & pointhead(p(:),x(j),y(i),alpha,kbg,arcq,arcr1,arcX1,arcY1, &
             & arcth1,arcth2,.false.,.false.),NUMFT)
        end if

     end do
     write(*,'(/)',advance='no')
  end do

  ! set to false to get contour maps
  lslice = .false.

  if (lslice) then

     ! indicies to make slices at (pseudo evenly spaced)
     is(1) = 1;  is(5) = NUMX
     is(2:4) = nint([(i*(NUMX-1.0)/4.0, i=1,3 )])

     ! write out results as x-slices for plotting in 2D
     do jj = 1,5
        write(chslice,'(i1)') jj
        open(unit=22, file=trim(bout)//'x'//chslice//'.out', &
             & status='replace', action='write')
        write (22,*) '# line source benchmarking -- slices along constant x'

        do i = 1, NUMT
           write (22,*) '# t=',t(i)
           write (22,*) '# x-slice',jj,'at x=',x(is(jj))
           write (22,*) '# Y, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
           write (22,995) (y(j),bh(is(jj),j,i), bvx(is(jj),j,i), bvy(is(jj),j,i),&
                & bh(is(jj),j,i)-head(is(jj),j,i), bvx(is(jj),j,i)-velx(is(jj),j,i), &
                & bvy(is(jj),j,i)-vely(is(jj),j,i),j=1,NUMY)
           write (22,'(//)')
        end do
        write(22,*) '# EOF'
        close(22)
     end do
     
     is(5) = NUMY
     is(2:4) = nint([(i*(NUMY-1.0)/4.0, i=1,3 )])

     ! write out results as y-slices for plotting in 2D
     do jj = 1,5
        write(chslice,'(i1)') jj
        open(unit=22, file=trim(bout)//'y'//chslice//'.out', &
             & status='replace', action='write')
        write (22,*) '# line source benchmarking -- slices along constant y'

        do i = 1, NUMT
           write (22,*) '# t=',t(i)
           write (22,*) '# y-slice',jj,'at y=',y(is(jj))
           write (22,*) '# X, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
           write (22,995) (x(k),bh(k,is(jj),i),bvx(k,is(jj),i),bvy(k,is(jj),i),&
                & bh(k,is(jj),i)-head(k,is(jj),i),bvx(k,is(jj),i)-velx(k,is(jj),i),&
                & bvy(k,is(jj),i)-vely(k,is(jj),i),k=1,NUMX)
           write (22,'(//)')
        end do
        write(22,*) '# EOF'
        close(22)
     end do

  else
     ! write out all results in x-y-z style, for gnuplot surface plotting
     open(unit=22, file='arc_bench.out', status='replace', action='write')
     write (22,*) '# passive arc source benchmark'
     
     do i = 1, NUMT
        write (22,*) '# t=',t(i)
        write (22,*) '# X, Y, bh, bvx, bvy, bh-head, bvx-velx, bvy-vely'
        do j=1,NUMY
           do k=1,NUMX
              write (22,992) x(k),y(j),bh(k,j,i),bvx(k,j,i),bvy(k,j,i),&
                   & bh(k,j,i)-head(k,j,i),bvx(k,j,i)-velx(k,j,i),&
                   & bvy(k,j,i)-vely(k,j,i)
           end do
        end do
        write (22,'(//)')
     end do
     
     write(22,*) '# EOF'
     close(22)
  end if

992 format(2(1X,ES13.5),6(1X,ES22.14e3))
995 format(1X,ES13.5,6(1X,ES22.14e3))
  
  ! call subroutine to write output to file
  call writeResults(head,velx,vely,x,y,t,1,BGoutfname,Presult)

contains

  !! consolidate calculation of laplace parameter vector
  subroutine invlap_setup(p,tee,M,alpha,tol,t)
    use constants, only : DP
    
    complex(DP), intent(out), dimension(2*M+1) :: p
    real(DP), intent(out) :: tee
    integer :: np,i
    integer, intent(in) :: M
    real(DP), dimension(0:2*M) :: run
    real(DP), intent(in) :: alpha,tol,t

    run = real((/ (i,i=0,2*M) /),DP)
    np = 2*M+1
    tee = 2.0_DP*t

    p(1:np) = cmplx(alpha - log(tol)/(2.0_DP*tee), run*PI/tee, DP)

    ! check for potential problems
    if(minval(abs(p)) < alpha) then
       stop 'decrease alpha parameter for de Hoog routine'
    end if
    
  end subroutine invlap_setup

end program arc_source_test
