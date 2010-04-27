program line_ellipse_test
  ! $Id: line_ellipse_test.f90,v 1.12 2010/04/08 01:33:33 klkuhlm Exp klkuhlm $

  ! constants.f90
  use constants, only : DP, PI, RTWO, RONE, TWOPI

  ! file_ops.f90
  use file_ops, only : writeResults

  ! invlap.f90
  use inverse_Laplace_Transform, only : invlap => dehoog_invlap, deHoog_pvalues

  ! mathieu_line.f90
  use mathieu_line, only : pointhead, mathieulinehead, line_external
  use mathieu_functions, only : mathieu
  
  implicit none 

  type(mathieu), allocatable :: mf(:,:)
  type(line_external) :: lin
  integer :: i, j, k, npts, m, output,ierr

  !! inverse Laplace transform tolerance, alpha, and tee=T is scaling parameter
  real(DP) :: tol, dhal, tee

  integer :: numx, numy, numt, numft
  real(DP), allocatable :: t(:)  ! NUMT
  real(DP), allocatable :: head(:,:,:),velx(:,:,:),vely(:,:,:), bh(:,:,:)
  real(DP), allocatable :: bvx(:,:,:),bvy(:,:,:), Presult(:,:,:) ! NUMX,NUMY,NUMT
  real(DP), allocatable :: x(:) ! NUMX
  real(DP), allocatable :: y(:) ! NUMY

  character(128) :: BGoutfname, bout='bench_'
  logical :: hydrograph, point

  complex(DP), allocatable :: p(:) ! 2*numft+1

  open(unit=23,file='line.in',action='read',status='old')
  read(23,*) lin%match, lin%flux, point, lin%cartesian, output, hydrograph

  read(23,*) numx, numy, numt
  allocate(t(numt),x(numx),y(numy),stat=ierr)
  if (ierr /= 0) then
     write(*,'(A,3(I0,1X))') 'error allocating x,y,t',numx,numy,numt
     stop 1001
  end if
  
  allocate(head(numx,numy,numt),Presult(numx,numy,numt),stat=ierr)
  if (ierr /= 0) then
     write(*,'(A,3(I0,1X))') 'error allocating head,Presult',numx,numy,numt
     stop 1004
  end if

  if (.not. hydrograph) then
     allocate(velx(numx,numy,numt),vely(numx,numy,numt),stat=ierr)
     if (ierr /= 0) then
        write(*,'(A,3(I0,1X))') 'error allocating velx,vely',numx,numy,numt
        stop 1005
     end if
     
  end if
  
  if (point) then
     allocate(bh(numx,numy,numt),bvx(numx,numy,numt),bvy(numx,numy,numt),stat=ierr)
     if (ierr /= 0) then
        write(*,'(A,3(I0,1X))') 'error allocating bh,bvx,bvy',numx,numy,numt
        stop 1002
     end if   
  end if
  
  read(23,*) numft, tol, dhal
  allocate(p(2*numft+1),mf(2*numft+1,2),stat=ierr)
  if (ierr /= 0) then
     write(*,'(A,(I0,1X))') 'error allocating p,lin%mf',numft
     stop 1003
  end if

  read(23,*) lin%q, npts
  read(23,*) lin%gamma
  read(23,*) lin%k, lin%Ek
  read(23,*) lin%Ss, lin%ESs
  read(23,*) x(1:NUMX)
  read(23,*) y(1:NUMY)
  read(23,*) t(1:NUMT)
  read(23,*) lin%sf  !!, y0??
  read(23,*) lin%eta0, lin%n, lin%ms
  read(23,*) lin%leak, lin%k2, lin%ss2, lin%b2, lin%b
  read(23,*) lin%unconfined, lin%Sy, lin%kz
  read(23,*) BGoutFname
  close(23)

  bout(7:8+len_trim(BGoutFname)) = BGoutFname

  write(*,'(5(A,L1),A,I2)') 'match=',lin%match,' flux=',lin%flux,' point=',point, &
       &' rectangular=',lin%cartesian ,' hydrograph=',hydrograph,' output=',output

  ! some minor error checking
  if(lin%match .and. lin%eta0 < epsilon(1.0)) then
     stop 'need finite radius for matching'
  elseif(.not. lin%match .and. lin%eta0 > epsilon(1.0)) then
     lin%eta0 = 0.0_DP
     print *, 'resetting radius to zero for line source'
  end if
  if(.not. lin%match .and. abs(lin%gamma) > epsilon(0.0_DP)) then
     lin%gamma = 0.0_DP
     print *, 'resetting recharge reate to zero if not matching'
  end if
  
  lin%alpha = lin%K/lin%Ss   ! outer alpha
  lin%Ealpha = lin%Ek/lin%ESs
  Presult = 0.0

  lin%well = cmplx(1.0,1.0,DP)

  do m = 1,NUMT
     write(*,'(A,I3,1X,ES10.2)') 'time',m,t(m)
     tee = 2*t(m)
     p(:) = deHoog_pvalues(tee,DHAL,TOL,NUMFT)
     lin%first = .true.

     do i = 1,NUMY
        write(*,'(A2,1x,f7.4,1X)') 'y:',y(i)
        do j = 1,NUMX
!!$           write(*,'(A2,1x,f7.4,1X)') 'x:',x(j)

           if (.not. lin%cartesian) then
              lin%WW = cmplx(0.0_DP, TWOPI/real(NUMX,DP)*real(j-1,DP))
              y = 0.0
              x(j) = aimag(lin%WW)
           end if
         
!!$           write(*,'(A)',ADVANCE='no') '='
           ! line source potential
           lin%calc_head = .true.
           head(j,i,m) = invlap(DHAL,TOL,t(m),tee, mathieulinehead(p(:),mf,x(j),y(i),lin),NUMFT)
           lin%first = .false.

           if(point .and. lin%flux .and. (.not. lin%match)) then
              bh(j,i,m) = invlap(DHAL,TOL,t(m),tee, pointhead(p(:),x(j),y(i),lin,npts),NUMFT)
           end if

           if (.not. hydrograph) then
              ! line source x-flux
              lin%calc_velx = .true.
              lin%calc_head = .false.
              velx(j,i,m) = invlap(DHAL,TOL,t(m),tee, mathieulinehead(p(:),mf,x(j),y(i),lin),NUMFT)
              if(point .and. lin%flux .and. (.not. lin%match)) then        
                 bvx(j,i,m) =  invlap(DHAL,TOL,t(m),tee, pointhead(p(:),x(j),y(i),lin,npts),NUMFT)
              end if
              
              ! line source y-flux
              lin%calc_velx = .false.
              vely(j,i,m) = invlap(DHAL,TOL,t(m),tee, mathieulinehead(p(:),mf,x(j),y(i),lin),NUMFT)
              if(point .and. lin%flux .and. (.not. lin%match)) then
                 bvy(j,i,m) = invlap(DHAL,TOL,t(m),tee, pointhead(p(:),x(j),y(i),lin,npts),NUMFT)
              end if
           end if
        end do
        write(*,'(/)',advance='no')
     end do
  end do

  if(point .and. (.not. hydrograph)) then
     ! write out all results in x-y-z style, for gnuplot surface plotting
     open(unit=22, file=bout, status='replace', action='write')
     write (22,*) '# instantaneous line source benchmark'
     
     do i = 1, NUMT
        write (22,*) '# t=',t(i)
        write (22,*) '# X, Y, point_head, point_velx, point_vely, MF_head, MF_velx, MF_vely'
        do j=1,NUMY
           do k=1,NUMX
              write (22,992) x(k),y(j),bh(k,j,i),bvx(k,j,i),bvy(k,j,i),&
                   & head(k,j,i),velx(k,j,i),vely(k,j,i)
           end do
        end do
        write (22,'(//)')
     end do
     
     write(22,*) '# EOF'
     close(22)
  end if
     
992 format(2(1X,ES13.5),6(1X,ES22.14e3))
  
  ! call subroutine to write output to file
  if (hydrograph) then
     call writeResults(head,head,head,x,y,t,output,BGoutfname,Presult)
  else
     call writeResults(head,velx,vely,x,y,t,output,BGoutfname,Presult)
  end if
  
  deallocate(p,mf,t,x,y,head,presult)
  if (.not. hydrograph) deallocate(velx,vely)
  if (point) deallocate(bh,bvx,bvy)

end program line_ellipse_test
