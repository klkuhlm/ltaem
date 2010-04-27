program area_source_test
  ! $Id: area_sources_test.f90,v 1.2 2007/07/11 01:50:27 kris Exp kris $

  use constants, only : DP, PI
  use inverse_Laplace_Transform
  use complex_bessel, only : cbesk
  implicit none 

  
  integer :: i, j, k, nz,ierr, np
  real(DP), parameter :: DHTOL = 1.0D-13, DHALPHA = 0.0D+0 
  integer, parameter :: NUMX = 40, NUMY = NUMX, DHM = 30
  real(DP) :: r, Q, K1, Ss1, alpha1, tee, t, gamx, gamy, rw
  real(DP) :: xmin, ymin, delx, dely
  real(DP), dimension(NUMX) :: x
  real(DP), dimension(NUMY) :: y
  real(DP), dimension(NUMX,NUMY) :: ft
  complex(DP), dimension(2*DHM+1,NUMX,NUMY) ::  fp
  complex(DP), dimension(0:1) :: besk
  complex(DP) :: src
  logical :: inst

  complex(DP), dimension(2*DHM+1) :: p, kappa1

  xmin = -2.0_DP
  ymin = xmin
  delx = 1.0D-1
  dely = delx

  x(1:NUMX) = (/(xmin + delx*i,i=0,NUMX-1)/)
  y(1:NUMX) = (/(ymin + dely*i,i=0,NUMY-1)/)

  open(unit=23,file='area_source.in',action='read',status='old')
  read(23,*) t               !! time
  read(23,*) Q, rw               !! pumping rate, borehole diameter
  read(23,*) K1, Ss1         !! aquifer properties
  read(23,*) gamx, gamy      !! area source strengths
  read(23,*) inst            !! instantaneous or continuous source
  close(23)

  alpha1 = K1/Ss1

  np = size(p)
  call invlap_setup(p(1:np),tee,DHM,DHALPHA,DHTOL,t)
  kappa1(1:np) = sqrt(p(1:np)/alpha1)

  do i=1,NUMX
     do j=1,NUMY
        r = sqrt(x(i)**2 + y(j)**2)
        if (r < rw) r = rw

        do k=1,np
           !! normal source
           call cbesk(r*kappa1(k), 0.0_DP, 1, 1, besk(0), nz, ierr)
           call cbesk(rw*kappa1(k), 1.0_DP, 1, 1, besk(1), nz, ierr)
           if(ierr /= 0 .and. ierr /= 3) then
              print *, 'BF error 1', ierr, 'x:',x(i),'y:',y(j)
           end if
           
           src = (gamx*x(i) + gamy*y(j))/p(k)

           if (.not. inst) then
              src = src/p(k)
           end if

           fp(k,i,j) = Q*besk(0)/(rw*kappa1(k)*2.0_DP*PI*K1*p(k)*besk(1)) - src

        end do       
        ft(i,j) = deHoog_invlap(DHALPHA,DHTOL,t,tee,fp(1:np,i,j),DHM)
     end do
  end do
     

  !! write dimensionless results to table
  open(unit=44,file='area_source.out',action='write',status='replace')
  write(44,'(A)') '# area source effects'
  write(44,'(3(A,ES10.3))') '# r=',r,' Q=',q, ' t=',t
  write(44,'(3(A,ES10.3))') '# K1=',K1,' Ss1=',Ss1,' alpha1=',alpha1
  write(44,'(2(A,ES10.3),A,L1)') '# gamx=',gamx,' gamy=',gamy,' inst=',inst
  write(44,'(2(A,ES10.3),A,I3)') '# deHoog: alpha=',dhAlpha,' tol=',dhTol,' M=',dhM
  write(44,'(A)') '# '
  write(44,'(A)') '#      x             y             soln '

  do i=1,NUMX
     do j=1,NUMY
        write(44,'(3(ES15.7E3,1X))') x(i),y(j),ft(i,j)
     end do
  end do

contains

  !! consolidate calculation of laplace parameter vector
  subroutine invlap_setup(p,tee,M,alpha,tol,t)
    
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
  
end program area_source_test
