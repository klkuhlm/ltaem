program double_circle
  ! $Id: source_terms_test.f90,v 1.4 2007/09/22 16:34:53 kris Exp kris $
  
  use constants, only : DP, PI
  use inverse_Laplace_Transform
  use complex_bessel, only : cbesk
  use expint_approx, only : expint
  implicit none 

  character(30) :: outfname
  integer :: nr, ws, ierr, nz, i,j,k
  real(DP), parameter :: DHTOL = 1.0D-13, DHALPHA = 0.0_DP, ORD = 10
  integer, parameter :: DHM = 100
  real(DP) :: t,td, rw, r1, r2, Q, tee, 
  real(DP), dimension(3) :: K,Ss, alpha
  real(DP), allocatable :: rob(:)
  real(DP), allocatable :: theis(:), circ(:)
  complex(DP), dimension(2*DHM+1) ::  fp
  complex(DP), dimension(2*DHM+1) :: p
  complex(DP), dimension(2*DHM+1,3) :: kappa
  complex(DP), dimension(2*DHM+1,ORD,3) :: bk,bi,dbk,dbi
  np = size(p)

  open(unit=23,file='sources.in',action='read',status='old')
  read(23,*) rw, r1, r2
  read(23,*) Q, ws               !! pumping rate
  read(23,*) K(1:3)
  read(23,*) Ss(1:3)
  read(23,*) nr
  allocate(rob(nr))
  read(23,*) rob(1:nr)
  read(23,*) t
  close(23)

  call invlap_setup(p(1:np),tee,DHM,DHALPHA,DHTOL,t)
     
  alpha = K/Ss
  do i=1,3
     kappa(1:np,i) = sqrt(p(1:np)/alpha(i))
  end do
  
     do j=1,np
        
        do k=1,3
           call cbesk(r*kappa(j,k), 0.0_DP, 1, ORD, bk(j,1:ORD,k), nz, ierr)
           if(ierr /= 0 .and. ierr /= 3) then
              print *, 'besk error ',k, ierr
           end if
           call cbesi(r*kappa(j,k), 0.0_DP, 1, ORD, bk(j,1:ORD,k), nz, ierr)
           if(ierr /= 0 .and. ierr /= 3) then
              print *, 'besk error ',k, ierr
           end if
           
        end do
        

     end do
     

     fp(1:np,0:NTYPES) = fp(1:np,0:NTYPES)*Q/ &
          & (2.0_DP*PI*K1*spread(p(1:np),dim=2,ncopies=NTYPES+1))

!!$     !! van Everdingen & Hurst solution for constant head source
!!$     fp(1:np,NTYPES) = 2.0_DP*PI*K1/(Q*fp(1:np,0)*p(1:np)**2)

     do k=0,NTYPES
        ft(i,k) = deHoog_invlap(DHALPHA,DHTOL,t(i),tee,fp(1:np,k),DHM)
     end do
     
  end do

  !! compute theis solutions directly (no Laplace transform)
  !! <<these are dimensionless>>
  do i = 1, numt
     theis(i,1) = expint(r**2*Ss1/(4.0_DP*K1*t(i)))   ! early Theis (Ss)
     theis(i,2) = expint(r**2*Sy/(4.0_DP*K1*t(i)))  ! late Theis (Sy)
  end do
!!$
  !! write dimensional results to table
  open(unit=44,file=trim(outfname)//'-dim',action='write',status='replace')
  write(44,'(A)') '# source terms'
  write(44,'(2(A,ES10.3))') '# r=',r,' Q=',q
  write(44,'(3(A,ES10.3))') '# K1=',K1,' Ss1=',Ss1,' alpha1=',alpha1
  write(44,'(4(A,ES10.3))') '# K2=',K2,' Ss2=',Ss2,' alpha2=',alpha2,' b=',b2
  write(44,'(2(A,ES10.3))') '# Kz1=',Kz1,' Sy=',Sy
  write(44,'(1(A,ES10.3))') '# tau=',tau
  write(44,'(2(A,ES10.3),A,I3)') '# deHoog: alpha=',dhAlpha,' tol=',dhTol,' M=',dhM
  write(44,'(A)') '# '
  write(44,'(A)') '#     t            orig            leak1        &
       &  leak2         unconf          wave          leak3         Theis_Ss       Theis_Sy'
  do i=1,NUMT
     write(44,'(10(ES18.10E3,1X))') t(i),(ft(i,k),k=0,NTYPES), theis(i,1), theis(i,2)
  end do
  

  !! write dimensionless results to table
  open(unit=44,file=trim(outfname),action='write',status='replace')
  write(44,'(A)') '# dimensionless source terms'
  write(44,'(2(A,ES10.3))') '# r=',r,' Q=',q
  write(44,'(3(A,ES10.3))') '# K1=',K1,' Ss1=',Ss1,' alpha1=',alpha1
  write(44,'(4(A,ES10.3))') '# K2=',K2,' Ss2=',Ss2,' alpha2=',alpha2,' b=',b2
  write(44,'(2(A,ES10.3))') '# Kz1=',Kz1,' Sy=',Sy
  write(44,'(1(A,ES10.3))') '# tau=',tau
  write(44,'(2(A,ES10.3),A,I3)') '# deHoog: alpha=',dhAlpha,' tol=',dhTol,' M=',dhM
  write(44,'(A)') '# '
  write(44,'(A)') '#     t              orig              leak1          &
       &  leak2           unconf            wave            leak3            herrera         Theis_Ss       Theis_Sy'

  Hc = 4.0_DP*PI*K1/Q
  td(1:NUMT) = K1*t(1:NUMT)/(Ss1*r**2)

  do i = 1, numt
     theis(i,1) = expint(0.25_DP/td(i))   ! early Theis (Ss)
     theis(i,2) = expint(Sy*r**2/(4.0_DP*Kz1*t(i)))  ! late Theis (Sy)
  end do

  do i=1,NUMT
     write(44,'(10(ES15.7E3,1X))') td(i), (ft(i,j)*Hc,j=0,NTYPES), theis(i,1), theis(i,2)
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


end program double_circle
