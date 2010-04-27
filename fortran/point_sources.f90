program source_terms_test
  ! $Id: source_terms_test.f90,v 1.4 2007/09/22 16:34:53 kris Exp kris $
  
  use constants, only : DP, PI
  use inverse_Laplace_Transform
  use complex_bessel, only : cbesk
  use expint_approx, only : expint
  implicit none 

!!$  interface tanh
!!$     module procedure tanh_double_complex
!!$  end interface
!!$  
!!$  interface coth
!!$     module procedure coth_double_complex
!!$  end interface
  
  character(30) :: outfname
  integer :: i, j, k, nz,ierr, np, kk
  real(DP), parameter :: DHTOL = 1.0D-13, DHALPHA = 0.0_DP
  integer, parameter :: NTERMS = 1500, MAXITER = 10000
  integer, parameter :: NUMT = 201, DHM = 100, NTYPES = 1, MINLOGT = -8, MAXLOGT = 8
  real(DP) :: r, Q, K1, Ss1, b2, K2, Ss2, Kz1,Sy, tau, alpha1, alpha2, tee, dellogt, &
       & Hc, sigma, delta, x, boulton, rw 
  real(DP), dimension(NTERMS) :: root, guess, gamma
  real(DP), dimension(NUMT) :: t, td
  real(DP), dimension(NUMT,0:NTYPES) :: ft
  real(DP), dimension(NUMt,2) :: theis
  complex(DP), dimension(2*DHM+1,0:NTYPES) ::  fp
  complex(DP) :: arg, kernel
  complex(DP), dimension(0:2) :: bk

  complex(DP), dimension(2*DHM+1) :: p, kappa1, kappa2 
  np = size(p)

  !! compute evenly-spaced log t
  dellogt = real(MAXLOGT - MINLOGT + 1,DP)/real(NUMT-1,DP)
  t(1:NUMT) = 10.0_DP**((/(MINLOGT + i*dellogt, i=0,NUMT-1)/))

  open(unit=23,file='sources.in',action='read',status='old')
  read(23,*) r, rw
  read(23,*) Q               !! pumping rate
  read(23,*) K1, Ss1           !! aquifer properties
  read(23,*) b2, K2, Ss2      !! aquitard properties
  read(23,*) Kz1,Sy           !! unconfined properties
  read(23,*) tau             !! time constant
  read(23,*) outfname
  close(23)

  if(r < rw) then
     print *, 'observation point inside well, moving to r_w'
     r = rw
  end if

  sigma = sqrt(Ss1/Sy)
  alpha1 = K1/Ss1
  alpha2 = K2/Ss2
  boulton = 3.0_DP*Kz1/Sy  !! Boulton's alpha interms of physical parameters
 
  do i=1,NUMT
     !! compute new vector of p for each value of t (optimal results but slow)
     call invlap_setup(p(1:np),tee,DHM,DHALPHA,DHTOL,t(i))
     
     kappa1(1:np) = sqrt(p(1:np)/alpha1)
     kappa2(1:np) = sqrt(p(1:np)/alpha2)
     
     do j=1,np
        
        !! normal finite radius source
        call cbesk(r*kappa1(j), 0.0_DP, 1, 1, bk(0), nz, ierr)
        call cbesk(rw*kappa1(j), 1.0_DP, 1, 2, bk(1), nz, ierr)
        call cbesk(rw*kappa1(j), 0.0_DP, 1, 1, bk(2), nz, ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 1', ierr
        end if
        
        fp(j,0) = bk(0)/(bk(1)*kappa1(j)*rw)

        !! wellbore storage solution
        fp(j,1) = bk(0)/(bk(1)*kappa1(j)*rw*(1.0_DP + 0.5_DP*rw*p(j)/kappa1(j)*bk(2)/bk(1)))


        !! leaky source (no drawdown at top of aquitard)
        arg = p(j)/alpha1 + kappa2(j)*K2/(K1)*coth(kappa2(j)*b2)
        call cbesk(r*sqrt(arg), 0.0_DP, 1, 1, fp(j,1), nz, ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 1', ierr
        end if

        !! leaky source (impermeable at top of aquitard)
        arg = p(j)/alpha1 + kappa2(j)*K2/(K1)*tanh(kappa2(j)*b2)
        call cbesk(r*sqrt(arg), 0.0_DP, 1, 1, fp(j,2), nz, ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 2', ierr
        end if

        !! Boulton kernel for unconfined source
        arg = p(j)/alpha1 + Sy*p(j)*boulton/(K1*(boulton + p(j)))
        call cbesk(r*sqrt(arg),0.0_DP,1,1,fp(j,3),nz,ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 3', ierr
        end if

        !! Boulton unconfined source (Herrera infinite sum Kernel)
        !! guess is halfway between asymptotes of cot()
        
        !! kernel is not a function of p
        if(i == 1 .and. j == 1) then
           guess(2:NTERMS) = PI*(real((/(k, k=1,NTERMS-1)/)) + 0.5_DP)/sigma
           guess(1) = 1.7D0
           
           !! first root is hard to find with NR, 
           !! use TS approximation for tangent and re-arrange
           x = guess(1)
           NR1: do kk = 1,MAXITER
              delta = (x + (sigma - 1.0_DP/sigma)*(x*sigma + (x*sigma)**3/3.0_DP + &
                   & 2.0_DP*(sigma*x)**5/15.0_DP) + 17.0_DP*(x*sigma)**7/315.0_DP)/ &
                   & (1.0_DP - (1.0_DP/sigma - sigma)*(sigma + x**2*sigma**3 + &
                   & 2.0_DP*x**4*sigma**5/3.0_DP + 17.0_DP*x**6*sigma**7/45.0_DP))
                 x = x - delta
                 if (abs(delta) <= 1.0D-10) then
                    root(1) = x
                    exit NR1
                 end if
                 if(kk == MAXITER) print *, '1 failed to converge'
              end do NR1

           do k = 2, NTERMS
              x = guess(k)
              NR: do kk = 1,MAXITER
                 delta = (1.0_DP/tan(x*sigma) + (sigma - 1.0_DP/sigma)/x)/&
                      & (sigma/(sin(sigma*x)**2) + (sigma + 1.0_DP/sigma)/x**2)
                 x = x + delta
                 if (abs(delta) <= spacing(x)) then
                    root(k) = x
                    exit NR
                 end if
                 if(kk == MAXITER) print *, k,'failed to converge'
              end do NR
           end do
           
           !! for debugging
           open(unit=666,file='debug_herrera_roots')
           write(666,*) '# sigma:',sigma
           write(666,*) '# n,  guess,  root, f(root)'
           write(666,*) '# extra part of Herrera kernel:', &
                & sum(root(:)**2/(root(:) - 1.0_DP + sigma**2))
           do k=1,NTERMS
              write(666,*) k,guess(k),root(k), &
                   & root(k)/tan(root(k)*sigma)+(sigma-1.0_DP/sigma)
           end do
           close(666)
           !! assume b=1
           gamma(1:NTERMS) = Kz1*root(1:NTERMS)**2/Sy
        end if

        kernel = 2.0_DP*sum(gamma(1:NTERMS)/((root(1:NTERMS)**2 - 1.0_DP + sigma**2)* &
             & (p(j) + gamma(1:NTERMS))))
        
        arg = p(j)/alpha1 + p(j)*Sy/K1*kernel
        call cbesk(r*sqrt(arg),0.0_DP,1,1,fp(j,6),nz,ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 6', ierr
        end if
!!$
        !! wave eqn
        call cbesk(r*sqrt(p(j)/alpha1*(1.0_DP + tau*p(j))) ,0.0_DP,1,1,bk(0),nz,ierr)
        call cbesk(rw*sqrt(p(j)/alpha1*(1.0_DP + tau*p(j))) ,1.0_DP,1,1,bk(1),nz,ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 4', ierr
        end if
        
        fp(j,4) = bk(0)/(sqrt(p(j)/alpha1*(1.0_DP + tau*p(j)))*bk(1)*rw)
!!$
        !! leaky source (b -> \infty)
        arg = p(j)/alpha1 + kappa2(j)*K2/(K1)
        call cbesk(r*sqrt(arg), 0.0_DP, 1, 1, fp(j,5), nz, ierr)
        if(ierr /= 0 .and. ierr /= 3) then
           print *, 'BF error 2', ierr
        end if
!!$
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

!!$  !########################################
!!$  elemental function tanh_double_complex(z) result(fz)
!!$    !! definition in terms of real & imag components from
!!$    !! Baker, "Less Complex Elementary Functions"
!!$    complex(DP), intent(in) :: z
!!$    complex(DP) :: fz
!!$    real(DP) :: x, y
!!$
!!$    x = real(z)
!!$    y = aimag(z)
!!$
!!$    fz = cmplx(tanh(2.0_DP*x)/(1.0_DP + cos(2.0_DP*y)/cosh(2.0_DP*x)), &
!!$         & sin(2.0_DP*y)/(cosh(2.0_DP*x) + cos(2.0_DP*y)))
!!$
!!$  end function tanh_double_complex
!!$  
!!$  elemental function coth_double_complex(z) result(fz)
!!$    complex(DP), intent(in) :: z
!!$    complex(DP) :: fz
!!$    
!!$    fz = 1.0_DP/tanh_double_complex(z)
!!$
!!$  end function coth_double_complex
!!$  
end program source_terms_test
