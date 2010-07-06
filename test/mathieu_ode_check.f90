!! $Id: mathieu_ode_check.f90,v 1.2 2007/07/28 19:16:00 kris Exp kris $
program mathieu_ode_check

  use constants, only : DP, PI, TWOPI, RTWO

#ifdef NORMAL
  use mathieu_functions_norm
#else
  use mathieu_functions
#endif

  use shared_mathieu, only : A,B,q,mcn
  use mcn_matrix_method, only : mcn_eigenvalues
  implicit none
  
  integer, parameter :: NUMQ = 10, NUMETA = 10, NUMPSI = 10, NUMORD = 5, MS = 45
  complex(DP), dimension(NUMQ) :: qq
  complex(DP), dimension(-1:1,0:NUMORD) :: fe, fo, Dfe, Dfo
  real(DP), dimension(NUMETA) :: eta
  real(DP), dimension(NUMPSI) :: psi
  
  character(2) :: chq,chord,chms
  character(32) :: fmt1,fmt2,fmt3, fmt4
  real(DP) :: h, sigma, step, lopsi,hipsi
  real(DP) :: mineta, maxeta, deleta, delpsi, sign
  integer :: i,j,m,r,mm
  integer, dimension(0:NUMORD) :: vi
  complex(DP), dimension(0:NUMORD) :: rese, reso, resDe, resDo
  
  vi = (/ (i,i=0,NUMORD) /)
  h = 1.0D-6
  sigma = 1.25D+1
  step = 5.0D-2

#ifdef DNORMAL
  qq(1:NUMQ) = -cmplx(sigma,(/(step*i, i=0,NUMQ-1)/)) 
  sign = +1.0_DP
#else
  qq(1:NUMQ) = cmplx(sigma,(/(step*i, i=0,NUMQ-1)/)) 
  sign = -1.0_DP
#endif


  write(chord,'(I2.2)') NUMORD+1
  fmt1 = '(ES13.6E1,2('//chord//'(2(1X,ES14.7E2))))'
  fmt2 = '(2('//chord//'(2(1X,ES14.7E2))))'

  !! numeta points distributed logarithmically from mineta to maxeta
  mineta = 1.0D-3
  maxeta = 2.5D+0
  deleta = (log10(maxeta) - log10(mineta))/real(NUMETA-1,DP)
!!$  eta(1) = 0.0D+0
  eta(1:NUMETA) = 1.0D+1**(/(log10(mineta)+ real(i,DP)*deleta,i=0,NUMETA-1)/)

  !! numpsi point distributed uniformly from -pi to pi
  delpsi = TWOPI/real(NUMPSI-1,DP)
  psi(1:NUMETA) = -PI + (/(delpsi*real(i,DP),i=0,NUMETA-1)/)

  allocate(A(1:MS,0:MS-1,1:2),B(1:MS,0:MS-1,1:2))

  do i=1,NUMQ
     write(chq,'(I2.2)') i
     open(unit=77,file='mathieu_ode_check_'//chq//'.out',action='write',status='replace')
     write(77,*) '# q=',qq(i), ' h=',h
     write(77,*) '# radial Mathieu functions, orders 0-',NUMORD
     write(77,'(A)') '#    eta       -----------Ie0--------------  -----------Ie1-------------- &
          & -----------Ie2--------------  -----------Ie3--------------  -----------Ie4-------------- &
          & -----------Ie5--------------  -----------Io0--------------  -----------Io1-------------- &
          & -----------Io2--------------  -----------Io3--------------  -----------Io4-------------- &
          & -----------Io5--------------  -----------Ke0--------------  -----------Ke1-------------- &
          & -----------Ke2--------------  -----------Ke3--------------  -----------Ke4-------------- &
          & -----------Ke5--------------  -----------Ko0--------------  -----------Ko1-------------- &
          & -----------Ko2--------------  -----------Ko3--------------  -----------Ko4-------------- &
          & -----------Ko5--------------'

     open(unit=33,file='Dmathieu_ode_check_'//chq//'.out',action='write',status='replace')
     write(33,*) '# q=',qq(i), ' h=',h
     write(33,*) '# derivatives of radial Mathieu functions, orders 0-',NUMORD
     write(33,'(A)') '#    eta       ----------DIe0--------------  ----------DIe1-------------- &
          & ----------DIe2--------------  ----------DIe3--------------  ----------DIe4-------------- &
          & ----------DIe5--------------  ----------DIo0--------------  ----------DIo1-------------- &
          & ----------DIo2--------------  ----------DIo3--------------  ----------DIo4-------------- &
          & ----------DIo5--------------  ----------DKe0--------------  ----------DKe1-------------- &
          & ----------DKe2--------------  ----------DKe3--------------  ----------DKe4-------------- &
          & ----------DKe5--------------  ----------DKo0--------------  ----------DKo1-------------- &
          & ----------DKo2--------------  ----------DKo3--------------  ----------DKo4-------------- &
          & ----------DKo5--------------'

     q = qq(i)
     call mcn_eigenvalues(qq(i),A(:,:,:),B(:,:,:),1)

     !! write out mathieu characteristic numbers
     open(unit=88,file='mathieu_mcn_check_'//chq//'.out')
     write(88,*) '# q=',qq(i)
     write(88,'(A)') '# A even (ce_2n)'
     do j=1,MS
        write(88,'(2(ES15.7E3,1X))') mcn(j)
     end do
     write(88,'(//A)') '# A odd (se_2n+1)'
     do j=MS+1,2*MS
        write(88,'(2(ES15.7E3,1X))') mcn(j)
     end do
     write(88,'(//A)') '# B even (se_2n+2)'
     do j=2*MS+1,3*MS
        write(88,'(2(ES15.7E3,1X))') mcn(j)
     end do     
     write(88,'(//A)') '# B odd (ce_2n+1)'
     do j=3*MS+1,4*MS
        write(88,'(2(ES15.7E3,1X))') mcn(j)
     end do     
     close(88)
    
     !! write out vectors of Mathieu coefficients 
     open(unit=99,file='mathieu_coefficients_'//chq//'.out')
     write(99,*) '# q=',qq(i)
     write(chms,'(I2.2)') MS
     fmt3 = '('//chms//'(2(ES15.7E3,1X)))'

     write(99,'(A)') '# A even (ce_2n)'
     do j=1,MS
        write(99,fmt3) A(1:MS,j-1,1)
     end do
     write(99,'(//A)') '# A odd (se_2n+1)'
     do j=1,MS
        write(99,fmt3) A(1:MS,j-1,2)
     end do
     write(99,'(//A)') '# B even (se_2n+2)'
     do j=1,MS
        write(99,fmt3) B(1:MS,j-1,1)
     end do
     write(99,'(//A)') '# B odd (ce_2n+1)'
     do j=1,MS
        write(99,fmt3) B(1:MS,j-1,2)
     end do
     close(99)

     !! values of first kind radial mathieu functions at end of interval
     open(unit=56,file='mathieu_fcn_vals_'//chq//'.out')
     write(56,*) '# q=',qq(i), ' h=',h
     write(56,'(A)') '#    eta       -----------Ie0--------------  -----------Ie1-------------- &
          & -----------Ie2--------------  -----------Ie3--------------  -----------Ie4-------------- &
          & -----------Ie5--------------  -----------Io0--------------  -----------Io1-------------- &
          & -----------Io2--------------  -----------Io3--------------  -----------Io4-------------- &
          & -----------Io5--------------'

     !! loop over radial locations
     do j=1,NUMETA
        
        !! radial functions (mod. RME)
        
        do r=1,2
           if(r==1) then
              !! even/odd first kind (Ie/Io)
              call mmatIeIo(vi,eta(j)-h,fe(-1,0:NUMORD),fo(-1,0:NUMORD))
              call mmatIeIo(vi,eta(j),  fe(0,0:NUMORD), fo(0,0:NUMORD))
              call mmatIeIo(vi,eta(j)+h,fe(+1,0:NUMORD),fo(+1,0:NUMORD))
              call mmatDIeDIo(vi,eta(j)-h,Dfe(-1,0:NUMORD),Dfo(-1,0:NUMORD))
              call mmatDIeDIo(vi,eta(j)+h,Dfe(+1,0:NUMORD),Dfo(+1,0:NUMORD))

              fmt4 = '(ES13.6E1,2('//chord//'(2(1X,ES14.7E2))))'

              if(j>=NUMETA-1) then  !! print values of function for debugging
                 write(56,fmt4) eta(j)-h,fe(-1,:),fo(-1,:)
                 write(56,fmt4) eta(j)  ,fe( 0,:),fo( 0,:)
                 write(56,fmt4) eta(j)+h,fe(+1,:),fo(+1,:)
                 write(56,fmt4) eta(j)-h,Dfe(-1,:),Dfo(-1,:)
                 write(56,fmt4) eta(j)+h,Dfe(+1,:),Dfo(+1,:)
              end if
              
           else
              !! even/odd second kind (Ke/Ko)
              call mmatKeKo(vi,eta(j)-h,fe(-1,0:NUMORD),fo(-1,0:NUMORD))
              call mmatKeKo(vi,eta(j),  fe(0,0:NUMORD), fo(0,0:NUMORD))
              call mmatKeKo(vi,eta(j)+h,fe(+1,0:NUMORD),fo(+1,0:NUMORD))
              call mmatDKeDKo(vi,eta(j)-h,Dfe(-1,0:NUMORD),Dfo(-1,0:NUMORD))
              call mmatDKeDKo(vi,eta(j)+h,Dfe(+1,0:NUMORD),Dfo(+1,0:NUMORD))
           end if
           
           do m=0,NUMORD
              
              !! even radial functions (Ie/Ke)
              if(mod(m,2) == 0) then  !! even order
                 mm = 1 + m/2  !! A even
              else  !! odd order
                 mm = 3*MS+1 + (m-1)/2  !! B odd
              end if
              
          
              !! compute residual in ODE using functions
              !! second derivative approximated using FD, O(h**2)
              rese(m) = (fe(1,m) - RTWO*fe(0,m) + fe(-1,m))/h**2 - &
                   & (mcn(mm) - sign*RTWO*qq(i)*cosh(RTWO*eta(j)))*fe(0,m)               
              
              !! use first derivative of analytic deriv of fcn
              !! first derivative approximated using central FD, O(h**2)
              resDe(m) = (Dfe(1,m) - Dfe(-1,m))/(RTWO*h) - &
                   & (mcn(mm) - sign*RTWO*qq(i)*cosh(RTWO*eta(j)))*fe(0,m)

              if(m > 0) then
                 
                 !! odd radial functions (Io/Ko)
                 if(mod(m,2) == 0) then  !! even order
                    mm = 2*MS+1 + (m-2)/2  !! B even
                 else  !! odd order
                    mm = MS+1 + (m-1)/2   !! A odd
                 end if
                 
                 reso(m) = (fo(1,m) - RTWO*fo(0,m) + fo(-1,m))/h**2 - &
                      & (mcn(mm) - sign*RTWO*qq(i)*cosh(RTWO*eta(j)))*fo(0,m) 
                 
                 resDo(m) = (Dfo(1,m) - Dfo(-1,m))/(RTWO*h) - &
                      & (mcn(mm) - sign*RTWO*qq(i)*cosh(RTWO*eta(j)))*fo(0,m)
              else
                 reso(0) = cmplx(0.0,0.0,DP)
                 resDo(0) = cmplx(0.0,0.0,DP)
              end if
              
           end do

           if(r==1) then
              write(77,fmt1,advance='no') eta(j),rese(0:NUMORD),reso(0:NUMORD)
              write(33,fmt1,advance='no') eta(j),resDe(0:NUMORD),resDo(0:NUMORD)
           else
              write(77,fmt2) rese(0:NUMORD),reso(0:NUMORD)
              write(33,fmt2) resDe(0:NUMORD),resDo(0:NUMORD)
           end if
        end do
     end do
     
     close(56)

     write(77,'(/)')  !! break between radial & angular functions
     write(77,*) '# angular Mathieu functions, orders 0-',NUMORD
     write(77,'(A)') '#    psi       -----------ce0--------------  -----------ce1-------------- &
          & -----------ce2--------------  -----------ce3--------------  -----------ce4-------------- &
          & -----------ce5--------------  -----------se0--------------  -----------se1-------------- &
          & -----------se2--------------  -----------se3--------------  -----------se4-------------- &
          & -----------se5--------------'
     write(33,'(/)')
     write(33,*) '# derivatives of angular Mathieu functions, orders 0-',NUMORD
     write(33,'(A)') '#    psi       ----------Dce0--------------  ----------Dce1-------------- &
          & ----------Dce2--------------  ----------Dce3--------------  ----------Dce4-------------- &
          & ----------Dce5--------------  ----------Dse0--------------  ----------Dse1-------------- &
          & ----------Dse2--------------  ----------Dse3--------------  ----------Dse4-------------- &
          & ----------Dse5--------------'

     do j=1,NUMPSI
        !! angular functions (mod. AME)
        
        !! keep argument in -pi < psi < pi range
        if((psi(j) - h) < -PI) then
           lopsi = (psi(j) - h) + TWOPI
           hipsi = psi(j) + h
        elseif((psi(j) + h) > PI) then
           hipsi = (psi(j) + h) - TWOPI
           lopsi = psi(j) - h
        else
           lopsi = psi(j) - h
           hipsi = psi(j) + h
        end if

        !! only first kind functions here (2nd kind are non-periodic)
        fe(-1,0:NUMORD) = mmatce(vi,lopsi)
        fe(0,0:NUMORD) =  mmatce(vi,psi(j))
        fe(+1,0:NUMORD) = mmatce(vi,hipsi)

        Dfe(-1,0:NUMORD) = mmatDce(vi,lopsi)
        Dfe(+1,0:NUMORD) = mmatDce(vi,hipsi)       

        fo(-1,0:NUMORD) = mmatse(vi,lopsi)
        fo(0,0:NUMORD) =  mmatse(vi,psi(j))
        fo(+1,0:NUMORD) = mmatse(vi,hipsi)

        Dfo(-1,0:NUMORD) = mmatDse(vi,lopsi)
        Dfo(+1,0:NUMORD) = mmatDse(vi,hipsi)

        do m=0,NUMORD
           
           !! even angular function (ce)
           if(mod(m,2) == 0) then  !! even order
              mm = 1 + m/2  !! A even
           else  !! odd order
              mm = 3*MS+1 + (m-1)/2  !! B odd
           end if
           
           !! compute residual in ODE
           rese(m) = (fe(1,m) - RTWO*fe(0,m) + fe(-1,m))/h**2 + &
                & (mcn(mm) - sign*RTWO*qq(i)*cos(RTWO*psi(j)))*fe(0,m)                 

           resDe(m) = (Dfe(1,m) - Dfe(-1,m))/(RTWO*h) + &
                & (mcn(mm) - sign*RTWO*qq(i)*cos(RTWO*psi(j)))*fe(0,m) 
           
           if(m > 0) then
              
              !! odd angular functions (se)
              if(mod(m,2) == 0) then  !! even order
                 mm = 2*MS+1 + (m-2)/2  !! B even
              else  !! odd order
                 mm = MS+1 + (m-1)/2   !! A odd
              end if
              
              reso(m) = (fo(1,m) - RTWO*fo(0,m) + fo(-1,m))/h**2 + &
                   & (mcn(mm) - sign*RTWO*qq(i)*cos(RTWO*psi(j)))*fo(0,m)

              resDo(m) = (Dfo(1,m) - Dfo(-1,m))/(RTWO*h) + &
                   & (mcn(mm) - sign*RTWO*qq(i)*cos(RTWO*psi(j)))*fo(0,m) 
           else
              reso(0) = cmplx(0.0,0.0,DP)
              resDo(0) = cmplx(0.0,0.0,DP)
           end if
           
        end do
        
        write(77,fmt1) psi(j),rese(0:NUMORD),reso(0:NUMORD)
        write(33,fmt1) psi(j),resDe(0:NUMORD),resDo(0:NUMORD)
        
     end do
     close(77)
     close(33)
  end do
  

end program mathieu_ode_check
