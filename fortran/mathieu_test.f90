! $Id: mathieu_test.f90,v 1.7 2010/02/26 15:54:33 klkuhlm Exp klkuhlm $

program mathieu_test

  use constants, only : DP, PI
  use mathieu_functions
  use utility, only : linspace, logspace

  implicit none
  integer, parameter :: LENE = 1000, LENP = 1000, MAXORD = 8, NUMQ = 15
  character(3) :: chnum, chqi
  character(43) :: fmt
  character(13), parameter :: fmt1 = '(ES14.5E3,1x,'
  character(25), parameter :: fmt2 = '(ES14.5E3,1x,ES14.5E3,1x)'
  character(3), dimension(12) :: type
  real(DP), dimension(NUMQ) :: qr, qi
  type(mathieu) :: mv
  real(DP), dimension(LENP) :: etaM
  real(DP), dimension(LENE) :: psiM
  complex(DP), dimension(0:MAXORD,1:LENP,1:4) :: aout
  complex(DP), dimension(0:MAXORD,1:LENE,5:12) :: rout
  complex(DP) :: q
  integer :: i,j,k,qq,minord
  integer, dimension(0:MAXORD), parameter :: vi = [(i, i=0,MAXORD)]

  aout = -999.0; rout = -999.0

  type = ['ce ','Dce','se ','Dse','Ke ','DKe','Ko ','DKo','Ie ','DIe','Io ','DIo']
  write(chnum,'(I2.2)') MAXORD+1
  fmt = fmt1//chnum//fmt2//')'
  write(chnum,'(I2.2)') (MAXORD+1)**2/2+1

  etaM = (/0.0_DP, logspace(-3.0,log10(2.0),LENE-1)/)
  psiM = linspace(0.0_DP,PI,LENP)

  qr = linspace(0.0_DP, 40.0_DP, NUMQ)
  qr(1) = 1.0D-3
  qi = qr

  chqi(1:1) = '_'
  
  do qq = 1, NUMQ
     write(chqi(2:3),'(I2.2)') qq

     q = cmplx(qr(qq),qi(qq),DP)
     
     mv = mathieu_init(q=q, MM=max(16,2*MAXORD,int(2.75*abs(q))))

     print *, 'q:',mv%q,'MS:',mv%M

     ! compute angular functions
     !========================================
     aout(0:,1:LENP,1) = ce( mv,vi(0:),psiM(1:LENP))
     aout(0:,1:LENP,2) = Dce(mv,vi(0:),psiM(1:LENP))
     aout(1:,1:LENP,3) = se( mv,vi(1:),psiM(1:LENP))
     aout(1:,1:LENP,4) = Dse(mv,vi(1:),psiM(1:LENP))

     ! output modified angular mathieu functions for plotting
     !========================================
     do k = 1, 4
        if (k < 3) then
           minord = 0
        else
           minord = 1
        end if
        open(unit=20,file=trim(type(k))//chqi//'_bench.out',status='replace',action='write')
        write(20,'(2(A,ES11.3),A)') '# Mathieu functions for q= (', &
             & real(mv%q),',',aimag(mv%q),')'
        write(20,*) '# ',type(k)
        do j = 1, LENP
           write(20,fmt) psiM(j),(real(aout(i,j,k)),aimag(aout(i,j,k)),i=minord,MAXORD)
        end do
        close(20)
     end do

     ! compute radial functions
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     rout(0:,1:LENE,5) =  Ke( mv,vi(0:),etaM(1:LENE))
     rout(1:,1:LENE,7) =  Ko( mv,vi(1:),etaM(1:LENE))
     rout(0:,1:LENE,6) =  DKe(mv,vi(0:),etaM(1:LENE))
     rout(1:,1:LENE,8) =  DKo(mv,vi(1:),etaM(1:LENE))
     rout(0:,1:LENE,9 ) = Ie( mv,vi(0:),etaM(1:LENE))
     rout(1:,1:LENE,11) = Io( mv,vi(1:),etaM(1:LENE))
     rout(0:,1:LENE,10) = DIe(mv,vi(0:),etaM(1:LENE))
     rout(1:,1:LENE,12) = DIo(mv,vi(1:),etaM(1:LENE))
     
     ! output first and second kind modified radial mathieu functions for plotting
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     do k = 5, 12
        if (k==5 .or. k==6 .or. k==9 .or. k==10) then
           minord = 0
        else
           minord = 1
        end if               
        open(unit=20,file=trim(type(k))//chqi//'_bench.out',status='replace',action='write')
        write(20,'(2(A,ES11.3),A)') '# Mathieu functions for q= (', &
             & real(mv%q),',',aimag(mv%q),')'
        write(20,*) '# ',type(k)
        do j = 1, LENE
           write(20,fmt) etaM(j),(real(rout(i,j,k)),aimag(rout(i,j,k)),i=minord,MAXORD)
        end do
        close(20)
     end do

     ! output NORMALIZED first and second kind modified radial mathieu functions for plotting
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !! K-type radial functions normalized by their value at eta=0
     do k = 5, 8
        if (k==5 .or. k==6) then
           minord = 0
        else
           minord = 1
        end if               
        open(unit=20,file=trim(type(k))//chqi//'_norm_bench.out',&
             & status='replace',action='write')
        write(20,'(2(A,ES12.4),A)') '# Mathieu functions for q= (', &
             & real(mv%q),',',aimag(mv%q),')'
        write(20,*) '# ',type(k)
        do j = 1, LENE
           write(20,fmt) etaM(j),(real(rout(i,j,k)/rout(i,1,k)),&
                & aimag(rout(i,j,k)/rout(i,1,k)),i=minord,MAXORD)
        end do
        close(20)
     end do

     ! output line source normalized functions
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !! even order even K-Mathieu functions normalized by their derivative at eta=0
     open(unit=20,file='flux_source_bench'//chqi//'.out',status='replace',action='write')
     write(20,'(2(A,ES12.4),A)') '# Mathieu functions for q= (', &
          & real(mv%q),',',aimag(mv%q),')'
     write(20,*) '# line source functions Ke_{2n}(eta)/DKe_{2n}(0)'
     do j = 1, LENE
        write(20,fmt) etaM(j),(real(rout(2*i,j,5)/rout(2*i,1,6)),&
             & aimag(rout(2*i,j,5)/rout(2*i,1,5)),i=0,floor(MAXORD/2.0))
     end do
     close(20)
  end do

end program mathieu_test

