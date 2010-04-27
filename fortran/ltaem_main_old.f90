! $Id: ltaem_main.f90,v 1.22 2006/06/09 16:01:32 kris Exp $
program ltaem_main
  use constants, only     : DP, PI, RTWO, SMALL,CZERO
  use calc_shared_data, only : CalcCoeff, CalcX, CalcY, Presult
  use element_specs, only : BGinfname,BGcalc,BGparticle,BGnumt,BGt,PARnum,PARdt,PARti,PARtf,&
       & BGx,BGy,Bgt,INVm,INVtol,INValpha,PARint,CIn,CIm,CInum,CImatchTol,BGCoeffFName,BGcontour,&
       & BGnumx,BGnumy,INVsmooth,BGoutput,BGoutfname
  use file_ops, only : readinput, writeresults
  use error_handler, only : fileerror
  use inverse_Laplace_Transform
  use matching
  use calc_routines
  use circular_geometry, only : DistanceAngleCalcs
  use particle_integrate

  implicit none
  integer :: i, j, part
  integer, allocatable :: parnumdt(:)
  real(DP), allocatable :: head(:,:,:), velx(:,:,:), vely(:,:,:)
  complex(DP), allocatable :: coeff(:,:,:,:) 
  complex(DP), allocatable :: Gm(:,:,:)

  real(DP), allocatable :: logt(:)
  integer, allocatable :: nt(:), run(:)
  complex(DP), allocatable :: s(:,:)
  real(DP) :: tee
  integer :: ilogt, iminlogt, imaxlogt, lo, hi, ierr
  character(128) :: subname = 'MainProgram (ltaem_main)'

  real(DP), parameter :: MOST = 0.99_DP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! either specify here or ask for at prompt
  BGinfname = 'input'

  ! read in locations, options & data from input file
  call readInput(BGinfname)

!!$  BGparticle = .false.
!!$  BGcontour = .true.

  ! nudge times on 'edge' of logcycle down a tiny bit to increase accuracy
  where ((nint(BGt(1:BGnumt)) - BGt(1:BGnumt)) < SMALL)
     BGt(1:BGnumt) = BGt(1:BGnumt) - SMALL
  end where

  allocate(run(1:2*INVm+1),stat=ierr)
  if (ierr /= 0) stop 'error allocating: RUN'
  run = (/ (i, i=0,2*INVm) /)

  ! independent of most choices
  call DistanceAngleCalcs()

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (BGcalc) then

     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     if (BGparticle) then   ! particle tracking
        
        ! need to do something here about particles starting at t=0
        where (PARti(1:PARnum) < 1.0D-5)
           PARti(1:PARnum) = 1.0D-5
        end where
        
        iminlogt = minval(  floor(log10(PARti(1:PARnum))))
        imaxlogt = maxval(ceiling(log10(PARtf(1:PARnum))))

        ! to make it possible for particles / contours to share code...
        imaxlogt = imaxlogt + 1

        allocate(s(iminlogt:imaxlogt-1,2*INVm+1), nt(iminlogt:imaxlogt-1), stat=ierr)
        if (ierr /= 0) stop 'error allocating: S, NT'

        do ilogt = iminlogt, imaxlogt-1
           print *, ilogt
           nt(ilogt) = 1
           tee = min(10.0_DP**(ilogt + MOST),real(imaxlogt-1,DP))*RTWO
           s(ilogt,:) = cmplx(INValpha - log(INVtol)/(RTWO*tee), PI*run/tee, DP)
        end do
        
        deallocate(run, stat=ierr)
        if (ierr /= 0) stop 'error deallocating: RUN'
        
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     else   ! contour maps / hydrographs

        allocate(logt(1:BGnumt),stat=ierr)
        if (ierr /= 0) stop 'error allocating: LOGT'

        lo = 1
        logt(1:BGnumt) = log10(BGt(1:BGnumt))
        iminlogt =   floor(minval(logt(1:BGnumt)))
        imaxlogt = ceiling(maxval(logt(1:BGnumt)))

        allocate(s(iminlogt:imaxlogt-1,2*INVm+1), nt(iminlogt:imaxlogt-1), stat=ierr)
        if (ierr /= 0) stop 'error allocating: S, NT'

        do ilogt = iminlogt, imaxlogt-1
           nt(ilogt) = count(logt >= real(ilogt,DP) .and. logt < real(ilogt+1,DP))
           tee = maxval(BGt(lo:lo+nt(ilogt)-1))*RTWO
           s(ilogt,:) = cmplx(INValpha - log(INVtol)/(RTWO*tee), PI*run/tee, DP)
           lo = lo + nt(ilogt)
        end do
        deallocate(logt,run, stat=ierr)
        if (ierr /= 0) stop 'error deallocating: S, NT'
     end if

     allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
          &  Gm(4*CIn+2, 2*CIm, CInum), stat=ierr)
     if (ierr /= 0) stop 'error allocating: COEFF, GM, VN'
     
     coeff = CZERO
     
     ! calculate coefficients for each value of Laplace parameter
     ! ** common between particle tracking and contours/hydrographs **

     do ilogt = iminlogt,imaxlogt-1
        write (*,*) 'log t=',ilogt
        do j = 1,2*INVm+1
           if (nt(ilogt) > 0) then
              print *, j, s(ilogt,j)
              call circInverse(s(ilogt,j),Gm(1:4*CIn+2, 1:2*CIm, 1:CInum))
              call matchIter(Gm(1:4*CIn+2, 1:2*CIm, 1:CInum), &
                   & coeff(j,ilogt,0:4*CIn+1,1:CInum), s(ilogt,j),CImatchTol)
              ! use results from last values of s as a first guess for next s
              if (j < 2*INVm+1) then
                 coeff(j+1,ilogt,0:4*CIn+1,1:CInum) = coeff(j,ilogt,0:4*CIn+1,1:CInum)
              end if
           end if
        end do
     end do


     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! save coefficient matrix in binary form to file (more accurate?)
     open(UNIT=77, FILE=BGCoeffFName, STATUS='REPLACE', ACTION='WRITE',&
          &  FORM='UNFORMATTED', IOSTAT=ierr)
     if (ierr /= 0) call fileError(BGCoeffFName,ierr,subname,0)
     write(77) iminlogt,imaxlogt,coeff,nt,s
     close(77)

     deallocate(Gm, stat=ierr)
     if (ierr /= 0) stop 'error deallocating: GM'
     print *, 'matching finished'

  else ! do not re-calculate coefficients
     open(UNIT=77, FILE=BGCoeffFName, STATUS='OLD', ACTION='READ',&
          &  FORM='UNFORMATTED', IOSTAT=ierr)
     if (ierr /= 0) then
        ! go back and recalc if no restart file
        print *, 'error opening restart file, recalculating...'
        BGcalc = .true.
        goto 111
     end if
     
     ! read two numbers needed to allocate arrays
     read(77) iminlogt,imaxlogt
     allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
            & s(iminlogt:imaxlogt-1,2*INVm+1), nt(iminlogt:imaxlogt-1),stat=ierr)
     if (ierr /= 0) stop 'error allocating: COEFF, S, NT'

     ! read needed data calculated in previous run from file
     rewind(77)
     read(77)  iminlogt,imaxlogt,coeff,nt,s
     close(77)
     print *, 'matching results read from file'

  end if ! re-calculate coefficients

  call freeMatchMem()
  allocate(calcCoeff(2*INVm+1,0:4*CIn+1,CInum), stat=ierr)
  if (ierr /= 0) stop 'error allocating: CALCCOEFF'

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(BGparticle) then ! integrate along particle path
     allocate(parnumdt(PARnum), stat=ierr)
     if (ierr /= 0) stop 'error allocating: PARNUMDT'

     forall(i = 1:PARnum)
        parnumdt(i) = ceiling((PARtf(i) - PARti(i))/PARdt)
     end forall

     allocate(Presult(0:maxval(parnumdt),5,PARnum), stat=ierr)
     if (ierr /= 0) stop 'error allocating: PRESULT'

     Presult = 0.0_DP

     print *, 'dt(i):',parnumdt

     print *, 'coefficient matrix bounds'
     print *, 'upper: ',ubound(coeff)
     print *, 'lower: ',lbound(coeff)

     ! cycle over particles, integrating each
     do part = 1, PARnum
        select case (PARint(part))
        case (1)
           call rungekuttamerson(coeff,lbound(coeff,dim=2),part)
        case (2)
           call rungekutta(coeff,lbound(coeff,dim=2),part)
        case (4)
           call fwdEuler(coeff,lbound(coeff,dim=2),part)
        case default
           print *, 'invalid integration code', part, PARint(part)
        end select
     end do

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(BGcontour) then ! contour output (x,y locations are independent; e.g. outerproduct)

     allocate(head(BGnumx,BGnumy,BGnumt), &
          &   velx(BGnumx,BGnumy,BGnumt), &
          &   vely(BGnumx,BGnumy,BGnumt), stat=ierr)
     if (ierr /= 0) stop 'error allocating: HEAD, VELX, VELY'

     head = 2.0_DP

     do ilogt = iminlogt,imaxlogt-1
        print *, 'logcycle:',ilogt
        do j = 1,BGnumx
           write (*,'(1A)',advance='NO') '#'
           do i = 1,BGnumy
              lo = 1 + sum(nt(iminlogt:ilogt-1))
              hi = sum(nt(iminlogt:ilogt))

              ! pass x, y and coefficients through module variables
              calcX = BGx(j)
              calcY = BGy(i)
              calcCoeff = coeff(:,ilogt,:,:)

              head(j,i,lo:hi) = invlap(headCalc,BGt(lo:hi),INValpha,INVtol,INVm,INVsmooth)
              ! velx and vely share most of the calculations (velx _must_ be called first)
              velx(j,i,lo:hi) = invlap(velxCalc,BGt(lo:hi),INValpha,INVtol,INVm,INVsmooth)
              vely(j,i,lo:hi) = invlap(velyCalc,BGt(lo:hi),INValpha,INVtol,INVm,INVsmooth)
           end do
        end do
     end do
  else ! hydrograph output (x,y locations are in pairs; e.g. inner product)

     allocate(head(BGnumx,1,BGnumt), &
          &   velx(BGnumx,1,BGnumt), &
          &   vely(BGnumx,1,BGnumt), stat=ierr)
     if (ierr /= 0) stop 'error allocating: HEAD, VELX, VELY'
     
     do i = 1,BGnumx
        print *, 'location:',BGx(i),BGy(i)
        do ilogt = iminlogt,imaxlogt-1
           lo = 1 + sum(nt(iminlogt:ilogt-1))
           hi = sum(nt(iminlogt:ilogt))

           ! pass x, y and coefficients through module variables
           calcX = BGx(i)
           calcY = BGy(i)
           calcCoeff = coeff(:,ilogt,:,:)

           ! don't need second dimension of results matricies
           head(i,1,lo:hi) = invlap(headCalc,BGt(lo:hi),INValpha,INVtol,INVm)
           velx(i,1,lo:hi) = invlap(velxCalc,BGt(lo:hi),INValpha,INVtol,INVm)
           vely(i,1,lo:hi) = invlap(velyCalc,BGt(lo:hi),INValpha,INVtol,INVm)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file

  deallocate(coeff,calcCoeff,s,nt, stat=ierr)
  if (ierr /= 0) stop 'error deallocating: COEFF, CALCCOEFF, S, NT'
  call freeCalcMem()

  if (BGparticle) then 
     ! make sure correct output if doing particle tracking
     BGoutput = 4
  else
     ! reset back to gnuplot output if not particle tracking
     if (BGoutput == 4) BGoutput = 1
  end if

  ! call subroutine to write output to file
  call writeResults(head,velx,vely,BGx,BGy,BGt,BGoutput,BGoutfname,Presult)

end program ltaem_main

