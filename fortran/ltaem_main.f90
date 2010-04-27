! $Id: ltaem_main.f90,v 1.7 2008/12/10 02:46:44 kris Exp kris $
program ltaem_main
  use constants, only     : DP, PI, RTWO, SMALL,CZERO
  use calc_shared_data, only : Presult
  use element_specs, only : BGinfname,BGcalc,BGparticle,BGnumt,BGt,PARnum,PARdt,PARti,PARtf,&
       & BGx,BGy,Bgt,INVm,INVtol,INValpha,PARint,CIn,CInum,BGCoeffFName,BGcontour,&
       & BGnumx,BGnumy,BGoutput,BGoutfname, CIm, CIMatchTol, WLstorCoeff
  use file_ops, only : readinput, writeresults
  use error_handler, only : fileerror
  use inverse_Laplace_Transform, only : invlap => deHoog_invlap
!!$  use matching_new, only : circInverse_matrix
  use matching_old
  use calc_routines
  use circular_geometry, only : DistanceAngleCalcs
  use particle_integrate

  implicit none
  integer :: i, j, part, tnp
  integer, dimension(4) :: ic
  integer, allocatable :: parnumdt(:)
  real(DP), allocatable :: head(:,:,:), velx(:,:,:), vely(:,:,:)
  complex(DP), allocatable :: headp(:),velxp(:),velyp(:)
  complex(DP), allocatable :: coeff(:,:,:,:) 
  complex(DP), allocatable :: Gm(:,:,:)

  real(DP), allocatable :: logt(:), tee(:)
  integer, allocatable :: nt(:), run(:)
  complex(DP), allocatable :: s(:,:)
  integer :: ilogt, iminlogt, imaxlogt, lot, hit, lop, hip, ierr, lo
  character(128) :: subname = 'MainProgram (ltaem_main)'
  character(20) :: tmpfname

  real(DP), parameter :: MOST = 0.99_DP  ! 'most' of a log-cycle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! either specify here or ask for at prompt
  BGinfname = 'input'

  ! read in locations, options & data from input file
  call readInput(BGinfname)

  ! nudge times on 'edge' of logcycle down a tiny bit to increase accuracy
  where ((nint(BGt(1:BGnumt)) - BGt(1:BGnumt)) < SMALL)
     BGt(1:BGnumt) = BGt(1:BGnumt) - SMALL
  end where

  allocate(run(1:2*INVm+1))
  forall(i=0:2*INVm) run(i+1)=i

  ! independent of most choices
  call DistanceAngleCalcs()

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (BGcalc) then

     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     if (BGparticle) then   ! particle tracking
        
        ! need to do something better here about particles starting at t=0?
        do i=1, PARnum
           if (PARti(i) < 1.0D-5) then
              PARti(1:PARnum) = 1.0D-5
              write(*,'(A,I3,A)') '^^^^ start time for particle ',i,' reset to 1.0E-5 ^^^^'
           end if
        end do
        
        iminlogt = floor(  minval(log10(PARti(1:PARnum))))
        imaxlogt = ceiling(maxval(log10(PARtf(1:PARnum))))

        allocate(s(2*INVm+1,iminlogt:imaxlogt), nt(iminlogt:imaxlogt), &
             & tee(iminlogt:imaxlogt))

        nt(iminlogt:imaxlogt) = 1

        do ilogt = iminlogt, imaxlogt
           tee(ilogt) = min(10.0_DP**(ilogt + MOST),maxval(PARtf(1:PARnum)))*RTWO
           s(:,ilogt) = cmplx(INValpha - log(INVtol)/(RTWO*tee(ilogt)), PI*run/tee(ilogt), DP)
        end do

        ! to make it possible for particles / contours to share code...
        imaxlogt = imaxlogt + 1

        
        deallocate(run)
        
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     else   ! contour maps / hydrographs
        allocate(logt(1:BGnumt))

        lo = 1
        logt(1:BGnumt) = log10(BGt(1:BGnumt))
        iminlogt =   floor(minval(logt(1:BGnumt)))
        imaxlogt = ceiling(maxval(logt(1:BGnumt)))

        allocate(s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1))

        do ilogt = iminlogt, imaxlogt-1
           nt(ilogt) = count(logt >= real(ilogt,DP) .and. logt < real(ilogt+1,DP))
!!$           print *, 'ilogt:',ilogt
!!$           print *, BGt(lo:lo+nt(ilogt)-1)
           tee(ilogt) = maxval(BGt(lo:lo+nt(ilogt)-1))*RTWO
           s(:,ilogt) = cmplx(INValpha - log(INVtol)/(RTWO*tee(ilogt)), PI*run/tee(ilogt), DP)
           lo = lo + nt(ilogt)
        end do
        deallocate(logt,run)
     end if
     

     if(CInum > 0) then
        allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
             & Gm(1:4*CIn+2, 1:2*CIm, 1:CInum))
        
        coeff = CZERO
     
        ! calculate coefficients for each value of Laplace parameter
        ! ** common between particle tracking and contours/hydrographs **

        do ilogt = iminlogt,imaxlogt-1
           write(*,'(A,I3)') 'log t=',ilogt

           do j = 1,2*INVm+1
              if (nt(ilogt) > 0) then
                 write(*,'(I4,1X,2(A,ES10.3),A)') j, '(',real(s(j,ilogt)),',',aimag(s(j,ilogt)),')'
                 ! compute generalized inverse of A matrix
                 call circInverse_old(s(j,ilogt),Gm)  
!!$               ! use direct approach to solve once
!!$              call circInverse_matrix(s(j,ilogt),coeff(j,ilogt,0:4*CIn+1,1:CInum))  

                 if (j > 1) then  ! copy previous as initial guess
                    coeff(j,ilogt,0:4*CIn+1,1:CInum) = coeff(j-1,ilogt,0:4*CIn+1,1:CInum)
                 end if
                 call matchIter(Gm, coeff(j,ilogt,0:4*CIn+1,1:CInum), s(j,ilogt),CImatchTol)
              end if
           end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! save coefficient matrix in binary form to file (more accurate?)
        open(UNIT=77, FILE=BGCoeffFName, STATUS='REPLACE', ACTION='WRITE',&
             &  FORM='UNFORMATTED', IOSTAT=ierr)
        if (ierr /= 0) call fileError(BGCoeffFName,ierr,subname,0)
        write(77) iminlogt,imaxlogt,coeff,nt,s,tee,WLstorCoeff
        close(77)
        write(*,'(A)') 'matching finished'

     end if

  else ! do not re-calculate coefficients (this only makes sense if CInum > 0)

     if(CInum > 0) then
        open(UNIT=77, FILE=BGCoeffFName, STATUS='OLD', ACTION='READ',&
             &  FORM='UNFORMATTED', IOSTAT=ierr)
        if (ierr /= 0) then
           ! go back and recalc if no restart file
           write(*,'(A)') 'error opening restart file, recalculating...'
           BGcalc = .true.
           goto 111
        end if

        ! read two numbers needed to allocate arrays
        read(77) iminlogt,imaxlogt
        allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*CIn+1, CInum), &
             & s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1))

        ! read needed data calculated in previous run from file
        rewind(77)
        read(77)  iminlogt,imaxlogt,coeff,nt,s,tee,WLstorCoeff
        close(77)
        write(*,'(A)') 'matching results read from file'
     end if

  end if ! re-calculate coefficients

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(BGparticle) then ! integrate along particle path

     write(*,'(A)') 'compute solution for tracking particles'
     allocate(parnumdt(PARnum))

     parnumdt(:) = ceiling((PARtf(:) - PARti(:))/PARdt)

     allocate(Presult(0:maxval(parnumdt),5,PARnum))
     Presult = 0.0_DP

     ! cycle over particles, integrating each
     ! re-shape matrix s into a vector
     tmpfname = 'calc_part_    .debug'

     do part = 1, PARnum
        write(tmpfname(11:14),'(I4.4)') part
        open(unit=77,file=tmpfname,action='write',status='replace')

        select case (PARint(part))
        case (1)
           call rungekuttamerson(s,tee,coeff,lbound(coeff,dim=2),part)
        case (2)
                 call rungekutta(s,tee,coeff,lbound(coeff,dim=2),part)
        case (4)
                   call fwdEuler(s,tee,coeff,lbound(coeff,dim=2),part)
        case default
           write(*,'(A,I3,1X,I3)') 'invalid integration code', part, PARint(part)
        end select
        close(77)
     end do

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(BGcontour) then ! contour output (x,y locations are independent; e.g. outerproduct)

     write(*,'(A)') 'compute solution for plotting contours'
     allocate(head(BGnumx,BGnumy,BGnumt), headp(size(s)), &
          &   velx(BGnumx,BGnumy,BGnumt), velXp(size(s)), &
          &   vely(BGnumx,BGnumy,BGnumt), velYp(size(s)))

     ic = shape(coeff)  !! integer vector of matrix size in each dimension
     tnp = ic(1)*ic(2)  !! total number of Laplace parameters, over all times

     do j = 1,BGnumx
        write (*,'(A,ES14.6E1)') 'x: ',BGx(j)
        do i = 1,BGnumy

           !! compute f(p) for all values of p at this location 
           headp(1:tnp) = headCalc(reshape(s,[tnp]),BGx(j),BGy(i),reshape(coeff,[tnp,ic(3:4)]))
           velxp(1:tnp) = velxCalc(reshape(s,[tnp]),BGx(j),BGy(i),reshape(coeff,[tnp,ic(3:4)]))
           ! velx and vely share most of the calculations (velx _must_ be called first)
           velyp(1:tnp) = velyCalc(reshape(s,[tnp]),BGx(j),BGy(i))

           !! invert solutions one log-cycle of t at a time
           do ilogt = iminlogt,imaxlogt-1
              
              !! group of times in current log cycle
              lot = 1 + sum(nt(iminlogt:ilogt-1))
              hit = sum(nt(iminlogt:ilogt))
              
              !! group of Laplace parameters corresponding to this logcycle
              lop = (ilogt - iminlogt)*ic(1) + 1
              hip = lop + ic(1) - 1

              head(j,i,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),headp(lop:hip),INVm)
              velx(j,i,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),velxp(lop:hip),INVm)
              vely(j,i,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),velyp(lop:hip),INVm)
           end do
        end do
     end do
     
  else ! hydrograph output (x,y locations are in pairs; e.g. inner product)

     write(*,'(A)') 'compute solution for plotting hydrograph'

     allocate(head(BGnumx,1,BGnumt), headp(size(s)), &
          &   velx(BGnumx,1,BGnumt), velXp(size(s)), &
          &   vely(BGnumx,1,BGnumt), velYp(size(s)))
     
     do i = 1,BGnumx
        write(*,'(A,2(1X,ES14.7E1))') 'location:',BGx(i),BGy(i)

        if (allocated(coeff)) then
           ic = shape(coeff)
           tnp = ic(1)*ic(2)   
        else
           tnp = size(s)
           ic(1:4) = [0,0,0,0]
        end if
        

        headp(1:tnp) = headCalc(reshape(s,[tnp]),BGx(i),BGy(i),reshape(coeff,[tnp,ic(3:4)]))
        velxp(1:tnp) = velxCalc(reshape(s,[tnp]),BGx(i),BGy(i),reshape(coeff,[tnp,ic(3:4)]))
        velyp(1:tnp) = velyCalc(reshape(s,[tnp]),BGx(i),BGy(i))
        
        do ilogt = iminlogt,imaxlogt-1
           lot = 1 + sum(nt(iminlogt:ilogt-1))
           hit = sum(nt(iminlogt:ilogt))

           lop = (ilogt - iminlogt)*ic(1) + 1
           hip = lop + ic(1) - 1

!!$           print *, lot,hit,ilogt

           ! don't need second dimension of results matricies
           head(i,1,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),headp(lop:hip),INVm)
           velx(i,1,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),velxp(lop:hip),INVm)
           vely(i,1,lot:hit) = invlap(INValpha,INVtol,BGt(lot:hit),tee(ilogt),velyp(lop:hip),INVm)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file

!!$  
!!$  deallocate(coeff,s,nt)
!!$  call freeCalcMem()

  if (BGparticle) then 
     ! make sure correct output if doing particle tracking
     if (BGoutput /= 4 .and. BGoutput /= 5) then
        BGoutput = 4
     end if
  else
     ! reset back to gnuplot output if not particle tracking
     if (BGoutput == 4) BGoutput = 1
  end if

  ! call subroutine to write output to file
  call writeResults(head,velx,vely,BGx,BGy,BGt,BGoutput,BGoutfname,Presult)


end program ltaem_main

