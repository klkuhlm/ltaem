! $Id: ltaem_el_main.f90,v 1.10 2007/09/19 21:11:10 kris Exp kris $
program ltaem_el_main
  use constants, only     : DP, PI, RTWO, SMALL,CZERO
  use calc_shared_data, only : Presult
  use element_specs, only : BGinfname,BGcalc,BGparticle,BGnumt,BGt,&
       & BGx,BGy,Bgt,INVm,INVtol,INValpha,EIn,EInum,&
       & BGCoeffFName,BGcontour,BGnumx,BGnumy,BGoutput,BGoutfname,EIm,&
       & av,EIInclUp,EIms,EIf
  use file_ops, only : readEllipseinput, writeresults
  use error_handler, only : fileerror
  use inverse_Laplace_Transform, only : deHoog_invlap
  use matching_ellipse_only
  use calc_ellipse_routines
  use elliptical_geometry, only : ElDistanceAngleCalcs
  use shared_mathieu, only : A,B
  use mcn_matrix_method, only : mcn_eigenvalues

  implicit none
  integer :: i, j, tnp,k,el
  integer, dimension(4) :: ic
  integer, dimension(6) :: mc
  real(DP), allocatable :: head(:,:,:), velx(:,:,:), vely(:,:,:)
  complex(DP), allocatable :: headp(:),velp(:,:)
  complex(DP), allocatable :: coeff(:,:,:,:) 
  complex(DP), allocatable :: AA(:,:,:,:,:,:),BB(:,:,:,:,:,:),qq(:,:,:)
  complex(DP), allocatable :: Am(:,:,:)

  real(DP), allocatable :: logt(:), tee(:)
  integer, allocatable :: nt(:), run(:)
  complex(DP), allocatable :: s(:,:)
  integer :: ilogt, iminlogt, imaxlogt, lot, hit, lop, hip, ierr, lo
  character(128) :: subname = 'MainProgram (ltaem_main)'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! either specify here or ask for at prompt
  BGinfname = 'input.ellipse'

  ! read in locations, options & data from input file
  call readEllipseInput(BGinfname)

  ! nudge times on 'edge' of logcycle down a tiny bit to increase accuracy
  where ((nint(BGt(1:BGnumt)) - BGt(1:BGnumt)) < SMALL)
     BGt(1:BGnumt) = BGt(1:BGnumt) - SMALL
  end where

  allocate(run(1:2*INVm+1),stat=ierr)
  if (ierr /= 0) stop 'error allocating: RUN'
  run = (/ (i, i=0,2*INVm) /)

  ! independent of most choices
  call ElDistanceAngleCalcs()

111 continue ! come back here if error opening restart file

  ! calculate the values of p needed for inverse LT
  if (BGcalc) then

        allocate(logt(1:BGnumt),stat=ierr)
        if (ierr /= 0) stop 'error allocating: LOGT'

        lo = 1
        logt(1:BGnumt) = log10(BGt(1:BGnumt))
        iminlogt =   floor(minval(logt(1:BGnumt)))
        imaxlogt = ceiling(maxval(logt(1:BGnumt)))

        allocate(s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
             & tee(iminlogt:imaxlogt-1), stat=ierr)
        if (ierr /= 0) stop 'error allocating: S, NT, TEE'

        do ilogt = iminlogt, imaxlogt-1
           nt(ilogt) = count(logt >= real(ilogt,DP) .and. logt < real(ilogt+1,DP))
           print *, BGt(lo:lo+nt(ilogt)-1)
           tee(ilogt) = maxval(BGt(lo:lo+nt(ilogt)-1))*RTWO
           s(:,ilogt) = cmplx(INValpha - log(INVtol)/(RTWO*tee(ilogt)), &
                & PI*run/tee(ilogt), DP)
           lo = lo + nt(ilogt)
        end do
        deallocate(logt,run, stat=ierr)
        if (ierr /= 0) stop 'error deallocating: LOGT, RUN'

     allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*EIn+1, EInum),&
          & Am(1:2*EIm, 1:4*EIn+2, 1:EInum), stat=ierr)
     if (ierr /= 0) stop 'error allocating: COEFF, Am'
     
     coeff = CZERO

     !! compute Mathieu coefficients for each ellipse and value of p
     
     allocate(AA(EIms,0:EIms-1,2,2*EInum,2*INVm+1,iminlogt:imaxlogt-1),&
            & BB(EIms,0:EIms-1,2,2*EInum,2*INVm+1,iminlogt:imaxlogt-1),&
            & qq(2*EInum,2*INVm+1,iminlogt:imaxlogt-1))
     
     !! mathieu parameter (assuming all elements have different f,
     !! which is worst case scenario)
     do k=1,EInum
        !! inside q
        qq(k,:,:) = EIf(k)**2*s(:,:)/(4.0_DP*av(k))
        !! outside q 
        qq(k+EInum,:,:) = EIf(k)**2*s(:,:)/(4.0_DP*av(EIInclUp(k)))
     end do

     !! compute Mathieu coefficients for each value of Laplace parameter
     do ilogt = iminlogt,imaxlogt-1
        do j = 1,2*INVm+1
           do el = 1,2*EInum
              call mcn_eigenvalues(qq(el,j,ilogt),&
                   & AA(:,:,:,el,j,ilogt),BB(:,:,:,el,j,ilogt),1)
           end do
        end do
     end do

#ifdef DEBUG
     open(unit=78,file='debug_AA_BB.out')
     write(78,*) '# debug output, all Mathieu coefficients (long)'
     do k=1,EInum
        do ilogt = iminlogt,imaxlogt-1
           do i=1,2*INVm+1
              write(78,*) 'inside q: element',k,' alpha(k):',av(k), qq(k,i,ilogt)
              write(78,*) 'outside q: element',k,' parent:',EIInclUp(k),&
                   & ' alpha(k):',av(EIInclUp(k)), qq(k+EInum,i,ilogt)
              write(78,*) k,' A even; p=',s(i,ilogt)
              do j = 1,EIms
                 write(78,*) AA(j,:,1,k,i,ilogt)
              end do
              write(78,*) 'element:',k,' A odd; p=',s(i,ilogt)
              do j = 1,EIms
                 write(78,*) AA(j,:,2,k,i,ilogt)
              end do
              write(78,*) 'element:',k,' B even; p=',s(i,ilogt)              
              do j = 1,EIms
                 write(78,*) BB(j,:,1,k,i,ilogt)
              end do
              write(78,*) 'element:',k,' B odd; p=',s(i,ilogt)
              do j = 1,EIms
                 write(78,*) BB(j,:,2,k,i,ilogt)
              end do 
           end do
        end do
     end do
     close(78)
#endif

     !! allocate shared "local" version of Mathieu coefficients
     allocate(A(1:EIms,0:EIms-1,1:2),B(1:EIms,0:EIms-1,1:2))

#ifdef DEBUG
     open(unit=67,file='debug_Gm.out')
     write(67,*) '# debugging Gm generalized inverse matrix (long)'
     open(unit=55,file='debug_coeff.out')
     write(55,*) '# debugging LT-AEM coefficients'
#endif

     do ilogt = iminlogt,imaxlogt-1
        write (*,*) 'log t=',ilogt
        do j = 1,2*INVm+1
           if (nt(ilogt) > 0) then
              print *, j, s(j,ilogt)
              ! compute Am matrix to use in Least-Squares
              call EllipseInverse_old(s(j,ilogt),Am(1:2*EIm, 1:4*EIn+2, 1:EInum),&
                   & AA(:,:,:,1:EInum,j,ilogt),BB(:,:,:,1:EInum,j,ilogt),qq(1:EInum,j,ilogt))

#ifdef DEBUG
              write(67,*) 's=',s(j,ilogt)
              do i=1,EInum
                 write(67,*) 'Am(), element=',i
                 do k=1,2*EIm
                    write(67,*) Am(k,:,i)
                 end do
              end do
#endif
              
              call matchIter(Am( 1:2*EIm, 1:4*EIn+2, 1:EInum), &                   
                   & coeff(j,ilogt,0:4*EIn+1,1:EInum), s(j,ilogt),1.0D-11,&
                   & AA(:,:,:,1:EInum,j,ilogt),BB(:,:,:,1:EInum,j,ilogt),qq(1:EInum,j,ilogt))

#ifdef DEBUG
              write(55,*) 's=',s(j,ilogt)
              do i=1,EInum
                 write(55,*) 'coeff element=',i
                 write(55,*) 'a:', coeff(j,ilogt,0:EIn,i)
                 write(55,*) 'b:', coeff(j,ilogt,EIn+1:2*EIn,i)
                 write(55,*) 'c:', coeff(j,ilogt,2*EIn+1:3*EIn+1,i)
                 write(55,*) 'd:', coeff(j,ilogt,3*EIn+2:4*EIn+1,i)
              end do
#endif
              
           end if
        end do
     end do

#ifdef DEBUG
     close(67)
     close(55)
#endif

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! save coefficient matrix in binary form to file (more accurate?)
     open(UNIT=77, FILE=BGCoeffFName, STATUS='REPLACE', ACTION='WRITE',&
          &  FORM='UNFORMATTED', IOSTAT=ierr)
     if (ierr /= 0) call fileError(BGCoeffFName,ierr,subname,0)
     write(77) iminlogt,imaxlogt,coeff,nt,s,tee
     close(77)

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
     allocate(coeff(2*INVm+1, iminlogt:imaxlogt-1, 0:4*EIn+1, EInum), &
            & s(2*INVm+1,iminlogt:imaxlogt-1), nt(iminlogt:imaxlogt-1), &
            & tee(iminlogt:imaxlogt-1),stat=ierr)
     if (ierr /= 0) stop 'error allocating: COEFF, S, NT'

     ! read needed data calculated in previous run from file
     rewind(77)
     read(77)  iminlogt,imaxlogt,coeff,nt,s,tee
     close(77)
     print *, 'matching results read from file'

  end if ! re-calculate coefficients

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(BGcontour) then ! contour output (x,y locations are independent; e.g. outerproduct)

!!$     print *, 'debugging 1'

     allocate(head(BGnumx,BGnumy,BGnumt), headp(size(coeff,dim=1)), &
          &   velx(BGnumx,BGnumy,BGnumt), velp(size(coeff,dim=1),2), &
          &   vely(BGnumx,BGnumy,BGnumt), stat=ierr)
     if (ierr /= 0) stop 'error allocating: HEAD, VELX, VELY'

     head = 2.0_DP
     ic = shape(coeff)
     tnp = ic(1)*ic(2)  !! total number of Laplace parameters, over all times

     mc = shape(AA)

!!$     print *, 'debugging 2', mc

     do j = 1,BGnumx
        write (*,*) 'x:',BGx(j)
        do i = 1,BGnumy

           !! compute f(p) for all values of p at this location 
           headp(1:tnp) = headCalc(reshape(s,[tnp]),BGx(j),BGy(i),&
                & reshape(coeff,[tnp,ic(3:4)]),reshape(AA,[mc(1:4),tnp]),&
                & reshape(BB,[mc(1:4),tnp]),reshape(qq,[mc(4),tnp]))
           velp(1:tnp,1:2) = velCalc(reshape(s,[tnp]),BGx(j),BGy(i),&
                & reshape(coeff,[tnp,ic(3:4)]),reshape(AA,[mc(1:4),tnp]),&
                & reshape(BB,[mc(1:4),tnp]),reshape(qq,[mc(4),tnp]))
           
           !! invert solutions one log-cycle of t at a time
           do ilogt = iminlogt,imaxlogt-1
              
              !! group of times in current log cycle
              lot = 1 + sum(nt(iminlogt:ilogt-1))
              hit = sum(nt(iminlogt:ilogt))
              
              !! group of Laplace parameters corresponding to this logcycle
              lop = ilogt - iminlogt + 1
              hip = lop + ic(1) - 1

              head(j,i,lot:hit) = deHoog_invlap(INValpha,INVtol,BGt(lot:hit),&
                   & tee(ilogt),headp(lop:hip),INVm)
              velx(j,i,lot:hit) = deHoog_invlap(INValpha,INVtol,BGt(lot:hit),&
                   & tee(ilogt),velp(lop:hip,1),INVm)
              vely(j,i,lot:hit) = deHoog_invlap(INValpha,INVtol,BGt(lot:hit),&
                   & tee(ilogt),velp(lop:hip,2),INVm)
           end do
        end do
     end do
  else 

     ! only head makes sense for a hydrograph output
     ! hydrograph output (x,y locations are in pairs)

     allocate(head(BGnumx,1,BGnumt), stat=ierr)
     if (ierr /= 0) stop 'error allocating: HEAD'
     
     do i = 1,BGnumx
        print *, 'location:',BGx(i),BGy(i)
        
        headp(1:tnp) = headCalc(reshape(s,[tnp]),BGx(i),BGy(i),&
                & reshape(coeff,[tnp,ic(3:4)]),reshape(AA,[mc(1:4),tnp]),&
                & reshape(BB,[mc(1:4),tnp]),reshape(qq,[mc(4),tnp]))
        
        do ilogt = iminlogt,imaxlogt-1
           lot = 1 + sum(nt(iminlogt:ilogt-1))
           hit = sum(nt(iminlogt:ilogt))

           lop = ilogt - iminlogt + 1
           hip = lop + ic(1) - 1

           ! don't need second dimension of results matricies
           head(i,1,lot:hit) = deHoog_invlap(INValpha,INVtol,BGt(lot:hit),&
                & tee(ilogt),headp(lop:hip),INVm)
        end do
     end do
  end if

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! cleanup memory and write output to file

  deallocate(coeff,s,nt, stat=ierr)
  if (ierr /= 0) stop 'error deallocating: COEFF, S, NT'

  if (BGparticle) then 
     ! make sure correct output if doing particle tracking
     BGoutput = 4
  else
     ! reset back to gnuplot output if not particle tracking
     if (BGoutput == 4) BGoutput = 1
  end if

  ! call subroutine to write output to file
  call writeResults(head,velx,vely,BGx,BGy,BGt,BGoutput,BGoutfname,Presult)

end program ltaem_el_main

