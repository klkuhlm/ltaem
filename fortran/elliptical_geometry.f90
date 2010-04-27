 ! $Id: elliptical_geometry.f90,v 1.9 2007/10/01 00:01:01 kris Exp kris $
module elliptical_geometry
  implicit none

  private
  public :: ElDistanceAngleCalcs

contains

  ! ##################################################
  ! initializes / calculates geometry - all global variables
  subroutine ElDistanceAngleCalcs()
    use constants, only: DP, PI, TWOPI, EYE
    use element_specs, only : EIm,EInum,WLnum,EIeta,EIf,EIx,EIy,EItheta,EIWellBg,EIWellIn,&
         & WLq,WLr,WLx,WLy,EIInclIn,EIInclUp,EIInclBg,EIWellIn,EIWellBg,EIWellUp,EICalcIn
    use shared_matching_el_data, only :  EIPcm, EIXcm, EIYcm, EIXom, EIYom, EIRwm, EIPwm, EIEgm, EIPgm
    use file_ops, only : writeGeometry

!!$    character(50) :: fmt

    ! INTERNAL VARIABLES
    ! x/y distances from wells to points on circumference of elements
    real(DP), dimension(EIm,WLnum,EInum) :: EIXwm, EIYwm
!!$    ! x/y distances from center of other circles to circumference of elements
!!$    real(DP), dimension(EIm,EInum,EInum) :: EIXgm, EIYgm
    real(DP) :: rm 
    integer :: incl, well, other, M, ni, nw, i
    complex(DP), dimension(EIm,EInum) :: omega, z
    complex(DP), dimension(EInum) :: z0
    complex(DP), dimension(EIm) :: wtmp !!, ztmp

    M = EIm; ni = EInum; nw = WLnum;

    ! variables shared through ELEMENT_SPECS module -- initialized and filled here
    allocate(EIInclIn(0:ni,ni), EIInclUp(ni), EIInclBg(ni,ni), &
         & EIWellIn(0:ni,nw), EIWellBg(ni,nw), EIWellUp(nw), EICalcIn(0:ni), EIPcm(1:M))

    ! vector of eqi-spaced (in terms of psi) locations on perimeter of ellipse
    rM = real(M,DP)
    EIPcm = (/( -PI + TWOPI/rM * real(i,DP), i= 0,M-1)/)

    call EllipticalElementHierarchy()

    ! elliptical element related geometry
    allocate(EIXcm(M,ni), EIYcm(M,ni), EIXom(M,ni), &
         & EIYom(M,ni), EIRwm(M,nw,ni), EIPwm(M,nw,ni), &
         & EIEgm(M,ni,ni), EIPgm(M,ni,ni))

    !! translation to each local coordinate
    !! complex elliptical coordinates of all matching pts
    z0(1:ni) = cmplx(EIx,EIy,DP)
    forall (incl=1:ni)
       omega(1:M,incl) = cmplx(spread(EIeta(incl),dim=1,ncopies=M),EIPcm(1:M),DP)
       z(1:M,incl) = EIf(incl)*ccosh(omega(1:M,incl))*exp(EYE*EItheta(incl)) + z0(incl)
    end forall

    EIXom(1:M,1:ni) = real(z(:,:))
    EIYom(1:M,1:ni) = aimag(z(:,:))
    
    do incl = 1,ni
       ! x,y,r & theta from each well to the points on the circumference of each element
       forall (well = 1:nw, EIWellBg(incl,well) .or. EIWellIn(incl,well))
          EIXwm(1:M,well,incl) = EIXom(1:M,incl) - WLx(well)
          EIYwm(1:M,well,incl) = EIYom(1:M,incl) - WLy(well)
          EIRwm(1:M,well,incl) = abs(cmplx(EIXwm(1:M,well,incl),EIYwm(1:M,well,incl),DP))
          !! theta used to project well flux effects onto Cartesian
          EIPwm(1:M,well,incl) = atan2(EIYwm(1:M,well,incl),EIXwm(1:M,well,incl))
       end forall

       ! x,y,eta & psi in local coords of 'other (aka source)' ellipse of the points
       !  on the circumference of 'incl (aka target)' ellipse
       do other = 1,ni
          if( EIInclBg(other,incl) .or. EIInclIn(other,incl) .or. EIInclIn(incl,other)) then
!!$             EIXgm(1:M,other,incl) = EIXom(1:M,incl) - EIx(other)
!!$             EIYgm(1:M,other,incl) = EIYom(1:M,incl) - EIy(other)
!!$             ztmp(1:M) = cmplx(EIXgm(1:M,other,incl),EIYgm(1:M,other,incl),DP)
          
             wtmp(1:M) = acosh((z(1:M,incl) - z0(other))*exp(-EYE*EItheta(other))/EIf(other))
             EIEgm(1:M,other,incl) = real(wtmp(1:M))
             EIPgm(1:M,other,incl) = aimag(wtmp(1:M))
          end if
       end do
    end do

!!$    open(unit=22,file='debug.geom',action='write',status='replace')
!!$    fmt = '(A,   (F8.4,1X))'
!!$    write(fmt(4:6),'(I3.3)') M
!!$    write(22,fmt) 'EIPcm',EIpcm(1:M)
!!$    do incl = 1,ni
!!$       write(22,'(A,I2,A,F8.4,A,2(F8.4,1X),A,F8.4)') 'element:',incl, ' eta_0:',EIeta(incl), &
!!$            & ' z0:',z0(incl), ' angle:',EItheta(incl)
!!$       write(22,fmt) 'EIXom',EIXom(1:M,incl)
!!$       write(22,fmt) 'EIYom',EIYom(1:M,incl)
!!$       write(22,'(A)') 'one well'
!!$       write(22,fmt) 'EIRwm',EIRwm(1:M,1,incl)
!!$       write(22,fmt) 'EIPwm',EIPwm(1:M,1,incl)
!!$       if(incl == 1) then
!!$          other = 2
!!$       else
!!$          other = 1
!!$       end if       
!!$       write(22,'(A,I3)') 'other element', other
!!$       write(22,fmt) 'EIXgm',EIXgm(1:M,other,incl)
!!$       write(22,fmt) 'EIYgm',EIYgm(1:M,other,incl)
!!$       write(22,fmt) 'EIEgm',EIEgm(1:M,other,incl)
!!$       write(22,fmt) 'EIPgm',EIPgm(1:M,other,incl)
!!$    end do
!!$    
!!$    close(22)

    ! create listing of points on circumference of circles for plotting
    call writeGeometry(EIXom,EIYom,WLx,WLy,Wlr,WLq)    
  end subroutine ElDistanceAngleCalcs


  !##################################################
  !##################################################
  subroutine EllipticalElementHierarchy()
    use constants, only : DP, SMALL, RZERO
    use element_specs, only : EInum, WLnum, EIWellIn, EICalcIn, EIWellBg, EIInclIn, &
         & EIInclIn,EIInclUp,EIInclBg,EIWellUp

    integer :: nw, ni, line
    nw = WLnum; ni = EInum
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ! heirarchy of circular elements
    ! a temporary fix until this is calculated from inclusion specs

    open(unit=56, file='el_match.input', status='old', action='read')
    open(unit=57, file='echo_el_element_hierarchy', status='replace', action='write')

    do line = 0, ni
       read(56,*) EIInclIn(line,1:ni)
       write(57,*) EIInclIn(line,1:ni),'EIInclIn(',line,',1:',ni,')'
    end do

    read(56,*) EIInclUp(1:ni)
    write(57,*) EIInclUp(1:ni),'EIInclUp(1:',ni,')'

    do line = 1, ni
       read(56,*) EIInclBg(line,1:ni)
       write(57,*) EIInclBg(line,1:ni), 'EIInclBg(',line,',1:',ni,')'
    end do

    read(56,*) EICalcIn(1:ni)
    write (57,*) EICalcIn(1:ni), 'EICalcIn(0:',ni,')'
    EICalcIn(0) = .true.

    if (.not. EICalcIn(0)) then
       write (57,*) '*****************************************************************************'
       write (57,*) ''
       write (57,*) 'will not calculate background: all level-1 elements must be specified BC.'
       write (57,*) ''
       write (57,*) '*****************************************************************************'
       print *, '*****************************************************************************'
       print *, ''
       print *, 'will not calculate background: all level-1 elements must be specified BC.'
       print *, ''
       print *, '*****************************************************************************'
    end if

    do line = 0, ni
       read(56,*) EIWellIn(line,1:nw)
       write(57,*) EIWellIn(line,1:nw), 'EIWellIn(',line,',1:',nw,')'
    end do

    do line = 1, ni
       read(56,*)  EIWellBg(line,1:nw)
       write(57,*) EIWellBg(line,1:nw),  'EIWellBg(',line,',1:',nw,')'
    end do

    read(56,*) EIWellUp(1:nw)
    write(57,*) EIWellUp(1:nw), 'EIWellUp(1:',nw,')'

    close(56)
    close(57)

  end subroutine EllipticalElementHierarchy
  
  pure elemental function ccosh(z) result(f)
    use constants, only : DP

    complex(DP), intent(in) :: z
    complex(DP) :: f
    real(DP) :: x,y
    
    x = real(z)
    y = aimag(z)

    f = cmplx(cosh(x)*cos(y),sinh(x)*sin(y),DP)

  end function ccosh

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  elemental function acosh(z) result(f)
    use constants, only : DP, RONE

    complex(DP), intent(in) :: z
    complex(DP) :: f

    if(real(z) > 0.0) then
       f = log(z + sqrt(z**2 - RONE))
    else
       f = log(z - sqrt(z**2 - RONE))
    end if

  end function acosh
  
end module  elliptical_geometry

