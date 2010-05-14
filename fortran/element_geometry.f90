! this module computes the geometry related to the element distributions

module geometry
  implicit none

  private
  public :: DistanceAngleCalcs

contains

  ! ##################################################
  ! initializes / calculates geometry - all global variables
  subroutine DistanceAngleCalcs(cm)
    use constants, only: DP, PI
    use element_specs, only :  circle, matching
    
    use file_ops, only : writeGeometry

    type(matching) :: cm

    ! INTERNAL VARIABLES
    ! x/y distances from wells to points on circumference of elements
    real(DP), dimension(CIm,WLnum,CInum) :: CIXwm, CIYwm
    ! x/y distances from center of other circles to circumference of elements
    real(DP), dimension(CIm,CInum,CInum) :: CIXgm, CIYgm
    real(DP) :: rm 
    integer :: incl, well, other, M, ni, nw, phi !, line

    M = CIm; ni = CInum; nw = WLnum;

    ! variables shared through ELEMENT_SPECS module -- initialized and filled here
    allocate(CIInclIn(0:ni,ni), CIInclUp(ni), CIInclBg(ni,ni), &
         & CIWellIn(0:ni,nw), CIWellBg(ni,nw), CIWellUp(nw), CIPcm(1:M))  !! CICalcIn(0:ni), 

    ! vector of eqi-spaced locations on perimeter of circle
    forall ()
       forall (phi = 1:M)
          CIPcm(phi) = -PI + 2.0*PI/M * real(phi-1,DP)
       end forall
    end forall
    
    call CircularElementHierarchy()

    ! circular inclusion related geometry
    allocate(CIXcm(M,ni), CIYcm(M,ni), CIXom(M,ni), &
         & CIYom(M,ni), CIRwm(M,nw,ni), CIPwm(M,nw,ni), &
         & CIRgm(M,ni,ni), CIPgm(M,ni,ni))

    forall (incl = 1:ni)
       ! x,y components from center of element to points on circumference
       CIXcm(1:M,incl) = cos(CIPcm(1:M))*CIr(incl)
       CIYcm(1:M,incl) = sin(CIPcm(1:M))*CIr(incl)

       ! now from Cartesian origin to point on circumference of element
       CIXom(1:M,incl) = CIXcm(1:M,incl) + CIx(incl)
       CIYom(1:M,incl) = CIYcm(1:M,incl) + CIy(incl)

       ! x,y,r & theta from each well to the points on the circumference of each element
       forall (well = 1:nw, CIWellBg(incl,well) .or. CIWellIn(incl,well))
          CIXwm(1:M,well,incl) = CIXom(1:M,incl) - WLx(well)
          CIYwm(1:M,well,incl) = CIYom(1:M,incl) - WLy(well)
          CIRwm(1:M,well,incl) = sqrt(CIXwm(1:M,well,incl)**2 + CIYwm(1:M,well,incl)**2)
          CIPwm(1:M,well,incl) = atan2(CIYwm(1:M,well,incl),CIXwm(1:M,well,incl))
       end forall

       ! x,y,r & theta from center of each element to the points on the circumference of BG and parent elements
       forall (other = 1:ni, CIInclBg(other,incl) .or. CIInclIn(other,incl) .or. CIInclIn(incl,other))
          CIXgm(1:M,other,incl) = CIXom(1:M,incl) - CIx(other)
          CIYgm(1:M,other,incl) = CIYom(1:M,incl) - CIy(other)
          CIRgm(1:M,other,incl) = sqrt(CIXgm(1:M,other,incl)**2 + CIYgm(1:M,other,incl)**2)
          CIPgm(1:M,other,incl) = atan2(CIYgm(1:M,other,incl),CIXgm(1:M,other,incl))
       end forall
    end forall

    ! create listing of points on circumference of circles for plotting
    call writeGeometry(CIXom,CIYom,WLx,WLy,Wlr,WLq)    
  end subroutine DistanceAngleCalcs


  !##################################################
  !##################################################
  subroutine CircularElementHierarchy()
    use constants, only : DP, SMALL, RZERO
    use element_specs, only : CInum, WLnum, CIWellIn,  CIWellBg, CIInclIn, &
         & WLx,WLy,CIx,CIy,CIr,CIInclIn,CIInclUp,CIInclBg,CIWellUp

!!$    character(40) :: fmt
    integer, dimension(CInum) :: CIlevel
    real(DP), dimension(CInum,WLnum) :: Rwg
    real(DP), dimension(CInum,CInum) :: Rgg
    logical, dimension(CInum,CInum) ::  upper, nondiag
    integer :: nw, ni, line, i, j, numw, numc, level
    integer, dimension(CInum) :: iwells, ichild
    logical, dimension(CInum) :: lwells

!!    open(unit=89,file='debug.geom',action='write',status='replace')

    !##################################################
    !  first run at determining hierarchy of geometry
    !##################################################

    nw = WLnum; ni = CInum
    CIWellIn = .false.; CIWellBg = .false.; CIInclIn = .false.; CIInclBg = .false.
    !! inner = .false.; outer = .false.
    CIlevel = huge(1)

    nondiag = .true. ! logical mask for non-diagonal elements
    forall(i = 1:ni) nondiag(i,i) = .false.

    upper = .false. ! logical mask for upper triangular (not including diag) elements
    forall(i = 1:ni, j = 1:ni, j > i)  upper(i,j) = .true.

    !########################################
    ! determine which wells are in which inclusions;
    ! multiply nested wells will appear multiply.

    ! radius from each well to center of each inclusion
    Rwg(1:ni,1:nw) = sqrt( &
         &(spread(CIx(:),dim=2,ncopies=nw) - spread(WLx(:),dim=1,ncopies=ni))**2 +&
         &(spread(CIy(:),dim=2,ncopies=nw) - spread(WLy(:),dim=1,ncopies=ni))**2 )

    ! wells that are inside circles
    where (Rwg(1:ni,1:nw) < spread(CIr(1:ni),dim=2,ncopies=nw)) 
       CIWellIn(1:ni,1:nw) = .true.
    end where

    ! if a well is not in any circle, it is in the background (element 0)
    where (.not. any(CIWellIn(1:ni,1:nw),dim=1))
       CIWellIn(0,1:nw) = .true.
    end where
    
    ! well within double-precision of circular element boundary
    if (any(abs(Rwg(:,:) - spread(CIr(:),dim=2,ncopies=nw)) < SMALL)) then
       stop 'GEOMETRY ERROR: well _on_ boundary of circular element'
    end if

    ! determine if any inclusions are nested; 
    ! multiply nested inclusions will appear twice
    ! Rgg => Rgg(incl,other) -- diagonal is meaningless: sqrt(0)
    Rgg(1:ni,1:ni) = sqrt( &
         &(spread(CIx(:),dim=2,ncopies=ni) - spread(CIx(:),dim=1,ncopies=ni))**2 +&
         &(spread(CIy(:),dim=2,ncopies=ni) - spread(CIy(:),dim=1,ncopies=ni))**2 )

    ! discr produces a symmetric matrix, but only need upper triangle
    if(any( discr(spread(CIx(1:ni),dim=2,ncopies=ni),&
         &  spread(CIy(1:ni),dim=2,ncopies=ni),&
         &  spread(CIr(1:ni),dim=2,ncopies=ni),&
         &  spread(CIx(1:ni),dim=1,ncopies=ni),&
         &  spread(CIy(1:ni),dim=1,ncopies=ni),&
         &  spread(CIr(1:ni),dim=1,ncopies=ni)) >= -SMALL .and. upper)) then
       stop 'GEOMETRY ERROR: intersecting circles'
    end if

    ! first dimension -> parent; second dimension -> child
    ! which circular elements are in which (guaranteed no intersecting elements now)
    where (Rgg(:,:) < spread(CIr(:),dim=2,ncopies=ni) .and. nondiag)
       ! this results in a symmetric matrix, 
       ! but the end product is definitely not symmetric
       CIInclIn(1:ni,:) = .true.
    end where

    do i=1,ni
       do j=1,ni
          if (CIInclIn(i,j) .and. CIr(j) > CIr(i)) then
             CIInclIn(i,j) = .false.
          end if
       end do
    end do
    
    ! determine which elements are in the background (level 0)
    ! background elements (columns) have all false entries
    where (.not. any(CIInclIn(1:ni,:),dim=1))
       CIlevel(1:ni) = 0
       CIInclUp(1:ni) = 0
       CIInclIn(0,1:ni) = .true.
    end where
   
    !##################################################
    !  remove nested duplicates
    !##################################################

    ! take care of inclusions listed as being interior
    ! to more than one inclusions (an inclusion nested
    ! inside an inclusion, which is nested again, etc.)
    ! only want the "innermost" pair to remain

    LEV: do level = 0, ni-1 ! at most ni levels deep
       do i = 1, ni
          if (CIlevel(i) == level) then
             call clear_grandkids(CIInclIn(:,:),i)
             forall(j = 1:ni, CIInclIn(i,j))
                CIlevel(j) = level + 1
                CIInclUp(j) = i
             end forall
          end if
       end do
       if(all(CIlevel <= level)) exit LEV
    end do LEV


    !##################################################
    !  determine neighbors (background)
    !##################################################

    forall(i = 1:ni, j = 1:ni, CIInclUp(j) == CIInclUp(i))
       CIInclBg(i,j) = .true.
    end forall
    forall(i = 1:ni) CIInclBg(i,i) = .false.

    !##################################################
    !  remove nested and determine parent inclusion for nested wells
    !##################################################
    CIWellUp = 0
    do i = 1, nw
       lwells = CIWellIn(1:ni,i) !! don't include 0 (background)
       ! does this well appear in at least one inclusion?
       if(count(lwells) > 0) then
          call log2intl(lwells, iwells, numw) ! list of circles
          if(numw == 1) then
             CIWellUp(i) = iwells(1)
          else
             do j = 1, numw
                call log2intl(CIInclIn(iwells(j),:),ichild,numc)
                ! if numc == 0, then there are no circular elements
                ! inside the first one, so it is the "innermost" one

                ! looking for inclusion which does not have another inclusion
                ! from the list as a child
                if (numc > 0) then
                   ! this element has children
                   CIWellIn(iwells(j),i) = .false.
                else
                   ! this element has no children (it must
                   !  be the most immedeate parent of the well, then)
                   CIWellUp(i) = iwells(j)
                end if
             end do
          end if
       end if
    end do

!!$    close(89)


    !## if a well and inclusion have the same parent, the well
    !   must be in the background of the inclusion

    forall(i = 1:ni, j = 1:nw, CIInclUp(i) == CIWellUp(j))
       CIWellBg(i,j) = .true.
    end forall

    open(unit=57, file='echo_element_hierarchy',status='replace',action='write')
    
    do line = 0, ni
       write(57,*) CIInclIn(line,1:ni),'CIInclIn(',line,',1:',ni,')'
    end do

    write(57,*) CIInclUp(1:ni),'CIInclUp(1:',ni,')'

    do line = 1, ni
       write(57,*) CIInclBg(line,1:ni), 'CIInclBg(',line,',1:',ni,')'
    end do

    do line = 0, ni
       write(57,*) CIWellIn(line,1:nw), 'CIWellIn(',line,',1:',nw,')'
    end do

    do line = 1, ni
       write(57,*) CIWellBg(line,1:nw),  'CIWellBg(',line,',1:',nw,')'
    end do

    write(57,*) CIWellUp(1:nw), 'CIWellUp(1:',nw,')'
    close(57)

    !*******************************************************
    !*******************************************************
    !*******************************************************


  contains

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! this is part of the discriminant of the quadratic equation
    ! from solving for the intersection of two circles
    ! if discr > 0 two real roots (circle intersect)
    ! if discr < 0 complex conjugate roots (no intersection)
    ! if discr == 0 one double real root (circles touch at 1 pt)
    elemental function discr(a,b,r,c,d,w) result(z)
      use constants, only : DP
      ! (a,b) are (x,y) of center and r is radius for circle 1
      ! (c,d) are (x,y) of center and w is radius for circle 2 
      real(DP), intent(in) :: a,b,r, c,d,w
      real(DP) :: z

      ! actual discriminat would be divided by (b-d)**2, but
      ! that goes to zero when circles have same y-coordinate
      ! it doesn't change the classification, so it is left out

      z = -((a**2 - r**2 - 2*a*c + c**2 + (b-d)**2)**2 - &
           &  2*(r**2 + (a-c)**2 + (b-d)**2)*w**2 + w**4)

    end function discr

    !##################################################
    ! -- conversion routine --
    ! return an integer vector of the true entries
    ! in the logical vector, padded with trailing zeros
    pure function log2int(l) result(v)
      logical, intent(in), dimension(:) :: l
      integer, dimension(size(l)) :: v
      integer :: i,j,n

      n = size(l)
      j = 0; v = 0;

      do i = 1,n
         if(l(i)) then
            j = j + 1
            v(j) = i
         end if
      end do

    end function log2int

    !##################################################
    ! -- conversion routine --
    ! return an integer vector of the true entries
    ! in the logical vector, padded with trailing zeros
    pure subroutine log2intl(l,v,s)
      logical, intent(in), dimension(:) :: l
      integer, intent(out), dimension(size(l)) :: v
      integer, intent(out) :: s
      integer :: i

      s = 0; v = 0;

      do i = 1,size(l)
         if(l(i)) then
            s = s + 1
            v(s) = i
         end if
      end do

    end subroutine log2intl

    !##################################################
    ! clears doubly-listed "grandkids" of specified parent
    ! element only using CIInclIn array. Grandkids are elements
    ! which are inside a larger element, but they have some
    ! intermediate element between them.
    pure subroutine clear_grandkids(IN,par)

      ! EXTERNAL VARIABLES
      ! CIInclIn which is modified
      logical, intent(inout), dimension(:,:) :: IN
      ! element which is having its grandkids set to false
      integer, intent(in) :: par

      logical, dimension(size(IN,dim=2)) :: lkid,union
      integer, dimension(size(IN,dim=2)) :: ikid
      integer :: i

      union = .false.

      lkid = IN(par,:)
      ikid = log2int(lkid)

      do i = 1, count(lkid)
         union = union .and. IN(ikid(i),:)
      end do

      where (IN(par,:) .and. union)  IN(par,:) = .false.

    end subroutine clear_grandkids
  end subroutine CircularElementHierarchy
end module  circular_geometry

