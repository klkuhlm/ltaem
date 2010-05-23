! this module computes the geometry related to the element distributions

module geometry
  implicit none

  private
  public :: DistanceAngleCalcs

contains

  ! ##################################################
  ! initializes / calculates geometry - all global variables
  subroutine DistanceAngleCalcs(c,e,bg,dom,sol)
    use constants, only: DP, PI
    use type_definitions
    use file_ops, only : writeGeometry

    type(domain), intent(inout) :: dom
    type(circle), intent(in), dimension(:) :: c
    type(ellipse), intent(in), dimension(:) :: e
    type(element), intent(in) :: bg
    type(solution), intent(in) :: sol

    ! INTERNAL VARIABLES
    ! x/y distances from wells to points on circumference of elements
    real(DP), dimension(CIm,WLnum,CInum) :: CIXwm, CIYwm
    ! x/y distances from center of other circles to circumference of elements
    real(DP), dimension(CIm,CInum,CInum) :: CIXgm, CIYgm
    real(DP) :: rm 
    integer :: i, ne, nc, ntot, par
    integer :: incl, well, other, M, ni, nw, phi 

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    allocate(dom%InclIn(0:ntot,ntot), dom%InclUp(ntot), dom%InclBg(ntot,ntot))

    ! vector of eqi-spaced locations on perimeter of circle and ellipse
    do i=1,nc
       allocate(c(i)%Pcm(c(i)%M))
       forall (phi = 1:c(i)%M)
          c(i)%Pcm(phi) = -PI + 2.0*PI/c(i)%M*real(phi-1,DP)
       end forall
    end do
    do i=1,ne
       allocate(e(i)%Pcm(e(i)%M))
       forall (phi = 1:e(i)%M)
          e(i)%Pcm(phi) = -PI + 2.0*PI/e(i)%M*real(phi-1,DP)
       end forall
    end do

    call ElementHierarchy(dom,sol)

    ! setup pointers to parent elements
    bg%parent => null()  ! background has no parent
    do i=1,nc
       par = dom%InclUp(i) 
       select case(par)
       case(0)
          ! circle has background as parent
          c(i)%parent => bg
       case(1:nc)
          ! circle has another circle as parent
          c(i)%parent => c(par)%element
       case(nc+1:ntot)   
          ! circle has ellipse as parent
          c(i)%parent => e(par)%element
       case default
          write(*,'(A,(1X,I0))') 'error in parent element index',par,i
          stop 200
       end select
    end do
    do i=1,ne
       par = dom%InclUp(nc+i) 
       select case(par)
       case(0)
          ! ellipse has background as parent
          e(i)%parent => bg
       case(1:nc)
          ! ellipse has circle as parent
          e(i)%parent => c(par)%element
       case(nc+1:ntot) 
          ! ellipse has another ellipse as parent
          e(i)%parent => e(par)%element
       case default
          write(*,'(A,(1X,I0))') 'error in parent element index',par,i
          stop 201
       end select
    end do

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
  subroutine ElementHierarchy(dom,sol)
    use constants, only : DP
    use type_definitions, only : domain, solution

    type(domain), intent(inout) :: dom
    type(solution), intent(in) :: sol
    integer :: nc,ne,ntot, line

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    ! later I will write code to do this automatically

    open(unit=75, file=sol%elemHfName, status='old', action='read')
    open(unit=57, file=trim(sol%elemHfName)//'.echo',status='replace',action='write')

    do line = 0,ntot
       read(75,*) dom%InclIn(line,1:ntot)
    end do
    do line = 0,ntot
       write(57,*) dom%InclIn(line,1:ntot),'InclIn(',line,',1:',ni,')'
    end do

    read(75,*) dom%InclUp(1:ntot)
    write(57,*) InclUp(1:ni),'InclUp(1:',ni,')'

    do line = 1,ntot
       read(75,*) dom%InclBg(line,1:ntot)
    end do    
    do line = 1,ntot
       write(57,*) dom%InclBg(line,1:ni), 'InclBg(',line,',1:',ni,')'
    end do

    close(75)
    close(57)
  end subroutine ElementHierarchy

end module  geometry

