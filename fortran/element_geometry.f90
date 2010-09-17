! this module computes the geometry related to the element distributions

module geometry
  implicit none

  private
  public :: DistanceAngleCalcs

contains

  ! ##################################################
  ! initializes / calculates geometry
  subroutine DistanceAngleCalcs(c,e,bg,dom,sol)
    use constants, only: DP, PI, EYE
    use type_definitions, only : domain, circle, ellipse, element, solution, matching
    use file_ops, only : writeGeometry, readElementHierarchy
    use utility, only : ccosh, cacosh

    type(domain), intent(inout) :: dom
    type(circle),  target, intent(inout), dimension(:) :: c
    type(ellipse), target, intent(inout), dimension(:) :: e
    type(element), target, intent(inout) :: bg
    type(solution), intent(in) :: sol
    type(matching), pointer :: other => null()

    integer :: i, j, ne, nc, ntot, par, M
    complex(DP), allocatable :: z(:)
#ifdef DEBUG
    integer :: k
#endif

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    allocate(dom%InclIn(0:ntot,ntot), dom%InclUp(ntot), dom%InclBg(ntot,ntot))
    bg%id = 0

    ! vector of eqi-spaced locations on perimeter of circle and ellipse
    ! each element can have a different number of matching locations
    ! starting at -PI (-x axis) continuing around the circle CCW.
    do i=1,nc
       allocate(c(i)%Pcm(c(i)%M))
       forall(j = 1:c(i)%M)
          c(i)%Pcm(j) = -PI + 2.0*PI/c(i)%M*real(j-1,DP)
       end forall
       
       c(i)%id = i ! global ID
    end do
    do i=1,ne
       allocate(e(i)%Pcm(e(i)%M))
       forall (j = 1:e(i)%M)
          e(i)%Pcm(j) = -PI + 2.0*PI/e(i)%M*real(j-1,DP)
       end forall
       
       e(i)%id = i+nc ! global ID
    end do

    ! setup pointers to parent elements
    bg%parent => null()  ! background has no parent
    do i=1,nc
       par = dom%InclUp(i) 
       if (par==0) then
          ! circle has background as parent
          c(i)%parent => bg
       elseif (par <= nc) then
          ! circle has another circle as parent
          c(i)%parent => c(par)%element
       elseif (par <= ntot) then
          ! circle has ellipse as parent
          c(i)%parent => e(par)%element
       else
          write(*,'(A,(1X,I0))') 'error in parent element index',par,i
          stop 200
       end if
    end do
    do i=1,ne
       par = dom%InclUp(nc+i) 
       if (par == 0) then
          ! ellipse has background as parent
          e(i)%parent => bg
       elseif (par <= nc) then
          ! ellipse has circle as parent
          e(i)%parent => c(par)%element
       elseif (par <= ntot) then 
          ! ellipse has another ellipse as parent
          e(i)%parent => e(par)%element
       else
          write(*,'(A,(1X,I0))') 'error in parent element index',par,i
          stop 201
       end if
    end do

#ifdef DEBUG
    open(unit=101,file='geom_self.debug',action='write',status='replace')
#endif

    ! circular element self-geometry
    do i = 1,nc
       M = c(i)%M
       allocate(c(i)%Zcm(M), c(i)%Zom(M), c(i)%G(ntot))

       ! x,y components from center of element to points on circumference
       if (M > 1) then
          c(i)%Zcm(1:M) = c(i)%r*exp(c(i)%Pcm(1:M)*EYE)
       else
          ! when only one matching point move to center of element
          c(i)%Zcm(1) = cmplx(0.0,0.0,DP)
       end if

       ! x,y from Cartesian origin to point on circumference of element
       c(i)%Zom(1:M) = c(i)%Zcm(:) + cmplx(c(i)%x,c(i)%y,DP)

#ifdef DEBUG
       write(101,*) '# elem',i
       do j=1,M
          write(101,*) 0.0,0.0,real(c(i)%Zom(j)),aimag(c(i)%Zom(j))
       end do
       write(101,'(/)')
#endif

    end do

    ! elliptical element self-geometry
    do i = 1,ne
       M = e(i)%M
       allocate(e(i)%Zcm(M), e(i)%Zom(M), e(i)%G(ntot), z(M))

       ! local elliptical coordinates (r is eta)
       z(1:M) = e(i)%f*ccosh(cmplx(e(i)%r,e(i)%Pcm(1:M),DP))
       ! z is local Cartesian coordinate, with +x parallel to semi-focal axis

       ! x,y components from center of element to points on circumference
       ! account for rotation of local elliptical coordinates
       if (M > 1) then
          e(i)%Zcm(1:M) = z(:)*exp(EYE*e(i)%theta)
       else
          ! when only one matching location move to center of line between foci
          e(i)%Zcm(1) = cmplx(0.0,0.0,DP)
       end if
       deallocate(z)

       ! x,y from Cartesian origin to point on circumference of element
       e(i)%Zom(1:M) = e(i)%Zcm(:) + cmplx(e(i)%x,e(i)%y,DP)

#ifdef DEBUG
       write(101,*) '# elem',i+nc
       do j=1,M
          write(101,*) 0.0,0.0,real(e(i)%Zom(j)),aimag(e(i)%Zom(j))
       end do
       write(101,'(/)')
#endif

    end do

#ifdef DEBUG    
    close(101)
    open(unit=202,file='geom_other.debug',action='write',status='replace')
#endif

    ! compute radial distances and angles to points on the circumferece of other elements
    ! from this element (cross-geometry), in terms of the current circle's or ellipse's
    ! coordinate system.
    do i = 1,nc
       ! this element a circle
       do j = 1,ntot
          if (i /= j) then
             if (j <= nc) then
                other => c(j)%matching    ! other element a circle             
             else
                other => e(j-nc)%matching ! other element an ellipse
             end if
             M = other%M    

             allocate(c(i)%G(j)%Zgm(M), c(i)%G(j)%Rgm(M), c(i)%G(j)%Pgm(M))

             c(i)%G(j)%Zgm(1:M) = other%Zom(1:M) - cmplx(c(i)%x,c(i)%y,DP)
             c(i)%G(j)%Rgm(1:M) = abs(c(i)%G(j)%Zgm(1:M)) ! r
             c(i)%G(j)%Pgm(1:M) = atan2(aimag(c(i)%G(j)%Zgm(1:M)), &
                                    & real(c(i)%G(j)%Zgm(1:M))) ! theta
             other => null()

#ifdef DEBUG
             write(202,*) '# src:',i,' tgt:',j
             do k=1,M
                write(202,*) c(i)%x,c(i)%y,real(c(i)%G(j)%Zgm(k)),aimag(c(i)%G(j)%Zgm(k))
             end do
             write(202,'(/)')
#endif
          end if
       end do
    end do
    do i = 1,ne
       ! this element an ellipse
       do j = 1,ntot
          if (i+nc /= j) then
             if (j <= nc) then
                other => c(j)%matching    ! other element a circle
             else
                other => e(j-nc)%matching ! other element an ellipse
             end if
             M = other%M    

             allocate(e(i)%G(j)%Zgm(M), e(i)%G(j)%Rgm(M), e(i)%G(j)%Pgm(M), z(M))

             e(i)%G(j)%Zgm(1:M) = other%Zom(1:M) - cmplx(e(i)%x,e(i)%y,DP)
             z(1:M) = cacosh( e(i)%G(j)%Zgm(1:M)*exp(-EYE*e(i)%theta)/e(i)%f )
             e(i)%G(j)%Rgm(1:M) =  real(z(1:M)) ! eta
             e(i)%G(j)%Pgm(1:M) = aimag(z(1:M)) ! psi
             
             deallocate(z)
             other => null()

#ifdef DEBUG
             write(202,*) '# src:',i+nc,' tgt:',j
             do k=1,M
                write(202,*) e(i)%x,e(i)%y,real(e(i)%G(j)%Zgm(k)),aimag(e(i)%G(j)%Zgm(k))
             end do
             write(202,'(/)')
          
             write(202,*) '# src:',i+nc,' tgt:',j, '*converted from elliptcial coords*'
             do k=1,M
                write(202,*) e(i)%x,e(i)%y, &
                     & real(exp(EYE*e(i)%theta)*e(i)%f*ccosh(cmplx(e(i)%G(j)%Rgm(k),e(i)%G(j)%Pgm(k),DP))),&
                     & aimag(exp(EYE*e(i)%theta)*e(i)%f*ccosh(cmplx(e(i)%G(j)%Rgm(k),e(i)%G(j)%Pgm(k),DP)))
             end do
             write(202,'(/)')
#endif
          end if
       end do
    end do
    
    call ReadElementHierarchy(dom,sol)

#ifdef DEBUG
    close(202)
#endif

    ! create listing of points on circumference of circles for plotting
    call writeGeometry(c,e,sol)    
  end subroutine DistanceAngleCalcs

  !##################################################
  subroutine ComputeElementHierarchy(dom,c,e)
    use constants, only : DP, EYE
    use type_definitions, only : domain, solution, circle, ellipse
    use utility, only : cacosh

    type(domain), intent(inout) :: dom
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    integer :: nc, ne, ntot

    real(DP), allocatable :: Rcg(:,:), Eeg(:,:), R(:)
    complex(DP), allocatable :: Z(:,:)
    logical, allocatable :: nondiag(:,:), upper(:,:)
    integer, allocatable :: iv(:), nest(:)
    integer :: i, j

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    if (ntot > 1) then

       allocate(nondiag(ntot,ntot), upper(ntot,ntot))
       ! logical mask for non-diagonal elements
       nondiag = .true. 
       forall(i=1:ntot) nondiag(i,i) = .false.

       ! logical mask for upper triangular (not including diag) elements
       upper = .false. 
       forall(i=1:ntot, j=1:ntot, j>i)  upper(i,j) = .true.

       dom%InclIn = .false. !! dom%InclIn(0:ntot,1:ntot) logical
       dom%InclUp = huge(1) !! dom%InclUp(1:ntot)        integer
       dom%InclBg = .false. !! dom%InclBg(1:ntot,1:ntot) logical

       ! ## step 1 ####################
       ! determine what circular element each circular + elliptical element falls inside 
       ! possibly multiple elements, if multiply nested.  Determine if elements intersect.

       ! check centers of elements (rows = circles, columns = all elements)
       if (nc > 0) then
          allocate(Rcg(nc,ntot))
          Rcg(1:nc,1:nc) =      abs(spread(cmplx(c%x,c%y,DP),1,nc) - &
                                  & spread(cmplx(c%x,c%y,DP),2,nc))
          Rcg(1:nc,nc+1:ntot) = abs(spread(cmplx(e%x,e%y,DP),1,nc) - &
                                  & spread(cmplx(c%x,c%y,DP),2,ne))
          
          ! nondiag handles zero distance-to-self case
          where (Rcg(1:nc,1:ntot) < spread(c(1:nc)%r,2,ntot) .and. nondiag)
             ! Is center of target element (2nd dim)
             ! inside radius of source element (1st dim)?
             dom%InclIn(1:nc,1:ntot) = .true.
          end where

          ! two cases where one element center is inside another
          !  1) elements intersect (2 pts) or touch (1 point) = BAD
          !  2) smaller element inside larger element         = OK

          ! check circle-on-circle intersection
          do i = 1, nc
             do j = 1, nc
                if (i /= j) then
                   ! all points on circumferencce must be either inside or outside
                   if (any(abs(c(i)%G(j)%Rgm(:) - c(i)%r) < epsilon(1.0)) .or. &
                        & (any(c(i)%G(j)%Rgm(:) < c(i)%r) .and. &
                        &  any(c(i)%G(j)%Rgm(:) > c(i)%r))) then
                      write(*,*) 'INTERSECTING CIRCLES: ',i,j
                      stop 400
                   end if
                end if
             end do
          end do

          ! check circle-on-ellipse interesction
          do i = 1, nc
             do j = 1, ne
                if (any(abs(c(i)%G(nc+j)%Rgm(:) - c(i)%r) < epsilon(1.0)) .or. &
                     & (any(c(i)%G(nc+j)%Rgm(:) < c(i)%r) .and. &
                     &  any(c(i)%G(nc+j)%Rgm(:) > c(i)%r))) then
                   write(*,*) 'INTERSECTING CIRCLE & ELLIPSE: ',i,j
                   stop 401
                end if
             end do
          end do
       end if
       
       ! ## step 2 ####################
       ! determine what elliptical element each circ + ellip element falls inside ...

       ! check centers of elements (rows = ellipses, columns = all elements)
       if (ne > 0) then
          allocate(Eeg(ne,ntot),Z(ne,ntot))
          Z(1:ne,1:nc) =      spread(cmplx(c%x,c%y,DP),1,ne) - &
                            & spread(cmplx(e%x,e%y,DP),2,nc)
          Eeg(1:ne,1:nc) =      real(cacosh(Z(1:ne,1:nc)* &
                                  & spread(exp(-EYE*e%theta)/e%f,2,nc)))

          Z(1:ne,nc+1:ntot) = spread(cmplx(e%x,e%y,DP),1,ne) - &
                            & spread(cmplx(e%x,e%y,DP),2,ne)
          Eeg(1:ne,nc+1:ntot) = real(cacosh(Z(1:ne,nc+1:ntot)* &
                                  & spread(exp(-EYE*e%theta)/e%f,2,ne)))
          deallocate(Z)

          where (Eeg(1:ne,1:ntot) < spread(e%r,2,ntot) .and. nondiag)
             dom%InclIn(nc+1:ntot,1:ntot) = .true.
          end where

          ! check ellipse-on-circle intersection
          do i = 1, ne
             do j = 1, nc
                if (any(abs(e(i)%G(j)%Rgm(:) - e(i)%r) < epsilon(1.0)) .or. &
                     & (any(e(i)%G(j)%Rgm(:) < e(i)%r) .and. &
                     &  any(e(i)%G(j)%Rgm(:) > e(i)%r))) then
                   write(*,*) 'INTERSECTING ELLIPSE & CIRCLE: ',i,j
                   stop 402
                end if
             end do
          end do

          ! check ellipse-on-ellipse interesction
          do i = 1, ne
             do j = 1, ne
                if (i /= j) then
                   if (any(abs(e(i)%G(nc+j)%Rgm(:) - e(i)%r) < epsilon(1.0)) .or. &
                        & (any(e(i)%G(nc+j)%Rgm(:) < e(i)%r) .and. &
                        &  any(e(i)%G(nc+j)%Rgm(:) > e(i)%r))) then
                      write(*,*) 'INTERSECTING ELLIPSES: ',i,j
                      stop 403
                   end if
                end if
             end do
          end do
       end if
              
       ! #1# remove bigger elements from inside 
       ! smaller elements, for the case they are concentric
       allocate(R(1:ntot), nest(1:ntot))
       R(1:nc) = c(1:nc)%r
       R(nc+1:ntot) = e(1:ne)%r  ! TODO : double-check that circ/ellip can be compared validly  

       forall (i=1:ntot, j=1:ntot, i/=j, &
            & dom%InclIn(i,j).and.dom%InclIn(j,i), R(i)>R(j))
          dom%InclIn(i,j) = .false.
       end forall

       ! number of True values in each column 
       nest(1:ntot) = count(dom%InclIn(1:ntot,1:ntot),dim=1)  

       ! #2# if an element is not inside any other elements
       !  (whole column is false) it must be in the background
       where (all(.not. dom%InclIn(1:ntot,:),dim=1))
          ! zero row indicating element is inside element zero
          dom%InclIn(0,1:ntot) = .true.
          dom%InclUp(1:ntot) = 0  
       end where

       ! #3# if an element is only listed as being inside one element
       ! (one element in column is true), it is simple to determine parent
       iv = [(i,i=1,ntot)]
       do i=1,ntot
          if (count(dom%InclIn(1:ntot,i)) == 1) then
             ! only one true in that column; use it to mask sum of an integer vector
             dom%InclUp(i) = sum(iv,mask=dom%InclIn(1:ntot,i))
          end if
       end do
       
       ! #4# handle wells listed as being inside multiple elements
       do i=1,ntot
          nest = 
          if ( > 1) then
             nest = 


    else
       ! special case of only one element, no matching -- it is in background
       dom%InclIn(0:1,1) = [.true.,.false.]
       dom%InclUp(1) = 0
       dom%InclBg(1,1) = .false.       
    end if
  end subroutine ComputeElementHierarchy
end module  geometry

