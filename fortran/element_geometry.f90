!
! Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

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
    use file_ops, only : writeGeometry
    use geomConv, only : c2xyA, e2xyA, xy2cA, xy2eA 

    type(domain), intent(inout) :: dom
    type(circle),  target, intent(inout), dimension(:) :: c
    type(ellipse), target, intent(inout), dimension(:) :: e
    type(element), target, intent(inout) :: bg
    type(solution), intent(in) :: sol
    type(matching), pointer :: other => null()

    integer :: i, j, ne, nc, ntot, par, M
    complex(DP), allocatable :: Zgm(:)

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

    ! circular element self-geometry
    do i = 1,nc
       M = c(i)%M
       allocate(c(i)%Zom(M), c(i)%G(ntot))
       
       if (M > 1) then 
          ! x,y from Cartesian origin to point on circumference of element
          c(i)%Zom(1:M) = c2xyA(cmplx(c(i)%r,c(i)%Pcm(1:M),DP),c(i)) 
       else
          ! when only one matching point move to center of element
          c(i)%Zom(1) = c(i)%z
       end if      
    end do

    ! elliptical element self-geometry
    do i = 1,ne
       M = e(i)%M
       allocate(e(i)%Zom(M), e(i)%G(ntot))

       if (M > 1) then
          ! x,y from Cartesian origin to point on circumference of element
          e(i)%Zom(1:M) = e2xyA(cmplx(e(i)%r,e(i)%Pcm(1:M),DP),e(i))
       else
          ! when only one matching location move to center of line between foci
          e(i)%Zom(1) = e(i)%z
       end if
    end do

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

             allocate(Zgm(M), c(i)%G(j)%Rgm(M), c(i)%G(j)%Pgm(M))

             Zgm(1:M) = xy2cA(other%Zom(1:M),c(i))
             c(i)%G(j)%Rgm(1:M) =   abs(Zgm(1:M)) ! r
             c(i)%G(j)%Pgm(1:M) = aimag(Zgm(1:M)) ! theta

             deallocate(Zgm)
             other => null()
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

             allocate(Zgm(M), e(i)%G(j)%Rgm(M), e(i)%G(j)%Pgm(M))

             Zgm(1:M) = xy2eA(other%Zom(1:M),e(i))
             e(i)%G(j)%Rgm(1:M) =  real(Zgm(1:M)) ! eta
             e(i)%G(j)%Pgm(1:M) = aimag(Zgm(1:M)) ! psi
             
             deallocate(Zgm)
             other => null()
          end if
       end do
    end do
    
    call ComputeElementHierarchy(dom,sol,c,e)

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

    ! create listing of points on circumference of circles for plotting
    call writeGeometry(c,e,sol)    
  end subroutine DistanceAngleCalcs

  !##################################################
  subroutine ComputeElementHierarchy(dom,sol,c,e)
    use constants, only : DP, EYE
    use type_definitions, only : domain, solution, circle, ellipse
    use utility, only : cacosh

    type(domain), intent(inout) :: dom
    type(solution), intent(in) :: sol
    type(circle),  dimension(:), intent(in) :: c
    type(ellipse), dimension(:), intent(in) :: e
    integer :: nc, ne, ntot, ierr

    real(DP), allocatable :: Rcg(:,:), Eeg(:,:), R(:)
    complex(DP), allocatable :: Z(:,:)
    ! single-byte integer representation of logical
    logical(1), allocatable :: nondiag(:,:)  
    integer, allocatable :: iv(:), nest(:), val(:)
    integer :: i, j, n, parent
    character(4) :: chint

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    if (ntot > 1) then

       allocate(nondiag(ntot,ntot))
       ! logical mask for non-diagonal elements
       nondiag = .true. 
       forall(i=1:ntot) nondiag(i,i) = .false.

       dom%InclIn = .false. !! dom%InclIn(0:ntot,1:ntot) logical
       dom%InclUp = -huge(1) !! dom%InclUp(1:ntot)       integer
       dom%InclBg = .false. !! dom%InclBg(1:ntot,1:ntot) logical

       ! ## step 1 ####################
       ! determine what circular element each circular + elliptical element falls inside 
       ! possibly multiple elements, if multiply nested.  Determine if elements intersect.

       ! check centers of elements (rows = circles, columns = all elements)
       if (nc > 0) then
          allocate(Rcg(nc,ntot))
          Rcg(1:nc,1:nc) =      abs(spread(c%z,1,nc) - spread(c%z,2,nc))
          Rcg(1:nc,nc+1:ntot) = abs(spread(e%z,1,nc) - spread(c%z,2,ne))
                    
          ! nondiag handles zero distance-to-self case
          where (Rcg(1:nc,1:ntot) < spread(c(1:nc)%r,2,ntot) .and. nondiag(1:nc,1:ntot))
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
                   if (any(abs(c(i)%G(j)%Rgm(:) - c(i)%r) < spacing(1.0)) .or. &
                        & (any(c(i)%G(j)%Rgm(:) < c(i)%r) .and. &
                        &  any(c(i)%G(j)%Rgm(:) > c(i)%r))) then
                      write(*,*) 'ERROR: INTERSECTING CIRCLES: ',i,j
                      stop 400
                   end if
                end if
             end do
          end do

          ! check circle-on-ellipse interesction
          do i = 1, nc
             do j = 1, ne
                if (any(abs(c(i)%G(nc+j)%Rgm(:) - c(i)%r) < spacing(1.0)) .or. &
                     & (any(c(i)%G(nc+j)%Rgm(:) < c(i)%r) .and. &
                     &  any(c(i)%G(nc+j)%Rgm(:) > c(i)%r))) then
                   write(*,*) 'ERROR: INTERSECTING CIRCLE & ELLIPSE: ',i,j
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
          Z(1:ne,1:nc) =      spread(c%z,1,ne) - spread(e%z,2,nc)
          Eeg(1:ne,1:nc) =      real(cacosh(Z(1:ne,1:nc)*spread(exp(-EYE*e%theta)/e%f,2,nc)))

          Z(1:ne,nc+1:ntot) = spread(e%z,1,ne) - spread(e%z,2,ne)
          Eeg(1:ne,nc+1:ntot) = real(cacosh(Z(1:ne,nc+1:ntot)*spread(exp(-EYE*e%theta)/e%f,2,ne)))
          deallocate(Z)

          where (Eeg(1:ne,1:ntot) < spread(e%r,2,ntot) .and. nondiag(1:ne,1:ntot))
             dom%InclIn(nc+1:ntot,1:ntot) = .true.
          end where

          ! check ellipse-on-circle intersection
          do i = 1, ne
             do j = 1, nc
                if (any(abs(e(i)%G(j)%Rgm(:) - e(i)%r) < spacing(1.0)) .or. &
                     & (any(e(i)%G(j)%Rgm(:) < e(i)%r) .and. &
                     &  any(e(i)%G(j)%Rgm(:) > e(i)%r))) then
                   write(*,*) 'ERROR: INTERSECTING ELLIPSE & CIRCLE: ',i,j
                   stop 402
                end if
             end do
          end do

          ! check ellipse-on-ellipse interesction
          do i = 1, ne
             do j = 1, ne
                if (i /= j) then
                   if (any(abs(e(i)%G(nc+j)%Rgm(:) - e(i)%r) < spacing(1.0)) .or. &
                        & (any(e(i)%G(nc+j)%Rgm(:) < e(i)%r) .and. &
                        &  any(e(i)%G(nc+j)%Rgm(:) > e(i)%r))) then
                      write(*,*) 'ERROR: INTERSECTING ELLIPSES: ',i,j
                      stop 403
                   end if
                end if
             end do
          end do
       end if
              
       ! ## step 3 ####################
       ! adjust InclIn, fill in InclUp and InclBg

       ! ## 3.1 remove bigger elements from inside 
       ! smaller elements, for the case they are concentric
       allocate(R(ntot), nest(ntot), iv(ntot))
       forall (i=1:ntot) iv(i) = i
       R(1:nc) = c(1:nc)%r
       R(nc+1:ntot) = e(1:ne)%r  ! TODO : double-check that circ/ellip can be compared validly  

       forall (i=1:ntot, j=1:ntot, i/=j .and. dom%InclIn(i,j).and.dom%InclIn(j,i) .and. R(i)<R(j))
          dom%InclIn(i,j) = .false.
       end forall

       deallocate(R)

       ! number of True values in each column 
       nest(1:ntot) = count(dom%InclIn(1:ntot,1:ntot),dim=1)  

       ! ## 3.2 if a column has no T entries, it is in background (level 1 of tree)
       where (nest(1:ntot) == 0)
          dom%InclIn(0,1:ntot) = .true. 
          dom%InclUp(1:ntot) = 0  
       end where

       ! ## 3.3 if a column has one T entry, it is inside a background element (level 2 of tree)
       if (any(nest == 1)) then
          do i = 1,ntot
             if (nest(i) == 1) then
                ! pack returns vector of one value, use sum to scalarize
                dom%InclUp(i) = sum(pack(iv, mask=dom%InclIn(1:ntot,i)))
             end if
          end do
       
          do n = 2,maxval(nest)
             do i = 1,ntot
                if (nest(i) == n) then
                   allocate(val(n))
                   val(1:n) = pack(iv, mask=dom%inclIn(1:ntot,i))

                   ! of these multiple potential parents, choose one
                   ! that has the highest number of T values in its column
                   ! set rest of the column to F

                   dom%InclIn(1:ntot,i) = .false.
                   ! use sum() to scalarize results of maxloc()
                   parent = val(sum(maxloc(nest(val(1:n))))) 
                   deallocate(val)

                   dom%InclIn(parent,i) = .true.
                   dom%InclUp(i) = parent
                end if
             end do
          end do
       end if
       
       deallocate(nest)
       allocate(nest(0:ntot))

       ! ## 3.4 elements in same row of InclIn are in the background together
       nest(0:ntot) = count(dom%InclIn(0:ntot,1:ntot),dim=2)

       ! set InclBg array up
       do i=0,ntot
          n = nest(i)
          if (n > 1) then
             allocate(val(n))
             val(1:n) = pack(iv, mask=dom%InclIn(i,1:ntot))
             
             do j=1,n
                dom%InclBg(val(1:n),val(j)) = .true.
             end do
             deallocate(val)
          end if
       end do
       
       ! an element isn't in its own background (by convention)
       forall (i=1:ntot)
          dom%InclBg(i,i) = .false.
       end forall

       deallocate(nest,iv)

    else
       ! special case of only one element, no matching -- it is in background
       dom%InclIn(0:1,1) = [.true.,.false.]
       dom%InclUp(1) = 0
       dom%InclBg(1,1) = .false.       
    end if

    ! echo results of heirarchy calcs to file
    open(unit=57, file=sol%elemhfname, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       ! non-fatal error
       write(*,'(A)') 'ElementHierarcy: ERROR opening file '//trim(sol%elemHfName)// &
            &'.echo for writing element hierarchy'
    else
       ! add a file variable to set Emacs to auto-revert mode
       write(57,'(A)') '-*-auto-revert-*-'  

       write(chint,'(I4.4)') ntot
       do i = 0,ntot
          write(57,'('//chint//'(L1,1X),2(A,I0),A)')  &
               & dom%InclIn(i,1:ntot),'InclIn(',i,',1:',ntot,')'
       end do

       write(57,'('//chint//'(I0,1X),A,I0,A)') &
            & dom%InclUp(1:ntot),'InclUp(1:',ntot,')'

       do i = 1,ntot
          write(57,'('//chint//'(L1,1X),2(A,I0),A)') &
               & dom%InclBg(i,1:ntot), 'InclBg(',i,',1:',ntot,')'
       end do
    end if
    close(57)

  end subroutine ComputeElementHierarchy
end module  geometry

