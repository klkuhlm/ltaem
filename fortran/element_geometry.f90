! this module computes the geometry related to the element distributions

module geometry
  implicit none

  private
  public :: DistanceAngleCalcs

contains

  ! ##################################################
  ! initializes / calculates geometry - all global variables
  subroutine DistanceAngleCalcs(c,e,bg,dom,sol)
    use constants, only: DP, PI, EYE
    use type_definitions, only : domain, circle, ellipse, element, solution, matching
    use file_ops, only : writeGeometry
    use utility, only : ccosh, cacosh

    type(domain), intent(inout) :: dom
    type(circle),  target, intent(inout), dimension(:) :: c
    type(ellipse), target, intent(inout), dimension(:) :: e
    type(element), target, intent(inout) :: bg
    type(solution), intent(in) :: sol
    type(matching), pointer :: other => null()

    integer :: i, j, k, ne, nc, ntot, par, M, ierr
    complex(DP), allocatable :: z(:)

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    allocate(dom%InclIn(0:ntot,ntot), dom%InclUp(ntot), dom%InclBg(ntot,ntot),stat=ierr)
    if (ierr /= 0) stop 'error allocating hierarchy matrices, dom%InclIn, dom%inclUp, dom%InclBg'

    ! vector of eqi-spaced locations on perimeter of circle and ellipse
    ! each element can have a different number of matching locations
    do i=1,nc
       allocate(c(i)%Pcm(c(i)%M),stat=ierr)
       if (ierr /= 0) stop 'error allocating angular matching location vector, c(i)%Pcm(1:M)'
       forall(j = 1:c(i)%M)
          c(i)%Pcm(j) = -PI + 2.0*PI/c(i)%M*real(j-1,DP)
       end forall
       
       c(i)%id = i ! global ID
    end do
    do i=1,ne
       allocate(e(i)%Pcm(e(i)%M),stat=ierr)
       if (ierr /= 0) stop 'error allocating angular matching location vector, e(i)%Pcm(1:M)'
       forall (j = 1:e(i)%M)
          e(i)%Pcm(j) = -PI + 2.0*PI/e(i)%M*real(j-1,DP)
       end forall
       
       e(i)%id = i+nc ! global ID
    end do

    call ElementHierarchy(dom,sol)

    ! setup pointers to parent elements
    bg%parent => null()  ! background has no parent
    do i=1,nc
       print *, 'geom circle parent ptrs: i=',i
       par = dom%InclUp(i) 
       if (par==0) then
          ! circle has background as parent
          c(i)%parent => bg
          print *, 'parent is bg',par
       elseif (par <= nc) then
          ! circle has another circle as parent
          c(i)%parent => c(par)%element
          print *, 'parent is circle',par
       elseif (par <= ntot) then
          ! circle has ellipse as parent
          c(i)%parent => e(par)%element
          print *, 'parent is ellipse',par
       else
          write(*,'(A,(1X,I0))') 'error in parent element index',par,i
          stop 200
       end if
    end do
    do i=1,ne
       print *, 'geom ellipse parent ptrs: i=',i
       par = dom%InclUp(nc+i) 
       if (par == 0) then
          ! ellipse has background as parent
          e(i)%parent => bg
          print *, 'parent is bg',par
       elseif (par <= nc) then
          ! ellipse has circle as parent
          e(i)%parent => c(par)%element
          print *, 'parent is circle',par
       elseif (par <= ntot) then 
          ! ellipse has another ellipse as parent
          e(i)%parent => e(par)%element
          print *, 'parent is ellipse',par
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
       allocate(c(i)%Zcm(M), c(i)%Zom(M), c(i)%G(ntot), stat=ierr)
       if (ierr /= 0) stop 'error allocating geometry matrices c(i)%{Zcm,Zom},c(i)%G(ntot)'

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
       allocate(e(i)%Zcm(M), e(i)%Zom(M), e(i)%G(ntot), z(M), stat=ierr)
       if (ierr /= 0) stop 'error allocating geometry matrices e(i)%{Zcm,Zom},e(i)%G(ntot), z'

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
       deallocate(z,stat=ierr)
       if (ierr /= 0) stop 'element_geometry.f90 error deallocating memory, z1'

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
          if (dom%InclBg(j,i) .or. dom%InclIn(j,i) .or. dom%InclIn(i,j)) then
             if (j <= nc) then
                other => c(j)%matching    ! other element a circle             
             else
                other => e(j-nc)%matching ! other element an ellipse
             end if
             M = other%M    

             allocate(c(i)%G(j)%Zgm(M), c(i)%G(j)%Rgm(M), c(i)%G(j)%Pgm(M),stat=ierr)
             if (ierr /= 0) stop 'error allocating geometry arrays c(i)%G(j)%{Zgm,Rgm,Pgm}'

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
          if (dom%InclBg(j,i+nc) .or. dom%InclIn(j,i+nc) .or. dom%InclIn(i+nc,j)) then
             if (j <= nc) then
                other => c(j)%matching    ! other element a circle
             else
                other => e(j-nc)%matching ! other element an ellipse
             end if
             M = other%M    

             allocate(e(i)%G(j)%Zgm(M),e(i)%G(j)%Rgm(M),e(i)%G(j)%Pgm(M),z(M),stat=ierr)
             if (ierr /= 0) stop 'error allocating geometry arrays e(i)%G(j)%{Zgm,Rgm,Pgm}'

             e(i)%G(j)%Zgm(1:M) = other%Zom(1:M) - cmplx(e(i)%x,e(i)%y,DP)
             z(1:M) = cacosh( e(i)%G(j)%Zgm(1:M)*exp(-EYE*e(i)%theta)/e(i)%f )
             e(i)%G(j)%Rgm(1:M) = real(z)  ! eta
             e(i)%G(j)%Pgm(1:M) = aimag(z) ! psi

             deallocate(z,stat=ierr)
             if (ierr /= 0) stop 'element_geometry.f90 error deallocating memory, z'
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

#ifdef DEBUG
    close(202)
#endif

    ! create listing of points on circumference of circles for plotting
    call writeGeometry(c,e,sol)    
  end subroutine DistanceAngleCalcs

  !##################################################
  subroutine ElementHierarchy(dom,sol)
    use constants, only : DP
    use type_definitions, only : domain, solution

    type(domain), intent(inout) :: dom
    type(solution), intent(in) :: sol
    integer :: nc,ne,ntot, line, ierr
    character(4) :: chint
    
    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    ! later I will write code to do this automatically
    open(UNIT=75, FILE=sol%elemHfName, STATUS='old', ACTION='read', IOSTAT=ierr)
    if (ierr /= 0) then
       write(*,'(A)') 'ElementHierarchy: ERROR opening file '//sol%elemHFName// &
            & ' for reading element hierarchy'
    end if
    open(UNIT=57, FILE=trim(sol%elemHfName)//'.echo', STATUS='replace', ACTION='write', IOSTAT=ierr)
    if (ierr /= 0) then
       write(*,'(A)') 'ElementHierarcy: ERROR opening file '//trim(sol%elemHfName)// &
            &'.echo for writing element hierarchy'
    else
       ! add a file variable to set Emacs to auto-revert mode
       write(57,'(A)') '-*-auto-revert-*-'  
    end if

    write(chint,'(I4.4)') ntot
    do line = 0,ntot
       read(75,*) dom%InclIn(line,1:ntot)
    end do
    do line = 0,ntot
       write(57,'('//chint//'(L1,1X),2(A,I0),A)')  &
            & dom%InclIn(line,1:ntot),'InclIn(',line,',1:',ntot,')'
    end do

    read(75,*) dom%InclUp(1:ntot)
    write(57,'('//chint//'(I0,1X),A,I0,A)') &
         & dom%InclUp(1:ntot),'InclUp(1:',ntot,')'

    do line = 1,ntot
       read(75,*) dom%InclBg(line,1:ntot)
    end do    
    do line = 1,ntot
       write(57,'('//chint//'(L1,1X),2(A,I0),A)') &
            & dom%InclBg(line,1:ntot), 'InclBg(',line,',1:',ntot,')'
    end do

    close(75)
    close(57)
  end subroutine ElementHierarchy

end module  geometry

