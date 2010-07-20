! this module computes the geometry related to the element distributions

module unsat_geometry
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

    integer :: i, j, ne, nc, ntot, par, M, ierr
    complex(DP), allocatable :: z(:)
#ifdef DEBUG
    integer :: k
#endif

    nc = dom%num(1)
    ne = dom%num(2)
    ntot = sum(dom%num)

    allocate(dom%InclIn(0:ntot,ntot), dom%InclUp(ntot), dom%InclBg(ntot,ntot),stat=ierr)
    if (ierr /= 0) stop 'error allocating hierarchy matrices, dom%InclIn, dom%inclUp, dom%InclBg'

    bg%id = 0

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


    ! setup pointers to parent elements
    bg%parent => null()  ! background has no parent
    do i=1,nc
       ! every circle has background as parent
       c(i)%parent => bg
    end do

    do i=1,ne
       ! ellipse has background as parent
       e(i)%parent => bg
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
          if (j /= i) then
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
          if (j /= i) then
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
             e(i)%G(j)%Rgm(1:M) =  real(z(1:M)) ! eta
             e(i)%G(j)%Pgm(1:M) = aimag(z(1:M)) ! psi

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

end module  unsat_geometry

