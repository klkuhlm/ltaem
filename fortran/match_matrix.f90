module matching
  implicit none

  private 
  public :: CircInverse_matrix,freeMatchMem

contains

  !##################################################
  subroutine Inverse_matrix(p,coeff)
    use constants, only: DP,PI,TWOPI,CZERO,CONE,CTWO,RONE,RZERO
    use element_specs, only : CIn,CIm,CInum,av,kv,CImatch,CIInclUp, &
         & CIibnd,CIr,CICalcIn,CIinclIn,CIInclBg
    use bessel_functions_deriv, only : besselk_val, besselk_val_and_deriv, besseli_val, besseli_val_and_deriv
    use bessel_functions
    use shared_matching_data, only : CIPcm,CIRgm,CIPgm,row,col

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
    ! scalar Laplace parameter
    complex(DP), intent(in) :: p

    ! generalized Fourier coefficients - the results
    complex(DP), dimension(0:4*CIn+1,CInum), intent(out) :: coeff

    ! INTERNAL VARIABLES
    ! orthogonal eigenvectors from separation of variables
    complex(DP), allocatable :: Amain(:,:)
    complex(DP), dimension(1:2*CIm, 1:4*CIn+2, CInum, CInum) :: Am

    ! Jacobians for computing flux effects between elements (intermediate)
    ! Jrr is deriv wrt radius 1 -> radius 2
    ! Jtr is deriv wrt theta 1 -> radius 2
    real(DP), dimension(1:CIm) :: Jrr, Jtr
    integer, save :: counter = 0

    ! vectors of Bessel functions to minimize calls to Bessel fcn subroutine
    complex(DP), dimension(CIm,0:CIn,CInum) :: besskRcm, bessiRcm
    complex(DP), dimension(CIm,-1:CIn+1,CInum,CInum) :: besskRgm, bessiRgm
    
    ! vectors of sin/cos to save calls in Jacobian calcs
    real(DP), dimension(CIm) :: cosPM, sinPM, cosGm, sinGm

    ! ratios of Bessel functions, for simplifying code
    complex(DP), dimension(1:CIm,0:CIn) :: BessRat1,BessRat2
    
    ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
    ! vectors of Bessel functions to minimize calls to Bessel fcn subroutine
    complex(DP), dimension(0:CIn) :: bessk, bessi, Dbessk, Dbessi
    ! cos/sin(theta) -- where theta is spaced evenly along circumference
    real(DP), dimension(1:CIm,0:CIn) :: cosnCIPcm, sinnCIPcm

    ! vector of trig functions to minimize calls to built-in functions
    complex(DP), dimension(CIm,0:CIn) :: cosnPgm, sinnPgm

    ! counters
    integer :: ni, e,t,s,par, k, M, N,i, colmax
    real(DP), dimension(0:CIn) :: rk

    !! fill the row/col vectors with the needed indicies, based on 
    !! the type of each element

    M = CIm; N = CIn; ni = CInum;
    q = sqrt(p/av)
    counter = counter + 1

    !! paste locally indexed sub-matricies into big global matrix
    if(.not. allocated(row)) then
       allocate(row(2,ni),col(2,ni))
       row(1,1) = 1; col(1,1) = 1
       if(CIibnd(1) == 0) then
          row(2,1) = 2*M
          col(2,1) = 4*N+2
       else
          row(2,1) = M
          col(2,1) = 2*N+1
       end if

       do i = 2,ni
          row(1,i) = row(2,i-1) + 1
          col(1,i) = col(2,i-1) + 1
          if(CIibnd(i) == 0) then
             row(2,i) = 2*M-1 + row(1,i)
             col(2,i) = 4*N+1 + col(1,i)
          else
             row(2,i) = M-1 + row(1,i)
             col(2,i) = 2*N + col(1,i)
          end if
       end do
    end if

    allocate(Amain(1:row(2,ni),1:col(2,ni)))
    Amain = CZERO

    sinnCIPcm(:,0) = RZERO;
    Am = CZERO

    rk(0:N) = real((/ (k, k=0,N) /),DP)

    !******************** Am matrix setup ********************
    ! calcualte coefficient matrix (Am) which depends on the geometry 
    ! of the problem and p

    ! diagonal sub-matricies (effect of an element at its own boundary)
    ELEM: do e = 1,ni
       if (CImatch(e)) then
          par = CIInclUp(e)  ! parent

          select case (CIibnd(e))  
          case (0)
             !^^^^^^^^^^^^^^^^^^^^ head & flux matching ^^^^^^^^^^^^^^^^^^^^
             ! pre-compute bessels
             call besselk_val_and_deriv(CIr(e)*q(par),N+1,bessk(0:N),Dbessk(0:N))
             Dbessk(0:N) = Dbessk*q(par)
             call besseli_val_and_deriv(CIr(e)*q(e),N+1,bessi(0:N),Dbessi(0:N))
             Dbessi(0:N) = Dbessi*q(e)

             cosnCIPcm(1:M,0:N) = cos(outer_prod(CIPcm(1:M),rk(0:N)))
             sinnCIPcm(1:M,1:N) = sin(outer_prod(CIPcm(1:M),rk(1:N)))

             ! head (PHI/k) matching
             Am(1:M, 1:N+1, e,e)       =  cosnCIPcm(1:M,0:N)/kv(par) ! a_n
             Am(1:M, N+2:2*N+1, e,e)   =  sinnCIPcm(1:M,1:N)/kv(par) ! b_n
             Am(1:M, 2*N+2:3*N+2, e,e) = -cosnCIPcm(1:M,0:N)/kv(e) ! c_n
             Am(1:M, 3*N+3:4*N+2, e,e) = -sinnCIPcm(1:M,1:N)/kv(e) ! d_n

             ! flux (dPHI/dr) matching
             ! a_n
             Am(M+1:2*M, 1:N+1, e,e) = &
                  &  spread(Dbessk(0:N)/bessk(0:N),1,M)*cosnCIPcm(1:M,0:N)
             ! b_n
             Am(M+1:2*M, N+2:2*N+1, e,e) = &
                  &  spread(Dbessk(1:N)/bessk(1:N),1,M)*sinnCIPcm(1:M,1:N)
             ! c_n
             Am(M+1:2*M, 2*N+2:3*N+2, e,e) = &
                  & -spread(Dbessi(0:N)/bessi(0:N),1,M)*cosnCIPcm(1:M,0:N)
             ! d_n
             Am(M+1:2*M, 3*N+3:4*N+2, e,e) = &
                  & -spread(Dbessi(1:N)/bessi(1:N),1,M)*sinnCIPcm(1:M,1:N)

             Amain(row(1,e):row(2,e),col(1,e):col(2,e)) = Am(1:2*M,1:4*N+2,e,e)

          case (-1)
             !^^^^^^^^^^^^^^^^^ specified head ^^^^^^^^^^^^^^^^^^^^
             ! head (PHI/k) matching
             Am(1:M, 1:N+1, e,e) =     cos(outer_prod(CIPcm(1:M),rk(0:N)))/kv(par) ! a_n
             Am(1:M, N+2:2*N+1, e,e) = sin(outer_prod(CIPcm(1:M),rk(1:N)))/kv(par) ! b_n

             !>>>>> calculate inside the specified head inclusion <<<<<
             if (CIcalcin(e)) then

                Am(1:M, 2*N+2:3*N+2, e,e) = -cos(outer_prod(CIPcm(1:M),rk(0:N)))/kv(e) ! c_n
                Am(1:M, 3*N+3:4*N+2, e,e) = -sin(outer_prod(CIPcm(1:M),rk(1:N)))/kv(e) ! d_n

                colmax = 4*N+2
             else
                colmax = 2*N+1
             end if

             Amain(row(1,e):row(2,e),col(1,e):col(2,e)) = Am(1:M,1:colmax,e,e)

          case (+1)

             !^^^^^^^^^^^^^^^^^^^^ specified flux ^^^^^^^^^^^^^^^^^^^^
             call besselk_val_and_deriv(CIr(e)*q(par),N+1,bessk(0:N),Dbessk(0:N))
             Dbessk(0:N) = Dbessk*q(par)

             Am(M+1:2*M, 1:N+1, e,e) = spread(Dbessk(0:N)/bessk(0:N),1,M)* &
                  & cos(outer_prod(CIPcm(1:M),rk(0:N))) ! a_n

             Am(M+1:2*M, N+2:2*N+1, e,e) = spread(Dbessk(1:N)/bessk(1:N),1,M)*&
                  & sin(outer_prod(CIPcm(1:M),rk(1:N))) ! b_n

             !>>>>> calculate inside the specified flux inclusion <<<<<
             if (CIcalcin(e)) then

                call besseli_val_and_deriv(CIr(e)*q(e),N+1,bessi(0:N),Dbessi(0:N))
                Dbessi(0:N) = Dbessi*q(e)

                Am(M+1:2*M, 2*N+2:3*N+2, e,e) = -spread(Dbessi(0:N)/bessi(0:N),1,M)*&
                     & cos(outer_prod(CIPcm(1:M),rk(0:N))) ! c_n
                Am(M+1:2*M, 3*N+3:4*N+2, e,e) = -spread(Dbessi(1:N)/bessi(1:N),1,M)*&
                     & sin(outer_prod(CIPcm(1:M),rk(1:N))) ! d_n

                colmax = 4*N+2
             else
                colmax = 2*N+1
             end if

             Amain(row(1,e):row(2,e),col(1,e):col(2,e)) = Am(M+1:2*M,1:colmax,e,e)

          end select 
       end if
    end do ELEM

    !! off-diagonal sub-matricies 
    !! effect of source element on circumference of target element
    !! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    TARG: do t = 1,ni  ! "target" element
       SRC: do s = 1,ni  ! "source" element

          ! all three cases lumped together into a more general framework
          ! to minimize code duplication

          ! @@@@@ source in outside background of target @@@@@ (s & t have same parent)
          !!               -or-
          ! @@@@@ source in inside background target @@@@@@  (t is parent of s)
          !!               -or-
          ! @@@@@ target element has source element as parent @@@@ 
          if (CIInclBg(s,t) .or. (CIInclIn(t,s) .and. CIcalcin(t)) .or. &
               & (CIInclUp(t) == s)) then

             if(t == CIInclUp(s))then  
                par = t
             elseif(s == CIInclUp(t)) then
                par = s
             else
                ! t & s have same parent (some other element)
                par = CIInclUp(t) 
             end if

             if(par /= s) then
                ! pre-compute Bessel-K
                besskRcm(1:M,0:N,s) = spread(besselk(CIr(s)*q(par),0,N+1),1,M)
                besskRgm(1:M,0:N+1,s,t) = besselk(CIRgm(1:M,s,t)*q(par),0,N+2)
                besskRgm(1:M,-1,s,t) = besskRgm(1:M,+1,s,t)

                ! ratio of deriv Kn(r) to Kn(r_0)
                BessRat1(1:M,0:N) = -q(par)*&
                     & (besskRgm(1:M,-1:N-1,s,t) + besskRgm(1:M,1:N+1,s,t))/&
                     & (CTWO*besskRcm(1:M,0:N,s))
                ! ratio of Kn(r) to Kn(r_0)
                BessRat2(1:M,0:N) = besskRgm(1:M,0:N,s,t)/besskRcm(1:M,0:N,s)
             else
                ! pre-compute Bessel-I
                bessiRcm(1:M,0:N,par) = spread(besseli(CIr(par)*q(par),0,N+1),1,M)
                bessiRgm(1:M,0:N+1,par,t) = besseli(CIRgm(:,par,t)*q(par),0,N+2)
                bessiRgm(1:M,-1,par,t) = bessiRgm(:,+1,par,t)

                ! ratio of deriv In(r) to In(r_0)
                BessRat1(1:M,0:N) = q(par)*&
                     & (bessiRgm(:,-1:N-1,par,t) + bessiRgm(1:M,1:N+1,par,t))/&
                     & (CTWO*bessiRcm(1:M,0:N,par))
                ! ratio of In(r) to In(r_0)
                BessRat2(1:M,0:N) = bessiRgm(1:M,0:N,par,t)/bessiRcm(1:M,0:N,par)
             end if

             ! pre-compute sin/cos
             cosnPgm(1:M,0:N) = cos(outer_prod(CIPgm(1:M,s,t),rk(0:N)))
             sinnPgm(1:M,0:N) = sin(outer_prod(CIPgm(1:M,s,t),rk(0:N)))

             !##########################################
             ! head effects on outside of target element
             !##########################################
             if (CIibnd(t) /= +1) then ! not specified flux 
                ! potential effects of source inlcusions
                Am(1:M, 1:N+1, t,s) =     BessRat2(1:M,0:N)*cosnPgm(1:M,0:N)
                Am(1:M, N+2:2*N+1, t,s) = BessRat2(1:M,1:N)*sinnPgm(1:M,1:N)
             end if

             !##########################################
             ! flux effects on outside of target element
             !##########################################
             if (CIibnd(t) /= -1) then ! not constant head 

                cosPm(1:M) = cos(CIPcm(:))
                sinPm(1:M) = sin(CIPcm(:))
                cosGm(1:M) = cos(CIPgm(:,s,t))
                sinGm(1:M) = sin(CIPgm(:,s,t))

                ! deriv wrt radius of source to deriv wrt r of target
                Jrr(1:M) = cosGm(:)*cosPm(:) + sinGm(:)*sinPm(:)

                ! deriv wrt theta of source to deriv wrt r of target
                Jtr(1:M) = -sinGm(:)/CIRgm(:,s,t)*cosPm(:) + &
                        &   cosGm(:)/CIRgm(:,s,t)*sinPm(:)

                !a_n
                Am(M+1:2*M, 1:N+1, t,s) = &
                     &  spread(Jrr(1:M),2,N+1)*cosnPgm(:,0:N)*BessRat1(:,0:N) - &
                     & outer_prod(Jtr(1:M),rk(0:N))*sinnPgm(1:M,0:N)*BessRat2(:,0:N) 

                !b_n
                Am(M+1:2*M, N+2:2*N+1, t,s) = &
                     &  spread(Jrr(1:M),2,N)*sinnPgm(:,1:N)*BessRat1(:,1:N) + &
                     & outer_prod(Jtr(1:M),rk(1:N))*cosnPgm(:,1:N)*BessRat2(:,1:N) 
             end if

             ! apply sign and 1/k to head portions of sub-matricies not on diagonals
             if(CIMatch(t)) then
                if(par /= s) then
                   Am(1:M, 1:2*N+1, t,s) = -Am(1:M, 1:2*N+1, t,s)
                endif
                Am(1:M, 1:2*N+1, t,s) = Am(1:M, 1:2*N+1, t,s)/kv(par)
             end if

             if(par == s) then
                select case(CIibnd(t))
                case(0) ! matching: head + flux effects on outside of target (a_n and b_n)
                   Amain(row(1,t):row(2,t),col(2,s)-(2*N):col(2,s)) = Am(1:2*M,1:2*N+1,t,s)
                case(+1) ! spec. flux: flux effects on outside of target
                   Amain(row(1,t):row(2,t),col(2,s)-(2*N):col(2,s)) = Am(1:M,1:2*N+1,t,s)
                case(-1) ! spec. head: head effects on outside of target
                   Amain(row(1,t):row(2,t),col(2,s)-(2*N):col(2,s)) = Am(1:M,1:2*N+1,t,s)
                end select
             else
                select case(CIibnd(t))
                case(0) ! matching:  head + flux effects on inside of target (c_n and d_n)
                   Amain(row(1,t):row(2,t),col(1,s):col(1,s)+2*N) = Am(1:2*M,1:2*N+1,t,s)
                case(+1) ! spec. flux: flux effects on inside of target
                   Amain(row(1,t):row(2,t),col(1,s):col(1,s)+2*N) = Am(1:M,1:2*N+1,t,s)
                case(-1) ! spec. head: head effects on inside of target
                   Amain(row(1,t):row(2,t),col(1,s):col(1,s)+2*N) = Am(1:M,1:2*N+1,t,s)
                end select
             end if
          end if
       end do SRC
    end do TARG

    call computeCircularStrengths(Amain,p,coeff,counter)
  end subroutine CircInverse_matrix


  !##################################################
  ! this routine takes the already-computed  inverse and 
  ! computes the RHS vector of knowns, then computes x
  subroutine computeCircularStrengths(Amain,p,coeff,counter)
    use shared_matching_data, only : row,col
    use constants, only : DP, TWOPI, CZERO, RZERO, RONE
    use wells
    use element_specs, only : Cin,CInum,CIm, av,kv,sv,WLnum, &
         & CIr,CIInclUp,CIibnd,CIarea,CICalcIn,CIspec,CIWellIn,CIWellBg
    use bessel_functions

    !! solve over-determined system via least-squares
    interface
       subroutine ZGELS(TRANSA, M, N, NRHS, A, LDA, B, LDB, WORK, &
            & LDWORK, INFO)
         integer, intent(in) :: M, N, NRHS, LDA, LDB, LDWORK
         character(LEN=1), intent(in) :: TRANSA
         complex(8), intent(inout), dimension(LDWORK) :: WORK
         complex(8), intent(inout), dimension(LDA,N) :: A
         complex(8), intent(inout), dimension(LDB,NRHS) ::  B
         integer, intent(out) :: INFO
       end subroutine ZGELS
    end interface

    integer, intent(in) :: counter
    complex(DP), dimension(:,:), intent(inout) :: Amain !! coefficient matrix 
    complex(DP), intent(in) :: p     ! scalar Laplace parameter

    ! generalized Fourier coefficients (previous value passed in as first guess)
    complex(DP), dimension(0:4*CIn+1,CInum), intent(inout) :: coeff

    complex(DP), dimension(size(Amain,dim=1)) :: bm   ! rhs: specified BC - well effects
    !! this contains the solution vector after call to ZGELS

    ! results from subroutines for passive circular elements
    complex(DP), dimension(CIm,CInum) :: PotWellIn,FluxWellIn,PotWellBg,FluxWellBg

    ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
    integer :: N, M, i, w, ni, nw, par, bigM, bigN

    ! things only needed for LAPACK routine
    integer, parameter :: iwork = 10250
    complex(DP), dimension(iwork) :: WORK
    integer :: IERR

    N = CIn; M = CIm; ni = CInum; nw = WLnum

    q(1:ni) = sqrt(p/av(1:ni))

    if (nw > 0) then
       call wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    end if

    bm = CZERO
    do i = 1,ni

       ! put specified BC onto RHS
       select case(CIibnd(i))
       case(-1) ! specified head
          bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) + CIspec(i)*CircTimeBdry(i,p)
       case(+1) ! specified flux
          bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) + &
               & CIspec(i)*CircTimeBdry(i,p)/(TWOPI*CIr(i))
       case(0)  ! matching: area flux on inside of element (no flux effect)
          bm(row(1,i)+M:row(2,i)) = bm(row(1,i)+M:row(2,i)) - &
               & CircTimeArea(i,p)*CIarea(i)*sv(i)/q(i)**2
       end select

       ! put well effects onto RHS
       do w = 1,nw
          if(CIWellIn(i,w)) then ! well is in inside background of element
             select case(CIibnd(i))
             case(0) ! matching element
                bm(row(1,i):row(1,i)+M-1) = bm(row(1,i):row(1,i)+M-1) + &
                     & PotWellIn(:,i)/kv(i)
                bm(row(1,i)+M:row(2,i))   = bm(row(1,i)+M:row(2,i)) + FluxWellIn(:,i)
             case(-1) ! specified head element
                bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) + PotWellIn(:,i)/kv(i)
             case(+1) ! specified flux element
                bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) + FluxWellIn(:,i)
             end select
          elseif(CIWellBg(i,w)) then  ! well is in outside background of element
             par = CIInclUp(i)
             select case(CIibnd(i))
             case(0) ! matching
                bm(row(1,i):row(1,i)+M-1) = bm(row(1,i):row(1,i)+M-1) - &
                     & PotWellBg(:,i)/kv(par)
                bm(row(1,i)+M:row(2,i)) = bm(row(1,i)+M:row(2,i)) - FluxWellBg(:,i)
             case(-1) ! specified head element
                bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) - PotWellBg(:,i)/kv(par)
             case(+1) ! specified flux element
                bm(row(1,i):row(2,i)) = bm(row(1,i):row(2,i)) - FluxWellBg(:,i)
             end select
          end if
       end do

       ! put passive circular element effects onto RHS
       ! <<<< TODO >>>>

    end do

    bigM = size(Amain,dim=1)
    bigN = size(Amain,dim=2)
    !! solve least-squares problem
    call ZGELS('N',bigM,bigN,1,Amain,bigM,bm(:),bigM,WORK,IWORK,ierr)
    if(ierr /= 0) print *, 'ZGELS error:',ierr,' p:',p
       
    ! put results into coeff variable
    do i = 1,ni
       select case(CIibnd(i))
       case(0) ! matching
          coeff(0:N,i) =         bm(col(1,i):col(1,i)+N) !a_n
          coeff(N+1:2*N,i) =     bm(col(1,i)+N+1:col(1,i)+2*N) !b_n
          coeff(2*N+1:3*N+1,i) = bm(col(1,i)+2*N+1:col(1,i)+3*N+1) !c_n
          coeff(3*N+2:4*N+1,i) = bm(col(1,i)+3*N+2:col(1,i)+4*N+1) !d_n
       case(-1,+1) ! spec head/flux
          coeff(0:N,i) =         bm(col(1,i):col(1,i)+N) !a_n
          coeff(N+1:2*N,i) =     bm(col(1,i)+N+1:col(1,i)+2*N) !b_n
          if(CIcalcin(i)) then
             coeff(2*N+1:3*N,i) =   bm(col(1,i)+2*N+1:col(1,i)+3*N+1) !c_n
             coeff(3*N+1:4*N+1,i) = bm(col(1,i)+3*N+2:col(2,i)) !d_n
          end if
       end select
    end do

  end subroutine computeCircularStrengths

  !##################################################
  ! calcuale effects of wells (passive elements)
  ! this subroutine called once for each value of p before beginning iteration
  subroutine wellEffects(p,PotWellIn,FluxWellIn,PotWellBg,FluxWellBg)
    use shared_matching_data, only : CIPcm, CIRwm, CIPwm
    use constants, only: DP, CZERO
    use element_specs, only : CIm,CInum,WLnum,av,CImatch,CIInclUp,CIWellIn,CIWellBg
    use wells ! generalized well-effect functions

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
    ! scalar Laplace parameter
    complex(DP), intent(in) :: p
    ! potential and flux due to wells inside and in background of each circular element
    complex(DP), dimension(CIm,CInum), intent(out) :: PotWellIn,FluxWellIn,&
         & PotWellBg,FluxWellBg

    ! INTERNAL VARIABLES
    ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
    ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(CIm) :: dPot_dRw, dPot_dX, dPot_dY
    integer :: M, well, inc, par, ni, nw

    M = CIm; ni = CInum; nw = WLnum;
    q = sqrt(p/av)

    PotWellIn = CZERO; FluxWellIn = CZERO;
    PotWellBg = CZERO; FluxWellBg = CZERO;

    !******************** well effects (head and flux) ********************
    do inc = 1,ni
       if (CImatch(inc)) then
          par = CIInclUp(inc);
          do well = 1,nw 
             !^^^^^^^^^^^^^^^^^^^^ well inside this inclusion ^^^^^^^^^^^^^^^^^^^^
             if (CIWellIn(inc,well)) then

                ! potential due to wells
                PotWellIn(1:M,inc) = PotWellIn(1:M,inc) + &
                     &  wellHead(well,p,q(inc),CIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,q(inc),CIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(CIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local polar coords and add to total
                FluxWellIn(1:M,inc) = FluxWellIn(1:M,inc) + &
                     &  dPot_dX(1:M)*cos(CIPcm(1:M)) + dPot_dY(1:M)*sin(CIPcm(1:M))

                !^^^^^^^^^^^^^^^^^^^^ well in bg of this inclusion ^^^^^^^^^^^^^^^^^^^^
             elseif (CIWellBg(inc,well)) then

                ! potential due to wells
                PotWellBg(1:M,inc) = PotWellBg(1:M,inc) + &
                     &  wellHead(well,p,q(par),CIRwm(:,well,inc))

                ! radial flux due to wells
                dPot_dRw(1:M) = wellFlux(well,p,q(par),CIRwm(:,well,inc))

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPwm(:,well,inc))*dPot_dRw(1:M)
                dPot_dY(1:M) = sin(CIPwm(:,well,inc))*dPot_dRw(1:M)

                ! project onto local polar coords and add to total
                FluxWellBg(1:M,inc) = FluxWellBg(1:M,inc) + &
                     &  dPot_dX(1:M)*cos(CIPcm(1:M)) + dPot_dY(1:M)*sin(CIPcm(1:M))
             else
                ! if well is inside a different element, do nothing
             end if
          end do
       end if
    end do
  end subroutine wellEffects

  !########################################
  ! free memory related to points along circular inclusions 
  ! (not needed after coeff are calculated)
  subroutine freematchmem()
    use shared_matching_data
    deallocate(CIXcm, CIYcm, CIXom, CIYom, CIRwm, CIPwm, CIRgm, CIPgm, CIPcm,row,col)
  end subroutine freematchmem

end module matching_new
