module circ_passive_elements
! $Id: circ_passive_elements.f90,v 1.1 2006/06/02 01:07:48 kris Exp kris $

private
public CircFluxElementEffects, CircHeadElementEffects

contains

  !##################################################
  ! effects of passive circular elements with specified flux
  ! this subroutine called once for each value of p before beginning iteration
  subroutine CircFluxElementEffects(p,PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg,coeff)
    use shared_matching_data, only : CIRgm, CIPgm, CIPcm
    use constants, only: DP, RONE, RTWO, CZERO, TWOPI
    use bessel_functions, only : besselk, besseli
    use wells, only : circTimeBdry
    use element_specs, only : CIm,CIn,CInum,av,CIibnd,CIInclUp,CIInclBg,CIInclIn,CImatch,&
         & CIr,CIspec

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter    
    complex(DP), intent(in) :: p
      ! potential and flux due to ring elements inside and in background of circular elements
    complex(DP), dimension(CIm,CInum), intent(inout) :: PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg
      ! generalized Fourier coefficients for passive circular elements
    complex(DP), dimension(0:4*CIn+1,CInum), intent(inout) :: coeff 

    ! INTERNAL VARIABLES
      ! vectors of Bessel functions to minimize calls to Bessel fcn routine
    complex(DP), dimension(CIm,0:1) :: besskRgm, bessiRgm
    complex(DP), dimension(0:1) :: besskRcm, bessiRcm
      ! sqrt(p/alpha)    
    complex(DP), dimension(0:CInum) :: q
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(CIm) :: dPot_dRg, dPot_dX, dPot_dY
    integer :: M, incl, ni, other, parent, N

    M = CIm; ni = CInum; N = CIn
    q = sqrt(p/av)

    !^^^^^^^^ specified circular element effects ^^^^^^^^
    do incl = 1,ni 
       if (CIibnd(incl) == -2) then ! loop over specified elemental flux inclusions

          parent = CIInclUp(incl)

          do other = 1,ni  ! calculate effects spec elem has on other (matching incl)

             ! ***** other matching inclusion in __background__ of this one *****
             if (CIInclBg(incl,other) .and. CImatch(other)) then
                
                ! pre-calculate K Bessel functions
                besskRgm(1:M,0:1) = besselk(CIRgm(:,incl,other)*q(parent),0,2) !vector
                besskRcm(0:1) = besselk(CIr(incl)*q(parent),0,2) !scalar

                ! this is currently just a ring source which is constant in theta
                coeff(0,incl) = -CIspec(incl)*CircTimeBdry(incl,p)/&
                     & (TWOPI*CIr(incl)*q(parent))*besskRcm(0)/besskRcm(1)  ! a_0
                coeff(1:2*N,incl) = CZERO  ! a_1:n and b_n

                ! potential due to circ elements
                PotCElmBg(1:M,other) = PotCelmBg(:,other) - &
                     & coeff(0,incl)*besskRgm(:,0)/besskRcm(0)

                ! radial flux due to circ elements
                dPot_dRg(1:M) = -coeff(0,incl)*q(parent)*besskRgm(:,1)/besskRcm(0)

                ! azimuthal flux == 0, since uniform specified head is symmetric wrt theta

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPgm(:,incl,other))*dPot_dRg(:)
                dPot_dY(1:M) = sin(CIPgm(:,incl,other))*dPot_dRg(:)
                
                ! project onto local polar coords and add to cum total
                FluxCElmBg(1:M,other) = FluxCElmBg(:,incl) + &
                     &  dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
                
             ! ***** other inclusion __inside__ of this one *****
             elseif (CIInclIn(incl,other) .and. CImatch(other)) then

                ! pre-calculate I Bessel functions
                bessiRgm(1:M,0:1) = besseli(CIRgm(:,incl,other)*q(incl),0,2)
                bessiRcm(0:1) = besseli(CIr(incl)*q(incl),0,2)

                coeff(2*N+1,incl) = CIspec(incl)*CircTimeBdry(incl,p)/&
                     &(TWOPI*CIr(incl)*q(incl))*bessiRcm(0)/bessiRcm(1) ! c_0
                coeff(2*N+2:4*N+1,incl) = CZERO ! c_1:n and d_n

                ! potential due to circ elements
                PotCElmIn(1:M,other) = PotCElmIn(:,other) + &
                     & coeff(2*N+1,incl)*bessiRgm(:,0)/bessiRcm(0)

                ! radial flux due to circ elements
                dPot_dRg(1:M) = coeff(2*N+1,incl)*q(incl)/bessiRcm(0)*bessiRgm(:,1)

                ! azimuthal flux == 0, since uniform specified head is symmetric wrt theta
                
                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPgm(:,incl,other))*dPot_dRg(:)
                dPot_dY(1:M) = sin(CIPgm(:,incl,other))*dPot_dRg(:)
                
                ! project onto local polar coords and add to cum total
                FluxCElmIn(1:M,other) = FluxCElmIn(:,incl) + &
                     &  dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
             end if
          end do
       end if   
    end do
  end subroutine CircFluxElementEffects
  

  !##################################################
  ! effects of passive circular elements with specified head
  ! this subroutine called once for each value of p before beginning iteration
  subroutine CircHeadElementEffects(p,PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg,coeff)
    use shared_matching_data, only : CIRgm, CIPgm, CIPcm
    use constants, only: DP, CZERO
    use wells, only : circTimeBdry
    use bessel_functions, only : besselk, besseli
    use element_specs, only : CIibnd,CIm,CInum,CIn,CIInclBg,CIInclIn,CImatch,av,kv,CIInclUp,CIr,&
         & CIspec

    ! ARGUMENTS EXTERNAL TO SUBROUTINE
      ! scalar Laplace parameter    
    complex(DP), intent(in) :: p
      ! potential and flux due to ring elements inside and in background of circular elements
    complex(DP), dimension(CIm,CInum), intent(inout) :: PotCElmIn,FluxCElmIn,PotCElmBg,FluxCElmBg
      ! generalized Fourier coefficients for passive circular elements
    complex(DP), dimension(0:4*CIn+1,CInum), intent(inout) :: coeff 

    ! INTERNAL VARIABLES
      ! vectors of Bessel functions to minimize calls to Bessel fcn routine
    complex(DP), dimension(CIm,0:1) :: besskRgm, bessiRgm
    complex(DP) :: besskRcm0, bessiRcm0
      ! sqrt(p/alpha)
    complex(DP), dimension(0:CInum) :: q
      ! flux projected onto cartesian coords (intermediate result)
    complex(DP), dimension(CIm) :: dPot_dRg, dPot_dX, dPot_dY
    integer :: M, incl, ni, other, parent, N

    M = CIm; ni = CInum; N = CIn
    q = sqrt(p/av)

    !^^^^^^^^ specified head circular element effects ^^^^^^^^
    do incl = 1,ni 
       if (CIibnd(incl) == +2) then ! loop over specified elemental head inclusions

          parent = CIInclUp(incl)

          do other = 1,ni  ! calculate effects spec elem has on other (matching incl)

             ! ***** other matching inclusion in __background__ of this one *****
             if (CIInclBg(incl,other) .and. CImatch(other)) then

                besskRgm(1:M,0:1) = besselk(CIRgm(:,incl,other)*q(parent),0,2)
                besskRcm0 = besselk(CIr(incl)*q(parent),0)

                coeff(0,incl) = CIspec(incl)*kv(parent)*CircTimeBdry(incl,p)  ! a_0
                coeff(1:2*N,incl) = CZERO  ! a_1:n and b_n

                ! potential due to circ elements
                PotCElmBg(1:M,other) = PotCelmBg(:,other) + &
                     & coeff(0,incl)*besskRgm(:,0)/besskRcm0
                
                ! radial flux due to circ elements
                dPot_dRg(1:M) = -q(parent)*coeff(0,incl)*besskRgm(:,1)/besskRcm0

                ! azimuthal flux == 0, since uniform specified head is symmetric wrt theta

                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPgm(:,incl,other))*dPot_dRg(:)
                dPot_dY(1:M) = sin(CIPgm(:,incl,other))*dPot_dRg(:)
                
                ! project onto local polar coords and add to cum total
                FluxCElmBg(1:M,other) = FluxCElmBg(:,incl) + &
                     &  dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
                
             ! ***** other inclusion __inside__ of this one *****
             elseif (CIInclIn(incl,other) .and. CImatch(other)) then

                bessiRgm(1:M,0:1) = besseli(CIRgm(:,incl,other)*q(incl),0,2)
                bessiRcm0 = besseli(CIr(incl)*q(incl),0)

                coeff(2*N+1,incl) = CIspec(incl)*kv(incl)*CircTimeBdry(incl,p) ! c_0
                coeff(2*N+2:4*N+1,incl) = CZERO ! c_1:n and d_n

                ! potential due to circ elements
                PotCElmIn(1:M,other) = PotCElmIn(:,other) + &
                     & coeff(2*N+1,incl)*bessiRgm(:,0)/bessiRcm0

                ! radial flux due to circ elements
                dPot_dRg(1:M) = q(incl)*coeff(2*N+1,incl)*bessiRgm(:,1)/besskRcm0

                ! azimuthal flux == 0, since uniform specified head is symmetric wrt theta
                
                ! project onto global Cartesian coords
                dPot_dX(1:M) = cos(CIPgm(:,incl,other))*dPot_dRg(:)
                dPot_dY(1:M) = sin(CIPgm(:,incl,other))*dPot_dRg(:)
                
                ! project onto local polar coords and add to cum total
                FluxCElmIn(1:M,other) = FluxCElmIn(:,incl) + &
                     &  dPot_dX(:)*cos(CIPcm(:)) + dPot_dY(:)*sin(CIPcm(:))
             end if
          end do
       end if   
    end do
  end subroutine CircHeadElementEffects


end module circ_passive_elements
