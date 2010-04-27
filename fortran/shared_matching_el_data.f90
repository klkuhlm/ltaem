! $Id: shared_matching_el_data.f90,v 1.3 2007/07/29 04:33:18 kris Exp kris $
module shared_matching_el_data
  use constants, only : DP
  implicit none

  public

  ! things which relate to positions on circumference of ellipse
  ! which only need to be seen within matching routines

  ! geometry: distances and angles
  real(DP), save, allocatable :: EIPcm(:), EIXcm(:,:), EIYcm(:,:), EIXom(:,:), EIYom(:,:)
  real(DP), save, allocatable :: EIRwm(:,:,:), EIPwm(:,:,:), EIEgm(:,:,:), EIPgm(:,:,:)

  integer, save, allocatable :: row(:,:), col(:,:)

end module shared_matching_el_data
