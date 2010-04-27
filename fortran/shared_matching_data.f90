! $Id: shared_matching_data.f90,v 1.1 2006/06/02 01:14:07 kris Exp $
module shared_matching_data
  use constants, only : DP
  implicit none

  public

  ! things which relate to positions on circumference of circle
  ! which only need to be seen within matching routines

  ! geometry: distances and angles
  real(DP), save, allocatable :: CIPcm(:), CIXcm(:,:), CIYcm(:,:), CIXom(:,:), CIYom(:,:)
  real(DP), save, allocatable :: CIRwm(:,:,:), CIPwm(:,:,:), CIRgm(:,:,:), CIPgm(:,:,:)

  integer, save, allocatable :: row(:,:), col(:,:)

end module shared_matching_data
