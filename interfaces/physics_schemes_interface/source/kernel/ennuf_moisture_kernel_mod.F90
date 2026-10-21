!-------------------------------------------------------------------------------
! Proportionally partition ENNUF total-water increments into moisture species.
!-------------------------------------------------------------------------------
module ennuf_moisture_kernel_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use argument_mod, only: arg_type, GH_FIELD, GH_REAL, GH_LOGICAL, GH_SCALAR, &
                         GH_READ, GH_READWRITE, CELL_COLUMN
  use constants_mod, only: r_def, i_def, l_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod, only: kernel_type
  implicit none
  private
  public :: ennuf_moisture_kernel_type, ennuf_moisture_code, partition_water_increment

  type, extends(kernel_type) :: ennuf_moisture_kernel_type
    type(arg_type) :: meta_args(10) = (/ &
      arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! vapour
      arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! liquid
      arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! snow/large ice
      arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! small ice
      arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta), & ! input vapour
      arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta), & ! input liquid
      arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta), & ! input snow/large ice
      arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta), & ! input small ice
      arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta), & ! total-water increment
      arg_type(GH_SCALAR, GH_LOGICAL, GH_READ) /) ! CASIM includes small ice
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: code => ennuf_moisture_code
  end type ennuf_moisture_kernel_type

contains

  subroutine partition_water_increment(initial, current, increment, updated)
    real(r_def), intent(in) :: initial(:), current(:), increment
    real(r_def), intent(out) :: updated(size(current))
    real(r_def) :: fractions(size(initial)), applied, total
    integer :: species

    if (size(initial) /= size(current)) error stop 'ENNUF moisture species size mismatch'
    if (.not. all(ieee_is_finite(initial)) .or. .not. all(ieee_is_finite(current)) .or. &
        .not. ieee_is_finite(increment)) error stop 'Non-finite ENNUF moisture values'
    if (any(initial < 0.0_r_def) .or. any(current < 0.0_r_def)) then
      error stop 'ENNUF moisture partition requires non-negative water species'
    end if
    updated = current
    total = sum(initial)
    if (total <= 0.0_r_def .or. increment == 0.0_r_def) return
    fractions = initial / total
    applied = increment
    if (applied < 0.0_r_def) then
      ! Other physics may have consumed a component since the inference snapshot.
      do species = 1, size(initial)
        if (fractions(species) > 0.0_r_def) then
          applied = max(applied, -max(current(species), 0.0_r_def) / fractions(species))
        end if
      end do
    end if
    updated = max(current + applied * fractions, 0.0_r_def)
  end subroutine partition_water_increment

  subroutine ennuf_moisture_code(nlayers, vapour, liquid, snow, ice, &
                                initial_vapour, initial_liquid, initial_snow, initial_ice, &
                                increment, include_ice, ndf, undf, map)
    integer(i_def), intent(in) :: nlayers, ndf, undf, map(ndf)
    real(r_def), intent(inout) :: vapour(undf), liquid(undf), snow(undf), ice(undf)
    real(r_def), intent(in) :: initial_vapour(undf), initial_liquid(undf)
    real(r_def), intent(in) :: initial_snow(undf), initial_ice(undf), increment(undf)
    logical(l_def), intent(in) :: include_ice
    real(r_def) :: initial(4), current(4), updated(4)
    integer(i_def) :: level, index

    do level = 1, nlayers
      index = map(1) + level
      if (increment(index) == 0.0_r_def) cycle
      initial = [initial_vapour(index), initial_liquid(index), initial_snow(index), 0.0_r_def]
      current = [vapour(index), liquid(index), snow(index), 0.0_r_def]
      if (include_ice) then
        initial(4) = initial_ice(index)
        current(4) = ice(index)
      end if
      call partition_water_increment(initial, current, increment(index), updated)
      vapour(index) = updated(1)
      liquid(index) = updated(2)
      snow(index) = updated(3)
      if (include_ice) ice(index) = updated(4)
    end do
  end subroutine ennuf_moisture_code
end module ennuf_moisture_kernel_mod
