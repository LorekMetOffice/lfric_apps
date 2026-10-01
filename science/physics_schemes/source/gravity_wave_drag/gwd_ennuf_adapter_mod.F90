!-----------------------------------------------------------------------------
! ENNUF adapter for column residual prediction in gravity_wave_drag.
! This is a small integration layer around the generated ENNUF model modules.
!-----------------------------------------------------------------------------
module gwd_ennuf_adapter_mod
  use, intrinsic :: iso_fortran_env, only: real32, real64

  use log_mod, only: log_event, LOG_LEVEL_INFO
  use gwd_ennuf_scaling_mod, only: ennuf_standardise_input

  use theta_global_flatdense_bs512_20km_pressure_orogonly_mod, only: &
    theta_global_flatdense_bs512_20km_pressure_orogonly
  use theta_us2_flatdense_bs512_20km_pressure_orogonly_mod, only: &
      theta_us2_flatdense_bs512_20km_pressure_orogonly
  use qt_global_flatdense_bs512_20km_pressure_orogonly_mod, only: &
    qt_global_flatdense_bs512_20km_pressure_orogonly
  use qt_us2_flatdense_bs512_20km_pressure_orogonly_mod, only: &
      qt_us2_flatdense_bs512_20km_pressure_orogonly
  use x_wind_global_flatdense_bs512_20km_pressure_orogonly_mod, only: &
    x_wind_global_flatdense_bs512_20km_pressure_orogonly
  use x_wind_us2_flatdense_bs512_20km_pressure_orogonly_mod, only: &
      x_wind_us2_flatdense_bs512_20km_pressure_orogonly
  use y_wind_global_flatdense_bs512_20km_pressure_orogonly_mod, only: &
    y_wind_global_flatdense_bs512_20km_pressure_orogonly
  use y_wind_us2_flatdense_bs512_20km_pressure_orogonly_mod, only: &
      y_wind_us2_flatdense_bs512_20km_pressure_orogonly

  implicit none

  private
  public :: thermodynamic_correction, build_ennuf_input, predict_column_residuals
  public :: set_nudge_options

  integer, parameter :: ml_profile_len = 53
  integer, parameter :: ml_input_len = 268
  real(real64), parameter :: cp = 1005.0_real64
  real(real64), parameter :: lc = 2.501e6_real64
  real(real64), parameter :: lf = 0.334e6_real64
  real(real64), parameter :: p0 = 1.0e5_real64
  real(real64), parameter :: kappa = 0.286_real64
  real(real64), parameter :: two_pi = 6.28318530717958647692528676655900577_real64

  logical, save :: nudge_apply_globally = .true.
  logical, save :: nudge_use_theta = .true.
  logical, save :: nudge_use_qt = .true.
  logical, save :: nudge_use_u = .true.
  logical, save :: nudge_use_v = .true.
  real(real64), save :: nudge_lat_min = -90.0_real64
  real(real64), save :: nudge_lat_max = 90.0_real64
  real(real64), save :: nudge_lon_min = -180.0_real64
  real(real64), save :: nudge_lon_max = 180.0_real64

contains

  subroutine set_nudge_options(apply_globally, use_theta, use_qt, use_u, use_v, &
                               lat_min, lat_max, lon_min, lon_max)
    logical, intent(in), optional :: apply_globally
    logical, intent(in), optional :: use_theta, use_qt, use_u, use_v
    real(real64), intent(in), optional :: lat_min, lat_max, lon_min, lon_max

    if (present(apply_globally)) nudge_apply_globally = apply_globally
    if (present(use_theta)) nudge_use_theta = use_theta
    if (present(use_qt)) nudge_use_qt = use_qt
    if (present(use_u)) nudge_use_u = use_u
    if (present(use_v)) nudge_use_v = use_v
    if (present(lat_min)) nudge_lat_min = lat_min
    if (present(lat_max)) nudge_lat_max = lat_max
    if (present(lon_min)) nudge_lon_min = lon_min
    if (present(lon_max)) nudge_lon_max = lon_max
  end subroutine set_nudge_options

  pure function thermodynamic_correction(theta, qcl, qcf) result(theta_li)
    real(real64), intent(in) :: theta(ml_profile_len)
    real(real64), intent(in) :: qcl(ml_profile_len)
    real(real64), intent(in) :: qcf(ml_profile_len)
    real(real64) :: theta_li(ml_profile_len)

    theta_li = theta - (lc / cp) * qcl - ((lc + lf) / cp) * qcf
  end function thermodynamic_correction

  pure function pressure_from_exner(exner) result(pressure)
    real(real64), intent(in) :: exner(ml_profile_len)
    real(real64) :: pressure(ml_profile_len)

    pressure = p0 * exner ** kappa
  end function pressure_from_exner

  pure function compute_qt(qv, qcl, qcf) result(qt)
    real(real64), intent(in) :: qv(ml_profile_len)
    real(real64), intent(in) :: qcl(ml_profile_len)
    real(real64), intent(in) :: qcf(ml_profile_len)
    real(real64) :: qt(ml_profile_len)

    qt = qv + qcl + qcf
  end function compute_qt

  subroutine build_ennuf_input(theta, qt, u, v, pressure, orography, seconds_of_day, combined_input)
    real(real64), intent(in) :: theta(ml_profile_len)
    real(real64), intent(in) :: qt(ml_profile_len)
    real(real64), intent(in) :: u(ml_profile_len)
    real(real64), intent(in) :: v(ml_profile_len)
    real(real64), intent(in) :: pressure(ml_profile_len)
    real(real64), intent(in) :: orography
    real(real64), intent(in) :: seconds_of_day
    real(real32), intent(out) :: combined_input(ml_input_len)

    real(real64) :: time_angle

    time_angle = two_pi * modulo(seconds_of_day, 86400.0_real64) / 86400.0_real64

    combined_input(1:53) = real(theta, real32)
    combined_input(54:106) = real(qt, real32)
    combined_input(107:159) = real(u, real32)
    combined_input(160:212) = real(v, real32)
    combined_input(213:265) = real(pressure, real32)
    combined_input(266) = real(orography, real32)
    combined_input(267) = real(sin(time_angle), real32)
    combined_input(268) = real(cos(time_angle), real32)
  end subroutine build_ennuf_input

  subroutine predict_column_residuals(theta, qv, qcl, qcf, u, v, pressure, orography, &
                                      seconds_of_day, theta_residual, qt_residual, &
                                      u_residual, v_residual,                   &
                                      use_theta, use_qt, use_u, use_v,          &
                                      apply_globally, latitude, longitude,      &
                                      lat_min, lat_max, lon_min, lon_max)
    real(real64), intent(in) :: theta(ml_profile_len)
    real(real64), intent(in) :: qv(ml_profile_len)
    real(real64), intent(in) :: qcl(ml_profile_len)
    real(real64), intent(in) :: qcf(ml_profile_len)
    real(real64), intent(in) :: u(ml_profile_len)
    real(real64), intent(in) :: v(ml_profile_len)
    real(real64), intent(in) :: pressure(ml_profile_len)
    real(real64), intent(in) :: orography
    real(real64), intent(in) :: seconds_of_day
    real(real64), intent(out) :: theta_residual(ml_profile_len)
    real(real64), intent(out) :: qt_residual(ml_profile_len)
    real(real64), intent(out) :: u_residual(ml_profile_len)
    real(real64), intent(out) :: v_residual(ml_profile_len)
    logical, intent(in), optional :: use_theta, use_qt, use_u, use_v
    logical, intent(in), optional :: apply_globally
    real(real64), intent(in), optional :: latitude, longitude
    real(real64), intent(in), optional :: lat_min, lat_max, lon_min, lon_max

    real(real64) :: theta_li(ml_profile_len)
    real(real64) :: qt(ml_profile_len)
    real(real32) :: combined_input(ml_input_len)
    real(real64) :: combined_input_dbl(ml_input_len)
    real(real64) :: combined_input_global_scaled(ml_input_len)
    real(real64) :: combined_input_us2_scaled(ml_input_len)
    real(real32) :: combined_input_global_for_model(ml_input_len)
    real(real32) :: combined_input_us2_for_model(ml_input_len)
    real(real32) :: theta_global_scaled(1, ml_profile_len)
    real(real32) :: theta_us2_scaled(1, ml_profile_len)
    real(real32) :: qt_global_scaled(1, ml_profile_len)
    real(real32) :: qt_us2_scaled(1, ml_profile_len)
    real(real32) :: u_global_scaled(1, ml_profile_len)
    real(real32) :: u_us2_scaled(1, ml_profile_len)
    real(real32) :: v_global_scaled(1, ml_profile_len)
    real(real32) :: v_us2_scaled(1, ml_profile_len)
    logical :: allow_apply
    logical :: in_region
    logical :: select_theta
    logical :: select_qt
    logical :: select_u
    logical :: select_v
    character(len=256) :: nudge_region_msg
    real(real64) :: active_lat_min
    real(real64) :: active_lat_max
    real(real64) :: active_lon_min
    real(real64) :: active_lon_max

    call log_event('ENNUF adapter module loaded: gwd_ennuf_adapter_mod', LOG_LEVEL_INFO)

    theta_li = thermodynamic_correction(theta, qcl, qcf)
    qt = compute_qt(qv, qcl, qcf)

    select_theta = nudge_use_theta
    select_qt = nudge_use_qt
    select_u = nudge_use_u
    select_v = nudge_use_v
    if (present(use_theta)) select_theta = use_theta
    if (present(use_qt)) select_qt = use_qt
    if (present(use_u)) select_u = use_u
    if (present(use_v)) select_v = use_v

    allow_apply = nudge_apply_globally
    if (present(apply_globally)) allow_apply = apply_globally

    in_region = .true.
    if (.not. allow_apply) then
      in_region = .false.
      if (present(latitude) .and. present(longitude)) then
        if (latitude >= merge(lat_min, nudge_lat_min, present(lat_min)) .and. &
            latitude <= merge(lat_max, nudge_lat_max, present(lat_max)) .and. &
            longitude >= merge(lon_min, nudge_lon_min, present(lon_min)) .and. &
            longitude <= merge(lon_max, nudge_lon_max, present(lon_max))) then
          in_region = .true.
        end if
      end if
    end if

    call build_ennuf_input(theta_li, qt, u, v, pressure, orography, seconds_of_day, combined_input)
    combined_input_dbl = real(combined_input, real64)

    call ennuf_standardise_input(combined_input_dbl, 'global', combined_input_global_scaled)
    call ennuf_standardise_input(combined_input_dbl, 'us2', combined_input_us2_scaled)
    combined_input_global_for_model = real(combined_input_global_scaled, real32)
    combined_input_us2_for_model = real(combined_input_us2_scaled, real32)

    call theta_global_flatdense_bs512_20km_pressure_orogonly(combined_input_global_for_model, theta_global_scaled)
    call theta_us2_flatdense_bs512_20km_pressure_orogonly(combined_input_us2_for_model, theta_us2_scaled)
    call qt_global_flatdense_bs512_20km_pressure_orogonly(combined_input_global_for_model, qt_global_scaled)
    call qt_us2_flatdense_bs512_20km_pressure_orogonly(combined_input_us2_for_model, qt_us2_scaled)
    call x_wind_global_flatdense_bs512_20km_pressure_orogonly(combined_input_global_for_model, u_global_scaled)
    call x_wind_us2_flatdense_bs512_20km_pressure_orogonly(combined_input_us2_for_model, u_us2_scaled)
    call y_wind_global_flatdense_bs512_20km_pressure_orogonly(combined_input_global_for_model, v_global_scaled)
    call y_wind_us2_flatdense_bs512_20km_pressure_orogonly(combined_input_us2_for_model, v_us2_scaled)

    if (in_region) then
      theta_residual = 0.0_real64
      qt_residual = 0.0_real64
      u_residual = 0.0_real64
      v_residual = 0.0_real64

      if (present(lat_min) .or. present(lat_max) .or. present(lon_min) .or. present(lon_max)) then
        active_lat_min = merge(lat_min, nudge_lat_min, present(lat_min))
        active_lat_max = merge(lat_max, nudge_lat_max, present(lat_max))
        active_lon_min = merge(lon_min, nudge_lon_min, present(lon_min))
        active_lon_max = merge(lon_max, nudge_lon_max, present(lon_max))
      else
        active_lat_min = nudge_lat_min
        active_lat_max = nudge_lat_max
        active_lon_min = nudge_lon_min
        active_lon_max = nudge_lon_max
      end if

      if (select_theta .or. select_qt .or. select_u .or. select_v) then
        if (select_theta .and. .not. select_qt .and. .not. select_u .and. .not. select_v) then
          call log_event('hecate nudge - theta=on qt=off u=off v=off', LOG_LEVEL_INFO)
        else if (.not. select_theta .and. select_qt .and. .not. select_u .and. .not. select_v) then
          call log_event('hecate nudge - theta=off qt=on u=off v=off', LOG_LEVEL_INFO)
        else if (.not. select_theta .and. .not. select_qt .and. select_u .and. .not. select_v) then
          call log_event('hecate nudge - theta=off qt=off u=on v=off', LOG_LEVEL_INFO)
        else if (.not. select_theta .and. .not. select_qt .and. .not. select_u .and. select_v) then
          call log_event('hecate nudge - theta=off qt=off u=off v=on', LOG_LEVEL_INFO)
        else
          call log_event(                                                  &
            'hecate nudge - theta=' // merge('on', 'off', select_theta) // &
            ' qt=' // merge('on', 'off', select_qt) //                     &
            ' u=' // merge('on', 'off', select_u) //                       &
            ' v=' // merge('on', 'off', select_v), LOG_LEVEL_INFO)
        end if

        write(nudge_region_msg, '(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)')        &
          'hecate nudge region: lat=[', active_lat_min, ',', active_lat_max, &
          '] lon=[', active_lon_min, ',', active_lon_max, ']'
        call log_event(trim(nudge_region_msg), LOG_LEVEL_INFO)
      end if

      if (select_theta) then
        call log_event('ENNUF inference running: applying theta nudging', LOG_LEVEL_INFO)
        theta_residual = real(theta_us2_scaled(1, :) - theta_global_scaled(1, :), real64)
      end if
      if (select_qt) then
        call log_event('ENNUF inference running: applying qt nudging', LOG_LEVEL_INFO)
        qt_residual = real(qt_us2_scaled(1, :) - qt_global_scaled(1, :), real64)
      end if
      if (select_u) then
        call log_event('ENNUF inference running: applying u nudging', LOG_LEVEL_INFO)
        u_residual = real(u_us2_scaled(1, :) - u_global_scaled(1, :), real64)
      end if
      if (select_v) then
        call log_event('ENNUF inference running: applying v nudging', LOG_LEVEL_INFO)
        v_residual = real(v_us2_scaled(1, :) - v_global_scaled(1, :), real64)
      end if
    else
      theta_residual = 0.0_real64
      qt_residual = 0.0_real64
      u_residual = 0.0_real64
      v_residual = 0.0_real64
    end if
  end subroutine predict_column_residuals

end module gwd_ennuf_adapter_mod
