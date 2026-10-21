!-----------------------------------------------------------------------------
! ENNUF adapter for column residual prediction in gravity_wave_drag.
! This is a small integration layer around the generated ENNUF model modules.
! Inputs match train_hecate_ants.py: four interleaved atmospheric channels,
! mean orography, signed x/y slopes, orography std, and cyclic time features.
! Return physical US2-minus-global tendencies after separate output unscaling.
! Qt is an input only; only theta and wind residuals are predicted.
! Only enabled targets are evaluated, and only within the selected region.
! The algorithm resets and reports rank-local inference counts around each drag call.
!-----------------------------------------------------------------------------
module gwd_ennuf_adapter_mod
  use, intrinsic :: iso_fortran_env, only: real32, real64

  use log_mod, only: log_event, LOG_LEVEL_INFO
  use gwd_ennuf_scaling_mod, only: ennuf_standardise_input, ennuf_destandardise_output, &
                                 ennuf_n_levels, ennuf_n_input

  use hecate_global_theta_dummy_mod, only: hecate_global_theta_dummy
  use hecate_us2_theta_dummy_mod, only: hecate_us2_theta_dummy
  use hecate_global_x_wind_dummy_mod, only: hecate_global_x_wind_dummy
  use hecate_us2_x_wind_dummy_mod, only: hecate_us2_x_wind_dummy
  use hecate_global_y_wind_dummy_mod, only: hecate_global_y_wind_dummy
  use hecate_us2_y_wind_dummy_mod, only: hecate_us2_y_wind_dummy

  implicit none

  private
  public :: thermodynamic_correction, build_ennuf_input, predict_column_residuals
  public :: set_nudge_options
  public :: reset_ennuf_inference_counts, log_ennuf_inference_summary
  public :: log_hecate_start, log_hecate_nudging_done

  integer, parameter :: ml_profile_len = ennuf_n_levels
  integer, parameter :: ml_input_len = ennuf_n_input
  real(real64), parameter :: cp = 1005.0_real64
  real(real64), parameter :: lc = 2.501e6_real64
  real(real64), parameter :: lf = 0.334e6_real64
  real(real64), parameter :: two_pi = 6.28318530717958647692528676655900577_real64

  logical, save :: nudge_apply_globally = .true.
  logical, save :: nudge_use_theta = .true.
  logical, save :: nudge_use_u = .true.
  logical, save :: nudge_use_v = .true.
  real(real64), save :: nudge_lat_min = -90.0_real64
  real(real64), save :: nudge_lat_max = 90.0_real64
  real(real64), save :: nudge_lon_min = -180.0_real64
  real(real64), save :: nudge_lon_max = 180.0_real64
  integer, save :: inference_columns(3) = 0

contains

  subroutine log_hecate_start()
    call log_event('HECATE starting: nudging switches theta=' // merge('on ', 'off', nudge_use_theta) // &
      ' qt=input-only' // &
      ' u=' // merge('on ', 'off', nudge_use_u) // &
      ' v=' // merge('on ', 'off', nudge_use_v), LOG_LEVEL_INFO)
  end subroutine log_hecate_start

  subroutine log_hecate_nudging_done()
    if (any(inference_columns > 0)) then
      call log_event('HECATE nudging done: local corrections added to theta/wind increments; moisture unchanged', &
                     LOG_LEVEL_INFO)
    else
      call log_event('HECATE nudging done: no local columns ran inference; no HECATE corrections applied', &
                     LOG_LEVEL_INFO)
    end if
  end subroutine log_hecate_nudging_done

  subroutine reset_ennuf_inference_counts()
    inference_columns = 0
  end subroutine reset_ennuf_inference_counts

  subroutine log_ennuf_inference_summary(seconds_of_day)
    real(real64), intent(in) :: seconds_of_day
    character(len=512) :: summary

    ! Called outside the OpenMP column loop, after all local segments finish.
    write(summary, '(A,F12.1,A,L1,A,3(L1,1X),A,4(F9.2,1X),A,3(I0,1X))') &
      'HECATE inference complete: time_seconds=', seconds_of_day, &
      ' global=', nudge_apply_globally, &
      ' enabled(theta u v)=', nudge_use_theta, nudge_use_u, nudge_use_v, &
      ' bounds(lat_min lat_max lon_min lon_max)=', &
      nudge_lat_min, nudge_lat_max, nudge_lon_min, nudge_lon_max, &
      ' columns with inference/residuals (theta u v)=', inference_columns
    call log_event(trim(summary), LOG_LEVEL_INFO)
  end subroutine log_ennuf_inference_summary

  subroutine set_nudge_options(apply_globally, use_theta, use_u, use_v, &
                               lat_min, lat_max, lon_min, lon_max)
    logical, intent(in), optional :: apply_globally
    logical, intent(in), optional :: use_theta, use_u, use_v
    real(real64), intent(in), optional :: lat_min, lat_max, lon_min, lon_max

    if (present(apply_globally)) nudge_apply_globally = apply_globally
    if (present(use_theta)) nudge_use_theta = use_theta
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

    theta_li = theta - ((lc / cp) * qcl + ((lc + lf) / cp) * qcf)
  end function thermodynamic_correction

  pure function compute_qt(qv, qcl, qcf) result(qt)
    real(real64), intent(in) :: qv(ml_profile_len)
    real(real64), intent(in) :: qcl(ml_profile_len)
    real(real64), intent(in) :: qcf(ml_profile_len)
    real(real64) :: qt(ml_profile_len)

    qt = qv + qcl + qcf
  end function compute_qt

  subroutine build_ennuf_input(theta, qt, u, v, orography, slope_x, slope_y, orog_std, &
                              seconds_of_day, combined_input)
    real(real64), intent(in) :: theta(ml_profile_len)
    real(real64), intent(in) :: qt(ml_profile_len)
    real(real64), intent(in) :: u(ml_profile_len)
    real(real64), intent(in) :: v(ml_profile_len)
    real(real64), intent(in) :: orography, slope_x, slope_y, orog_std
    real(real64), intent(in) :: seconds_of_day
    real(real32), intent(out) :: combined_input(ml_input_len)

    real(real64) :: time_angle

    time_angle = two_pi * modulo(seconds_of_day, 86400.0_real64) / 86400.0_real64

    ! Match np.stack(channels, axis=-1).reshape(...) in training.
    combined_input(1:4*ml_profile_len:4) = real(theta, real32)
    combined_input(2:4*ml_profile_len:4) = real(qt, real32)
    combined_input(3:4*ml_profile_len:4) = real(u, real32)
    combined_input(4:4*ml_profile_len:4) = real(v, real32)
    combined_input(4*ml_profile_len+1:ml_input_len) = &
      real([orography, slope_x, slope_y, orog_std, sin(time_angle), cos(time_angle)], real32)
  end subroutine build_ennuf_input

  subroutine predict_column_residuals(theta, qv, qcl, qcf, u, v, orography, &
                                      slope_x, slope_y, orog_std,                &
                                      seconds_of_day, theta_residual, &
                                      u_residual, v_residual,                   &
                                      use_theta, use_u, use_v,                  &
                                      apply_globally, latitude, longitude,      &
                                      lat_min, lat_max, lon_min, lon_max)
    real(real64), intent(in) :: theta(ml_profile_len)
    real(real64), intent(in) :: qv(ml_profile_len)
    real(real64), intent(in) :: qcl(ml_profile_len)
    real(real64), intent(in) :: qcf(ml_profile_len)
    real(real64), intent(in) :: u(ml_profile_len)
    real(real64), intent(in) :: v(ml_profile_len)
    real(real64), intent(in) :: orography, slope_x, slope_y, orog_std
    real(real64), intent(in) :: seconds_of_day
    real(real64), intent(out) :: theta_residual(ml_profile_len)
    real(real64), intent(out) :: u_residual(ml_profile_len)
    real(real64), intent(out) :: v_residual(ml_profile_len)
    logical, intent(in), optional :: use_theta, use_u, use_v
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
    real(real32) :: u_global_scaled(1, ml_profile_len)
    real(real32) :: u_us2_scaled(1, ml_profile_len)
    real(real32) :: v_global_scaled(1, ml_profile_len)
    real(real32) :: v_us2_scaled(1, ml_profile_len)
    real(real64) :: prediction_global(ml_profile_len), prediction_us2(ml_profile_len)
    logical :: allow_apply
    logical :: in_region
    logical :: select_theta
    logical :: select_u
    logical :: select_v
    real(real64) :: active_lat_min
    real(real64) :: active_lat_max
    real(real64) :: active_lon_min
    real(real64) :: active_lon_max

    theta_residual = 0.0_real64
    u_residual = 0.0_real64
    v_residual = 0.0_real64

    select_theta = nudge_use_theta
    select_u = nudge_use_u
    select_v = nudge_use_v
    if (present(use_theta)) select_theta = use_theta
    if (present(use_u)) select_u = use_u
    if (present(use_v)) select_v = use_v

    if (.not. (select_theta .or. select_u .or. select_v)) return

    allow_apply = nudge_apply_globally
    if (present(apply_globally)) allow_apply = apply_globally

    active_lat_min = nudge_lat_min
    active_lat_max = nudge_lat_max
    active_lon_min = nudge_lon_min
    active_lon_max = nudge_lon_max
    if (present(lat_min)) active_lat_min = lat_min
    if (present(lat_max)) active_lat_max = lat_max
    if (present(lon_min)) active_lon_min = lon_min
    if (present(lon_max)) active_lon_max = lon_max

    in_region = .true.
    if (.not. allow_apply) then
      in_region = .false.
      if (present(latitude) .and. present(longitude)) then
        if (latitude >= active_lat_min .and. latitude <= active_lat_max .and. &
            longitude >= active_lon_min .and. longitude <= active_lon_max) then
          in_region = .true.
        end if
      end if
    end if

    if (.not. in_region) return

    theta_li = thermodynamic_correction(theta, qcl, qcf)
    qt = compute_qt(qv, qcl, qcf)

    call build_ennuf_input(theta_li, qt, u, v, orography, slope_x, slope_y, orog_std, &
                          seconds_of_day, combined_input)
    combined_input_dbl = real(combined_input, real64)

    call ennuf_standardise_input(combined_input_dbl, 'global', combined_input_global_scaled)
    call ennuf_standardise_input(combined_input_dbl, 'us2', combined_input_us2_scaled)
    combined_input_global_for_model = real(combined_input_global_scaled, real32)
    combined_input_us2_for_model = real(combined_input_us2_scaled, real32)

    if (in_region) then
      theta_residual = 0.0_real64
      u_residual = 0.0_real64
      v_residual = 0.0_real64

      if (select_theta) then
        call hecate_global_theta_dummy(combined_input_global_for_model, theta_global_scaled)
        call hecate_us2_theta_dummy(combined_input_us2_for_model, theta_us2_scaled)
        call ennuf_destandardise_output(real(theta_global_scaled(1,:), real64), 'global', 'theta', prediction_global)
        call ennuf_destandardise_output(real(theta_us2_scaled(1,:), real64), 'us2', 'theta', prediction_us2)
        theta_residual = prediction_us2 - prediction_global
        !$omp atomic update
        inference_columns(1) = inference_columns(1) + 1
      end if
      if (select_u) then
        call hecate_global_x_wind_dummy(combined_input_global_for_model, u_global_scaled)
        call hecate_us2_x_wind_dummy(combined_input_us2_for_model, u_us2_scaled)
        call ennuf_destandardise_output(real(u_global_scaled(1,:), real64), 'global', 'x_wind', prediction_global)
        call ennuf_destandardise_output(real(u_us2_scaled(1,:), real64), 'us2', 'x_wind', prediction_us2)
        u_residual = prediction_us2 - prediction_global
        !$omp atomic update
        inference_columns(2) = inference_columns(2) + 1
      end if
      if (select_v) then
        call hecate_global_y_wind_dummy(combined_input_global_for_model, v_global_scaled)
        call hecate_us2_y_wind_dummy(combined_input_us2_for_model, v_us2_scaled)
        call ennuf_destandardise_output(real(v_global_scaled(1,:), real64), 'global', 'y_wind', prediction_global)
        call ennuf_destandardise_output(real(v_us2_scaled(1,:), real64), 'us2', 'y_wind', prediction_us2)
        v_residual = prediction_us2 - prediction_global
        !$omp atomic update
        inference_columns(3) = inference_columns(3) + 1
      end if
    else
      theta_residual = 0.0_real64
      u_residual = 0.0_real64
      v_residual = 0.0_real64
    end if
  end subroutine predict_column_residuals

end module gwd_ennuf_adapter_mod
