module run_star_extras
   use star_lib
   use star_def
   use const_def, only: dp, ln10
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_get_iso_id
   use eos_def, only: nv, i_frac_HELM, i_frac_FreeEOS, i_frac_WDM_AV, num_eos_d_dxa_results
   use eos_def, only: i_mu, i_eta
   use eos_lib, only: eosDT_get
   use freeeos_runtime_adapter, only: runtime_adapter, runtime_identity, runtime_initialize, runtime_query, &
      runtime_hybrid_policy
   use runtime_sha256, only: file_sha256
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: pin_count=20
   real(dp), parameter :: pin_tolerance=5d-13
   integer :: pin_ids(pin_count), query_unit, event_unit, matrix_unit, solver_unit
   integer :: query_count=0, model_tag=0, attempt_tag=0, matrix_count=0
   integer :: active_star_id=0, abundance_mode=0
   integer :: replayed_attempt=0
   type(runtime_adapter) :: adapter
   integer :: component_unit, component_calls=0, component_errors=0, hybrid_checks=0
   logical :: candidate_enabled=.false.
   logical :: force_retry=.false., retry_injected=.false.
   real(dp) :: pin_xa(pin_count)
   character(len=20) :: query_phase='initialization'
contains
   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: model
      character(len=32) :: isotope, mode
      integer :: pin_unit, position
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      active_star_id=id
      call get_environment_variable('TASK052_BACKEND', mode)
      if (trim(mode) /= 'BASELINE' .and. trim(mode) /= 'CANDIDATE') error stop 'TASK052 unknown backend'
      candidate_enabled=trim(mode) == 'CANDIDATE'
      abundance_mode=0
      call get_environment_variable('TASK052_FORCE_RETRY', mode)
      force_retry=trim(mode) == '1'
      open(newunit=pin_unit, file='uniform_xa.data', status='old', action='read')
      do position=1,pin_count
         read(pin_unit, *) isotope, pin_xa(position)
         pin_ids(position)=chem_get_iso_id(trim(isotope))
         if (pin_ids(position) <= 0) error stop 'missing pinned isotope'
      end do
      close(pin_unit)
      call initialize_backend()
      open(newunit=component_unit,file='component_queries.tsv',status='new',action='write')
      write(component_unit,'(a)') '# call model attempt phase logRho logT X Z abar zbar status hybrid_ok xa_in_network_order'
      open(newunit=query_unit, file='queries.tsv', status='new', action='write')
      open(newunit=event_unit, file='events.tsv', status='new', action='write')
      open(newunit=solver_unit, file='iterations.tsv', status='new', action='write')
      open(newunit=matrix_unit, file='matrices.bin', status='new', access='stream', form='unformatted')
      write(query_unit,'(a)') '# query model attempt phase logRho logT HELM OPAL_SCVH FreeEOS PC Skye CMS ideal WDM_AV max_relative_composition_error'
      write(event_unit,'(a)') '# phase model attempt retries age_yr dt_seconds solver_iterations max_relative_composition_error'
      write(solver_unit,'(a)') '# model attempt iter passed correction_norm max_correction residual_norm max_residual'
      model% eos_rq% okay_to_convert_ierr_to_skip=.false.
      model% eos_rq% use_other_FreeEOS_results=candidate_enabled
      model% eos_rq% other_FreeEOS_results => candidate_results
      model% eos_rq% other_eos_results => guarded_legacy_results
      model% extras_startup => extras_startup
      model% extras_start_step => extras_start_step
      model% extras_check_model => extras_check_model
      model% extras_finish_step => extras_finish_step
      model% extras_after_evolve => extras_after_evolve
      model% other_before_struct_burn_mix => before_attempt
      model% other_solver_monitor => solver_monitor
      model% other_after_solver_setmatrix => audit_matrix
   end subroutine extras_controls


   subroutine initialize_backend()
      type(runtime_identity) :: identity
      character(len=4096) :: cache_path, pattern_path, header
      character(len=64) :: identity_sha
      character(len=8) :: element
      integer :: identity_unit, pattern_unit, position, ierr, backend_unit
      real(dp) :: metal_sum, original_pin(pin_count), reference_error
      original_pin=pin_xa
      call get_environment_variable('TASK052_CACHE',cache_path)
      call get_environment_variable('TASK052_PATTERN',pattern_path)
      call file_sha256('solar.identity',identity_sha,ierr)
      if (ierr /= 0 .or. identity_sha /= &
          'ac3eef30bba5389bb8c04f9adf56a22ce053e4df47260f27e0bb00fd63f03ea3') &
         error stop 'TASK052 identity fixture authentication failed'
      open(newunit=identity_unit,file='solar.identity',status='old',action='read')
      read(identity_unit,*) identity%cache_sha256
      read(identity_unit,*) identity%table%nq,identity%table%nt
      read(identity_unit,*) identity%table%hydrogen,identity%table%metals,identity%table%crad,identity%table%split
      read(identity_unit,*) identity%table%source_sha256
      read(identity_unit,*) identity%table%preprocessing_sha256,identity%table%radius,identity%table%route
      read(identity_unit,*) identity%table%u_min,identity%table%du,identity%table%tau_min,identity%table%dtau
      read(identity_unit,*) identity%minimum_logT,identity%pattern_sha256
      close(identity_unit)
      call runtime_initialize(adapter,trim(cache_path),identity,trim(pattern_path),ierr)
      if (ierr /= 0) then
         write(*,*) 'TASK052_ADAPTER_INIT_ERROR',ierr
         error stop 'TASK052 adapter initialization failed'
      end if
      open(newunit=pattern_unit,file=trim(pattern_path),status='old',action='read')
      read(pattern_unit,'(a)') header
      pin_xa(1:2)=[0.7_dp,0.28_dp]
      do position=3,pin_count
         read(pattern_unit,*) element,pin_xa(position)
      end do
      close(pattern_unit)
      metal_sum=sum(pin_xa(3:))
      pin_xa(3:)=0.02_dp*pin_xa(3:)/metal_sum
      reference_error=maxval(abs(original_pin-pin_xa)/pin_xa)
      if (reference_error > pin_tolerance) error stop 'TASK052 seed pattern differs from047'
      open(newunit=backend_unit,file='backend.txt',status='new',action='write')
      write(backend_unit,*) 'candidate_enabled',candidate_enabled
      write(backend_unit,*) 'cache_sha256',identity%cache_sha256
      write(backend_unit,*) 'identity_sha256',identity_sha
      write(backend_unit,*) 'pattern_sha256',identity%pattern_sha256
      write(backend_unit,*) 'minimum_logT',identity%minimum_logT
      write(backend_unit,*) 'seed_reference_relative_error',reference_error
      write(backend_unit,*) 'hybrid_policy',runtime_hybrid_policy
      write(backend_unit,*) 'pin_ids',pin_ids
      write(backend_unit,*) 'pin_xa',pin_xa
      close(backend_unit)
   end subroutine initialize_backend

   subroutine candidate_results(handle,metals,hydrogen,abar,zbar,species,chem_id,net_iso,xa, &
         rho,logRho,temperature,logT,values,deriv_rho,deriv_temp,deriv_xa,status)
      integer, intent(in) :: handle,species
      integer, pointer :: chem_id(:),net_iso(:)
      real(dp), intent(in) :: metals,hydrogen,abar,zbar,xa(:),rho,logRho,temperature,logT
      real(dp), intent(inout) :: values(:),deriv_rho(:),deriv_temp(:),deriv_xa(:,:)
      integer, intent(out) :: status
      real(dp) :: legacy_abundance(size(deriv_xa,1),size(deriv_xa,2))
      real(dp) :: legacy_aux(3,3)
      logical :: hybrid_ok
      legacy_abundance=deriv_xa
      legacy_aux(:,1)=values(i_mu:i_eta)
      legacy_aux(:,2)=deriv_rho(i_mu:i_eta)
      legacy_aux(:,3)=deriv_temp(i_mu:i_eta)
      call runtime_query(adapter,species,chem_id,net_iso,xa,hydrogen,metals,abar,zbar, &
         rho,logRho,temperature,logT,values,deriv_rho,deriv_temp,deriv_xa,status)
      hybrid_ok=.false.
      if (status == 0) then
         hybrid_ok=all(deriv_xa == legacy_abundance) .and. all(ieee_is_finite(deriv_xa))
         hybrid_ok=hybrid_ok .and. all(values(i_mu:i_eta) == legacy_aux(:,1))
         hybrid_ok=hybrid_ok .and. all(deriv_rho(i_mu:i_eta) == legacy_aux(:,2))
         hybrid_ok=hybrid_ok .and. all(deriv_temp(i_mu:i_eta) == legacy_aux(:,3))
      end if
      !$omp critical(task052_components)
      component_calls=component_calls+1
      if (status /= 0) component_errors=component_errors+1
      if (hybrid_ok) hybrid_checks=hybrid_checks+1
      write(component_unit,*) component_calls,model_tag,attempt_tag,trim(query_phase), &
         logRho,logT,hydrogen,metals,abar,zbar,status,hybrid_ok,xa
      if (status /= 0 .or. .not. hybrid_ok) then
         flush(component_unit)
         flush(query_unit)
         flush(event_unit)
      end if
      !$omp end critical(task052_components)
      if (status /= 0) error stop 'TASK052 candidate query rejected; see component_queries.tsv'
      if (.not. hybrid_ok) error stop 'TASK052 legacy auxiliary/abundance policy violated'
   end subroutine candidate_results

   subroutine check_vector(species, chem_id, net_iso, xa, difference)
      integer, intent(in) :: species, chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(out) :: difference
      integer :: position, mapped
      if (species /= pin_count .or. size(chem_id) /= pin_count .or. size(xa) /= pin_count) &
         error stop 'TASK052 wrong species count'
      if (size(net_iso) /= num_chem_isos) error stop 'TASK052 wrong mapping length'
      if (.not. all(ieee_is_finite(xa)) .or. any(xa <= 0)) error stop 'TASK052 invalid abundances'
      if (abs(sum(xa)-1d0) > pin_tolerance) error stop 'TASK052 nonnormalized abundances'
      if (any(net_iso < 0) .or. any(net_iso > species)) error stop 'TASK052 invalid mapping'
      if (count(net_iso /= 0) /= pin_count) error stop 'TASK052 extra or missing isotope'
      difference=0
      do position=1,pin_count
         mapped=net_iso(pin_ids(position))
         if (mapped < 1 .or. mapped > species) error stop 'TASK052 isotope absent'
         if (chem_id(mapped) /= pin_ids(position)) error stop 'TASK052 inconsistent isotope mapping'
         difference=max(difference,abs(xa(mapped)-pin_xa(position))/pin_xa(position))
      end do
      if (difference > pin_tolerance) error stop 'TASK052 composition drift'
   end subroutine check_vector

   subroutine guarded_legacy_results(handle, species, chem_id, net_iso, xa, rho, logRho, temperature, logT, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)
      integer, intent(in) :: handle, species
      integer, pointer :: chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:), rho, logRho, temperature, logT
      real(dp), intent(inout) :: res(:), d_dlnd(:), d_dlnT(:), d_dxa(:,:)
      integer, intent(out) :: ierr
      real(dp) :: difference
      integer :: species_index, row
      ierr=0
      call check_vector(species, chem_id, net_iso, xa, difference)
      if (.not. all(ieee_is_finite(res))) error stop 'TASK052 nonfinite EOS result'
      if (res(i_frac_WDM_AV) /= 0) error stop 'TASK052 WDM_AV active'
      if (abs(sum(res(i_frac_HELM:i_frac_WDM_AV))-1d0) > 1d-12) error stop 'TASK052 invalid route fractions'
      if (size(d_dxa,1) /= nv) error stop 'TASK052 incomplete abundance rows'
      if (abundance_mode /= 0) then
         do species_index=1,species
            do row=1,size(d_dxa,1)
               d_dxa(row,species_index)=abundance_mode*1d6*(row+species_index)
            end do
         end do
      end if
      !$omp critical(task049_queries)
      query_count=query_count+1
      write(query_unit,*) query_count, model_tag, attempt_tag, trim(query_phase), logRho, logT, &
         res(i_frac_HELM:i_frac_WDM_AV), difference
      !$omp end critical(task049_queries)
   end subroutine guarded_legacy_results

   subroutine check_model_state(model, phase, previous)
      type(star_info), pointer :: model
      character(len=*), intent(in) :: phase
      logical, intent(in) :: previous
      integer :: zone
      real(dp) :: difference, maximum
      if (model% dxdt_nuc_factor /= 0 .or. model% mix_factor /= 0 .or. model% eps_nuc_factor /= 0) &
         error stop 'TASK052 burning or mixing enabled'
      if (model% do_element_diffusion .or. model% do_conv_premix .or. model% okay_to_remesh) &
         error stop 'TASK052 abundance transport enabled'
      if (model% rotation_flag .or. model% mass_change /= 0) error stop 'TASK052 mass/rotation change enabled'
      if (model% include_composition_in_eps_grav .or. model% fix_d_eos_dxa_partials) &
         error stop 'TASK052 composition EOS terms enabled'
      if (trim(model% energy_eqn_option) /= 'eps_grav' .or. .not. model% use_time_centered_eps_grav) &
         error stop 'TASK052 energy or time centering changed'
      maximum=0
      do zone=1,model% nz
         call check_vector(model% species, model% chem_id, model% net_iso, model% xa(:,zone), difference)
         maximum=max(maximum,difference)
         if (previous) then
            call check_vector(model% species, model% chem_id, model% net_iso, model% xa_start(:,zone), difference)
            maximum=max(maximum,difference)
         end if
      end do
      write(event_unit,*) trim(phase), model% model_number, attempt_tag, model% num_retries, &
         model% star_age, model% dt, model% num_solver_iterations, maximum
      flush(event_unit)
   end subroutine check_model_state

   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type(star_info), pointer :: model
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      call check_model_state(model, 'startup', .false.)
   end subroutine extras_startup

   integer function extras_start_step(id) result(status)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      model_tag=model% model_number
      query_phase='step_setup'
      call check_model_state(model, 'start_step', .false.)
      status=keep_going
   end function extras_start_step

   subroutine before_attempt(id, dt, status)
      integer, intent(in) :: id
      real(dp), intent(in) :: dt
      integer, intent(out) :: status
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      attempt_tag=attempt_tag+1
      model_tag=model% model_number
      call check_model_state(model, 'attempt', .false.)
      query_phase='trial_or_retry'
      status=keep_going
   end subroutine before_attempt

   subroutine replay_previous(model)
      type(star_info), pointer :: model
      integer :: zone, ierr
      real(dp) :: res(nv), density_partials(nv), temperature_partials(nv)
      real(dp) :: abundance_partials(num_eos_d_dxa_results,pin_count)
      query_phase='previous_replay'
      do zone=1,model% nz
         call eosDT_get(model% eos_handle, model% species, model% chem_id, model% net_iso, model% xa_start(:,zone), &
            model% rho_start(zone), log10(model% rho_start(zone)), &
            exp(model% lnT_start(zone)), model% lnT_start(zone)/ln10, &
            res, density_partials, temperature_partials, abundance_partials, ierr)
         if (ierr /= 0) error stop 'TASK052 previous-state EOS query failed'
      end do
      query_phase='trial_or_retry'
      replayed_attempt=attempt_tag
   end subroutine replay_previous

   subroutine audit_matrix(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer :: hydro_count, zones
      type(star_info), pointer :: model
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      if (model% do_burn .or. model% do_mix) error stop 'TASK052 abundance equations active'
      hydro_count=model% nvar_hydro
      zones=model% nz
      if (size(model% dblk,1) /= hydro_count) error stop 'TASK052 solver includes nonstructural variables'
      call check_model_state(model, 'matrix', .true.)
      if (replayed_attempt /= attempt_tag) call replay_previous(model)
      matrix_count=matrix_count+1
      write(matrix_unit) model% model_number, attempt_tag, hydro_count, zones
      write(matrix_unit) model% equ(1:hydro_count,1:zones)
      write(matrix_unit) model% dblk(1:hydro_count,1:hydro_count,1:zones)
      write(matrix_unit) model% lblk(1:hydro_count,1:hydro_count,1:zones)
      write(matrix_unit) model% ublk(1:hydro_count,1:hydro_count,1:zones)
      flush(matrix_unit)
   end subroutine audit_matrix

   subroutine solver_monitor(id, iter, passed_tol_tests, correction_norm, max_correction, residual_norm, max_residual, ierr)
      integer, intent(in) :: id, iter
      logical, intent(in) :: passed_tol_tests
      real(dp), intent(in) :: correction_norm, max_correction, residual_norm, max_residual
      integer, intent(out) :: ierr
      ierr=0
      write(solver_unit,*) model_tag, attempt_tag, iter, passed_tol_tests, &
         correction_norm, max_correction, residual_norm, max_residual
      flush(solver_unit)
   end subroutine solver_monitor

   integer function extras_check_model(id) result(status)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      call check_model_state(model, 'check_model', .true.)
      status=keep_going
      if (force_retry .and. .not. retry_injected .and. model% model_number == 2) then
         retry_injected=.true.
         model% retry_message='TASK052 deliberate retry audit'
         call check_model_state(model, 'forced_retry', .true.)
         status=retry
      end if
   end function extras_check_model

   integer function extras_finish_step(id) result(status)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      call check_model_state(model, 'accepted', .true.)
      model% need_to_save_profiles_now=.true.
      model% need_to_update_history_now=.true.
      flush(query_unit)
      status=keep_going
   end function extras_finish_step

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer :: summary_unit
      type(star_info), pointer :: model
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      call check_model_state(model, 'after_evolve', .false.)
      open(newunit=summary_unit,file='completion.txt',status='new')
      write(summary_unit,*) 'model_number',model% model_number
      write(summary_unit,*) 'termination_code',model% termination_code
      write(summary_unit,*) 'queries',query_count
      write(summary_unit,*) 'attempts',attempt_tag
      write(summary_unit,*) 'matrix_evaluations',matrix_count
      write(summary_unit,*) 'component_calls',component_calls
      write(summary_unit,*) 'component_errors',component_errors
      write(summary_unit,*) 'hybrid_checks',hybrid_checks
      close(summary_unit)
      close(query_unit)
      close(event_unit)
      close(solver_unit)
      close(matrix_unit)
      close(component_unit)
   end subroutine extras_after_evolve
end module run_star_extras
