module run_star_extras
   use star_lib
   use star_def
   use const_def, only: dp
   use chem_lib, only: chem_get_iso_id
   use freeeos_family, only: freeeos_family_get_counters
   use iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer :: event_unit, iteration_unit, attempt_count=0, matrix_count=0
   integer :: hydrogen_id, helium_id, helium3_id
   integer :: model_tag=0, attempt_tag=0
   real(dp) :: initial_hydrogen=0, initial_helium=0, initial_helium3=0
   logical :: saw_burn=.false., saw_mix=.false., saw_abundance_equations=.false.
contains
   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: model
      call star_ptr(id, model, ierr)
      if (ierr /= 0) return
      hydrogen_id=chem_get_iso_id('h1')
      helium_id=chem_get_iso_id('he4')
      helium3_id=chem_get_iso_id('he3')
      open(newunit=event_unit,file='events.tsv',status='new',action='write')
      open(newunit=iteration_unit,file='iterations.tsv',status='new',action='write')
      write(event_unit,'(a)') '# phase model attempt retries age_yr dt_s zones species do_burn do_mix h1_center he4_center he3_center eps_nuc_center solver_iterations'
      write(iteration_unit,'(a)') '# model attempt iteration passed correction_norm max_correction residual_norm max_residual'
      model% extras_startup => extras_startup
      model% extras_start_step => extras_start_step
      model% extras_finish_step => extras_finish_step
      model% extras_after_evolve => extras_after_evolve
      model% other_before_struct_burn_mix => before_attempt
      model% other_solver_monitor => solver_monitor
      model% other_after_solver_setmatrix => audit_matrix
   end subroutine extras_controls

   subroutine isotope_positions(model, hydrogen, helium, helium3)
      type(star_info), pointer :: model
      integer, intent(out) :: hydrogen, helium, helium3
      hydrogen=model% net_iso(hydrogen_id)
      helium=model% net_iso(helium_id)
      helium3=model% net_iso(helium3_id)
      if (hydrogen <= 0 .or. helium <= 0 .or. helium3 <= 0) error stop 'TASK064 burning network lacks H/He isotopes'
   end subroutine isotope_positions

   subroutine record_state(model, phase)
      type(star_info), pointer :: model
      character(len=*), intent(in) :: phase
      integer :: hydrogen, helium, helium3
      call isotope_positions(model,hydrogen,helium,helium3)
      if (.not. all(ieee_is_finite(model% xa(:,model% nz)))) error stop 'TASK064 nonfinite composition'
      write(event_unit,*) trim(phase),model% model_number,attempt_tag,model% num_retries, &
         model% star_age,model% dt,model% nz,model% species,model% do_burn,model% do_mix, &
         model% xa(hydrogen,model% nz),model% xa(helium,model% nz), &
         model% xa(helium3,model% nz),model% eps_nuc(model% nz),model% num_solver_iterations
      flush(event_unit)
   end subroutine record_state

   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      integer :: hydrogen, helium, helium3
      type(star_info), pointer :: model
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      call isotope_positions(model,hydrogen,helium,helium3)
      initial_hydrogen=model% xa(hydrogen,model% nz)
      initial_helium=model% xa(helium,model% nz)
      initial_helium3=model% xa(helium3,model% nz)
      call record_state(model,'startup')
   end subroutine extras_startup

   integer function extras_start_step(id) result(status)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      model_tag=model% model_number
      call record_state(model,'start_step')
      status=keep_going
   end function extras_start_step

   subroutine before_attempt(id, dt, status)
      integer, intent(in) :: id
      real(dp), intent(in) :: dt
      integer, intent(out) :: status
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      attempt_count=attempt_count+1
      attempt_tag=attempt_count
      model_tag=model% model_number
      call record_state(model,'attempt')
      status=keep_going
   end subroutine before_attempt

   subroutine audit_matrix(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: model
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      matrix_count=matrix_count+1
      saw_burn=saw_burn .or. model% do_burn
      saw_mix=saw_mix .or. model% do_mix
      saw_abundance_equations=saw_abundance_equations .or. size(model% dblk,1) > model% nvar_hydro
      if (matrix_count == 1) then
         write(event_unit,*) 'matrix',model% model_number,attempt_tag,model% nvar_hydro, &
            size(model% dblk,1),model% do_burn,model% do_mix
         flush(event_unit)
      end if
   end subroutine audit_matrix

   subroutine solver_monitor(id, iter, passed_tol_tests, correction_norm, max_correction, residual_norm, max_residual, ierr)
      integer, intent(in) :: id,iter
      logical, intent(in) :: passed_tol_tests
      real(dp), intent(in) :: correction_norm,max_correction,residual_norm,max_residual
      integer, intent(out) :: ierr
      ierr=0
      write(iteration_unit,*) model_tag,attempt_tag,iter,passed_tol_tests, &
         correction_norm,max_correction,residual_norm,max_residual
      flush(iteration_unit)
   end subroutine solver_monitor

   integer function extras_finish_step(id) result(status)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: model
      status=terminate
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      call record_state(model,'accepted')
      model% need_to_save_profiles_now=.true.
      model% need_to_update_history_now=.true.
      status=keep_going
   end function extras_finish_step

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer :: summary_unit, hydrogen, helium, helium3
      integer(int64) :: native_successes, native_failures, native_loads
      type(star_info), pointer :: model
      call star_ptr(id,model,ierr)
      if (ierr /= 0) return
      call isotope_positions(model,hydrogen,helium,helium3)
      call freeeos_family_get_counters(native_successes,native_failures,native_loads)
      call record_state(model,'after_evolve')
      open(newunit=summary_unit,file='completion.txt',status='new',action='write')
      write(summary_unit,*) 'model_number',model% model_number
      write(summary_unit,*) 'termination_code',model% termination_code
      write(summary_unit,*) 'attempts',attempt_count
      write(summary_unit,*) 'matrix_evaluations',matrix_count
      write(summary_unit,*) 'saw_burn',saw_burn
      write(summary_unit,*) 'saw_mix',saw_mix
      write(summary_unit,*) 'saw_abundance_equations',saw_abundance_equations
      write(summary_unit,*) 'native_successes',native_successes
      write(summary_unit,*) 'native_failures',native_failures
      write(summary_unit,*) 'native_loads',native_loads
      write(summary_unit,*) 'delta_center_h1',model% xa(hydrogen,model% nz)-initial_hydrogen
      write(summary_unit,*) 'delta_center_he4',model% xa(helium,model% nz)-initial_helium
      write(summary_unit,*) 'delta_center_he3',model% xa(helium3,model% nz)-initial_helium3
      close(summary_unit)
      close(event_unit)
      close(iteration_unit)
   end subroutine extras_after_evolve
end module run_star_extras
