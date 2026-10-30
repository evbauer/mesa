program candidate_driver
   use iso_fortran_env, only: input_unit, output_unit, iostat_end
   use bridge_support
   use eos_lib, only: eos_ptr, free_eos_handle, eos_shutdown
   use freeeos_runtime_adapter, only: runtime_adapter, runtime_identity, runtime_initialize, runtime_query
   implicit none
   type(EoS_General_Info), pointer :: controls
   type(runtime_adapter) :: adapter
   type(runtime_identity) :: identity
   character(len=4096) :: mesa_directory, inlist, massfile, cache_path, identity_path, record
   character(len=32) :: scenario, route
   real(dp) :: log_rho, log_temp, res(nv), density_partials(nv), temperature_partials(nv)
   real(dp) :: abundance_partials(num_eos_d_dxa_results,pinned_species), full_abundance(nv,pinned_species)
   real(dp) :: saved_res(nv), saved_rho(nv), saved_temp(nv), saved_composition(4), altered_xa(pinned_species)
   integer :: handle, ierr, io_status, stream, calls, retained_hybrid, row_index, identity_unit

   open(newunit=stream,file='/dev/stdout',action='write')
   open(unit=output_unit,file='/dev/stderr',action='write')
   if (command_argument_count() /= 5) error stop 'usage: candidate_driver MESA INLIST MASSFILE CACHE IDENTITY'
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,inlist)
   call get_command_argument(3,massfile)
   call get_command_argument(4,cache_path)
   call get_command_argument(5,identity_path)
   call initialize_eos(trim(mesa_directory),trim(massfile),ierr)
   if (ierr /= 0) error stop 'EOS initialization failed'
   open(newunit=identity_unit,file=trim(identity_path),status='old',action='read')
   read(identity_unit,*) identity%cache_sha256
   read(identity_unit,*) identity%table%nq,identity%table%nt
   read(identity_unit,*) identity%table%hydrogen,identity%table%metals,identity%table%crad,identity%table%split
   read(identity_unit,*) identity%table%source_sha256
   read(identity_unit,*) identity%table%preprocessing_sha256,identity%table%radius,identity%table%route
   read(identity_unit,*) identity%table%u_min,identity%table%du,identity%table%tau_min,identity%table%dtau
   read(identity_unit,*) identity%minimum_logT,identity%pattern_sha256
   close(identity_unit)
   call runtime_initialize(adapter,trim(cache_path),identity,trim(massfile),ierr)
   write(output_unit,'(a,1x,i0)') 'ADAPTER_INIT',ierr
   if (ierr /= 0) stop 2
   handle=alloc_eos_handle_using_inlist(trim(inlist),ierr)
   if (ierr /= 0) error stop 'handle allocation failed'
   call eos_ptr(handle,controls,ierr)
   if (ierr /= 0 .or. controls%use_WDM_AV) error stop 'invalid EOS handle'
   controls%okay_to_convert_ierr_to_skip=.false.
   controls%other_FreeEOS_results=>potential_hook
   row_index=0
   do
      read(input_unit,'(a)',iostat=io_status) record
      if (io_status == iostat_end) exit
      if (io_status /= 0) error stop 'stdin failed'
      if (len_trim(record) == 0) cycle
      read(record,*,iostat=io_status) scenario,route,log_rho,log_temp
      if (io_status /= 0) error stop 'bad query'
      if (scenario /= 'off' .and. scenario /= 'candidate' .and. scenario /= 'direct' .and. &
          scenario /= 'composition_fault' .and. scenario /= 'input_fault') error stop 'unknown scenario'
      if (scenario == 'direct' .and. route /= 'component') error stop 'direct mode requires component route'
      controls%use_other_FreeEOS_results=scenario /= 'off'
      calls=0
      retained_hybrid=1
      call query_eos(trim(route),handle,log_rho,log_temp,res,density_partials,temperature_partials,abundance_partials,ierr)
      if (scenario == 'direct' .and. ierr == 0) then
         if (calls /= 1) error stop 'missing captured component'
         call runtime_query(adapter,pinned_species,pinned_chem_id,pinned_net_iso,pinned_xa, &
            saved_composition(1),saved_composition(2),saved_composition(3),saved_composition(4), &
            10.0_dp**log_rho,log_rho,10.0_dp**log_temp,log_temp, &
            saved_res,saved_rho,saved_temp,full_abundance,ierr)
         res=saved_res
         density_partials=saved_rho
         temperature_partials=saved_temp
         abundance_partials=full_abundance(1:num_eos_d_dxa_results,:)
      end if
      row_index=row_index+1
      write(stream,'(4(i0,1x),*(es25.16e3,1x))') row_index,ierr,calls,retained_hybrid, &
         res,density_partials,temperature_partials,abundance_partials
      flush(stream)
   end do
   call free_eos_handle(handle)
   call eos_shutdown()
contains
   subroutine potential_hook(eos_handle,metals,hydrogen,abar,zbar,species,chem_id,net_iso,xa, &
         rho,logRho,temperature,logT,values,deriv_rho,deriv_temp,deriv_xa,status)
      integer, intent(in) :: eos_handle,species
      integer, pointer :: chem_id(:),net_iso(:)
      real(dp), intent(in) :: metals,hydrogen,abar,zbar,xa(:),rho,logRho,temperature,logT
      real(dp), intent(inout) :: values(:),deriv_rho(:),deriv_temp(:),deriv_xa(:,:)
      integer, intent(out) :: status
      real(dp) :: trial_rho
      calls=calls+1
      full_abundance=deriv_xa
      saved_res=values
      saved_rho=deriv_rho
      saved_temp=deriv_temp
      saved_composition=[hydrogen,metals,abar,zbar]
      status=0
      if (scenario == 'direct') return
      altered_xa=xa
      trial_rho=rho
      if (scenario == 'composition_fault') then
         altered_xa(3)=xa(4)
         altered_xa(4)=xa(3)
      end if
      if (scenario == 'input_fault') trial_rho=rho*1.001_dp
      call runtime_query(adapter,species,chem_id,net_iso,altered_xa,hydrogen,metals,abar,zbar, &
         trial_rho,logRho,temperature,logT,values,deriv_rho,deriv_temp,deriv_xa,status)
      if (any(deriv_xa /= full_abundance) .or. .not. all(ieee_is_finite(deriv_xa))) retained_hybrid=0
      if (any(values(i_mu:i_eta) /= saved_res(i_mu:i_eta))) retained_hybrid=0
      if (any(deriv_rho(i_mu:i_eta) /= saved_rho(i_mu:i_eta))) retained_hybrid=0
      if (any(deriv_temp(i_mu:i_eta) /= saved_temp(i_mu:i_eta))) retained_hybrid=0
   end subroutine potential_hook
end program candidate_driver
