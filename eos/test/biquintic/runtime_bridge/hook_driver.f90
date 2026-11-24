program hook_driver
   use iso_fortran_env, only: input_unit, output_unit, iostat_end
   use bridge_support
   use eos_lib, only: eos_ptr, free_eos_handle, eos_shutdown
   implicit none
   type(EoS_General_Info), pointer :: controls
   character(len=4096) :: mesa_directory, inlist, massfile, record
   character(len=32) :: scenario, route
   real(dp) :: log_rho, log_temp, res(nv), density_partials(nv), temperature_partials(nv)
   real(dp) :: abundance_partials(num_eos_d_dxa_results,pinned_species)
   integer :: handle, ierr, io_status, stream, calls, finite_hybrid, row_index

   open(newunit=stream,file='/dev/stdout',action='write')
   open(unit=output_unit,file='/dev/stderr',action='write')
   if (command_argument_count() /= 3) error stop 'usage: hook_driver MESA INLIST MASSFILE'
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,inlist)
   call get_command_argument(3,massfile)
   call initialize_eos(trim(mesa_directory),trim(massfile),ierr)
   if (ierr /= 0) error stop 'EOS initialization failed'
   handle=alloc_eos_handle_using_inlist(trim(inlist),ierr)
   if (ierr /= 0) error stop 'handle allocation failed'
   call eos_ptr(handle,controls,ierr)
   if (ierr /= 0 .or. controls%use_WDM_AV) error stop 'invalid EOS handle'
   row_index=0
   do
      read(input_unit,'(a)',iostat=io_status) record
      if (io_status == iostat_end) exit
      if (io_status /= 0) error stop 'stdin failed'
      if (len_trim(record) == 0) cycle
      read(record,*,iostat=io_status) scenario, route, log_rho, log_temp
      if (io_status /= 0) error stop 'bad query'
      controls%use_other_FreeEOS_results = trim(scenario) /= 'off'
      controls%okay_to_convert_ierr_to_skip = trim(scenario) == 'unsafe' .or. trim(scenario) == 'off'
      controls%other_FreeEOS_results => diagnostic_hook
      if (trim(scenario) == 'missing') nullify(controls%other_FreeEOS_results)
      calls=0
      finite_hybrid=1
      call query_eos(trim(route),handle,log_rho,log_temp,res,density_partials,temperature_partials,abundance_partials,ierr)
      row_index=row_index+1
      write(stream,'(4(i0,1x),*(es25.16e3,1x))') row_index,ierr,calls,finite_hybrid, &
         res,density_partials,temperature_partials,abundance_partials
      flush(stream)
   end do
   call free_eos_handle(handle)
   call eos_shutdown()
contains
   subroutine diagnostic_hook(eos_handle,metals,hydrogen,abar,zbar,species,chem_id,net_iso,xa, &
         rho,logRho,temperature,logT,values,deriv_rho,deriv_temp,deriv_xa,status)
      integer, intent(in) :: eos_handle,species
      integer, pointer :: chem_id(:),net_iso(:)
      real(dp), intent(in) :: metals,hydrogen,abar,zbar,xa(:),rho,logRho,temperature,logT
      real(dp), intent(inout) :: values(:),deriv_rho(:),deriv_temp(:),deriv_xa(:,:)
      integer, intent(out) :: status
      calls=calls+1
      status=0
      if (.not. all(ieee_is_finite(deriv_xa))) then
         finite_hybrid=0
         status=-7998
         return
      end if
      if (trim(scenario) == 'failure') then
         status=-7999
      else if (trim(scenario) == 'inject') then
         values(i_Cv)=2.0_dp*values(i_Cv)
         deriv_rho(i_Cv)=2.0_dp*deriv_rho(i_Cv)
         deriv_temp(i_Cv)=2.0_dp*deriv_temp(i_Cv)
      end if
   end subroutine diagnostic_hook
end program hook_driver
