program native_driver
   use iso_fortran_env, only: input_unit, output_unit, iostat_end, int64
   use bridge_support
   use eos_lib, only: free_eos_handle, eos_shutdown
   use freeeos_family, only: freeeos_family_get_counters
   implicit none
   character(len=4096) :: mesa_directory, inlist, massfile, record
   character(len=32) :: route
   real(dp) :: log_rho, log_temp, hydrogen, metals, metal_pattern(pinned_species-2)
   real(dp) :: original_xa(pinned_species)
   real(dp) :: res(nv), density_partials(nv), temperature_partials(nv)
   real(dp) :: abundance_partials(num_eos_d_dxa_results,pinned_species)
   integer :: handle, ierr, io_status, stream, row_index
   integer(int64) :: successes, failures, loads

   open(newunit=stream,file='/dev/stdout',action='write')
   open(unit=output_unit,file='/dev/stderr',action='write')
   if (command_argument_count() /= 3) error stop 'usage: native_driver MESA INLIST MASSFILE'
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,inlist)
   call get_command_argument(3,massfile)
   call initialize_eos(trim(mesa_directory),trim(massfile),ierr)
   if (ierr /= 0) error stop 'EOS initialization failed'
   handle=alloc_eos_handle_using_inlist(trim(inlist),ierr)
   if (ierr /= 0) error stop 'handle allocation failed'
   metal_pattern=pinned_xa(3:)/sum(pinned_xa(3:))
   original_xa=pinned_xa
   row_index=0
   do
      read(input_unit,'(a)',iostat=io_status) record
      if (io_status == iostat_end) exit
      if (io_status /= 0) error stop 'stdin failed'
      if (len_trim(record) == 0) cycle
      read(record,*,iostat=io_status) route,log_rho,log_temp,hydrogen,metals
      if (io_status /= 0) error stop 'bad query'
      pinned_xa(1)=hydrogen
      pinned_xa(2)=1-hydrogen-metals
      pinned_xa(3:)=metals*metal_pattern
      if (hydrogen == .70_dp .and. metals == .02_dp) pinned_xa=original_xa
      call query_eos(trim(route),handle,log_rho,log_temp,res,density_partials,temperature_partials,abundance_partials,ierr)
      call freeeos_family_get_counters(successes,failures,loads)
      row_index=row_index+1
      write(stream,'(5(i0,1x),*(es25.16e3,1x))') row_index,ierr,successes,failures,loads, &
         res,density_partials,temperature_partials,abundance_partials
      flush(stream)
   end do
   call free_eos_handle(handle)
   call eos_shutdown()
end program native_driver
