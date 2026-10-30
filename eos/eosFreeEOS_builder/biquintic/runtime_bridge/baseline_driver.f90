program baseline_driver
   use iso_fortran_env, only: input_unit,output_unit,iostat_end
   use bridge_support
   use eos_lib, only: free_eos_handle,eos_shutdown
   implicit none
   character(len=4096) :: mesa_directory,inlist,massfile,record
   character(len=32) :: route
   real(dp) :: log_rho,log_temp,res(nv),density_partials(nv),temperature_partials(nv)
   real(dp) :: abundance_partials(num_eos_d_dxa_results,pinned_species)
   integer :: handle,ierr,io_status,stream,row_index

   open(newunit=stream,file='/dev/stdout',action='write')
   open(unit=output_unit,file='/dev/stderr',action='write')
   if (command_argument_count() /= 3) error stop 'usage: baseline_driver MESA INLIST MASSFILE'
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,inlist)
   call get_command_argument(3,massfile)
   call initialize_eos(trim(mesa_directory),trim(massfile),ierr)
   if (ierr /= 0) error stop 'EOS initialization failed'
   handle=alloc_eos_handle_using_inlist(trim(inlist),ierr)
   if (ierr /= 0) error stop 'handle allocation failed'
   row_index=0
   do
      read(input_unit,'(a)',iostat=io_status) record
      if (io_status == iostat_end) exit
      if (io_status /= 0) error stop 'stdin failed'
      if (len_trim(record) == 0) cycle
      read(record,*,iostat=io_status) route,log_rho,log_temp
      if (io_status /= 0) error stop 'bad query'
      call query_eos(trim(route),handle,log_rho,log_temp,res,density_partials,temperature_partials,abundance_partials,ierr)
      row_index=row_index+1
      write(stream,'(4(i0,1x),*(es25.16e3,1x))') row_index,ierr,0,1,res,density_partials,temperature_partials,abundance_partials
      flush(stream)
   end do
   call free_eos_handle(handle)
   call eos_shutdown()
end program baseline_driver
