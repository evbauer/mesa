program dual_handle_driver
   use iso_fortran_env, only: int64
   use bridge_support
   use eos_lib, only: free_eos_handle, eos_shutdown
   use freeeos_family, only: freeeos_family_get_counters
   implicit none
   character(len=4096) :: mesa_directory, configured_inlist, ordinary_inlist, massfile
   real(dp) :: res(nv), density_partials(nv), temperature_partials(nv)
   real(dp) :: abundance_partials(num_eos_d_dxa_results,pinned_species)
   real(dp) :: metal_pattern(pinned_species-2)
   integer :: configured_handle, ordinary_handle, handle, status, iteration
   integer(int64) :: successes, failures, loads

   if (command_argument_count() /= 4) error stop 'arguments'
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,configured_inlist)
   call get_command_argument(3,ordinary_inlist)
   call get_command_argument(4,massfile)
   call initialize_eos(trim(mesa_directory),trim(massfile),status)
   if (status /= 0) error stop 'EOS initialization'
   configured_handle=alloc_eos_handle_using_inlist(trim(configured_inlist),status)
   if (status /= 0) error stop 'configured handle'
   ordinary_handle=alloc_eos_handle_using_inlist(trim(ordinary_inlist),status)
   if (status /= 0) error stop 'ordinary handle'
   metal_pattern=pinned_xa(3:)/sum(pinned_xa(3:))
   pinned_xa(1)=0
   pinned_xa(2)=0.15_dp
   pinned_xa(3:)=0.85_dp*metal_pattern
   do iteration=1,6
      if (mod(iteration,2) == 1) then
         handle=configured_handle
      else
         handle=ordinary_handle
      end if
      call query_eos('component',handle,-7.0_dp,5.0_dp, &
         res,density_partials,temperature_partials,abundance_partials,status)
      write(*,'(i0,1x,i0,1x,es25.16e3)') iteration,status,res(i_lnPgas)
   end do
   call freeeos_family_get_counters(successes,failures,loads)
   write(*,'(a,3(1x,i0))') 'COUNTS',successes,failures,loads
   call free_eos_handle(configured_handle)
   call free_eos_handle(ordinary_handle)
   call eos_shutdown()
end program dual_handle_driver
