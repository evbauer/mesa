program control_test
   use bridge_support
   use eos_lib, only: eos_get_control_namelist, eos_set_control_namelist, free_eos_handle, eos_shutdown
   implicit none
   character(len=4096) :: mesa_directory,massfile
   character(len=64) :: value
   integer :: handle,status
   logical :: enabled
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,massfile)
   call initialize_eos(trim(mesa_directory),trim(massfile),status)
   if(status /= 0) error stop 'initialize EOS'
   handle=alloc_eos_handle_using_inlist('',status)
   if(status /= 0) error stop 'default handle'
   if(eos_handles(handle)%FreeEOS_use_HELM_interpolation) error stop 'default must be false'
   call eos_get_control_namelist(handle,'FreeEOS_use_HELM_interpolation',value,status)
   if(status /= 0) error stop 'get default control'
   read(value,*) enabled
   if(enabled) error stop 'default namelist value'
   call eos_set_control_namelist(handle,'FreeEOS_use_HELM_interpolation','.true.',status)
   if(status /= 0 .or. .not.eos_handles(handle)%FreeEOS_use_HELM_interpolation) error stop 'enable'
   call eos_get_control_namelist(handle,'FreeEOS_use_HELM_interpolation',value,status)
   if(status /= 0) error stop 'get enabled control'
   read(value,*) enabled
   if(.not.enabled) error stop 'enabled namelist value'
   call eos_set_control_namelist(handle,'FreeEOS_use_HELM_interpolation','.false.',status)
   if(status /= 0 .or. eos_handles(handle)%FreeEOS_use_HELM_interpolation) error stop 'disable'
   call eos_get_control_namelist(handle,'FreeEOS_use_biquintic_interpolation',value,status)
   if(status == 0) error stop 'old control remains registered'
   call free_eos_handle(handle)
   call eos_shutdown()
   print *, 'PASS renamed control: false default, public getter/setter, and old name unregistered'
end program control_test
