program test_loader
   use const_def, only: dp
   use freeeos_table, only: biquintic_table
   use freeeos_supplement, only: load_source_supplement
   implicit none
   type(biquintic_table) :: table
   character(len=2048) :: source_path,supplement_path,argument
   character(len=64) :: source_sha,supplement_sha
   real(dp) :: hydrogen,metals
   integer :: family_id,expected_status,status
   call get_command_argument(1,source_path)
   call get_command_argument(2,supplement_path)
   call get_command_argument(3,argument); read(argument,*) hydrogen
   call get_command_argument(4,argument); read(argument,*) metals
   call get_command_argument(5,argument); read(argument,*) family_id
   call get_command_argument(6,source_sha)
   call get_command_argument(7,supplement_sha)
   call get_command_argument(8,argument); read(argument,*) expected_status
   call load_source_supplement(trim(source_path),trim(supplement_path),hydrogen,metals,family_id, &
      source_sha,supplement_sha,table,status)
   if (status /= expected_status) then
      write(*,*) 'expected/status',expected_status,status
      error stop 'loader status'
   end if
   write(*,'(a,i0)') 'PASS source+supplement load status ',status
end program test_loader
