program table_query
   use const_def, only: dp
   use freeeos_table
   use auto_diff_real_2var_order3_module
   use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit, iostat_end
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   type(biquintic_table) :: table
   type(table_metadata) :: metadata
   type(potential_partials) :: partials
   type(thermodynamic_state) :: state
   type(auto_diff_real_2var_order3) :: x_ad, t_ad, potential_ad
   character(len=4096) :: path
   character(len=128) :: tokens(2), argument, expected_sha
   real(dp) :: expected_x, expected_z, x, t, nan
   integer :: status, count, io_status, real_status, ad_status, state_status, state_index

   if (command_argument_count() < 3 .or. command_argument_count() > 4) then
      write(error_unit,'(a)') 'Usage: table_query CACHE EXPECTED_X EXPECTED_Z [EXPECTED_SOURCE_SHA256]'
      stop 2
   end if
   call get_command_argument(1,path,status=io_status)
   if (io_status /= 0) call input_error('cache path too long')
   call get_command_argument(2,argument,status=io_status)
   if (io_status /= 0) call input_error('invalid expected X')
   if (verify(trim(argument),'0123456789+-.eEdD') /= 0) call input_error('invalid expected X')
   read(argument,*,iostat=io_status) expected_x
   if (io_status /= 0) call input_error('invalid expected X')
   call get_command_argument(3,argument,status=io_status)
   if (io_status /= 0) call input_error('invalid expected Z')
   if (verify(trim(argument),'0123456789+-.eEdD') /= 0) call input_error('invalid expected Z')
   read(argument,*,iostat=io_status) expected_z
   if (io_status /= 0) call input_error('invalid expected Z')
   if (command_argument_count() == 4) then
      call get_command_argument(4,expected_sha,status=io_status)
      if (io_status /= 0) call input_error('invalid expected SHA')
      call load_table(trim(path),expected_x,expected_z,table,status,trim(expected_sha))
   else
      call load_table(trim(path),expected_x,expected_z,table,status)
   end if
   if (status /= table_ok) then
      write(error_unit,'(a,i0)') 'Table load rejected; status=',status
      stop 1
   end if
   metadata=get_metadata(table)
   write(output_unit,'(a)',advance='no') 'x t status_real status_ad status_state'// &
      ' f f_x f_t f_xx f_xt f_tt f_xxx f_xxt f_xtt f_ttt'// &
      ' ad_f ad_f_x ad_f_t ad_f_xx ad_f_xt ad_f_tt ad_f_xxx ad_f_xxt ad_f_xtt ad_f_ttt'
   do state_index=1,state_count
      write(output_unit,'(a)',advance='no') ' '//trim(state_names(state_index))// &
         ' '//trim(state_names(state_index))//'_x '//trim(state_names(state_index))//'_t'
   end do
   write(output_unit,*)
   flush(output_unit)
   nan=ieee_value(0.0_dp,ieee_quiet_nan)
   do
      call read_tokens(input_unit,tokens,count,status)
      if (status == iostat_end) exit
      if (status /= 0 .or. count /= 2) call input_error('expected exactly two physical coordinates')
      if (verify(trim(tokens(1)),'0123456789+-.eEdD') /= 0 .or. &
          verify(trim(tokens(2)),'0123456789+-.eEdD') /= 0) call input_error('invalid numeric query')
      read(tokens(1),*,iostat=io_status) x
      if (io_status /= 0) call input_error('invalid x')
      read(tokens(2),*,iostat=io_status) t
      if (io_status /= 0) call input_error('invalid t')
      call evaluate_table_real(table,x,t,partials,real_status)
      x_ad=x
      x_ad%d1val2=1.0_dp
      t_ad=t
      t_ad%d1val1=1.0_dp
      call evaluate_table_ad(table,x_ad,t_ad,potential_ad,ad_status)
      state%value=nan
      state%dlnrho=nan
      state%dlnT=nan
      state_status=real_status
      if (real_status == table_ok) call derive_state(partials,x,t,metadata%crad,state,state_status)
      write(output_unit,'(2(es25.16e3,1x),3(i0,1x),20(es25.16e3,1x))',advance='no') &
         x,t,real_status,ad_status,state_status,partial_vector(partials),ad_vector(potential_ad)
      do state_index=1,state_count
         write(output_unit,'(3(es25.16e3,1x))',advance='no') &
            state%value(state_index),state%dlnrho(state_index),state%dlnT(state_index)
      end do
      write(output_unit,*)
      flush(output_unit)
   end do

contains

   subroutine input_error(reason)
      character(len=*), intent(in) :: reason
      write(error_unit,'(a)') 'Invalid driver input: '//reason
      stop 2
   end subroutine input_error

end program table_query
