program cell_driver
   use const_def, only: dp
   use auto_diff_real_2var_order3_module
   use freeeos_biquintic
   use, intrinsic :: iso_fortran_env, only: error_unit, input_unit, output_unit, iostat_end, iostat_eor
   implicit none
   real(dp) :: nodal(0:2,0:2,2,2), bounds(4), query(2), corner_values(9)
   type(freeeos_partials) :: partials
   type(auto_diff_real_2var_order3) :: u_ad, tau_ad, f_ad
   integer :: status, corner_u, corner_tau, cell_number, record_number
   logical :: at_eof

   cell_number = 0
   do
      cell_number = cell_number + 1
      record_number = 1
      call read_record(bounds, .true., at_eof)
      if (at_eof) exit
      record_number = 2
      call read_record(query, .false., at_eof)
      do corner_tau = 1, 2
         do corner_u = 1, 2
            record_number = record_number + 1
            call read_record(corner_values, .false., at_eof)
            nodal(:,:,corner_u,corner_tau) = reshape(corner_values, [3, 3])
         end do
      end do
      call evaluate_cell_real(nodal, bounds(1), bounds(2), bounds(3), bounds(4), &
         query(1), query(2), partials, status)
      call check_status(status)
      u_ad = query(1)
      u_ad%d1val2 = 1.0_dp
      tau_ad = query(2)
      tau_ad%d1val1 = 1.0_dp
      call evaluate_cell_ad(nodal, bounds(1), bounds(2), bounds(3), bounds(4), u_ad, tau_ad, f_ad, status)
      call check_status(status)
      write(output_unit, '(a)') 'f f_u f_tau f_uu f_u_tau f_tau_tau f_uuu f_uu_tau f_u_tau_tau f_tau_tau_tau'
      write(output_unit, '(10(es25.16e3,1x))') partials%f, partials%f_u, partials%f_tau, partials%f_uu, &
         partials%f_u_tau, partials%f_tau_tau, partials%f_uuu, partials%f_uu_tau, &
         partials%f_u_tau_tau, partials%f_tau_tau_tau
      write(output_unit, '(10(es25.16e3,1x))') f_ad%val, f_ad%d1val2, f_ad%d1val1, f_ad%d2val2, &
         f_ad%d1val1_d1val2, f_ad%d2val1, f_ad%d3val2, f_ad%d1val1_d2val2, &
         f_ad%d2val1_d1val2, f_ad%d3val1
      flush(output_unit)
   end do

contains

   subroutine read_record(values, allow_eof, reached_eof)
      real(dp), intent(out) :: values(:)
      logical, intent(in) :: allow_eof
      logical, intent(out) :: reached_eof
      character(len=4096) :: record
      integer :: read_status, record_length, position, token_start, token_count

      reached_eof = .false.
      read(input_unit, '(a)', advance='no', size=record_length, iostat=read_status) record
      if (read_status == iostat_end .and. record_length == 0 .and. allow_eof) then
         reached_eof = .true.
         return
      end if
      if (read_status /= iostat_eor) call fail_input('missing or oversized record')
      token_count = 0
      position = 1
      do while (position <= record_length)
         if (is_separator(record(position:position))) then
            position = position + 1
            cycle
         end if
         token_start = position
         do while (position <= record_length)
            if (is_separator(record(position:position))) exit
            position = position + 1
         end do
         token_count = token_count + 1
         if (token_count > size(values)) call fail_input('too many values')
         if (scan(record(token_start:position-1), ',/*') /= 0) call fail_input('invalid numeric token')
         read(record(token_start:position-1), *, iostat=read_status) values(token_count)
         if (read_status /= 0) call fail_input('invalid numeric token')
      end do
      if (token_count /= size(values)) call fail_input('incomplete record')
   end subroutine read_record

   logical function is_separator(character_value)
      character, intent(in) :: character_value

      is_separator = index(' '//achar(9)//achar(13), character_value) /= 0
   end function is_separator

   subroutine fail_input(reason)
      character(len=*), intent(in) :: reason

      write(error_unit, '(a,i0,a,i0,a,a)') 'Invalid or incomplete cell input; cell=', cell_number, &
         '; record=', record_number, '; ', reason
      stop 2
   end subroutine fail_input

   subroutine check_status(cell_status)
      integer, intent(in) :: cell_status

      if (cell_status == cell_ok) return
      write(error_unit, '(a,i0,a,i0)') 'Cell rejected; cell=', cell_number, '; status=', cell_status
      stop 1
   end subroutine check_status

end program cell_driver
