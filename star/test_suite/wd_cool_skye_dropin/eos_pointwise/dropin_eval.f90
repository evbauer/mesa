! WP8 (scratch branch WDM_EOS_dropin only): point-wise eosDT evaluation along
! stellar profiles, to compare the Skye drop-in (WDM_AV_replaces_Skye) with the
! released cascade at IDENTICAL (rho, T, X).  Not part of the test case proper.
!
! usage: dropin_eval <eos inlist> <input file> <output file>
!   input : first line n; then n rows "logRho logT x(1..11)" with the isotopes
!           h1 he3 he4 c12 n14 o16 o18 ne20 ne22 mg24 si28 (the case's net)
!   output: header with result names; one row per point: logRho logT ierr res(1:nv)
! Prints the drop-in call/decline counters and the mean wall time per eosDT_get.
program dropin_eval

   use eos_def
   use eos_lib, only: eos_init, eos_shutdown, eos_ptr, free_eos_handle, &
      alloc_eos_handle_using_inlist, eosDT_get
   use eos_wdm_av_def, only: wdm_av_dropin_calls, wdm_av_dropin_declines, &
      WDM_AV_replaces_Skye_flag
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use const_lib, only: const_init
   use const_def, only: dp, i8
   use math_lib
   use utils_lib, only: set_nan, mesa_error

   implicit none

   character(len=256), parameter :: mesa_root = '/Users/bauer39/work/mesas/evbauer/mesa'
   integer, parameter :: nsp = 11
   character(len=8), parameter :: iso(nsp) = [character(len=8) :: &
      'h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'o18', 'ne20', 'ne22', 'mg24', 'si28']

   character(len=512) :: inlist, infile, outfile
   integer :: handle, ierr, n, k, j, iin, iout
   integer(i8) :: t0, t1, rate
   type (EoS_General_Info), pointer :: rq
   integer, pointer :: net_iso(:), chem_id(:)
   real(dp) :: xa(nsp), logRho, logT
   real(dp), allocatable :: d_dxa(:,:)
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
   real(dp) :: elapsed

   if (command_argument_count() /= 3) call mesa_error(__FILE__, __LINE__, 'usage: dropin_eval inlist in out')
   call get_command_argument(1, inlist)
   call get_command_argument(2, infile)
   call get_command_argument(3, outfile)

   call const_init(trim(mesa_root), ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'const_init failed')
   call math_init()
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'chem_init failed')
   call eos_init(' ', .true., ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_init failed')
   handle = alloc_eos_handle_using_inlist(trim(inlist), ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'alloc_eos_handle_using_inlist failed')
   call eos_ptr(handle, rq, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_ptr failed')

   write(*,'(a,a)') ' inlist                 ', trim(inlist)
   write(*,'(a,l2)') ' WDM_AV_replaces_Skye   ', WDM_AV_replaces_Skye_flag
   write(*,'(a,l2)') ' use_WDM_AV             ', rq% use_WDM_AV
   write(*,'(a,a)') ' WDM_AV_table_set       ', trim(rq% WDM_AV_table_set)
   write(*,'(a,a)') ' WDM_AV_mixing_rule     ', trim(rq% WDM_AV_mixing_rule)
   write(*,'(a,2f8.2)') ' Gamma_max all/any      ', rq% Gamma_max_for_all_WDM_AV, rq% Gamma_max_for_any_WDM_AV
   write(*,'(a,a)') ' Skye_solid_mixing_rule ', trim(rq% Skye_solid_mixing_rule)

   allocate(net_iso(num_chem_isos), chem_id(nsp), d_dxa(num_eos_d_dxa_results, nsp))
   net_iso(:) = 0
   do j = 1, nsp
      chem_id(j) = chem_get_iso_id(trim(iso(j)))
      if (chem_id(j) <= 0) call mesa_error(__FILE__, __LINE__, 'bad iso')
      net_iso(chem_id(j)) = j
   end do

   open(newunit=iin, file=trim(infile), action='read', status='old')
   read(iin,*) n
   open(newunit=iout, file=trim(outfile), action='write', status='replace')
   write(iout,'(a)', advance='no') 'logRho logT ierr'
   do j = 1, num_eos_basic_results
      write(iout,'(1x,a)', advance='no') trim(eosDT_result_names(j))
   end do
   write(iout,'(a)') ''

   elapsed = 0d0
   do k = 1, n
      read(iin,*) logRho, logT, xa
      call set_nan(res); call set_nan(d_dlnd); call set_nan(d_dlnT)
      ierr = 0
      call system_clock(t0, rate)
      call eosDT_get(handle, nsp, chem_id, net_iso, xa, &
                     exp10(logRho), logRho, exp10(logT), logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
      call system_clock(t1)
      elapsed = elapsed + real(t1 - t0, dp)/real(rate, dp)
      write(iout,'(2(es24.16e3,1x),i4,999(1x,es24.16e3))') logRho, logT, ierr, res
   end do
   close(iin)
   close(iout)

   write(*,'(a,i0)') ' points                 ', n
   write(*,'(a,es12.4)') ' mean seconds per call  ', elapsed/max(1, n)
   write(*,'(a,i0)') ' dropin calls           ', wdm_av_dropin_calls
   write(*,'(a,i0)') ' dropin declines        ', wdm_av_dropin_declines

   call free_eos_handle(handle)
   call eos_shutdown()

end program dropin_eval
