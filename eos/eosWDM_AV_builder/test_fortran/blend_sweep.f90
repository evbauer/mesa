! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

program blend_sweep

   ! Blend-edge sweeps of the integrated WDM_AV eosDT component, evaluated
   ! through the FULL cascade (eosDT_get, not eosDT_get_component).  As the
   ! sweep coordinate crosses the WDM_AV blend box, Do_Blend ramps the result
   ! between WDM_AV and whichever lower cascade level covers the point (Skye,
   ! FreeEOS, HELM, ideal, ...).  The i_frac_* columns record the mix.
   !
   ! For two compositions (C/O 50/50 and pure c12) and two table sets ('skye'
   ! or 'helm', chosen by the first command-line argument), five 1-D cuts at
   ! step 0.005 dex are written:
   !
   !   rho_right_T6.0 : logT=6.0 fixed, logRho 3.0 -> 6.0  (right edge, into Skye)
   !   rho_right_T6.8 : logT=6.8 fixed, logRho 3.0 -> 6.0  (right edge, into Skye)
   !   T_top_R2.0     : logRho=2.0 fixed, logT 6.0 -> 8.0  (top edge, into FreeEOS)
   !   T_top_R4.0     : logRho=4.0 fixed, logT 6.0 -> 8.0  (top edge, into Skye)
   !   rho_left_T6.8  : logT=6.8 fixed, logRho -4.0 -> -1.5 (left edge)
   !
   ! Each row carries the sweep coordinate, logRho, logT, a flag, the res
   ! entries lnPgas,lnE,lnS,chiRho,chiT,Cv,Cp,grad_ad,gamma1, the d_dlnd and
   ! d_dlnT of lnPgas,lnE,lnS,Cv,grad_ad, and all eight i_frac_* values.
   !
   ! A 2-D frac map on the 10-points-per-decade reference grid (logRho -8..10,
   ! logT 3..9.5) is also written for C/O with each set.
   !
   ! Output: data/eosWDM_AV_data/blend/<set>/<mixture>_<cut>.txt and
   !         data/eosWDM_AV_data/blend/<set>/fracmap_c12_o16.txt
   !
   ! Radiation: WDM_AV and Skye add radiation internally (pack_for_export), so
   ! lnE,lnS,Cv,chiRho,chiT INCLUDE radiation; lnPgas is gas-only.  The Python
   ! analysis (validation/blend.py) is internal-consistency only (returned
   ! partials vs finite differences, value steps across the ramp), so no
   ! radiation handling is needed here.

   use eos_def
   use eos_lib, only: eos_init, eos_shutdown, eos_ptr, free_eos_handle, &
      alloc_eos_handle_using_inlist, eosDT_get
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use const_lib, only: const_init
   use const_def, only: dp, mesa_dir
   use math_lib
   use utils_lib, only: set_nan, is_bad, mesa_error

   implicit none

   ! locations under MESA_DIR, set after const_init
   character(len=512) :: inlist_path, data_root

   integer, parameter :: max_species = 2

   ! --- the five 1-D cuts -------------------------------------------------
   integer, parameter :: n_cut = 5
   character(len=16) :: cut_name(n_cut) = [character(len=16) :: &
      'rho_right_T6.0', 'rho_right_T6.8', 'T_top_R2.0', 'T_top_R4.0', 'rho_left_T6.8']
   ! axis: 1 => vary logRho at fixed logT ; 2 => vary logT at fixed logRho
   integer :: cut_axis(n_cut) = [1, 1, 2, 2, 1]
   real(dp) :: cut_fixed(n_cut) = [6.0d0, 6.8d0, 2.0d0, 4.0d0, 6.8d0]
   real(dp) :: cut_lo(n_cut) = [3.0d0, 3.0d0, 6.0d0, 6.0d0, -4.0d0]
   real(dp) :: cut_hi(n_cut) = [6.0d0, 6.0d0, 8.0d0, 8.0d0, -1.5d0]
   real(dp), parameter :: cut_step = 0.005d0

   integer :: handle, ierr, imix, icut
   type (EoS_General_Info), pointer :: rq
   character(len=32) :: table_set, rule_arg, set_label
   character(len=512) :: out_dir
   integer :: nargs

   ! composition + eos work arrays
   integer, pointer :: net_iso(:), chem_id(:)
   real(dp) :: xa(max_species)
   real(dp), allocatable :: d_dxa(:,:)
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT

   ! ------------------------------------------------------------------
   ! table set from the command line (default: whatever the inlist sets)
   nargs = command_argument_count()
   if (nargs >= 1) then
      call get_command_argument(1, table_set)
   else
      table_set = ''
   end if

   ! optional second argument: the mixing rule ('pressure' | 'electron_density'),
   ! overriding WDM_AV_mixing_rule from the inlist; '<set>_ne' output dir for the
   ! electron-density rule so it does not clobber the additive-volume sweeps.
   if (nargs >= 2) then
      call get_command_argument(2, rule_arg)
   else
      rule_arg = ''
   end if

   ! ------------------------------------------------------------------
   call const_init('', ierr)  ! MESA_DIR from the environment
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'const_init failed')
   inlist_path = trim(mesa_dir)//'/eos/eosWDM_AV_builder/test_fortran/inlist_blend_sweep'
   data_root = trim(mesa_dir)//'/data/eosWDM_AV_data'
   call math_init()
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'chem_init failed')
   call eos_init(' ', .true., ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_init failed')

   handle = alloc_eos_handle_using_inlist(trim(inlist_path), ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'alloc_eos_handle_using_inlist failed')

   call eos_ptr(handle, rq, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_ptr failed')

   ! force the WDM_AV path on and select the requested table set
   rq% use_WDM_AV = .true.
   if (len_trim(table_set) > 0) rq% WDM_AV_table_set = trim(table_set)
   if (len_trim(rule_arg) > 0) rq% WDM_AV_mixing_rule = trim(rule_arg)
   table_set = trim(rq% WDM_AV_table_set)

   set_label = trim(table_set)
   if (trim(rq% WDM_AV_mixing_rule) == 'electron_density') set_label = trim(table_set)//'_ne'
   out_dir = trim(data_root)//'/blend/'//trim(set_label)
   call execute_command_line('mkdir -p '//trim(out_dir))

   write(*,'(a)') '================ blend_sweep ================'
   write(*,'(a,a)') ' table set            ', trim(table_set)
   write(*,'(a,a)') ' WDM_AV_mixing_rule    ', trim(rq% WDM_AV_mixing_rule)
   write(*,'(a,a)') ' WDM_AV_mixing_entropy ', trim(rq% WDM_AV_mixing_entropy)
   write(*,'(a,a)') ' output dir           ', trim(out_dir)
   write(*,'(a)') ' blend box (logRho any/all/all/any, logT any/all/all/any):'
   write(*,'(4(a,f7.3))') '   logRho ', rq% logRho_min_for_any_WDM_AV, ' /', &
      rq% logRho_min_for_all_WDM_AV, ' /', rq% logRho_max_for_all_WDM_AV, ' /', &
      rq% logRho_max_for_any_WDM_AV
   write(*,'(4(a,f7.3))') '   logT   ', rq% logT_min_for_any_WDM_AV, ' /', &
      rq% logT_min_for_all_WDM_AV, ' /', rq% logT_max_for_all_WDM_AV, ' /', &
      rq% logT_max_for_any_WDM_AV
   write(*,'(a)') '============================================='

   allocate(net_iso(num_chem_isos), chem_id(max_species), &
            d_dxa(num_eos_d_dxa_results, max_species), stat=ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'allocate failed')

   do imix = 1, 2
      do icut = 1, n_cut
         call do_one_cut(imix, icut)
      end do
      if (imix == 1) call do_fracmap('c12_o16', 1)
   end do

   deallocate(net_iso, chem_id, d_dxa)
   call free_eos_handle(handle)
   call eos_shutdown()
   write(*,'(a)') 'blend_sweep: done.'

contains

   ! set composition for mixture index (1 = C/O 50/50, 2 = pure c12).
   ! returns the number of species and the mixture stem.
   subroutine set_mixture(imix, nsp, stem)
      integer, intent(in) :: imix
      integer, intent(out) :: nsp
      character(len=*), intent(out) :: stem
      integer :: cid_c12, cid_o16

      cid_c12 = chem_get_iso_id('c12')
      if (cid_c12 <= 0) call mesa_error(__FILE__, __LINE__, 'no c12')

      net_iso(:) = 0
      xa = 0d0
      if (imix == 1) then
         cid_o16 = chem_get_iso_id('o16')
         if (cid_o16 <= 0) call mesa_error(__FILE__, __LINE__, 'no o16')
         nsp = 2
         chem_id(1) = cid_c12; chem_id(2) = cid_o16
         net_iso(cid_c12) = 1; net_iso(cid_o16) = 2
         xa(1) = 0.5d0; xa(2) = 0.5d0
         stem = 'c12_o16'
      else
         nsp = 1
         chem_id(1) = cid_c12
         net_iso(cid_c12) = 1
         xa(1) = 1d0
         stem = 'c12'
      end if
   end subroutine set_mixture


   subroutine do_one_cut(imix, icut)
      integer, intent(in) :: imix, icut
      integer :: nsp, npts, k, iou, n_ok, n_bad
      character(len=32) :: stem
      character(len=512) :: out_file
      real(dp) :: coord, logRho, logT, Rho, T

      call set_mixture(imix, nsp, stem)

      npts = nint((cut_hi(icut) - cut_lo(icut)) / cut_step) + 1

      out_file = trim(out_dir)//'/'//trim(stem)//'_'//trim(cut_name(icut))//'.txt'
      open(newunit=iou, file=trim(out_file), action='write', status='replace')
      call write_cut_header(iou, stem, icut, nsp)

      n_ok = 0; n_bad = 0
      do k = 1, npts
         coord = cut_lo(icut) + cut_step * real(k-1, dp)
         if (cut_axis(icut) == 1) then
            logRho = coord
            logT = cut_fixed(icut)
         else
            logRho = cut_fixed(icut)
            logT = coord
         end if
         Rho = exp10(logRho)
         T = exp10(logT)

         call set_nan(res); call set_nan(d_dlnd); call set_nan(d_dlnT)
         ierr = 0
         call eosDT_get(handle, nsp, chem_id, net_iso, xa, &
                        Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)

         if (ierr /= 0 .or. is_bad(res(i_lnPgas))) then
            n_bad = n_bad + 1
            call write_cut_row(iou, coord, logRho, logT, 1, .true.)
            ierr = 0
         else
            n_ok = n_ok + 1
            call write_cut_row(iou, coord, logRho, logT, 0, .false.)
         end if
      end do

      close(iou)
      write(*,'(a,2x,i0,a,2x,i0,a)') &
         '  '//trim(stem)//'_'//trim(cut_name(icut)), n_ok, ' ok', n_bad, ' bad'
   end subroutine do_one_cut


   ! write one data row of a cut.  res/d_dlnd/d_dlnT are module-level and hold
   ! the current evaluation.  bad => write the coordinate + flag and NaNs.
   subroutine write_cut_row(iou, coord, logRho, logT, flag, bad)
      integer, intent(in) :: iou, flag
      real(dp), intent(in) :: coord, logRho, logT
      logical, intent(in) :: bad
      real(dp) :: nanv
      integer :: j

      if (bad) then
         call set_nan(nanv)
         write(iou,'(3(es24.16e3,1x),i4,1x,99(es24.16e3,1x))') &
            coord, logRho, logT, flag, (nanv, j=1,27)
      else
         write(iou,'(3(es24.16e3,1x),i4,1x,99(es24.16e3,1x))') &
            coord, logRho, logT, flag, &
            res(i_lnPgas), res(i_lnE), res(i_lnS), res(i_chiRho), res(i_chiT), &
            res(i_Cv), res(i_Cp), res(i_grad_ad), res(i_gamma1), &
            d_dlnd(i_lnPgas), d_dlnT(i_lnPgas), &
            d_dlnd(i_lnE), d_dlnT(i_lnE), &
            d_dlnd(i_lnS), d_dlnT(i_lnS), &
            d_dlnd(i_Cv), d_dlnT(i_Cv), &
            d_dlnd(i_grad_ad), d_dlnT(i_grad_ad), &
            res(i_frac_HELM), res(i_frac_OPAL_SCVH), res(i_frac_FreeEOS), &
            res(i_frac_PC), res(i_frac_Skye), res(i_frac_CMS), &
            res(i_frac_ideal), res(i_frac_WDM_AV)
      end if
   end subroutine write_cut_row


   subroutine write_cut_header(iou, stem, icut, nsp)
      integer, intent(in) :: iou, icut, nsp
      character(len=*), intent(in) :: stem
      character(len=16) :: axis_str
      if (cut_axis(icut) == 1) then
         axis_str = 'logRho'
      else
         axis_str = 'logT'
      end if
      write(iou,'(a)') '# WDM_AV blend-edge sweep (full cascade, eosDT_get), version 1'
      write(iou,'(a,a,a)') 'table_set            "', trim(table_set), '"'
      write(iou,'(a,a,a)') 'mixing_entropy       "', trim(rq% WDM_AV_mixing_entropy), '"'
      write(iou,'(a,a,a)') 'mixture              "', trim(stem), '"'
      write(iou,'(a,i0)')  'nspecies             ', nsp
      write(iou,'(a,a,a)') 'cut                  "', trim(cut_name(icut)), '"'
      write(iou,'(a,a)')   'sweep_axis           ', trim(axis_str)
      write(iou,'(a,f10.4)') 'fixed_value          ', cut_fixed(icut)
      write(iou,'(a,f10.4,a,f10.4,a,f10.4)') &
         'coord_min            ', cut_lo(icut), '  coord_max  ', cut_hi(icut), &
         '  step  ', cut_step
      write(iou,'(a,4(1x,f7.3))') 'box_logRho_any_all_all_any', &
         rq% logRho_min_for_any_WDM_AV, rq% logRho_min_for_all_WDM_AV, &
         rq% logRho_max_for_all_WDM_AV, rq% logRho_max_for_any_WDM_AV
      write(iou,'(a,4(1x,f7.3))') 'box_logT_any_all_all_any  ', &
         rq% logT_min_for_any_WDM_AV, rq% logT_min_for_all_WDM_AV, &
         rq% logT_max_for_all_WDM_AV, rq% logT_max_for_any_WDM_AV
      write(iou,'(a)') 'radiation_note       "lnE,lnS,Cv,chiRho,chiT include radiation '// &
         '(pack_for_export); lnPgas is gas-only"'
      write(iou,'(a)') 'columns              coord logRho logT flag lnPgas lnE lnS chiRho chiT '// &
         'Cv Cp grad_ad gamma1 dlnPgas_dlnRho dlnPgas_dlnT dlnE_dlnRho dlnE_dlnT '// &
         'dlnS_dlnRho dlnS_dlnT dCv_dlnRho dCv_dlnT dgrad_ad_dlnRho dgrad_ad_dlnT '// &
         'frac_HELM frac_OPAL_SCVH frac_FreeEOS frac_PC frac_Skye frac_CMS frac_ideal frac_WDM_AV'
      write(iou,'(a)') '# flag: 0 valid, 1 bad (ierr or NaN result); NaN data columns when flag /= 0'
   end subroutine write_cut_header


   ! 2-D blend-fraction map on the 10-points-per-decade reference grid.
   subroutine do_fracmap(stem_in, imix)
      character(len=*), intent(in) :: stem_in
      integer, intent(in) :: imix
      real(dp), parameter :: lr_min = -8d0, lr_max = 10d0, lt_min = 3d0, lt_max = 9.5d0
      integer, parameter :: nR = 181, nT = 66  ! 10 pts/decade inclusive
      integer :: nsp, iou, jR, jT, n_ok, n_bad
      character(len=32) :: stem
      character(len=512) :: out_file
      real(dp) :: logRho, logT, Rho, T

      call set_mixture(imix, nsp, stem)

      out_file = trim(out_dir)//'/fracmap_'//trim(stem_in)//'.txt'
      open(newunit=iou, file=trim(out_file), action='write', status='replace')
      write(iou,'(a)') '# WDM_AV blend-fraction map (full cascade, eosDT_get), version 1'
      write(iou,'(a,a,a)') 'table_set            "', trim(table_set), '"'
      write(iou,'(a,a,a)') 'mixture              "', trim(stem), '"'
      write(iou,'(a,i0,a,f8.3,a,f8.3)') &
         'nlogRho              ', nR, '     logRho_min  ', lr_min, '    logRho_max  ', lr_max
      write(iou,'(a,i0,a,f8.3,a,f8.3)') &
         'nlogT                ', nT, '     logT_min    ', lt_min, '    logT_max    ', lt_max
      write(iou,'(a)') 'columns              logRho logT flag frac_HELM frac_OPAL_SCVH '// &
         'frac_FreeEOS frac_PC frac_Skye frac_CMS frac_ideal frac_WDM_AV'
      write(iou,'(a)') '# rows: logT outer loop, logRho inner loop'

      n_ok = 0; n_bad = 0
      do jT = 1, nT
         logT = lt_min + (lt_max - lt_min) * real(jT-1, dp) / real(nT-1, dp)
         T = exp10(logT)
         do jR = 1, nR
            logRho = lr_min + (lr_max - lr_min) * real(jR-1, dp) / real(nR-1, dp)
            Rho = exp10(logRho)
            call set_nan(res); call set_nan(d_dlnd); call set_nan(d_dlnT)
            ierr = 0
            call eosDT_get(handle, nsp, chem_id, net_iso, xa, &
                           Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
            if (ierr /= 0 .or. is_bad(res(i_lnPgas))) then
               n_bad = n_bad + 1
               write(iou,'(2(es24.16e3,1x),i4,1x,8(es24.16e3,1x))') &
                  logRho, logT, 1, &
                  res(i_frac_HELM), res(i_frac_OPAL_SCVH), res(i_frac_FreeEOS), &
                  res(i_frac_PC), res(i_frac_Skye), res(i_frac_CMS), &
                  res(i_frac_ideal), res(i_frac_WDM_AV)
               ierr = 0
            else
               n_ok = n_ok + 1
               write(iou,'(2(es24.16e3,1x),i4,1x,8(es24.16e3,1x))') &
                  logRho, logT, 0, &
                  res(i_frac_HELM), res(i_frac_OPAL_SCVH), res(i_frac_FreeEOS), &
                  res(i_frac_PC), res(i_frac_Skye), res(i_frac_CMS), &
                  res(i_frac_ideal), res(i_frac_WDM_AV)
            end if
         end do
      end do

      close(iou)
      write(*,'(a,2x,i0,a,2x,i0,a)') '  fracmap_'//trim(stem_in), n_ok, ' ok', n_bad, ' bad'
   end subroutine do_fracmap

end program blend_sweep
