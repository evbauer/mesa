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

program mix_vs_reference

   ! Mixture-level validation of the integrated WDM_AV eosDT component.
   !
   ! For each of the five reference mixtures (same compositions as
   ! make_mesa_tables), this program reads the composition and grid from the
   ! reference file header, then at every reference grid node calls the
   ! integrated path
   !
   !    eosDT_get_component(handle, i_eos_WDM_AV, ...)
   !
   ! and writes the full WDM_AV result vector (plus the first partials that the
   ! reference files carry) to
   !
   !    data/eosWDM_AV_data/validation_mixtures/<set>/<mixture>.txt
   !
   ! The table set ('skye' or 'helm') is chosen by the first command-line
   ! argument and overrides WDM_AV_table_set from inlist_mix_vs_reference.
   !
   ! Radiation: WDM_AV adds radiation internally (via pack_for_export, like
   ! Skye), so the written lnE, lnS, Cv, chiRho, chiT INCLUDE radiation for both
   ! table sets.  The Python comparison (validation/mixtures.py) strips radiation
   ! analytically for the HELM reference (which is gas-only) and leaves the Skye
   ! comparison alone.  The flag column records 0 for a valid result and a
   ! non-zero code for a declined / failed node (the result columns are NaN).

   use eos_def
   use eos_lib, only: eos_init, eos_shutdown, eos_ptr, free_eos_handle, &
      alloc_eos_handle_using_inlist, eosDT_get_component
   use chem_def, only: chem_isos, num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use const_lib, only: const_init
   use const_def, only: dp, mesa_dir
   use math_lib
   use utils_lib, only: set_nan, is_bad, mesa_error

   implicit none

   character(len=256), parameter :: mesa_root = &
      '/Users/bauer39/work/mesas/evbauer/mesa'
   character(len=512), parameter :: inlist_path = &
      '/Users/bauer39/work/mesas/evbauer/mesa/eos/eosWDM_AV_builder/test_fortran/inlist_mix_vs_reference'
   character(len=512), parameter :: data_root = &
      '/Users/bauer39/work/mesas/evbauer/mesa/data/eosWDM_AV_data'

   integer, parameter :: max_species = 8
   integer, parameter :: n_mix = 5

   ! the five reference mixtures (file stems, shared by both reference sets)
   character(len=32) :: mix_files(n_mix) = [character(len=32) :: &
      'mix_c12_o16', 'mix_o16_ne20_mg24', 'mix_he4_c12', 'mix_h1_he4', &
      'mix_c12_o16_ne22']

   integer :: handle, ierr, i, nargs
   type (EoS_General_Info), pointer :: rq
   character(len=32) :: table_set
   character(len=512) :: ref_dir, out_dir

   ! composition + eos work arrays
   integer, pointer :: net_iso(:), chem_id(:)
   real(dp) :: xa(max_species)
   real(dp), allocatable :: d_dxa(:,:)
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
   real(dp) :: nanv

   call set_nan(nanv)

   ! ------------------------------------------------------------------
   ! table set from the command line (default: whatever the inlist sets)
   nargs = command_argument_count()
   if (nargs >= 1) then
      call get_command_argument(1, table_set)
   else
      table_set = ''
   end if

   ! ------------------------------------------------------------------
   ! initialise MESA infrastructure
   call const_init(trim(mesa_root), ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'const_init failed')
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
   table_set = trim(rq% WDM_AV_table_set)

   if (trim(table_set) == 'skye') then
      ref_dir = trim(data_root)//'/reference_skye_mixtures'
   else if (trim(table_set) == 'helm') then
      ref_dir = trim(data_root)//'/reference_helm_mixtures'
   else
      write(*,'(a)') 'unknown table set "'//trim(table_set)//'" (expected skye or helm)'
      call mesa_error(__FILE__, __LINE__, 'bad table set')
   end if
   out_dir = trim(data_root)//'/validation_mixtures/'//trim(table_set)
   call execute_command_line('mkdir -p '//trim(out_dir))

   write(*,'(a)') '================ mix_vs_reference ================'
   write(*,'(a,a)') ' table set            ', trim(table_set)
   write(*,'(a,a)') ' WDM_AV_mixing_entropy ', trim(rq% WDM_AV_mixing_entropy)
   write(*,'(a,a)') ' reference dir        ', trim(ref_dir)
   write(*,'(a,a)') ' output dir           ', trim(out_dir)
   write(*,'(a)') '=================================================='

   allocate(net_iso(num_chem_isos), chem_id(max_species), &
            d_dxa(num_eos_d_dxa_results, max_species), stat=ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'allocate failed')

   do i = 1, n_mix
      call do_one_mixture(trim(mix_files(i)))
   end do

   deallocate(net_iso, chem_id, d_dxa)
   call free_eos_handle(handle)
   call eos_shutdown()

   write(*,'(a)') 'mix_vs_reference: done.'

contains

   ! --------------------------------------------------------------------------
   ! one mixture: read header (composition + grid) from the reference file,
   ! evaluate WDM_AV on that grid, write the output table.
   ! --------------------------------------------------------------------------
   subroutine do_one_mixture(stem)
      character(len=*), intent(in) :: stem
      character(len=512) :: ref_file, out_file
      integer :: nsp, cids(max_species)
      real(dp) :: fracs(max_species)
      character(len=32) :: names(max_species)
      integer :: nlogRho, nlogT
      real(dp) :: logRho_min, logRho_max, logT_min, logT_max
      character(len=256) :: comp_str
      integer :: iou, jT, jRho, k, flag
      integer :: n_ok, n_skip, n_err
      real(dp) :: logRho, logT, Rho, T

      ref_file = trim(ref_dir)//'/'//trim(stem)//'.txt'
      out_file = trim(out_dir)//'/'//trim(stem)//'.txt'

      call read_reference_header(trim(ref_file), comp_str, &
         nlogRho, logRho_min, logRho_max, nlogT, logT_min, logT_max, ierr)
      if (ierr /= 0) then
         write(*,'(a)') '  SKIP mixture (cannot read header): '//trim(ref_file)
         return
      end if

      call parse_composition(comp_str, nsp, names, fracs)
      if (nsp < 1) then
         write(*,'(a)') '  SKIP mixture (empty composition): '//trim(stem)
         return
      end if

      ! build chem_id / net_iso / xa
      net_iso(:) = 0
      xa = 0d0
      do k = 1, nsp
         cids(k) = chem_get_iso_id(trim(names(k)))
         if (cids(k) <= 0) then
            write(*,'(a)') '  SKIP mixture, unknown isotope: '//trim(names(k))
            return
         end if
         chem_id(k) = cids(k)
         net_iso(cids(k)) = k
         xa(k) = fracs(k)
      end do

      write(*,'(a)') '--- '//trim(stem)//'  ['//trim(comp_str)//'] ---'

      open(newunit=iou, file=trim(out_file), action='write', status='replace')
      call write_output_header(iou, comp_str, nlogRho, logRho_min, logRho_max, &
                               nlogT, logT_min, logT_max)

      n_ok = 0; n_skip = 0; n_err = 0
      do jT = 1, nlogT
         logT = logT_min + (logT_max - logT_min) * real(jT-1, dp) / real(max(nlogT-1, 1), dp)
         T = exp10(logT)
         do jRho = 1, nlogRho
            logRho = logRho_min + (logRho_max - logRho_min) * real(jRho-1, dp) / real(max(nlogRho-1, 1), dp)
            Rho = exp10(logRho)

            ! pre-fill NaN so a silent skip (ierr=0, res untouched) is detected
            call set_nan(res)
            call set_nan(d_dlnd)
            call set_nan(d_dlnT)

            call eosDT_get_component(handle, i_eos_WDM_AV, nsp, chem_id, net_iso, xa, &
                                     Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)

            flag = 0
            if (ierr /= 0) then
               flag = ierr
               n_err = n_err + 1
               ierr = 0
            else if (is_bad(res(i_lnPgas)) .or. &
                     (res(i_lnPgas) == 0d0 .and. res(i_lnE) == 0d0 .and. &
                      res(i_lnS) == 0d0)) then
               ! declined: either res was left at the pre-filled NaN (component
               ! untouched), or it is the all-zero sentinel that eos_wdm_av_get
               ! leaves (it does "res = 0d0" up front) when its own solver
               ! declines (no-converge / bad-domain / out-of-table) and the
               ! integrated path swallows that ierr into "skip".
               flag = -1
               n_skip = n_skip + 1
            else
               n_ok = n_ok + 1
            end if

            if (flag /= 0) then
               write(iou,'(2(es24.16e3,1x),i4,1x,99(es24.16e3,1x))') &
                  logRho, logT, flag, (nan_val(), k=1,23)
            else
               write(iou,'(2(es24.16e3,1x),i4,1x,99(es24.16e3,1x))') &
                  logRho, logT, flag, &
                  res(i_lnPgas), res(i_lnE), res(i_lnS), res(i_chiRho), res(i_chiT), &
                  res(i_Cv), res(i_Cp), res(i_grad_ad), res(i_gamma1), res(i_gamma3), &
                  res(i_eta), res(i_mu), res(i_lnfree_e), &
                  d_dlnd(i_lnPgas), d_dlnT(i_lnPgas), &
                  d_dlnd(i_lnE), d_dlnT(i_lnE), &
                  d_dlnd(i_lnS), d_dlnT(i_lnS), &
                  d_dlnd(i_Cv), d_dlnT(i_Cv), &
                  d_dlnd(i_grad_ad), d_dlnT(i_grad_ad)
            end if
         end do
      end do

      close(iou)
      write(*,'(a,1x,i0)') '   valid nodes   ', n_ok
      write(*,'(a,1x,i0)') '   declined nodes', n_skip
      write(*,'(a,1x,i0)') '   error nodes   ', n_err
   end subroutine do_one_mixture


   real(dp) function nan_val()
      call set_nan(nan_val)
   end function nan_val


   subroutine write_output_header(iou, comp_str, nlogRho, logRho_min, logRho_max, &
                                  nlogT, logT_min, logT_max)
      integer, intent(in) :: iou, nlogRho, nlogT
      character(len=*), intent(in) :: comp_str
      real(dp), intent(in) :: logRho_min, logRho_max, logT_min, logT_max
      write(iou,'(a)') '# WDM_AV mixture validation output, version 1'
      write(iou,'(a,a)') 'component            i_eos_WDM_AV (eosDT_get_component, integrated path)'
      write(iou,'(a,a,a)') 'table_set            "', trim(table_set), '"'
      write(iou,'(a,a,a)') 'mixing_entropy       "', trim(rq% WDM_AV_mixing_entropy), '"'
      write(iou,'(a,a,a)') 'composition          "', trim(comp_str), '"'
      write(iou,'(a)') 'radiation_note       "raw WDM_AV result vector; radiation INCLUDED in '// &
         'lnE,lnS,Cv,chiRho,chiT (pack_for_export, like Skye); lnPgas is gas-only."'
      write(iou,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogRho              ', nlogRho, '     logRho_min  ', logRho_min, '    logRho_max  ', logRho_max
      write(iou,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogT                ', nlogT, '     logT_min    ', logT_min, '    logT_max    ', logT_max
      write(iou,'(a)') 'columns              logRho logT flag lnPgas lnE lnS chiRho chiT Cv Cp '// &
         'grad_ad gamma1 gamma3 eta mu lnfree_e dlnPgas_dlnRho dlnPgas_dlnT dlnE_dlnRho dlnE_dlnT '// &
         'dlnS_dlnRho dlnS_dlnT dCv_dlnRho dCv_dlnT dgrad_ad_dlnRho dgrad_ad_dlnT'
      write(iou,'(a)') '# flag: 0 valid, -1 declined (out of table / solver skip; '// &
         'res was NaN or the all-zero sentinel), else ierr; NaN result columns when flag /= 0'
      write(iou,'(a)') '# rows: logT outer loop, logRho inner loop'
   end subroutine write_output_header


   ! --------------------------------------------------------------------------
   ! read the reference-mixture header: composition string + grid definition.
   ! header lines have the form  'keyword   value...'  (no leading '#').
   ! --------------------------------------------------------------------------
   subroutine read_reference_header(fname, comp_str, nlogRho, logRho_min, logRho_max, &
                                    nlogT, logT_min, logT_max, ierr)
      character(len=*), intent(in) :: fname
      character(len=*), intent(out) :: comp_str
      integer, intent(out) :: nlogRho, nlogT
      real(dp), intent(out) :: logRho_min, logRho_max, logT_min, logT_max
      integer, intent(out) :: ierr
      integer :: u, ios, p, q
      character(len=1024) :: line
      character(len=64) :: key
      logical :: have_comp, have_rho, have_t

      ierr = 0
      comp_str = ''
      nlogRho = 0; nlogT = 0
      logRho_min = 0d0; logRho_max = 0d0; logT_min = 0d0; logT_max = 0d0
      have_comp = .false.; have_rho = .false.; have_t = .false.

      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) then; ierr = -1; return; end if

      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         ! stop at the first data row (begins with a number or sign)
         if (scan(line(1:1), '0123456789+-.') > 0) exit
         read(line,*,iostat=ios) key
         if (ios /= 0) cycle
         select case (trim(key))
         case ('composition')
            p = index(line, '"')
            q = index(line, '"', back=.true.)
            if (p > 0 .and. q > p) then
               comp_str = line(p+1:q-1)
               have_comp = .true.
            end if
         case ('nlogRho')
            call parse_grid_line(line, nlogRho, logRho_min, logRho_max)
            have_rho = .true.
         case ('nlogT')
            call parse_grid_line(line, nlogT, logT_min, logT_max)
            have_t = .true.
         end select
      end do
      close(u)

      if (.not. (have_comp .and. have_rho .and. have_t)) ierr = -1
   end subroutine read_reference_header


   ! parse e.g.  'nlogRho   181   logRho_min  -8.0000  logRho_max  10.0000'
   subroutine parse_grid_line(line, n, vmin, vmax)
      character(len=*), intent(in) :: line
      integer, intent(out) :: n
      real(dp), intent(out) :: vmin, vmax
      character(len=64) :: key, k2, k3
      read(line,*) key, n, k2, vmin, k3, vmax
   end subroutine parse_grid_line


   ! parse e.g.  'c12=0.50 o16=0.50'  ->  names, fracs (normalised)
   subroutine parse_composition(comp_str, nsp, names, fracs)
      character(len=*), intent(in) :: comp_str
      integer, intent(out) :: nsp
      character(len=*), intent(out) :: names(:)
      real(dp), intent(out) :: fracs(:)
      character(len=256) :: s, tok
      integer :: pos, sp, eq
      real(dp) :: fsum

      nsp = 0
      names = ''
      fracs = 0d0
      s = adjustl(comp_str)
      do
         s = adjustl(s)
         if (len_trim(s) == 0) exit
         sp = index(s, ' ')
         if (sp == 0) then
            tok = s
            s = ''
         else
            tok = s(1:sp-1)
            s = s(sp+1:)
         end if
         if (len_trim(tok) == 0) cycle
         eq = index(tok, '=')
         if (eq <= 0) cycle
         nsp = nsp + 1
         names(nsp) = tok(1:eq-1)
         read(tok(eq+1:), *) fracs(nsp)
      end do

      fsum = sum(fracs(1:nsp))
      if (fsum > 0d0) fracs(1:nsp) = fracs(1:nsp)/fsum
   end subroutine parse_composition

end program mix_vs_reference
