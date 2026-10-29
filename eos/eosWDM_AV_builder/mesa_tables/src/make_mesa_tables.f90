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

program make_mesa_tables

   ! Writes WDM_AV source-format ASCII tables (FORMAT.md Section 1) for a list of
   ! pure isotopes, evaluated with MESA's analytic EOS components (HELM and Skye),
   ! plus reference mixture data (full MESA result vector) for later validation of
   ! the additive-volume mixture against Skye/HELM run directly on the mixture.
   !
   ! All output is gas-only (radiation stripped). HELM honours include_radiation
   ! (= .false. here); Skye's component output ignores include_radiation, so for
   ! Skye the radiation contributions to E, S and Cv are subtracted analytically:
   !    e_rad = crad*T^4/rho,  s_rad = 4*crad*T^3/(3*rho),  cv_rad = 4*crad*T^3/rho.
   ! The gas pressure (res(i_lnPgas)) already excludes radiation in both components.

   use eos_def
   use eos_lib, only: eos_init, eos_shutdown, eos_ptr, free_eos_handle, &
      alloc_eos_handle_using_inlist, eosDT_get_component
   use chem_def
   use chem_lib, only: chem_init, chem_get_iso_id, basic_composition_info
   use const_lib, only: const_init
   use const_def, only: dp, mesa_dir, crad, kerg, amu, avo, arg_not_provided
   use math_lib
   use utils_lib, only: set_nan, is_bad

   implicit none

   integer, parameter :: max_species = 8
   integer, parameter :: max_iso = 64

   ! namelist controls
   real(dp) :: logRho_min, logRho_max, logT_min, logT_max
   integer :: pts_per_decade, pts_per_decade_mix
   integer :: n_iso
   character(len=32) :: iso_names(max_iso)
   character(len=256) :: output_root
   ! run-selection controls: skip the (unchanged) isotope source tables and/or
   ! restrict the mixture loop to a sub-range, so a run can regenerate only the
   ! new mixtures without overwriting existing isotope / mixture files.
   logical :: skip_isotope_tables
   integer :: mix_first, mix_last

   namelist /table_controls/ &
      logRho_min, logRho_max, logT_min, logT_max, &
      pts_per_decade, pts_per_decade_mix, n_iso, iso_names, output_root, &
      skip_isotope_tables, mix_first, mix_last

   ! eos handle and runtime options
   integer :: handle
   type (EoS_General_Info), pointer :: rq
   logical :: helm_has_rad  ! HELM include_radiation setting (should be .false.)
   integer :: ierr

   ! grids
   integer :: nlogRho, nlogT, nlogRho_mix, nlogT_mix

   ! work arrays for a single eos call
   integer, pointer :: net_iso(:), chem_id(:)
   real(dp) :: xa(max_species)
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
   real(dp), allocatable :: d_dxa(:,:)

   character(len=256) :: version_string, dir_helm, dir_skye, dir_mix_skye, dir_mix_helm
   real(dp) :: nan_val

   integer :: i

   ! mixture definitions (hardcoded per WP4; cases 1-5 WP4, cases 6-11 WP4c)
   integer, parameter :: n_mix = 11

   include 'formats'

   ! ------------------------------------------------------------------
   ! defaults
   logRho_min = -8d0
   logRho_max = 10d0
   logT_min = 3d0
   logT_max = 9.5d0
   pts_per_decade = 20
   pts_per_decade_mix = 10
   output_root = ''
   skip_isotope_tables = .false.
   mix_first = 1
   mix_last = n_mix
   iso_names = ''
   n_iso = 7
   iso_names(1) = 'h1'
   iso_names(2) = 'he4'
   iso_names(3) = 'c12'
   iso_names(4) = 'o16'
   iso_names(5) = 'ne20'
   iso_names(6) = 'ne22'
   iso_names(7) = 'mg24'

   ierr = 0
   call set_nan(nan_val)

   ! ------------------------------------------------------------------
   ! initialise MESA infrastructure (empty dir => use $MESA_DIR)
   call const_init('', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'const_init failed')

   call math_init()

   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'chem_init failed')

   call eos_init(' ', .true., ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_init failed')

   ! read our own controls
   call read_table_controls()

   handle = alloc_eos_handle_using_inlist('inlist_mesa_tables', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'alloc_eos_handle_using_inlist failed')

   call eos_ptr(handle, rq, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'eos_ptr failed')

   helm_has_rad = rq% include_radiation

   allocate(net_iso(num_chem_isos), chem_id(max_species), d_dxa(num_eos_d_dxa_results, max_species), stat=ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'allocate failed')

   ! grid sizes
   nlogRho = nint((logRho_max - logRho_min) * pts_per_decade) + 1
   nlogT = nint((logT_max - logT_min) * pts_per_decade) + 1
   nlogRho_mix = nint((logRho_max - logRho_min) * pts_per_decade_mix) + 1
   nlogT_mix = nint((logT_max - logT_min) * pts_per_decade_mix) + 1

   call get_version_string(version_string)

   call setup_output_dirs()

   write(*,*)
   write(*,'(a)') '================ make_mesa_tables ================'
   write(*,1) 'logRho_min', logRho_min
   write(*,1) 'logRho_max', logRho_max
   write(*,1) 'logT_min', logT_min
   write(*,1) 'logT_max', logT_max
   write(*,2) 'nlogRho (iso)', nlogRho
   write(*,2) 'nlogT (iso)', nlogT
   write(*,2) 'nlogRho (mix)', nlogRho_mix
   write(*,2) 'nlogT (mix)', nlogT_mix
   write(*,'(a,l1)') ' HELM include_radiation = ', helm_has_rad
   write(*,'(a,l1)') ' Skye_use_ion_offsets   = ', rq% Skye_use_ion_offsets
   write(*,'(a)') ' NOTE: Skye ignores include_radiation; its E,S,Cv are radiation-stripped analytically.'
   write(*,'(a)') '================================================='
   write(*,*)

   ! ------------------------------------------------------------------
   ! pure-isotope source tables (skipped when regenerating only mixtures)
   if (skip_isotope_tables) then
      write(*,'(a)') ' skip_isotope_tables=T: leaving existing isotope source tables untouched.'
   else
      do i = 1, n_iso
         call write_isotope_table(i_eos_HELM, iso_names(i), dir_helm)
         call write_isotope_table(i_eos_Skye, iso_names(i), dir_skye)
      end do
   end if

   ! ------------------------------------------------------------------
   ! reference mixture tables (mix_first..mix_last; default 1..n_mix)
   write(*,'(a,i0,a,i0)') ' generating reference mixtures ', mix_first, ' .. ', mix_last
   do i = mix_first, mix_last
      call write_mixture_table(i_eos_Skye, i, dir_mix_skye)
      call write_mixture_table(i_eos_HELM, i, dir_mix_helm)
   end do

   call free_eos_handle(handle)
   call eos_shutdown

   deallocate(net_iso, chem_id, d_dxa)

   write(*,*)
   write(*,'(a)') 'make_mesa_tables: done.'

contains

   subroutine read_table_controls()
      integer :: iounit, ios
      logical :: have_file
      inquire(file='inlist_mesa_tables', exist=have_file)
      if (.not. have_file) then
         write(*,'(a)') 'inlist_mesa_tables not found; using built-in defaults.'
         return
      end if
      open(newunit=iounit, file='inlist_mesa_tables', action='read', status='old', iostat=ios)
      if (ios /= 0) then
         write(*,'(a)') 'could not open inlist_mesa_tables; using built-in defaults.'
         return
      end if
      read(iounit, nml=table_controls, iostat=ios)
      close(iounit)
      if (ios /= 0) then
         write(*,'(a)') 'could not read &table_controls; using built-in defaults.'
      end if
   end subroutine read_table_controls


   subroutine get_version_string(str)
      character(len=*), intent(out) :: str
      integer :: iounit, ios
      str = 'unknown'
      open(newunit=iounit, file=trim(mesa_dir)//'/data/version_number', &
         action='read', status='old', iostat=ios)
      if (ios /= 0) return
      read(iounit, '(a)', iostat=ios) str
      close(iounit)
      if (ios /= 0) str = 'unknown'
   end subroutine get_version_string


   subroutine setup_output_dirs()
      character(len=256) :: root
      if (len_trim(output_root) > 0) then
         root = trim(output_root)
      else
         root = trim(mesa_dir)//'/data/eosWDM_AV_data'
      end if
      dir_helm = trim(root)//'/source_helm'
      dir_skye = trim(root)//'/source_skye'
      dir_mix_skye = trim(root)//'/reference_skye_mixtures'
      dir_mix_helm = trim(root)//'/reference_helm_mixtures'
      call execute_command_line('mkdir -p '//trim(dir_helm)//' '//trim(dir_skye)//' '// &
         trim(dir_mix_skye)//' '//trim(dir_mix_helm))
   end subroutine setup_output_dirs


   ! returns whether the given component's raw output includes radiation
   logical function component_includes_radiation(which_eos)
      integer, intent(in) :: which_eos
      if (which_eos == i_eos_Skye) then
         component_includes_radiation = .true.  ! Skye always adds radiation to E,S,Cv,P_total
      else
         component_includes_radiation = helm_has_rad  ! HELM honours the control
      end if
   end function component_includes_radiation


   character(len=16) function component_name(which_eos)
      integer, intent(in) :: which_eos
      if (which_eos == i_eos_Skye) then
         component_name = 'Skye'
      else if (which_eos == i_eos_HELM) then
         component_name = 'HELM'
      else
         component_name = 'unknown'
      end if
   end function component_name


   ! ----------------------------------------------------------------
   ! pure isotope table
   subroutine write_isotope_table(which_eos, iso_name, dir)
      integer, intent(in) :: which_eos
      character(len=*), intent(in) :: iso_name, dir

      integer :: cid, Z
      real(dp) :: A
      integer :: iounit, jT, jRho
      real(dp) :: logRho, logT, Rho, T
      real(dp) :: P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta
      logical :: comp_rad
      character(len=16) :: cname

      ! sanity-check accumulators
      integer :: n_fail, n_nonmono, n_Sneg
      real(dp) :: fail_logRho_min, fail_logRho_max, fail_logT_min, fail_logT_max
      real(dp) :: prev_P
      logical :: row_has_prev
      real(dp) :: max_viol

      cid = chem_get_iso_id(iso_name)
      if (cid <= 0) then
         write(*,'(a)') 'skipping unknown isotope: '//trim(iso_name)
         return
      end if
      Z = chem_isos% Z(cid)
      A = chem_isos% W(cid)
      cname = component_name(which_eos)
      comp_rad = component_includes_radiation(which_eos)

      ! single-species composition
      net_iso(:) = 0
      net_iso(cid) = 1
      chem_id(1) = cid
      xa = 0d0
      xa(1) = 1d0

      n_fail = 0; n_nonmono = 0; n_Sneg = 0
      max_viol = 0d0
      fail_logRho_min = 1d99; fail_logRho_max = -1d99
      fail_logT_min = 1d99; fail_logT_max = -1d99

      open(newunit=iounit, file=trim(dir)//'/eosWDM_AV_'//trim(iso_name)//'.txt', &
         action='write', status='replace')

      call write_isotope_header(iounit, iso_name, Z, A, which_eos)

      do jT = 1, nlogT
         logT = logT_min + (logT_max - logT_min) * real(jT-1, dp) / real(max(nlogT-1, 1), dp)
         T = exp10(logT)
         row_has_prev = .false.
         prev_P = 0d0
         do jRho = 1, nlogRho
            logRho = logRho_min + (logRho_max - logRho_min) * real(jRho-1, dp) / real(max(nlogRho-1, 1), dp)
            Rho = exp10(logRho)

            call eos_one_node(which_eos, 1, Rho, logRho, T, logT, res, d_dlnd, d_dlnT, ierr)

            if (ierr /= 0) then
               P = nan_val; E = nan_val; F = nan_val; S = nan_val
               Zbar = nan_val; dPdRho_T = nan_val; dPdT_Rho = nan_val
               Cv = nan_val; eta = nan_val
               n_fail = n_fail + 1
               fail_logRho_min = min(fail_logRho_min, logRho)
               fail_logRho_max = max(fail_logRho_max, logRho)
               fail_logT_min = min(fail_logT_min, logT)
               fail_logT_max = max(fail_logT_max, logT)
               row_has_prev = .false.  ! only test monotonicity between adjacent valid nodes
            else
               call extract_isotope_columns(which_eos, comp_rad, Z, Rho, T, &
                  res, d_dlnd, d_dlnT, P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta)
               if (S <= 0d0) n_Sneg = n_Sneg + 1
               if (row_has_prev) then
                  if (P <= prev_P) then
                     n_nonmono = n_nonmono + 1
                     max_viol = max(max_viol, (prev_P - P) / prev_P)
                  end if
               end if
               prev_P = P
               row_has_prev = .true.
            end if

            write(iounit,'(99(es24.16e3,1x))') logRho, logT, P, E, F, S, Zbar, &
               dPdRho_T, dPdT_Rho, Cv, eta
         end do
      end do

      close(iounit)

      ! report per-table diagnostics
      write(*,'(a)') '--- '//trim(cname)//'  '//trim(iso_name)//' ---'
      write(*,'(a,i0)') '   failed nodes ', n_fail
      if (n_fail > 0) then
         write(*,'(a,4(f8.3))') '    fail bbox [logRho_min,max logT_min,max]: ', &
            fail_logRho_min, fail_logRho_max, fail_logT_min, fail_logT_max
      end if
      write(*,'(a,i0)') '   S<=0 nodes ', n_Sneg
      write(*,'(a,i0,a,es9.2)') '   non-monotonic P-in-rho steps (adjacent valid) ', n_nonmono, &
         '   max frac decrease ', max_viol
      call ideal_gas_check(which_eos, comp_rad, cid, iso_name)
      call fd_check(which_eos, comp_rad, cid, iso_name)

   end subroutine write_isotope_table


   subroutine extract_isotope_columns(which_eos, comp_rad, Z, Rho, T, &
         res, d_dlnd, d_dlnT, P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta)
      integer, intent(in) :: which_eos, Z
      logical, intent(in) :: comp_rad
      real(dp), intent(in) :: Rho, T
      real(dp), intent(in), dimension(:) :: res, d_dlnd, d_dlnT
      real(dp), intent(out) :: P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta
      real(dp) :: e_rad, s_rad, cv_rad

      ! gas pressure (always excludes radiation in HELM & Skye)
      P = exp(res(i_lnPgas))

      E = exp(res(i_lnE))
      S = exp(res(i_lnS))
      Cv = res(i_Cv)

      if (comp_rad) then
         ! strip radiation analytically (needed for Skye; for HELM only if flag on)
         e_rad = crad * pow4(T) / Rho
         s_rad = 4d0 * crad * pow3(T) / (3d0 * Rho)
         cv_rad = 4d0 * crad * pow3(T) / Rho
         E = E - e_rad
         S = S - s_rad
         Cv = Cv - cv_rad
      end if

      F = E - T * S
      Zbar = real(Z, dp)  ! full ionization assumed by both components
      eta = res(i_eta)

      ! gas-only pressure derivatives from d(ln Pgas); exact for both components.
      ! (For HELM with radiation off this equals chiRho*P/rho; for Skye chiRho is
      !  the TOTAL-pressure susceptibility, so we must use d(ln Pgas) instead.)
      dPdRho_T = P * d_dlnd(i_lnPgas) / Rho
      dPdT_Rho = P * d_dlnT(i_lnPgas) / T
   end subroutine extract_isotope_columns


   subroutine write_isotope_header(iounit, iso_name, Z, A, which_eos)
      integer, intent(in) :: iounit, Z, which_eos
      character(len=*), intent(in) :: iso_name
      real(dp), intent(in) :: A
      character(len=16) :: cname
      character(len=512) :: energy_zero, entropy_zero
      cname = component_name(which_eos)

      if (which_eos == i_eos_HELM) then
         energy_zero = 'HELM: electron/positron kinetic energy (no electron rest mass) '// &
            '+ ideal-ion thermal energy + Coulomb corrections (HELM defaults); '// &
            'no ionization offset; radiation excluded (include_radiation=F).'
         entropy_zero = 'HELM electron/positron + ideal-ion (Sackur-Tetrode) + Coulomb entropy; '// &
            'radiation entropy excluded.'
      else
         if (rq% Skye_use_ion_offsets) then
            energy_zero = 'Skye: ideal-ion + electron/positron + Coulomb free energy; '// &
               'Skye_use_ion_offsets=T adds the total electron binding (full-ionization) energy '// &
               'per gram so the fully-ionized ground state is the zero of energy (FreeEOS convention); '// &
               'radiation subtracted analytically (e_rad=crad*T^4/rho).'
         else
            energy_zero = 'Skye: ideal-ion + electron/positron + Coulomb free energy; '// &
               'Skye_use_ion_offsets=F (no ionization offset applied); '// &
               'radiation subtracted analytically (e_rad=crad*T^4/rho).'
         end if
         entropy_zero = 'Skye ideal-ion + electron/positron + Coulomb entropy; '// &
            'radiation entropy subtracted analytically (s_rad=4*crad*T^3/(3*rho)).'
      end if

      write(iounit,'(a)') '# WDM_AV source table, version 1'
      write(iounit,'(a,a)')         'isotope              ', trim(iso_name)
      write(iounit,'(a,i0)')        'Z                    ', Z
      write(iounit,'(a,f14.8)')     'A                    ', A
      write(iounit,'(a)')           'grid                 log10 uniform'
      write(iounit,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogRho              ', nlogRho, '     logRho_min  ', logRho_min, '    logRho_max  ', logRho_max
      write(iounit,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogT                ', nlogT, '     logT_min    ', logT_min, '    logT_max    ', logT_max
      write(iounit,'(a)')           'includes_radiation   F'
      write(iounit,'(a,a,a)')       'energy_zero          "', trim(energy_zero), '"'
      write(iounit,'(a,a,a)')       'entropy_zero         "', trim(entropy_zero), '"'
      write(iounit,'(a)')           'constants            "MESA const_def (CODATA 2018)"'
      write(iounit,'(a,a,a,a,a)')   'model                "MESA ', trim(cname), ' ', trim(version_string), '"'
      write(iounit,'(a)')           'tolerance            analytic'
      write(iounit,'(a)')           'columns              logRho logT P E F S Zbar dPdRho_T dPdT_Rho Cv eta'
      write(iounit,'(a)')           '# rows: logT outer loop, logRho inner loop; NaN for undefined entries'
   end subroutine write_isotope_header


   ! ----------------------------------------------------------------
   ! mixture table (full result vector + first partials)
   subroutine write_mixture_table(which_eos, imix, dir)
      integer, intent(in) :: which_eos, imix
      character(len=*), intent(in) :: dir

      integer :: nsp
      character(len=32) :: names(max_species)
      real(dp) :: fracs(max_species)
      integer :: cids(max_species)
      integer :: iounit, jT, jRho, k
      real(dp) :: logRho, logT, Rho, T
      integer :: n_fail
      character(len=16) :: cname
      character(len=256) :: comp_str, fname

      cname = component_name(which_eos)
      call get_mixture(imix, nsp, names, fracs, fname, comp_str)

      ! composition setup
      net_iso(:) = 0
      xa = 0d0
      do k = 1, nsp
         cids(k) = chem_get_iso_id(names(k))
         if (cids(k) <= 0) then
            write(*,'(a)') 'mixture skipped, unknown isotope: '//trim(names(k))
            return
         end if
         chem_id(k) = cids(k)
         net_iso(cids(k)) = k
         xa(k) = fracs(k)
      end do

      n_fail = 0

      open(newunit=iounit, file=trim(dir)//'/'//trim(fname)//'.txt', &
         action='write', status='replace')
      call write_mixture_header(iounit, which_eos, comp_str)

      do jT = 1, nlogT_mix
         logT = logT_min + (logT_max - logT_min) * real(jT-1, dp) / real(max(nlogT_mix-1, 1), dp)
         T = exp10(logT)
         do jRho = 1, nlogRho_mix
            logRho = logRho_min + (logRho_max - logRho_min) * real(jRho-1, dp) / real(max(nlogRho_mix-1, 1), dp)
            Rho = exp10(logRho)

            call eos_one_node(which_eos, nsp, Rho, logRho, T, logT, res, d_dlnd, d_dlnT, ierr)

            if (ierr /= 0) then
               n_fail = n_fail + 1
               write(iounit,'(99(es24.16e3,1x))') logRho, logT, (nan_val, k=1,24)
            else
               write(iounit,'(99(es24.16e3,1x))') logRho, logT, &
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

      close(iounit)
      write(*,'(a)') '--- '//trim(cname)//'  mixture '//trim(fname)//' ---'
      write(*,'(a,i0)') '   failed nodes ', n_fail

   end subroutine write_mixture_table


   subroutine write_mixture_header(iounit, which_eos, comp_str)
      integer, intent(in) :: iounit, which_eos
      character(len=*), intent(in) :: comp_str
      character(len=16) :: cname
      character(len=512) :: radnote
      cname = component_name(which_eos)
      if (which_eos == i_eos_Skye) then
         radnote = 'raw MESA Skye component output; Skye ignores include_radiation, '// &
            'so lnE,lnS,Cv,chiRho,chiT include radiation (lnPgas is gas-only).'
      else
         radnote = 'raw MESA HELM component output with include_radiation=F (gas-only).'
      end if
      write(iounit,'(a)') '# WDM_AV reference mixture table, version 1'
      write(iounit,'(a,a,a,a,a)') 'model                "MESA ', trim(cname), ' ', trim(version_string), '"'
      write(iounit,'(a,a,a)') 'composition          "', trim(comp_str), '"'
      write(iounit,'(a,a,a)') 'radiation_note       "', trim(radnote), '"'
      write(iounit,'(a)') 'constants            "MESA const_def (CODATA 2018)"'
      write(iounit,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogRho              ', nlogRho_mix, '     logRho_min  ', logRho_min, '    logRho_max  ', logRho_max
      write(iounit,'(a,i0,a,f10.4,a,f10.4)') &
         'nlogT                ', nlogT_mix, '     logT_min    ', logT_min, '    logT_max    ', logT_max
      write(iounit,'(a)') 'columns              logRho logT lnPgas lnE lnS chiRho chiT Cv Cp '// &
         'grad_ad gamma1 gamma3 eta mu lnfree_e dlnPgas_dlnRho dlnPgas_dlnT dlnE_dlnRho dlnE_dlnT '// &
         'dlnS_dlnRho dlnS_dlnT dCv_dlnRho dCv_dlnT dgrad_ad_dlnRho dgrad_ad_dlnT'
      write(iounit,'(a)') '# rows: logT outer loop, logRho inner loop; NaN for undefined entries'
   end subroutine write_mixture_header


   subroutine get_mixture(imix, nsp, names, fracs, fname, comp_str)
      integer, intent(in) :: imix
      integer, intent(out) :: nsp
      character(len=*), intent(out) :: names(:)
      real(dp), intent(out) :: fracs(:)
      character(len=*), intent(out) :: fname, comp_str
      names = ''
      fracs = 0d0
      select case (imix)
      case (1)
         nsp = 2
         names(1) = 'c12';  fracs(1) = 0.50d0
         names(2) = 'o16';  fracs(2) = 0.50d0
         fname = 'mix_c12_o16'
         comp_str = 'c12=0.50 o16=0.50'
      case (2)
         nsp = 3
         names(1) = 'o16';  fracs(1) = 0.55d0
         names(2) = 'ne20'; fracs(2) = 0.35d0
         names(3) = 'mg24'; fracs(3) = 0.10d0
         fname = 'mix_o16_ne20_mg24'
         comp_str = 'o16=0.55 ne20=0.35 mg24=0.10'
      case (3)
         nsp = 2
         names(1) = 'he4';  fracs(1) = 0.50d0
         names(2) = 'c12';  fracs(2) = 0.50d0
         fname = 'mix_he4_c12'
         comp_str = 'he4=0.50 c12=0.50'
      case (4)
         nsp = 2
         names(1) = 'h1';   fracs(1) = 0.70d0
         names(2) = 'he4';  fracs(2) = 0.30d0
         fname = 'mix_h1_he4'
         comp_str = 'h1=0.70 he4=0.30'
      case (5)
         nsp = 3
         names(1) = 'c12';  fracs(1) = 0.49d0
         names(2) = 'o16';  fracs(2) = 0.49d0
         names(3) = 'ne22'; fracs(3) = 0.02d0
         fname = 'mix_c12_o16_ne22'
         comp_str = 'c12=0.49 o16=0.49 ne22=0.02'
      case (6)
         nsp = 2
         names(1) = 'c12';  fracs(1) = 0.90d0
         names(2) = 'o16';  fracs(2) = 0.10d0
         fname = 'mix_c12_o16_90_10'
         comp_str = 'c12=0.90 o16=0.10'
      case (7)
         nsp = 2
         names(1) = 'c12';  fracs(1) = 0.99d0
         names(2) = 'o16';  fracs(2) = 0.01d0
         fname = 'mix_c12_o16_99_01'
         comp_str = 'c12=0.99 o16=0.01'
      case (8)
         nsp = 2
         names(1) = 'he4';  fracs(1) = 0.99d0
         names(2) = 'c12';  fracs(2) = 0.01d0
         fname = 'mix_he4_c12_99_01'
         comp_str = 'he4=0.99 c12=0.01'
      case (9)
         nsp = 2
         names(1) = 'he4';  fracs(1) = 0.10d0
         names(2) = 'c12';  fracs(2) = 0.90d0
         fname = 'mix_he4_c12_10_90'
         comp_str = 'he4=0.10 c12=0.90'
      case (10)
         nsp = 3
         names(1) = 'c12';  fracs(1) = 0.45d0
         names(2) = 'o16';  fracs(2) = 0.45d0
         names(3) = 'ne20'; fracs(3) = 0.10d0
         fname = 'mix_c12_o16_ne20'
         comp_str = 'c12=0.45 o16=0.45 ne20=0.10'
      case (11)
         nsp = 3
         names(1) = 'he4';  fracs(1) = 0.20d0
         names(2) = 'c12';  fracs(2) = 0.40d0
         names(3) = 'o16';  fracs(3) = 0.40d0
         fname = 'mix_he4_c12_o16'
         comp_str = 'he4=0.20 c12=0.40 o16=0.40'
      case default
         nsp = 0
         fname = 'mix_unknown'
         comp_str = 'unknown'
      end select
   end subroutine get_mixture


   ! ----------------------------------------------------------------
   ! single eos node (sets res/d_dlnd/d_dlnT to NaN on failure)
   subroutine eos_one_node(which_eos, nsp, Rho, logRho, T, logT, res, d_dlnd, d_dlnT, ierr)
      integer, intent(in) :: which_eos, nsp
      real(dp), intent(in) :: Rho, logRho, T, logT
      real(dp), intent(inout), dimension(:) :: res, d_dlnd, d_dlnT
      integer, intent(out) :: ierr
      ! pre-fill with NaN so that an off_table "skip" (ierr=0, res untouched)
      ! cannot leave stale values in the output.
      call set_nan(res)
      call set_nan(d_dlnd)
      call set_nan(d_dlnT)
      call eosDT_get_component( &
         handle, which_eos, nsp, chem_id, net_iso, xa, &
         Rho, logRho, T, logT, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)
      if (ierr == 0) then
         if (is_bad(res(i_lnPgas))) ierr = -1  ! skipped / off_table with no result
      end if
      if (ierr /= 0) then
         call set_nan(res)
         call set_nan(d_dlnd)
         call set_nan(d_dlnT)
      end if
   end subroutine eos_one_node


   ! ----------------------------------------------------------------
   ! sanity checks
   subroutine ideal_gas_check(which_eos, comp_rad, cid, iso_name)
      integer, intent(in) :: which_eos, cid
      logical, intent(in) :: comp_rad
      character(len=*), intent(in) :: iso_name
      real(dp) :: logRho_t, logT_t, Rho, T
      real(dp) :: P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta
      real(dp) :: bc_X, bc_Y, bc_Z, bc_abar, bc_zbar, bc_z2bar, bc_z53bar, bc_ye, bc_mass_corr, bc_sumx
      real(dp) :: P_ideal, rel_err
      integer :: Z

      ! low density, high T, non-degenerate, fully ionized
      logRho_t = -6d0
      logT_t = 7d0
      Rho = exp10(logRho_t)
      T = exp10(logT_t)
      Z = chem_isos% Z(cid)

      call eos_one_node(which_eos, 1, Rho, logRho_t, T, logT_t, res, d_dlnd, d_dlnT, ierr)
      if (ierr /= 0) then
         write(*,'(a)') '   ideal-gas check: node failed'
         return
      end if
      call extract_isotope_columns(which_eos, comp_rad, Z, Rho, T, &
         res, d_dlnd, d_dlnT, P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta)

      call basic_composition_info(1, chem_id, xa, bc_X, bc_Y, bc_Z, &
         bc_abar, bc_zbar, bc_z2bar, bc_z53bar, bc_ye, bc_mass_corr, bc_sumx)

      ! P_ideal = (ions + electrons) kT = (1+zbar) * rho * avo/abar * kT
      P_ideal = (1d0 + bc_zbar) * Rho * avo * kerg * T / bc_abar
      rel_err = abs(P / P_ideal - 1d0)
      write(*,'(a,es12.4,a,es12.4,a,es10.2)') &
         '   ideal-gas check @ logRho=-6 logT=7: P=', P, ' P_ideal=', P_ideal, ' rel_err=', rel_err
   end subroutine ideal_gas_check


   subroutine fd_check(which_eos, comp_rad, cid, iso_name)
      integer, intent(in) :: which_eos, cid
      logical, intent(in) :: comp_rad
      character(len=*), intent(in) :: iso_name
      integer, parameter :: npt = 3
      real(dp) :: test_logRho(npt), test_logT(npt)
      integer :: ip, Z
      real(dp) :: logRho, logT, Rho, T, dlnRho
      real(dp) :: Pm, Pp, Rm, Rp
      real(dp) :: P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta
      real(dp) :: fd, rel

      test_logRho = [-4d0, 0d0, 4d0]
      test_logT = [5d0, 6d0, 7d0]
      Z = chem_isos% Z(cid)
      dlnRho = 1d-4

      do ip = 1, npt
         logRho = test_logRho(ip)
         logT = test_logT(ip)
         Rho = exp10(logRho)
         T = exp10(logT)

         call eos_one_node(which_eos, 1, Rho, logRho, T, logT, res, d_dlnd, d_dlnT, ierr)
         if (ierr /= 0) cycle
         call extract_isotope_columns(which_eos, comp_rad, Z, Rho, T, &
            res, d_dlnd, d_dlnT, P, E, F, S, Zbar, dPdRho_T, dPdT_Rho, Cv, eta)

         ! central difference of gas pressure w.r.t. rho at fixed T
         Rp = Rho * exp(dlnRho)
         Rm = Rho * exp(-dlnRho)
         call eos_one_node(which_eos, 1, Rp, log10(Rp), T, logT, res, d_dlnd, d_dlnT, ierr)
         if (ierr /= 0) cycle
         Pp = exp(res(i_lnPgas))
         call eos_one_node(which_eos, 1, Rm, log10(Rm), T, logT, res, d_dlnd, d_dlnT, ierr)
         if (ierr /= 0) cycle
         Pm = exp(res(i_lnPgas))

         fd = (Pp - Pm) / (Rp - Rm)
         rel = abs(dPdRho_T / fd - 1d0)
         write(*,'(a,f6.1,a,f6.1,a,es10.2,a,es10.2,a,es9.2)') &
            '   FD dP/drho @ logRho=', logRho, ' logT=', logT, &
            ': analytic=', dPdRho_T, ' fd=', fd, ' rel=', rel
      end do
   end subroutine fd_check

end program make_mesa_tables
