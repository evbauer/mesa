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

module eos_wdm_av_eval
   ! eosDT-cascade wrapper for the WDM_AV (warm dense matter, additive volume)
   ! component (WP5).  Mirrors skye.f90's get_Skye_for_eosdt / Get_Skye_alfa_simple:
   !
   !   get_WDM_AV_for_eosdt : build the per-species table list for the mixture from
   !                          the controls, call eos_wdm_av_get, and do the frac
   !                          bookkeeping.
   !   Get_WDM_AV_alfa      : simple Rho/T ramp blend (trapezoidal windows between
   !                          the any/all limits), alfa = 1 outside the box (not
   !                          used), 0 inside (fully used).  Do_Blend applies the
   !                          quintic smoothstep.

   use const_def, only: dp, mesa_data_dir, qe, avo, kerg, four_thirds_pi, one_third
   use eos_def
   use eos_wdm_av_def, only: eos_wdm_av_table_index
   use eos_wdm_av_load, only: eos_wdm_av_load_isotope
   use eos_wdm_av_mix, only: eos_wdm_av_options, eos_wdm_av_state, eos_wdm_av_get, &
                             wdm_av_rule_pressure, wdm_av_rule_electron_density
   use chem_def, only: chem_isos, num_chem_isos, iso_name_length
   use auto_diff
   use math_lib
   use utils_lib, only: mesa_error

   implicit none

   private
   public :: get_WDM_AV_for_eosdt, Get_WDM_AV_alfa

contains

   subroutine get_WDM_AV_for_eosdt(handle, dbg, Z, X, abar, zbar, species, chem_id, net_iso, xa, &
                                   rho, logRho, T, logT, remaining_fraction, &
                                   res, d_dlnd, d_dlnT, d_dxa, skip, ierr)
      integer, intent(in) :: handle
      logical, intent(in) :: dbg
      real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
      integer, intent(in) :: species
      integer, pointer :: chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: rho, logRho, T, logT
      real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
      real(dp), intent(inout), dimension(nv, species) :: d_dxa
      logical, intent(out) :: skip
      integer, intent(out) :: ierr

      type (EoS_General_Info), pointer :: rq
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: state
      character(len=512) :: table_dir
      character(len=iso_name_length) :: iso_name
      integer :: nsp, j, idx, cid
      integer :: tbl_idx(species)
      real(dp) :: A_sp(species), Z_sp(species), X_sp(species)
      real(dp) :: covered_mass, uncovered_mass, mass_limit, Xsum

      include 'formats'

      ierr = 0
      skip = .false.
      rq => eos_handles(handle)

      call wdm_av_resolve_table_dir(rq, table_dir)
      mass_limit = rq% mass_fraction_limit_for_WDM_AV

      ! build the per-species table list
      nsp = 0
      covered_mass = 0d0
      uncovered_mass = 0d0
      do j = 1, species
         if (xa(j) <= mass_limit) cycle   ! trace species dropped and X renormalised
         cid = chem_id(j)
         iso_name = chem_isos% name(cid)

         idx = 0
         if (wdm_av_file_exists(table_dir, iso_name)) then
            call eos_wdm_av_load_isotope(trim(table_dir), trim(iso_name), idx, ierr)
            if (ierr /= 0) then
               idx = 0
               ierr = 0
            end if
         end if

         if (idx == 0 .and. rq% WDM_AV_rescale_missing_isotopes) then
            call wdm_av_find_same_Z(table_dir, chem_isos% Z(cid), chem_isos% W(cid), idx)
         end if

         if (idx > 0) then
            nsp = nsp + 1
            tbl_idx(nsp) = idx
            A_sp(nsp) = chem_isos% W(cid)
            Z_sp(nsp) = real(chem_isos% Z(cid), dp)
            X_sp(nsp) = xa(j)
            covered_mass = covered_mass + xa(j)
         else
            uncovered_mass = uncovered_mass + xa(j)
         end if
      end do

      ! a non-trivial uncovered mass means the tables cannot represent the mixture
      if (uncovered_mass > mass_limit .or. nsp < 1) then
         skip = .true.
         return
      end if

      ! renormalise mass fractions over the covered species
      Xsum = covered_mass
      if (Xsum <= 0d0) then
         skip = .true.
         return
      end if
      do j = 1, nsp
         X_sp(j) = X_sp(j)/Xsum
      end do

      ! fill solver options from the controls
      opts% newton_tol = rq% WDM_AV_newton_tol
      opts% max_iter = rq% WDM_AV_max_iter
      opts% mixing_entropy = wdm_av_mixing_entropy_code(rq% WDM_AV_mixing_entropy)
      opts% mixing_rule = wdm_av_mixing_rule_code(rq% WDM_AV_mixing_rule)
      opts% debug_crosscheck = rq% WDM_AV_debug

      call eos_wdm_av_get(opts, T, rho, nsp, tbl_idx(1:nsp), A_sp(1:nsp), Z_sp(1:nsp), X_sp(1:nsp), &
                          res, d_dlnd, d_dlnT, state, ierr)

      if (ierr /= 0) then
         if (rq% WDM_AV_debug) then
            if (dbg) write(*,2) 'WDM_AV get failed, ierr', ierr
            return
         end if
         ! consistent with okay_to_convert_ierr_to_skip: decline the point
         skip = .true.
         ierr = 0
         return
      end if

      ! composition derivatives not provided (Decision D17)
      d_dxa = 0d0

      ! zero all frac components, then mark this one
      res(i_frac:i_frac+num_eos_frac_results-1) = 0d0
      d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0d0
      d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0d0
      res(i_frac_WDM_AV) = 1d0

   end subroutine get_WDM_AV_for_eosdt


   subroutine Get_WDM_AV_alfa(rq, logRho, logT, Z, abar, zbar, alfa, d_alfa_dlogT, d_alfa_dlogRho, ierr)
      ! Simple Rho/T ramp blend in the style of Get_Skye_alfa_simple.  alfa is the
      ! fraction handed to the lower cascade levels: 1 outside the "any" box (not
      ! used), 0 inside the "all" box (fully used), linear ramps between.
      !
      ! When Gamma_max_for_any_WDM_AV > 0 the in-box weight is additionally
      ! multiplied by a window in log10(Gamma), where Gamma is the mixture ion
      ! coupling parameter under full ionization computed as the PC blend does
      ! (Gamma = zbar**(5/3) * Gamma_e; see Get_PC_alfa / pc_eos GAMImean=Z53*GAME).
      ! The window is 1 for Gamma <= Gamma_max_for_all_WDM_AV, ramps to 0 at
      ! Gamma_max_for_any_WDM_AV, and 0 above, keeping WDM_AV clear of the
      ! crystallized regime and handing over to Skye before the melting line (D23).
      ! A value <= 0 for Gamma_max_for_any_WDM_AV disables the limit and reproduces
      ! the present behaviour bit for bit.
      type (EoS_General_Info), pointer :: rq
      real(dp), intent(in) :: logRho, logT, Z, abar, zbar
      real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
      integer, intent(out) :: ierr

      type(auto_diff_real_2var_order1) :: logRho_auto, logT_auto
      type(auto_diff_real_2var_order1) :: wRho, wT, wGamma, logGamma_auto, blend
      real(dp) :: logGamma0

      ierr = 0

      ! logRho is val1
      logRho_auto% val = logRho
      logRho_auto% d1val1 = 1d0
      logRho_auto% d1val2 = 0d0

      ! logT is val2
      logT_auto% val = logT
      logT_auto% d1val1 = 0d0
      logT_auto% d1val2 = 1d0

      wRho = wdm_av_window(logRho_auto, rq% logRho_min_for_any_WDM_AV, rq% logRho_min_for_all_WDM_AV, &
                           rq% logRho_max_for_all_WDM_AV, rq% logRho_max_for_any_WDM_AV)
      wT = wdm_av_window(logT_auto, rq% logT_min_for_any_WDM_AV, rq% logT_min_for_all_WDM_AV, &
                         rq% logT_max_for_all_WDM_AV, rq% logT_max_for_any_WDM_AV)

      ! in-box fraction wRho*wT, optionally cut by the coupling-parameter window
      blend = wRho*wT

      if (rq% Gamma_max_for_any_WDM_AV > 0d0) then
         ! mixture ion coupling parameter under full ionization, PC convention:
         !   Gamma_e = Ge0 * rho**(1/3) / T,  Ge0 = qe^2 (4/3 pi avo zbar/abar)^(1/3)/kerg
         !   Gamma   = zbar**(5/3) * Gamma_e
         ! => log10 Gamma = logGamma0 + logRho/3 - logT, with
         !    d logGamma/d logRho = 1/3,  d logGamma/d logT = -1 at fixed composition.
         logGamma0 = (5d0/3d0)*log10(zbar) &
                     + log10(qe*qe*pow(four_thirds_pi*avo*zbar/abar, one_third)/kerg)
         logGamma_auto% val = logGamma0 + logRho/3d0 - logT
         logGamma_auto% d1val1 = one_third   ! d/dlogRho
         logGamma_auto% d1val2 = -1d0        ! d/dlogT
         ! window: 1 for log Gamma <= log(Gamma_all), ramp to 0 at log(Gamma_any);
         ! lower edges set far below so only the upper ramp acts.
         wGamma = wdm_av_window(logGamma_auto, -1d99, -1d98, &
                                log10(rq% Gamma_max_for_all_WDM_AV), &
                                log10(rq% Gamma_max_for_any_WDM_AV))
         blend = blend*wGamma
      end if

      ! alfa handed to lower levels is the complement of the in-box fraction
      blend = 1d0 - blend

      alfa = blend% val
      d_alfa_dlogRho = blend% d1val1
      d_alfa_dlogT = blend% d1val2

   end subroutine Get_WDM_AV_alfa


   ! trapezoidal window: 0 for x < a_lo, ramps to 1 at b_lo, 1 on [b_lo,b_hi],
   ! ramps back to 0 at a_hi, 0 for x > a_hi
   function wdm_av_window(x, a_lo, b_lo, b_hi, a_hi) result(w)
      type(auto_diff_real_2var_order1), intent(in) :: x
      real(dp), intent(in) :: a_lo, b_lo, b_hi, a_hi
      type(auto_diff_real_2var_order1) :: w

      if (x < a_lo) then
         w = 0d0
      else if (x <= b_lo) then
         w = (x - a_lo)/(b_lo - a_lo)
      else if (x <= b_hi) then
         w = 1d0
      else if (x <= a_hi) then
         w = (a_hi - x)/(a_hi - b_hi)
      else
         w = 0d0
      end if
   end function wdm_av_window


   ! directory holding the per-isotope tables, including the table-set subdir
   subroutine wdm_av_resolve_table_dir(rq, table_dir)
      type (EoS_General_Info), pointer :: rq
      character(len=*), intent(out) :: table_dir
      if (len_trim(rq% WDM_AV_table_dir) > 0) then
         table_dir = trim(rq% WDM_AV_table_dir)
      else
         table_dir = trim(mesa_data_dir)//'/eosWDM_AV_data'
      end if
      if (len_trim(rq% WDM_AV_table_set) > 0) then
         table_dir = trim(table_dir)//'/'//trim(rq% WDM_AV_table_set)
      end if
   end subroutine wdm_av_resolve_table_dir


   ! true iff the processed HDF5 table file for isotope exists (checked before
   ! calling the loader so missing-file diagnostics are not printed)
   logical function wdm_av_file_exists(table_dir, isotope) result(ex)
      character(len=*), intent(in) :: table_dir, isotope
      character(len=512) :: fname
      fname = trim(table_dir)//'/eosWDM_AV_'//trim(isotope)//'.h5'
      inquire(file=trim(fname), exist=ex)
   end function wdm_av_file_exists


   ! same-Z nearest-A fallback (Decision D6): find the loadable table with Z == Zint
   ! whose tabulated atomic weight is closest to A_target, and load it.  Returns
   ! its store index, or 0 if none is available.
   subroutine wdm_av_find_same_Z(table_dir, Zint, A_target, idx)
      character(len=*), intent(in) :: table_dir
      integer, intent(in) :: Zint
      real(dp), intent(in) :: A_target
      integer, intent(out) :: idx

      character(len=iso_name_length) :: best_name, nm
      real(dp) :: best_dA, dA
      integer :: k, ierr_load

      idx = 0
      best_dA = huge(1d0)
      best_name = ''
      do k = 1, num_chem_isos
         if (chem_isos% Z(k) /= Zint) cycle
         nm = chem_isos% name(k)
         if (.not. wdm_av_file_exists(table_dir, nm)) cycle
         dA = abs(chem_isos% W(k) - A_target)
         if (dA < best_dA) then
            best_dA = dA
            best_name = nm
         end if
      end do

      if (len_trim(best_name) == 0) return
      call eos_wdm_av_load_isotope(trim(table_dir), trim(best_name), idx, ierr_load)
      if (ierr_load /= 0) idx = 0
   end subroutine wdm_av_find_same_Z


   ! map the mixing-entropy control string to the integer code used by eos_wdm_av_mix
   integer function wdm_av_mixing_entropy_code(mode) result(code)
      character(len=*), intent(in) :: mode
      select case (trim(adjustl(mode)))
      case ('none')
         code = 0
      case ('ions')
         code = 1
      case ('ions_electrons')
         code = 2
      case default
         code = 2
      end select
   end function wdm_av_mixing_entropy_code


   integer function wdm_av_mixing_rule_code(mode) result(code)
      ! Decision D27: 'pressure' (additive volume, default) | 'electron_density'
      character(len=*), intent(in) :: mode
      select case (trim(adjustl(mode)))
      case ('pressure')
         code = wdm_av_rule_pressure
      case ('electron_density')
         code = wdm_av_rule_electron_density
      case default
         code = wdm_av_rule_pressure
      end select
   end function wdm_av_mixing_rule_code

end module eos_wdm_av_eval
