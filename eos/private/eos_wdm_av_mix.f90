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

module eos_wdm_av_mix
   ! Additive-volume (AV) mixture equation of state (WP3).  Given per-species
   ! (rho,T) free-energy tables and a composition, this module returns the full
   ! MESA eosDT result vector and its first partials, via
   !
   !   (1) a deterministic real Newton solve of the AV constraint
   !       P_i(rho_i,T) = P* for all i, sum_i X_i/rho_i = 1/rho  (Sec. 4.2 notes);
   !   (2) a derivative pass of TWO flattened Newton steps in auto_diff arithmetic
   !       seeded at the converged real solution, which delivers the partials
   !       through third order (order-doubling lemma, Sec. 4.3 notes);
   !   (3) the mixture free energy F_mix = sum_i X_i f_i - T*Delta s_mix with the
   !       ideal mixing entropy evaluated with the ACTUAL AV partial volumes
   !       (Decision D8, Sec. 5 notes);
   !   (4) electron auxiliaries under full ionization (Decision D9);
   !   (5) Skye's pack_for_export to build the result vector (radiation added
   !       there analytically).
   !
   ! The converged state is returned in eos_wdm_av_state for a future composition
   ! derivative routine (Decision D17); d_dxa is not computed here.

   use const_def, only: dp, kerg, amu, avo, pi, me, planck_h
   use eos_def
   use eos_wdm_av_def
   use eos_wdm_av_interp, only: eos_wdm_av_eval_real, &
                                eos_wdm_av_eval_species_real, eos_wdm_av_eval_species_ad
   use skye_thermodynamics, only: pack_for_export, thermodynamics_from_free_energy, &
                                  compute_derived_quantities
   use skye_ideal, only: compute_ideal_ele, compute_xne
   use auto_diff
   use math_lib
   use utils_lib, only: is_bad

   implicit none
   private

   public :: eos_wdm_av_options
   public :: eos_wdm_av_state
   public :: eos_wdm_av_get
   public :: wdm_av_ierr_no_converge
   public :: wdm_av_ierr_bad_domain
   public :: wdm_av_rule_pressure, wdm_av_rule_electron_density

   ! mix-level error codes (distinct from the def-level codes)
   integer, parameter :: wdm_av_ierr_no_converge = -111
   integer, parameter :: wdm_av_ierr_bad_domain  = -112

   ! mixing-entropy modes
   integer, parameter :: wdm_av_mix_none = 0
   integer, parameter :: wdm_av_mix_ions = 1
   integer, parameter :: wdm_av_mix_ions_electrons = 2

   ! mixing rule (Decision D27): how the species' partial densities are chosen
   integer, parameter :: wdm_av_rule_pressure = 0          ! additive volume at common P, T
   integer, parameter :: wdm_av_rule_electron_density = 1  ! common mean electron density, rho_i = rho (Z_i/A_i)/y_e

   ! run-time options
   type eos_wdm_av_options
      real(dp) :: newton_tol = 1d-12
      integer :: max_iter = 50
      integer :: mixing_entropy = wdm_av_mix_ions_electrons   ! 0 none, 1 ions, 2 ions+electrons
      integer :: mixing_rule = wdm_av_rule_pressure            ! 0 common pressure (AV), 1 common electron density
      logical :: debug_crosscheck = .false.
      logical :: include_radiation = .true.    ! for the tests; pack always adds radiation
      integer :: deriv_steps = 2               ! hidden: # auto_diff Newton steps (default 2)
   end type eos_wdm_av_options

   ! converged AV state (hook for composition derivatives, Decision D17)
   type eos_wdm_av_state
      real(dp) :: P = 0d0, lnP = 0d0
      integer :: nsp = 0
      real(dp), allocatable :: rho_i(:), xi_i(:), f_i(:), g_i(:), h_i(:), e_i(:), &
                               s_i(:), v_i(:), P_rho_i(:), P_T_i(:), chi_rho_i(:)
      real(dp) :: G_PP = 0d0, G_PT = 0d0, G_TT = 0d0
      real(dp) :: P_rho_mix = 0d0, P_T_mix = 0d0, cv_mix = 0d0, cp_mix = 0d0
      integer :: n_iter = 0
      logical :: converged = .false.
      real(dp) :: delta_s_mix = 0d0
      real(dp) :: residual = 0d0
   end type eos_wdm_av_state

contains

   !----------------------------------------------------------------------------
   ! public entry
   !----------------------------------------------------------------------------

   subroutine eos_wdm_av_get(opts, T, rho, nsp, tbl_idx, A, Z, X, res, d_dlnd, d_dlnT, state, ierr)
      type(eos_wdm_av_options), intent(in) :: opts
      real(dp), intent(in) :: T, rho
      integer, intent(in) :: nsp
      integer, intent(in) :: tbl_idx(nsp)
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp)
      real(dp), intent(out) :: res(nv), d_dlnd(nv), d_dlnT(nv)
      type(eos_wdm_av_state), intent(out) :: state
      integer, intent(out) :: ierr

      real(dp) :: Xn(nsp), xi(nsp), chi(nsp), lnc(nsp)
      real(dp) :: tau, v_target, L, Pstar, resid, Xsum
      integer :: i, n_iter

      ! auto_diff workspace
      type(auto_diff_real_2var_order3) :: temp, den, tau_ad, v_ad, L_ad
      type(auto_diff_real_2var_order3) :: xi_ad(nsp), f_ad(nsp), rho_ad(nsp)
      type(auto_diff_real_2var_order3) :: s_fx, s_ft, s_fxx
      type(auto_diff_real_2var_order3) :: F_AV, F_mix, ds_mix
      type(auto_diff_real_2var_order3) :: xnefer, etaele
      type(auto_diff_real_2var_order3) :: F_zero, phase, latent_ddlnT, latent_ddlnRho

      real(dp) :: ye, abar, zbar, ytot1

      ierr = 0
      res = 0d0; d_dlnd = 0d0; d_dlnT = 0d0
      Pstar = 0d0; L = 0d0   ! set by av_solve_real (AV path); unused by the ne path

      if (nsp < 1) then
         ierr = wdm_av_ierr_bad_domain
         return
      end if

      ! defensive normalisation of the mass fractions
      Xsum = sum(X(1:nsp))
      if (Xsum <= 0d0) then
         ierr = wdm_av_ierr_bad_domain
         return
      end if
      do i = 1, nsp
         Xn(i) = X(i)/Xsum
      end do

      tau = log(T)
      v_target = 1d0/rho

      ! -- (1) choose the species partial densities rho_i (real) --------------
      ! The mixing rule decides how the per-species densities are placed; the
      ! species evaluation, mixing entropy, electron auxiliaries and the pack to
      ! MESA results below are shared between the two rules.
      select case (opts% mixing_rule)

      case (wdm_av_rule_electron_density)
         ! Common mean electron density (Decision D27).  There is NO solve: for a
         ! fully ionized plasma the species sit at a common electron density n_e,
         ! so species i occupies rho_i = n_e A_i m_u / Z_i = rho * y_e * A_i / Z_i
         ! with y_e = sum_i X_i Z_i/A_i.  These are fixed multiples of rho, c_i rho,
         ! so (a) sum_i X_i/rho_i = 1/rho identically (no volume equation) and
         ! (b) every n_e,i = Z_i rho_i/(A_i m_u) equals n_e (the electron mixing
         ! term vanishes identically).
         ! NOTE on the formula: rho_i scales as A_i/Z_i (NOT Z_i/A_i).  The share
         ! of specific volume is V_i/V = X_i rho/rho_i = (X_i Z_i/A_i)/y_e, i.e.
         ! the species' share of the electrons, as it must be for one common n_e.
         ! HOOK: the Zbar-dependent generalisation of D27 (partial ionization)
         ! would replace Z_i by a tabulated Zbar_i(rho_i,T) here and reinstate a
         ! Newton solve of the structure of av_solve_real (ln n_e in place of
         ! ln P*); not implemented -- D27 is open and no Zbar tables exist.
         call ne_densities_real(nsp, A, Z, Xn, rho, xi, lnc, ierr)
         if (ierr /= 0) return
         state% n_iter = 0
         state% residual = 0d0
         state% converged = .true.

      case default
         ! -- additive volume at common pressure (the original path) -----------
         call av_solve_real(opts, tau, v_target, nsp, tbl_idx, A, Z, Xn, L, xi, chi, &
                            Pstar, n_iter, resid, ierr)
         state% n_iter = n_iter
         state% residual = resid
         state% converged = (ierr == 0)
         if (ierr /= 0) return

      end select

      ! validity-mask gate (Skye tables flag unusable cold-corner cells); shared
      do i = 1, nsp
         if (.not. cell_usable(wdm_av_tables(tbl_idx(i)), A(i), xi(i), tau)) then
            ierr = wdm_av_ierr_out_of_table
            return
         end if
      end do

      ! -- fill the converged-state block (real) ------------------------------
      ! For AV this also builds the hand-derived Gibbs cross-check (G_PP,...); for
      ! the common-n_e rule the pressures are NOT equal, so the Gibbs cross-check
      ! is meaningless and skipped -- fill_state_ne records the mixture P (the
      ! volume-weighted sum) and the per-species block only.
      select case (opts% mixing_rule)
      case (wdm_av_rule_electron_density)
         call fill_state_ne(nsp, tbl_idx, A, Xn, xi, tau, T, rho, state, ierr)
      case default
         call fill_state(nsp, tbl_idx, A, Xn, xi, tau, T, rho, Pstar, L, state, ierr)
      end select
      if (ierr /= 0) return

      ! -- (2) auto_diff derivative seeds (shared) ----------------------------
      temp = T
      temp% d1val1 = 1d0
      den = rho
      den% d1val2 = 1d0
      tau_ad = log(temp)
      v_ad = 1d0/den

      ! -- obtain xi_ad(i) = ln rho_i as auto_diff, by the mixing rule --------
      select case (opts% mixing_rule)

      case (wdm_av_rule_electron_density)
         ! rho_i = c_i rho with constant c_i, so ln rho_i = ln rho + ln c_i and
         ! d(ln rho_i)/d(ln rho) = 1 exactly: a single auto_diff evaluation, no
         ! Newton steps (the order-doubling pass of the AV path is not needed).
         do i = 1, nsp
            xi_ad(i) = log(den) + lnc(i)
         end do

      case default
         L_ad = L                      ! seed: constant at converged value
         do i = 1, nsp
            xi_ad(i) = xi(i)           ! seed: constant at converged value
         end do
         call av_newton_pass_ad(opts% deriv_steps, tau_ad, v_ad, nsp, tbl_idx, A, Xn, L_ad, xi_ad, ierr)
         if (ierr /= 0) return

      end select

      ! final species free energies and partial densities at the converged xi_ad
      do i = 1, nsp
         call eos_wdm_av_eval_species_ad(wdm_av_tables(tbl_idx(i)), A(i), xi_ad(i), tau_ad, &
                                         f_ad(i), s_fx, s_ft, s_fxx, ierr)
         if (ierr /= 0) return
         rho_ad(i) = exp(xi_ad(i))
      end do

      ! -- (3) mixture free energy + mixing term ------------------------------
      F_AV = 0d0
      do i = 1, nsp
         F_AV = F_AV + Xn(i)*f_ad(i)
      end do

      call mixing_entropy_ad(opts% mixing_entropy, nsp, A, Z, Xn, den, rho_ad, ds_mix)
      F_mix = F_AV - temp*ds_mix
      state% delta_s_mix = ds_mix% val

      ! -- (4) electron auxiliaries under full ionization (Decision D9) -------
      ytot1 = sum(Xn(1:nsp)/A(1:nsp))         ! = 1/abar
      abar = 1d0/ytot1
      ye = sum(Xn(1:nsp)*Z(1:nsp)/A(1:nsp))
      zbar = abar*ye

      xnefer = compute_xne(den, ytot1, zbar)
      call electron_eta(temp, den, ytot1, ye, zbar, xnefer, etaele)

      ! -- (5) pack into the MESA result vector -------------------------------
      F_zero = 0d0
      phase = 0d0
      latent_ddlnT = 0d0
      latent_ddlnRho = 0d0

      call pack_for_export(F_mix, F_zero, F_zero, F_zero, temp, den, xnefer, etaele, &
                           abar, zbar, phase, latent_ddlnT, latent_ddlnRho, &
                           res, d_dlnd, d_dlnT, ierr)
      if (ierr /= 0) return

      ! res(i_mu) is already abar/(1+zbar) from pack_for_export (full ionization).
      ! The eosdt wrapper will set the i_frac_WDM_AV slot later; leave the
      ! fraction block (i_frac : i_frac+num_eos_frac_results-1) at zero here.

      ! -- (6) optional debug cross-check (hand-derived Gibbs vs auto_diff) ----
      ! The Gibbs cross-check assumes a common pressure (G_PP,... are built from
      ! the equal-pressure solve), so it is only meaningful for the AV rule.
      if (opts% debug_crosscheck .and. opts% mixing_rule == wdm_av_rule_pressure) then
         call crosscheck_print(F_AV, temp, den, state)
      end if

   end subroutine eos_wdm_av_get


   !----------------------------------------------------------------------------
   ! (1) real AV solve: bracketed flattened Newton on L = ln P*, with per-species
   ! 1-D inversions ln P_i(xi_i)=L (keeps r_i ~ 0 so the bracket is exact).
   !----------------------------------------------------------------------------

   subroutine av_solve_real(opts, tau, v_target, nsp, tbl_idx, A, Z, X, L, xi, chi, &
                            Pstar, n_iter, resid, ierr)
      type(eos_wdm_av_options), intent(in) :: opts
      real(dp), intent(in) :: tau, v_target
      integer, intent(in) :: nsp, tbl_idx(nsp)
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp)
      real(dp), intent(out) :: L, xi(nsp), chi(nsp), Pstar, resid
      integer, intent(out) :: n_iter, ierr

      real(dp) :: Tk, rho, P0, sumB, L_lo, L_hi, V, S, dL, L_try, V_try
      real(dp) :: xi_try(nsp), chi_try(nsp)
      real(dp) :: tol, logv, prev_resid
      integer :: i, it, nexp, n_stall
      logical :: have_lo, have_hi, converged

      ! Convergence floors for the outer AV solve (documented at the loop below).
      ! loose_tol: the residual below which the constraint is "satisfied to the
      !   attainable round-off floor" on real (interpolated) tables.  The inner
      !   1-D inversions converge each r_i = L - ln P_i only to ~1e-11, which
      !   caps |ln V - ln v| (and the Newton step |dL|) at a few 1e-12 -- above
      !   newton_tol = 1e-12.  loose_tol is set well above that floor.
      ! step_floor: a Newton step in L this small changes nothing (full
      !   convergence; analytic tables reach this via newton_tol first).
      real(dp), parameter :: loose_tol = 1d-9
      real(dp), parameter :: step_floor = 1d-13
      integer, parameter :: n_stall_accept = 3

      ierr = 0
      n_iter = 0
      resid = 0d0
      Pstar = 0d0
      tol = opts% newton_tol
      Tk = exp(tau)
      rho = 1d0/v_target
      logv = log(v_target)

      ! ideal fully ionized initial pressure
      sumB = 0d0
      do i = 1, nsp
         sumB = sumB + X(i)*(1d0 + Z(i))/A(i)
      end do
      P0 = rho*kerg*Tk*sumB/amu
      if (P0 <= 0d0) then
         ierr = wdm_av_ierr_bad_domain
         return
      end if
      L = log(P0)

      ! per-species initial density guess (ideal) then invert ln P_i = L
      do i = 1, nsp
         xi(i) = log(max(1d-300, P0*A(i)*amu/((1d0 + Z(i))*kerg*Tk)))
      end do
      call av_volume_at_L(L, tau, nsp, tbl_idx, A, X, xi, chi, V, ierr)
      if (ierr /= 0) return

      ! single species: rho_1 = rho exactly; evaluate L and chi directly
      if (nsp == 1) then
         xi(1) = -logv                 ! ln rho
         call eval_lnP_chi(wdm_av_tables(tbl_idx(1)), A(1), xi(1), tau, L, chi(1), ierr)
         if (ierr /= 0) return
         Pstar = exp(L)
         n_iter = 1
         resid = 0d0
         return
      end if

      ! Safeguarded flattened Newton on L = ln P* with an on-the-fly bracket
      ! (V is strictly decreasing in L).  Three ways to accept a solution:
      !   (A) resid = |ln V - ln v| < newton_tol    -- full convergence; the
      !       analytic/ideal tables (r_i ~ 0 to round-off) always reach this.
      !   (B) |dL| < step_floor                      -- the Newton update to L
      !       (and hence every xi_i = (r_i + dL)/chi_i) is negligible.
      !   (C) resid < loose_tol and it has stopped decreasing for n_stall_accept
      !       consecutive steps -- the real-table round-off floor.  On
      !       interpolated tables the inner inversions leave r_i ~ 1e-11, so
      !       resid and |dL| bottom out at a few 1e-12 and the plain
      !       resid < newton_tol / |dL| < newton_tol tests of the previous
      !       version were unreachable: the iteration limit-cycled to max_iter
      !       and returned no_converge.  Accepting the stalled floor fixes that.
      ! The test is deterministic (no state saved across calls, Decision D15);
      ! the bracket + bisection safeguard plus the stall counter guarantee the
      ! loop terminates and cannot limit-cycle.
      have_lo = .false.; have_hi = .false.
      L_lo = L; L_hi = L
      xi_try = xi; chi_try = chi
      converged = .false.
      prev_resid = huge(1d0)
      n_stall = 0
      do it = 1, opts% max_iter
         n_iter = it
         resid = abs(log(V) - logv)
         if (V > v_target) then
            L_lo = L; have_lo = .true.
         else
            L_hi = L; have_hi = .true.
         end if
         ! (A) full convergence
         if (resid < tol) then
            converged = .true.; exit
         end if

         ! flattened Newton step (r_i ~ 0 after inner inversion):
         !   delta_L = (V - v)/S, S = sum_i w_i/chi_i, w_i = X_i/rho_i
         S = 0d0
         do i = 1, nsp
            S = S + X(i)*exp(-xi(i))/chi(i)
         end do
         if (S <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return
         end if
         dL = (V - v_target)/S

         ! (B) negligible Newton step: nothing left to improve
         if (abs(dL) < step_floor) then
            converged = .true.; exit
         end if

         ! (C) stalled iteration.  The flattened Newton converges quadratically,
         ! so n_stall_accept consecutive steps with no residual decrease mean it
         ! has hit a floor.  Two outcomes:
         !   - resid below loose_tol: the real-table round-off floor (the inner
         !     1-D inversions pin r_i at ~1e-11, which caps |lnV - lnv| and |dL|
         !     at a few 1e-12).  Accept the solution.
         !   - resid still large with a one-sided bracket (V never crossed v):
         !     the volume constraint cannot be met within the tabulated density
         !     range -- a partial density is pinned at a table edge (e.g. the
         !     logRho=10 corner).  Decline as out_of_table (D10) instead of
         !     running to max_iter and returning a spurious no_converge.
         if (resid >= 0.9d0*prev_resid) then
            n_stall = n_stall + 1
         else
            n_stall = 0
         end if
         if (n_stall >= n_stall_accept) then
            if (resid < loose_tol) then
               converged = .true.; exit
            else
               ! stalled well above the floor: either a one-sided bracket (the
               ! root needs a density beyond the table, e.g. the logRho=10
               ! corner) or a collapsed two-sided bracket straddling a jump in
               ! V(L) (a non-monotonic inner inversion, so the AV constraint has
               ! no solution in the attainable range).  Either way no usable AV
               ! state exists here -> decline as out_of_table (D10) rather than
               ! iterating to max_iter and returning no_converge.
               ierr = wdm_av_ierr_out_of_table; return
            end if
         end if

         L_try = L + dL
         if (have_lo .and. have_hi) then
            if (L_try <= L_lo .or. L_try >= L_hi) L_try = 0.5d0*(L_lo + L_hi)
         end if

         ! commit the step, backtracking if it lands outside any species' table
         nexp = 0
         do
            call av_volume_at_L(L_try, tau, nsp, tbl_idx, A, X, xi_try, chi_try, V_try, ierr)
            if (ierr == 0) exit
            nexp = nexp + 1
            if (nexp > 60) then
               ierr = wdm_av_ierr_out_of_table; return
            end if
            L_try = 0.5d0*(L + L_try)
            xi_try = xi
         end do
         prev_resid = resid
         L = L_try
         V = V_try
         xi = xi_try; chi = chi_try
      end do

      if (.not. converged) then
         ierr = wdm_av_ierr_no_converge
         return
      end if
      Pstar = exp(L)

   end subroutine av_solve_real


   ! for a given L, invert each species to ln rho_i and return V = sum X_i/rho_i
   subroutine av_volume_at_L(L, tau, nsp, tbl_idx, A, X, xi, chi, V, ierr)
      real(dp), intent(in) :: L, tau
      integer, intent(in) :: nsp, tbl_idx(nsp)
      real(dp), intent(in) :: A(nsp), X(nsp)
      real(dp), intent(inout) :: xi(nsp)
      real(dp), intent(out) :: chi(nsp), V
      integer, intent(out) :: ierr
      integer :: i
      ierr = 0
      V = 0d0
      do i = 1, nsp
         call invert_species(wdm_av_tables(tbl_idx(i)), A(i), tau, L, xi(i), chi(i), ierr)
         if (ierr /= 0) return
         V = V + X(i)*exp(-xi(i))
      end do
   end subroutine av_volume_at_L


   ! 1-D Newton inversion of ln P_i(xi,tau) = Ltar for xi (warm started).
   subroutine invert_species(tbl, A_sp, tau, Ltar, xi, chi, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_sp, tau, Ltar
      real(dp), intent(inout) :: xi
      real(dp), intent(out) :: chi
      integer, intent(out) :: ierr
      type(eos_wdm_av_derivs) :: d
      real(dp) :: lnr, xi_lo, xi_hi, lnP, g, dxi, fx, xi_prev
      integer :: it
      integer, parameter :: maxit = 60
      real(dp), parameter :: step_cap = 4d0

      ierr = 0
      ! allowed xi window for this species (eval shifts by lnr = ln(A_tbl/A_sp))
      lnr = log(tbl% A/A_sp)
      xi_lo = tbl% xi_min - lnr
      xi_hi = tbl% xi_min + (tbl% nxi - 1)*tbl% dxi - lnr
      if (xi < xi_lo) xi = xi_lo
      if (xi > xi_hi) xi = xi_hi

      do it = 1, maxit
         call eos_wdm_av_eval_species_real(tbl, A_sp, xi, tau, d, ierr)
         if (ierr /= 0) return
         fx = d% fx
         ! NaN/Inf nodal data: the containing cell has invalid (Skye-failed)
         ! corners, so this species' state is outside the usable table.  Report
         ! out_of_table (D10) rather than burning the iteration budget and
         ! returning a spurious no_converge.
         if (is_bad(fx) .or. is_bad(d% fxx)) then
            ierr = wdm_av_ierr_out_of_table; return
         end if
         if (fx <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return      ! P_i <= 0
         end if
         chi = 1d0 + d% fxx/fx
         if (chi <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return      ! dP_i/drho_i <= 0 (non-monotone)
         end if
         lnP = xi + log(fx)                 ! ln(rho_i f_xi) = ln P_i
         g = lnP - Ltar
         ! inner tolerance kept above the biquintic evaluation round-off floor
         ! (~1e-12 at these log magnitudes) so the 1-D Newton terminates cleanly
         if (abs(g) < 1d-11) return
         dxi = -g/chi
         if (dxi > step_cap) dxi = step_cap
         if (dxi < -step_cap) dxi = -step_cap
         xi_prev = xi
         xi = xi + dxi
         if (xi < xi_lo) xi = xi_lo
         if (xi > xi_hi) xi = xi_hi
         ! pinned at a table edge with the Newton step pointing further out: the
         ! requested pressure L is not attainable within this species' tabulated
         ! density range (e.g. a partial density below/above the table, or the
         ! flat-pressure pair corner where chi->0).  Decline as out_of_table (D10)
         ! instead of iterating to maxit and returning no_converge.
         if (xi == xi_prev) then
            ierr = wdm_av_ierr_out_of_table; return
         end if
      end do
      ! not converged in the inner loop: report as non-convergence for this point
      ierr = wdm_av_ierr_no_converge
   end subroutine invert_species


   ! evaluate ln P_i and chi_rho,i at a given xi (no inversion)
   subroutine eval_lnP_chi(tbl, A_sp, xi, tau, lnP, chi, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_sp, xi, tau
      real(dp), intent(out) :: lnP, chi
      integer, intent(out) :: ierr
      type(eos_wdm_av_derivs) :: d
      ierr = 0
      call eos_wdm_av_eval_species_real(tbl, A_sp, xi, tau, d, ierr)
      if (ierr /= 0) return
      if (d% fx <= 0d0) then
         ierr = wdm_av_ierr_bad_domain; return
      end if
      chi = 1d0 + d% fxx/d% fx
      if (chi <= 0d0) then
         ierr = wdm_av_ierr_bad_domain; return
      end if
      lnP = xi + log(d% fx)
   end subroutine eval_lnP_chi


   !----------------------------------------------------------------------------
   ! (2) derivative pass: apply the flattened Newton update in auto_diff nsteps
   ! times, recomputing the Jacobian (chi_i, w_i) in auto_diff each step.
   !----------------------------------------------------------------------------

   subroutine av_newton_pass_ad(nsteps, tau_ad, v_ad, nsp, tbl_idx, A, X, L_ad, xi_ad, ierr)
      integer, intent(in) :: nsteps, nsp, tbl_idx(nsp)
      type(auto_diff_real_2var_order3), intent(in) :: tau_ad, v_ad
      real(dp), intent(in) :: A(nsp), X(nsp)
      type(auto_diff_real_2var_order3), intent(inout) :: L_ad, xi_ad(nsp)
      integer, intent(out) :: ierr

      type(auto_diff_real_2var_order3) :: f_ad, fx_ad, ft_ad, fxx_ad
      type(auto_diff_real_2var_order3) :: P_i, V, num, den_s, dL
      type(auto_diff_real_2var_order3) :: wq(nsp), rq(nsp), chiq(nsp)
      integer :: step, i

      ierr = 0
      do step = 1, nsteps
         V = 0d0
         den_s = 0d0
         do i = 1, nsp
            call eos_wdm_av_eval_species_ad(wdm_av_tables(tbl_idx(i)), A(i), xi_ad(i), tau_ad, &
                                            f_ad, fx_ad, ft_ad, fxx_ad, ierr)
            if (ierr /= 0) return
            P_i = exp(xi_ad(i))*fx_ad
            chiq(i) = 1d0 + fxx_ad/fx_ad
            rq(i) = L_ad - log(P_i)
            wq(i) = X(i)*exp(-xi_ad(i))         ! w_i = X_i/rho_i
            V = V + wq(i)
            den_s = den_s + wq(i)/chiq(i)
         end do
         ! delta_L = -[(v - V) + sum_i w_i r_i/chi_i] / sum_i w_i/chi_i
         num = (v_ad - V)
         do i = 1, nsp
            num = num + wq(i)*rq(i)/chiq(i)
         end do
         dL = -num/den_s
         L_ad = L_ad + dL
         do i = 1, nsp
            xi_ad(i) = xi_ad(i) + (rq(i) + dL)/chiq(i)
         end do
      end do
   end subroutine av_newton_pass_ad


   !----------------------------------------------------------------------------
   ! (3) ideal mixing entropy with actual AV partial volumes (Decision D8)
   !----------------------------------------------------------------------------

   subroutine mixing_entropy_ad(mode, nsp, A, Z, X, den, rho_ad, ds)
      integer, intent(in) :: mode, nsp
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp)
      type(auto_diff_real_2var_order3), intent(in) :: den, rho_ad(nsp)
      type(auto_diff_real_2var_order3), intent(out) :: ds

      type(auto_diff_real_2var_order3) :: ds_ion, ds_e, ne_sum, ne_i, frac
      real(dp) :: kmu, ye
      integer :: i

      ds = 0d0
      if (mode == wdm_av_mix_none) return
      kmu = kerg/amu

      ! ion term: Delta s_ion = -(k/mu) sum_i (X_i/A_i) ln(X_i rho/rho_i)
      ds_ion = 0d0
      do i = 1, nsp
         if (X(i) <= 0d0) cycle
         frac = X(i)*den/rho_ad(i)
         ds_ion = ds_ion - (X(i)/A(i))*log(frac)
      end do
      ds_ion = kmu*ds_ion
      ds = ds_ion

      if (mode == wdm_av_mix_ions_electrons) then
         ! electron term: Delta s_e = (k/mu) sum_i (X_i Z_i/A_i) ln(n_e,i/n_e)
         ! n_e,i = Z_i rho_i/(A_i mu),  n_e = (rho/mu) sum_i X_i Z_i/A_i
         ye = sum(X(1:nsp)*Z(1:nsp)/A(1:nsp))
         if (ye > 0d0) then
            ne_sum = (den/amu)*ye
            ds_e = 0d0
            do i = 1, nsp
               if (X(i) <= 0d0 .or. Z(i) <= 0d0) cycle
               ne_i = Z(i)*rho_ad(i)/(A(i)*amu)
               ds_e = ds_e + (X(i)*Z(i)/A(i))*log(ne_i/ne_sum)
            end do
            ds = ds + kmu*ds_e
         end if
      end if
   end subroutine mixing_entropy_ad


   !----------------------------------------------------------------------------
   ! (4) electron degeneracy parameter eta.  Production path uses HELM
   ! (compute_ideal_ele, as Skye does); if the HELM table is not loaded, fall
   ! back to the nondegenerate ideal-gas expression eta = ln(n_e/(2 n_Q,e)).
   !----------------------------------------------------------------------------

   subroutine electron_eta(temp, den, ytot1, ye, zbar, xnefer, etaele)
      type(auto_diff_real_2var_order3), intent(in) :: temp, den, xnefer
      real(dp), intent(in) :: ytot1, ye, zbar
      type(auto_diff_real_2var_order3), intent(out) :: etaele

      type(auto_diff_real_2var_order3) :: F_ele, xne_tmp, nQ, base
      real(dp) :: din_val, logtemp_val, logden_val
      integer :: ierr

      if (associated(eos_ht)) then
         din_val = ye*den% val
         logtemp_val = log10(temp% val)
         logden_val = log10(den% val)
         call compute_ideal_ele(temp% val, den% val, din_val, logtemp_val, logden_val, &
                                zbar, ytot1, ye, eos_ht, F_ele, etaele, xne_tmp, ierr)
         if (ierr == 0) return
      end if

      ! nondegenerate fallback: n_Q,e = (2 pi m_e k T / h^2)^{3/2}
      base = (2d0*pi*me*kerg/pow2(planck_h))*temp
      nQ = pow(base, 1.5d0)
      etaele = log(xnefer/(2d0*nQ))
   end subroutine electron_eta


   !----------------------------------------------------------------------------
   ! (6) fill the converged state (real, hand-derived Gibbs-space cross-check
   ! values for F_AV without the mixing term)
   !----------------------------------------------------------------------------

   subroutine fill_state(nsp, tbl_idx, A, X, xi, tau, T, rho, Pstar, L, state, ierr)
      integer, intent(in) :: nsp, tbl_idx(nsp)
      real(dp), intent(in) :: A(nsp), X(nsp), xi(nsp), tau, T, rho, Pstar, L
      type(eos_wdm_av_state), intent(inout) :: state
      integer, intent(out) :: ierr

      type(eos_wdm_av_derivs) :: d
      real(dp) :: rho_i, P_rho, P_T, cV, cP, gPP, gPT, gTT
      integer :: i

      ierr = 0
      state% P = Pstar
      state% lnP = L
      state% nsp = nsp
      if (allocated(state% rho_i)) deallocate(state% rho_i, state% xi_i, state% f_i, &
         state% g_i, state% h_i, state% e_i, state% s_i, state% v_i, &
         state% P_rho_i, state% P_T_i, state% chi_rho_i)
      allocate(state% rho_i(nsp), state% xi_i(nsp), state% f_i(nsp), state% g_i(nsp), &
               state% h_i(nsp), state% e_i(nsp), state% s_i(nsp), state% v_i(nsp), &
               state% P_rho_i(nsp), state% P_T_i(nsp), state% chi_rho_i(nsp))

      state% G_PP = 0d0; state% G_PT = 0d0; state% G_TT = 0d0
      do i = 1, nsp
         call eos_wdm_av_eval_species_real(wdm_av_tables(tbl_idx(i)), A(i), xi(i), tau, d, ierr)
         if (ierr /= 0) return
         rho_i = exp(xi(i))
         if (d% fx <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return
         end if
         P_rho = d% fx + d% fxx                   ! dP/drho = f_xi + f_xixi
         if (P_rho <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return
         end if
         P_T = rho_i*d% fxt/T                      ! dP/dT at fixed rho
         cV = (d% ft - d% ftt)/T
         cP = cV + T*P_T*P_T/(rho_i*rho_i*P_rho)

         state% rho_i(i) = rho_i
         state% xi_i(i) = xi(i)
         state% f_i(i) = d% f
         state% e_i(i) = d% f - d% ft             ! e = f - f_tau
         state% s_i(i) = -d% ft/T
         state% v_i(i) = 1d0/rho_i
         state% P_rho_i(i) = P_rho
         state% P_T_i(i) = P_T
         state% chi_rho_i(i) = 1d0 + d% fxx/d% fx
         state% g_i(i) = d% f + Pstar/rho_i
         state% h_i(i) = state% e_i(i) + Pstar/rho_i

         gPP = -1d0/(rho_i*rho_i*P_rho)
         gPT = P_T/(rho_i*rho_i*P_rho)
         gTT = -cP/T
         state% G_PP = state% G_PP + X(i)*gPP
         state% G_PT = state% G_PT + X(i)*gPT
         state% G_TT = state% G_TT + X(i)*gTT
      end do

      ! mixture Gibbs-space back-conversion (F_AV without the mixing term)
      state% P_rho_mix = -1d0/(rho*rho*state% G_PP)
      state% P_T_mix = -state% G_PT/state% G_PP
      state% cv_mix = -T*(state% G_TT - state% G_PT*state% G_PT/state% G_PP)
      state% cp_mix = -T*state% G_TT

   end subroutine fill_state


   !----------------------------------------------------------------------------
   ! common-electron-density rule (Decision D27): place each species at the
   ! density where its electron density equals the common mean n_e.  For full
   ! ionization rho_i = n_e A_i m_u / Z_i = rho * y_e * A_i / Z_i, with
   ! y_e = sum_i X_i Z_i/A_i.  Returns xi(i) = ln rho_i and lnc(i) = ln(rho_i/rho)
   ! (the constant multiplier, for the auto_diff seed).  No solve.
   !----------------------------------------------------------------------------

   subroutine ne_densities_real(nsp, A, Z, X, rho, xi, lnc, ierr)
      integer, intent(in) :: nsp
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp), rho
      real(dp), intent(out) :: xi(nsp), lnc(nsp)
      integer, intent(out) :: ierr
      real(dp) :: ye, logrho
      integer :: i

      ierr = 0
      ye = sum(X(1:nsp)*Z(1:nsp)/A(1:nsp))
      if (ye <= 0d0) then
         ierr = wdm_av_ierr_bad_domain       ! no electrons: no common n_e exists
         return
      end if
      logrho = log(rho)
      do i = 1, nsp
         if (Z(i) <= 0d0) then
            ! a neutral species cannot match a finite n_e at any density
            ! (the opposite weakness of additive volume; D27 writeup)
            ierr = wdm_av_ierr_bad_domain
            return
         end if
         lnc(i) = log(ye*A(i)/Z(i))          ! ln(rho_i/rho) = ln(y_e A_i/Z_i)
         xi(i) = logrho + lnc(i)
      end do
   end subroutine ne_densities_real


   !----------------------------------------------------------------------------
   ! fill the converged state for the common-n_e rule.  Unlike fill_state, the
   ! species pressures are NOT equal, so there is no single P* and the Gibbs-space
   ! cross-check (G_PP,...) does not apply; those are left zero.  The mixture
   ! pressure is the volume-weighted sum P = sum_i (V_i/V) P_i, which (because the
   ! rho_i are fixed multiples of rho) is exactly rho^2 dF_AV/drho -- the same P
   ! the auto_diff pass returns.
   !----------------------------------------------------------------------------

   subroutine fill_state_ne(nsp, tbl_idx, A, X, xi, tau, T, rho, state, ierr)
      integer, intent(in) :: nsp, tbl_idx(nsp)
      real(dp), intent(in) :: A(nsp), X(nsp), xi(nsp), tau, T, rho
      type(eos_wdm_av_state), intent(inout) :: state
      integer, intent(out) :: ierr

      type(eos_wdm_av_derivs) :: d
      real(dp) :: rho_i, P_i, P_rho, P_T, cV, Pmix, wvol
      integer :: i

      ierr = 0
      state% nsp = nsp
      if (allocated(state% rho_i)) deallocate(state% rho_i, state% xi_i, state% f_i, &
         state% g_i, state% h_i, state% e_i, state% s_i, state% v_i, &
         state% P_rho_i, state% P_T_i, state% chi_rho_i)
      allocate(state% rho_i(nsp), state% xi_i(nsp), state% f_i(nsp), state% g_i(nsp), &
               state% h_i(nsp), state% e_i(nsp), state% s_i(nsp), state% v_i(nsp), &
               state% P_rho_i(nsp), state% P_T_i(nsp), state% chi_rho_i(nsp))

      ! cross-check quantities are not defined for this rule
      state% G_PP = 0d0; state% G_PT = 0d0; state% G_TT = 0d0
      state% P_rho_mix = 0d0; state% P_T_mix = 0d0
      state% cv_mix = 0d0; state% cp_mix = 0d0

      Pmix = 0d0
      do i = 1, nsp
         call eos_wdm_av_eval_species_real(wdm_av_tables(tbl_idx(i)), A(i), xi(i), tau, d, ierr)
         if (ierr /= 0) return
         rho_i = exp(xi(i))
         if (d% fx <= 0d0) then
            ierr = wdm_av_ierr_bad_domain; return      ! P_i <= 0
         end if
         P_i = rho_i*d% fx                         ! P_i = rho_i f_xi
         P_rho = d% fx + d% fxx                     ! dP_i/drho_i
         P_T = rho_i*d% fxt/T                        ! dP_i/dT at fixed rho_i
         cV = (d% ft - d% ftt)/T

         state% rho_i(i) = rho_i
         state% xi_i(i) = xi(i)
         state% f_i(i) = d% f
         state% e_i(i) = d% f - d% ft               ! e = f - f_tau
         state% s_i(i) = -d% ft/T
         state% v_i(i) = 1d0/rho_i
         state% P_rho_i(i) = P_rho
         state% P_T_i(i) = P_T
         state% chi_rho_i(i) = 1d0 + d% fxx/d% fx
         state% g_i(i) = d% f + P_i/rho_i           ! each species at its own P_i
         state% h_i(i) = state% e_i(i) + P_i/rho_i

         wvol = X(i)*rho/rho_i                       ! V_i/V (sum_i wvol = 1)
         Pmix = Pmix + wvol*P_i
      end do

      state% P = Pmix
      state% lnP = log(Pmix)

   end subroutine fill_state_ne


   !----------------------------------------------------------------------------
   ! debug cross-check: hand-derived Gibbs values vs auto_diff for F_AV alone
   !----------------------------------------------------------------------------

   subroutine crosscheck_print(F_AV, temp, den, state)
      type(auto_diff_real_2var_order3), intent(in) :: F_AV, temp, den
      type(eos_wdm_av_state), intent(in) :: state
      type(auto_diff_real_2var_order3) :: s, e, p, cv, cp, chit, chid, gam1, gam2, gam3, nabad, cs
      real(dp) :: P_rho_ad, P_T_ad, cv_ad, cp_ad
      real(dp) :: r_prho, r_pt, r_cv, r_cp

      call thermodynamics_from_free_energy(F_AV, temp, den, s, e, p)
      call compute_derived_quantities(temp, den, s, e, p, cv, cp, chit, chid, &
                                      gam1, gam2, gam3, nabad, cs)
      P_rho_ad = chid% val * p% val / den% val
      P_T_ad = chit% val * p% val / temp% val
      cv_ad = cv% val
      cp_ad = cp% val

      r_prho = reldiff(state% P_rho_mix, P_rho_ad)
      r_pt   = reldiff(state% P_T_mix, P_T_ad)
      r_cv   = reldiff(state% cv_mix, cv_ad)
      r_cp   = reldiff(state% cp_mix, cp_ad)
      write(*,'(a,4es12.3)') ' WDM_AV crosscheck relerr [P_rho,P_T,cV,cP] = ', &
         r_prho, r_pt, r_cv, r_cp
   end subroutine crosscheck_print

   real(dp) function reldiff(a, b) result(r)
      real(dp), intent(in) :: a, b
      r = abs(a - b)/max(abs(b), 1d-99)
   end function reldiff


   !----------------------------------------------------------------------------
   ! validity-mask helper (replicates find_cell indexing; eos_wdm_av_interp does
   ! not expose a validity check).  Returns .true. when the table has no mask or
   ! all four corners of the containing cell are valid.
   !----------------------------------------------------------------------------

   logical function cell_usable(tbl, A_sp, xi, tau) result(ok)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_sp, xi, tau
      real(dp) :: xi_eval, lnr
      integer :: ixi, itau
      ok = .true.
      if (.not. allocated(tbl% valid)) return
      lnr = log(tbl% A/A_sp)
      xi_eval = xi + lnr
      ixi = floor((xi_eval - tbl% xi_min)/tbl% dxi) + 1
      ixi = max(1, min(ixi, tbl% nxi - 1))
      itau = floor((tau - tbl% tau_min)/tbl% dtau) + 1
      itau = max(1, min(itau, tbl% ntau - 1))
      if (tbl% valid(ixi,   itau)   == 0 .or. tbl% valid(ixi+1, itau)   == 0 .or. &
          tbl% valid(ixi,   itau+1) == 0 .or. tbl% valid(ixi+1, itau+1) == 0) ok = .false.
   end function cell_usable

end module eos_wdm_av_mix
