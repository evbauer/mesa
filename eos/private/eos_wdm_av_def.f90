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

module eos_wdm_av_def
   ! Data structures for the WDM_AV (warm dense matter, additive volume)
   ! per-species tabulated equation of state.  The processed-table layout this
   ! module mirrors is documented in eos/eosWDM_AV_builder/FORMAT.md.

   use const_def, only: dp, i8

   implicit none
   public

   ! ------------------------------------------------------------------------
   ! TEST HOOK, scratch branch only, never merge into WDM_EOS.
   ! Module-level flag enabling the "WDM_AV replaces Skye" drop-in experiment
   ! (brief WP8 / Decision D31).  When true, the Skye cascade level evaluates
   ! WDM_AV on the Skye-derived per-isotope tables in Skye's place, blended into
   ! Skye through the Gamma window (see eos_wdm_av_eval::get_WDM_AV_as_Skye_for_eosdt).
   ! It is stored as a module variable, set from the &eos namelist in
   ! eos_ctrls_io, rather than a field of EoS_General_Info in eos_def.f90, so that
   ! this test hook requires no full ./install (a type change in eos_def forces a
   ! rebuild of star/kap/net/star_data).  Consequence: the flag is GLOBAL, not
   ! per-handle; this is acceptable for a single-handle test and never ships.
   logical :: WDM_AV_replaces_Skye_flag = .false.

   ! Drop-in call/decline counters (brief WP8/T1).  Updated under the named OMP
   ! critical section wdm_av_dropin_counters.  A "decline" is a point where the
   ! Gamma window asks WDM_AV to serve (w > 0) but WDM_AV returns skip/ierr, so
   ! Skye serves it alone; w == 0 (crystallizing core handed to Skye by design)
   ! is NOT a decline.
   integer(i8) :: wdm_av_dropin_calls = 0_i8
   integer(i8) :: wdm_av_dropin_declines = 0_i8
   ! ------------------------------------------------------------------------

   ! format version understood by this reader (see FORMAT.md)
   integer, parameter :: wdm_av_format_version = 1

   ! maximum number of per-isotope tables held in the module-level store
   integer, parameter :: max_wdm_av_tables = 64

   ! interpolation orders
   integer, parameter :: wdm_av_order_bicubic = 3
   integer, parameter :: wdm_av_order_biquintic = 5

   ! error codes (negative so they never collide with a normal ierr==0)
   integer, parameter :: wdm_av_ierr_out_of_table = -101
   integer, parameter :: wdm_av_ierr_bad_format = -102
   integer, parameter :: wdm_av_ierr_not_loaded = -103

   ! one processed table for a single isotope.  All nodal data are the Helmholtz
   ! free energy per gram (erg/g) and its derivatives with respect to the log
   ! coordinates xi = ln(rho) and tau = ln(T), on a uniform (nxi, ntau) grid.
   type eos_wdm_av_table
      character(len=8) :: isotope = ''
      integer :: Z = 0
      real(dp) :: A = 0d0
      integer :: order = wdm_av_order_biquintic   ! 3 (bicubic) or 5 (biquintic)
      integer :: nxi = 0, ntau = 0
      real(dp) :: xi_min = 0d0, dxi = 0d0
      real(dp) :: tau_min = 0d0, dtau = 0d0
      real(dp) :: e_offset = 0d0, f_offset = 0d0
      logical :: loaded = .false.
      logical :: has_Zbar = .false.
      logical :: has_eta = .false.
      ! nodal free energy and derivatives, all dimensioned (nxi, ntau):
      !   f      = f
      !   fx     = df/dxi
      !   ft     = df/dtau
      !   fxx    = d2f/dxi2
      !   ftt    = d2f/dtau2
      !   fxt    = d2f/dxi dtau
      !   fxxt   = d3f/dxi2 dtau
      !   fxtt   = d3f/dxi dtau2
      !   fxxtt  = d4f/dxi2 dtau2
      ! (for order 3 only f, fx, ft, fxt are used in interpolation)
      real(dp), allocatable :: f(:,:), fx(:,:), ft(:,:), fxx(:,:), ftt(:,:), fxt(:,:), &
                               fxxt(:,:), fxtt(:,:), fxxtt(:,:)
      ! optional bicubic nodal data for the mean ionization (allocated only when present)
      real(dp), allocatable :: Zbar(:,:), Zbar_x(:,:), Zbar_t(:,:), Zbar_xt(:,:)
      ! optional bicubic nodal data for the degeneracy parameter
      real(dp), allocatable :: eta(:,:), eta_x(:,:), eta_t(:,:), eta_xt(:,:)
      ! optional validity mask (1 = usable node, 0 = unusable)
      integer(1), allocatable :: valid(:,:)
   end type eos_wdm_av_table

   ! free energy and its derivatives through third order at one (xi, tau) point,
   ! with respect to xi = ln(rho) and tau = ln(T)
   type eos_wdm_av_derivs
      real(dp) :: f = 0d0
      real(dp) :: fx = 0d0, ft = 0d0
      real(dp) :: fxx = 0d0, fxt = 0d0, ftt = 0d0
      real(dp) :: fxxx = 0d0, fxxt = 0d0, fxtt = 0d0, fttt = 0d0
   end type eos_wdm_av_derivs

   ! module-level table store (read-only after lazy load under OMP critical)
   type(eos_wdm_av_table), allocatable, target :: wdm_av_tables(:)
   integer :: num_wdm_av_tables = 0

contains

   ! return the store index of a loaded table for the given isotope, or 0 if absent
   integer function eos_wdm_av_table_index(isotope_name) result(idx)
      character(len=*), intent(in) :: isotope_name
      integer :: i
      idx = 0
      if (.not. allocated(wdm_av_tables)) return
      do i = 1, num_wdm_av_tables
         if (wdm_av_tables(i)% loaded) then
            if (trim(wdm_av_tables(i)% isotope) == trim(isotope_name)) then
               idx = i
               return
            end if
         end if
      end do
   end function eos_wdm_av_table_index

end module eos_wdm_av_def
