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

module eos_wdm_av_interp
   ! Per-species Hermite interpolation of the WDM_AV free-energy tables in the
   ! log coordinates xi = ln(rho), tau = ln(T).  Both a bicubic (order 3) and a
   ! biquintic (order 5) local Hermite scheme are supported, reusing the HELM
   ! basis polynomials and tensor-product sums in helm_polynomials.  The weight
   ! construction (including the sign flips on the mirrored first-derivative
   ! terms) follows helm_electron_positron.dek, specialized to a uniform grid.
   !
   ! Correspondence with the HELM machinery (there d = density, t = temperature;
   ! here d = xi, t = tau):
   !   fi( 1: 4) f       fi( 5: 8) f_tau    fi( 9:12) f_tautau
   !   fi(13:16) f_xi    fi(17:20) f_xixi   fi(21:24) f_xitau
   !   fi(25:28) f_xixitau   fi(29:32) f_xitautau   fi(33:36) f_xixitautau
   ! corner order within each group: (ixi,itau) (ixi+1,itau) (ixi,itau+1) (ixi+1,itau+1).

   use const_def, only: dp, kerg, amu
   use eos_wdm_av_def
   use helm_polynomials
   use auto_diff
   use math_lib

   implicit none
   private

   public :: eos_wdm_av_eval_real
   public :: eos_wdm_av_eval_ad
   public :: eos_wdm_av_eval_species_real
   public :: eos_wdm_av_eval_species_ad

   ! out-of-table tolerance, as a fraction of the grid spacing
   real(dp), parameter :: edge_tol = 1d-6
   ! relative tolerance for treating A_species as equal to the tabulated A
   real(dp), parameter :: same_A_rtol = 1d-14

contains

   !----------------------------------------------------------------------------
   ! real evaluator
   !----------------------------------------------------------------------------

   subroutine eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: xi, tau
      type(eos_wdm_av_derivs), intent(out) :: d
      integer, intent(out) :: ierr

      integer :: ixi, itau
      real(dp) :: tx, ut
      real(dp) :: fi(36)
      real(dp) :: vx(6), d1x(6), d2x(6), d3x(6)
      real(dp) :: vt(6), d1t(6), d2t(6), d3t(6)
      real(dp) :: cvx(4), cd1x(4), cd2x(4), cd3x(4)
      real(dp) :: cvt(4), cd1t(4), cd2t(4), cd3t(4)

      if (.not. tbl% loaded) then
         ierr = wdm_av_ierr_not_loaded
         return
      end if

      call find_cell(tbl, xi, tau, ixi, itau, tx, ut, ierr)
      if (ierr /= 0) return

      if (tbl% order == wdm_av_order_biquintic) then
         call qwts(tx, tbl% dxi, vx, d1x, d2x, d3x)
         call qwts(ut, tbl% dtau, vt, d1t, d2t, d3t)
         call fill_fi5(tbl, ixi, itau, fi)
         d% f    = hq(fi, vt,  vx)
         d% fx   = hq(fi, vt,  d1x)
         d% ft   = hq(fi, d1t, vx)
         d% fxx  = hq(fi, vt,  d2x)
         d% fxt  = hq(fi, d1t, d1x)
         d% ftt  = hq(fi, d2t, vx)
         d% fxxx = hq(fi, vt,  d3x)
         d% fxxt = hq(fi, d1t, d2x)
         d% fxtt = hq(fi, d2t, d1x)
         d% fttt = hq(fi, d3t, vx)
      else
         call cwts(tx, tbl% dxi, cvx, cd1x, cd2x, cd3x)
         call cwts(ut, tbl% dtau, cvt, cd1t, cd2t, cd3t)
         call fill_fi3(tbl, ixi, itau, fi)
         d% f    = hc(fi, cvt,  cvx)
         d% fx   = hc(fi, cvt,  cd1x)
         d% ft   = hc(fi, cd1t, cvx)
         d% fxx  = hc(fi, cvt,  cd2x)
         d% fxt  = hc(fi, cd1t, cd1x)
         d% ftt  = hc(fi, cd2t, cvx)
         d% fxxx = hc(fi, cvt,  cd3x)
         d% fxxt = hc(fi, cd1t, cd2x)
         d% fxtt = hc(fi, cd2t, cd1x)
         d% fttt = hc(fi, cd3t, cvx)
      end if

      ierr = 0
   end subroutine eos_wdm_av_eval_real


   !----------------------------------------------------------------------------
   ! auto_diff evaluator: returns f, fx = df/dxi, ft = df/dtau, fxx = d2f/dxi2,
   ! each carrying correct derivatives through third order with respect to
   ! whatever xi and tau depend on.
   !----------------------------------------------------------------------------

   subroutine eos_wdm_av_eval_ad(tbl, xi, tau, f, fx, ft, fxx, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      type(auto_diff_real_2var_order3), intent(in) :: xi, tau
      type(auto_diff_real_2var_order3), intent(out) :: f, fx, ft, fxx
      integer, intent(out) :: ierr

      integer :: ixi, itau
      real(dp) :: tx0, ut0, xi0, tau0
      real(dp) :: fi(36)
      type(auto_diff_real_2var_order3) :: tx, ut
      type(auto_diff_real_2var_order3) :: vx(6), d1x(6), d2x(6)
      type(auto_diff_real_2var_order3) :: vt(6), d1t(6), d2t(6)
      type(auto_diff_real_2var_order3) :: cvx(4), cd1x(4), cd2x(4)
      type(auto_diff_real_2var_order3) :: cvt(4), cd1t(4), cd2t(4)

      if (.not. tbl% loaded) then
         ierr = wdm_av_ierr_not_loaded
         return
      end if

      call find_cell(tbl, xi% val, tau% val, ixi, itau, tx0, ut0, ierr)
      if (ierr /= 0) return

      xi0 = tbl% xi_min + (ixi-1)*tbl% dxi
      tau0 = tbl% tau_min + (itau-1)*tbl% dtau
      tx = (xi - xi0)/tbl% dxi
      ut = (tau - tau0)/tbl% dtau

      if (tbl% order == wdm_av_order_biquintic) then
         call qwts_ad(tx, tbl% dxi, vx, d1x, d2x)
         call qwts_ad(ut, tbl% dtau, vt, d1t, d2t)
         call fill_fi5(tbl, ixi, itau, fi)
         f   = hq_ad(fi, vt,  vx)
         fx  = hq_ad(fi, vt,  d1x)
         ft  = hq_ad(fi, d1t, vx)
         fxx = hq_ad(fi, vt,  d2x)
      else
         call cwts_ad(tx, tbl% dxi, cvx, cd1x, cd2x)
         call cwts_ad(ut, tbl% dtau, cvt, cd1t, cd2t)
         call fill_fi3(tbl, ixi, itau, fi)
         f   = hc_ad(fi, cvt,  cvx)
         fx  = hc_ad(fi, cvt,  cd1x)
         ft  = hc_ad(fi, cd1t, cvx)
         fxx = hc_ad(fi, cvt,  cd2x)
      end if

      ierr = 0
   end subroutine eos_wdm_av_eval_ad


   !----------------------------------------------------------------------------
   ! isotope rescaling (Decision D6): evaluate the table for a species of mass
   ! A_species using the number-density-preserving rescaling
   !   f_A'(rho,T) = (A/A') f_A(rho A/A', T) - (3/2) k T/(A' m_u) ln(A'/A)
   ! In log coordinates: shift xi by ln(r), r = A_tbl/A_species, scale all
   ! derivatives by r, and add the mass term -c T (linear in T) to f and to the
   ! pure-tau derivatives ft, ftt, fttt.
   !----------------------------------------------------------------------------

   subroutine eos_wdm_av_eval_species_real(tbl, A_species, xi, tau, d, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_species, xi, tau
      type(eos_wdm_av_derivs), intent(out) :: d
      integer, intent(out) :: ierr

      type(eos_wdm_av_derivs) :: dt
      real(dp) :: r, lnr, c, cT

      if (abs(A_species - tbl% A) <= same_A_rtol*tbl% A) then
         call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
         return
      end if

      r = tbl% A / A_species
      lnr = log(r)
      c = 1.5d0 * kerg/(A_species*amu) * log(A_species/tbl% A)

      call eos_wdm_av_eval_real(tbl, xi + lnr, tau, dt, ierr)
      if (ierr /= 0) return

      cT = c * exp(tau)

      d% f    = r*dt% f    - cT
      d% fx   = r*dt% fx
      d% ft   = r*dt% ft   - cT
      d% fxx  = r*dt% fxx
      d% fxt  = r*dt% fxt
      d% ftt  = r*dt% ftt  - cT
      d% fxxx = r*dt% fxxx
      d% fxxt = r*dt% fxxt
      d% fxtt = r*dt% fxtt
      d% fttt = r*dt% fttt - cT

      ierr = 0
   end subroutine eos_wdm_av_eval_species_real


   subroutine eos_wdm_av_eval_species_ad(tbl, A_species, xi, tau, f, fx, ft, fxx, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_species
      type(auto_diff_real_2var_order3), intent(in) :: xi, tau
      type(auto_diff_real_2var_order3), intent(out) :: f, fx, ft, fxx
      integer, intent(out) :: ierr

      real(dp) :: r, lnr, c
      type(auto_diff_real_2var_order3) :: cT

      if (abs(A_species - tbl% A) <= same_A_rtol*tbl% A) then
         call eos_wdm_av_eval_ad(tbl, xi, tau, f, fx, ft, fxx, ierr)
         return
      end if

      r = tbl% A / A_species
      lnr = log(r)
      c = 1.5d0 * kerg/(A_species*amu) * log(A_species/tbl% A)

      call eos_wdm_av_eval_ad(tbl, xi + lnr, tau, f, fx, ft, fxx, ierr)
      if (ierr /= 0) return

      cT = c * exp(tau)

      f   = r*f  - cT
      fx  = r*fx
      ft  = r*ft - cT
      fxx = r*fxx

      ierr = 0
   end subroutine eos_wdm_av_eval_species_ad


   !----------------------------------------------------------------------------
   ! cell lookup
   !----------------------------------------------------------------------------

   subroutine find_cell(tbl, xi, tau, ixi, itau, tx, ut, ierr)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: xi, tau
      integer, intent(out) :: ixi, itau
      real(dp), intent(out) :: tx, ut
      integer, intent(out) :: ierr
      real(dp) :: xi_max, tau_max

      ierr = 0
      xi_max = tbl% xi_min + (tbl% nxi-1)*tbl% dxi
      tau_max = tbl% tau_min + (tbl% ntau-1)*tbl% dtau

      if (xi  < tbl% xi_min  - edge_tol*tbl% dxi  .or. xi  > xi_max  + edge_tol*tbl% dxi  .or. &
          tau < tbl% tau_min - edge_tol*tbl% dtau .or. tau > tau_max + edge_tol*tbl% dtau) then
         ierr = wdm_av_ierr_out_of_table
         return
      end if

      ixi = floor((xi - tbl% xi_min)/tbl% dxi) + 1
      ixi = max(1, min(ixi, tbl% nxi-1))
      itau = floor((tau - tbl% tau_min)/tbl% dtau) + 1
      itau = max(1, min(itau, tbl% ntau-1))

      tx = (xi  - (tbl% xi_min  + (ixi -1)*tbl% dxi)) / tbl% dxi
      ut = (tau - (tbl% tau_min + (itau-1)*tbl% dtau)) / tbl% dtau
   end subroutine find_cell


   !----------------------------------------------------------------------------
   ! fill the fi(36) nodal array for biquintic / fi(1:16) for bicubic
   !----------------------------------------------------------------------------

   subroutine set_corners(arr, ixi, itau, fi, base)
      real(dp), intent(in) :: arr(:,:)
      integer, intent(in) :: ixi, itau, base
      real(dp), intent(inout) :: fi(36)
      fi(base+1) = arr(ixi,   itau)
      fi(base+2) = arr(ixi+1, itau)
      fi(base+3) = arr(ixi,   itau+1)
      fi(base+4) = arr(ixi+1, itau+1)
   end subroutine set_corners

   subroutine fill_fi5(tbl, ixi, itau, fi)
      type(eos_wdm_av_table), intent(in) :: tbl
      integer, intent(in) :: ixi, itau
      real(dp), intent(out) :: fi(36)
      call set_corners(tbl% f,     ixi, itau, fi,  0)
      call set_corners(tbl% ft,    ixi, itau, fi,  4)
      call set_corners(tbl% ftt,   ixi, itau, fi,  8)
      call set_corners(tbl% fx,    ixi, itau, fi, 12)
      call set_corners(tbl% fxx,   ixi, itau, fi, 16)
      call set_corners(tbl% fxt,   ixi, itau, fi, 20)
      call set_corners(tbl% fxxt,  ixi, itau, fi, 24)
      call set_corners(tbl% fxtt,  ixi, itau, fi, 28)
      call set_corners(tbl% fxxtt, ixi, itau, fi, 32)
   end subroutine fill_fi5

   subroutine fill_fi3(tbl, ixi, itau, fi)
      type(eos_wdm_av_table), intent(in) :: tbl
      integer, intent(in) :: ixi, itau
      real(dp), intent(out) :: fi(36)
      fi = 0d0
      call set_corners(tbl% f,   ixi, itau, fi,  0)
      call set_corners(tbl% ft,  ixi, itau, fi,  4)
      call set_corners(tbl% fx,  ixi, itau, fi,  8)
      call set_corners(tbl% fxt, ixi, itau, fi, 12)
   end subroutine fill_fi3


   !----------------------------------------------------------------------------
   ! real quintic / cubic weight 6- / 4-vectors on a uniform grid of spacing h.
   ! Order: (w0, w1, w2, w0m, w1m, w2m) for quintic; (w0, w1, w0m, w1m) cubic.
   !----------------------------------------------------------------------------

   subroutine qwts(t, h, v, d1, d2, d3)
      real(dp), intent(in) :: t, h
      real(dp), intent(out) :: v(6), d1(6), d2(6), d3(6)
      real(dp) :: mt, hi1, hi2, hi3, h2
      mt = 1d0 - t
      hi1 = 1d0/h
      hi2 = hi1*hi1
      hi3 = hi2*hi1
      h2 = h*h
      ! value weights
      v(1) = psi0(t);        v(2) = psi1(t)*h;        v(3) = psi2(t)*h2
      v(4) = psi0(mt);       v(5) = -psi1(mt)*h;      v(6) = psi2(mt)*h2
      ! first derivative wrt the physical coordinate
      d1(1) = dpsi0(t)*hi1;  d1(2) = dpsi1(t);        d1(3) = dpsi2(t)*h
      d1(4) = -dpsi0(mt)*hi1; d1(5) = dpsi1(mt);      d1(6) = -dpsi2(mt)*h
      ! second derivative
      d2(1) = ddpsi0(t)*hi2; d2(2) = ddpsi1(t)*hi1;   d2(3) = ddpsi2(t)
      d2(4) = ddpsi0(mt)*hi2; d2(5) = -ddpsi1(mt)*hi1; d2(6) = ddpsi2(mt)
      ! third derivative
      d3(1) = dddpsi0(t)*hi3; d3(2) = dddpsi1(t)*hi2; d3(3) = dddpsi2(t)*hi1
      d3(4) = -dddpsi0(mt)*hi3; d3(5) = dddpsi1(mt)*hi2; d3(6) = -dddpsi2(mt)*hi1
   end subroutine qwts

   subroutine cwts(t, h, v, d1, d2, d3)
      real(dp), intent(in) :: t, h
      real(dp), intent(out) :: v(4), d1(4), d2(4), d3(4)
      real(dp) :: mt, hi1, hi2, hi3
      mt = 1d0 - t
      hi1 = 1d0/h
      hi2 = hi1*hi1
      hi3 = hi2*hi1
      v(1) = xpsi0(t);        v(2) = xpsi1(t)*h;      v(3) = xpsi0(mt);        v(4) = -xpsi1(mt)*h
      d1(1) = xdpsi0(t)*hi1;  d1(2) = xdpsi1(t);      d1(3) = -xdpsi0(mt)*hi1; d1(4) = xdpsi1(mt)
      d2(1) = xddpsi0(t)*hi2; d2(2) = xddpsi1(t)*hi1; d2(3) = xddpsi0(mt)*hi2; d2(4) = -xddpsi1(mt)*hi1
      d3(1) = xdddpsi0(t)*hi3; d3(2) = xdddpsi1(t)*hi2; d3(3) = -xdddpsi0(mt)*hi3; d3(4) = xdddpsi1(mt)*hi2
   end subroutine cwts


   !----------------------------------------------------------------------------
   ! real tensor-product sums (thin wrappers on helm_polynomials h5 / h3)
   !----------------------------------------------------------------------------

   real(dp) function hq(fi, wt, wd) result(r)
      real(dp), intent(in) :: fi(36), wt(6), wd(6)
      r = h5(1, 1, fi, &
             wt(1), wt(2), wt(3), wt(4), wt(5), wt(6), &
             wd(1), wd(2), wd(3), wd(4), wd(5), wd(6))
   end function hq

   real(dp) function hc(fi, wt, wd) result(r)
      real(dp), intent(in) :: fi(36), wt(4), wd(4)
      r = h3(1, 1, fi, &
             wt(1), wt(2), wt(3), wt(4), &
             wd(1), wd(2), wd(3), wd(4))
   end function hc


   !----------------------------------------------------------------------------
   ! auto_diff basis polynomials (expanded forms, using powN from math_lib)
   !----------------------------------------------------------------------------

   function ad_psi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -6d0*pow5(z) + 15d0*pow4(z) - 10d0*pow3(z) + 1d0
   end function ad_psi0

   function ad_psi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -3d0*pow5(z) + 8d0*pow4(z) - 6d0*pow3(z) + z
   end function ad_psi1

   function ad_psi2(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -0.5d0*pow5(z) + 1.5d0*pow4(z) - 1.5d0*pow3(z) + 0.5d0*pow2(z)
   end function ad_psi2

   function ad_dpsi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -30d0*pow4(z) + 60d0*pow3(z) - 30d0*pow2(z)
   end function ad_dpsi0

   function ad_dpsi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -15d0*pow4(z) + 32d0*pow3(z) - 18d0*pow2(z) + 1d0
   end function ad_dpsi1

   function ad_dpsi2(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -2.5d0*pow4(z) + 6d0*pow3(z) - 4.5d0*pow2(z) + z
   end function ad_dpsi2

   function ad_ddpsi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -120d0*pow3(z) + 180d0*pow2(z) - 60d0*z
   end function ad_ddpsi0

   function ad_ddpsi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -60d0*pow3(z) + 96d0*pow2(z) - 36d0*z
   end function ad_ddpsi1

   function ad_ddpsi2(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = -10d0*pow3(z) + 18d0*pow2(z) - 9d0*z + 1d0
   end function ad_ddpsi2

   function ad_xpsi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = 2d0*pow3(z) - 3d0*pow2(z) + 1d0
   end function ad_xpsi0

   function ad_xpsi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = pow3(z) - 2d0*pow2(z) + z
   end function ad_xpsi1

   function ad_xdpsi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = 6d0*pow2(z) - 6d0*z
   end function ad_xdpsi0

   function ad_xdpsi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = 3d0*pow2(z) - 4d0*z + 1d0
   end function ad_xdpsi1

   function ad_xddpsi0(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = 12d0*z - 6d0
   end function ad_xddpsi0

   function ad_xddpsi1(z) result(r)
      type(auto_diff_real_2var_order3), intent(in) :: z
      type(auto_diff_real_2var_order3) :: r
      r = 6d0*z - 4d0
   end function ad_xddpsi1


   !----------------------------------------------------------------------------
   ! auto_diff quintic / cubic weight vectors
   !----------------------------------------------------------------------------

   subroutine qwts_ad(t, h, v, d1, d2)
      type(auto_diff_real_2var_order3), intent(in) :: t
      real(dp), intent(in) :: h
      type(auto_diff_real_2var_order3), intent(out) :: v(6), d1(6), d2(6)
      type(auto_diff_real_2var_order3) :: mt
      real(dp) :: hi1, hi2, h2
      mt = 1d0 - t
      hi1 = 1d0/h
      hi2 = hi1*hi1
      h2 = h*h
      v(1) = ad_psi0(t);          v(2) = h*ad_psi1(t);        v(3) = h2*ad_psi2(t)
      v(4) = ad_psi0(mt);         v(5) = (-h)*ad_psi1(mt);    v(6) = h2*ad_psi2(mt)
      d1(1) = hi1*ad_dpsi0(t);    d1(2) = ad_dpsi1(t);        d1(3) = h*ad_dpsi2(t)
      d1(4) = (-hi1)*ad_dpsi0(mt); d1(5) = ad_dpsi1(mt);      d1(6) = (-h)*ad_dpsi2(mt)
      d2(1) = hi2*ad_ddpsi0(t);   d2(2) = hi1*ad_ddpsi1(t);   d2(3) = ad_ddpsi2(t)
      d2(4) = hi2*ad_ddpsi0(mt);  d2(5) = (-hi1)*ad_ddpsi1(mt); d2(6) = ad_ddpsi2(mt)
   end subroutine qwts_ad

   subroutine cwts_ad(t, h, v, d1, d2)
      type(auto_diff_real_2var_order3), intent(in) :: t
      real(dp), intent(in) :: h
      type(auto_diff_real_2var_order3), intent(out) :: v(4), d1(4), d2(4)
      type(auto_diff_real_2var_order3) :: mt
      real(dp) :: hi1, hi2
      mt = 1d0 - t
      hi1 = 1d0/h
      hi2 = hi1*hi1
      v(1) = ad_xpsi0(t);          v(2) = h*ad_xpsi1(t);       v(3) = ad_xpsi0(mt);          v(4) = (-h)*ad_xpsi1(mt)
      d1(1) = hi1*ad_xdpsi0(t);    d1(2) = ad_xdpsi1(t);       d1(3) = (-hi1)*ad_xdpsi0(mt); d1(4) = ad_xdpsi1(mt)
      d2(1) = hi2*ad_xddpsi0(t);   d2(2) = hi1*ad_xddpsi1(t);  d2(3) = hi2*ad_xddpsi0(mt);   d2(4) = (-hi1)*ad_xddpsi1(mt)
   end subroutine cwts_ad


   !----------------------------------------------------------------------------
   ! auto_diff tensor-product sums (mirror helm_polynomials h5 / h3)
   !----------------------------------------------------------------------------

   function hq_ad(fi, wt, wd) result(r)
      real(dp), intent(in) :: fi(36)
      type(auto_diff_real_2var_order3), intent(in) :: wt(6), wd(6)
      type(auto_diff_real_2var_order3) :: r
      ! wd = (w0d,w1d,w2d,w0md,w1md,w2md); wt = (w0t,w1t,w2t,w0mt,w1mt,w2mt)
      r = fi(1) *(wd(1)*wt(1)) + fi(2) *(wd(4)*wt(1)) &
        + fi(3) *(wd(1)*wt(4)) + fi(4) *(wd(4)*wt(4)) &
        + fi(5) *(wd(1)*wt(2)) + fi(6) *(wd(4)*wt(2)) &
        + fi(7) *(wd(1)*wt(5)) + fi(8) *(wd(4)*wt(5)) &
        + fi(9) *(wd(1)*wt(3)) + fi(10)*(wd(4)*wt(3)) &
        + fi(11)*(wd(1)*wt(6)) + fi(12)*(wd(4)*wt(6)) &
        + fi(13)*(wd(2)*wt(1)) + fi(14)*(wd(5)*wt(1)) &
        + fi(15)*(wd(2)*wt(4)) + fi(16)*(wd(5)*wt(4)) &
        + fi(17)*(wd(3)*wt(1)) + fi(18)*(wd(6)*wt(1)) &
        + fi(19)*(wd(3)*wt(4)) + fi(20)*(wd(6)*wt(4)) &
        + fi(21)*(wd(2)*wt(2)) + fi(22)*(wd(5)*wt(2)) &
        + fi(23)*(wd(2)*wt(5)) + fi(24)*(wd(5)*wt(5)) &
        + fi(25)*(wd(3)*wt(2)) + fi(26)*(wd(6)*wt(2)) &
        + fi(27)*(wd(3)*wt(5)) + fi(28)*(wd(6)*wt(5)) &
        + fi(29)*(wd(2)*wt(3)) + fi(30)*(wd(5)*wt(3)) &
        + fi(31)*(wd(2)*wt(6)) + fi(32)*(wd(5)*wt(6)) &
        + fi(33)*(wd(3)*wt(3)) + fi(34)*(wd(6)*wt(3)) &
        + fi(35)*(wd(3)*wt(6)) + fi(36)*(wd(6)*wt(6))
   end function hq_ad

   function hc_ad(fi, wt, wd) result(r)
      real(dp), intent(in) :: fi(36)
      type(auto_diff_real_2var_order3), intent(in) :: wt(4), wd(4)
      type(auto_diff_real_2var_order3) :: r
      ! wd = (w0d,w1d,w0md,w1md); wt = (w0t,w1t,w0mt,w1mt)
      r = fi(1) *(wd(1)*wt(1)) + fi(2) *(wd(3)*wt(1)) &
        + fi(3) *(wd(1)*wt(3)) + fi(4) *(wd(3)*wt(3)) &
        + fi(5) *(wd(1)*wt(2)) + fi(6) *(wd(3)*wt(2)) &
        + fi(7) *(wd(1)*wt(4)) + fi(8) *(wd(3)*wt(4)) &
        + fi(9) *(wd(2)*wt(1)) + fi(10)*(wd(4)*wt(1)) &
        + fi(11)*(wd(2)*wt(3)) + fi(12)*(wd(4)*wt(3)) &
        + fi(13)*(wd(2)*wt(2)) + fi(14)*(wd(4)*wt(2)) &
        + fi(15)*(wd(2)*wt(4)) + fi(16)*(wd(4)*wt(4))
   end function hc_ad

end module eos_wdm_av_interp
