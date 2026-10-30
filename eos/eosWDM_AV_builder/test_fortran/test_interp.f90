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

program test_interp
   ! Unit tests for the WDM_AV per-species Hermite interpolant (WP2).
   use const_def, only: dp, kerg, amu, me, planck_h, pi
   use math_lib
   use auto_diff
   use helm_polynomials
   use eos_wdm_av_def
   use eos_wdm_av_interp
   use eos_wdm_av_load
   use utils_lib, only: mesa_error

   implicit none

   integer :: npass, nfail

   ! shared grid (logRho in [-10,6], logT in [3,9], 10 pts/decade; natural logs)
   integer, parameter :: NXI = 161, NTAU = 61
   real(dp) :: ln10, xi_min, dxi, tau_min, dtau

   ! generic test-function parameters (two exponentials, nontrivial in both axes)
   real(dp), parameter :: gc1 = 1.3d0, gp1 = 0.13d0, gq1 = 0.17d0
   real(dp), parameter :: gc2 = 0.8d0, gp2 = -0.09d0, gq2 = 0.11d0

   npass = 0
   nfail = 0
   ln10 = log(10d0)
   xi_min = ln10*(-10d0)
   dxi = ln10*0.1d0
   tau_min = ln10*3d0
   dtau = ln10*0.1d0

   write(*,'(a)') '============================================================'
   write(*,'(a)') ' WDM_AV interpolant unit tests (WP2)'
   write(*,'(a)') '============================================================'

   call test_cardinality()
   call test_generic(5)
   call test_generic(3)
   call test_edge_continuity()
   call test_physics(5)
   call test_physics(3)
   call test_isotope()
   call test_hdf5()

   write(*,'(a)') '------------------------------------------------------------'
   write(*,'(a,i0,a,i0,a)') ' SUMMARY: ', npass, ' passed, ', nfail, ' failed'
   write(*,'(a)') '============================================================'
   if (nfail > 0) call mesa_error(__FILE__, __LINE__, 'test_interp: failures')

contains

   !----------------------------------------------------------------------------
   ! reporting helpers
   !----------------------------------------------------------------------------

   subroutine check(name, val, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: val, tol
      if (val <= tol .and. val == val) then
         write(*,'(a,a40,a,es11.3,a,es9.2)') 'PASS  ', name, '  ', val, '  tol ', tol
         npass = npass + 1
      else
         write(*,'(a,a40,a,es11.3,a,es9.2)') 'FAIL  ', name, '  ', val, '  tol ', tol
         nfail = nfail + 1
      end if
   end subroutine check

   subroutine check_lt(name, a, b)
      ! pass if a < b (used for order5 vs order3 accuracy ordering)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: a, b
      if (a < b) then
         write(*,'(a,a40,a,es11.3,a,es11.3)') 'PASS  ', name, '  ', a, '  <  ', b
         npass = npass + 1
      else
         write(*,'(a,a40,a,es11.3,a,es11.3)') 'FAIL  ', name, '  ', a, '  >= ', b
         nfail = nfail + 1
      end if
   end subroutine check_lt

   real(dp) function relerr(a, b)
      real(dp), intent(in) :: a, b
      relerr = abs(a - b)/max(abs(b), 1d-300)
   end function relerr

   !----------------------------------------------------------------------------
   ! generic test function: g = c1 exp(p1 xi + q1 tau) + c2 exp(p2 xi + q2 tau)
   !----------------------------------------------------------------------------

   real(dp) function ipow(x, n) result(r)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      integer :: k
      r = 1d0
      do k = 1, n
         r = r*x
      end do
   end function ipow

   real(dp) function gder(a, b, xi, tau) result(r)
      ! d^(a+b) g / dxi^a dtau^b
      integer, intent(in) :: a, b
      real(dp), intent(in) :: xi, tau
      r = gc1*ipow(gp1,a)*ipow(gq1,b)*exp(gp1*xi + gq1*tau) &
        + gc2*ipow(gp2,a)*ipow(gq2,b)*exp(gp2*xi + gq2*tau)
   end function gder

   subroutine gen_nodal(xi, tau, n9)
      real(dp), intent(in) :: xi, tau
      real(dp), intent(out) :: n9(9)
      ! order: f, fx, ft, fxx, ftt, fxt, fxxt, fxtt, fxxtt
      n9(1) = gder(0,0,xi,tau)
      n9(2) = gder(1,0,xi,tau)
      n9(3) = gder(0,1,xi,tau)
      n9(4) = gder(2,0,xi,tau)
      n9(5) = gder(0,2,xi,tau)
      n9(6) = gder(1,1,xi,tau)
      n9(7) = gder(2,1,xi,tau)
      n9(8) = gder(1,2,xi,tau)
      n9(9) = gder(2,2,xi,tau)
   end subroutine gen_nodal

   subroutine gen_truth(xi, tau, d)
      real(dp), intent(in) :: xi, tau
      type(eos_wdm_av_derivs), intent(out) :: d
      d% f    = gder(0,0,xi,tau)
      d% fx   = gder(1,0,xi,tau)
      d% ft   = gder(0,1,xi,tau)
      d% fxx  = gder(2,0,xi,tau)
      d% fxt  = gder(1,1,xi,tau)
      d% ftt  = gder(0,2,xi,tau)
      d% fxxx = gder(3,0,xi,tau)
      d% fxxt = gder(2,1,xi,tau)
      d% fxtt = gder(1,2,xi,tau)
      d% fttt = gder(0,3,xi,tau)
   end subroutine gen_truth

   !----------------------------------------------------------------------------
   ! physics model: fully ionized ideal ions + ideal nondegenerate electrons
   !----------------------------------------------------------------------------

   subroutine phys_consts(A, Z, K, B, A0)
      real(dp), intent(in) :: A
      integer, intent(in) :: Z
      real(dp), intent(out) :: K, B, A0
      real(dp) :: Ci, Ce, a_i, a_e
      K = kerg/amu
      B = (1d0 + Z)/A
      Ci = pow(2d0*pi*(A*amu)*kerg/pow2(planck_h), 1.5d0)
      Ce = pow(2d0*pi*me*kerg/pow2(planck_h), 1.5d0)
      a_i = log(Ci) + log(A*amu) + 1d0
      a_e = log(2d0*Ce) - log(dble(Z)) + log(A*amu) + 1d0
      A0 = (a_i + Z*a_e)/A
   end subroutine phys_consts

   subroutine phys_nodal(A, Z, xi, tau, n9)
      real(dp), intent(in) :: A, xi, tau
      integer, intent(in) :: Z
      real(dp), intent(out) :: n9(9)
      real(dp) :: K, B, A0, E, P1
      call phys_consts(A, Z, K, B, A0)
      E = exp(tau)
      P1 = A0 + 1.5d0*B*tau - B*xi
      n9(1) = -K*E*P1               ! f
      n9(2) = K*B*E                 ! fx
      n9(3) = -K*E*(P1 + 1.5d0*B)   ! ft
      n9(4) = 0d0                   ! fxx
      n9(5) = -K*E*(P1 + 3d0*B)     ! ftt
      n9(6) = K*B*E                 ! fxt
      n9(7) = 0d0                   ! fxxt
      n9(8) = K*B*E                 ! fxtt
      n9(9) = 0d0                   ! fxxtt
   end subroutine phys_nodal

   subroutine phys_truth(A, Z, xi, tau, d)
      real(dp), intent(in) :: A, xi, tau
      integer, intent(in) :: Z
      type(eos_wdm_av_derivs), intent(out) :: d
      real(dp) :: K, B, A0, E, P1
      call phys_consts(A, Z, K, B, A0)
      E = exp(tau)
      P1 = A0 + 1.5d0*B*tau - B*xi
      d% f    = -K*E*P1
      d% fx   = K*B*E
      d% ft   = -K*E*(P1 + 1.5d0*B)
      d% fxx  = 0d0
      d% fxt  = K*B*E
      d% ftt  = -K*E*(P1 + 3d0*B)
      d% fxxx = 0d0
      d% fxxt = 0d0
      d% fxtt = K*B*E
      d% fttt = -K*E*(P1 + 4.5d0*B)
   end subroutine phys_truth

   !----------------------------------------------------------------------------
   ! table construction
   !----------------------------------------------------------------------------

   subroutine alloc_table(tbl, order, A, Z)
      type(eos_wdm_av_table), intent(inout) :: tbl
      integer, intent(in) :: order, Z
      real(dp), intent(in) :: A
      tbl% isotope = 'test'
      tbl% Z = Z
      tbl% A = A
      tbl% order = order
      tbl% nxi = NXI
      tbl% ntau = NTAU
      tbl% xi_min = xi_min
      tbl% dxi = dxi
      tbl% tau_min = tau_min
      tbl% dtau = dtau
      allocate(tbl% f(NXI,NTAU), tbl% fx(NXI,NTAU), tbl% ft(NXI,NTAU))
      allocate(tbl% fxx(NXI,NTAU), tbl% ftt(NXI,NTAU), tbl% fxt(NXI,NTAU))
      allocate(tbl% fxxt(NXI,NTAU), tbl% fxtt(NXI,NTAU), tbl% fxxtt(NXI,NTAU))
      tbl% loaded = .true.
   end subroutine alloc_table

   subroutine store_nodal(tbl, i, j, n9)
      type(eos_wdm_av_table), intent(inout) :: tbl
      integer, intent(in) :: i, j
      real(dp), intent(in) :: n9(9)
      tbl% f(i,j)     = n9(1)
      tbl% fx(i,j)    = n9(2)
      tbl% ft(i,j)    = n9(3)
      tbl% fxx(i,j)   = n9(4)
      tbl% ftt(i,j)   = n9(5)
      tbl% fxt(i,j)   = n9(6)
      tbl% fxxt(i,j)  = n9(7)
      tbl% fxtt(i,j)  = n9(8)
      tbl% fxxtt(i,j) = n9(9)
   end subroutine store_nodal

   subroutine build_generic(tbl, order)
      type(eos_wdm_av_table), intent(inout) :: tbl
      integer, intent(in) :: order
      integer :: i, j
      real(dp) :: xi, tau, n9(9)
      call alloc_table(tbl, order, 4d0, 2)
      do j = 1, NTAU
         tau = tau_min + (j-1)*dtau
         do i = 1, NXI
            xi = xi_min + (i-1)*dxi
            call gen_nodal(xi, tau, n9)
            call store_nodal(tbl, i, j, n9)
         end do
      end do
   end subroutine build_generic

   subroutine build_physics(tbl, order, A, Z)
      type(eos_wdm_av_table), intent(inout) :: tbl
      integer, intent(in) :: order, Z
      real(dp), intent(in) :: A
      integer :: i, j
      real(dp) :: xi, tau, n9(9)
      call alloc_table(tbl, order, A, Z)
      do j = 1, NTAU
         tau = tau_min + (j-1)*dtau
         do i = 1, NXI
            xi = xi_min + (i-1)*dxi
            call phys_nodal(A, Z, xi, tau, n9)
            call store_nodal(tbl, i, j, n9)
         end do
      end do
   end subroutine build_physics

   !----------------------------------------------------------------------------
   ! (a) basis cardinality
   !----------------------------------------------------------------------------

   subroutine test_cardinality()
      real(dp) :: e
      write(*,'(a)') '--- (a) basis cardinality (helm_polynomials) ---'
      ! quintic value basis psi0: psi0(0)=1, all others 0 at the two nodes
      e = 0d0
      e = max(e, abs(psi0(0d0)-1d0), abs(psi0(1d0)), abs(dpsi0(0d0)), abs(dpsi0(1d0)))
      e = max(e, abs(ddpsi0(0d0)), abs(ddpsi0(1d0)))
      call check('psi0 cardinality', e, 1d-12)
      ! psi1: first-derivative basis, dpsi1(0)=1
      e = 0d0
      e = max(e, abs(psi1(0d0)), abs(psi1(1d0)), abs(dpsi1(0d0)-1d0), abs(dpsi1(1d0)))
      e = max(e, abs(ddpsi1(0d0)), abs(ddpsi1(1d0)))
      call check('psi1 cardinality', e, 1d-12)
      ! psi2: second-derivative basis, ddpsi2(0)=1
      e = 0d0
      e = max(e, abs(psi2(0d0)), abs(psi2(1d0)), abs(dpsi2(0d0)), abs(dpsi2(1d0)))
      e = max(e, abs(ddpsi2(0d0)-1d0), abs(ddpsi2(1d0)))
      call check('psi2 cardinality', e, 1d-12)
      ! cubic value basis
      e = 0d0
      e = max(e, abs(xpsi0(0d0)-1d0), abs(xpsi0(1d0)), abs(xdpsi0(0d0)), abs(xdpsi0(1d0)))
      call check('xpsi0 cardinality', e, 1d-12)
      ! cubic first-derivative basis
      e = 0d0
      e = max(e, abs(xpsi1(0d0)), abs(xpsi1(1d0)), abs(xdpsi1(0d0)-1d0), abs(xdpsi1(1d0)))
      call check('xpsi1 cardinality', e, 1d-12)
   end subroutine test_cardinality

   !----------------------------------------------------------------------------
   ! generic-function tests: node reproduction, held-out accuracy, FD
   ! consistency, auto_diff vs real, edge continuity (both axes).
   !----------------------------------------------------------------------------

   subroutine test_generic(order)
      integer, intent(in) :: order
      type(eos_wdm_av_table) :: tbl
      type(eos_wdm_av_derivs) :: d, t
      integer :: i, j, ierr
      real(dp) :: xi, tau, emax
      character(len=8) :: tag

      write(tag,'(a,i0)') 'order ', order
      write(*,'(a,a,a)') '--- (b) generic function, ', trim(tag), ' ---'
      call build_generic(tbl, order)

      ! node reproduction: f and nodal derivatives exact at interior nodes
      emax = 0d0
      do j = 2, NTAU-1, 7
         do i = 2, NXI-1, 11
            xi = xi_min + (i-1)*dxi
            tau = tau_min + (j-1)*dtau
            call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
            call gen_truth(xi, tau, t)
            emax = max(emax, relerr(d% f, t% f), relerr(d% fx, t% fx), &
                       relerr(d% ft, t% ft), relerr(d% fxt, t% fxt))
            if (order == 5) emax = max(emax, relerr(d% fxx, t% fxx), relerr(d% ftt, t% ftt))
         end do
      end do
      call check('node reproduction (f,derivs)', emax, 1d-9)

      ! held-out mid-cell accuracy in f
      emax = 0d0
      do j = 2, NTAU-2, 5
         do i = 2, NXI-2, 9
            xi = xi_min + (i-1+0.5d0)*dxi
            tau = tau_min + (j-1+0.5d0)*dtau
            call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
            call gen_truth(xi, tau, t)
            emax = max(emax, relerr(d% f, t% f))
         end do
      end do
      if (order == 5) then
         call check('held-out f accuracy', emax, 1d-6)
      else
         call check('held-out f accuracy', emax, 1d-3)
      end if

      ! FD consistency of all returned derivatives, within a single cell
      call fd_consistency_generic(tbl, order)

      ! auto_diff vs real
      call ad_vs_real_generic(tbl)

      deallocate(tbl% f, tbl% fx, tbl% ft, tbl% fxx, tbl% ftt, tbl% fxt, &
                 tbl% fxxt, tbl% fxtt, tbl% fxxtt)
   end subroutine test_generic

   subroutine fd_consistency_generic(tbl, order)
      type(eos_wdm_av_table), intent(in) :: tbl
      integer, intent(in) :: order
      type(eos_wdm_av_derivs) :: d
      real(dp) :: xi0, tau0, h, fv(-2:2,-2:2)
      real(dp) :: fx, ft, fxx, ftt, fxt, fxxx, fttt, fxxt, fxtt
      real(dp) :: emax
      integer :: a, b, ierr
      ! a point at a cell center well inside the grid
      xi0 = xi_min + (40-1+0.5d0)*dxi
      tau0 = tau_min + (20-1+0.5d0)*dtau
      h = 0.02d0
      do a = -2, 2
         do b = -2, 2
            call eos_wdm_av_eval_real(tbl, xi0 + a*h, tau0 + b*h, d, ierr)
            fv(a,b) = d% f
         end do
      end do
      call eos_wdm_av_eval_real(tbl, xi0, tau0, d, ierr)
      fx  = (fv(1,0) - fv(-1,0))/(2d0*h)
      ft  = (fv(0,1) - fv(0,-1))/(2d0*h)
      fxx = (fv(1,0) - 2d0*fv(0,0) + fv(-1,0))/(h*h)
      ftt = (fv(0,1) - 2d0*fv(0,0) + fv(0,-1))/(h*h)
      fxt = (fv(1,1) - fv(1,-1) - fv(-1,1) + fv(-1,-1))/(4d0*h*h)
      fxxx = (fv(2,0) - 2d0*fv(1,0) + 2d0*fv(-1,0) - fv(-2,0))/(2d0*h*h*h)
      fttt = (fv(0,2) - 2d0*fv(0,1) + 2d0*fv(0,-1) - fv(0,-2))/(2d0*h*h*h)
      fxxt = ((fv(1,1) - 2d0*fv(0,1) + fv(-1,1)) - (fv(1,-1) - 2d0*fv(0,-1) + fv(-1,-1)))/(2d0*h*h*h)
      fxtt = ((fv(1,1) - 2d0*fv(1,0) + fv(1,-1)) - (fv(-1,1) - 2d0*fv(-1,0) + fv(-1,-1)))/(2d0*h*h*h)
      emax = 0d0
      emax = max(emax, relerr(d% fx, fx), relerr(d% ft, ft))
      emax = max(emax, relerr(d% fxx, fxx), relerr(d% ftt, ftt), relerr(d% fxt, fxt))
      emax = max(emax, relerr(d% fxxx, fxxx), relerr(d% fttt, fttt))
      emax = max(emax, relerr(d% fxxt, fxxt), relerr(d% fxtt, fxtt))
      call check('FD consistency (1st-3rd derivs)', emax, 1d-3)
   end subroutine fd_consistency_generic

   subroutine ad_vs_real_generic(tbl)
      type(eos_wdm_av_table), intent(in) :: tbl
      type(eos_wdm_av_derivs) :: d
      type(auto_diff_real_2var_order3) :: xi_ad, tau_ad, f_ad, fx_ad, ft_ad, fxx_ad
      real(dp) :: xi0, tau0, emax
      integer :: ierr
      xi0 = xi_min + (40-1+0.37d0)*dxi
      tau0 = tau_min + (20-1+0.61d0)*dtau
      call eos_wdm_av_eval_real(tbl, xi0, tau0, d, ierr)
      xi_ad = xi0
      xi_ad% d1val2 = 1d0     ! var2 = xi
      tau_ad = tau0
      tau_ad% d1val1 = 1d0    ! var1 = tau
      call eos_wdm_av_eval_ad(tbl, xi_ad, tau_ad, f_ad, fx_ad, ft_ad, fxx_ad, ierr)
      emax = 0d0
      ! f output carries all derivatives through third order
      emax = max(emax, relerr(f_ad% val, d% f))
      emax = max(emax, relerr(f_ad% d1val2, d% fx))
      emax = max(emax, relerr(f_ad% d1val1, d% ft))
      emax = max(emax, relerr(f_ad% d2val2, d% fxx))
      emax = max(emax, relerr(f_ad% d1val1_d1val2, d% fxt))
      emax = max(emax, relerr(f_ad% d2val1, d% ftt))
      emax = max(emax, relerr(f_ad% d3val2, d% fxxx))
      emax = max(emax, relerr(f_ad% d1val1_d2val2, d% fxxt))   ! d/dtau d2/dxi2
      emax = max(emax, relerr(f_ad% d2val1_d1val2, d% fxtt))   ! d2/dtau2 d/dxi
      emax = max(emax, relerr(f_ad% d3val1, d% fttt))
      call check('auto_diff f vs real (all 3rd order)', emax, 1d-8)
      ! fx output vs real third derivatives
      emax = 0d0
      emax = max(emax, relerr(fx_ad% val, d% fx))
      emax = max(emax, relerr(fx_ad% d1val2, d% fxx))
      emax = max(emax, relerr(fx_ad% d1val1, d% fxt))
      emax = max(emax, relerr(fx_ad% d2val2, d% fxxx))
      emax = max(emax, relerr(fx_ad% d1val1_d1val2, d% fxxt))
      emax = max(emax, relerr(fx_ad% d2val1, d% fxtt))
      call check('auto_diff fx vs real', emax, 1d-8)
      ! fxx output vs real third derivatives
      emax = 0d0
      emax = max(emax, relerr(fxx_ad% val, d% fxx))
      emax = max(emax, relerr(fxx_ad% d1val2, d% fxxx))
      emax = max(emax, relerr(fxx_ad% d1val1, d% fxxt))
      call check('auto_diff fxx vs real', emax, 1d-8)
   end subroutine ad_vs_real_generic

   subroutine test_edge_continuity()
      ! second derivatives: biquintic (order 5) is C2 so they are continuous
      ! across cell edges; bicubic (order 3) is only C1 so they jump.
      type(eos_wdm_av_table) :: t5, t3
      real(dp) :: j5x, j5t, j3x, j3t
      write(*,'(a)') '--- edge continuity of second derivatives ---'
      call build_generic(t5, 5)
      call build_generic(t3, 3)
      call edge_jumps(t5, j5x, j5t)
      call edge_jumps(t3, j3x, j3t)
      call check('order5 fxx continuous across xi-edge', j5x, 1d-7)
      call check('order5 ftt continuous across tau-edge', j5t, 1d-7)
      call check_lt('order3 fxx jumps (>> order5)', max(50d0*j5x, 1d-9), j3x)
      call check_lt('order3 ftt jumps (>> order5)', max(50d0*j5t, 1d-9), j3t)
      write(*,'(a,2es11.3,a,2es11.3)') '  jumps order5 (fxx,ftt)=', j5x, j5t, &
         '   order3=', j3x, j3t
      deallocate(t5% f, t5% fx, t5% ft, t5% fxx, t5% ftt, t5% fxt, t5% fxxt, t5% fxtt, t5% fxxtt)
      deallocate(t3% f, t3% fx, t3% ft, t3% fxx, t3% ftt, t3% fxt, t3% fxxt, t3% fxtt, t3% fxxtt)
   end subroutine test_edge_continuity

   subroutine edge_jumps(tbl, jump_fxx, jump_ftt)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(out) :: jump_fxx, jump_ftt
      type(eos_wdm_av_derivs) :: dm, dp_
      real(dp) :: xi_node, tau_node, xi_mid, tau_mid, delta
      integer :: inode, jnode, ierr
      inode = 60
      jnode = 30
      xi_node = xi_min + (inode-1)*dxi
      tau_node = tau_min + (jnode-1)*dtau
      xi_mid = xi_min + (inode-1+0.5d0)*dxi
      tau_mid = tau_min + (jnode-1+0.5d0)*dtau
      delta = 1d-9      ! fraction of the spacing; small enough to isolate the jump
      call eos_wdm_av_eval_real(tbl, xi_node - delta*dxi, tau_mid, dm, ierr)
      call eos_wdm_av_eval_real(tbl, xi_node + delta*dxi, tau_mid, dp_, ierr)
      jump_fxx = abs(dp_% fxx - dm% fxx)/max(abs(dm% fxx), abs(dp_% fxx), 1d-99)
      call eos_wdm_av_eval_real(tbl, xi_mid, tau_node - delta*dtau, dm, ierr)
      call eos_wdm_av_eval_real(tbl, xi_mid, tau_node + delta*dtau, dp_, ierr)
      jump_ftt = abs(dp_% ftt - dm% ftt)/max(abs(dm% ftt), abs(dp_% ftt), 1d-99)
   end subroutine edge_jumps

   !----------------------------------------------------------------------------
   ! physics-model tests: thermodynamic quantities held out of the table
   !----------------------------------------------------------------------------

   subroutine test_physics(order)
      integer, intent(in) :: order
      type(eos_wdm_av_table) :: tbl
      type(eos_wdm_av_derivs) :: d, t
      integer :: i, j, ierr
      real(dp) :: xi, tau, rho, Tk
      real(dp) :: emax_f, emax_P, emax_s, emax_cv, emax_chr, emax_cht
      real(dp), parameter :: A = 4.0026d0
      integer, parameter :: Z = 2

      write(*,'(a,i0,a)') '--- (b) physics model he4, order ', order, ' ---'
      call build_physics(tbl, order, A, Z)

      emax_f = 0d0; emax_P = 0d0; emax_s = 0d0; emax_cv = 0d0; emax_chr = 0d0; emax_cht = 0d0
      do j = 2, NTAU-2, 5
         do i = 2, NXI-2, 9
            xi = xi_min + (i-1+0.5d0)*dxi
            tau = tau_min + (j-1+0.5d0)*dtau
            rho = exp(xi); Tk = exp(tau)
            call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
            call phys_truth(A, Z, xi, tau, t)
            ! f crosses zero (f = -K T P1), so scale its error by the thermal
            ! energy floor K*T to avoid a meaningless blow-up near the zero line
            emax_f = max(emax_f, abs(d% f - t% f)/(abs(t% f) + (kerg/amu)*Tk))
            ! P = rho*fx
            emax_P = max(emax_P, relerr(rho*d% fx, rho*t% fx))
            ! s = -ft/T
            emax_s = max(emax_s, relerr(-d% ft/Tk, -t% ft/Tk))
            ! c_V = (ft - ftt)/T  (ft and ftt are large and nearly equal here, so
            ! this is a genuine cancellation; see Decision D16 on carrying s)
            emax_cv = max(emax_cv, relerr((d% ft - d% ftt)/Tk, (t% ft - t% ftt)/Tk))
            ! chi_rho = 1 + fxx/fx
            emax_chr = max(emax_chr, relerr(1d0 + d% fxx/d% fx, 1d0 + t% fxx/t% fx))
            ! chi_T = fxt/fx
            emax_cht = max(emax_cht, relerr(d% fxt/d% fx, t% fxt/t% fx))
         end do
      end do
      call check('held-out f   error (/kT scale)', emax_f, merge(1d-5,1d-2,order==5))
      call check('held-out P   error', emax_P, merge(1d-6,1d-3,order==5))
      call check('held-out s   error', emax_s, merge(1d-6,1d-3,order==5))
      call check('held-out cV  error (cancellation)', emax_cv, merge(1d-4,1d-1,order==5))
      call check('held-out chi_rho error', emax_chr, 1d-9)
      call check('held-out chi_T   error', emax_cht, merge(1d-6,1d-3,order==5))

      ! node reproduction to roundoff
      emax_f = 0d0
      do j = 2, NTAU-1, 7
         do i = 2, NXI-1, 11
            xi = xi_min + (i-1)*dxi
            tau = tau_min + (j-1)*dtau
            call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
            call phys_truth(A, Z, xi, tau, t)
            emax_f = max(emax_f, relerr(d% f, t% f), relerr(d% fx, t% fx), &
                         relerr(d% ft, t% ft))
         end do
      end do
      call check('node reproduction', emax_f, 1d-10)

      deallocate(tbl% f, tbl% fx, tbl% ft, tbl% fxx, tbl% ftt, tbl% fxt, &
                 tbl% fxxt, tbl% fxtt, tbl% fxxtt)
   end subroutine test_physics

   !----------------------------------------------------------------------------
   ! (c) isotope rescaling: ne20 table evaluated for ne22 vs the ne22 model
   !----------------------------------------------------------------------------

   subroutine test_isotope()
      type(eos_wdm_av_table) :: tbl
      type(eos_wdm_av_derivs) :: d, t
      integer :: i, j, ierr
      real(dp) :: xi, tau, rho, Tk, r, lnr, emax, emax_node
      real(dp), parameter :: A_ne20 = 19.9924401762d0
      real(dp), parameter :: A_ne22 = 21.99138511d0
      integer, parameter :: Z = 10

      write(*,'(a)') '--- (c) isotope rescaling ne20 -> ne22 ---'
      call build_physics(tbl, 5, A_ne20, Z)

      r = A_ne20/A_ne22
      lnr = log(r)

      ! held-out interior points: rescaled vs analytic ne22
      emax = 0d0
      do j = 10, NTAU-10, 7
         do i = 30, NXI-30, 11
            xi = xi_min + (i-1+0.5d0)*dxi
            tau = tau_min + (j-1+0.5d0)*dtau
            rho = exp(xi); Tk = exp(tau)
            call eos_wdm_av_eval_species_real(tbl, A_ne22, xi, tau, d, ierr)
            if (ierr /= 0) cycle
            call phys_truth(A_ne22, Z, xi, tau, t)
            ! f scaled by thermal energy (f crosses zero); P, s, chi_T are clean
            emax = max(emax, abs(d% f - t% f)/(abs(t% f) + (kerg/amu)*Tk))
            emax = max(emax, relerr(rho*d% fx, rho*t% fx))
            emax = max(emax, relerr(-d% ft/Tk, -t% ft/Tk))
            emax = max(emax, relerr(d% fxt/d% fx, t% fxt/t% fx))
         end do
      end do
      call check('rescaled interior f,P,s,chiT', emax, 1d-5)

      ! at points that map to table nodes (xi + lnr = node), result is exact
      emax_node = 0d0
      do j = 15, NTAU-15, 9
         do i = 50, NXI-15, 13
            xi = (xi_min + (i-1)*dxi) - lnr    ! so xi + lnr hits node i
            tau = tau_min + (j-1)*dtau         ! tau at a node
            rho = exp(xi); Tk = exp(tau)
            call eos_wdm_av_eval_species_real(tbl, A_ne22, xi, tau, d, ierr)
            if (ierr /= 0) cycle
            call phys_truth(A_ne22, Z, xi, tau, t)
            emax_node = max(emax_node, relerr(d% f, t% f))
            emax_node = max(emax_node, relerr(rho*d% fx, rho*t% fx))
            emax_node = max(emax_node, relerr(d% f - d% ft, t% f - t% ft))  ! e = f - ft
         end do
      end do
      call check('rescaled at mapped nodes (P,e exact)', emax_node, 1d-10)

      ! auto_diff species rescaling agrees with real
      call isotope_ad(tbl, A_ne22)

      deallocate(tbl% f, tbl% fx, tbl% ft, tbl% fxx, tbl% ftt, tbl% fxt, &
                 tbl% fxxt, tbl% fxtt, tbl% fxxtt)
   end subroutine test_isotope

   subroutine isotope_ad(tbl, A_species)
      type(eos_wdm_av_table), intent(in) :: tbl
      real(dp), intent(in) :: A_species
      type(eos_wdm_av_derivs) :: d
      type(auto_diff_real_2var_order3) :: xi_ad, tau_ad, f_ad, fx_ad, ft_ad, fxx_ad
      real(dp) :: xi0, tau0, emax, scale
      integer :: ierr
      xi0 = xi_min + (70-1+0.4d0)*dxi
      tau0 = tau_min + (25-1+0.3d0)*dtau
      call eos_wdm_av_eval_species_real(tbl, A_species, xi0, tau0, d, ierr)
      xi_ad = xi0; xi_ad% d1val2 = 1d0
      tau_ad = tau0; tau_ad% d1val1 = 1d0
      call eos_wdm_av_eval_species_ad(tbl, A_species, xi_ad, tau_ad, f_ad, fx_ad, ft_ad, fxx_ad, ierr)
      ! AD and real do nearly identical arithmetic, so compare absolute
      ! differences against the overall magnitude scale (ft crosses zero, so a
      ! pointwise relative error would be meaningless there)
      scale = max(abs(d% f), abs(d% fx), abs(d% ft), abs(d% ftt), abs(d% fxt), 1d-99)
      emax = 0d0
      emax = max(emax, abs(f_ad% val - d% f))
      emax = max(emax, abs(f_ad% d1val2 - d% fx))
      emax = max(emax, abs(f_ad% d1val1 - d% ft))
      emax = max(emax, abs(f_ad% d2val2 - d% fxx))
      emax = max(emax, abs(f_ad% d2val1 - d% ftt))
      emax = max(emax, abs(f_ad% d1val1_d1val2 - d% fxt))
      emax = max(emax, abs(ft_ad% val - d% ft))
      emax = max(emax, abs(fx_ad% val - d% fx))
      call check('isotope auto_diff vs real (/scale)', emax/scale, 1d-10)
   end subroutine isotope_ad

   !----------------------------------------------------------------------------
   ! (d) HDF5 loader
   !----------------------------------------------------------------------------

   subroutine test_hdf5()
      type(eos_wdm_av_table) :: tbl
      type(eos_wdm_av_derivs) :: d, t
      integer :: ierr, idx, i, j
      real(dp) :: xi, tau, rho, Tk, emax
      logical :: real_exists
      character(len=256) :: real_file, phys_file, mesa_dir_env
      real(dp), parameter :: A = 4.0026d0
      integer, parameter :: Z = 2

      write(*,'(a)') '--- (d) HDF5 loader ---'

      ! always generate the synthetic files (needed for the bad-version and
      ! lazy-load sub-tests, and as a fallback for the physics comparison)
      call execute_command_line('python3 -I make_test_h5.py scratch/test_data', wait=.true., exitstat=ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'FAIL  make_test_h5.py failed to run'
         nfail = nfail + 1
         return
      end if

      ! prefer the builder-produced file for the physics comparison, if present
      call get_environment_variable('MESA_DIR', mesa_dir_env)
      real_file = trim(mesa_dir_env)//'/data/eosWDM_AV_data/ideal/eosWDM_AV_he4.h5'
      inquire(file=trim(real_file), exist=real_exists)
      if (real_exists) then
         phys_file = trim(real_file)
         write(*,'(a)') '  physics comparison uses builder-produced file'
      else
         phys_file = 'scratch/test_data/eosWDM_AV_he4.h5'
         write(*,'(a)') '  physics comparison uses synthetic make_test_h5.py file'
      end if

      ! read the valid file
      call eos_wdm_av_read_table(trim(phys_file), tbl, ierr)
      call check('read valid file ierr==0', dble(abs(ierr)), 0.5d0)
      if (ierr == 0) then
         write(*,'(a,i0,a,f9.4,a,i0)') '  attrs: Z=', tbl% Z, ' A=', tbl% A, ' order=', tbl% order
         write(*,'(a,2i5,a,4es12.4)') '  grid: nxi,ntau=', tbl% nxi, tbl% ntau, &
            ' xi_min,dxi,tau_min,dtau=', tbl% xi_min, tbl% dxi, tbl% tau_min, tbl% dtau
         ! compare interpolated f, P, s at interior points vs analytic model
         emax = 0d0
         do j = 10, tbl% ntau-10, 13
            do i = 20, tbl% nxi-20, 23
               xi = tbl% xi_min + (i-1+0.5d0)*tbl% dxi
               tau = tbl% tau_min + (j-1+0.5d0)*tbl% dtau
               rho = exp(xi); Tk = exp(tau)
               call eos_wdm_av_eval_real(tbl, xi, tau, d, ierr)
               call phys_truth(A, Z, xi, tau, t)
               emax = max(emax, abs(d% f - t% f)/(abs(t% f) + (kerg/amu)*Tk))
               emax = max(emax, relerr(rho*d% fx, rho*t% fx))
               emax = max(emax, relerr(-d% ft/Tk, -t% ft/Tk))
            end do
         end do
         call check('loaded table f,P,s vs analytic', emax, 1d-5)
         deallocate(tbl% f, tbl% fx, tbl% ft, tbl% fxx, tbl% ftt, tbl% fxt, &
                    tbl% fxxt, tbl% fxtt, tbl% fxxtt)
      end if

      ! wrong format_version must be refused
      call eos_wdm_av_read_table('scratch/test_data/eosWDM_AV_badver.h5', tbl, ierr)
      call check('bad format_version refused (ierr/=0)', merge(0d0,1d0,ierr/=0), 0.5d0)

      ! lazy-load path through the module store
      call eos_wdm_av_load_isotope('scratch/test_data', 'he4', idx, ierr)
      call check('load_isotope returns index', merge(0d0,1d0,(idx>0 .and. ierr==0)), 0.5d0)
      ! second call returns the same cached index
      call eos_wdm_av_load_isotope('scratch/test_data', 'he4', idx, ierr)
      call check('load_isotope cached', merge(0d0,1d0,(idx>0 .and. ierr==0)), 0.5d0)
      call eos_wdm_av_free_tables()
   end subroutine test_hdf5

end program test_interp
