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

program test_mix
   ! Unit tests for the WDM_AV additive-volume mixture module (WP3).
   use const_def, only: dp, kerg, amu, avo, crad, pi, me, planck_h
   use math_lib
   use auto_diff
   use eos_def
   use eos_wdm_av_def
   use eos_wdm_av_load
   use eos_wdm_av_interp, only: eos_wdm_av_eval_real
   use eos_wdm_av_mix
   use const_lib, only: const_init
   use chem_lib, only: chem_init
   use eos_lib, only: eos_init, eos_shutdown
   use utils_lib, only: mesa_error

   implicit none

   integer :: npass, nfail
   character(len=256), parameter :: ideal_dir = &
      '/Users/bauer39/work/mesas/evbauer/mesa/data/eosWDM_AV_data/ideal'
   character(len=256), parameter :: skye_dir = &
      '/Users/bauer39/work/mesas/evbauer/mesa/data/eosWDM_AV_data/skye'
   character(len=256), parameter :: skye_src = &
      '/Users/bauer39/work/mesas/evbauer/mesa/data/eosWDM_AV_data/source_skye/eosWDM_AV_c12.txt'

   real(dp) :: kmu
   integer :: info

   npass = 0
   nfail = 0
   kmu = kerg/amu

   call do_init()

   write(*,'(a)') '============================================================'
   write(*,'(a)') ' WDM_AV mixture unit tests (WP3)'
   write(*,'(a)') '============================================================'

   call load_ideal_tables()

   call test_convergence()          ! (a)
   call test_gas_vs_analytic()      ! (b)
   call test_identities()           ! (c)
   call test_fd_partials()          ! (d)
   call test_crosscheck()           ! (e)
   call test_first_law()            ! (f)
   call test_robustness()           ! (g)
   call test_timing()               ! (h)
   call test_electron_density_rule() ! (i) common-n_e mixing rule (D27)

   call eos_wdm_av_free_tables()
   call test_skye_c12()             ! realistic-physics chain

   write(*,'(a)') '------------------------------------------------------------'
   write(*,'(a,i0,a,i0,a)') ' SUMMARY: ', npass, ' passed, ', nfail, ' failed'
   write(*,'(a)') '============================================================'

   call eos_shutdown()
   if (nfail > 0) call mesa_error(__FILE__, __LINE__, 'test_mix: failures')

contains

   !----------------------------------------------------------------------------
   ! init / reporting helpers
   !----------------------------------------------------------------------------

   subroutine do_init()
      call const_init('/Users/bauer39/work/mesas/evbauer/mesa', info)
      if (info /= 0) call mesa_error(__FILE__, __LINE__, 'const_init')
      call math_init()
      call chem_init('isotopes.data', info)
      if (info /= 0) call mesa_error(__FILE__, __LINE__, 'chem_init')
      call eos_init(' ', .true., info)
      if (info /= 0) call mesa_error(__FILE__, __LINE__, 'eos_init')
   end subroutine do_init

   subroutine check(name, val, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: val, tol
      if (val <= tol .and. val == val) then
         write(*,'(a,a44,a,es11.3,a,es9.2)') 'PASS  ', name, '  ', val, '  tol ', tol
         npass = npass + 1
      else
         write(*,'(a,a44,a,es11.3,a,es9.2)') 'FAIL  ', name, '  ', val, '  tol ', tol
         nfail = nfail + 1
      end if
   end subroutine check

   subroutine check_true(name, ok)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      if (ok) then
         write(*,'(a,a44)') 'PASS  ', name
         npass = npass + 1
      else
         write(*,'(a,a44)') 'FAIL  ', name
         nfail = nfail + 1
      end if
   end subroutine check_true

   real(dp) function relerr(a, b)
      real(dp), intent(in) :: a, b
      relerr = abs(a - b)/max(abs(b), 1d-300)
   end function relerr

   !----------------------------------------------------------------------------
   ! table loading and composition setup
   !----------------------------------------------------------------------------

   subroutine load_ideal_tables()
      character(len=8) :: names(7)
      integer :: i, idx, ierr
      names = [character(len=8) :: 'h1','he4','c12','o16','ne20','ne22','mg24']
      do i = 1, 7
         call eos_wdm_av_load_isotope(trim(ideal_dir), trim(names(i)), idx, ierr)
         if (ierr /= 0 .or. idx <= 0) then
            write(*,'(a,a)') 'FAIL  could not load ideal table ', trim(names(i))
            call mesa_error(__FILE__, __LINE__, 'load_ideal_tables')
         end if
      end do
      write(*,'(a,i0,a)') '  loaded ', num_wdm_av_tables, ' ideal isotope tables'
   end subroutine load_ideal_tables

   integer function iso_idx(name) result(idx)
      character(len=*), intent(in) :: name
      idx = eos_wdm_av_table_index(name)
      if (idx <= 0) call mesa_error(__FILE__, __LINE__, 'iso_idx: '//trim(name))
   end function iso_idx

   ! build composition arrays from isotope names and mass fractions (A,Z from tables)
   subroutine make_comp(names, Xin, nsp, idxs, A, Z, X)
      character(len=*), intent(in) :: names(:)
      real(dp), intent(in) :: Xin(:)
      integer, intent(out) :: nsp
      integer, allocatable, intent(out) :: idxs(:)
      real(dp), allocatable, intent(out) :: A(:), Z(:), X(:)
      integer :: i
      nsp = size(names)
      allocate(idxs(nsp), A(nsp), Z(nsp), X(nsp))
      do i = 1, nsp
         idxs(i) = iso_idx(trim(names(i)))
         A(i) = wdm_av_tables(idxs(i))% A
         Z(i) = dble(wdm_av_tables(idxs(i))% Z)
         X(i) = Xin(i)
      end do
      X = X/sum(X)
   end subroutine make_comp

   !----------------------------------------------------------------------------
   ! analytic ideal-mixture reference formulas
   !----------------------------------------------------------------------------

   real(dp) function ideal_Pgas(rho, T, nsp, A, Z, X) result(P)
      real(dp), intent(in) :: rho, T, A(:), Z(:), X(:)
      integer, intent(in) :: nsp
      P = rho*kmu*T*sum(X(1:nsp)*(1d0 + Z(1:nsp))/A(1:nsp))
   end function ideal_Pgas

   real(dp) function ideal_egas(rho, T, nsp, A, Z, X) result(e)
      real(dp), intent(in) :: rho, T, A(:), Z(:), X(:)
      integer, intent(in) :: nsp
      e = 1.5d0*kmu*T*sum(X(1:nsp)*(1d0 + Z(1:nsp))/A(1:nsp))
   end function ideal_egas

   real(dp) function ideal_cv(nsp, A, Z, X) result(cv)
      real(dp), intent(in) :: A(:), Z(:), X(:)
      integer, intent(in) :: nsp
      cv = 1.5d0*kmu*sum(X(1:nsp)*(1d0 + Z(1:nsp))/A(1:nsp))
   end function ideal_cv

   ! exact Sackur-Tetrode mixture entropy (ions, g=1, + pooled electrons, g=2)
   real(dp) function sackur_s(rho, T, nsp, A, Z, X) result(s)
      real(dp), intent(in) :: rho, T, A(:), Z(:), X(:)
      integer, intent(in) :: nsp
      real(dp) :: ni, zi, ze, ne, acc
      integer :: i
      acc = 0d0
      ze = 2d0*pow(2d0*pi*me*kerg*T/pow2(planck_h), 1.5d0)
      do i = 1, nsp
         if (X(i) <= 0d0) cycle
         ni = X(i)*rho/(A(i)*amu)
         zi = pow(2d0*pi*A(i)*amu*kerg*T/pow2(planck_h), 1.5d0)
         acc = acc + ni*kerg*(log(zi/ni) + 2.5d0)
      end do
      ne = (rho/amu)*sum(X(1:nsp)*Z(1:nsp)/A(1:nsp))
      acc = acc + ne*kerg*(log(ze/ne) + 2.5d0)
      s = acc/rho
   end function sackur_s

   ! analytic ideal mixing entropy (ions + electrons) with ideal partial volumes
   real(dp) function ideal_dsmix(nsp, A, Z, X) result(ds)
      real(dp), intent(in) :: A(:), Z(:), X(:)
      integer, intent(in) :: nsp
      real(dp) :: sumB, ye, di, de
      integer :: i
      sumB = sum(X(1:nsp)*(1d0 + Z(1:nsp))/A(1:nsp))
      ye = sum(X(1:nsp)*Z(1:nsp)/A(1:nsp))
      di = 0d0; de = 0d0
      do i = 1, nsp
         if (X(i) <= 0d0) cycle
         di = di - (X(i)/A(i))*log(X(i)*(1d0 + Z(i))/(A(i)*sumB))
         if (Z(i) > 0d0 .and. ye > 0d0) &
            de = de + (X(i)*Z(i)/A(i))*log(Z(i)*sumB/((1d0 + Z(i))*ye))
      end do
      ds = kmu*(di + de)
   end function ideal_dsmix

   ! radiation pieces that pack_for_export adds
   real(dp) function rad_e(rho, T)
      real(dp), intent(in) :: rho, T
      rad_e = crad*pow4(T)/rho
   end function rad_e
   real(dp) function rad_s(rho, T)
      real(dp), intent(in) :: rho, T
      rad_s = 4d0*crad*pow3(T)/(3d0*rho)
   end function rad_s
   real(dp) function rad_p(T)
      real(dp), intent(in) :: T
      rad_p = crad*pow4(T)/3d0
   end function rad_p
   real(dp) function rad_cv(rho, T)
      real(dp), intent(in) :: rho, T
      rad_cv = 4d0*crad*pow3(T)/rho
   end function rad_cv

   !----------------------------------------------------------------------------
   ! (a) convergence
   !----------------------------------------------------------------------------

   subroutine test_convergence()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr, il, it, maxiter
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T, maxres
      real(dp) :: lrho(3) = [-8d0,-4d0,0d0], lt(3) = [4d0,6d0,8d0]
      integer :: npures

      write(*,'(a)') '--- (a) convergence ---'
      opts% mixing_entropy = 2
      maxiter = 0; maxres = 0d0
      call make_comp([character(len=8)::'o16','ne20','ne22','mg24'], &
                     [55d0,33d0,2d0,10d0], nsp, idxs, A, Z, X)
      do il = 1, 3
         do it = 1, 3
            rho = pow(10d0, lrho(il)); T = pow(10d0, lt(it))
            call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
            if (ierr /= 0) then
               write(*,'(a,2f7.2,a,i0,a,i0,a,es10.2)') '  o16/ne/mg at logRho,logT=', &
                  lrho(il), lt(it), ' ierr=', ierr, ' n_iter=', st% n_iter, ' resid=', st% residual
               nfail = nfail + 1
            end if
            maxiter = max(maxiter, st% n_iter)
            maxres = max(maxres, st% residual)
         end do
      end do
      call check('4-species residual < tol', maxres, 1d-12)
      write(*,'(a,i0)') '  4-species max n_iter = ', maxiter
      deallocate(idxs, A, Z, X)

      ! pure species: 1-2 iterations
      call make_comp([character(len=8)::'c12'], [1d0], nsp, idxs, A, Z, X)
      npures = 0
      do il = 1, 3
         do it = 1, 3
            rho = pow(10d0, lrho(il)); T = pow(10d0, lt(it))
            call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
            npures = max(npures, st% n_iter)
            if (ierr /= 0) nfail = nfail + 1
         end do
      end do
      call check_true('pure c12 converges in <= 2 iters', npures <= 2)
      write(*,'(a,i0)') '  pure c12 max n_iter = ', npures
      deallocate(idxs, A, Z, X)
   end subroutine test_convergence

   !----------------------------------------------------------------------------
   ! (b) gas quantities vs analytic ideal mixture
   !----------------------------------------------------------------------------

   subroutine test_gas_vs_analytic()
      write(*,'(a)') '--- (b) gas quantities vs analytic ideal mixture ---'
      call gas_case('pure c12          ', [character(len=8)::'c12'], [1d0])
      call gas_case('c12/o16 50/50     ', [character(len=8)::'c12','o16'], [50d0,50d0])
      call gas_case('he4/c12 50/50     ', [character(len=8)::'he4','c12'], [50d0,50d0])
      call gas_case('h1/he4 70/30      ', [character(len=8)::'h1','he4'], [70d0,30d0])
      call gas_case('o16/ne20/ne22/mg24', [character(len=8)::'o16','ne20','ne22','mg24'], &
                    [55d0,33d0,2d0,10d0])
      call gas_case('c12 + trace ne22  ', [character(len=8)::'c12','ne22'], [1d0,1d-6])
   end subroutine test_gas_vs_analytic

   subroutine gas_case(label, names, Xin)
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: names(:)
      real(dp), intent(in) :: Xin(:)
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv), res0(nv)
      integer :: nsp, ierr, ip
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T, Pgas, egas, sex, cv, Pt, Pr
      real(dp) :: chr_g, cht_g, cv_g, g3m1, g1_g, gad_g
      real(dp) :: e_lnP, e_e, e_s, e_cv, e_chr, e_cht, e_g1, e_gad, e_dsmix
      integer :: nok
      ! nondegenerate points (classical ideal entropy > 0); radiation is stripped
      ! analytically from every derived quantity, so it need not be negligible
      integer, parameter :: np = 6
      real(dp) :: lrho(np) = [-2d0,-1d0, 0d0, 0d0, 1d0, 2d0]
      real(dp) :: lt(np)   = [ 5d0, 6d0, 6d0, 7d0, 7d0, 7d0]

      call make_comp(names, Xin, nsp, idxs, A, Z, X)
      opts% mixing_entropy = 2
      nok = 0
      e_lnP=0d0; e_e=0d0; e_s=0d0; e_cv=0d0; e_chr=0d0; e_cht=0d0; e_g1=0d0; e_gad=0d0; e_dsmix=0d0
      do ip = 1, np
         rho = pow(10d0, lrho(ip)); T = pow(10d0, lt(ip))
         Pgas = ideal_Pgas(rho, T, nsp, A, Z, X)
         egas = ideal_egas(rho, T, nsp, A, Z, X)
         sex  = sackur_s(rho, T, nsp, A, Z, X)
         cv   = ideal_cv(nsp, A, Z, X)
         Pr = rad_p(T); Pt = Pgas + Pr

         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         nok = nok + 1

         ! strip radiation to recover gas-only quantities (dP/drho has no radiation
         ! term; dP/dT and Cv do; gamma1/grad_ad rebuilt from gas chi's and Cv)
         chr_g = res(i_chiRho)*Pt/Pgas
         cht_g = (res(i_chiT)*Pt - 4d0*Pr)/Pgas
         cv_g  = res(i_Cv) - rad_cv(rho,T)
         g3m1  = Pgas*cht_g/(rho*T*cv_g)
         g1_g  = chr_g + cht_g*g3m1
         gad_g = g3m1/g1_g

         e_lnP = max(e_lnP, relerr(res(i_lnPgas), log(Pgas)))
         e_e   = max(e_e,   relerr(exp(res(i_lnE)) - rad_e(rho,T), egas))
         e_s   = max(e_s,   relerr(exp(res(i_lnS)) - rad_s(rho,T), sex))
         e_cv  = max(e_cv,  relerr(cv_g, cv))
         e_chr = max(e_chr, relerr(chr_g, 1d0))
         e_cht = max(e_cht, relerr(cht_g, 1d0))
         e_g1  = max(e_g1,  relerr(g1_g, 5d0/3d0))
         e_gad = max(e_gad, relerr(gad_g, 0.4d0))

         ! mixing-entropy delta: s(mix=2) - s(mix=0) == analytic ideal mixing entropy
         opts% mixing_entropy = 0
         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res0, ddd, ddt, st, ierr)
         opts% mixing_entropy = 2
         if (ierr == 0) e_dsmix = max(e_dsmix, &
            abs((exp(res(i_lnS)) - exp(res0(i_lnS))) - ideal_dsmix(nsp,A,Z,X)) &
            / max(abs(ideal_dsmix(nsp,A,Z,X)), sex))
      end do
      write(*,'(a,a,a,i0,a)') '  ', label, '  (', nok, ' points)'
      call check_true('  points evaluated', nok >= 4)
      call check('  lnPgas vs analytic', e_lnP, 1d-6)
      call check('  e_gas  vs analytic', e_e, 1d-5)
      call check('  s_gas(mix2) vs Sackur-Tetrode', e_s, 1d-5)
      call check('  cV     vs analytic', e_cv, 1d-4)
      call check('  chiRho vs 1', e_chr, 1d-5)
      call check('  chiT   vs 1', e_cht, 1d-5)
      call check('  gamma1 vs 5/3', e_g1, 1d-4)
      call check('  grad_ad vs 0.4', e_gad, 1d-4)
      call check('  mixing delta = analytic', e_dsmix, 1d-5)
      deallocate(idxs, A, Z, X)
   end subroutine gas_case

   !----------------------------------------------------------------------------
   ! (c) identities: X_k=1 reproduces pure interpolant; split species
   !----------------------------------------------------------------------------

   subroutine test_identities()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      type(eos_wdm_av_derivs) :: d
      real(dp) :: res(nv), ddd(nv), ddt(nv), res2(nv)
      integer :: ierr, idxc, idxs1(1), idxs2(2)
      real(dp) :: A1(1), Z1(1), X1(1), A2(2), Z2(2), X2(2)
      real(dp) :: rho, T, xi, tau, Pp, ep, sp_, emax, a, b, dsmix_analytic, dsg

      write(*,'(a)') '--- (c) identities ---'
      opts% mixing_entropy = 1
      idxc = iso_idx('c12')
      rho = 1d2; T = 1d6
      xi = log(rho); tau = log(T)

      ! X_k = 1 reproduces the pure-species interpolant (gas-only)
      idxs1(1) = idxc; A1(1) = wdm_av_tables(idxc)% A; Z1(1) = dble(wdm_av_tables(idxc)% Z); X1(1) = 1d0
      call eos_wdm_av_get(opts, T, rho, 1, idxs1, A1, Z1, X1, res, ddd, ddt, st, ierr)
      call eos_wdm_av_eval_real(wdm_av_tables(idxc), xi, tau, d, ierr)
      Pp = rho*d% fx
      ep = d% f - d% ft
      sp_ = -d% ft/T
      emax = 0d0
      emax = max(emax, relerr(res(i_lnPgas), log(Pp)))
      emax = max(emax, relerr(exp(res(i_lnE)) - rad_e(rho,T), ep))
      emax = max(emax, relerr(exp(res(i_lnS)) - rad_s(rho,T), sp_))
      call check('X_k=1 reproduces pure interpolant', emax, 1d-9)

      ! split c12 into two identical copies (same table index), mixing_entropy=1
      a = 0.3d0; b = 0.7d0
      idxs2 = [idxc, idxc]
      A2 = wdm_av_tables(idxc)% A; Z2 = dble(wdm_av_tables(idxc)% Z); X2 = [a, b]
      call eos_wdm_av_get(opts, T, rho, 2, idxs2, A2, Z2, X2, res2, ddd, ddt, st, ierr)
      ! P and e unchanged vs pure; s increased by -(k/(A mu))(a ln a + b ln b)
      call check('split: Pgas unchanged', relerr(res2(i_lnPgas), res(i_lnPgas)), 1d-9)
      call check('split: e unchanged', relerr(exp(res2(i_lnE)), exp(res(i_lnE))), 1d-9)
      dsmix_analytic = -(kerg/(A2(1)*amu))*(a*log(a) + b*log(b))
      dsg = (exp(res2(i_lnS)) - rad_s(rho,T)) - (exp(res(i_lnS)) - rad_s(rho,T))
      call check('split: delta_s = -(k/Amu)(a ln a+b ln b)', relerr(dsg, dsmix_analytic), 1d-8)
      write(*,'(a,2es13.5)') '  split delta_s [code, analytic] = ', dsg, dsmix_analytic
   end subroutine test_identities

   !----------------------------------------------------------------------------
   ! (d) finite-difference verification of every d_dlnd/d_dlnT in res
   !----------------------------------------------------------------------------

   subroutine test_fd_partials()
      integer :: nsp
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: wA1, wB1, wA2, wB2
      write(*,'(a)') '--- (d) FD verification of res partials (THE derivative-pass test) ---'
      call make_comp([character(len=8)::'o16','ne20','ne22','mg24'], &
                     [55d0,33d0,2d0,10d0], nsp, idxs, A, Z, X)
      ! mid-cell point so the biquintic 3rd derivative is smooth across +-h
      call fd_at(nsp, idxs, A, Z, X, 1.06d0, 1.06d6, 2, wA2, wB2)
      call fd_at(nsp, idxs, A, Z, X, 1.06d0, 1.06d6, 1, wA1, wB1)
      write(*,'(a,2es11.3)') '  two-step  worst relerr [1st/2nd-deriv group, 3rd-deriv group] = ', wA2, wB2
      write(*,'(a,2es11.3)') '  one-step  worst relerr [1st/2nd-deriv group, 3rd-deriv group] = ', wA1, wB1
      call check('  all res partials FD-consistent (2 steps)', max(wA2, wB2), 1d-3)
      call check_true('  3rd-deriv partials much worse with 1 step (>=50x)', wB1 > 50d0*wB2)
      deallocate(idxs, A, Z, X)
   end subroutine test_fd_partials

   ! Richardson central-difference check of every res partial at (rho,T).
   ! Returns worst relative error over the "1st/2nd-derivative" group (wA) and the
   ! "3rd-derivative content" group {Cv,chiT,grad_ad} (wB).
   subroutine fd_at(nsp, idxs, A, Z, X, rho, T, nsteps, wA, wB, rule)
      integer, intent(in) :: nsp, idxs(nsp), nsteps
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp), rho, T
      real(dp), intent(out) :: wA, wB
      integer, intent(in), optional :: rule
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res0(nv), dd(nv), dt(nv)
      real(dp) :: fd_d, fd_t, er_d, er_t, atol
      integer :: ierr, k
      integer :: idx3(3)
      idx3 = [i_Cv, i_chiT, i_grad_ad]

      opts% mixing_entropy = 2
      opts% deriv_steps = nsteps
      opts% mixing_rule = 0
      if (present(rule)) opts% mixing_rule = rule

      call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res0, dd, dt, st, ierr)

      wA = 0d0; wB = 0d0
      do k = 1, i_latent_ddlnRho
         fd_d = fd_partial(opts, nsp, idxs, A, Z, X, rho, T, k, .true.)
         fd_t = fd_partial(opts, nsp, idxs, A, Z, X, rho, T, k, .false.)
         atol = 1d-6*max(abs(res0(k)), 1d0)     ! floor for near-zero partials
         er_d = abs(dd(k) - fd_d)/(abs(fd_d) + atol)
         er_t = abs(dt(k) - fd_t)/(abs(fd_t) + atol)
         if (any(idx3 == k)) then
            wB = max(wB, er_d, er_t)
         else
            wA = max(wA, er_d, er_t)
         end if
      end do
   end subroutine fd_at

   ! Richardson-extrapolated central difference of res(k) w.r.t. lnRho or lnT
   real(dp) function fd_partial(opts, nsp, idxs, A, Z, X, rho, T, k, wrt_rho) result(fd)
      type(eos_wdm_av_options), intent(in) :: opts
      integer, intent(in) :: nsp, idxs(nsp), k
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp), rho, T
      logical, intent(in) :: wrt_rho
      real(dp) :: f1, f2
      f1 = cdiff_res(opts, nsp, idxs, A, Z, X, rho, T, k, wrt_rho, 2d-3)
      f2 = cdiff_res(opts, nsp, idxs, A, Z, X, rho, T, k, wrt_rho, 1d-3)
      fd = (4d0*f2 - f1)/3d0          ! Richardson (error O(h^4))
   end function fd_partial

   real(dp) function cdiff_res(opts, nsp, idxs, A, Z, X, rho, T, k, wrt_rho, h) result(c)
      type(eos_wdm_av_options), intent(in) :: opts
      integer, intent(in) :: nsp, idxs(nsp), k
      real(dp), intent(in) :: A(nsp), Z(nsp), X(nsp), rho, T, h
      logical, intent(in) :: wrt_rho
      type(eos_wdm_av_state) :: st
      real(dp) :: rp(nv), rm(nv), dd(nv), dt(nv)
      integer :: ierr
      if (wrt_rho) then
         call eos_wdm_av_get(opts, T, rho*exp(h), nsp, idxs, A, Z, X, rp, dd, dt, st, ierr)
         call eos_wdm_av_get(opts, T, rho*exp(-h), nsp, idxs, A, Z, X, rm, dd, dt, st, ierr)
      else
         call eos_wdm_av_get(opts, T*exp(h), rho, nsp, idxs, A, Z, X, rp, dd, dt, st, ierr)
         call eos_wdm_av_get(opts, T*exp(-h), rho, nsp, idxs, A, Z, X, rm, dd, dt, st, ierr)
      end if
      c = (rp(k) - rm(k))/(2d0*h)
   end function cdiff_res

   !----------------------------------------------------------------------------
   ! (e) hand-derived Gibbs cross-check vs auto_diff (via res) and vs analytic
   !----------------------------------------------------------------------------

   subroutine test_crosscheck()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T, Pg, Pt, P_rho_res, P_T_res, cv_res, cp_res, e1, e2, e3, e4

      write(*,'(a)') '--- (e) hand-derived Gibbs vs auto_diff (F_AV, mixing off) ---'
      call make_comp([character(len=8)::'c12','o16'], [50d0,50d0], nsp, idxs, A, Z, X)
      opts% mixing_entropy = 0          ! F_AV only
      opts% debug_crosscheck = .true.   ! also prints internal hand-vs-autodiff relerr
      rho = 1d-2; T = 1d6                ! nondegenerate (classical entropy > 0)
      call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)

      ! extract GAS quantities from res (pack adds radiation): dP/drho has no
      ! radiation term (Prad indep of rho); dP/dT and Cv do.
      Pg = ideal_Pgas(rho, T, nsp, A, Z, X)
      Pt = Pg + rad_p(T)
      P_rho_res = res(i_chiRho)*Pt/rho                       ! = dP/drho (gas = total)
      P_T_res = res(i_chiT)*Pt/T - 4d0*rad_p(T)/T            ! subtract dPrad/dT
      cv_res = res(i_Cv) - rad_cv(rho,T)
      cp_res = cv_res + T*P_T_res*P_T_res/(rho*rho*P_rho_res)
      e1 = relerr(st% P_rho_mix, P_rho_res)
      e2 = relerr(st% P_T_mix, P_T_res)
      e3 = relerr(st% cv_mix, cv_res)
      e4 = relerr(st% cp_mix, cp_res)
      call check('cross-check P_rho (hand vs auto_diff)', e1, 1d-4)
      call check('cross-check P_T   (hand vs auto_diff)', e2, 1d-4)
      call check('cross-check c_V   (hand vs auto_diff)', e3, 1d-4)
      call check('cross-check c_P   (hand vs auto_diff)', e4, 1d-4)
      deallocate(idxs, A, Z, X)
   end subroutine test_crosscheck

   !----------------------------------------------------------------------------
   ! (f) first law de = T ds + (P/rho^2) drho
   !----------------------------------------------------------------------------

   subroutine test_first_law()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr, k
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho0, T0, h, worst
      real(dp) :: r1, r2, t1, t2, e_a, e_b, s_a, s_b, Pm, Tm, rhom, de, rhs
      real(dp) :: lrh, lth

      write(*,'(a)') '--- (f) first law de = T ds + (P/rho^2) drho ---'
      call make_comp([character(len=8)::'he4','c12'], [50d0,50d0], nsp, idxs, A, Z, X)
      opts% mixing_entropy = 2
      h = 1d-4
      worst = 0d0
      call random_reset()
      do k = 1, 20
         rho0 = pow(10d0, rand_range(-2d0, 4d0))
         T0 = pow(10d0, rand_range(4d0, 7d0))
         lrh = rand_range(-1d0, 1d0)*h
         lth = rand_range(-1d0, 1d0)*h
         r1 = rho0*exp(-lrh); r2 = rho0*exp(lrh)
         t1 = T0*exp(-lth); t2 = T0*exp(lth)
         call eos_wdm_av_get(opts, t1, r1, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         e_a = exp(res(i_lnE)); s_a = exp(res(i_lnS))
         call eos_wdm_av_get(opts, t2, r2, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         e_b = exp(res(i_lnE)); s_b = exp(res(i_lnS))
         ! midpoint state
         rhom = rho0; Tm = T0
         call eos_wdm_av_get(opts, Tm, rhom, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         Pm = exp(res(i_lnPgas)) + rad_p(Tm)
         de = e_b - e_a
         rhs = Tm*(s_b - s_a) + (Pm/(rhom*rhom))*(r2 - r1)
         worst = max(worst, abs(de - rhs)/max(abs(de), abs(rhs), 1d-99))
      end do
      call check('first law residual (O(h^2))', worst, 1d-6)
      deallocate(idxs, A, Z, X)
   end subroutine test_first_law

   !----------------------------------------------------------------------------
   ! (g) robustness
   !----------------------------------------------------------------------------

   subroutine test_robustness()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr, k, nfailpts, maxiter
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T

      write(*,'(a)') '--- (g) robustness ---'
      opts% mixing_entropy = 2

      ! near-table-edge point: expect a clean ierr, not a crash
      call make_comp([character(len=8)::'h1','he4'], [70d0,30d0], nsp, idxs, A, Z, X)
      rho = pow(10d0, 7d0); T = pow(10d0, 3.0d0)   ! above logRho max (6) for a partial density
      call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
      call check_true('near-edge point returns clean ierr', ierr /= 0)
      write(*,'(a,i0)') '    edge-case ierr = ', ierr
      deallocate(idxs, A, Z, X)

      ! extreme ratio: X = 1e-8 for a trace species
      call make_comp([character(len=8)::'c12','mg24'], [1d0,1d-8], nsp, idxs, A, Z, X)
      rho = 1d0; T = 1d6
      call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
      call check_true('extreme ratio (X=1e-8) converges', ierr == 0 .and. st% converged)
      deallocate(idxs, A, Z, X)

      ! random scan of 200 (rho,T,X) points
      call make_comp([character(len=8)::'o16','ne20','ne22','mg24'], &
                     [1d0,1d0,1d0,1d0], nsp, idxs, A, Z, X)
      nfailpts = 0; maxiter = 0
      call random_reset()
      do k = 1, 200
         rho = pow(10d0, rand_range(-6d0, 5d0))
         T = pow(10d0, rand_range(3.5d0, 8.5d0))
         X(1) = rand_range(0.1d0, 1d0); X(2) = rand_range(0.1d0, 1d0)
         X(3) = rand_range(0.1d0, 1d0); X(4) = rand_range(0.1d0, 1d0)
         X = X/sum(X)
         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) then
            nfailpts = nfailpts + 1
         else
            maxiter = max(maxiter, st% n_iter)
         end if
      end do
      write(*,'(a,i0,a,i0)') '    random scan: ', nfailpts, &
         ' point(s) declined (out-of-table ok); max n_iter = ', maxiter
      call check_true('random scan: no crashes, converged points iterate finitely', maxiter > 0)
      deallocate(idxs, A, Z, X)
   end subroutine test_robustness

   !----------------------------------------------------------------------------
   ! (h) timing
   !----------------------------------------------------------------------------

   subroutine test_timing()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr, k, ncall
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T, t0, t1, us_full, us_solve

      write(*,'(a)') '--- (h) timing (4-species) ---'
      call make_comp([character(len=8)::'o16','ne20','ne22','mg24'], &
                     [55d0,33d0,2d0,10d0], nsp, idxs, A, Z, X)
      rho = 1d0; T = 1d6
      ncall = 20000

      opts% mixing_entropy = 2; opts% deriv_steps = 0
      call cpu_time(t0)
      do k = 1, ncall
         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
      end do
      call cpu_time(t1)
      us_solve = (t1 - t0)/ncall*1d6

      opts% deriv_steps = 2
      call cpu_time(t0)
      do k = 1, ncall
         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
      end do
      call cpu_time(t1)
      us_full = (t1 - t0)/ncall*1d6

      write(*,'(a,f9.3,a)') '    real solve + pack (deriv_steps=0): ', us_solve, ' us/call'
      write(*,'(a,f9.3,a)') '    full (deriv_steps=2):              ', us_full, ' us/call'
      call check_true('timing completed', us_full > 0d0)
      deallocate(idxs, A, Z, X)
   end subroutine test_timing

   !----------------------------------------------------------------------------
   ! (i) common-electron-density mixing rule (Decision D27).  For fully ionized
   ! species rho_i = rho * y_e * A_i/Z_i (= n_e A_i m_u/Z_i), no solve.  Proves:
   !   - the volume identity sum X_i/rho_i = 1/rho holds to round-off;
   !   - the rule is EXACT for the ideal fully ionized tables (P, e, s to the
   !     same tolerances as the AV rows), once the mixing entropy is included;
   !   - the auto_diff mixture P equals the volume-weighted sum of species P_i;
   !   - every returned partial is FD-consistent (no Newton pass -> interpolant
   !     FD floor, better than the AV two-step pass);
   !   - the X_k=1 and split-species identities (for equal Z/A the ne densities
   !     coincide with the mixture rho, so these match the AV values);
   !   - the first law.
   !----------------------------------------------------------------------------

   subroutine test_electron_density_rule()
      integer :: nsp
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: wA, wB
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st, st2
      type(eos_wdm_av_derivs) :: d
      real(dp) :: res(nv), ddd(nv), ddt(nv), res2(nv)
      integer :: ierr, idxc, idxs1(1), idxs2(2)
      real(dp) :: A1(1), Z1(1), X1(1), A2(2), Z2(2), X2(2)
      real(dp) :: rho, T, xi, tau, Pp, ep, sp_, emax, fa, fb, dsmix_analytic, dsg
      real(dp) :: rho0, T0, h, worst, r1, r2, t1, t2, e_a, e_b, s_a, s_b
      real(dp) :: Pm, Tm, rhom, de, rhs, lrh, lth
      integer :: k

      write(*,'(a)') '--- (i) common-electron-density mixing rule (D27) ---'

      ! ideal-gas exactness, volume identity, auto_diff P = volume-weighted sum
      call ne_case('pure c12          ', [character(len=8)::'c12'], [1d0])
      call ne_case('c12/o16 50/50     ', [character(len=8)::'c12','o16'], [50d0,50d0])
      call ne_case('he4/c12 50/50     ', [character(len=8)::'he4','c12'], [50d0,50d0])
      call ne_case('h1/he4 70/30      ', [character(len=8)::'h1','he4'], [70d0,30d0])
      call ne_case('o16/ne20/ne22/mg24', [character(len=8)::'o16','ne20','ne22','mg24'], &
                   [55d0,33d0,2d0,10d0])

      ! FD verification of every res partial (same machinery as (d), rule=1)
      call make_comp([character(len=8)::'o16','ne20','ne22','mg24'], &
                     [55d0,33d0,2d0,10d0], nsp, idxs, A, Z, X)
      call fd_at(nsp, idxs, A, Z, X, 1.06d0, 1.06d6, 2, wA, wB, rule=1)
      write(*,'(a,2es11.3)') '  ne-rule worst relerr [1st/2nd-deriv group, 3rd-deriv group] = ', wA, wB
      call check('  all res partials FD-consistent (ne rule)', max(wA, wB), 1d-5)
      deallocate(idxs, A, Z, X)

      ! identities (mixing_entropy=1 so only the ion term enters)
      opts% mixing_entropy = 1
      opts% mixing_rule = wdm_av_rule_electron_density
      idxc = iso_idx('c12')
      rho = 1d2; T = 1d6
      xi = log(rho); tau = log(T)

      ! X_k = 1 reproduces the pure-species interpolant (rho_1 = rho for ne rule)
      idxs1(1) = idxc; A1(1) = wdm_av_tables(idxc)% A; Z1(1) = dble(wdm_av_tables(idxc)% Z); X1(1) = 1d0
      call eos_wdm_av_get(opts, T, rho, 1, idxs1, A1, Z1, X1, res, ddd, ddt, st, ierr)
      call eos_wdm_av_eval_real(wdm_av_tables(idxc), xi, tau, d, ierr)
      Pp = rho*d% fx
      ep = d% f - d% ft
      sp_ = -d% ft/T
      emax = 0d0
      emax = max(emax, relerr(res(i_lnPgas), log(Pp)))
      emax = max(emax, relerr(exp(res(i_lnE)) - rad_e(rho,T), ep))
      emax = max(emax, relerr(exp(res(i_lnS)) - rad_s(rho,T), sp_))
      call check('ne: X_k=1 reproduces pure interpolant', emax, 1d-9)

      ! split c12 into two identical copies: P,e unchanged; s up by the ion mixing term
      fa = 0.3d0; fb = 0.7d0
      idxs2 = [idxc, idxc]
      A2 = wdm_av_tables(idxc)% A; Z2 = dble(wdm_av_tables(idxc)% Z); X2 = [fa, fb]
      call eos_wdm_av_get(opts, T, rho, 2, idxs2, A2, Z2, X2, res2, ddd, ddt, st2, ierr)
      call check('ne: split Pgas unchanged', relerr(res2(i_lnPgas), res(i_lnPgas)), 1d-9)
      call check('ne: split e unchanged', relerr(exp(res2(i_lnE)), exp(res(i_lnE))), 1d-9)
      dsmix_analytic = -(kerg/(A2(1)*amu))*(fa*log(fa) + fb*log(fb))
      dsg = (exp(res2(i_lnS)) - rad_s(rho,T)) - (exp(res(i_lnS)) - rad_s(rho,T))
      call check('ne: split delta_s = -(k/Amu)(a ln a+b ln b)', relerr(dsg, dsmix_analytic), 1d-8)

      ! first law de = T ds + (P/rho^2) drho (ne rule)
      call make_comp([character(len=8)::'he4','c12'], [50d0,50d0], nsp, idxs, A, Z, X)
      opts% mixing_entropy = 2
      opts% mixing_rule = wdm_av_rule_electron_density
      h = 1d-4
      worst = 0d0
      call random_reset()
      do k = 1, 20
         rho0 = pow(10d0, rand_range(-2d0, 4d0))
         T0 = pow(10d0, rand_range(4d0, 7d0))
         lrh = rand_range(-1d0, 1d0)*h
         lth = rand_range(-1d0, 1d0)*h
         r1 = rho0*exp(-lrh); r2 = rho0*exp(lrh)
         t1 = T0*exp(-lth); t2 = T0*exp(lth)
         call eos_wdm_av_get(opts, t1, r1, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         e_a = exp(res(i_lnE)); s_a = exp(res(i_lnS))
         call eos_wdm_av_get(opts, t2, r2, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         e_b = exp(res(i_lnE)); s_b = exp(res(i_lnS))
         rhom = rho0; Tm = T0
         call eos_wdm_av_get(opts, Tm, rhom, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         Pm = exp(res(i_lnPgas)) + rad_p(Tm)
         de = e_b - e_a
         rhs = Tm*(s_b - s_a) + (Pm/(rhom*rhom))*(r2 - r1)
         worst = max(worst, abs(de - rhs)/max(abs(de), abs(rhs), 1d-99))
      end do
      call check('ne: first law residual (O(h^2))', worst, 1d-6)
      deallocate(idxs, A, Z, X)
   end subroutine test_electron_density_rule

   ! ideal-gas / identity checks for one mixture under the common-n_e rule
   subroutine ne_case(label, names, Xin)
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: names(:)
      real(dp), intent(in) :: Xin(:)
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: nsp, ierr, ip, i
      integer, allocatable :: idxs(:)
      real(dp), allocatable :: A(:), Z(:), X(:)
      real(dp) :: rho, T, Pgas, egas, sex, cv, Pt, Pr, vsum
      real(dp) :: chr_g, cht_g, cv_g, g3m1, g1_g, gad_g
      real(dp) :: e_lnP, e_e, e_s, e_cv, e_chr, e_cht, e_g1, e_gad, e_vid, e_pw
      integer :: nok
      integer, parameter :: np = 6
      real(dp) :: lrho(np) = [-2d0,-1d0, 0d0, 0d0, 1d0, 2d0]
      real(dp) :: lt(np)   = [ 5d0, 6d0, 6d0, 7d0, 7d0, 7d0]

      call make_comp(names, Xin, nsp, idxs, A, Z, X)
      opts% mixing_entropy = 2
      opts% mixing_rule = wdm_av_rule_electron_density
      nok = 0
      e_lnP=0d0; e_e=0d0; e_s=0d0; e_cv=0d0; e_chr=0d0; e_cht=0d0; e_g1=0d0; e_gad=0d0
      e_vid=0d0; e_pw=0d0
      do ip = 1, np
         rho = pow(10d0, lrho(ip)); T = pow(10d0, lt(ip))
         Pgas = ideal_Pgas(rho, T, nsp, A, Z, X)
         egas = ideal_egas(rho, T, nsp, A, Z, X)
         sex  = sackur_s(rho, T, nsp, A, Z, X)
         cv   = ideal_cv(nsp, A, Z, X)
         Pr = rad_p(T); Pt = Pgas + Pr

         call eos_wdm_av_get(opts, T, rho, nsp, idxs, A, Z, X, res, ddd, ddt, st, ierr)
         if (ierr /= 0) cycle
         nok = nok + 1

         ! volume identity sum X_i/rho_i = 1/rho (fixed-multiple densities)
         vsum = 0d0
         do i = 1, nsp
            vsum = vsum + X(i)/st% rho_i(i)
         end do
         e_vid = max(e_vid, abs(vsum - 1d0/rho)*rho)
         ! auto_diff mixture P (gas-only lnPgas) vs volume-weighted sum (st% P)
         e_pw = max(e_pw, relerr(exp(res(i_lnPgas)), st% P))

         chr_g = res(i_chiRho)*Pt/Pgas
         cht_g = (res(i_chiT)*Pt - 4d0*Pr)/Pgas
         cv_g  = res(i_Cv) - rad_cv(rho,T)
         g3m1  = Pgas*cht_g/(rho*T*cv_g)
         g1_g  = chr_g + cht_g*g3m1
         gad_g = g3m1/g1_g

         e_lnP = max(e_lnP, relerr(res(i_lnPgas), log(Pgas)))
         e_e   = max(e_e,   relerr(exp(res(i_lnE)) - rad_e(rho,T), egas))
         e_s   = max(e_s,   relerr(exp(res(i_lnS)) - rad_s(rho,T), sex))
         e_cv  = max(e_cv,  relerr(cv_g, cv))
         e_chr = max(e_chr, relerr(chr_g, 1d0))
         e_cht = max(e_cht, relerr(cht_g, 1d0))
         e_g1  = max(e_g1,  relerr(g1_g, 5d0/3d0))
         e_gad = max(e_gad, relerr(gad_g, 0.4d0))
      end do
      write(*,'(a,a,a,i0,a)') '  ', label, '  (', nok, ' points)'
      call check_true('  points evaluated', nok >= 4)
      call check('  volume identity sum X_i/rho_i=1/rho', e_vid, 1d-14)
      call check('  auto_diff P = volume-weighted sum', e_pw, 1d-12)
      call check('  lnPgas vs analytic (exact ideal)', e_lnP, 1d-6)
      call check('  e_gas  vs analytic (exact ideal)', e_e, 1d-5)
      call check('  s_gas  vs Sackur-Tetrode (exact ideal)', e_s, 1d-5)
      call check('  cV     vs analytic', e_cv, 1d-4)
      call check('  chiRho vs 1', e_chr, 1d-5)
      call check('  chiT   vs 1', e_cht, 1d-5)
      call check('  gamma1 vs 5/3', e_g1, 1d-4)
      call check('  grad_ad vs 0.4', e_gad, 1d-4)
      deallocate(idxs, A, Z, X)
   end subroutine ne_case

   !----------------------------------------------------------------------------
   ! Skye c12 realistic-physics chain
   !----------------------------------------------------------------------------

   subroutine test_skye_c12()
      type(eos_wdm_av_options) :: opts
      type(eos_wdm_av_state) :: st
      type(eos_wdm_av_derivs) :: d
      real(dp) :: res(nv), ddd(nv), ddt(nv)
      integer :: ierr, idx, i
      integer :: idxs(1)
      real(dp) :: A1(1), Z1(1), X1(1)
      real(dp) :: rho, T, xi, tau
      real(dp) :: lrho(3) = [2d0,4d0,6d0]
      real(dp) :: P_r, S_r, dPdR_r, dPdT_r, Cv_r, cp_r, chiR_r, chiT_r, gad_r
      real(dp) :: e_lnP, e_lnS, e_chr, e_cht, e_cv, e_gad
      real(dp) :: lrn, ltn

      write(*,'(a)') '--- Skye c12 realistic-physics chain ---'
      call eos_wdm_av_load_isotope(trim(skye_dir), 'c12', idx, ierr)
      if (ierr /= 0 .or. idx <= 0) then
         write(*,'(a)') '  SKIP: could not load Skye c12 table'
         return
      end if
      idxs(1) = idx; A1(1) = wdm_av_tables(idx)% A; Z1(1) = dble(wdm_av_tables(idx)% Z); X1(1) = 1d0
      opts% mixing_entropy = 2

      ! (i) off-node: compare to nearest node of the source table
      e_lnP=0d0; e_lnS=0d0; e_chr=0d0; e_cht=0d0; e_cv=0d0; e_gad=0d0
      do i = 1, 3
         rho = pow(10d0, lrho(i)); T = pow(10d0, 6.5d0)
         call eos_wdm_av_get(opts, T, rho, 1, idxs, A1, Z1, X1, res, ddd, ddt, st, ierr)
         if (ierr /= 0) then
            write(*,'(a,2f7.2,a,i0)') '  get failed at logRho,logT=', lrho(i), 6.5d0, ' ierr=', ierr
            nfail = nfail + 1; cycle
         end if
         call skye_ref(lrho(i), 6.5d0, P_r, S_r, dPdR_r, dPdT_r, Cv_r, ierr)
         if (ierr /= 0) cycle
         chiR_r = dPdR_r*rho/P_r; chiT_r = dPdT_r*T/P_r
         cp_r = Cv_r + T*dPdT_r*dPdT_r/(rho*rho*dPdR_r)
         gad_r = P_r*chiT_r/(rho*T*cp_r*chiR_r)
         e_lnP = max(e_lnP, relerr(res(i_lnPgas), log(P_r)))
         e_lnS = max(e_lnS, relerr(exp(res(i_lnS)) - rad_s(rho,T), S_r))
         e_chr = max(e_chr, relerr(res(i_chiRho), chiR_r))
         e_cht = max(e_cht, relerr(res(i_chiT), chiT_r))
         e_cv  = max(e_cv,  relerr(res(i_Cv) - rad_cv(rho,T), Cv_r))
         e_gad = max(e_gad, relerr(res(i_grad_ad), gad_r))
      end do
      write(*,'(a)') '  off-node vs nearest source node (loose, includes interp error):'
      call check('  Skye lnPgas', e_lnP, 3d-3)
      call check('  Skye lnS', e_lnS, 3d-3)
      call check('  Skye chiRho', e_chr, 2d-2)
      call check('  Skye chiT', e_cht, 2d-2)
      call check('  Skye Cv', e_cv, 2d-2)
      call check('  Skye grad_ad', e_gad, 2d-2)

      ! (ii) exactly at a node: P and S reproduced near roundoff by the interpolant
      lrn = 4d0; ltn = 7d0
      rho = pow(10d0, lrn); T = pow(10d0, ltn)
      xi = log(rho); tau = log(T)
      call eos_wdm_av_eval_real(wdm_av_tables(idx), xi, tau, d, ierr)
      call skye_ref(lrn, ltn, P_r, S_r, dPdR_r, dPdT_r, Cv_r, ierr)
      if (ierr == 0) then
         call check('  Skye at-node P (interp vs source)', relerr(rho*d% fx, P_r), 1d-4)
         call check('  Skye at-node S (interp vs source)', relerr(-d% ft/T, S_r), 1d-4)
      end if
   end subroutine test_skye_c12

   ! read the nearest-node reference row from the Skye source ASCII table
   subroutine skye_ref(lrho_q, lt_q, P, S, dPdR, dPdT, Cv, ierr)
      real(dp), intent(in) :: lrho_q, lt_q
      real(dp), intent(out) :: P, S, dPdR, dPdT, Cv
      integer, intent(out) :: ierr
      integer :: u, ios
      character(len=512) :: line
      real(dp) :: col(11), best, dist, lr, lt
      logical :: found
      ierr = 0
      found = .false.; best = 1d99
      P=0d0; S=0d0; dPdR=0d0; dPdT=0d0; Cv=0d0
      open(newunit=u, file=trim(skye_src), status='old', action='read', iostat=ios)
      if (ios /= 0) then; ierr = -1; return; end if
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read(line,*,iostat=ios) col
         if (ios /= 0) cycle
         if (col(3) /= col(3)) cycle      ! skip NaN rows
         lr = col(1); lt = col(2)
         dist = (lr - lrho_q)**2 + (lt - lt_q)**2
         if (dist < best) then
            best = dist; found = .true.
            ! columns: logRho logT P E F S Zbar dPdRho_T dPdT_Rho Cv eta
            P = col(3); S = col(6); dPdR = col(8); dPdT = col(9); Cv = col(10)
         end if
      end do
      close(u)
      if (.not. found) ierr = -1
   end subroutine skye_ref

   !----------------------------------------------------------------------------
   ! deterministic pseudo-random helpers
   !----------------------------------------------------------------------------

   subroutine random_reset()
      integer :: n
      integer, allocatable :: seed(:)
      call random_seed(size=n)
      allocate(seed(n))
      seed = 20261026
      call random_seed(put=seed)
      deallocate(seed)
   end subroutine random_reset

   real(dp) function rand_range(lo, hi) result(r)
      real(dp), intent(in) :: lo, hi
      real(dp) :: u
      call random_number(u)
      r = lo + (hi - lo)*u
   end function rand_range

end program test_mix
