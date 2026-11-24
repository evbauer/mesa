program direct_freeeos
   use iso_fortran_env, only: real64, iostat_end
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none

   integer, parameter :: dp = real64, neps = 20
   integer, parameter :: ifopt = 3, ifmod = 1, ifion = -2, kif = 2
   real(dp), parameter :: hydrogen = 0.70_dp, metals = 0.02_dp
   real(dp), parameter :: ln10 = 2.3025850929940455_dp
   real(dp), parameter :: atom_wgt(neps) = [ &
      1.007825_dp, 4.0026_dp, 12.0111_dp, 14.0067_dp, 15.9994_dp, &
      20.179_dp, 22.98977_dp, 24.305_dp, 26.9815_dp, 28.086_dp, &
      30.9738_dp, 32.06_dp, 35.453_dp, 39.948_dp, 40.08_dp, &
      47.9_dp, 51.996_dp, 54.938_dp, 55.847_dp, 58.71_dp ]
   character(len=2), parameter :: elements(neps) = [character(len=2) :: &
      'H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', &
      'P', 'S', 'Cl', 'Ar', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe', 'Ni']
   character(len=4096) :: mass_path, query_path, output_name, record
   character(len=2) :: element
   integer :: mass_unit, query_unit, output_unit, status, element_index, iter, queries
   real(dp) :: mass_frac(neps), eps(neps), requested_logrho, requested_logt, nan
   real(dp) :: logf, temperature, rho, logrho, logt, pressure_value, logp
   real(dp) :: cf, cp, qf, qp, sf, st, grada, rtp, rmue, fh2, fhe2, fhe3
   real(dp) :: xmu1, xmu3, eta, gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e
   real(dp) :: degeneracy(3), pressure(3), density(3), energy(3), enthalpy(3), entropy(3)
   external :: free_eos

   if (command_argument_count() /= 3) then
      error stop 'usage: ./build/direct_freeeos MASS_LIST QUERIES OUTPUT_BASENAME (under runs/)'
   end if
   call get_command_argument(1, mass_path)
   call get_command_argument(2, query_path)
   call get_command_argument(3, output_name)
   if (len_trim(output_name) == 0 .or. index(output_name, '/') /= 0 .or. &
       trim(output_name) == '.' .or. trim(output_name) == '..') error stop 'invalid output basename'

   open(newunit=mass_unit, file=trim(mass_path), status='old', action='read', iostat=status)
   if (status /= 0) error stop 'cannot open mass list'
   read(mass_unit, '(a)', iostat=status) record
   if (status /= 0) error stop 'cannot read mass list header'
   mass_frac(1) = hydrogen
   mass_frac(2) = 1.0_dp - hydrogen - metals
   do element_index = 3, neps
      read(mass_unit, *, iostat=status) element, mass_frac(element_index)
      if (status /= 0) error stop 'cannot read metal abundance'
      if (element /= elements(element_index)) error stop 'incorrect metal order'
   end do
   close(mass_unit)
   if (.not. all(ieee_is_finite(mass_frac))) error stop 'nonfinite composition'
   if (any(mass_frac < 0.0_dp)) error stop 'negative composition'
   if (abs(sum(mass_frac(3:neps)) - 1.0_dp) > 1.0e-12_dp) error stop 'metal pattern is not normalized'
   mass_frac(3:neps) = metals*mass_frac(3:neps)
   eps = mass_frac/atom_wgt

   open(newunit=query_unit, file=trim(query_path), status='old', action='read', iostat=status)
   if (status /= 0) error stop 'cannot open queries'
   open(newunit=output_unit, file='runs/'//trim(output_name), status='new', action='write', iostat=status)
   if (status /= 0) error stop 'cannot create output; run from reference/ with runs/ present and a new basename'
   write(output_unit, '(a)') '# Direct-call scaffold; FreeEOS release and success semantics require external verification.'
   write(output_unit, '(a,4(i0,1x))') '# ifopt ifmod ifion kif: ', ifopt, ifmod, ifion, kif
   write(output_unit, '(a,20(es25.16e3,1x))') '# mass_frac H He C N O Ne Na Mg Al Si P S Cl Ar Ca Ti Cr Mn Fe Ni: ', mass_frac
   write(output_unit, '(a,20(es25.16e3,1x))') '# eps: ', eps
   write(output_unit, '(a)') &
      '# log10rho_request log10T_request lnrho_return T rho P lnP '// &
      'degeneracy1 degeneracy2 degeneracy3 pressure1 pressure2 pressure3 '// &
      'density1 density2 density3 energy1 energy2 energy3 enthalpy1 enthalpy2 enthalpy3 '// &
      'entropy1 entropy2 entropy3 Cp grad_ad gamma1 gamma3 xmu1 xmu3 eta iter'
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   queries = 0
   do
      read(query_unit, '(a)', iostat=status) record
      if (status == iostat_end) exit
      if (status /= 0) error stop 'query read failed'
      record = adjustl(record)
      if (len_trim(record) == 0) cycle
      if (record(1:1) == '#') cycle
      read(record, *, iostat=status) requested_logrho, requested_logt
      if (status /= 0) error stop 'queries must contain log10rho log10T'
      if (.not. ieee_is_finite(requested_logrho) .or. &
          .not. ieee_is_finite(requested_logt)) error stop 'nonfinite query'
      logrho = requested_logrho*ln10
      logt = requested_logt*ln10
      temperature = nan
      rho = nan
      pressure_value = nan
      logp = nan
      degeneracy = nan
      pressure = nan
      density = nan
      energy = nan
      enthalpy = nan
      entropy = nan
      cp = nan
      grada = nan
      gamma1 = nan
      gamma3 = nan
      xmu1 = nan
      xmu3 = nan
      eta = nan
      iter = -huge(iter)
      call free_eos(ifopt, ifmod, ifion, kif, eps, neps, &
         logrho, logt, logf, temperature, rho, logrho, pressure_value, logp, &
         cf, cp, qf, qp, sf, st, grada, rtp, rmue, fh2, &
         fhe2, fhe3, xmu1, xmu3, eta, gamma1, gamma2, gamma3, &
         h2rat, h2plusrat, lambda, gamma_e, degeneracy, &
         pressure, density, energy, enthalpy, entropy, iter)
      write(output_unit, '(32(es25.16e3,1x),i0)') &
         requested_logrho, requested_logt, logrho, temperature, rho, pressure_value, logp, &
         degeneracy, pressure, density, energy, enthalpy, entropy, &
         cp, grada, gamma1, gamma3, xmu1, xmu3, eta, iter
      queries = queries + 1
   end do
   close(query_unit)
   close(output_unit)
   if (queries == 0) error stop 'no queries supplied'
   write(*, '(a,i0,a)') 'Wrote ', queries, ' raw calls; this is not a verified fresh-source accuracy result.'
end program direct_freeeos
