program helm_probe
   use iso_fortran_env, only: input_unit, output_unit, error_unit, iostat_end
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use const_def, only: dp, crad, kerg, avo, use_mesa_temp_cache
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_lib, only: chem_init, chem_get_iso_id, basic_composition_info
   use eos_def
   use eos_lib, only: eos_init, eos_shutdown, helmeos2_eval
   implicit none

   integer, parameter :: species = 20
   character(len=4), parameter :: isotope_names(species) = [character(len=4) :: &
      'h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'na23', 'mg24', 'al27', 'si28', &
      'p31', 's32', 'cl35', 'ar40', 'ca40', 'ti48', 'cr52', 'mn55', 'fe56', 'ni58']
   character(len=2), parameter :: element_names(species) = [character(len=2) :: &
      'H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', &
      'P', 'S', 'Cl', 'Ar', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe', 'Ni']
   integer :: chem_ids(species), ierr, status, stream_unit, data_unit, element_index, axis_index
   real(dp) :: mass_fractions(species), metal_sum, xh, xhe, metal, abar, zbar, z2bar, z53bar
   real(dp) :: ye, mass_correction, sumx, density, temperature, log_density, log_temperature
   real(dp) :: results(num_helm_results), electron(3), ions(3), radiation(3), total(3), coulomb(3)
   character(len=4096) :: mesa_directory, massfile, output_directory, record, cache_directory
   character(len=2) :: element
   logical :: off_table

   if (command_argument_count() /= 3) error stop 'usage: helm_probe MESA_DIR MASSFILE OUTPUT_DIR'
   call get_command_argument(1, mesa_directory)
   call get_command_argument(2, massfile)
   call get_command_argument(3, output_directory)
   call get_environment_variable('MESA_CACHES_DIR', cache_directory, status=status)
   if (status /= 0 .or. len_trim(cache_directory) == 0) error stop 'private cache root required'
   call get_environment_variable('MESA_TEMP_CACHES_DIR', cache_directory, status=status)
   if (status /= 0 .or. len_trim(cache_directory) == 0) error stop 'private temporary cache required'
   open(newunit=stream_unit, file='/dev/stdout', action='write')
   open(unit=output_unit, file='/dev/stderr', action='write')
   call const_init(trim(mesa_directory), ierr)
   if (ierr /= 0 .or. .not. use_mesa_temp_cache) error stop 'constant/cache initialization failed'
   call math_init()
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) error stop 'chemistry initialization failed'
   do element_index = 1, species
      chem_ids(element_index) = chem_get_iso_id(trim(isotope_names(element_index)))
      if (chem_ids(element_index) <= 0) error stop 'isotope unavailable'
   end do
   open(newunit=data_unit, file=trim(massfile), status='old', action='read')
   read(data_unit, '(a)') record
   mass_fractions(1:2) = [0.70_dp, 0.28_dp]
   do element_index = 3, species
      read(data_unit, *) element, mass_fractions(element_index)
      if (element /= element_names(element_index)) error stop 'metal order mismatch'
   end do
   close(data_unit)
   metal_sum = sum(mass_fractions(3:species))
   if (abs(metal_sum - 1.0_dp) > 1e-12_dp) error stop 'metal normalization mismatch'
   mass_fractions(3:species) = 0.02_dp * mass_fractions(3:species) / metal_sum
   call basic_composition_info(species, chem_ids, mass_fractions, xh, xhe, metal, &
      abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
   open(newunit=data_unit, file=trim(output_directory)//'/composition.tsv', status='replace')
   write(data_unit, '(a)') 'X Y Z abar zbar Ye crad kerg avo'
   write(data_unit, '(9(es26.17e3,1x))') xh, xhe, metal, abar, zbar, ye, crad, kerg, avo
   close(data_unit)
   call eos_init(' ', .true., ierr)
   if (ierr /= 0) error stop 'EOS initialization failed'
   open(newunit=data_unit, file=trim(output_directory)//'/density_axis.txt', status='replace')
   do axis_index = 1, eos_ht%imax
      write(data_unit, '(es26.17e3)') eos_ht%d(axis_index)
   end do
   close(data_unit)
   open(newunit=data_unit, file=trim(output_directory)//'/temperature_axis.txt', status='replace')
   do axis_index = 1, eos_ht%jmax
      write(data_unit, '(es26.17e3)') eos_ht%t(axis_index)
   end do
   close(data_unit)
   write(stream_unit, '(a)') 'rho T ierr ep_Cv ep_Cv_x ep_Cv_t ion_Cv ion_Cv_x ion_Cv_t '// &
      'rad_Cv rad_Cv_x rad_Cv_t coul_Cv coul_Cv_x coul_Cv_t total_Cv total_Cv_x total_Cv_t Gamma'
   flush(stream_unit)
   do
      read(input_unit, '(a)', iostat=status) record
      if (status == iostat_end) exit
      if (status /= 0) error stop 'input read failed'
      read(record, *, iostat=status) density, temperature
      if (status /= 0) error stop 'invalid query'
      if (.not. ieee_is_finite(density) .or. .not. ieee_is_finite(temperature)) error stop 'nonfinite query'
      if (min(density, temperature) <= 0) error stop 'nonpositive query'
      log_density = log10(density)
      log_temperature = log10(temperature)
      results = 0
      call helmeos2_eval(temperature, log_temperature, density, log_density, abar, zbar, &
         -1d99, -1d99, results, .false., .true., .true., off_table, ierr)
      if (off_table .or. ierr /= 0) error stop 'HELM query failed'
      electron = [results(h_deept), density*results(h_deepdt), temperature*results(h_deeptt)]
      ions = [results(h_deiont), density*results(h_deiondt), temperature*results(h_deiontt)]
      radiation = [results(h_deradt), density*results(h_deraddt), temperature*results(h_deradtt)]
      total = [results(h_cv), density*results(h_dcvdd), temperature*results(h_dcvdt)]
      coulomb = total - electron - ions - radiation
      if (.not. all(ieee_is_finite(total))) error stop 'nonfinite HELM result'
      write(stream_unit, '(2(es26.17e3,1x),i0,1x,16(es26.17e3,1x))') &
         density, temperature, ierr, electron, ions, radiation, coulomb, total, results(h_plasg)
      flush(stream_unit)
   end do
   call eos_shutdown()
   close(stream_unit)
end program helm_probe
