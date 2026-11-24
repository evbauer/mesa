program mesa_probe
   use iso_fortran_env, only: input_unit, output_unit, error_unit, iostat_end
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use const_def, only: dp, use_mesa_temp_cache
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use eos_def
   use eos_lib, only: eos_init, alloc_eos_handle_using_inlist, eos_ptr, eosDT_get, free_eos_handle, eos_shutdown
   implicit none

   integer, parameter :: species = 20
   character(len=4), parameter :: isotope_names(species) = [character(len=4) :: &
      'h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'na23', 'mg24', 'al27', 'si28', &
      'p31', 's32', 'cl35', 'ar40', 'ca40', 'ti48', 'cr52', 'mn55', 'fe56', 'ni58']
   character(len=2), parameter :: element_names(species) = [character(len=2) :: &
      'H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', &
      'P', 'S', 'Cl', 'Ar', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe', 'Ni']
   integer, pointer :: chem_id(:), net_iso(:)
   type(EoS_General_Info), pointer :: controls
   real(dp) :: xa(species), metal_sum, nan, logrho, logt, rho, temperature
   real(dp) :: results(num_eos_basic_results), density_partials(num_eos_basic_results)
   real(dp) :: temperature_partials(num_eos_basic_results), abundance_partials(num_eos_d_dxa_results, species)
   real(dp) :: fractions(num_eos_frac_results), cv_fields(3)
   integer :: handle, ierr, status, stream_unit, mass_unit, species_index, extra_status
   character(len=4096) :: mesa_directory, inlist, massfile, record, cache_directory, temp_directory, extra
   character(len=2) :: element

   open(newunit=stream_unit, file='/dev/stdout', action='write', iostat=status)
   if (status /= 0) error stop 'cannot open output stream'
   open(unit=output_unit, file='/dev/stderr', action='write', iostat=status)
   if (status /= 0) error stop 'cannot redirect MESA diagnostics'
   if (command_argument_count() /= 3) error stop 'usage: mesa_probe MESA_DIR INLIST MASSFILE'
   call get_command_argument(1, mesa_directory)
   call get_command_argument(2, inlist)
   call get_command_argument(3, massfile)
   call get_environment_variable('MESA_CACHES_DIR', cache_directory, status=status)
   if (status /= 0 .or. len_trim(cache_directory) == 0) error stop 'report-owned MESA_CACHES_DIR required'
   call get_environment_variable('MESA_TEMP_CACHES_DIR', temp_directory, status=status)
   if (status /= 0 .or. len_trim(temp_directory) == 0) error stop 'report-owned MESA_TEMP_CACHES_DIR required'

   call const_init(trim(mesa_directory), ierr)
   if (ierr /= 0) error stop 'const_init failed'
   if (.not. use_mesa_temp_cache) error stop 'temporary cache safety must remain enabled'
   call math_init()
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) error stop 'chem_init failed'
   allocate(chem_id(species), net_iso(num_chem_isos))
   net_iso = 0
   do species_index = 1, species
      chem_id(species_index) = chem_get_iso_id(trim(isotope_names(species_index)))
      if (chem_id(species_index) <= 0 .or. chem_id(species_index) > num_chem_isos) then
         write(error_unit, '(a,a)') 'Unsupported isotope: ', trim(isotope_names(species_index))
         error stop '20-isotope composition unsupported'
      end if
      net_iso(chem_id(species_index)) = species_index
   end do
   open(newunit=mass_unit, file=trim(massfile), status='old', action='read', iostat=status)
   if (status /= 0) error stop 'cannot open metal pattern'
   read(mass_unit, '(a)', iostat=status) record
   if (status /= 0) error stop 'cannot read metal pattern header'
   xa(1) = 0.70_dp
   xa(2) = 0.28_dp
   do species_index = 3, species
      read(mass_unit, *, iostat=status) element, xa(species_index)
      if (status /= 0) error stop 'cannot read metal pattern'
      if (element /= element_names(species_index)) error stop 'incorrect metal order'
   end do
   close(mass_unit)
   if (.not. all(ieee_is_finite(xa)) .or. any(xa < 0)) error stop 'invalid metal pattern'
   metal_sum = sum(xa(3:species))
   if (abs(metal_sum - 1.0_dp) > 1e-12_dp) error stop 'metal pattern must sum to one'
   xa(3:species) = 0.02_dp * xa(3:species) / metal_sum
   if (abs(sum(xa) - 1.0_dp) > 1e-12_dp) error stop 'composition must sum to one'
   write(error_unit, '(a,es25.16e3)') 'metal_pattern_sum_before_normalization=', metal_sum
   do species_index = 1, species
      write(error_unit, '(a,1x,i0,1x,es25.16e3)') trim(isotope_names(species_index)), &
         chem_id(species_index), xa(species_index)
   end do

   call eos_init(' ', .true., ierr)
   if (ierr /= 0) error stop 'eos_init failed'
   handle = alloc_eos_handle_using_inlist(trim(inlist), ierr)
   if (ierr /= 0) error stop 'EOS handle allocation failed'
   call eos_ptr(handle, controls, ierr)
   if (ierr /= 0) error stop 'EOS pointer failed'
   if (controls%use_WDM_AV) error stop 'use_WDM_AV must be false'

   write(stream_unit, '(a)') 'logRho logT ierr frac_HELM frac_OPAL_SCVH frac_FreeEOS frac_PC '// &
      'frac_Skye frac_CMS frac_ideal frac_WDM_AV Cv dCv_dlnRho dCv_dlnT'
   flush(stream_unit)
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   do
      read(input_unit, '(a)', iostat=status) record
      if (status == iostat_end) exit
      if (status /= 0) error stop 'stdin read failed'
      record = adjustl(record)
      if (len_trim(record) == 0) cycle
      if (record(1:1) == '#') cycle
      logrho = nan
      logt = nan
      fractions = nan
      cv_fields = nan
      read(record, *, iostat=status) logrho, logt
      read(record, *, iostat=extra_status) logrho, logt, extra
      ierr = 0
      if (status /= 0 .or. extra_status == 0) then
         ierr = -1001
      else if (.not. ieee_is_finite(logrho) .or. .not. ieee_is_finite(logt)) then
         ierr = -1002
      else if (min(logrho, logt) <= log10(tiny(1.0_dp)) + 1 .or. &
               max(logrho, logt) >= log10(huge(1.0_dp)) - 1) then
         ierr = -1003
      else
         rho = 10.0_dp**logrho
         temperature = 10.0_dp**logt
         results = nan
         density_partials = nan
         temperature_partials = nan
         abundance_partials = nan
         call eosDT_get(handle, species, chem_id, net_iso, xa, rho, logrho, temperature, logt, &
            results, density_partials, temperature_partials, abundance_partials, ierr)
         if (ierr == 0) then
            fractions = results(i_frac:i_frac + num_eos_frac_results - 1)
            cv_fields = [results(i_Cv), density_partials(i_Cv), temperature_partials(i_Cv)]
            if (.not. all(ieee_is_finite(results)) .or. .not. all(ieee_is_finite(cv_fields))) ierr = -1004
            if (ierr == 0) then
               if (any(fractions < -1e-10_dp) .or. any(fractions > 1 + 1e-10_dp) .or. &
                   abs(sum(fractions) - 1) > 1e-10_dp) ierr = -1005
            end if
         end if
      end if
      if (ierr /= 0) then
         fractions = nan
         cv_fields = nan
      end if
      write(stream_unit, '(2(es25.16e3,1x),i0,1x,11(es25.16e3,1x))') &
         logrho, logt, ierr, fractions, cv_fields
      flush(stream_unit)
   end do
   call free_eos_handle(handle)
   call eos_shutdown()
   deallocate(chem_id, net_iso)
   close(stream_unit)
end program mesa_probe
