module bridge_support
   use const_def, only: dp, use_mesa_temp_cache
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use eos_def
   use eos_lib, only: eos_init, eosDT_get, eosDT_get_component, alloc_eos_handle_using_inlist
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none

   integer, parameter :: pinned_species=20
   character(len=4), parameter :: isotope_names(pinned_species) = [character(len=4) :: &
      'h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'na23', 'mg24', 'al27', 'si28', &
      'p31', 's32', 'cl35', 'ar40', 'ca40', 'ti48', 'cr52', 'mn55', 'fe56', 'ni58']
   character(len=2), parameter :: element_names(pinned_species) = [character(len=2) :: &
      'H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', &
      'P', 'S', 'Cl', 'Ar', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe', 'Ni']
   integer, pointer :: pinned_chem_id(:), pinned_net_iso(:)
   real(dp) :: pinned_xa(pinned_species)

contains

   subroutine initialize_eos(mesa_directory, massfile, ierr)
      character(len=*), intent(in) :: mesa_directory, massfile
      integer, intent(out) :: ierr
      integer :: species_index, mass_unit, env_status
      character(len=4096) :: record, cache_path, temp_path
      character(len=2) :: element
      real(dp) :: metal_sum

      call get_environment_variable('MESA_CACHES_DIR', cache_path, status=env_status)
      if (env_status /= 0 .or. len_trim(cache_path) == 0) error stop 'private MESA cache required'
      call get_environment_variable('MESA_TEMP_CACHES_DIR', temp_path, status=env_status)
      if (env_status /= 0 .or. len_trim(temp_path) == 0) error stop 'private temporary MESA cache required'
      call const_init(trim(mesa_directory), ierr)
      if (ierr /= 0) return
      if (.not. use_mesa_temp_cache) error stop 'MESA temporary cache must remain enabled'
      call math_init()
      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) return
      allocate(pinned_chem_id(pinned_species), pinned_net_iso(num_chem_isos))
      pinned_net_iso = 0
      do species_index=1,pinned_species
         pinned_chem_id(species_index)=chem_get_iso_id(trim(isotope_names(species_index)))
         if (pinned_chem_id(species_index) < 1) error stop 'missing pinned isotope'
         pinned_net_iso(pinned_chem_id(species_index))=species_index
      end do
      open(newunit=mass_unit, file=trim(massfile), status='old', action='read', iostat=ierr)
      if (ierr /= 0) return
      read(mass_unit, '(a)', iostat=ierr) record
      if (ierr /= 0) error stop 'bad metal-pattern header'
      pinned_xa(1)=0.70_dp
      pinned_xa(2)=0.28_dp
      do species_index=3,pinned_species
         read(mass_unit, *, iostat=ierr) element, pinned_xa(species_index)
         if (ierr /= 0) error stop 'bad metal-pattern row'
         if (element /= element_names(species_index)) error stop 'unexpected metal-pattern order'
      end do
      close(mass_unit)
      if (.not. all(ieee_is_finite(pinned_xa)) .or. any(pinned_xa < 0)) error stop 'bad mass fractions'
      metal_sum=sum(pinned_xa(3:))
      if (abs(metal_sum-1.0_dp) > 1e-12_dp) error stop 'metal-pattern sum differs from one'
      pinned_xa(3:)=0.02_dp*pinned_xa(3:)/metal_sum
      if (abs(sum(pinned_xa)-1.0_dp) > 1e-12_dp) error stop 'mixture sum differs from one'
      call eos_init(' ', .true., ierr)
   end subroutine initialize_eos

   subroutine query_eos(route, handle, log_rho, log_temp, res, density_partials, temperature_partials, abundance_partials, ierr)
      character(len=*), intent(in) :: route
      integer, intent(in) :: handle
      real(dp), intent(in) :: log_rho, log_temp
      real(dp), intent(out) :: res(nv), density_partials(nv), temperature_partials(nv)
      real(dp), intent(out) :: abundance_partials(num_eos_d_dxa_results,pinned_species)
      integer, intent(out) :: ierr
      real(dp) :: rho, temperature, nan_value

      nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
      res=nan_value
      density_partials=nan_value
      temperature_partials=nan_value
      abundance_partials=nan_value
      rho=10.0_dp**log_rho
      temperature=10.0_dp**log_temp
      select case(trim(route))
      case('component')
         call eosDT_get_component(handle, i_eos_FreeEOS, pinned_species, pinned_chem_id, pinned_net_iso, pinned_xa, &
            rho, log_rho, temperature, log_temp, res, density_partials, temperature_partials, abundance_partials, ierr)
      case('final')
         call eosDT_get(handle, pinned_species, pinned_chem_id, pinned_net_iso, pinned_xa, &
            rho, log_rho, temperature, log_temp, res, density_partials, temperature_partials, abundance_partials, ierr)
      case default
         ierr=-7801
      end select
   end subroutine query_eos

end module bridge_support
