program mesa_benchmark
   use const_def, only: dp, ln10, use_mesa_temp_cache
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use eos_def
   use eos_lib, only: eos_init, alloc_eos_handle_using_inlist, eos_ptr, eosDT_get, free_eos_handle, eos_shutdown
   use benchmark_support
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: species=20, field_count=3*num_eos_basic_results+num_eos_d_dxa_results*species
   character(len=4), parameter :: isotope_names(species)=[character(len=4) :: &
      'h1','he4','c12','n14','o16','ne20','na23','mg24','al27','si28', &
      'p31','s32','cl35','ar40','ca40','ti48','cr52','mn55','fe56','ni58']
   character(len=2), parameter :: element_names(species)=[character(len=2) :: &
      'H','He','C','N','O','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','Ca','Ti','Cr','Mn','Fe','Ni']
   character(len=16), parameter :: names(1)=[character(len=16) :: 'default_eosDT']
   integer, parameter :: fields(1)=[field_count]
   integer, pointer :: chem_id(:), net_iso(:)
   type(EoS_General_Info), pointer :: controls
   real(dp) :: xa(species), metal_sum, values(field_count), scales(field_count,1), raw_sums(field_count)
   real(dp), allocatable :: linear_density(:), linear_temperature(:), log_density(:), log_temperature(:)
   real(dp) :: process_started, started, target, fraction_sums(num_eos_frac_results)
   integer :: handle, status, mass_unit, species_index, query, repeats, fraction_positive(num_eos_frac_results)
   character(len=4096) :: mesa_directory, inlist, massfile, query_path, record, cache_directory, temp_directory
   character(len=2) :: element

   process_started=wall_time()
   if (command_argument_count() /= 6) &
      error stop 'usage: mesa_benchmark MESA_DIR INLIST MASSFILE QUERIES REPEATS TARGET_SECONDS'
   call open_report()
   call get_command_argument(1,mesa_directory)
   call get_command_argument(2,inlist)
   call get_command_argument(3,massfile)
   call get_command_argument(4,query_path)
   call read_controls(5,6,repeats,target)
   call get_environment_variable('MESA_CACHES_DIR',cache_directory,status=status)
   if (status /= 0 .or. len_trim(cache_directory) == 0) error stop 'experiment-owned MESA_CACHES_DIR required'
   call get_environment_variable('MESA_TEMP_CACHES_DIR',temp_directory,status=status)
   if (status /= 0 .or. len_trim(temp_directory) == 0) error stop 'experiment-owned temporary caches required'
   call read_queries(trim(query_path))
   allocate(linear_density(size(coordinates,2)),linear_temperature(size(coordinates,2)), &
      log_density(size(coordinates,2)),log_temperature(size(coordinates,2)))
   log_density=coordinates(1,:)/ln10
   log_temperature=coordinates(2,:)/ln10
   linear_density=exp(coordinates(1,:))
   linear_temperature=exp(coordinates(2,:))
   if (.not. all(ieee_is_finite(linear_density)) .or. .not. all(ieee_is_finite(linear_temperature))) &
      error stop 'nonfinite linear MESA input'
   call memory_snapshot('before_load')
   started=wall_time()
   call const_init(trim(mesa_directory),status)
   if (status /= 0) error stop 'const_init failed'
   if (.not. use_mesa_temp_cache) error stop 'temporary cache safety must remain enabled'
   call math_init()
   call chem_init('isotopes.data',status)
   if (status /= 0) error stop 'chem_init failed'
   allocate(chem_id(species),net_iso(num_chem_isos))
   net_iso=0
   do species_index=1,species
      chem_id(species_index)=chem_get_iso_id(trim(isotope_names(species_index)))
      if (chem_id(species_index) <= 0 .or. chem_id(species_index) > num_chem_isos) error stop 'unsupported isotope'
      net_iso(chem_id(species_index))=species_index
   end do
   open(newunit=mass_unit,file=trim(massfile),status='old',action='read',iostat=status)
   if (status /= 0) error stop 'cannot open scaled-solar pattern'
   read(mass_unit,'(a)',iostat=status) record
   if (status /= 0) error stop 'cannot read pattern header'
   xa(1)=0.70_dp
   xa(2)=0.28_dp
   do species_index=3,species
      read(mass_unit,*,iostat=status) element,xa(species_index)
      if (status /= 0) error stop 'cannot read metal pattern'
      if (element /= element_names(species_index)) error stop 'incorrect element order'
   end do
   close(mass_unit)
   if (.not. all(ieee_is_finite(xa)) .or. any(xa < 0)) error stop 'invalid abundances'
   metal_sum=sum(xa(3:species))
   if (abs(metal_sum-1.0_dp) > 1e-12_dp) error stop 'metal pattern does not sum to one'
   xa(3:species)=0.02_dp*xa(3:species)/metal_sum
   if (abs(sum(xa)-1.0_dp) > 1e-12_dp) error stop 'composition does not sum to one'
   do species_index=1,species
      write(report_unit,*) 'isotope ',trim(isotope_names(species_index)),chem_id(species_index),xa(species_index)
   end do
   write(report_unit,*) 'phase_seconds',' const_chem_composition ',wall_time()-started
   started=wall_time()
   call eos_init(' ',.true.,status)
   if (status /= 0) error stop 'eos_init failed'
   handle=alloc_eos_handle_using_inlist(trim(inlist),status)
   if (status /= 0) error stop 'EOS allocation failed'
   call eos_ptr(handle,controls,status)
   if (status /= 0) error stop 'eos_ptr failed'
   if (controls%use_WDM_AV) error stop 'WDM_AV must be disabled as in experiment 014'
   write(report_unit,*) 'phase_seconds',' eos_init_handle ',wall_time()-started
   write(report_unit,*) 'mesa_output_counts',num_eos_basic_results,num_eos_d_dxa_results,species,field_count
   call memory_snapshot('after_load')
   scales=1
   raw_sums=0
   fraction_sums=0
   fraction_positive=0
   started=wall_time()
   do query=1,size(coordinates,2)
      call compute(1,query,values,status)
      if (status /= 0) then
         write(report_unit,*) 'rejected_query',query,status
         error stop 'default MESA rejected query from identical benchmark corpus'
      end if
      if (any(values(i_frac:i_frac+num_eos_frac_results-1) < -1e-10_dp) .or. &
          abs(sum(values(i_frac:i_frac+num_eos_frac_results-1))-1.0_dp) > 1e-10_dp) &
         error stop 'invalid MESA blend fractions'
      fraction_sums=fraction_sums+values(i_frac:i_frac+num_eos_frac_results-1)
      fraction_positive=fraction_positive+merge(1,0,values(i_frac:i_frac+num_eos_frac_results-1) > 0)
      scales(:,1)=max(scales(:,1),abs(values))
      raw_sums=raw_sums+values
      if (query == 1) write(report_unit,*) 'first_mesa',values
   end do
   write(report_unit,*) 'phase_seconds',' validation ',wall_time()-started
   write(report_unit,*) 'fraction_means',fraction_sums/size(coordinates,2)
   write(report_unit,*) 'fraction_positive_counts',fraction_positive
   write(report_unit,*) 'raw_sum_mesa',raw_sums
   write(report_unit,*) 'normalization ',trim(names(1)),scales(:,1)
   write(report_unit,*) 'phase_seconds',' startup_to_validated ',wall_time()-process_started
   call benchmark_modes(compute,names,fields,scales,repeats,target)
   call memory_snapshot('finished')
   call free_eos_handle(handle)
   call eos_shutdown()
   write(report_unit,*) 'phase_seconds',' total ',wall_time()-process_started
contains
   subroutine compute(mode, query, values, status)
      integer, intent(in) :: mode, query
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: status
      real(dp) :: results(num_eos_basic_results), density_partials(num_eos_basic_results)
      real(dp) :: temperature_partials(num_eos_basic_results), abundance_partials(num_eos_d_dxa_results,species)
      call eosDT_get(handle,species,chem_id,net_iso,xa,linear_density(query),log_density(query), &
         linear_temperature(query),log_temperature(query),results,density_partials,temperature_partials,abundance_partials,status)
      values=[results,density_partials,temperature_partials,reshape(abundance_partials,[num_eos_d_dxa_results*species])]
      if (status == 0 .and. .not. all(ieee_is_finite(values))) status=-1004
   end subroutine compute
end program mesa_benchmark
