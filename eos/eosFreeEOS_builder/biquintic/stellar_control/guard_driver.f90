program guard_driver
   use iso_fortran_env, only: output_unit
   use const_def, only: dp
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id
   use eos_def, only: nv, num_eos_d_dxa_results, EoS_General_Info, get_eos_ptr
   use eos_lib, only: eos_init, eosDT_get, alloc_eos_handle
   use run_star_extras, only: pin_count, pin_ids, pin_xa, query_unit, query_phase, guarded_legacy_results
   implicit none
   character(len=4096) :: mesa_directory, pattern_file, original_pattern, header
   character(len=32) :: mode, isotope
   integer :: ierr, position, pin_unit, species, handle
   integer, pointer :: ids(:), mapping(:)
   real(dp), allocatable :: abundances(:), abundance_partials(:,:)
   real(dp) :: res(nv), density_partials(nv), temperature_partials(nv), saved_mass
   real(dp) :: reference_masses(pin_count), metal_sum, relative_error
   type(EoS_General_Info), pointer :: request
   call get_environment_variable('MESA_DIR',mesa_directory)
   call get_command_argument(1,mode)
   call get_command_argument(2,pattern_file)
   call get_command_argument(3,original_pattern)
   call const_init(trim(mesa_directory),ierr)
   if (ierr /= 0) stop 2
   call math_init()
   call chem_init('isotopes.data',ierr)
   if (ierr /= 0) stop 2
   call eos_init(' ',.true.,ierr)
   if (ierr /= 0) stop 2
   open(newunit=pin_unit,file=trim(pattern_file),status='old')
   do position=1,pin_count
      read(pin_unit,*) isotope,pin_xa(position)
      pin_ids(position)=chem_get_iso_id(trim(isotope))
   end do
   close(pin_unit)
   open(newunit=pin_unit,file=trim(original_pattern),status='old')
   read(pin_unit,'(a)') header
   reference_masses(1:2)=[0.7_dp,0.28_dp]
   do position=3,pin_count
      read(pin_unit,*) isotope,reference_masses(position)
   end do
   close(pin_unit)
   metal_sum=sum(reference_masses(3:))
   reference_masses(3:)=0.02_dp*reference_masses(3:)/metal_sum
   relative_error=maxval(abs(pin_xa-reference_masses)/reference_masses)
   write(output_unit,'(a,es25.16)') 'TASK049_REFERENCE_RELATIVE_ERROR ',relative_error
   flush(output_unit)
   if (relative_error > 5d-13) error stop 'TASK049 normalization differs from 047 reference'
   species=pin_count
   if (trim(mode) == 'extra_zero') species=pin_count+1
   allocate(ids(species),mapping(num_chem_isos),abundances(species),abundance_partials(num_eos_d_dxa_results,species))
   ids(:pin_count)=pin_ids
   abundances(:pin_count)=pin_xa
   query_phase='guard_current'
   select case(trim(mode))
   case('ordered')
   case('reordered')
      ids(:pin_count)=pin_ids(pin_count:1:-1)
      abundances(:pin_count)=pin_xa(pin_count:1:-1)
   case('metal_swap','previous_swap','retry_swap')
      saved_mass=abundances(3)
      abundances(3)=abundances(5)
      abundances(5)=saved_mass
      if (trim(mode) == 'previous_swap') query_phase='previous_replay'
      if (trim(mode) == 'retry_swap') query_phase='trial_or_retry'
   case('normalization')
      abundances(1)=abundances(1)+1d-8
   case('extra_zero')
      ids(species)=chem_get_iso_id('he3')
      abundances(species)=0
   case default
      error stop 'unknown TASK049 guard case'
   end select
   mapping=0
   do position=1,species
      mapping(ids(position))=position
   end do
   handle=alloc_eos_handle(ierr)
   if (ierr /= 0) stop 2
   call get_eos_ptr(handle,request,ierr)
   if (ierr /= 0) stop 2
   request% use_WDM_AV=.false.
   request% okay_to_convert_ierr_to_skip=.false.
   request% use_other_eos_results=.true.
   request% other_eos_results=>guarded_legacy_results
   open(newunit=query_unit,file='guard_queries.tsv',status='new')
   call eosDT_get(handle,species,ids,mapping,abundances,10d0**(-4.2_dp),-4.2_dp,10d0**4.54_dp,4.54_dp, &
      res,density_partials,temperature_partials,abundance_partials,ierr)
   if (ierr /= 0) error stop 'EOS rejected guard query before final callback'
   close(query_unit)
   write(*,*) 'TASK049_GUARD_ACCEPTED ',trim(mode)
end program guard_driver
