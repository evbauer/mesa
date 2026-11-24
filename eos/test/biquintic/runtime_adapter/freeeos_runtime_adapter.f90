module freeeos_runtime_adapter
   use const_def, only: dp, ln10, crad
   use chem_def, only: num_chem_isos, chem_isos, chem_has_been_initialized
   use chem_lib, only: chem_get_iso_id, basic_composition_info
   use eos_def, only: nv, i_lnPgas, i_lnE, i_lnS, i_mu, i_eta, i_grad_ad, i_chiRho, i_chiT, &
      i_Cp, i_Cv, i_dE_dRho, i_dS_dT, i_dS_dRho, i_gamma1, i_gamma3, i_phase, i_frac_FreeEOS
   use freeeos_table
   use runtime_sha256, only: file_sha256
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: runtime_adapter, runtime_identity, runtime_initialize, runtime_query, runtime_map_partials
   public :: runtime_ok, runtime_uninitialized, runtime_identity_error, runtime_composition_error
   public :: runtime_input_error, runtime_legacy_error, runtime_policy_error, runtime_get_identity
   public :: runtime_species, runtime_composition_tolerance, runtime_hybrid_policy
   integer, parameter :: runtime_ok=0, runtime_uninitialized=100, runtime_identity_error=101
   integer, parameter :: runtime_composition_error=102, runtime_input_error=103, runtime_legacy_error=104
   integer, parameter :: runtime_policy_error=105, runtime_species=20
   real(dp), parameter :: runtime_composition_tolerance=5.0e-13_dp
   character(len=*), parameter :: runtime_hybrid_policy='legacy_aux_and_all27_abundance_rows_no_composition_consistency'
   character(len=64), parameter :: audited_pattern_sha= &
      'bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6'
   character(len=4), parameter :: isotope_names(runtime_species)=[character(len=4) :: &
      'h1','he4','c12','n14','o16','ne20','na23','mg24','al27','si28', &
      'p31','s32','cl35','ar40','ca40','ti48','cr52','mn55','fe56','ni58']
   character(len=2), parameter :: element_names(runtime_species)=[character(len=2) :: &
      'H','He','C','N','O','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','Ca','Ti','Cr','Mn','Fe','Ni']
   type :: runtime_identity
      type(table_metadata) :: table
      character(len=64) :: cache_sha256=''
      character(len=64) :: pattern_sha256=audited_pattern_sha
      real(dp) :: minimum_logT=3.2_dp
   end type runtime_identity
   type :: runtime_adapter
      private
      type(biquintic_table) :: table
      type(runtime_identity) :: identity
      integer :: ids(runtime_species)=0, mapping_size=0
      real(dp) :: masses(runtime_species)=0, summaries(4)=0
      real(dp) :: mass_numbers(runtime_species)=0, charges(runtime_species)=0
      logical :: initialized=.false.
   end type runtime_adapter
contains
   subroutine runtime_initialize(adapter, cache_path, identity, pattern_path, ierr)
      type(runtime_adapter), intent(out) :: adapter
      character(len=*), intent(in) :: cache_path, pattern_path
      type(runtime_identity), intent(in) :: identity
      integer, intent(out) :: ierr
      type(runtime_adapter) :: candidate
      type(table_metadata) :: actual
      character(len=64) :: digest
      integer :: hash_status
      ierr=runtime_policy_error
      if (nv /= 27) return
      if (.not. identity_policy_valid(identity)) return
      ierr=runtime_identity_error
      call file_sha256(cache_path,digest,hash_status)
      if (hash_status /= 0 .or. digest /= identity%cache_sha256) return
      call load_table(cache_path,0.7_dp,0.02_dp,candidate%table,ierr,identity%table%source_sha256)
      if (ierr /= table_ok) return
      actual=get_metadata(candidate%table)
      ierr=runtime_identity_error
      if (.not. metadata_equal(actual,identity%table)) return
      call file_sha256(cache_path,digest,hash_status)
      if (hash_status /= 0 .or. digest /= identity%cache_sha256) return
      call load_pattern(candidate,pattern_path,identity%pattern_sha256,ierr)
      if (ierr /= runtime_ok) return
      candidate%identity=identity
      candidate%initialized=.true.
      adapter=candidate
      ierr=runtime_ok
   end subroutine runtime_initialize

   logical function identity_policy_valid(identity) result(valid)
      type(runtime_identity), intent(in) :: identity
      valid=.false.
      if (.not. valid_sha(identity%cache_sha256) .or. .not. valid_sha(identity%table%source_sha256) .or. &
          .not. valid_sha(identity%table%preprocessing_sha256)) return
      if (identity%pattern_sha256 /= audited_pattern_sha) return
      if (.not. all(ieee_is_finite([identity%minimum_logT,identity%table%crad,identity%table%hydrogen, &
          identity%table%metals,identity%table%u_min,identity%table%du,identity%table%tau_min,identity%table%dtau]))) return
      if (identity%minimum_logT /= 3.2_dp) return
      if (identity%table%split /= 1 .or. identity%table%radius /= 2 .or. identity%table%route /= 'mean') return
      if (abs(identity%table%crad-crad) > 32*epsilon(1.0_dp)*crad) return
      if (abs(identity%table%hydrogen-0.7_dp) > 32*epsilon(1.0_dp)) return
      if (abs(identity%table%metals-0.02_dp) > 32*epsilon(1.0_dp)) return
      if (min(identity%table%nq,identity%table%nt) < 2) return
      if (min(identity%table%du,identity%table%dtau) <= 0) return
      valid=.true.
   end function identity_policy_valid

   logical function valid_sha(digest)
      character(len=*), intent(in) :: digest
      valid_sha=len_trim(digest) == 64 .and. verify(trim(digest),'0123456789abcdef') == 0
   end function valid_sha

   logical function metadata_equal(actual, expected) result(equal)
      type(table_metadata), intent(in) :: actual, expected
      equal=actual%nq == expected%nq .and. actual%nt == expected%nt .and. &
         actual%split == expected%split .and. actual%radius == expected%radius .and. actual%route == expected%route .and. &
         actual%source_sha256 == expected%source_sha256 .and. actual%preprocessing_sha256 == expected%preprocessing_sha256
      equal=equal .and. all([actual%hydrogen,actual%metals,actual%crad,actual%u_min,actual%du,actual%tau_min,actual%dtau] == &
         [expected%hydrogen,expected%metals,expected%crad,expected%u_min,expected%du,expected%tau_min,expected%dtau])
   end function metadata_equal

   subroutine load_pattern(adapter, path, expected_sha, ierr)
      type(runtime_adapter), intent(inout) :: adapter
      character(len=*), intent(in) :: path, expected_sha
      integer, intent(out) :: ierr
      character(len=64) :: digest
      character(len=128) :: record
      character(len=2) :: element
      integer :: unit, io_status, position
      real(dp) :: metal_sum, helium, z2bar, z53bar, ye, correction, total
      ierr=runtime_identity_error
      call file_sha256(path,digest,io_status)
      if (io_status /= 0 .or. digest /= expected_sha) return
      ierr=runtime_composition_error
      if (.not. chem_has_been_initialized .or. num_chem_isos < runtime_species) return
      adapter%mapping_size=num_chem_isos
      do position=1,runtime_species
         adapter%ids(position)=chem_get_iso_id(trim(isotope_names(position)))
         if (adapter%ids(position) <= 0 .or. adapter%ids(position) > num_chem_isos) return
         adapter%mass_numbers(position)=real(chem_isos%Z_plus_N(adapter%ids(position)),dp)
         adapter%charges(position)=real(chem_isos%Z(adapter%ids(position)),dp)
      end do
      open(newunit=unit,file=path,status='old',action='read',iostat=io_status)
      if (io_status /= 0) return
      read(unit,'(a)',iostat=io_status) record
      if (io_status /= 0) then
         close(unit)
         return
      end if
      adapter%masses(1:2)=[0.7_dp,0.28_dp]
      do position=3,runtime_species
         read(unit,*,iostat=io_status) element,adapter%masses(position)
         if (io_status /= 0) exit
         if (element /= element_names(position)) exit
      end do
      close(unit)
      if (position <= runtime_species) return
      if (.not. all(ieee_is_finite(adapter%masses)) .or. any(adapter%masses <= 0)) return
      metal_sum=sum(adapter%masses(3:))
      if (abs(metal_sum-1.0_dp) > 1e-12_dp) return
      adapter%masses(3:)=0.02_dp*adapter%masses(3:)/metal_sum
      call basic_composition_info(runtime_species,adapter%ids,adapter%masses,adapter%summaries(1),helium, &
         adapter%summaries(2),adapter%summaries(3),adapter%summaries(4),z2bar,z53bar,ye,correction,total)
      if (.not. all(ieee_is_finite(adapter%summaries))) return
      if (abs(total-1.0_dp) > runtime_composition_tolerance) return
      call file_sha256(path,digest,io_status)
      if (io_status /= 0 .or. digest /= expected_sha) then
         ierr=runtime_identity_error
         return
      end if
      ierr=runtime_ok
   end subroutine load_pattern

   subroutine runtime_get_identity(adapter, identity, ierr)
      type(runtime_adapter), intent(in) :: adapter
      type(runtime_identity), intent(out) :: identity
      integer, intent(out) :: ierr
      identity=runtime_identity()
      ierr=runtime_uninitialized
      if (.not. adapter%initialized) return
      identity=adapter%identity
      ierr=runtime_ok
   end subroutine runtime_get_identity

   logical function composition_matches(adapter, species, chem_id, net_iso, xa, hydrogen, metals, abar, zbar) result(matches)
      type(runtime_adapter), intent(in) :: adapter
      integer, intent(in) :: species, chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:), hydrogen, metals, abar, zbar
      integer :: position, canonical, mapping_index
      real(dp) :: ordered(runtime_species), actual_summary(4), number_fractions(runtime_species)
      matches=.false.
      if (species /= runtime_species) return
      if (size(chem_id) /= species .or. size(xa) /= species .or. size(net_iso) /= adapter%mapping_size) return
      if (.not. all(ieee_is_finite(xa)) .or. .not. all(ieee_is_finite([hydrogen,metals,abar,zbar]))) return
      if (any(xa <= 0) .or. abs(sum(xa)-1.0_dp) > runtime_composition_tolerance) return
      if (any(chem_id <= 0) .or. any(chem_id > adapter%mapping_size)) return
      if (any(net_iso < 0) .or. any(net_iso > species)) return
      if (count(net_iso /= 0) /= species) return
      do position=1,species
         if (net_iso(chem_id(position)) /= position) return
      end do
      do canonical=1,runtime_species
         mapping_index=net_iso(adapter%ids(canonical))
         if (mapping_index <= 0) return
         if (chem_id(mapping_index) /= adapter%ids(canonical)) return
         ordered(canonical)=xa(mapping_index)
      end do
      if (any(abs(ordered-adapter%masses) > runtime_composition_tolerance*adapter%masses)) return
      number_fractions=ordered/adapter%mass_numbers
      actual_summary(1)=ordered(1)
      actual_summary(2)=1-(ordered(1)+ordered(2))
      actual_summary(3)=sum(ordered)/sum(number_fractions)
      actual_summary(4)=sum(number_fractions*adapter%charges)*actual_summary(3)
      if (.not. all(ieee_is_finite(actual_summary))) return
      if (any(abs(actual_summary-[hydrogen,metals,abar,zbar]) > &
          runtime_composition_tolerance*max(abs(actual_summary),tiny(1.0_dp)))) return
      matches=.true.
   end function composition_matches

   subroutine runtime_query(adapter, species, chem_id, net_iso, xa, X, Z, abar, zbar, &
         rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
      type(runtime_adapter), intent(in) :: adapter
      integer, intent(in) :: species, chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:), X, Z, abar, zbar, rho, logRho, T, logT
      real(dp), intent(inout) :: res(:), d_dlnd(:), d_dlnT(:), d_dxa(:,:)
      integer, intent(out) :: ierr
      type(potential_partials) :: partials
      real(dp) :: lnrho, lnT, native_u, high_u, high_t
      ierr=runtime_uninitialized
      if (.not. adapter%initialized) return
      ierr=runtime_composition_error
      if (.not. composition_matches(adapter,species,chem_id,net_iso,xa,X,Z,abar,zbar)) return
      ierr=runtime_input_error
      if (.not. all(ieee_is_finite([rho,logRho,T,logT]))) return
      if (rho <= 0 .or. T <= 0) return
      lnrho=log(rho)
      lnT=log(T)
      if (abs(lnrho/ln10-logRho) > 64*epsilon(1.0_dp)*max(1.0_dp,abs(logRho))) return
      if (abs(lnT/ln10-logT) > 64*epsilon(1.0_dp)*max(1.0_dp,abs(logT))) return
      ierr=table_out_of_bounds
      if (logT < adapter%identity%minimum_logT) return
      native_u=lnrho-2*lnT+12*ln10
      high_u=adapter%identity%table%u_min+(adapter%identity%table%nq-1)*adapter%identity%table%du
      high_t=adapter%identity%table%tau_min+(adapter%identity%table%nt-1)*adapter%identity%table%dtau
      if (native_u < adapter%identity%table%u_min .or. native_u > high_u) return
      if (lnT < adapter%identity%table%tau_min .or. lnT > high_t) return
      call evaluate_table_real(adapter%table,lnrho,lnT,partials,ierr)
      if (ierr /= table_ok) return
      call runtime_map_partials(partials,rho,T,adapter%identity%table%crad,res,d_dlnd,d_dlnT,d_dxa,ierr)
   end subroutine runtime_query

   subroutine runtime_map_partials(partials, rho, temperature, radiation_constant, res, d_dlnd, d_dlnT, d_dxa, ierr)
      type(potential_partials), intent(in) :: partials
      real(dp), intent(in) :: rho, temperature, radiation_constant
      real(dp), intent(inout) :: res(:), d_dlnd(:), d_dlnT(:), d_dxa(:,:)
      integer, intent(out) :: ierr
      type(thermodynamic_state) :: state
      real(dp) :: values(nv), density_partials(nv), temperature_partials(nv)
      integer, parameter :: destination(7)=[i_grad_ad,i_chiRho,i_chiT,i_Cp,i_Cv,i_gamma1,i_gamma3]
      integer, parameter :: origin(7)=[state_grad_ad,state_chirho,state_chit,state_cp,state_cv,state_gamma1,state_gamma3]
      integer, parameter :: log_origin(3)=[state_pgas,state_energy,state_entropy]
      integer :: position, slot, source
      ierr=runtime_input_error
      if (size(res) /= nv .or. size(d_dlnd) /= nv .or. size(d_dlnT) /= nv) return
      if (size(d_dxa,1) /= nv .or. size(d_dxa,2) /= runtime_species) return
      if (.not. all(ieee_is_finite([rho,temperature,radiation_constant]))) return
      if (rho <= 0 .or. temperature <= 0 .or. radiation_constant < 0) return
      ierr=runtime_legacy_error
      if (.not. all(ieee_is_finite(res(i_mu:i_eta)))) return
      if (.not. all(ieee_is_finite(d_dlnd(i_mu:i_eta)))) return
      if (.not. all(ieee_is_finite(d_dlnT(i_mu:i_eta)))) return
      if (.not. all(ieee_is_finite(d_dxa))) return
      call derive_state(partials,log(rho),log(temperature),radiation_constant,state,ierr)
      if (ierr /= table_ok) return
      values=0
      density_partials=0
      temperature_partials=0
      values(i_mu:i_eta)=res(i_mu:i_eta)
      density_partials(i_mu:i_eta)=d_dlnd(i_mu:i_eta)
      temperature_partials(i_mu:i_eta)=d_dlnT(i_mu:i_eta)
      do slot=i_lnPgas,i_lnS
         source=log_origin(slot)
         values(slot)=log(state%value(source))
         density_partials(slot)=state%dlnrho(source)/state%value(source)
         temperature_partials(slot)=state%dlnT(source)/state%value(source)
      end do
      do position=1,size(destination)
         slot=destination(position)
         source=origin(position)
         values(slot)=state%value(source)
         density_partials(slot)=state%dlnrho(source)
         temperature_partials(slot)=state%dlnT(source)
      end do
      values(i_dE_dRho)=(partials%f_x-partials%f_xt)/rho
      density_partials(i_dE_dRho)=(partials%f_xx-partials%f_xxt-partials%f_x+partials%f_xt)/rho
      temperature_partials(i_dE_dRho)=(partials%f_xt-partials%f_xtt)/rho
      values(i_dS_dT)=state%value(state_cv)/temperature
      density_partials(i_dS_dT)=state%dlnrho(state_cv)/temperature
      temperature_partials(i_dS_dT)=(state%dlnT(state_cv)-state%value(state_cv))/temperature
      values(i_dS_dRho)=(-partials%f_xt/rho)/temperature
      density_partials(i_dS_dRho)=((partials%f_xt-partials%f_xxt)/rho)/temperature
      temperature_partials(i_dS_dRho)=((partials%f_xt-partials%f_xtt)/rho)/temperature
      values(i_phase:nv)=0
      values(i_frac_FreeEOS)=1
      ierr=table_arithmetic_failure
      if (.not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(density_partials)) .or. &
          .not. all(ieee_is_finite(temperature_partials))) return
      res=values
      d_dlnd=density_partials
      d_dlnT=temperature_partials
      ierr=runtime_ok
   end subroutine runtime_map_partials
end module freeeos_runtime_adapter
