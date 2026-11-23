program test_co
   use const_def, only: dp
   use eos_def, only: eos_def_init, FreeEOS_XZ_struct, nv, i_lnPgas, i_lnE, i_lnS, i_Cv, i_mu, i_eta
   use chem_def, only: chem_isos
   use freeeos_family, only: freeeos_family_results, freeeos_family_get_counters, &
      family_unavailable, family_manifest_error, family_unsupported_suffix
   use freeeos_table, only: biquintic_table, potential_partials, thermodynamic_state, &
      load_table, evaluate_table_real, derive_state, partial_vector, table_ok, &
      state_pgas, state_energy, state_entropy, state_cv
   use iso_fortran_env, only: int64
   implicit none
   character(len=1024) :: directory, mode
   character(len=30) :: suffixes(15)
   real(dp) :: values(nv), density_derivatives(nv), temperature_derivatives(nv), abundance_derivatives(nv,3)
   real(dp) :: ordinary(nv), selected(nv), lower(nv), upper(nv), lower_abundance(nv,3)
   real(dp), parameter :: rho=1.0e-7_dp, temperature=1.0e5_dp
   real(dp), parameter :: step=1.0e-5_dp
   real(dp) :: metals
   integer :: status, case_index
   integer(int64) :: successes, failures, loads

   call get_command_argument(1,directory)
   call get_command_argument(2,mode)
   if (len_trim(directory) == 0) error stop 'directory'
   call eos_def_init
   allocate(chem_isos%Z(3))
   chem_isos%Z=[1,2,6]
   suffixes=''
   suffixes(13:15)='_CO'
   if (trim(mode) == 'missing' .or. trim(mode) == 'duplicate' .or. trim(mode) == 'unknown' .or. &
       trim(mode) == 'unavailable' .or. trim(mode) == 'absent_payload') then
      call query(0.0_dp,0.8_dp,values,abundance_derivatives,status)
      if (trim(mode) == 'missing' .or. trim(mode) == 'unavailable') &
         call require(status == family_unavailable,'unavailable variant')
      if (trim(mode) == 'duplicate' .or. trim(mode) == 'unknown') &
         call require(status == family_manifest_error,'malformed variant')
      if (trim(mode) == 'absent_payload') call require(status == 10,'absent payload')
      write(*,'(a,1x,a,1x,i0)') 'PASS strict failure',trim(mode),status
      stop
   end if
   if (trim(mode) == 'unsupported') then
      suffixes(13)='_C'
      call query(0.0_dp,0.8_dp,values,abundance_derivatives,status)
      call require(status == family_unsupported_suffix,'unsupported suffix')
      write(*,'(a,i0)') 'PASS unsupported suffix ',status
      stop
   end if
   if (trim(mode) == 'parallel') then
      call reference_state(0.95_dp,.true.,selected,status)
      call require(status == table_ok,'parallel reference')
!$omp parallel do num_threads(4) private(case_index,values,abundance_derivatives,status)
      do case_index=1,64
         call query(0.0_dp,0.95_dp,values,abundance_derivatives,status)
         if (status /= table_ok) error stop 'parallel query'
         call check(values(i_lnPgas),selected(i_lnPgas),2.0e-6_dp,'parallel pressure')
         call check(values(i_Cv),selected(i_Cv),2.0e-6_dp,'parallel Cv')
      end do
!$omp end parallel do
      call freeeos_family_get_counters(successes,failures,loads)
      call require(successes == 64 .and. failures == 0 .and. loads == 2,'parallel counts')
      write(*,'(a,3(1x,i0))') 'PASS parallel',successes,failures,loads
      stop
   end if
   do case_index=1,7
      select case(case_index)
      case(1)
         metals=0.02_dp
      case(2)
         metals=0.75_dp
      case(3)
         metals=0.8_dp
      case(4)
         metals=0.85_dp
      case(5)
         metals=0.9_dp
      case(6)
         metals=0.95_dp
      case(7)
         metals=1.0_dp
      end select
      call query(0.0_dp,metals,values,abundance_derivatives,status)
      call require(status == table_ok,'configured query')
      call reference_state(metals,.true.,selected,status)
      call require(status == table_ok,'configured reference')
      call check(values(i_lnPgas),selected(i_lnPgas),2.0e-6_dp,'configured pressure')
      call check(values(i_lnE),selected(i_lnE),2.0e-6_dp,'configured energy')
      call check(values(i_lnS),selected(i_lnS),2.0e-6_dp,'configured entropy')
      call check(values(i_Cv),selected(i_Cv),2.0e-6_dp,'configured Cv')
      if (metals < 0.75_dp) cycle
      suffixes=''
      call query(0.0_dp,metals,values,lower_abundance,status)
      call require(status == table_ok,'ordinary override')
      call reference_state(metals,.false.,ordinary,status)
      call require(status == table_ok,'ordinary reference')
      call check(values(i_lnPgas),ordinary(i_lnPgas),2.0e-6_dp,'ordinary pressure')
      call require(abs(selected(i_lnPgas)-ordinary(i_lnPgas)) > 1.0e-8_dp,'families differ')
      suffixes(13:15)='_CO'
      if (metals == 0.85_dp .or. metals == 0.95_dp) then
         call query(0.0_dp,metals-step,lower,lower_abundance,status)
         call require(status == table_ok,'lower composition')
         call query(0.0_dp,metals+step,upper,lower_abundance,status)
         call require(status == table_ok,'upper composition')
         call query(0.0_dp,metals,values,abundance_derivatives,status)
         call check(abundance_derivatives(i_lnPgas,3), &
            (upper(i_lnPgas)-lower(i_lnPgas))/(2*step),2.0e-4_dp,'Z derivative')
      end if
   end do
   suffixes(13:15)='_CO'
   call query(0.0_dp,0.95_dp,values,abundance_derivatives,status)
   call require(status == table_ok,'reselected CO')
   call reference_state(0.95_dp,.true.,selected,status)
   call check(values(i_lnPgas),selected(i_lnPgas),2.0e-6_dp,'alternating cache')
   suffixes(13)=''
   call query(0.0_dp,0.85_dp,values,abundance_derivatives,status)
   call require(status == table_ok,'mixed layer configuration')
   call reference_state(0.85_dp,.true.,selected,status)
   call require(status == table_ok,'mixed layer reference')
   call check(values(i_lnPgas),selected(i_lnPgas),2.0e-6_dp,'mixed layer pressure')
   suffixes(13)='_CO'
   call query(0.05_dp,0.85_dp,values,abundance_derivatives,status)
   call require(status == table_ok,'X derivative base')
   call query(0.05_dp-step,0.85_dp,lower,lower_abundance,status)
   call require(status == table_ok,'X derivative lower')
   call query(0.05_dp+step,0.85_dp,upper,lower_abundance,status)
   call require(status == table_ok,'X derivative upper')
   call check(abundance_derivatives(i_lnPgas,1), &
      (upper(i_lnPgas)-lower(i_lnPgas))/(2*step),2.0e-4_dp,'X derivative')
   call freeeos_family_get_counters(successes,failures,loads)
   write(*,'(a,3(1x,i0))') 'PASS configured/ordinary payloads and Z derivatives; counts', &
      successes,failures,loads
contains
   subroutine query(hydrogen,metal_fraction,result,composition,status)
      real(dp), intent(in) :: hydrogen,metal_fraction
      real(dp), intent(out) :: result(nv),composition(nv,3)
      integer, intent(out) :: status
      real(dp) :: derivative_rho(nv),derivative_temperature(nv)
      result=17; derivative_rho=17; derivative_temperature=17; composition=17
      call freeeos_family_results(trim(directory),.true.,1.0e-6_dp,FreeEOS_XZ_struct, &
         hydrogen,metal_fraction,3,[1,2,3],rho,temperature,result,derivative_rho, &
         derivative_temperature,composition,status,suffixes)
      if (status == table_ok) then
         call require(all(result(i_mu:i_eta) == 17),'legacy auxiliary retention')
      end if
   end subroutine query

   subroutine reference_state(metal_fraction,use_co,result,status)
      real(dp), intent(in) :: metal_fraction
      logical, intent(in) :: use_co
      real(dp), intent(out) :: result(nv)
      integer, intent(out) :: status
      type(biquintic_table) :: source
      type(potential_partials) :: potential
      type(thermodynamic_state) :: state
      real(dp) :: combined(10), components(10), weight, lower_z, upper_z
      integer :: lower_layer, upper_layer, layer, table_status
      result=0
      lower_layer=FreeEOS_XZ_struct%nZs
      upper_layer=lower_layer
      do layer=2,FreeEOS_XZ_struct%nZs
         if (metal_fraction < FreeEOS_XZ_struct%Zs(layer)) then
            lower_layer=layer-1
            upper_layer=layer
            exit
         end if
      end do
      lower_z=FreeEOS_XZ_struct%Zs(lower_layer)
      upper_z=FreeEOS_XZ_struct%Zs(upper_layer)
      combined=0
      do layer=lower_layer,upper_layer
         if (layer == lower_layer) then
            if (upper_layer == lower_layer) then
               weight=1
            else
               weight=(upper_z-metal_fraction)/(upper_z-lower_z)
            end if
         else
            weight=(metal_fraction-lower_z)/(upper_z-lower_z)
         end if
         call load_selected(FreeEOS_XZ_struct%Zs(layer),use_co .and. &
            trim(suffixes(layer)) == '_CO',source,table_status)
         if (table_status /= table_ok) then
            status=table_status
            return
         end if
         call evaluate_table_real(source,log(rho),log(temperature),potential,table_status)
         if (table_status /= table_ok) then
            status=table_status
            return
         end if
         components=partial_vector(potential)
         combined=combined+weight*components
      end do
      potential=potential_partials(combined(1),combined(2),combined(3),combined(4),combined(5), &
         combined(6),combined(7),combined(8),combined(9),combined(10))
      call derive_state(potential,log(rho),log(temperature),7.56573325028000399e-15_dp,state,status)
      if (status /= table_ok) return
      result(i_lnPgas)=log(state%value(state_pgas))
      result(i_lnE)=log(state%value(state_energy))
      result(i_lnS)=log(state%value(state_entropy))
      result(i_Cv)=state%value(state_cv)
   end subroutine reference_state

   subroutine load_selected(metal_fraction,use_co,source,status)
      real(dp), intent(in) :: metal_fraction
      logical, intent(in) :: use_co
      type(biquintic_table), intent(out) :: source
      integer, intent(out) :: status
      character(len=256) :: filename, source_sha, payload_sha
      character(len=80) :: magic
      real(dp) :: hydrogen, metals_in_manifest
      integer :: unit, count, position, available, io_status
      status=family_unavailable
      open(newunit=unit,file=trim(directory)//'/family.manifest',status='old',action='read')
      read(unit,*) magic
      read(unit,*) count
      do position=1,count
         read(unit,*,iostat=io_status) hydrogen,metals_in_manifest,available,filename,source_sha,payload_sha
         if (io_status /= 0) exit
         if (abs(hydrogen) > 1.0e-12_dp .or. &
             abs(metals_in_manifest-metal_fraction) > 1.0e-12_dp) cycle
         if (index(trim(filename),'_CO.bq') > 0 .neqv. use_co) cycle
         if (available /= 1) exit
         call load_table(trim(directory)//'/'//trim(filename),hydrogen,metals_in_manifest, &
            source,status,trim(source_sha))
         exit
      end do
      close(unit)
   end subroutine load_selected

   subroutine check(actual,expected,tolerance,label)
      real(dp), intent(in) :: actual,expected,tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) <= tolerance*max(1.0_dp,abs(expected))) return
      write(*,'(a,1x,a,2es23.12)') 'FAIL comparison',label,actual,expected
      error stop 1
   end subroutine check

   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (condition) return
      write(*,'(a,1x,a)') 'FAIL',label
      error stop 1
   end subroutine require
end program test_co
