program runtime_driver
   use const_def, only: dp
   use const_lib, only: const_init
   use math_lib, only: math_init
   use chem_def, only: num_chem_isos
   use chem_lib, only: chem_init, chem_get_iso_id, basic_composition_info
   use eos_def, only: nv
   use freeeos_table, only: potential_partials, table_metadata
   use freeeos_runtime_adapter
   use runtime_sha256, only: file_sha256
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
   use omp_lib, only: omp_get_num_threads
   implicit none
   type(runtime_adapter) :: adapter
   type(runtime_identity) :: identity
   type(potential_partials) :: partials
   character(len=4096) :: mode, cache_path, identity_path, mesa_path, pattern_path, record, argument
   character(len=64) :: digest
   character(len=4), parameter :: names(20)=[character(len=4) :: &
      'h1','he4','c12','n14','o16','ne20','na23','mg24','al27','si28', &
      'p31','s32','cl35','ar40','ca40','ti48','cr52','mn55','fe56','ni58']
   character(len=2) :: element
   integer :: ierr, io_status, unit, position, alteration, failures, iteration, threads, extra_id, query_index
   integer, allocatable :: chem_id(:), net_iso(:), trial_ids(:), trial_net(:)
   real(dp) :: values(nv), density_partials(nv), temperature_partials(nv), abundance(nv,20)
   real(dp) :: reference(nv,3), original_abundance(nv,20), trial_xa(20), xa(20), hydrogen, helium, metals, abar, zbar
   real(dp) :: thread_reference(nv,3,4), thread_logrho(4), thread_logtemp(4), thread_rho(4), thread_temp(4)
   real(dp) :: z2bar, z53bar, ye, correction, total, rho, temperature, logrho, logt, rad, saved
   real(dp) :: vector(10), nan
   logical :: unchanged
   call get_command_argument(1,mode)
   nan=ieee_value(0.0_dp,ieee_quiet_nan)
   if (mode == 'hash') then
      call get_command_argument(2,argument)
      call file_sha256(trim(argument),digest,ierr)
      write(*,'(i0,1x,a)') ierr,digest
      stop
   end if
   if (mode == 'map') then
      do
         read(input_unit,*,iostat=io_status) rho,temperature,rad,vector
         if (io_status /= 0) exit
         partials=potential_partials(vector(1),vector(2),vector(3),vector(4),vector(5),vector(6), &
            vector(7),vector(8),vector(9),vector(10))
         call seed(values,density_partials,temperature_partials,abundance)
         call runtime_map_partials(partials,rho,temperature,rad,values,density_partials,temperature_partials,abundance,ierr)
         call emit(ierr,values,density_partials,temperature_partials,abundance)
      end do
      stop
   end if
   call get_command_argument(2,mesa_path)
   call get_command_argument(3,cache_path)
   call get_command_argument(4,identity_path)
   call get_command_argument(5,pattern_path)
   call const_init(trim(mesa_path),ierr)
   if (ierr /= 0) error stop 'const_init'
   call math_init()
   call chem_init('isotopes.data',ierr)
   if (ierr /= 0) error stop 'chem_init'
   open(newunit=unit,file=trim(identity_path),status='old',action='read')
   read(unit,*) identity%cache_sha256
   read(unit,*) identity%table%nq,identity%table%nt
   read(unit,*) identity%table%hydrogen,identity%table%metals,identity%table%crad,identity%table%split
   read(unit,*) identity%table%source_sha256
   read(unit,*) identity%table%preprocessing_sha256,identity%table%radius,identity%table%route
   read(unit,*) identity%table%u_min,identity%table%du,identity%table%tau_min,identity%table%dtau
   read(unit,*) identity%minimum_logT,identity%pattern_sha256
   close(unit)
   if (mode /= 'uninitialized') then
      call runtime_initialize(adapter,trim(cache_path),identity,trim(pattern_path),ierr)
      write(*,'(a,1x,i0)') 'INIT',ierr
      if (ierr /= 0 .or. mode == 'initialize') stop
   end if
   allocate(chem_id(20),net_iso(num_chem_isos),trial_ids(20),trial_net(num_chem_isos))
   net_iso=0
   do position=1,20
      chem_id(position)=chem_get_iso_id(trim(names(position)))
      net_iso(chem_id(position))=position
   end do
   xa(1:2)=[0.7_dp,0.28_dp]
   open(newunit=unit,file=trim(pattern_path),status='old',action='read')
   read(unit,'(a)') record
   do position=3,20
      read(unit,*) element,xa(position)
   end do
   close(unit)
   xa(3:)=0.02_dp*xa(3:)/sum(xa(3:))
   call basic_composition_info(20,chem_id,xa,hydrogen,helium,metals,abar,zbar,z2bar,z53bar,ye,correction,total)
   if (mode == 'threads') then
      read(input_unit,*) logrho,logt
      thread_logrho=logrho+[-0.07_dp,-0.02_dp,0.02_dp,0.07_dp]
      thread_logtemp=logt+[0.0_dp,0.003_dp,-0.003_dp,0.007_dp]
      thread_rho=10.0_dp**thread_logrho
      thread_temp=10.0_dp**thread_logtemp
      do query_index=1,4
         call seed(values,density_partials,temperature_partials,abundance)
         call runtime_query(adapter,20,chem_id,net_iso,xa,hydrogen,metals,abar,zbar,thread_rho(query_index), &
            thread_logrho(query_index),thread_temp(query_index),thread_logtemp(query_index), &
            values,density_partials,temperature_partials,abundance,ierr)
         if (ierr /= 0) error stop 'serial thread reference'
         thread_reference(:,1,query_index)=values
         thread_reference(:,2,query_index)=density_partials
         thread_reference(:,3,query_index)=temperature_partials
      end do
      original_abundance=abundance
      failures=0
      threads=0
      !$omp parallel default(shared) private(iteration,query_index,values,density_partials,temperature_partials,abundance,ierr) &
      !$omp reduction(+:failures) reduction(max:threads)
      threads=max(threads,omp_get_num_threads())
      !$omp do
      do iteration=1,1000
         query_index=1+mod(iteration-1,4)
         call seed(values,density_partials,temperature_partials,abundance)
         call runtime_query(adapter,20,chem_id,net_iso,xa,hydrogen,metals,abar,zbar,thread_rho(query_index), &
            thread_logrho(query_index),thread_temp(query_index),thread_logtemp(query_index), &
            values,density_partials,temperature_partials,abundance,ierr)
         if (ierr /= 0 .or. any(values /= thread_reference(:,1,query_index)) .or. &
             any(density_partials /= thread_reference(:,2,query_index)) .or. &
             any(temperature_partials /= thread_reference(:,3,query_index)) .or. &
             any(abundance /= original_abundance)) failures=failures+1
      end do
      !$omp end do
      !$omp end parallel
      write(*,'(a,3(1x,i0))') 'THREADS',threads,1000,failures
      if (failures /= 0) stop 1
      stop
   end if
   do
      read(input_unit,*,iostat=io_status) alteration,logrho,logt
      if (io_status /= 0) exit
      rho=10.0_dp**logrho
      temperature=10.0_dp**logt
      trial_ids=chem_id
      trial_net=net_iso
      trial_xa=xa
      call seed(values,density_partials,temperature_partials,abundance)
      call basic_composition_info(20,chem_id,xa,hydrogen,helium,metals,abar,zbar,z2bar,z53bar,ye,correction,total)
      select case(alteration)
      case(1)
         do position=1,20
            trial_ids(position)=chem_id(21-position)
            trial_xa(position)=xa(21-position)
            trial_net(trial_ids(position))=position
         end do
      case(2)
         saved=trial_xa(3)
         trial_xa(3)=trial_xa(4)
         trial_xa(4)=saved
      case(3)
         trial_ids(4)=trial_ids(3)
      case(4)
         trial_net(trial_ids(3))=4
      case(5)
         trial_xa(3)=nan
      case(6)
         trial_xa(3)=-trial_xa(3)
      case(7)
         trial_xa(3)=0
      case(8)
         trial_xa=trial_xa*1.00000001_dp
      case(9)
         hydrogen=hydrogen+1e-8_dp
      case(10)
         metals=metals+1e-8_dp
      case(11)
         abar=abar+1e-8_dp
      case(12)
         zbar=zbar+1e-8_dp
      case(13)
         rho=rho*1.001_dp
      case(14)
         temperature=temperature*1.001_dp
      case(15)
         rho=0
      case(16)
         logt=nan
      case(17)
         abundance(27,20)=nan
      case(18)
         values(4)=nan
      case(19)
         density_partials(5)=nan
      case(20)
         temperature_partials(6)=nan
      case(21)
         trial_ids(1)=num_chem_isos+1
      case(22)
         trial_net=0
      case(23)
         trial_net(1)=-1
      case(28)
         extra_id=chem_get_iso_id('he3')
         trial_net(extra_id)=21
      case(29)
         trial_ids(3)=chem_get_iso_id('c13')
         trial_net(chem_id(3))=0
         trial_net(trial_ids(3))=3
      case(30)
         extra_id=chem_get_iso_id('he3')
         trial_net(extra_id)=1
      case(31)
         hydrogen=nan
      case(32)
         abar=nan
      case(33)
         zbar=nan
      case(34)
         metals=nan
      end select
      reference(:,1)=values
      reference(:,2)=density_partials
      reference(:,3)=temperature_partials
      original_abundance=abundance
      select case(alteration)
      case(24)
         call runtime_query(adapter,19,trial_ids(:19),trial_net,trial_xa(:19),hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      case(25)
         call runtime_query(adapter,20,trial_ids,trial_net(:num_chem_isos-1),trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      case(26)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values(:26),density_partials,temperature_partials,abundance,ierr)
      case(27)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance(:2,:),ierr)
      case(28)
         call runtime_query(adapter,21,[trial_ids,extra_id],trial_net,[trial_xa,0.0_dp],hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      case(35)
         call runtime_query(adapter,20,trial_ids(:19),trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      case(36)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa(:19),hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      case(37)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials(:26),temperature_partials,abundance,ierr)
      case(38)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials(:26),abundance,ierr)
      case(39)
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance(:,:19),ierr)
      case default
         call runtime_query(adapter,20,trial_ids,trial_net,trial_xa,hydrogen,metals,abar,zbar, &
            rho,logrho,temperature,logt,values,density_partials,temperature_partials,abundance,ierr)
      end select
      unchanged=all_same(values,reference(:,1)) .and. all_same(density_partials,reference(:,2)) .and. &
         all_same(temperature_partials,reference(:,3)) .and. &
         all_same(reshape(abundance,[nv*20]),reshape(original_abundance,[nv*20]))
      write(*,'(a,1x,i0,1x,l1)') 'QUERY',alteration,unchanged
      call emit(ierr,values,density_partials,temperature_partials,abundance)
   end do
contains
   subroutine seed(values, density_partials, temperature_partials, abundance)
      real(dp), intent(out) :: values(nv), density_partials(nv), temperature_partials(nv), abundance(nv,20)
      integer :: slot, species_index
      do slot=1,nv
         values(slot)=real(slot,dp)/10
         density_partials(slot)=real(slot,dp)/100
         temperature_partials(slot)=-real(slot,dp)/1000
         do species_index=1,20
            abundance(slot,species_index)=real(100*slot+species_index,dp)
         end do
      end do
   end subroutine seed
   subroutine emit(ierr, values, density_partials, temperature_partials, abundance)
      integer, intent(in) :: ierr
      real(dp), intent(in) :: values(nv), density_partials(nv), temperature_partials(nv), abundance(nv,20)
      write(output_unit,'(i0,1x,*(es25.16e3,1x))') ierr,values,density_partials,temperature_partials,abundance
   end subroutine emit
   logical function all_same(left, right)
      real(dp), intent(in) :: left(:), right(:)
      all_same=all((left == right) .or. (.not. ieee_is_finite(left) .and. .not. ieee_is_finite(right)))
   end function all_same
end program runtime_driver
