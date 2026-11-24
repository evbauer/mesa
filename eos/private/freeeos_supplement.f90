module freeeos_supplement
   use const_def, only: dp
   use forum_m, only: hdf5io_t, OPEN_FILE_RO
   use freeeos_table, only: table_metadata, biquintic_table, assemble_table, &
      table_ok, table_io_error, table_format_error, table_mismatch
   use freeeos_legacy_nodes, only: reconstruct_legacy_nodes
   use eos_def, only: EosDT_XZ_Info
   use runtime_sha256, only: file_sha256
   use iso_fortran_env, only: int64
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: load_legacy_supplement, read_supplement
   character(len=64), parameter :: policy_sha='40b513e48b512e44ea2c020cd29c000b60128661260f38c0c0b7365c2ae10a7d'
   real(dp), parameter :: policy_crad=7.565733250280004e-15_dp
contains
   subroutine load_legacy_supplement(ep, supplement_path, x, z, family_id, source_sha, &
         supplement_sha, table, status)
      type(EosDT_XZ_Info), pointer, intent(in) :: ep
      character(len=*), intent(in) :: supplement_path, source_sha, supplement_sha
      real(dp), intent(in) :: x,z
      integer, intent(in) :: family_id
      type(biquintic_table), intent(out) :: table
      integer, intent(out) :: status
      type(table_metadata) :: metadata
      real(dp), allocatable :: mixed(:,:,:),nodes(:,:,:,:)
      logical, allocatable :: valid(:,:)
      call read_supplement(supplement_path,x,z,family_id,source_sha,supplement_sha,metadata,mixed,valid,status)
      if (status /= table_ok) return
      call reconstruct_legacy_nodes(ep,metadata,valid,nodes,status)
      if (status /= table_ok) return
      ! Preserve the legacy loader's dimensions. A full-source supplement may
      ! cover more Q nodes than that loader reads; use the common affine grid.
      metadata%nq=size(nodes,3)
      nodes(2,1,:,:)=mixed(1,1:metadata%nq,:)
      nodes(1,2,:,:)=mixed(2,1:metadata%nq,:)
      nodes(2,2,:,:)=mixed(3,1:metadata%nq,:)
      call assemble_table(metadata,nodes,valid(1:metadata%nq,:),table,status)
   end subroutine load_legacy_supplement

   subroutine read_supplement(path,x,z,family_id,source_sha,supplement_sha,metadata,mixed,valid,status)
      character(len=*), intent(in) :: path,source_sha,supplement_sha
      real(dp), intent(in) :: x,z
      integer, intent(in) :: family_id
      type(table_metadata), intent(out) :: metadata
      real(dp), allocatable, intent(out) :: mixed(:,:,:)
      logical, allocatable, intent(out) :: valid(:,:)
      integer, intent(out) :: status
      type(hdf5io_t) :: file_handle
      character(len=24), parameter :: int_names(9)=[character(len=24) :: &
         'format_version','source_version','family_id','nu','ntau','radiation_split', &
         'stencil_radius','twist_route','preprocessing_version']
      character(len=24), parameter :: real_names(8)=[character(len=24) :: &
         'X','Z','u_min','du','tau_min','dtau','radiation_constant','minimum_log10_T']
      character(len=24), parameter :: mixed_names(3)=[character(len=24) :: &
         'f_uu_tau','f_u_tau_tau','f_uu_tau_tau']
      character(len=64) :: supplement_digest
      integer :: integer_metadata(9),digest_bytes(32),hash_status,metadata_index,allocation_status
      integer, allocatable :: integer_mask(:,:)
      real(dp) :: real_metadata(8)
      real(dp), allocatable :: derivative_field(:,:)
      logical :: exists
      status=table_io_error
      inquire(file=path,exist=exists)
      if (.not. exists) return
      call file_sha256(path,supplement_digest,hash_status)
      if (hash_status /= 0) return
      status=table_mismatch
      if (supplement_digest /= supplement_sha) return
      file_handle=hdf5io_t(path,OPEN_FILE_RO)
      status=table_format_error
      do metadata_index=1,size(int_names)
         if (.not. scalar_attr(file_handle,trim(int_names(metadata_index)))) goto 900
         call file_handle%read_attr(trim(int_names(metadata_index)),integer_metadata(metadata_index))
      end do
      do metadata_index=1,size(real_names)
         if (.not. scalar_attr(file_handle,trim(real_names(metadata_index)))) goto 900
         call file_handle%read_attr(trim(real_names(metadata_index)),real_metadata(metadata_index))
      end do
      if (.not. all(ieee_is_finite(real_metadata))) goto 900
      if (integer_metadata(1) /= 1 .or. integer_metadata(2) /= 51 .or. integer_metadata(6) /= 1) goto 900
      if (integer_metadata(7) /= 2 .or. integer_metadata(8) /= 1 .or. integer_metadata(9) /= 1) goto 900
      if (min(integer_metadata(4),integer_metadata(5)) < 2 .or. max(integer_metadata(4),integer_metadata(5)) > 100000) goto 900
      if (int(integer_metadata(4),int64)*int(integer_metadata(5),int64) > 2000000) goto 900
      if (real_metadata(4) <= 0 .or. real_metadata(6) <= 0 .or. real_metadata(8) /= 3.2_dp) goto 900
      if (abs(real_metadata(7)-policy_crad) > 32*epsilon(1.0_dp)*policy_crad) goto 900
      status=table_mismatch
      if (integer_metadata(3) /= family_id) goto 900
      if (abs(real_metadata(1)-x) > 32*epsilon(1.0_dp) .or. abs(real_metadata(2)-z) > 32*epsilon(1.0_dp)) goto 900
      status=table_format_error
      if (.not. vector_dset(file_handle,'source_sha256',32)) goto 900
      call file_handle%read_dset('source_sha256',digest_bytes)
      if (any(digest_bytes < 0) .or. any(digest_bytes > 255)) goto 900
      metadata%source_sha256=bytes_hex(digest_bytes)
      if (.not. vector_dset(file_handle,'preprocessing_sha256',32)) goto 900
      call file_handle%read_dset('preprocessing_sha256',digest_bytes)
      if (any(digest_bytes < 0) .or. any(digest_bytes > 255)) goto 900
      metadata%preprocessing_sha256=bytes_hex(digest_bytes)
      status=table_mismatch
      if (metadata%source_sha256 /= source_sha .or. metadata%preprocessing_sha256 /= policy_sha) goto 900
      status=table_format_error
      do metadata_index=1,3
         if (.not. matrix_dset(file_handle,trim(mixed_names(metadata_index)),integer_metadata(4:5))) goto 900
      end do
      if (.not. matrix_dset(file_handle,'valid',integer_metadata(4:5))) goto 900
      allocate(derivative_field(integer_metadata(4),integer_metadata(5)),integer_mask(integer_metadata(4),integer_metadata(5)), &
         mixed(3,integer_metadata(4),integer_metadata(5)),valid(integer_metadata(4),integer_metadata(5)),stat=allocation_status)
      if (allocation_status /= 0) then
         status=table_io_error
         goto 900
      end if
      call file_handle%read_dset('valid',integer_mask)
      if (any(integer_mask /= 0 .and. integer_mask /= 1)) goto 900
      valid=integer_mask == 1
      if (.not. any(valid)) goto 900
      do metadata_index=1,3
         call file_handle%read_dset(trim(mixed_names(metadata_index)),derivative_field)
         if (.not. all(ieee_is_finite(derivative_field))) goto 900
         if (any(derivative_field /= 0 .and. .not. valid)) goto 900
         mixed(metadata_index,:,:)=derivative_field
      end do
      metadata%nq=integer_metadata(4); metadata%nt=integer_metadata(5); metadata%split=1; metadata%radius=2
      metadata%route='mean'; metadata%hydrogen=real_metadata(1); metadata%metals=real_metadata(2)
      metadata%u_min=real_metadata(3); metadata%du=real_metadata(4); metadata%tau_min=real_metadata(5); metadata%dtau=real_metadata(6)
      metadata%crad=real_metadata(7)
      status=table_ok
900   call file_handle%final()
   end subroutine read_supplement

   logical function scalar_attr(file_handle,name)
      type(hdf5io_t), intent(inout) :: file_handle
      character(len=*), intent(in) :: name
      integer(int64), allocatable :: dimensions(:)
      scalar_attr=.false.
      if (.not. file_handle%attr_exists(name)) return
      dimensions=file_handle%attr_shape(name)
      scalar_attr=size(dimensions) == 0
   end function scalar_attr

   logical function vector_dset(file_handle,name,length)
      type(hdf5io_t), intent(inout) :: file_handle
      character(len=*), intent(in) :: name
      integer, intent(in) :: length
      integer(int64), allocatable :: dimensions(:)
      vector_dset=.false.
      if (.not. file_handle%dset_exists(name)) return
      dimensions=file_handle%dset_shape(name)
      if (size(dimensions) /= 1) return
      vector_dset=dimensions(1) == length
   end function vector_dset

   logical function matrix_dset(file_handle,name,expected)
      type(hdf5io_t), intent(inout) :: file_handle
      character(len=*), intent(in) :: name
      integer, intent(in) :: expected(2)
      integer(int64), allocatable :: dimensions(:)
      matrix_dset=.false.
      if (.not. file_handle%dset_exists(name)) return
      dimensions=file_handle%dset_shape(name)
      if (size(dimensions) /= 2) return
      matrix_dset=all(dimensions == expected)
   end function matrix_dset

   function bytes_hex(digest_bytes) result(hex)
      integer, intent(in) :: digest_bytes(32)
      character(len=64) :: hex
      character(len=16), parameter :: digits='0123456789abcdef'
      integer :: metadata_index
      do metadata_index=1,32
         hex(2*metadata_index-1:2*metadata_index-1)=digits(digest_bytes(metadata_index)/16+1:digest_bytes(metadata_index)/16+1)
         hex(2*metadata_index:2*metadata_index)=digits(mod(digest_bytes(metadata_index),16)+1:mod(digest_bytes(metadata_index),16)+1)
      end do
   end function bytes_hex
end module freeeos_supplement
