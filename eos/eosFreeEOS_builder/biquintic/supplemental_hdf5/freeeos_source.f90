! Historical raw-source comparison utility; omitted from the production EOS build.
module freeeos_source
   use const_def, only: dp, ln10
   use freeeos_table, only: table_metadata, read_tokens, table_ok, table_io_error, table_format_error, table_mismatch
   use runtime_sha256, only: file_sha256
   use iso_fortran_env, only: iostat_end
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: reconstruct_source, cross_support
contains
   subroutine reconstruct_source(path, metadata, nodes, complete, status)
      character(len=*), intent(in) :: path
      type(table_metadata), intent(in) :: metadata
      real(dp), allocatable, intent(out) :: nodes(:,:,:,:)
      logical, allocatable, intent(out) :: complete(:,:)
      integer, intent(out) :: status
      character(len=128) :: source_tokens(22)
      character(len=64) :: source_digest
      real(dp) :: source_header(11), source_row(22), logq_value, logt_value, first_logq, first_logt, delta_logq, delta_logt, previous_logq
      real(dp) :: temperature, density, gas_pressure, radiation_pressure, pressure, energy, entropy, radiation_scale
      real(dp) :: potential_value, density_first, temperature_first, density_second, mixed_second, temperature_second, delta_density
      real(dp), allocatable :: logq_axis(:), logt_axis(:)
      logical, allocatable :: source_valid(:,:), first_support(:,:)
      logical :: expect_label, accepted, anomaly_qmax, anomaly_nq
      integer :: source_unit, io_status, token_count, token_status, index_q, index_temperature, field_index
      status=table_io_error
      call file_sha256(path,source_digest,io_status)
      if (io_status /= 0) return
      status=table_mismatch
      if (source_digest /= metadata%source_sha256) return
      status=table_io_error
      open(newunit=source_unit,file=path,status='old',action='read',iostat=io_status)
      if (io_status /= 0) return
      status=table_format_error
      call read_tokens(source_unit,source_tokens,token_count,token_status)
      if (token_status /= 0) goto 900
      call read_tokens(source_unit,source_tokens,token_count,token_status)
      if (token_status /= 0 .or. token_count /= 11) goto 900
      do field_index=1,11
         read(source_tokens(field_index),*,iostat=io_status) source_header(field_index)
         if (io_status /= 0) goto 900
      end do
      if (.not. all(ieee_is_finite(source_header))) goto 900
      if (source_header(1) /= 51 .or. source_header(7) <= 0 .or. source_header(11) <= 0) goto 900
      if (source_header(4) /= metadata%nt) goto 900
      anomaly_qmax=ends_with(path,'mesa-FreeEOS_02z60x.data')
      anomaly_nq=ends_with(path,'mesa-FreeEOS_04z50x.data')
      if (source_header(8) /= metadata%nq .and. .not. anomaly_nq) goto 900
      if (anomaly_nq .and. source_header(8) /= 445) goto 900
      if (anomaly_qmax .and. source_header(10) /= 3.69_dp) goto 900
      if (anomaly_nq .and. source_header(10) /= 3.23_dp) goto 900
      status=table_mismatch
      if (.not. close(source_header(2),metadata%hydrogen) .or. .not. close(source_header(3),metadata%metals)) goto 900
      status=table_format_error
      allocate(nodes(0:2,0:2,metadata%nq,metadata%nt),source_valid(metadata%nq,metadata%nt), &
         logq_axis(metadata%nq),logt_axis(metadata%nt),stat=io_status)
      if (io_status /= 0) then
         status=table_io_error
         goto 900
      end if
      nodes=0; source_valid=.false.; logq_axis=0; logt_axis=0
      index_q=0; index_temperature=0; expect_label=.false.; first_logq=0; first_logt=0; delta_logq=0; delta_logt=0; previous_logq=0
      do
         call read_tokens(source_unit,source_tokens,token_count,token_status)
         if (token_status == iostat_end) exit
         if (token_status /= 0) goto 900
         if (token_count == 0) cycle
         if (source_tokens(1) == 'logQ' .and. token_count == 7) then
            if (expect_label) goto 900
            if (index_q > 0 .and. index_temperature /= metadata%nt) goto 900
            expect_label=.true.
            cycle
         end if
         if (expect_label) then
            if (token_count /= 1 .or. index_q >= metadata%nq) goto 900
            read(source_tokens(1),*,iostat=io_status) logq_value
            if (io_status /= 0 .or. .not. ieee_is_finite(logq_value)) goto 900
            index_q=index_q+1; index_temperature=0; expect_label=.false.; logq_axis(index_q)=logq_value
            if (index_q == 1) first_logq=logq_value
            if (index_q == 2) delta_logq=logq_value-previous_logq
            if (index_q >= 2) then
               if (delta_logq <= 0 .or. abs(logq_value-previous_logq-delta_logq) > 1e-10_dp) goto 900
            end if
            previous_logq=logq_value
            cycle
         end if
         if (scan(source_tokens(1)(1:1),'+-0123456789') == 0) cycle
         if (token_count /= 22 .or. index_q == 0 .or. index_temperature >= metadata%nt) goto 900
         do field_index=1,22
            read(source_tokens(field_index),*,iostat=io_status) source_row(field_index)
            if (io_status /= 0) goto 900
         end do
         index_temperature=index_temperature+1; logt_value=source_row(1)
         if (.not. ieee_is_finite(logt_value)) goto 900
         if (index_q == 1) then
            logt_axis(index_temperature)=logt_value
            if (index_temperature == 1) first_logt=logt_value
            if (index_temperature == 2) delta_logt=logt_value-first_logt
            if (index_temperature >= 2) then
               if (delta_logt <= 0 .or. abs(logt_value-logt_axis(index_temperature-1)-delta_logt) > 1e-10_dp) goto 900
            end if
         else
            if (abs(logt_value-logt_axis(index_temperature)) > 1e-10_dp) goto 900
         end if
         temperature=10.0_dp**logt_axis(index_temperature)
         density=10.0_dp**(logq_value+2*logt_axis(index_temperature)-12)
         delta_density=source_row(19)-(logq_value+2*logt_axis(index_temperature)-12)
         if (ends_with(path,'mesa-FreeEOS_02z70x.data') .and. index_q == 521) then
            if (.not. ieee_is_finite(delta_density)) goto 900
            if (abs(delta_density+1e-4_dp) > 1e-10_dp) goto 900
         end if
         accepted=source_row(18) == 0 .and. logt_axis(index_temperature) >= 3.2_dp-1e-12_dp
         accepted=accepted .and. ieee_is_finite(delta_density) .and. abs(delta_density) < 1e-10_dp
         accepted=accepted .and. all(ieee_is_finite(source_row(2:11)))
         if (.not. accepted) cycle
         gas_pressure=10.0_dp**source_row(2); energy=10.0_dp**source_row(3); entropy=10.0_dp**source_row(4)
         radiation_pressure=metadata%crad*temperature**4/3; pressure=gas_pressure+radiation_pressure; radiation_scale=radiation_pressure/density
         potential_value=energy-temperature*entropy
         temperature_first=-temperature*entropy
         temperature_second=-temperature*(entropy+source_row(8))
         potential_value=potential_value+radiation_scale; density_first=gas_pressure/density; temperature_first=temperature_first+4*radiation_scale
         density_second=(pressure*source_row(5)-gas_pressure)/density
         mixed_second=(pressure*source_row(6)-4*radiation_pressure)/density
         temperature_second=temperature_second+16*radiation_scale
         accepted=pressure > 0 .and. energy > 0 .and. entropy > 0 .and. source_row(8) > 0 .and. source_row(5) > 0
         accepted=accepted .and. all(ieee_is_finite([potential_value,density_first,temperature_first,density_second,mixed_second,temperature_second]))
         if (.not. accepted) cycle
         nodes(0,0,index_q,index_temperature)=potential_value
         nodes(1,0,index_q,index_temperature)=density_first
         nodes(0,1,index_q,index_temperature)=temperature_first+2*density_first
         nodes(2,0,index_q,index_temperature)=density_second
         nodes(1,1,index_q,index_temperature)=mixed_second+2*density_second
         nodes(0,2,index_q,index_temperature)=temperature_second+4*mixed_second+4*density_second
         source_valid(index_q,index_temperature)=all(ieee_is_finite(nodes(:,:,index_q,index_temperature)))
      end do
      if (expect_label .or. index_q /= metadata%nq .or. index_temperature /= metadata%nt) goto 900
      if (abs(source_header(5)-first_logt) > 1e-10_dp .or. abs(source_header(6)-logt_axis(index_temperature)) > 1e-10_dp) goto 900
      if (abs(source_header(7)-delta_logt) > 1e-10_dp .or. abs(source_header(9)-first_logq) > 1e-10_dp) goto 900
      if (abs(source_header(11)-delta_logq) > 1e-10_dp) goto 900
      if (.not. anomaly_nq .and. .not. anomaly_qmax) then
         if (abs(source_header(10)-logq_axis(index_q)) > 1e-10_dp) goto 900
      end if
      if (anomaly_nq .or. anomaly_qmax) then
         if (index_q /= 527 .or. abs(logq_axis(index_q)-5.69_dp) > 1e-10_dp) goto 900
      end if
      status=table_mismatch
      if (.not. close(first_logq*ln10,metadata%u_min) .or. .not. close(delta_logq*ln10,metadata%du)) goto 900
      if (.not. close(first_logt*ln10,metadata%tau_min) .or. .not. close(delta_logt*ln10,metadata%dtau)) goto 900
      call cross_support(source_valid,first_support)
      call cross_support(first_support,complete)
      do index_temperature=1,metadata%nt
         do index_q=1,metadata%nq
            if (.not. complete(index_q,index_temperature)) nodes(:,:,index_q,index_temperature)=0
         end do
      end do
      status=table_ok
900   close(source_unit)
   end subroutine reconstruct_source

   subroutine cross_support(source_valid, support)
      logical, intent(in) :: source_valid(:,:)
      logical, allocatable, intent(out) :: support(:,:)
      integer :: index_q,index_temperature,nq,nt
      nq=size(source_valid,1); nt=size(source_valid,2)
      allocate(support(nq,nt)); support=.false.
      do index_temperature=3,nt-2
         do index_q=3,nq-2
            support(index_q,index_temperature)=all(source_valid(index_q-2:index_q+2,index_temperature)) .and. all(source_valid(index_q,index_temperature-2:index_temperature+2))
         end do
      end do
   end subroutine cross_support

   logical function close(a,b)
      real(dp), intent(in) :: a,b
      close=abs(a-b) <= 32*epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(b))
   end function close

   logical function ends_with(path,tail)
      character(len=*), intent(in) :: path,tail
      integer :: path_length,tail_length
      path_length=len_trim(path); tail_length=len_trim(tail); ends_with=.false.
      if (path_length < tail_length) return
      ends_with=path(path_length-tail_length+1:path_length) == tail(:tail_length)
   end function ends_with
end module freeeos_source
