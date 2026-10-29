program table_benchmark
   use const_def, only: dp
   use freeeos_table
   use auto_diff_real_2var_order3_module
   use benchmark_support
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: fields(4)=[10,10,43,43]
   character(len=16), parameter :: names(4)=[character(len=16) :: 'real_eval','ad_eval','real_state','ad_state']
   type(biquintic_table) :: table
   type(table_metadata) :: metadata
   real(dp) :: scales(43,4), real_values(43), ad_values(43), difference(43), maximum_difference(43)
   real(dp) :: raw_sums(43,2), started, process_started, target
   integer :: status, real_status, ad_status, query, repeats, field
   character(len=4096) :: cache_path, query_path
   character(len=64) :: expected_sha
   character(len=128) :: option
   logical :: preflight_only

   process_started=wall_time()
   if (command_argument_count() < 5 .or. command_argument_count() > 6) &
      error stop 'usage: table_benchmark CACHE QUERIES EXPECTED_SHA REPEATS TARGET_SECONDS [--preflight-only]'
   preflight_only=command_argument_count() == 6
   if (preflight_only) then
      call get_command_argument(6,option)
      if (trim(option) /= '--preflight-only') error stop 'unknown benchmark option'
   end if
   call open_report()
   call get_command_argument(1,cache_path)
   call get_command_argument(2,query_path)
   call get_command_argument(3,expected_sha)
   call read_controls(4,5,repeats,target)
   call read_queries(trim(query_path))
   call memory_snapshot('before_load')
   started=wall_time()
   call load_table(trim(cache_path),0.7_dp,0.02_dp,table,status,expected_sha)
   write(report_unit,*) 'phase_seconds',' cache_load ',wall_time()-started
   if (status /= table_ok) then
      write(report_unit,*) 'load_status',status
      error stop 'existing table rejected'
   end if
   metadata=get_metadata(table)
   write(report_unit,*) 'metadata_grid',metadata%nq,metadata%nt
   write(report_unit,*) 'metadata_composition',metadata%hydrogen,metadata%metals,metadata%crad,metadata%split
   write(report_unit,*) 'metadata_preprocessing',metadata%radius,' '//trim(metadata%route)//' '//metadata%preprocessing_sha256
   write(report_unit,*) 'metadata_source ',metadata%source_sha256
   write(report_unit,*) 'payload_bytes',int(metadata%nq,int64)*metadata%nt* &
      (9_int64*storage_size(0.0_dp)/8+storage_size(.true.)/8)+ &
      int(metadata%nq+metadata%nt,int64)*storage_size(0.0_dp)/8
   call memory_snapshot('after_load')
   scales=1
   maximum_difference=0
   raw_sums=0
   started=wall_time()
   do query=1,size(coordinates,2)
      call compute(3,query,real_values,real_status)
      call compute(4,query,ad_values,ad_status)
      if (real_status /= 0 .or. ad_status /= 0) then
         write(report_unit,*) 'rejected_query',query,real_status,ad_status
         error stop 'all four modes require the same successful physical-state corpus'
      end if
      if (.not. all(ieee_is_finite(real_values)) .or. .not. all(ieee_is_finite(ad_values))) &
         error stop 'nonfinite validation output'
      difference=abs(real_values-ad_values)/max(1.0_dp,abs(real_values),abs(ad_values))
      maximum_difference=max(maximum_difference,difference)
      scales(:,3)=max(scales(:,3),abs(real_values),abs(ad_values))
      raw_sums(:,1)=raw_sums(:,1)+real_values
      raw_sums(:,2)=raw_sums(:,2)+ad_values
      write(report_unit,*) 'validation',query,real_status,ad_status,maxval(difference(1:10)),maxval(difference(11:43))
      if (preflight_only) write(report_unit,*) 'trace',query,real_values,ad_values
      if (query == 1) then
         write(report_unit,*) 'first_real',real_values
         write(report_unit,*) 'first_ad',ad_values
      end if
   end do
   scales(:,4)=scales(:,3)
   scales(1:10,1)=scales(1:10,3)
   scales(1:10,2)=scales(1:10,3)
   write(report_unit,*) 'phase_seconds',' validation ',wall_time()-started
   write(report_unit,*) 'max_scaled_difference',maximum_difference
   write(report_unit,*) 'raw_sum_real',raw_sums(:,1)
   write(report_unit,*) 'raw_sum_ad',raw_sums(:,2)
   do field=1,4
      write(report_unit,*) 'normalization ',trim(names(field)),scales(1:fields(field),field)
   end do
   write(report_unit,*) 'phase_seconds',' startup_to_validated ',wall_time()-process_started
   if (preflight_only) then
      write(report_unit,*) 'preflight_only',1
   else
      call benchmark_modes(compute,names,fields,scales,repeats,target)
   end if
   call memory_snapshot('finished')
   write(report_unit,*) 'phase_seconds',' total ',wall_time()-process_started
contains
   subroutine compute(mode, query, values, status)
      integer, intent(in) :: mode, query
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: status
      type(potential_partials) :: partials
      type(thermodynamic_state) :: state
      type(auto_diff_real_2var_order3) :: x_ad, t_ad, potential_ad
      real(dp) :: vector(10), physical_x, physical_t
      physical_x=coordinates(1,query)
      physical_t=coordinates(2,query)
      if (mode == 1 .or. mode == 3) then
         call evaluate_table_real(table,physical_x,physical_t,partials,status)
         vector=partial_vector(partials)
      else
         x_ad=physical_x
         x_ad%d1val2=1.0_dp
         t_ad=physical_t
         t_ad%d1val1=1.0_dp
         call evaluate_table_ad(table,x_ad,t_ad,potential_ad,status)
         vector=ad_vector(potential_ad)
         if (mode == 4) partials=potential_partials(vector(1),vector(2),vector(3),vector(4), &
            vector(5),vector(6),vector(7),vector(8),vector(9),vector(10))
      end if
      values(1:10)=vector
      if (status /= 0 .or. mode <= 2) return
      call derive_state(partials,physical_x,physical_t,metadata%crad,state,status)
      values(11:43)=[state%value,state%dlnrho,state%dlnT]
   end subroutine compute
end program table_benchmark
