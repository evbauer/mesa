module benchmark_support
   use const_def, only: dp
   use iso_fortran_env, only: int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer :: report_unit
   real(dp), allocatable :: coordinates(:,:)
   abstract interface
      subroutine query_evaluator(mode, query, values, status)
         import dp
         integer, intent(in) :: mode, query
         real(dp), intent(out) :: values(:)
         integer, intent(out) :: status
      end subroutine query_evaluator
   end interface
contains
   subroutine open_report()
      integer :: status
      open(newunit=report_unit, file='/dev/stdout', action='write', iostat=status)
      if (status /= 0) error stop 'cannot open benchmark report'
      open(unit=output_unit, file='/dev/stderr', action='write', iostat=status)
      if (status /= 0) error stop 'cannot redirect library diagnostics'
      write(report_unit,*) 'timer_resolution_seconds', clock_resolution()
   end subroutine open_report

   real(dp) function wall_time()
      integer(int64) :: ticks, rate
      call system_clock(ticks, rate)
      wall_time=real(ticks,dp)/real(rate,dp)
   end function wall_time

   real(dp) function clock_resolution()
      integer(int64) :: rate
      call system_clock(count_rate=rate)
      clock_resolution=1.0_dp/real(rate,dp)
   end function clock_resolution

   subroutine read_queries(path)
      character(len=*), intent(in) :: path
      integer :: unit, status, query, count
      real(dp) :: started
      started=wall_time()
      open(newunit=unit,file=path,status='old',action='read',iostat=status)
      if (status /= 0) error stop 'cannot open saved queries'
      read(unit,*,iostat=status) count
      if (status /= 0) error stop 'missing query count'
      if (count < 1 .or. count > 100000) error stop 'query count outside bounded limit'
      allocate(coordinates(2,count))
      do query=1,count
         read(unit,*,iostat=status) coordinates(:,query)
         if (status /= 0) error stop 'missing query coordinates'
      end do
      close(unit)
      if (.not. all(ieee_is_finite(coordinates))) error stop 'nonfinite coordinates'
      write(report_unit,*) 'query_count',count
      write(report_unit,*) 'phase_seconds',' query_read ',wall_time()-started
   end subroutine read_queries

   subroutine memory_snapshot(phase)
      character(len=*), intent(in) :: phase
      character(len=256) :: record
      character(len=32) :: label
      integer :: unit, status, rss, peak, amount
      rss=-1
      peak=-1
      open(newunit=unit,file='/proc/self/status',status='old',action='read',iostat=status)
      if (status == 0) then
         do
            read(unit,'(a)',iostat=status) record
            if (status /= 0) exit
            if (index(record,'VmRSS:') == 1 .or. index(record,'VmHWM:') == 1) then
               read(record,*,iostat=status) label,amount
               if (status == 0 .and. label == 'VmRSS:') rss=amount
               if (status == 0 .and. label == 'VmHWM:') peak=amount
            end if
         end do
         close(unit)
      end if
      write(report_unit,*) 'memory_kib ',trim(phase),rss,peak
   end subroutine memory_snapshot

   subroutine run_batch(evaluate, mode, fields, scale, passes, elapsed, cpu_elapsed, checksum, failures)
      procedure(query_evaluator) :: evaluate
      integer, intent(in) :: mode, fields, passes
      real(dp), intent(in) :: scale(:)
      real(dp), intent(out) :: elapsed, cpu_elapsed, checksum
      integer, intent(out) :: failures
      real(dp) :: values(fields), started, cpu_started, cpu_finished
      integer :: pass, query, status
      checksum=0
      failures=0
      call cpu_time(cpu_started)
      started=wall_time()
      do pass=1,passes
         do query=1,size(coordinates,2)
            call evaluate(mode,query,values,status)
            if (status /= 0) then
               failures=failures+1
            else
               checksum=checksum+sum(values/scale(1:fields))
            end if
         end do
      end do
      elapsed=wall_time()-started
      call cpu_time(cpu_finished)
      cpu_elapsed=cpu_finished-cpu_started
      if (.not. ieee_is_finite(checksum)) error stop 'nonfinite benchmark checksum'
   end subroutine run_batch

   subroutine benchmark_modes(evaluate, names, fields, scales, repeats, target)
      procedure(query_evaluator) :: evaluate
      character(len=*), intent(in) :: names(:)
      integer, intent(in) :: fields(:), repeats
      real(dp), intent(in) :: scales(:,:), target
      integer :: passes(size(names)), mode, repeat, position, failures
      real(dp) :: elapsed, cpu_elapsed, checksum
      do mode=1,size(names)
         call run_batch(evaluate,mode,fields(mode),scales(:,mode),1,elapsed,cpu_elapsed,checksum,failures)
         write(report_unit,*) 'warmup ',trim(names(mode)),elapsed,cpu_elapsed,checksum,failures
         if (failures /= 0) error stop 'warmup failed'
         call run_batch(evaluate,mode,fields(mode),scales(:,mode),1,elapsed,cpu_elapsed,checksum,failures)
         passes(mode)=max(1,min(10000,ceiling(target/max(elapsed,clock_resolution()))))
         write(report_unit,*) 'calibration ',trim(names(mode)),passes(mode),elapsed
      end do
      call memory_snapshot('warmed')
      do repeat=1,repeats
         do position=1,size(names)
            mode=1+mod(position+repeat-2,size(names))
            call run_batch(evaluate,mode,fields(mode),scales(:,mode),passes(mode),elapsed,cpu_elapsed,checksum,failures)
            write(report_unit,*) 'timing ',trim(names(mode)),repeat,passes(mode),size(coordinates,2), &
               elapsed,cpu_elapsed,checksum,failures
            flush(report_unit)
            if (failures /= 0) error stop 'timed evaluation failed'
         end do
      end do
   end subroutine benchmark_modes

   subroutine read_controls(repeat_argument, target_argument, repeats, target)
      integer, intent(in) :: repeat_argument, target_argument
      integer, intent(out) :: repeats
      real(dp), intent(out) :: target
      character(len=128) :: argument
      integer :: status
      call get_command_argument(repeat_argument,argument)
      read(argument,*,iostat=status) repeats
      if (status /= 0) error stop 'invalid repeat count'
      call get_command_argument(target_argument,argument)
      read(argument,*,iostat=status) target
      if (status /= 0) error stop 'invalid target duration'
      if (repeats < 3 .or. repeats > 15) error stop 'repeats must be 3..15'
      if (.not. ieee_is_finite(target)) error stop 'nonfinite target duration'
      if (target < 0.001_dp .or. target > 1.0_dp) error stop 'target must be .001..1 seconds'
   end subroutine read_controls
end module benchmark_support
