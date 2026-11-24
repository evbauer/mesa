program test_reader
   use const_def, only: dp
   use freeeos_table, only: table_metadata, table_ok
   use freeeos_supplement, only: read_supplement
   use runtime_sha256, only: file_sha256
   implicit none
   character(len=1024) :: path,argument
   character(len=64) :: digest
   type(table_metadata) :: metadata
   real(dp), allocatable :: mixed(:,:,:)
   logical, allocatable :: valid(:,:)
   integer :: status,expected_status,io,iu,it,k
   logical :: accepted
   call get_command_argument(1,path)
   call get_command_argument(2,argument)
   read(argument,*) expected_status
   call file_sha256(trim(path),digest,io)
   if (io /= 0 .and. expected_status /= 10) error stop 'fixture hash failed'
   call read_supplement(trim(path),0.7_dp,0.02_dp,0,repeat('a',64),digest,metadata,mixed,valid,status)
   if (status /= expected_status) then
      write(*,*) 'expected/status',expected_status,status
      error stop 'reader status'
   end if
   if (status /= table_ok) then
      write(*,'(a,i0)') 'PASS recoverable reader error ',status
      stop
   end if
   if (metadata%nq /= 7 .or. metadata%nt /= 5) error stop 'fixture dimensions'
   do it=1,5
      do iu=1,7
         accepted=mod(iu+2*it,4) /= 0
         if (valid(iu,it) .neqv. accepted) error stop 'mask orientation'
         do k=1,3
            if (accepted) then
               if (mixed(k,iu,it) /= 100000*k+10*iu+it+0.125_dp) error stop 'float64 orientation'
            else
               if (mixed(k,iu,it) /= 0) error stop 'invalid placeholder'
            end if
         end do
      end do
   end do
   write(*,'(a)') 'PASS nonsquare 7x5 asymmetric float64/mask/numeric metadata'
end program test_reader
