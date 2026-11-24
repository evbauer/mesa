program family_tests
   use bridge_support
   use freeeos_family
   use freeeos_table, only: table_ok
   use omp_lib
   use iso_fortran_env, only:int64
   implicit none
   character(len=4096)::mesa_dir,inlist,massfile,directory
   type(EoS_General_Info),pointer::rq
   real(dp)::res(nv),dr(nv),dt(nv),da(nv,pinned_species),reference(nv,23,2)
   integer::handles(10),i,j,k,h,status,errors,tid
   integer(int64)::successes,failures,loads
   logical::linear
   call get_command_argument(1,mesa_dir);call get_command_argument(2,inlist);call get_command_argument(3,massfile)
   call get_command_argument(4,directory)
   call initialize_eos(trim(mesa_dir),trim(massfile),status)
   if(status/=0)error stop 'init'
   do i=1,10
      handles(i)=alloc_eos_handle_using_inlist(trim(inlist),status)
      if(status/=0)error stop 'handle'
      if(mod(i,2)==0)eos_handles(handles(i))%suffix_for_FreeEOS_Z=''
   enddo
   ! Missing context must fail before poisoning attempted/load state.
   res=0;dr=0;dt=0;da=0;res(i_mu:i_eta)=1
   rq=>eos_handles(handles(1))
   call freeeos_family_results(trim(directory),.false.,rq%tiny_fuzz,FreeEOS_XZ_struct, &
      0._dp,.95_dp,pinned_species,pinned_chem_id,1e-7_dp,1e5_dp,res,dr,dt,da,status,rq%suffix_for_FreeEOS_Z, &
      trim(mesa_dir)//'/data/eosFreeEOS_data')
   write(*,*)'MISSING_CONTEXT',status
   if(status/=family_legacy_error)error stop 'missing context'
   errors=0
   ! Real request context, concurrent first use including both suffixes.
!$omp parallel do num_threads(16) private(i,h,res,dr,dt,da,status) reduction(+:errors)
   do i=1,160
      h=handles(mod(i-1,10)+1)
      call query(h,.false.,0._dp,.95_dp,res,dr,dt,da,status)
      if(status/=0)errors=errors+1
   enddo
!$omp end parallel do
   if(errors/=0)error stop 'concurrent first use'
   call freeeos_family_get_counters(successes,failures,loads)
   write(*,*)'COLD_COUNTS',successes,failures,loads
   do j=1,2
      call query(handles(j),.false.,0._dp,.95_dp,res,dr,dt,da,status)
      if(status/=0)error stop 'reference'
      reference(:,1,j)=res;reference(:,2,j)=dr;reference(:,3,j)=dt;reference(:,4:23,j)=da
   enddo
   if(all(reference(:,:,1)==reference(:,:,2)))error stop 'suffix alias'
!$omp parallel do num_threads(16) private(i,j,h,res,dr,dt,da,status) reduction(+:errors)
   do i=1,320
      j=mod(i,2)+1;h=handles(j)
      call query(h,.false.,0._dp,.95_dp,res,dr,dt,da,status)
      if(status/=0.or.any(res/=reference(:,1,j)).or.any(dr/=reference(:,2,j)).or. &
         any(dt/=reference(:,3,j)).or.any(da/=reference(:,4:23,j)))errors=errors+1
   enddo
!$omp end parallel do
   if(errors/=0)error stop 'warm concurrent returns'
   write(*,*)'WARM_THREAD_EXACT',320
   do k=1,2
      linear=k==2
      do j=1,2
         do i=1,5
            select case(i)
            case(1);call query(handles(j),linear,.7_dp,.02_dp,res,dr,dt,da,status)
            case(2);call query(handles(j),linear,.65_dp,.03_dp,res,dr,dt,da,status)
            case(3);call query(handles(j),linear,.3_dp,.1_dp,res,dr,dt,da,status)
            case(4);call query(handles(j),linear,.05_dp,.85_dp,res,dr,dt,da,status)
            case(5);call query(handles(j),linear,0._dp,.95_dp,res,dr,dt,da,status)
            end select
            if(status/=0)error stop 'nodal stencil corpus'
            write(*,'(a,3(1x,i0),*(1x,es25.16e3))')'CORPUS',k,j,i,res,dr,dt,da
         enddo
      enddo
   enddo
   rq=>eos_handles(handles(1))
   res=0;dr=0;dt=0;da=0;res(i_mu:i_eta)=1
   call freeeos_family_results(trim(directory),.false.,rq%tiny_fuzz,FreeEOS_XZ_struct, &
      .7_dp,.02_dp,pinned_species,pinned_chem_id,1e30_dp,100._dp,res,dr,dt,da,status,rq%suffix_for_FreeEOS_Z, &
      trim(mesa_dir)//'/data/eosFreeEOS_data',legacy_request=rq)
   write(*,*)'OUTSIDE_SUPPORT',status
   if(status==0)error stop 'outside support'
   call freeeos_family_get_counters(successes,failures,loads)
   write(*,*)'FINAL_COUNTS',successes,failures,loads
   write(*,*)'PASS real handles, cold/warm threads, suffix identity, cubic/linear composition, support and context errors'
contains
   subroutine query(handle,linear,x,z,res,dr,dt,da,status)
      integer,intent(in)::handle
      logical,intent(in)::linear
      real(dp),intent(in)::x,z
      real(dp),intent(out)::res(nv),dr(nv),dt(nv),da(nv,pinned_species)
      integer,intent(out)::status
      type(EoS_General_Info),pointer::rq
      rq=>eos_handles(handle)
      res=0;dr=0;dt=0;da=0;res(i_mu:i_eta)=1
      call freeeos_family_results(trim(directory),linear,rq%tiny_fuzz,FreeEOS_XZ_struct, &
         x,z,pinned_species,pinned_chem_id,1e-7_dp,1e5_dp,res,dr,dt,da,status,rq%suffix_for_FreeEOS_Z, &
         trim(mesa_dir)//'/data/eosFreeEOS_data',legacy_request=rq)
   end subroutine query
end program family_tests
