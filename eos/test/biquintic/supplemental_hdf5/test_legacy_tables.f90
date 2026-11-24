program all_tables
   use bridge_support
   use eos_lib, only:eos_shutdown
   use eosDT_load_tables, only:load_single_eosDT_table_by_id
   use freeeos_supplement, only:load_legacy_supplement
   use freeeos_table
   implicit none
   character(len=4096)::mesa_dir,inlist,massfile,directory,filename
   character(len=64)::sha,payload
   character(len=80)::magic
   type(EoS_General_Info),pointer::rq
   type(EosDT_XZ_Info),pointer::ep
   type(biquintic_table)::table
   real(dp),allocatable::nodes(:,:,:,:)
   logical,allocatable::valid(:,:)
   real(dp)::x,z
   integer::handle,status,ix,iz,k,u,n,available,family_id,ordinary,co,total_valid
   call get_command_argument(1,mesa_dir);call get_command_argument(2,inlist);call get_command_argument(3,massfile)
   call get_command_argument(4,directory)
   call initialize_eos(trim(mesa_dir),trim(massfile),status)
   if(status/=0)error stop 'init'
   handle=alloc_eos_handle_using_inlist(trim(inlist),status)
   if(status/=0)error stop 'handle'
   rq=>eos_handles(handle)
   open(newunit=u,file=trim(directory)//'/supplements.manifest',status='old')
   read(u,*)magic;read(u,*)n
   if(n/=89)error stop 'coverage'
   ordinary=0;co=0;total_valid=0
   do k=1,n
      read(u,*)x,z,available,filename,sha,payload
      if(available/=1)error stop 'unavailable delivered entry'
      family_id=0;rq%suffix_for_FreeEOS_Z=''
      if(index(filename,'_CO.')>0)then
         family_id=1;rq%suffix_for_FreeEOS_Z='_CO'
      endif
      do iz=1,FreeEOS_XZ_struct%nZs
         if(abs(FreeEOS_XZ_struct%Zs(iz)-z)<1e-14_dp)exit
      enddo
      if(iz>FreeEOS_XZ_struct%nZs)error stop 'Z mapping'
      do ix=1,FreeEOS_XZ_struct%nXs_for_Z(iz)
         if(abs(FreeEOS_XZ_struct%Xs_for_Z(ix,iz)-x)<1e-14_dp)exit
      enddo
      if(ix>FreeEOS_XZ_struct%nXs_for_Z(iz))error stop 'X mapping'
      call load_single_eosDT_table_by_id(rq,eosdt_max_FreeEOS,ep,ix,iz,status)
      if(status/=0)error stop 'legacy load'
      call load_legacy_supplement(ep,trim(directory)//'/'//trim(filename),x,z,family_id,sha,payload,table,status)
      write(*,'(a,1x,a,4(1x,i0))') 'TABLE',trim(filename),status,ep%version,ep%num_logQs,ep%num_logTs
      if(status/=0)error stop 'adapter/supplement validation'
      call get_node_constraints(table,nodes,valid,status)
      if(status/=0.or.any(shape(nodes)/=[3,3,ep%num_logQs,ep%num_logTs]))error stop 'shape'
      if(.not.all(ieee_is_finite(nodes)))error stop 'nonfinite'
      if(.not.any(valid))error stop 'empty support'
      do iz=1,size(valid,2)
         do ix=1,size(valid,1)
            if(.not.valid(ix,iz).and.any(nodes(:,:,ix,iz)/=0))error stop 'invalid placeholders'
         enddo
      enddo
      if(family_id==0)ordinary=ordinary+1
      if(family_id==1)co=co+1
      total_valid=total_valid+count(valid)
      write(*,*) 'MASK',trim(filename),count(valid),size(valid)
   enddo
   close(u)
   call eos_shutdown()
   if(associated(FreeEOS_alternate_tables))error stop 'alternate shutdown'
   if(ordinary/=83.or.co/=6)error stop 'family counts'
   write(*,*)'PASS_ALL_TABLES',n,ordinary,co,total_valid
end program all_tables
