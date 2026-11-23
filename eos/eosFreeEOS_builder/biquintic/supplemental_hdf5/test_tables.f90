program test_tables
   use const_def, only: dp
   use freeeos_table
   use freeeos_supplement, only: load_source_supplement
   implicit none
   character(len=2048) :: source_dir,supplement_dir,reference_dir,filename,source_path,reference_path
   character(len=80) :: magic
   character(len=64) :: source_sha,supplement_sha
   type(biquintic_table) :: table,reference
   type(table_metadata) :: metadata
   type(potential_partials) :: potential,oracle
   type(thermodynamic_state) :: state,oracle_state
   real(dp), allocatable :: nodes(:,:,:,:),old_nodes(:,:,:,:)
   logical, allocatable :: valid(:,:),old_valid(:,:)
   real(dp) :: x,z,rho,tau,p(10),o(10),node_max(9),potential_max(10),state_max(11),partial_max(11)
   real(dp) :: delta(9),denominator(9),flat(9),old_flat(9)
   real(dp) :: node_abs(9),potential_abs(10),state_abs(11),partial_abs(11),scaled(11),actual_partial(11),old_partial(11)
   real(dp) :: worst_partial,partial_context(7),worst_node,node_context(6),relative
   integer :: partial_field,partial_axis,node_field,max_field,axis
   integer :: unit,entry_count,entry,available,status,old_status,iu,it,k,samples,accepted_states,nonphysical
   integer :: node_index,family_id,stride
   call get_command_argument(1,source_dir)
   call get_command_argument(2,supplement_dir)
   call get_command_argument(3,reference_dir)
   open(newunit=unit,file=trim(supplement_dir)//'/supplements.manifest',status='old',action='read')
   read(unit,*) magic
   if (magic /= 'FREEEOS_BIQUINTIC_SUPPLEMENTS_V1') error stop 'manifest'
   read(unit,*) entry_count
   do entry=1,entry_count
      read(unit,*) x,z,available,filename,source_sha,supplement_sha
      k=len_trim(filename)
      source_path=trim(source_dir)//'/'//filename(:k-6)//'.data'
      reference_path=trim(reference_dir)//'/'//filename(:k-3)
      family_id=0
      if (index(trim(filename),'_CO.bq.h5') > 0) family_id=1
      call load_source_supplement(trim(source_path),trim(supplement_dir)//'/'//trim(filename), &
         x,z,family_id,source_sha,supplement_sha,table,status)
      if (status /= table_ok) then
         write(*,*) 'load failure',trim(filename),status
         error stop 'source+HDF5 load'
      end if
      call load_table(trim(reference_path),x,z,reference,status,source_sha)
      if (status /= table_ok) error stop 'V1 load'
      call get_node_constraints(table,nodes,valid,status)
      call get_node_constraints(reference,old_nodes,old_valid,status)
      if (any(valid .neqv. old_valid)) error stop 'accepted masks changed'
      if (any(nodes(2,1,:,:) /= old_nodes(2,1,:,:))) error stop 'mixed21 changed'
      if (any(nodes(1,2,:,:) /= old_nodes(1,2,:,:))) error stop 'mixed12 changed'
      if (any(nodes(2,2,:,:) /= old_nodes(2,2,:,:))) error stop 'mixed22 changed'
      node_max=0; potential_max=0; state_max=0; partial_max=0; samples=0; accepted_states=0; nonphysical=0
      node_abs=0; potential_abs=0; state_abs=0; partial_abs=0
      worst_partial=-1; worst_node=-1; partial_context=0; node_context=0
      partial_field=0; partial_axis=0; node_field=0
      metadata=get_metadata(table)
      stride=max(1,count(valid)/64)
      node_index=0
      do it=1,metadata%nt
         do iu=1,metadata%nq
            if (.not. valid(iu,it)) cycle
            flat=reshape(nodes(:,:,iu,it),[9]); old_flat=reshape(old_nodes(:,:,iu,it),[9])
            node_max=max(node_max,abs(flat-old_flat)/max(1.0_dp,abs(old_flat)))
            node_abs=max(node_abs,abs(flat-old_flat))
            delta=abs(flat-old_flat)/max(1.0_dp,abs(old_flat))
            max_field=maxloc(delta,dim=1)
            if (delta(max_field) > worst_node) then
               worst_node=delta(max_field); node_field=max_field
               node_context=[real(iu,dp),real(it,dp),flat(max_field),old_flat(max_field), &
                  abs(flat(max_field)-old_flat(max_field)),worst_node]
            end if
            node_index=node_index+1
            if (iu >= metadata%nq .or. it >= metadata%nt) cycle
            if (.not. all(valid(iu:iu+1,it:it+1))) cycle
            if (mod(node_index,stride) /= 0) cycle
            tau=metadata%tau_min+(it-1+0.37_dp)*metadata%dtau
            rho=metadata%u_min+(iu-1+0.61_dp)*metadata%du+2*tau-12*log(10.0_dp)
            call evaluate_table_real(table,rho,tau,potential,status)
            call evaluate_table_real(reference,rho,tau,oracle,old_status)
            if (status /= old_status .or. status /= table_ok) error stop 'evaluator status'
            p=partial_vector(potential); o=partial_vector(oracle)
            potential_max=max(potential_max,abs(p-o)/max(1.0_dp,abs(o)))
            potential_abs=max(potential_abs,abs(p-o))
            samples=samples+1
            call derive_state(potential,rho,tau,metadata%crad,state,status)
            call derive_state(oracle,rho,tau,metadata%crad,oracle_state,old_status)
            if (status /= old_status) then
               write(*,*) 'state status differs',trim(filename),status,old_status,rho,tau
               error stop 'state status changed'
            end if
            if (status /= table_ok) then
               nonphysical=nonphysical+1
               cycle
            end if
            accepted_states=accepted_states+1
            state_max=max(state_max,abs(state%value-oracle_state%value)/max(1.0_dp,abs(oracle_state%value)))
            state_abs=max(state_abs,abs(state%value-oracle_state%value))
            partial_max=max(partial_max,abs(state%dlnrho-oracle_state%dlnrho)/max(1.0_dp,abs(oracle_state%dlnrho)))
            partial_max=max(partial_max,abs(state%dlnT-oracle_state%dlnT)/max(1.0_dp,abs(oracle_state%dlnT)))
            do axis=1,2
               if (axis == 1) then
                  actual_partial=state%dlnrho; old_partial=oracle_state%dlnrho
               else
                  actual_partial=state%dlnT; old_partial=oracle_state%dlnT
               end if
               partial_abs=max(partial_abs,abs(actual_partial-old_partial))
               scaled=abs(actual_partial-old_partial)/max(1.0_dp,abs(old_partial))
               max_field=maxloc(scaled,dim=1)
               if (scaled(max_field) > worst_partial) then
                  worst_partial=scaled(max_field); partial_field=max_field; partial_axis=axis
                  relative=huge(1.0_dp)
                  if (old_partial(max_field) /= 0) &
                     relative=abs(actual_partial(max_field)-old_partial(max_field))/abs(old_partial(max_field))
                  partial_context=[rho/log(10.0_dp),tau/log(10.0_dp),actual_partial(max_field), &
                     old_partial(max_field),abs(actual_partial(max_field)-old_partial(max_field)),worst_partial,relative]
               end if
            end do
         end do
      end do
      write(*,'(a,1x,a,3(1x,i0),*(1x,es18.10))') 'TABLE',trim(filename),samples,accepted_states,nonphysical, &
         node_max,potential_max,state_max,partial_max,node_abs,potential_abs,state_abs,partial_abs
      write(*,'(a,1x,a,2(1x,i0),*(1x,es25.16))') 'PARTIAL_OUTLIER',trim(filename),partial_field,partial_axis,partial_context
      write(*,'(a,1x,a,1x,i0,*(1x,es25.16))') 'NODE_OUTLIER',trim(filename),node_field,node_context
      flush(6)
   end do
   close(unit)
   write(*,'(a,i0)') 'PASS all source+HDF5/V1 masks and mixed constraints: ',entry_count
end program test_tables
