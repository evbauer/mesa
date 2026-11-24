! Historical comparison tests only. Never list this module in eos/Makefile.
module historical_source_test
   use const_def, only:dp
   use freeeos_source, only:reconstruct_source
   use freeeos_supplement, only:read_supplement
   use freeeos_table, only:table_metadata,biquintic_table,assemble_table,table_ok,table_mismatch
   implicit none
   private
   public :: load_historical_source_supplement
contains
   subroutine load_historical_source_supplement(source_path,supplement_path,x,z,family_id,source_sha, &
      supplement_sha,table,status)
      character(len=*),intent(in)::source_path,supplement_path,source_sha,supplement_sha
      real(dp),intent(in)::x,z
      integer,intent(in)::family_id
      type(biquintic_table),intent(out)::table
      integer,intent(out)::status
      type(table_metadata)::metadata
      real(dp),allocatable::mixed(:,:,:),nodes(:,:,:,:)
      logical,allocatable::valid(:,:),source_valid(:,:)
      call read_supplement(supplement_path,x,z,family_id,source_sha,supplement_sha,metadata,mixed,valid,status)
      if(status/=table_ok)return
      call reconstruct_source(source_path,metadata,nodes,source_valid,status)
      if(status/=table_ok)return
      status=table_mismatch
      if(any(valid.neqv.source_valid))return
      nodes(2,1,:,:)=mixed(1,:,:);nodes(1,2,:,:)=mixed(2,:,:);nodes(2,2,:,:)=mixed(3,:,:)
      call assemble_table(metadata,nodes,valid,table,status)
   end subroutine load_historical_source_supplement
end module historical_source_test
