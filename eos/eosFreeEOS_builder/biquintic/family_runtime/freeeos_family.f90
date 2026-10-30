module freeeos_family
   use const_def, only: dp
   use eos_def, only: DT_XZ_Info, nv, i_lnPgas, i_lnE, i_lnS, i_grad_ad, i_chiRho, i_chiT, &
      i_Cp, i_Cv, i_dE_dRho, i_dS_dT, i_dS_dRho, i_gamma1, i_gamma3, i_mu, i_eta
   use chem_def, only: chem_isos
   use freeeos_table, only: biquintic_table, table_metadata, potential_partials, thermodynamic_state, &
      load_table, get_metadata, evaluate_table_real, derive_state, partial_vector, table_ok, &
      state_pgas, state_energy, state_entropy, state_cv, state_cp, state_chirho, state_chit, &
      state_gamma1, state_gamma3, state_grad_ad
   use iso_fortran_env, only: int64, error_unit
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: freeeos_family_results, freeeos_family_get_counters
   public :: family_invalid_input, family_manifest_error, family_directory_changed, &
      family_grid_mismatch, family_unavailable, family_legacy_error

   integer, parameter :: family_invalid_input=301, family_manifest_error=302
   integer, parameter :: family_directory_changed=303, family_grid_mismatch=304
   integer, parameter :: family_unavailable=305, family_legacy_error=306
   integer, parameter :: max_entries=4096
   integer, parameter :: fields(13)=[i_lnPgas,i_lnE,i_lnS,i_grad_ad,i_chiRho,i_chiT, &
      i_Cp,i_Cv,i_dE_dRho,i_dS_dT,i_dS_dRho,i_gamma1,i_gamma3]
   type table_entry
      real(dp) :: hydrogen, metals
      integer :: available
      character(len=256) :: filename
      character(len=64) :: source_sha, payload_sha
      logical :: loaded=.false.
      logical :: attempted=.false.
      integer :: load_status=table_ok
      type(biquintic_table) :: table
   end type table_entry
   type(table_entry), allocatable, save :: entries(:)
   type(DT_XZ_Info), save :: saved_grid
   character(len=:), allocatable, save :: family_directory
   integer(int64), save :: successful_queries=0, failed_queries=0, table_loads=0

contains

   subroutine freeeos_family_get_counters(successes, failures, loads)
      integer(int64), intent(out) :: successes, failures, loads
!$omp atomic read
      successes=successful_queries
!$omp atomic read
      failures=failed_queries
!$omp atomic read
      loads=table_loads
   end subroutine freeeos_family_get_counters

   subroutine freeeos_family_results(directory, linear_x, tiny_fuzz, xz, &
         hydrogen, metals, species, chem_id, rho, temperature, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)
      character(len=*), intent(in) :: directory
      logical, intent(in) :: linear_x
      real(dp), intent(in) :: tiny_fuzz, hydrogen, metals, rho, temperature
      type(DT_XZ_Info), intent(in) :: xz
      integer, intent(in) :: species, chem_id(:)
      real(dp), intent(inout) :: res(nv), d_dlnd(nv), d_dlnT(nv)
      real(dp), intent(inout) :: d_dxa(nv,species)
      integer, intent(out) :: ierr
      real(dp) :: weights(4,2), dx_weights(4,2), dz_weights(4,2)
      real(dp) :: z_weights(2), z_derivatives(2)
      real(dp) :: combined(10), by_x(10), by_z(10), components(10)
      real(dp) :: reference(10)
      real(dp) :: mapped(nv), mapped_rho(nv), mapped_t(nv), mapped_x(nv), mapped_z(nv)
      real(dp) :: crad, xcoord, zcoord
      integer :: iz(2), ix(4,2), nz, nx(2), layer, point, entry_index, status, species_index
      type(potential_partials) :: potential
      type(thermodynamic_state) :: state
      type(table_metadata) :: metadata

      ierr=family_invalid_input
      iz=0; ix=0; nz=0; nx=0
      if (len_trim(directory) == 0 .or. species < 0 .or. size(chem_id) /= species) goto 900
      if (.not. all(ieee_is_finite([hydrogen,metals,rho,temperature,tiny_fuzz]))) goto 900
      if (rho <= 0 .or. temperature <= 0 .or. tiny_fuzz < 0) goto 900
      if (.not. associated(chem_isos%Z)) goto 900
      if (any(chem_id < 1) .or. any(chem_id > size(chem_isos%Z))) goto 900
      ierr=family_legacy_error
      if (.not. all(ieee_is_finite(res(i_mu:i_eta))) .or. &
          .not. all(ieee_is_finite(d_dlnd(i_mu:i_eta))) .or. &
          .not. all(ieee_is_finite(d_dlnT(i_mu:i_eta))) .or. &
          .not. all(ieee_is_finite(d_dxa(i_mu:i_eta,:)))) goto 900
      ierr=family_grid_mismatch
      if (xz%nZs < 1 .or. xz%nZs > size(xz%Zs)) goto 900
      if (any(xz%nXs_for_Z(1:xz%nZs) < 1) .or. &
          any(xz%nXs_for_Z(1:xz%nZs) > size(xz%Xs_for_Z,1))) goto 900
      call select_z(xz,metals,iz,nz,z_weights,z_derivatives,ierr)
      if (ierr /= table_ok) goto 900
      do layer=1,nz
         call select_x(xz,iz(layer),hydrogen,linear_x,tiny_fuzz, &
            ix(:,layer),nx(layer),weights(:,layer),dx_weights(:,layer),ierr)
         if (ierr /= table_ok) goto 900
      end do
      do layer=1,nz
         do point=1,nx(layer)
            dz_weights(point,layer)=z_derivatives(layer)*weights(point,layer)
            dx_weights(point,layer)=z_weights(layer)*dx_weights(point,layer)
            weights(point,layer)=z_weights(layer)*weights(point,layer)
         end do
      end do
      combined=0; by_x=0; by_z=0; crad=-1
      do layer=1,nz
         do point=1,nx(layer)
            xcoord=xz%Xs_for_Z(ix(point,layer),iz(layer))
            zcoord=xz%Zs(iz(layer))
!$omp critical(freeeos_family_cache)
            call ensure_entry(directory,xz,xcoord,zcoord,entry_index,status)
!$omp end critical(freeeos_family_cache)
            ierr=status
            if (ierr /= table_ok) goto 900
            metadata=get_metadata(entries(entry_index)%table)
            if (metadata%split /= 1 .or. metadata%radius /= 2 .or. metadata%route /= 'mean') then
               ierr=family_grid_mismatch
               goto 900
            end if
            if (crad < 0) then
               crad=metadata%crad
            else
               if (crad /= metadata%crad) then
                  ierr=family_grid_mismatch
                  goto 900
               end if
            end if
            call evaluate_table_real(entries(entry_index)%table,log(rho),log(temperature),potential,ierr)
            if (ierr /= table_ok) goto 900
            components=partial_vector(potential)
            if (layer == 1 .and. point == 1) then
               reference=components
               combined=reference
            end if
            combined=combined+weights(point,layer)*(components-reference)
            by_x=by_x+dx_weights(point,layer)*(components-reference)
            by_z=by_z+dz_weights(point,layer)*(components-reference)
         end do
      end do
      potential=potential_partials(combined(1),combined(2),combined(3),combined(4),combined(5), &
         combined(6),combined(7),combined(8),combined(9),combined(10))
      call derive_state(potential,log(rho),log(temperature),crad,state,ierr)
      if (ierr /= table_ok) goto 900
      call map_family_state(potential,state,by_x,by_z,rho,temperature, &
         mapped,mapped_rho,mapped_t,mapped_x,mapped_z,ierr)
      if (ierr /= table_ok) goto 900
      res(fields)=mapped(fields)
      d_dlnd(fields)=mapped_rho(fields)
      d_dlnT(fields)=mapped_t(fields)
      do species_index=1,species
         select case(chem_isos%Z(chem_id(species_index)))
         case(1)
            d_dxa(fields,species_index)=mapped_x(fields)
         case(2)
            d_dxa(fields,species_index)=0
         case default
            d_dxa(fields,species_index)=mapped_z(fields)
         end select
      end do
!$omp atomic update
      successful_queries=successful_queries+1
      return
900   call record_failure(ierr,hydrogen,metals,rho,temperature,iz,ix,nz,nx)
   end subroutine freeeos_family_results

   subroutine record_failure(ierr,hydrogen,metals,rho,temperature,iz,ix,nz,nx)
      integer, intent(in) :: ierr, iz(2), ix(4,2), nz, nx(2)
      real(dp), intent(in) :: hydrogen, metals, rho, temperature
      integer(int64) :: count
!$omp atomic capture
      failed_queries=failed_queries+1
      count=failed_queries
!$omp end atomic
      if (count > 5) return
!$omp critical(freeeos_family_diagnostic)
      write(error_unit,*) 'freeeos_family failure',ierr,' X Z rho T',hydrogen,metals,rho,temperature
      write(error_unit,*) 'freeeos_family stencil nZ iz nX ix',nz,iz,nx,ix
!$omp end critical(freeeos_family_diagnostic)
   end subroutine record_failure

   subroutine ensure_entry(directory,xz,hydrogen,metals,index_entry,status)
      character(len=*), intent(in) :: directory
      type(DT_XZ_Info), intent(in) :: xz
      real(dp), intent(in) :: hydrogen, metals
      integer, intent(out) :: index_entry, status
      integer :: position
      index_entry=0
      if (.not. allocated(family_directory)) then
         call read_manifest(directory,xz,status)
         if (status /= table_ok) return
      else
         status=family_directory_changed
         if (trim(directory) /= family_directory) return
         status=family_grid_mismatch
         if (.not. same_grid(xz,saved_grid)) return
      end if
      status=family_manifest_error
      do position=1,size(entries)
         if (is_variant(entries(position)%filename)) cycle
         if (close_coord(hydrogen,entries(position)%hydrogen) .and. &
             close_coord(metals,entries(position)%metals)) then
            index_entry=position
            exit
         end if
      end do
      if (index_entry == 0) return
      status=family_unavailable
      if (entries(index_entry)%available /= 1) return
      if (entries(index_entry)%loaded) then
         status=table_ok
         return
      end if
      if (entries(index_entry)%attempted) then
         status=entries(index_entry)%load_status
         return
      end if
      entries(index_entry)%attempted=.true.
      call load_table(trim(directory)//'/'//trim(entries(index_entry)%filename),hydrogen,metals, &
         entries(index_entry)%table,status,entries(index_entry)%source_sha)
      entries(index_entry)%load_status=status
      if (status /= table_ok) return
      entries(index_entry)%loaded=.true.
!$omp atomic update
      table_loads=table_loads+1
   end subroutine ensure_entry

   subroutine read_manifest(directory,xz,status)
      character(len=*), intent(in) :: directory
      type(DT_XZ_Info), intent(in) :: xz
      integer, intent(out) :: status
      integer :: unit, io_status, count, position, layer, point, matches, ordinary_count
      character(len=80) :: magic
      type(table_entry), allocatable :: candidate(:)
      status=family_manifest_error
      open(newunit=unit,file=trim(directory)//'/family.manifest',status='old',action='read',iostat=io_status)
      if (io_status /= 0) return
      read(unit,*,iostat=io_status) magic
      if (io_status /= 0) goto 900
      if (trim(magic) /= 'FREEEOS_BIQUINTIC_FAMILY_V1') goto 900
      read(unit,*,iostat=io_status) count
      if (io_status /= 0) goto 900
      if (count < 1 .or. count > max_entries) goto 900
      allocate(candidate(count),stat=io_status)
      if (io_status /= 0) goto 900
      do position=1,count
         read(unit,*,iostat=io_status) candidate(position)%hydrogen,candidate(position)%metals, &
            candidate(position)%available,candidate(position)%filename, &
            candidate(position)%source_sha,candidate(position)%payload_sha
         if (io_status /= 0) goto 900
         if (.not. all(ieee_is_finite([candidate(position)%hydrogen,candidate(position)%metals]))) goto 900
         if (candidate(position)%available /= 0 .and. candidate(position)%available /= 1) goto 900
         if (candidate(position)%available == 1) then
            if (trim(candidate(position)%filename) == '-' .or. &
                len_trim(candidate(position)%source_sha) /= 64 .or. &
                len_trim(candidate(position)%payload_sha) /= 64) goto 900
            if (index(candidate(position)%filename,'/') /= 0 .or. &
                index(candidate(position)%filename,'..') /= 0) goto 900
         end if
      end do
      ordinary_count=0
      do position=1,count
         if (is_variant(candidate(position)%filename)) cycle
         ordinary_count=ordinary_count+1
      end do
      if (ordinary_count /= sum(xz%nXs_for_Z(1:xz%nZs))) goto 900
      do layer=1,xz%nZs
         do point=1,xz%nXs_for_Z(layer)
            matches=0
            do position=1,count
               if (is_variant(candidate(position)%filename)) cycle
               if (close_coord(xz%Zs(layer),candidate(position)%metals) .and. &
                   close_coord(xz%Xs_for_Z(point,layer),candidate(position)%hydrogen)) matches=matches+1
            end do
            if (matches /= 1) goto 900
         end do
      end do
      call move_alloc(candidate,entries)
      family_directory=trim(directory)
      saved_grid=xz
      status=table_ok
900   close(unit)
   end subroutine read_manifest

   logical function same_grid(left,right)
      type(DT_XZ_Info), intent(in) :: left,right
      integer :: layer
      same_grid=.false.
      if (left%nZs /= right%nZs) return
      if (any(left%nXs_for_Z(1:left%nZs) /= right%nXs_for_Z(1:left%nZs))) return
      if (any(left%Zs(1:left%nZs) /= right%Zs(1:left%nZs))) return
      do layer=1,left%nZs
         if (any(left%Xs_for_Z(1:left%nXs_for_Z(layer),layer) /= &
             right%Xs_for_Z(1:left%nXs_for_Z(layer),layer))) return
      end do
      same_grid=.true.
   end function same_grid

   logical function close_coord(left,right)
      real(dp), intent(in) :: left,right
      close_coord=abs(left-right) <= 32*epsilon(1.0_dp)*max(1.0_dp,abs(left),abs(right))
   end function close_coord

   logical function is_variant(filename)
      character(len=*), intent(in) :: filename
      integer :: length
      length=len_trim(filename)
      is_variant=.false.
      if (length < 6) return
      is_variant=filename(length-5:length) == '_CO.bq'
   end function is_variant

   subroutine select_z(xz,metals,indices,count,weight,dweight,status)
      type(DT_XZ_Info), intent(in) :: xz
      real(dp), intent(in) :: metals
      integer, intent(out) :: indices(2),count,status
      real(dp), intent(out) :: weight(2),dweight(2)
      integer :: layer
      status=table_ok; indices=0; weight=0; dweight=0
      if (xz%nZs == 1 .or. metals >= xz%Zs(xz%nZs)) then
         count=1; indices(1)=xz%nZs; weight(1)=1
         return
      end if
      count=2
      do layer=2,xz%nZs
         if (metals < xz%Zs(layer)) exit
      end do
      indices=[layer-1,layer]
      status=family_grid_mismatch
      if (xz%Zs(layer) <= xz%Zs(layer-1)) return
      dweight(2)=1/(xz%Zs(layer)-xz%Zs(layer-1))
      dweight(1)=-dweight(2)
      weight(2)=(metals-xz%Zs(layer-1))*dweight(2)
      weight(1)=1-weight(2)
      status=table_ok
   end subroutine select_z

   subroutine select_x(xz,layer,hydrogen,linear_x,tiny_fuzz,indices,count,weight,dweight,status)
      type(DT_XZ_Info), intent(in) :: xz
      integer, intent(in) :: layer
      real(dp), intent(in) :: hydrogen,tiny_fuzz
      logical, intent(in) :: linear_x
      integer, intent(out) :: indices(4),count,status
      real(dp), intent(out) :: weight(4),dweight(4)
      real(dp) :: spacing, spacing2, spacing3, delta, denominator, coefficient, derivative
      integer :: number, point, lower, upper
      indices=0; weight=0; dweight=0; status=table_ok
      number=xz%nXs_for_Z(layer)
      status=family_grid_mismatch
      if (xz%Xs_for_Z(1,layer) /= 0) return
      do point=2,number
         if (xz%Xs_for_Z(point,layer) <= xz%Xs_for_Z(point-1,layer)) return
      end do
      status=table_ok
      if (hydrogen < tiny_fuzz .or. number == 1) then
         count=1; indices(1)=1; weight(1)=1
         return
      end if
      if (hydrogen >= xz%Xs_for_Z(number,layer)) then
         count=1; indices(1)=number; weight(1)=1
         return
      end if
      if (.not. linear_x .and. number >= 3) then
         if (hydrogen <= xz%Xs_for_Z(2,layer)) then
            lower=1; upper=3
         else if (hydrogen >= xz%Xs_for_Z(number-1,layer)) then
            lower=number-2; upper=number
         else
            lower=0; upper=0
            do point=3,number-1
               if (hydrogen <= xz%Xs_for_Z(point,layer)) then
                  lower=point-2; upper=point+1
                  exit
               end if
            end do
         end if
         if (lower < 1) then
            status=family_grid_mismatch
            return
         end if
         spacing=xz%Xs_for_Z(lower+1,layer)-xz%Xs_for_Z(lower,layer)
         spacing2=xz%Xs_for_Z(lower+2,layer)-xz%Xs_for_Z(lower+1,layer)
         if (abs(spacing-spacing2) <= tiny_fuzz) then
            if (upper-lower == 2) then
               count=3; indices(1:3)=[lower,lower+1,upper]
               delta=hydrogen-xz%Xs_for_Z(lower,layer)
               denominator=2*spacing*spacing
               weight(1:3)=[(2*spacing*spacing-3*spacing*delta+delta*delta)/denominator, &
                  2*(2*spacing-delta)*delta/denominator,delta*(delta-spacing)/denominator]
               dweight(1:3)=[(-3*spacing+2*delta)/denominator, &
                  2*(2*spacing-2*delta)/denominator,(2*delta-spacing)/denominator]
               return
            end if
            spacing3=xz%Xs_for_Z(upper,layer)-xz%Xs_for_Z(lower+2,layer)
            if (abs(spacing2-spacing3) <= tiny_fuzz) then
               count=4; indices=[lower,lower+1,lower+2,upper]
               coefficient=(hydrogen-xz%Xs_for_Z(lower+1,layer))/spacing
               coefficient=min(1.0_dp,max(0.0_dp,coefficient))
               weight=[-coefficient*(coefficient-1)**2/2, &
                  (2-coefficient**2*(5-3*coefficient))/2, &
                  coefficient*(1+coefficient*(4-3*coefficient))/2, &
                  coefficient**2*(coefficient-1)/2]
               derivative=1/spacing
               dweight=[-(3*coefficient**2-4*coefficient+1)*derivative/2, &
                  (9*coefficient**2-10*coefficient)*derivative/2, &
                  -(9*coefficient**2-8*coefficient-1)*derivative/2, &
                  coefficient*(3*coefficient-2)*derivative/2]
               return
            end if
         end if
      end if
      do point=2,number
         if (xz%Xs_for_Z(point,layer) >= hydrogen) exit
      end do
      count=2; indices(1:2)=[point-1,point]
      spacing=xz%Xs_for_Z(point,layer)-xz%Xs_for_Z(point-1,layer)
      weight(2)=(hydrogen-xz%Xs_for_Z(point-1,layer))/spacing
      weight(1)=1-weight(2)
      dweight(1:2)=[-1/spacing,1/spacing]
   end subroutine select_x

   subroutine map_family_state(potential,state,by_x,by_z,rho,temperature, &
         values,density_partials,temperature_partials,x_partials,z_partials,status)
      type(potential_partials), intent(in) :: potential
      type(thermodynamic_state), intent(in) :: state
      real(dp), intent(in) :: by_x(10),by_z(10),rho,temperature
      real(dp), intent(out) :: values(nv),density_partials(nv),temperature_partials(nv)
      real(dp), intent(out) :: x_partials(nv),z_partials(nv)
      integer, intent(out) :: status
      integer, parameter :: destinations(7)=[i_grad_ad,i_chiRho,i_chiT,i_Cp,i_Cv,i_gamma1,i_gamma3]
      integer, parameter :: origins(7)=[state_grad_ad,state_chirho,state_chit,state_cp,state_cv, &
         state_gamma1,state_gamma3]
      integer, parameter :: logs(3)=[state_pgas,state_energy,state_entropy]
      integer :: position,slot,source
      values=0; density_partials=0; temperature_partials=0; x_partials=0; z_partials=0
      do slot=i_lnPgas,i_lnS
         source=logs(slot)
         values(slot)=log(state%value(source))
         density_partials(slot)=state%dlnrho(source)/state%value(source)
         temperature_partials(slot)=state%dlnT(source)/state%value(source)
      end do
      do position=1,size(destinations)
         slot=destinations(position); source=origins(position)
         values(slot)=state%value(source)
         density_partials(slot)=state%dlnrho(source)
         temperature_partials(slot)=state%dlnT(source)
      end do
      values(i_dE_dRho)=(potential%f_x-potential%f_xt)/rho
      density_partials(i_dE_dRho)=(potential%f_xx-potential%f_xxt-potential%f_x+potential%f_xt)/rho
      temperature_partials(i_dE_dRho)=(potential%f_xt-potential%f_xtt)/rho
      values(i_dS_dT)=state%value(state_cv)/temperature
      density_partials(i_dS_dT)=state%dlnrho(state_cv)/temperature
      temperature_partials(i_dS_dT)=(state%dlnT(state_cv)-state%value(state_cv))/temperature
      values(i_dS_dRho)=-potential%f_xt/(rho*temperature)
      density_partials(i_dS_dRho)=(potential%f_xt-potential%f_xxt)/(rho*temperature)
      temperature_partials(i_dS_dRho)=(potential%f_xt-potential%f_xtt)/(rho*temperature)
      call composition_partials(by_x,potential,state,rho,temperature,x_partials)
      call composition_partials(by_z,potential,state,rho,temperature,z_partials)
      status=family_invalid_input
      if (.not. all(ieee_is_finite(values(fields))) .or. &
          .not. all(ieee_is_finite(density_partials(fields))) .or. &
          .not. all(ieee_is_finite(temperature_partials(fields))) .or. &
          .not. all(ieee_is_finite(x_partials(fields))) .or. &
          .not. all(ieee_is_finite(z_partials(fields)))) return
      status=table_ok
   end subroutine map_family_state

   subroutine composition_partials(derivative,potential,state,rho,temperature,result)
      real(dp), intent(in) :: derivative(10),rho,temperature
      type(potential_partials), intent(in) :: potential
      type(thermodynamic_state), intent(in) :: state
      real(dp), intent(out) :: result(nv)
      real(dp) :: dcv,dchir,dchit,dgamma3,dgamma1,dcp,dfx,dfxx,dfxt
      real(dp) :: chir,chit,cv,gamma3minus1,excess
      result=0
      dfx=derivative(2); dfxx=derivative(4); dfxt=derivative(5)
      cv=state%value(state_cv); chir=state%value(state_chirho); chit=state%value(state_chit)
      gamma3minus1=state%value(state_gamma3)-1
      dcv=(derivative(3)-derivative(6))/temperature
      dchir=(dfxx-(potential%f_xx/potential%f_x)*dfx)/potential%f_x
      dchit=(dfxt-chit*dfx)/potential%f_x
      dgamma3=dfxt/(temperature*cv)-gamma3minus1*dcv/cv
      dgamma1=dchir+dchit*gamma3minus1+chit*dgamma3
      excess=state%value(state_cp)-cv
      dcp=dcv+excess*(dfx/potential%f_x+2*dchit/chit-dchir/chir)
      result(i_lnPgas)=rho*dfx/state%value(state_pgas)
      result(i_lnE)=(derivative(1)-derivative(3))/state%value(state_energy)
      result(i_lnS)=(-derivative(3)/temperature)/state%value(state_entropy)
      result(i_grad_ad)=(dgamma3-state%value(state_grad_ad)*dgamma1)/state%value(state_gamma1)
      result(i_chiRho)=dchir
      result(i_chiT)=dchit
      result(i_Cp)=dcp
      result(i_Cv)=dcv
      result(i_dE_dRho)=(dfx-dfxt)/rho
      result(i_dS_dT)=dcv/temperature
      result(i_dS_dRho)=-dfxt/(rho*temperature)
      result(i_gamma1)=dgamma1
      result(i_gamma3)=dgamma3
   end subroutine composition_partials
end module freeeos_family
