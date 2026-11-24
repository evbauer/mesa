module freeeos_legacy_nodes
   use const_def, only: dp, ln10
   use math_lib, only: exp
   use eos_def, only: EosDT_XZ_Info, sz_per_eos_point, nv, i_lnPgas, i_lnE, i_lnS, i_Cv, i_chiRho, i_chiT
   use freeeos_table, only: table_metadata, table_ok, table_io_error, table_mismatch, table_nonphysical
   use ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: reconstruct_legacy_nodes
contains
   subroutine reconstruct_legacy_nodes(ep,metadata,valid,nodes,status)
      ! Slot 1 contains nodal physics; slots 2:4 contain legacy spline coefficients.
      ! Flags and density-coordinate source audits are encoded in the distributed
      ! complete-jet mask, since they are not columns of the legacy binary table.
      type(EosDT_XZ_Info), pointer, intent(in) :: ep
      type(table_metadata), intent(in) :: metadata
      logical, intent(in) :: valid(:,:)
      real(dp), allocatable, intent(out) :: nodes(:,:,:,:)
      integer, intent(out) :: status
      real(dp), pointer :: values(:,:,:,:)
      real(dp) :: logq,logt,temperature,density,gas_pressure,energy,entropy,cv,chirho,chit
      real(dp) :: radiation_pressure,pressure,radiation_scale,f,fx,ft,fxx,fxt,ftt
      integer :: iq,it,nq,allocation_status
      status=table_mismatch
      if (.not. associated(ep)) return
      if (.not. associated(ep%tbl1)) return
      if (.not. associated(ep%logQs) .or. .not. associated(ep%logTs)) return
      if (ep%version /= 51) return
      if (ep%num_logTs /= metadata%nt) return
      nq=min(ep%num_logQs,metadata%nq)
      if (nq < 2) return
      if (any(shape(valid) /= [metadata%nq,metadata%nt])) return
      if (size(ep%tbl1) /= sz_per_eos_point*nv*ep%num_logQs*ep%num_logTs) return
      if (.not. all(ieee_is_finite([ep%logQ_min,ep%logQ_max,ep%del_logQ, &
         ep%logT_min,ep%logT_max,ep%del_logT]))) return
      if (.not. coordinate_matches(ep%logQ_min,metadata%u_min/ln10)) return
      if (.not. coordinate_matches(ep%del_logQ,metadata%du/ln10)) return
      if (.not. coordinate_matches(ep%logT_min,metadata%tau_min/ln10)) return
      if (.not. coordinate_matches(ep%del_logT,metadata%dtau/ln10)) return
      ! Header endpoint arithmetic has decimal rounding; tolerance is in dex,
      ! matching the existing native grid audit, not an EOS-return tolerance.
      if (abs(ep%logT_max-(ep%logT_min+(metadata%nt-1)*ep%del_logT)) > 1e-10_dp) return
      if (ep%del_logQ <= 0 .or. ep%del_logT <= 0) return
      if (size(ep%logQs) /= ep%num_logQs .or. size(ep%logTs) /= metadata%nt) return
      if (.not. all(ieee_is_finite(ep%logQs)) .or. .not. all(ieee_is_finite(ep%logTs))) return
      ! Locate_logQ uses the loaded count, origin and spacing, not logQ_max.
      ! The loader stores logQ_max as its last abscissa even when that header
      ! field is wrong. Leave it unchanged and validate the interior grid.
      ! A shorter loaded table uses only the aligned prefix of the supplement.
      do iq=1,min(nq,ep%num_logQs-1)
         if (abs(ep%logQs(iq)-(ep%logQ_min+(iq-1)*ep%del_logQ)) > 1e-10_dp) return
      end do
      do it=1,metadata%nt
         if (abs(ep%logTs(it)-(ep%logT_min+(it-1)*ep%del_logT)) > 1e-10_dp) return
      end do
      values(1:sz_per_eos_point,1:nv,1:ep%num_logQs,1:metadata%nt) => ep%tbl1
      allocate(nodes(0:2,0:2,nq,metadata%nt),stat=allocation_status)
      status=table_io_error
      if (allocation_status /= 0) return
      nodes=0
      status=table_nonphysical
      do it=1,metadata%nt
         logt=ep%logT_min+(it-1)*ep%del_logT
         temperature=10.0_dp**logt
         do iq=1,nq
            if (.not. valid(iq,it)) cycle
            if (logt < 3.2_dp-1e-12_dp) return
            if (.not. all(ieee_is_finite(values(1,[i_lnPgas,i_lnE,i_lnS,i_Cv,i_chiRho,i_chiT],iq,it)))) return
            logq=ep%logQ_min+(iq-1)*ep%del_logQ
            density=10.0_dp**(logq+2*logt-12)
            gas_pressure=exp(values(1,i_lnPgas,iq,it))
            energy=exp(values(1,i_lnE,iq,it))
            entropy=exp(values(1,i_lnS,iq,it))
            cv=values(1,i_Cv,iq,it)
            chirho=values(1,i_chiRho,iq,it)
            chit=values(1,i_chiT,iq,it)
            radiation_pressure=metadata%crad*temperature**4/3
            pressure=gas_pressure+radiation_pressure
            radiation_scale=radiation_pressure/density
            if (.not. all(ieee_is_finite([density,temperature,pressure,energy,entropy,cv,chirho,chit]))) return
            if (min(density,temperature,pressure,energy,entropy,cv,chirho) <= 0) return
            f=energy-temperature*entropy
            ft=-temperature*entropy
            ftt=-temperature*(entropy+cv)
            f=f+radiation_scale
            fx=gas_pressure/density
            ft=ft+4*radiation_scale
            fxx=(pressure*chirho-gas_pressure)/density
            fxt=(pressure*chit-4*radiation_pressure)/density
            ftt=ftt+16*radiation_scale
            nodes(0,0,iq,it)=f
            nodes(1,0,iq,it)=fx
            nodes(0,1,iq,it)=ft+2*fx
            nodes(2,0,iq,it)=fxx
            nodes(1,1,iq,it)=fxt+2*fxx
            nodes(0,2,iq,it)=ftt+4*fxt+4*fxx
            if (.not. all(ieee_is_finite(nodes(:,:,iq,it)))) return
         end do
      end do
      status=table_ok
   end subroutine reconstruct_legacy_nodes

   ! The original source/header audit allows 1e-10 dex. In particular the
   ! printed row T origin differs from its header by about 9e-14 dex.
   logical function coordinate_matches(a,b)
      real(dp), intent(in) :: a,b
      coordinate_matches=abs(a-b) <= 1e-10_dp
   end function coordinate_matches
end module freeeos_legacy_nodes
