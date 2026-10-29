module freeeos_table
   use const_def, only: dp, ln10
   use auto_diff_real_2var_order3_module
   use freeeos_biquintic, only: freeeos_partials, evaluate_cell_real, evaluate_cell_ad, cell_ok
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use, intrinsic :: iso_fortran_env, only: int64, iostat_end, iostat_eor
   implicit none
   private
   public :: potential_partials, table_metadata, biquintic_table, thermodynamic_state
   public :: load_table, get_metadata, evaluate_table_real, evaluate_table_ad, derive_state
   public :: partial_vector, ad_vector, read_tokens
   public :: table_ok, table_io_error, table_format_error, table_mismatch, table_unloaded
   public :: table_invalid_input, table_out_of_bounds, table_invalid_cell, table_arithmetic_failure
   public :: table_nonphysical, state_count, state_names
   public :: state_pressure, state_pgas, state_energy, state_entropy, state_cv, state_cp
   public :: state_chirho, state_chit, state_gamma1, state_gamma3, state_grad_ad

   integer, parameter :: table_ok=0, table_io_error=10, table_format_error=11, table_mismatch=12
   integer, parameter :: table_unloaded=13, table_invalid_input=14, table_out_of_bounds=15
   integer, parameter :: table_invalid_cell=16, table_arithmetic_failure=17, table_nonphysical=18
   integer, parameter :: max_axis=100000, max_nodes=2000000, state_count=11
   integer, parameter :: state_pressure=1, state_pgas=2, state_energy=3, state_entropy=4
   integer, parameter :: state_cv=5, state_cp=6, state_chirho=7, state_chit=8
   integer, parameter :: state_gamma1=9, state_gamma3=10, state_grad_ad=11
   character(len=8), parameter :: state_names(state_count) = [character(len=8) :: &
      'P', 'Pgas', 'E', 'S', 'Cv', 'Cp', 'chiRho', 'chiT', 'gamma1', 'gamma3', 'grad_ad']
   real(dp), parameter :: coordinate_offset=12.0_dp*ln10

   type :: potential_partials
      real(dp) :: f, f_x, f_t, f_xx, f_xt, f_tt, f_xxx, f_xxt, f_xtt, f_ttt
   end type potential_partials

   type :: table_metadata
      integer :: nq=0, nt=0, split=0, radius=0
      real(dp) :: hydrogen=0.0_dp, metals=0.0_dp, crad=0.0_dp
      real(dp) :: u_min=0.0_dp, du=0.0_dp, tau_min=0.0_dp, dtau=0.0_dp
      character(len=64) :: source_sha256=''
      character(len=64) :: preprocessing_sha256=''
      character(len=6) :: route=''
   end type table_metadata

   type :: biquintic_table
      private
      type(table_metadata) :: metadata
      real(dp), allocatable :: nodes(:,:,:,:), u_axis(:), tau_axis(:)
      logical, allocatable :: valid(:,:)
   end type biquintic_table

   type :: thermodynamic_state
      real(dp) :: value(state_count), dlnrho(state_count), dlnT(state_count)
   end type thermodynamic_state

contains

   subroutine load_table(path, expected_x, expected_z, table, status, expected_sha)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: expected_x, expected_z
      type(biquintic_table), intent(out) :: table
      integer, intent(out) :: status
      character(len=*), intent(in), optional :: expected_sha
      type(biquintic_table) :: candidate
      type(table_metadata) :: metadata
      character(len=128) :: tokens(12)
      integer :: unit, io_status, token_status, count, index_u, index_tau, record_u, record_tau
      integer :: validity, derivative_u, derivative_tau, allocation_status
      real(dp) :: header(4), node_values(9)

      status = table_invalid_input
      if (.not. all(ieee_is_finite([expected_x, expected_z]))) return
      if (present(expected_sha)) then
         if (.not. valid_sha(expected_sha)) return
      end if
      status = table_io_error
      open(newunit=unit, file=path, status='old', action='read', iostat=io_status)
      if (io_status /= 0) return
      status = table_format_error
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 1) goto 900
      if (tokens(1) /= 'FREEEOS_BIQUINTIC_V1') goto 900
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 2) goto 900
      if (.not. parse_integer(tokens(1), metadata%nq)) goto 900
      if (.not. parse_integer(tokens(2), metadata%nt)) goto 900
      if (min(metadata%nq, metadata%nt) < 2 .or. max(metadata%nq, metadata%nt) > max_axis) goto 900
      if (int(metadata%nq,int64)*int(metadata%nt,int64) > max_nodes) goto 900
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 4) goto 900
      if (.not. parse_real(tokens(1), metadata%hydrogen)) goto 900
      if (.not. parse_real(tokens(2), metadata%metals)) goto 900
      if (.not. parse_real(tokens(3), metadata%crad)) goto 900
      if (.not. parse_integer(tokens(4), metadata%split)) goto 900
      if (metadata%hydrogen < 0 .or. metadata%metals < 0 .or. &
          metadata%hydrogen+metadata%metals > 1 .or. metadata%crad < 0) goto 900
      if (metadata%split /= 0 .and. metadata%split /= 1) goto 900
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 1) goto 900
      if (.not. valid_sha(trim(tokens(1)))) goto 900
      metadata%source_sha256 = lower_hex(tokens(1))
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 3) goto 900
      if (.not. valid_sha(trim(tokens(1)))) goto 900
      metadata%preprocessing_sha256=lower_hex(tokens(1))
      if (.not. parse_integer(tokens(2),metadata%radius)) goto 900
      if (metadata%radius < 1 .or. metadata%radius > 3) goto 900
      if (tokens(3) /= 'mean' .and. tokens(3) /= 'first' .and. tokens(3) /= 'second') goto 900
      metadata%route=trim(tokens(3))
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= 0 .or. count /= 4) goto 900
      do index_u=1,4
         if (.not. parse_real(tokens(index_u), header(index_u))) goto 900
      end do
      metadata%u_min=header(1)
      metadata%du=header(2)
      metadata%tau_min=header(3)
      metadata%dtau=header(4)
      if (metadata%du <= 0 .or. metadata%dtau <= 0) goto 900
      if (abs(metadata%hydrogen-expected_x) > 32*epsilon(1.0_dp) .or. &
          abs(metadata%metals-expected_z) > 32*epsilon(1.0_dp)) then
         status = table_mismatch
         goto 900
      end if
      if (present(expected_sha)) then
         if (metadata%source_sha256 /= lower_hex(expected_sha)) then
            status = table_mismatch
            goto 900
         end if
      end if
      allocate(candidate%u_axis(metadata%nq), candidate%tau_axis(metadata%nt), stat=allocation_status)
      if (allocation_status /= 0) then
         status = table_io_error
         goto 900
      end if
      do index_u=1,metadata%nq
         candidate%u_axis(index_u)=metadata%u_min+real(index_u-1,dp)*metadata%du
      end do
      do index_tau=1,metadata%nt
         candidate%tau_axis(index_tau)=metadata%tau_min+real(index_tau-1,dp)*metadata%dtau
      end do
      if (.not. axis_valid(candidate%u_axis) .or. .not. axis_valid(candidate%tau_axis)) goto 900
      allocate(candidate%nodes(0:2,0:2,metadata%nq,metadata%nt), &
         candidate%valid(metadata%nq,metadata%nt), stat=allocation_status)
      if (allocation_status /= 0) then
         status = table_io_error
         goto 900
      end if
      do index_tau=1,metadata%nt
         do index_u=1,metadata%nq
            call read_tokens(unit, tokens, count, token_status)
            if (token_status /= 0 .or. count /= 12) goto 900
            if (.not. parse_integer(tokens(1), record_u)) goto 900
            if (.not. parse_integer(tokens(2), record_tau)) goto 900
            if (.not. parse_integer(tokens(3), validity)) goto 900
            if (record_u /= index_u .or. record_tau /= index_tau) goto 900
            if (validity /= 0 .and. validity /= 1) goto 900
            do derivative_u=1,9
               if (.not. parse_real(tokens(derivative_u+3), node_values(derivative_u))) goto 900
            end do
            if (validity == 0 .and. any(node_values /= 0.0_dp)) goto 900
            candidate%valid(index_u,index_tau)=validity == 1
            do derivative_tau=0,2
               do derivative_u=0,2
                  candidate%nodes(derivative_u,derivative_tau,index_u,index_tau) = &
                     node_values(1+derivative_u+3*derivative_tau)
               end do
            end do
         end do
      end do
      call read_tokens(unit, tokens, count, token_status)
      if (token_status /= iostat_end) goto 900
      candidate%metadata=metadata
      call move_alloc(candidate%nodes, table%nodes)
      call move_alloc(candidate%valid, table%valid)
      call move_alloc(candidate%u_axis, table%u_axis)
      call move_alloc(candidate%tau_axis, table%tau_axis)
      table%metadata=metadata
      status=table_ok
900   close(unit)
   end subroutine load_table

   function get_metadata(table) result(metadata)
      type(biquintic_table), intent(in) :: table
      type(table_metadata) :: metadata
      metadata=table%metadata
   end function get_metadata

   logical function axis_valid(axis)
      real(dp), intent(in) :: axis(:)
      integer :: position
      real(dp) :: spacing
      axis_valid=.false.
      if (.not. all(ieee_is_finite(axis))) return
      do position=1,size(axis)-1
         spacing=axis(position+1)-axis(position)
         if (.not. ieee_is_finite(spacing)) return
         if (spacing < 8.0_dp*huge(1.0_dp)**(-1.0_dp/3.0_dp) .or. &
             spacing > 0.5_dp*sqrt(huge(1.0_dp))) return
         if (spacing <= &
             8*epsilon(1.0_dp)*max(abs(axis(position)),abs(axis(position+1)),tiny(1.0_dp))) return
      end do
      axis_valid=.true.
   end function axis_valid

   subroutine read_tokens(unit, tokens, count, status)
      integer, intent(in) :: unit
      character(len=*), intent(out) :: tokens(:)
      integer, intent(out) :: count, status
      character(len=4096) :: record
      integer :: length, position, start
      count=0
      tokens=''
      read(unit,'(a)',advance='no',size=length,iostat=status) record
      if (status == iostat_end .and. length == 0) return
      if (status /= iostat_eor) then
         status=table_format_error
         return
      end if
      status=table_format_error
      position=1
      do while (position <= length)
         if (index(' '//achar(9)//achar(13),record(position:position)) /= 0) then
            position=position+1
            cycle
         end if
         start=position
         do while (position <= length)
            if (index(' '//achar(9)//achar(13),record(position:position)) /= 0) exit
            position=position+1
         end do
         count=count+1
         if (count > size(tokens)) return
         if (position-start > len(tokens)) return
         tokens(count)=record(start:position-1)
      end do
      status=0
   end subroutine read_tokens

   logical function parse_integer(token, value)
      character(len=*), intent(in) :: token
      integer, intent(out) :: value
      integer :: status
      parse_integer=.false.
      if (len_trim(token) == 0) return
      if (verify(trim(token),'0123456789+-') /= 0) return
      read(token,*,iostat=status) value
      parse_integer=status == 0
   end function parse_integer

   logical function parse_real(token, value)
      character(len=*), intent(in) :: token
      real(dp), intent(out) :: value
      integer :: status
      parse_real=.false.
      if (len_trim(token) == 0) return
      if (verify(trim(token),'0123456789+-.eEdD') /= 0) return
      read(token,*,iostat=status) value
      if (status /= 0) return
      parse_real=ieee_is_finite(value)
   end function parse_real

   logical function valid_sha(sha)
      character(len=*), intent(in) :: sha
      valid_sha=len_trim(sha) == 64 .and. verify(trim(sha),'0123456789abcdefABCDEF') == 0
   end function valid_sha

   function lower_hex(sha) result(lower)
      character(len=*), intent(in) :: sha
      character(len=64) :: lower
      integer :: position, code
      lower=sha
      do position=1,64
         code=iachar(lower(position:position))
         if (code >= iachar('A') .and. code <= iachar('F')) lower(position:position)=achar(code+32)
      end do
   end function lower_hex

   subroutine locate_cell(table, x, t, u, index_u, index_tau, status)
      type(biquintic_table), intent(in) :: table
      real(dp), intent(in) :: x, t
      real(dp), intent(out) :: u
      integer, intent(out) :: index_u, index_tau, status
      real(dp) :: endpoint_tolerance
      status=table_unloaded
      if (.not. allocated(table%nodes)) return
      status=table_invalid_input
      if (.not. all(ieee_is_finite([x,t]))) return
      u=x-2*t+coordinate_offset
      if (.not. ieee_is_finite(u)) return
      status=table_out_of_bounds
      endpoint_tolerance=8*epsilon(1.0_dp)*max(1.0_dp,abs(x),abs(2*t),coordinate_offset, &
         abs(table%u_axis(1)),abs(table%u_axis(table%metadata%nq)))
      if (u < table%u_axis(1)) then
         if (table%u_axis(1)-u > endpoint_tolerance) return
         u=table%u_axis(1)
      else if (u > table%u_axis(table%metadata%nq)) then
         if (u-table%u_axis(table%metadata%nq) > endpoint_tolerance) return
         u=table%u_axis(table%metadata%nq)
      end if
      if (t < table%tau_axis(1) .or. t > table%tau_axis(table%metadata%nt)) return
      index_u=lower_cell(table%u_axis,u)
      index_tau=lower_cell(table%tau_axis,t)
      status=table_invalid_cell
      if (.not. all(table%valid(index_u:index_u+1,index_tau:index_tau+1))) return
      status=table_ok
   end subroutine locate_cell

   integer function lower_cell(axis, coordinate) result(cell)
      real(dp), intent(in) :: axis(:), coordinate
      integer :: lower, upper, middle
      lower=1
      upper=size(axis)
      do while (upper-lower > 1)
         middle=lower+(upper-lower)/2
         if (coordinate >= axis(middle)) then
            lower=middle
         else
            upper=middle
         end if
      end do
      cell=lower
   end function lower_cell

   subroutine evaluate_table_real(table, x, t, partials, status)
      type(biquintic_table), intent(in) :: table
      real(dp), intent(in) :: x, t
      type(potential_partials), intent(out) :: partials
      integer, intent(out) :: status
      type(freeeos_partials) :: native
      real(dp) :: u, radiation, values(10), nan, nodal(0:2,0:2,2,2)
      integer :: index_u, index_tau
      nan=ieee_value(0.0_dp,ieee_quiet_nan)
      partials=potential_partials(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan)
      call locate_cell(table,x,t,u,index_u,index_tau,status)
      if (status /= table_ok) return
      nodal=table%nodes(:,:,index_u:index_u+1,index_tau:index_tau+1)
      call evaluate_cell_real(nodal, &
         table%u_axis(index_u),table%u_axis(index_u+1),table%tau_axis(index_tau), &
         table%tau_axis(index_tau+1),u,t,native,status)
      if (status /= cell_ok) then
         status=table_arithmetic_failure
         return
      end if
      values=[native%f, native%f_u, native%f_tau-2*native%f_u, native%f_uu, &
         native%f_u_tau-2*native%f_uu, native%f_tau_tau-4*native%f_u_tau+4*native%f_uu, &
         native%f_uuu, native%f_uu_tau-2*native%f_uuu, &
         native%f_u_tau_tau-4*native%f_uu_tau+4*native%f_uuu, &
         native%f_tau_tau_tau-6*native%f_u_tau_tau+12*native%f_uu_tau-8*native%f_uuu]
      if (table%metadata%split == 1 .and. table%metadata%crad > 0) then
         radiation=-table%metadata%crad*exp(4*t-x)/3
         values=values+radiation*[1.0_dp,-1.0_dp,4.0_dp,1.0_dp,-4.0_dp,16.0_dp, &
            -1.0_dp,4.0_dp,-16.0_dp,64.0_dp]
      end if
      if (.not. all(ieee_is_finite(values))) then
         status=table_arithmetic_failure
         return
      end if
      partials=potential_partials(values(1),values(2),values(3),values(4),values(5), &
         values(6),values(7),values(8),values(9),values(10))
   end subroutine evaluate_table_real

   subroutine evaluate_table_ad(table, x, t, potential, status)
      type(biquintic_table), intent(in) :: table
      type(auto_diff_real_2var_order3), intent(in) :: x, t
      type(auto_diff_real_2var_order3), intent(out) :: potential
      integer, intent(out) :: status
      type(auto_diff_real_2var_order3) :: u_ad, candidate
      real(dp) :: u, nan, nodal(0:2,0:2,2,2)
      integer :: index_u, index_tau
      nan=ieee_value(0.0_dp,ieee_quiet_nan)
      potential=auto_diff_real_2var_order3(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan)
      call locate_cell(table,x%val,t%val,u,index_u,index_tau,status)
      if (status /= table_ok) return
      status=table_invalid_input
      if (.not. all(ieee_is_finite(ad_vector(x))) .or. .not. all(ieee_is_finite(ad_vector(t)))) return
      u_ad=x-2.0_dp*t+coordinate_offset
      u_ad%val=u
      nodal=table%nodes(:,:,index_u:index_u+1,index_tau:index_tau+1)
      call evaluate_cell_ad(nodal, &
         table%u_axis(index_u),table%u_axis(index_u+1),table%tau_axis(index_tau), &
         table%tau_axis(index_tau+1),u_ad,t,candidate,status)
      if (status /= cell_ok) then
         status=table_arithmetic_failure
         return
      end if
      if (table%metadata%split == 1 .and. table%metadata%crad > 0) &
         candidate=candidate-(table%metadata%crad/3.0_dp)*exp(4.0_dp*t-x)
      if (.not. all(ieee_is_finite(ad_vector(candidate)))) then
         status=table_arithmetic_failure
         return
      end if
      potential=candidate
   end subroutine evaluate_table_ad

   function partial_vector(partials) result(values)
      type(potential_partials), intent(in) :: partials
      real(dp) :: values(10)
      values=[partials%f,partials%f_x,partials%f_t,partials%f_xx,partials%f_xt,partials%f_tt, &
         partials%f_xxx,partials%f_xxt,partials%f_xtt,partials%f_ttt]
   end function partial_vector

   function ad_vector(potential) result(values)
      type(auto_diff_real_2var_order3), intent(in) :: potential
      real(dp) :: values(10)
      values=[potential%val,potential%d1val2,potential%d1val1,potential%d2val2, &
         potential%d1val1_d1val2,potential%d2val1,potential%d3val2, &
         potential%d1val1_d2val2,potential%d2val1_d1val2,potential%d3val1]
   end function ad_vector

   subroutine derive_state(partials, x, t, crad, state, status)
      type(potential_partials), intent(in) :: partials
      real(dp), intent(in) :: x, t, crad
      type(thermodynamic_state), intent(out) :: state
      integer, intent(out) :: status
      type(thermodynamic_state) :: candidate
      real(dp) :: density, temperature, radiation_pressure, nan
      real(dp) :: values(state_count), derivatives(state_count,2), excess, gamma3_minus1
      real(dp) :: fx, ft, fxx, fxt, ftt, fxxx, fxxt, fxtt, fttt
      integer :: axis
      nan=ieee_value(0.0_dp,ieee_quiet_nan)
      state%value=nan
      state%dlnrho=nan
      state%dlnT=nan
      status=table_invalid_input
      if (.not. all(ieee_is_finite([x,t,crad,partial_vector(partials)]))) return
      if (crad < 0) return
      density=exp(x)
      temperature=exp(t)
      radiation_pressure=0
      if (crad > 0) radiation_pressure=crad*exp(4*t)/3
      status=table_arithmetic_failure
      if (.not. all(ieee_is_finite([density,temperature,radiation_pressure]))) return
      if (density <= 0 .or. temperature <= 0) return
      fx=partials%f_x
      ft=partials%f_t
      fxx=partials%f_xx
      fxt=partials%f_xt
      ftt=partials%f_tt
      fxxx=partials%f_xxx
      fxxt=partials%f_xxt
      fxtt=partials%f_xtt
      fttt=partials%f_ttt
      values(state_pressure)=density*fx
      values(state_pgas)=values(state_pressure)-radiation_pressure
      values(state_energy)=partials%f-ft
      values(state_entropy)=-ft/temperature
      values(state_cv)=(ft-ftt)/temperature
      status=table_nonphysical
      if (fx <= 0 .or. any(values(1:5) <= 0)) return
      values(state_chirho)=1+fxx/fx
      values(state_chit)=fxt/fx
      if (values(state_chirho) <= 0 .or. values(state_chit) <= 0) return
      gamma3_minus1=fxt/(temperature*values(state_cv))
      values(state_gamma3)=1+gamma3_minus1
      values(state_gamma1)=values(state_chirho)+values(state_chit)*gamma3_minus1
      excess=fx*values(state_chit)**2/(temperature*values(state_chirho))
      values(state_cp)=values(state_cv)+excess
      if (.not. all(ieee_is_finite(values(1:10)))) then
         status=table_arithmetic_failure
         return
      end if
      if (any(values(1:10) <= 0) .or. gamma3_minus1 <= 0) return
      values(state_grad_ad)=gamma3_minus1/values(state_gamma1)
      derivatives(state_pressure,:)=[density*(fx+fxx),density*fxt]
      derivatives(state_pgas,:)=derivatives(state_pressure,:)-[0.0_dp,4*radiation_pressure]
      derivatives(state_energy,:)=[fx-fxt,ft-ftt]
      derivatives(state_entropy,:)=[-fxt/temperature,(ft-ftt)/temperature]
      derivatives(state_cv,:)=[(fxt-fxtt)/temperature,(2*ftt-fttt-ft)/temperature]
      derivatives(state_chirho,:)=[(fxxx-(fxx/fx)*fxx)/fx,(fxxt-(fxx/fx)*fxt)/fx]
      derivatives(state_chit,:)=[(fxxt-(fxt/fx)*fxx)/fx,(fxtt-(fxt/fx)*fxt)/fx]
      do axis=1,2
         derivatives(state_cp,axis)=derivatives(state_cv,axis)+excess*( &
            merge(fxx,fxt,axis == 1)/fx+2*derivatives(state_chit,axis)/values(state_chit) &
            -derivatives(state_chirho,axis)/values(state_chirho)-merge(0.0_dp,1.0_dp,axis == 1))
         derivatives(state_gamma3,axis)=merge(fxxt,fxtt,axis == 1)/(temperature*values(state_cv)) &
            -gamma3_minus1*(derivatives(state_cv,axis)/values(state_cv)+merge(0.0_dp,1.0_dp,axis == 1))
         derivatives(state_gamma1,axis)=derivatives(state_chirho,axis) &
            +derivatives(state_chit,axis)*gamma3_minus1+values(state_chit)*derivatives(state_gamma3,axis)
         derivatives(state_grad_ad,axis)=(derivatives(state_gamma3,axis) &
            -values(state_grad_ad)*derivatives(state_gamma1,axis))/values(state_gamma1)
      end do
      status=table_arithmetic_failure
      if (.not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(derivatives))) return
      status=table_nonphysical
      if (any(values <= 0)) return
      candidate%value=values
      candidate%dlnrho=derivatives(:,1)
      candidate%dlnT=derivatives(:,2)
      state=candidate
      status=table_ok
   end subroutine derive_state

end module freeeos_table
