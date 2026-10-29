program test_table
   use const_def, only: dp, ln10
   use freeeos_table
   use auto_diff_real_2var_order3_module
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp), parameter :: gas_constant=2.0e7_dp, radiation_constant=1.0e-6_dp
   character(len=64), parameter :: source_sha=repeat('a',64), processing_sha=repeat('b',64)
   type(biquintic_table) :: table, empty_table
   type(table_metadata) :: metadata
   type(potential_partials) :: partials, shifted_partials
   type(thermodynamic_state) :: state, lower_state, upper_state
   type(auto_diff_real_2var_order3) :: x_ad, t_ad, potential_ad, delta_x, delta_t, taylor
   character(len=4096) :: directory, path
   real(dp) :: u, t, x, values(10), expected(10), steps(2), shifted_x, shifted_t
   real(dp) :: density, temperature, pressure_rad, cv, gamma3_minus1, expected_state(state_count), radiation
   real(dp) :: nan, started, finished, accumulator
   integer :: status, failures, probe, axis, split, benchmark_index, evaluations

   failures=0
   evaluations=0
   nan=ieee_value(0.0_dp,ieee_quiet_nan)
   call get_command_argument(1,directory)
   if (len_trim(directory) == 0) directory='build'
   do split=0,1
      if (split == 1) then
         path=trim(directory)//'/ideal_rad.cache'
      else
         path=trim(directory)//'/total_rad.cache'
      end if
      call write_fixture(path,split,radiation_constant,.false.)
      call load_table(trim(path),0.7_dp,0.02_dp,table,status,source_sha)
      call require(status == table_ok,'load analytic fixture')
      metadata=get_metadata(table)
      call require(metadata%radius == 2 .and. metadata%route == 'mean','processing metadata')
      call require(metadata%preprocessing_sha256 == processing_sha,'processing SHA')
      do probe=1,5
         u=0.14_dp*real(probe,dp)
         t=5.0_dp+0.015_dp*real(probe,dp)
         x=u+2*t-12*ln10
         call evaluate_table_real(table,x,t,partials,status)
         call require(status == table_ok,'real fixture evaluation')
         values=partial_vector(partials)
         x_ad=x
         x_ad%d1val2=1.0_dp
         t_ad=t
         t_ad%d1val1=1.0_dp
         call evaluate_table_ad(table,x_ad,t_ad,potential_ad,status)
         call require(status == table_ok,'AD fixture evaluation')
         call close_vector(values,ad_vector(potential_ad),3.0e-7_dp,'independent real/AD')
         temperature=exp(t)
         radiation=-radiation_constant*exp(4*t-x)/3
         expected=gas_constant*temperature*[x-1.5_dp*t-5.0_dp,1.0_dp,x-1.5_dp*t-6.5_dp, &
            0.0_dp,1.0_dp,x-1.5_dp*t-8.0_dp,0.0_dp,0.0_dp,1.0_dp,x-1.5_dp*t-9.5_dp] &
            +radiation*[1.0_dp,-1.0_dp,4.0_dp,1.0_dp,-4.0_dp,16.0_dp,-1.0_dp,4.0_dp,-16.0_dp,64.0_dp]
         call close_vector(values,expected,merge(2.0e-4_dp,1.0e-3_dp,split == 1),'analytic physical partials')
         call derive_state(partials,x,t,radiation_constant,state,status)
         call require(status == table_ok,'derive positive fixture state')
         density=exp(x)
         pressure_rad=radiation_constant*exp(4*t)/3
         cv=1.5_dp*gas_constant+12*pressure_rad/(density*temperature)
         expected_state(state_pressure)=density*gas_constant*temperature+pressure_rad
         expected_state(state_pgas)=density*gas_constant*temperature
         expected_state(state_energy)=1.5_dp*gas_constant*temperature+3*pressure_rad/density
         expected_state(state_entropy)=gas_constant*(6.5_dp+1.5_dp*t-x)+4*pressure_rad/(density*temperature)
         expected_state(state_cv)=cv
         expected_state(state_chirho)=expected_state(state_pgas)/expected_state(state_pressure)
         expected_state(state_chit)=(expected_state(state_pgas)+4*pressure_rad)/expected_state(state_pressure)
         gamma3_minus1=(expected_state(state_pgas)+4*pressure_rad)/(density*temperature*cv)
         expected_state(state_gamma3)=1+gamma3_minus1
         expected_state(state_gamma1)=expected_state(state_chirho)+expected_state(state_chit)*gamma3_minus1
         expected_state(state_cp)=cv+expected_state(state_pressure)*expected_state(state_chit)**2 &
            /(density*temperature*expected_state(state_chirho))
         expected_state(state_grad_ad)=gamma3_minus1/expected_state(state_gamma1)
         call close_vector(state%value,expected_state,3.0e-5_dp,'analytic gas plus radiation state')
         steps=[2.0e-4_dp,2.0e-4_dp]
         do axis=1,2
            shifted_x=x-merge(steps(axis),0.0_dp,axis == 1)
            shifted_t=t-merge(0.0_dp,steps(axis),axis == 1)
            call evaluate_table_real(table,shifted_x,shifted_t,shifted_partials,status)
            call require(status == table_ok,'lower difference query')
            call derive_state(shifted_partials,shifted_x,shifted_t,radiation_constant,lower_state,status)
            call require(status == table_ok,'lower difference state')
            shifted_x=x+merge(steps(axis),0.0_dp,axis == 1)
            shifted_t=t+merge(0.0_dp,steps(axis),axis == 1)
            call evaluate_table_real(table,shifted_x,shifted_t,shifted_partials,status)
            call require(status == table_ok,'upper difference query')
            call derive_state(shifted_partials,shifted_x,shifted_t,radiation_constant,upper_state,status)
            call require(status == table_ok,'upper difference state')
            if (axis == 1) then
               call close_vector(state%dlnrho,(upper_state%value-lower_state%value)/(2*steps(axis)), &
                  8.0e-5_dp,'state log-density numerical derivatives')
            else
               call close_vector(state%dlnT,(upper_state%value-lower_state%value)/(2*steps(axis)), &
                  8.0e-5_dp,'state log-temperature numerical derivatives')
            end if
         end do
         evaluations=evaluations+1
      end do
   end do
   x_ad=x
   x_ad%d1val2=0.7_dp
   x_ad%d1val1=0.2_dp
   x_ad%d2val2=-0.3_dp
   x_ad%d1val1_d1val2=0.4_dp
   x_ad%d3val1=0.1_dp
   t_ad=t
   t_ad%d1val2=-0.1_dp
   t_ad%d1val1=0.8_dp
   t_ad%d2val1=0.2_dp
   t_ad%d1val1_d2val2=-0.2_dp
   delta_x=x_ad-x
   delta_t=t_ad-t
   taylor=partials%f+partials%f_x*delta_x+partials%f_t*delta_t &
      +0.5_dp*partials%f_xx*delta_x*delta_x+partials%f_xt*delta_x*delta_t &
      +0.5_dp*partials%f_tt*delta_t*delta_t &
      +(partials%f_xxx/6.0_dp)*delta_x*delta_x*delta_x &
      +0.5_dp*partials%f_xxt*delta_x*delta_x*delta_t &
      +0.5_dp*partials%f_xtt*delta_x*delta_t*delta_t &
      +(partials%f_ttt/6.0_dp)*delta_t*delta_t*delta_t
   call evaluate_table_ad(table,x_ad,t_ad,potential_ad,status)
   call require(status == table_ok,'nonlinear physical AD seeds')
   call close_vector(ad_vector(potential_ad),ad_vector(taylor),3.0e-7_dp,'nonlinear AD chain rules through three')
   do probe=0,2
      u=0.4_dp*real(probe,dp)
      do axis=0,1
         t=5.0_dp+0.1_dp*real(axis,dp)
         x=u+2*t-12*ln10
         call evaluate_table_real(table,x,t,partials,status)
         call require(status == table_ok,'grid boundary and two-cell face')
      end do
   end do
   call evaluate_table_real(table,-100.0_dp,5.05_dp,partials,status)
   call require(status == table_out_of_bounds,'reject out of bounds')
   call require(all(ieee_is_nan(partial_vector(partials))),'poison rejected real output')
   call evaluate_table_real(table,nan,5.05_dp,partials,status)
   call require(status == table_invalid_input,'reject NaN query')
   call evaluate_table_real(empty_table,x,t,partials,status)
   call require(status == table_unloaded,'reject unloaded table')
   x_ad=x
   x_ad%d3val2=ieee_value(0.0_dp,ieee_positive_inf)
   t_ad=t
   call evaluate_table_ad(table,x_ad,t_ad,potential_ad,status)
   call require(status == table_invalid_input,'reject nonfinite AD seed')
   call require(all(ieee_is_nan(ad_vector(potential_ad))),'poison rejected AD output')
   call load_table(trim(path),0.6_dp,0.02_dp,table,status)
   call require(status == table_mismatch,'reject composition mismatch')
   call evaluate_table_real(table,x,t,partials,status)
   call require(status == table_unloaded,'failed reload clears old table')
   call load_table(trim(path),0.7_dp,0.02_dp,table,status,processing_sha)
   call require(status == table_mismatch,'reject source SHA mismatch')
   path=trim(directory)//'/invalid_node.cache'
   call write_fixture(path,1,radiation_constant,.true.)
   call load_table(trim(path),0.7_dp,0.02_dp,table,status)
   call require(status == table_ok,'load invalid-node mask')
   t=5.05_dp
   x=0.2_dp+2*t-12*ln10
   call evaluate_table_real(table,x,t,partials,status)
   call require(status == table_invalid_cell,'first cell requires all four corners')
   x=0.6_dp+2*t-12*ln10
   call evaluate_table_real(table,x,t,partials,status)
   call require(status == table_ok,'second cell remains valid')
   path=trim(directory)//'/ideal_gas.cache'
   call write_fixture(path,1,0.0_dp,.false.)
   call load_table(trim(path),0.7_dp,0.02_dp,table,status)
   call require(status == table_ok,'load ideal gas without radiation')
   call evaluate_table_real(table,x,t,partials,status)
   call derive_state(partials,x,t,0.0_dp,state,status)
   call require(status == table_ok,'ideal gas state')
   call close_vector(state%value(7:11),[1.0_dp,1.0_dp,5.0_dp/3,5.0_dp/3,0.4_dp], &
      2.0e-6_dp,'ideal gas gamma3 is one plus minus1')
   partials%f_x=-1.0_dp
   call derive_state(partials,x,t,0.0_dp,state,status)
   call require(status == table_nonphysical,'reject negative required pressure')
   call require(all(ieee_is_nan(state%value)),'poison rejected state')
   partials=potential_partials(1.0_dp,1.0_dp,-1.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp)
   call derive_state(partials,0.0_dp,0.0_dp,0.0_dp,state,status)
   call require(status == table_nonphysical,'reject negative Cv without clipping')
   accumulator=0
   call cpu_time(started)
   do benchmark_index=1,100000
      call evaluate_table_real(table,x,t,partials,status)
      accumulator=accumulator+partials%f
   end do
   call cpu_time(finished)
   write(*,'(a,f10.5,a,es12.4)') '100000 real queries, no I/O: ',finished-started,' seconds; sum=',accumulator
   write(*,'(a,i0,a,i0)') 'Analytic real/AD pairs: ',evaluations,'; failures: ',failures
   if (failures /= 0) stop 1
   write(*,'(a)') 'PASS standalone table analytic, state derivative, AD, two-cell and rejection tests'

contains

   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (condition) return
      failures=failures+1
      write(*,'(a)') 'FAIL '//label
   end subroutine require

   subroutine close_vector(measured,reference,tolerance,label)
      real(dp), intent(in) :: measured(:), reference(:), tolerance
      character(len=*), intent(in) :: label
      real(dp) :: error
      error=maxval(abs(measured-reference)/max(1.0_dp,abs(reference)))
      if (error <= tolerance) return
      failures=failures+1
      write(*,'(a,a,a,es12.4,a,es12.4)') 'FAIL ',label,': ',error,' > ',tolerance
      write(*,'(a,i0,a,es16.8,a,es16.8)') 'component=', &
         maxloc(abs(measured-reference)/max(1.0_dp,abs(reference)),dim=1), &
         '; measured=',measured(maxloc(abs(measured-reference)/max(1.0_dp,abs(reference)),dim=1)), &
         '; reference=',reference(maxloc(abs(measured-reference)/max(1.0_dp,abs(reference)),dim=1))
   end subroutine close_vector

   subroutine write_fixture(filename,split_flag,crad,invalid_first)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: split_flag
      real(dp), intent(in) :: crad
      logical, intent(in) :: invalid_first
      real(dp) :: nodal(0:2,0:2), native_u, native_t, factor, bracket, rad
      integer :: unit, index_u, index_tau, derivative_u, derivative_tau, validity
      open(newunit=unit,file=filename,status='replace',action='write')
      write(unit,'(a)') 'FREEEOS_BIQUINTIC_V1'
      write(unit,'(a)') '3 2'
      write(unit,'(3(es25.16e3,1x),i0)') 0.7_dp,0.02_dp,crad,split_flag
      write(unit,'(a)') source_sha
      write(unit,'(a)') processing_sha//' 2 mean'
      write(unit,'(4(es25.16e3,1x))') 0.0_dp,0.4_dp,5.0_dp,0.1_dp
      do index_tau=1,2
         native_t=5.0_dp+real(index_tau-1,dp)*0.1_dp
         do index_u=1,3
            native_u=real(index_u-1,dp)*0.4_dp
            factor=gas_constant*exp(native_t)
            bracket=native_u+0.5_dp*native_t-12*ln10-5.0_dp
            nodal=0
            nodal(0,:)=[factor*bracket,factor*(bracket+0.5_dp),factor*(bracket+1.0_dp)]
            nodal(1,:)=factor
            if (split_flag == 0 .and. crad > 0) then
               rad=-crad*exp(2*native_t-native_u+12*ln10)/3
               do derivative_tau=0,2
                  do derivative_u=0,2
                     nodal(derivative_u,derivative_tau)=nodal(derivative_u,derivative_tau) &
                        +rad*real((-1)**derivative_u*2**derivative_tau,dp)
                  end do
               end do
            end if
            validity=1
            if (invalid_first .and. index_u == 1 .and. index_tau == 1) then
               validity=0
               nodal=0
            end if
            write(unit,'(3(i0,1x),9(es25.16e3,1x))') index_u,index_tau,validity,nodal
         end do
      end do
      close(unit)
   end subroutine write_fixture

end program test_table
