program jet_driver
   use const_def, only: dp
   use auto_diff_real_2var_order3_module
   use freeeos_biquintic, only: evaluate_cell_ad
   use iso_fortran_env, only: input_unit, output_unit, iostat_end
   implicit none
   real(dp) :: bounds(4), seed_u(10), seed_tau(10), corners(9), nodal(0:2,0:2,2,2)
   type(auto_diff_real_2var_order3) :: coordinate_u, coordinate_tau, value
   integer :: io_status, status, endpoint_u, endpoint_tau
   do
      read(input_unit,*,iostat=io_status) bounds
      if (io_status == iostat_end) exit
      if (io_status /= 0) error stop 'invalid bounds'
      read(input_unit,*) seed_u
      read(input_unit,*) seed_tau
      do endpoint_tau=1,2
         do endpoint_u=1,2
            read(input_unit,*) corners
            nodal(:,:,endpoint_u,endpoint_tau)=reshape(corners,[3,3])
         end do
      end do
      coordinate_u=seed(seed_u)
      coordinate_tau=seed(seed_tau)
      call evaluate_cell_ad(nodal,bounds(1),bounds(2),bounds(3),bounds(4), &
         coordinate_u,coordinate_tau,value,status)
      write(output_unit,'(i0,1x,10(es25.16e3,1x))') status,value%val,value%d1val2,value%d1val1, &
         value%d2val2,value%d1val1_d1val2,value%d2val1,value%d3val2,value%d1val1_d2val2, &
         value%d2val1_d1val2,value%d3val1
   end do
contains
   function seed(components) result(value)
      real(dp), intent(in) :: components(10)
      type(auto_diff_real_2var_order3) :: value
      value=0.0_dp
      value%val=components(1)
      value%d1val2=components(2)
      value%d1val1=components(3)
      value%d2val2=components(4)
      value%d1val1_d1val2=components(5)
      value%d2val1=components(6)
      value%d3val2=components(7)
      value%d1val1_d2val2=components(8)
      value%d2val1_d1val2=components(9)
      value%d3val1=components(10)
   end function seed
end program jet_driver
