module freeeos_biquintic
   use const_def, only: dp
   use auto_diff_real_2var_order3_module
   use helm_polynomials, only: psi0, psi1, psi2, dpsi0, dpsi1, dpsi2, &
      ddpsi0, ddpsi1, ddpsi2, dddpsi0, dddpsi1, dddpsi2, h5
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: freeeos_partials, evaluate_cell_real, evaluate_cell_ad
   public :: cell_ok, cell_invalid_input, cell_invalid_spacing, cell_out_of_bounds, cell_arithmetic_failure

   integer, parameter :: cell_ok = 0, cell_invalid_input = 1, cell_invalid_spacing = 2
   integer, parameter :: cell_out_of_bounds = 3, cell_arithmetic_failure = 4
   integer, parameter :: packed_u_orders(9) = [0, 0, 0, 1, 2, 1, 2, 1, 2]
   integer, parameter :: packed_tau_orders(9) = [0, 1, 2, 0, 0, 1, 1, 2, 2]

   type :: freeeos_partials
      real(dp) :: f, f_u, f_tau, f_uu, f_u_tau, f_tau_tau
      real(dp) :: f_uuu, f_uu_tau, f_u_tau_tau, f_tau_tau_tau
   end type freeeos_partials

contains

   subroutine evaluate_cell_real(nodal, u0, u1, t0, t1, u, tau, partials, status)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), u0, u1, t0, t1, u, tau
      type(freeeos_partials), intent(out) :: partials
      integer, intent(out) :: status
      real(dp) :: spacing_u, spacing_tau, weights_u(6,0:3), weights_tau(6,0:3), fi5(36), nan

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      partials = freeeos_partials(nan, nan, nan, nan, nan, nan, nan, nan, nan, nan)
      call validate_cell(nodal, u0, u1, t0, t1, u, tau, spacing_u, spacing_tau, status)
      if (status /= cell_ok) return
      call pack_fi5(nodal, fi5)
      call qwts((u-u0)/spacing_u, spacing_u, weights_u)
      call qwts((tau-t0)/spacing_tau, spacing_tau, weights_tau)
      partials%f = tensor_real(fi5, weights_u(:,0), weights_tau(:,0))
      partials%f_u = tensor_real(fi5, weights_u(:,1), weights_tau(:,0))
      partials%f_tau = tensor_real(fi5, weights_u(:,0), weights_tau(:,1))
      partials%f_uu = tensor_real(fi5, weights_u(:,2), weights_tau(:,0))
      partials%f_u_tau = tensor_real(fi5, weights_u(:,1), weights_tau(:,1))
      partials%f_tau_tau = tensor_real(fi5, weights_u(:,0), weights_tau(:,2))
      partials%f_uuu = tensor_real(fi5, weights_u(:,3), weights_tau(:,0))
      partials%f_uu_tau = tensor_real(fi5, weights_u(:,2), weights_tau(:,1))
      partials%f_u_tau_tau = tensor_real(fi5, weights_u(:,1), weights_tau(:,2))
      partials%f_tau_tau_tau = tensor_real(fi5, weights_u(:,0), weights_tau(:,3))
      if (.not. all(ieee_is_finite(partial_components(partials)))) then
         status = cell_arithmetic_failure
         partials = freeeos_partials(nan, nan, nan, nan, nan, nan, nan, nan, nan, nan)
      end if
   end subroutine evaluate_cell_real

   subroutine evaluate_cell_ad(nodal, u0, u1, t0, t1, u, tau, f, status)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), u0, u1, t0, t1
      type(auto_diff_real_2var_order3), intent(in) :: u, tau
      type(auto_diff_real_2var_order3), intent(out) :: f
      integer, intent(out) :: status
      real(dp) :: spacing_u, spacing_tau, fi5(36), nan
      type(auto_diff_real_2var_order3) :: weights_u(6), weights_tau(6)

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      f = auto_diff_real_2var_order3(nan, nan, nan, nan, nan, nan, nan, nan, nan, nan)
      call validate_cell(nodal, u0, u1, t0, t1, u%val, tau%val, spacing_u, spacing_tau, status)
      if (status /= cell_ok) return
      if (.not. all(ieee_is_finite(ad_components(u))) .or. &
          .not. all(ieee_is_finite(ad_components(tau)))) then
         status = cell_invalid_input
         return
      end if
      call pack_fi5(nodal, fi5)
      call qwts_ad((u-u0)/spacing_u, spacing_u, weights_u)
      call qwts_ad((tau-t0)/spacing_tau, spacing_tau, weights_tau)
      f = tensor_ad(fi5, weights_u, weights_tau)
      if (.not. all(ieee_is_finite(ad_components(f)))) then
         status = cell_arithmetic_failure
         f = auto_diff_real_2var_order3(nan, nan, nan, nan, nan, nan, nan, nan, nan, nan)
      end if
   end subroutine evaluate_cell_ad

   subroutine validate_cell(nodal, u0, u1, t0, t1, u, tau, spacing_u, spacing_tau, status)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), u0, u1, t0, t1, u, tau
      real(dp), intent(out) :: spacing_u, spacing_tau
      integer, intent(out) :: status

      status = cell_invalid_input
      spacing_u = 0.0_dp
      spacing_tau = 0.0_dp
      if (.not. all(ieee_is_finite([u0, u1, t0, t1, u, tau]))) return
      if (.not. all(ieee_is_finite(nodal))) return
      status = cell_invalid_spacing
      if (.not. valid_spacing(u0, u1, spacing_u)) return
      if (.not. valid_spacing(t0, t1, spacing_tau)) return
      status = cell_out_of_bounds
      if (u < u0 .or. u > u1 .or. tau < t0 .or. tau > t1) return
      status = cell_ok
   end subroutine validate_cell

   logical function valid_spacing(lower, upper, spacing)
      real(dp), intent(in) :: lower, upper
      real(dp), intent(out) :: spacing
      real(dp), parameter :: minimum_spacing = 8.0_dp * huge(1.0_dp)**(-1.0_dp/3.0_dp)
      real(dp), parameter :: maximum_spacing = 0.5_dp * sqrt(huge(1.0_dp))

      valid_spacing = .false.
      spacing = 0.0_dp
      if (upper <= lower) return
      if (lower < 0.0_dp) then
         if (upper > huge(1.0_dp) + lower) return
      end if
      spacing = upper - lower
      if (spacing < minimum_spacing .or. spacing > maximum_spacing) return
      valid_spacing = .true.
   end function valid_spacing

   subroutine pack_fi5(nodal, fi5)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2)
      real(dp), intent(out) :: fi5(36)
      integer :: block, corner_u, corner_tau, position

      do block = 1, 9
         do corner_tau = 1, 2
            do corner_u = 1, 2
               position = 4*(block-1) + 2*(corner_tau-1) + corner_u
               fi5(position) = nodal(packed_u_orders(block), packed_tau_orders(block), corner_u, corner_tau)
            end do
         end do
      end do
   end subroutine pack_fi5

   subroutine qwts(fraction, spacing, weights)
      real(dp), intent(in) :: fraction, spacing
      real(dp), intent(out) :: weights(6,0:3)
      real(dp) :: complement, inverse, inverse_squared, inverse_cubed, spacing_squared

      complement = 1.0_dp - fraction
      inverse = 1.0_dp/spacing
      inverse_squared = inverse*inverse
      inverse_cubed = inverse_squared*inverse
      spacing_squared = spacing*spacing
      weights(:,0) = [psi0(fraction), psi1(fraction)*spacing, psi2(fraction)*spacing_squared, &
         psi0(complement), -psi1(complement)*spacing, psi2(complement)*spacing_squared]
      weights(:,1) = [dpsi0(fraction)*inverse, dpsi1(fraction), dpsi2(fraction)*spacing, &
         -dpsi0(complement)*inverse, dpsi1(complement), -dpsi2(complement)*spacing]
      weights(:,2) = [ddpsi0(fraction)*inverse_squared, ddpsi1(fraction)*inverse, ddpsi2(fraction), &
         ddpsi0(complement)*inverse_squared, -ddpsi1(complement)*inverse, ddpsi2(complement)]
      weights(:,3) = [dddpsi0(fraction)*inverse_cubed, dddpsi1(fraction)*inverse_squared, &
         dddpsi2(fraction)*inverse, -dddpsi0(complement)*inverse_cubed, &
         dddpsi1(complement)*inverse_squared, -dddpsi2(complement)*inverse]
   end subroutine qwts

   real(dp) function tensor_real(fi5, weights_u, weights_tau) result(value)
      real(dp), intent(in) :: fi5(36), weights_u(6), weights_tau(6)

      value = h5(1, 1, fi5, &
         weights_tau(1), weights_tau(2), weights_tau(3), weights_tau(4), weights_tau(5), weights_tau(6), &
         weights_u(1), weights_u(2), weights_u(3), weights_u(4), weights_u(5), weights_u(6))
   end function tensor_real

   subroutine qwts_ad(fraction, spacing, weights)
      type(auto_diff_real_2var_order3), intent(in) :: fraction
      real(dp), intent(in) :: spacing
      type(auto_diff_real_2var_order3), intent(out) :: weights(6)
      type(auto_diff_real_2var_order3) :: complement

      complement = 1.0_dp - fraction
      weights(1) = psi0_ad(fraction)
      weights(2) = spacing*psi1_ad(fraction)
      weights(3) = (spacing*spacing)*psi2_ad(fraction)
      weights(4) = psi0_ad(complement)
      weights(5) = (-spacing)*psi1_ad(complement)
      weights(6) = (spacing*spacing)*psi2_ad(complement)
   end subroutine qwts_ad

   function psi0_ad(argument) result(value)
      type(auto_diff_real_2var_order3), intent(in) :: argument
      type(auto_diff_real_2var_order3) :: value

      value = argument*argument*argument*(argument*(-6.0_dp*argument + 15.0_dp) - 10.0_dp) + 1.0_dp
   end function psi0_ad

   function psi1_ad(argument) result(value)
      type(auto_diff_real_2var_order3), intent(in) :: argument
      type(auto_diff_real_2var_order3) :: value

      value = argument*(argument*argument*(argument*(-3.0_dp*argument + 8.0_dp) - 6.0_dp) + 1.0_dp)
   end function psi1_ad

   function psi2_ad(argument) result(value)
      type(auto_diff_real_2var_order3), intent(in) :: argument
      type(auto_diff_real_2var_order3) :: value

      value = 0.5_dp*argument*argument*(argument*(argument*(-argument + 3.0_dp) - 3.0_dp) + 1.0_dp)
   end function psi2_ad

   function tensor_ad(fi5, weights_u, weights_tau) result(value)
      real(dp), intent(in) :: fi5(36)
      type(auto_diff_real_2var_order3), intent(in) :: weights_u(6), weights_tau(6)
      type(auto_diff_real_2var_order3) :: value
      integer :: block, corner_u, corner_tau, position, index_u, index_tau

      value = 0.0_dp
      do block = 1, 9
         do corner_tau = 1, 2
            index_tau = packed_tau_orders(block) + 1 + 3*(corner_tau-1)
            do corner_u = 1, 2
               index_u = packed_u_orders(block) + 1 + 3*(corner_u-1)
               position = 4*(block-1) + 2*(corner_tau-1) + corner_u
               value = value + fi5(position)*(weights_u(index_u)*weights_tau(index_tau))
            end do
         end do
      end do
   end function tensor_ad

   pure function partial_components(partials) result(components)
      type(freeeos_partials), intent(in) :: partials
      real(dp) :: components(10)

      components = [partials%f, partials%f_u, partials%f_tau, partials%f_uu, partials%f_u_tau, &
         partials%f_tau_tau, partials%f_uuu, partials%f_uu_tau, partials%f_u_tau_tau, partials%f_tau_tau_tau]
   end function partial_components

   pure function ad_components(value) result(components)
      type(auto_diff_real_2var_order3), intent(in) :: value
      real(dp) :: components(10)

      components = [value%val, value%d1val2, value%d1val1, value%d2val2, value%d1val1_d1val2, &
         value%d2val1, value%d3val2, value%d1val1_d2val2, value%d2val1_d1val2, value%d3val1]
   end function ad_components

end module freeeos_biquintic
