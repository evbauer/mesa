program test_biquintic
   use const_def, only: dp, ln10
   use auto_diff_real_2var_order3_module
   use freeeos_biquintic
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   integer, parameter :: orders_u(10) = [0, 1, 0, 2, 1, 0, 3, 2, 1, 0]
   integer, parameter :: orders_tau(10) = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3]
   real(dp), parameter :: energy_scale = 2.5e12_dp
   real(dp), parameter :: probes(7) = [0.0_dp, 1.0e-10_dp, 0.17_dp, 0.5_dp, &
      0.89_dp, 1.0_dp-1.0e-10_dp, 1.0_dp]
   real(dp) :: coefficients(0:5,0:5)
   integer :: failures = 0, evaluations = 0, rejections = 0

   call test_tensor_polynomials()
   call test_endpoint_data()
   call test_real_ad_agreement()
   call test_exponential_convergence()
   call test_shared_faces()
   call test_sheared_coordinates()
   call test_nonlinear_ad_inputs()
   call test_invalid_inputs()
   write(*, '(a,i0,a,i0)') 'Successful real/AD cell pairs: ', evaluations, '; rejected calls: ', rejections
   if (failures /= 0) then
      write(*, '(a,i0)') 'FAIL assertions: ', failures
      stop 1
   end if
   write(*, '(a)') 'PASS all standalone FreeEOS biquintic tests'

contains

   subroutine report_error(label, measured, tolerance)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: measured, tolerance

      if (ieee_is_finite(measured) .and. measured <= tolerance) then
         write(*, '(a,a,a,es12.4,a,es12.4)') 'PASS ', label, ': ', measured, ' <= ', tolerance
      else
         failures = failures + 1
         write(*, '(a,a,a,es12.4,a,es12.4)') 'FAIL ', label, ': ', measured, ' > ', tolerance
      end if
   end subroutine report_error

   subroutine require_ok(status)
      integer, intent(in) :: status

      if (status == cell_ok) return
      write(*, '(a,i0)') 'FAIL unexpected evaluation status: ', status
      stop 1
   end subroutine require_ok

   function real_vector(partials) result(values)
      type(freeeos_partials), intent(in) :: partials
      real(dp) :: values(10)

      values = [partials%f, partials%f_u, partials%f_tau, partials%f_uu, partials%f_u_tau, &
         partials%f_tau_tau, partials%f_uuu, partials%f_uu_tau, partials%f_u_tau_tau, partials%f_tau_tau_tau]
   end function real_vector

   function ad_vector(value) result(values)
      type(auto_diff_real_2var_order3), intent(in) :: value
      real(dp) :: values(10)

      values = [value%val, value%d1val2, value%d1val1, value%d2val2, value%d1val1_d1val2, &
         value%d2val1, value%d3val2, value%d1val1_d2val2, value%d2val1_d1val2, value%d3val1]
   end function ad_vector

   real(dp) function cell_coordinate(lower, upper, fraction) result(coordinate)
      real(dp), intent(in) :: lower, upper, fraction

      if (fraction == 0.0_dp) then
         coordinate = lower
      else if (fraction == 1.0_dp) then
         coordinate = upper
      else
         coordinate = lower + fraction*(upper-lower)
      end if
   end function cell_coordinate

   subroutine evaluate_pair(nodal, bounds, u, tau, real_values, ad_values)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), bounds(4), u, tau
      real(dp), intent(out) :: real_values(10), ad_values(10)
      type(freeeos_partials) :: partials
      type(auto_diff_real_2var_order3) :: u_ad, tau_ad, f_ad
      integer :: status

      call evaluate_cell_real(nodal, bounds(1), bounds(2), bounds(3), bounds(4), u, tau, partials, status)
      call require_ok(status)
      u_ad = u
      u_ad%d1val2 = 1.0_dp
      tau_ad = tau
      tau_ad%d1val1 = 1.0_dp
      call evaluate_cell_ad(nodal, bounds(1), bounds(2), bounds(3), bounds(4), u_ad, tau_ad, f_ad, status)
      call require_ok(status)
      real_values = real_vector(partials)
      ad_values = ad_vector(f_ad)
      evaluations = evaluations + 1
   end subroutine evaluate_pair

   real(dp) function monomial_derivative(coordinate, power, order) result(value)
      real(dp), intent(in) :: coordinate
      integer, intent(in) :: power, order
      integer :: factor

      value = 0.0_dp
      if (order > power) return
      value = coordinate**(power-order)
      do factor = power-order+1, power
         value = value*real(factor, dp)
      end do
   end function monomial_derivative

   real(dp) function polynomial_derivative(u, tau, order_u, order_tau) result(value)
      real(dp), intent(in) :: u, tau
      integer, intent(in) :: order_u, order_tau
      integer :: power_u, power_tau

      value = 0.0_dp
      do power_tau = 0, 5
         do power_u = 0, 5
            value = value + coefficients(power_u,power_tau)* &
               monomial_derivative(u, power_u, order_u)*monomial_derivative(tau, power_tau, order_tau)
         end do
      end do
   end function polynomial_derivative

   real(dp) function exponential_derivative(u, tau, order_u, order_tau) result(value)
      real(dp), intent(in) :: u, tau
      integer, intent(in) :: order_u, order_tau

      value = energy_scale*(0.8_dp**order_u*(-1.1_dp)**order_tau*exp(0.8_dp*u-1.1_dp*tau) + &
         0.31_dp*(-0.45_dp)**order_u*0.65_dp**order_tau*exp(-0.45_dp*u+0.65_dp*tau))
   end function exponential_derivative

   real(dp) function physical_polynomial(x, tau, order_x, order_tau) result(value)
      real(dp), intent(in) :: x, tau
      integer, intent(in) :: order_x, order_tau
      integer :: power_x, power_tau

      value = 0.0_dp
      do power_tau = 0, 2
         do power_x = 0, 3
            value = value + energy_scale*(-1.0_dp)**(power_x+power_tau)/real(1+power_x+2*power_tau, dp)* &
               monomial_derivative(x, power_x, order_x)*monomial_derivative(tau-14.0_dp, power_tau, order_tau)
         end do
      end do
   end function physical_polynomial

   integer function binomial(power, selected) result(value)
      integer, intent(in) :: power, selected
      integer :: factor

      value = 1
      do factor = 1, selected
         value = value*(power-factor+1)/factor
      end do
   end function binomial

   real(dp) function sheared_derivative(u, tau, order_u, order_tau) result(value)
      real(dp), intent(in) :: u, tau
      integer, intent(in) :: order_u, order_tau
      integer :: selected
      real(dp) :: x

      x = u + 2.0_dp*tau - 12.0_dp*ln10
      value = 0.0_dp
      do selected = 0, order_tau
         value = value + real(binomial(order_tau, selected), dp)*2.0_dp**selected* &
            physical_polynomial(x, tau, order_u+selected, order_tau-selected)
      end do
   end function sheared_derivative

   real(dp) function analytic_derivative(family, u, tau, order_u, order_tau) result(value)
      integer, intent(in) :: family, order_u, order_tau
      real(dp), intent(in) :: u, tau

      select case (family)
      case (1)
         value = polynomial_derivative(u, tau, order_u, order_tau)
      case (2)
         value = exponential_derivative(u, tau, order_u, order_tau)
      case (3)
         value = sheared_derivative(u, tau, order_u, order_tau)
      case default
         stop 2
      end select
   end function analytic_derivative

   subroutine fill_nodes(family, bounds, nodal)
      integer, intent(in) :: family
      real(dp), intent(in) :: bounds(4)
      real(dp), intent(out) :: nodal(0:2,0:2,2,2)
      integer :: corner_u, corner_tau, order_u, order_tau

      do corner_tau = 1, 2
         do corner_u = 1, 2
            do order_tau = 0, 2
               do order_u = 0, 2
                  nodal(order_u,order_tau,corner_u,corner_tau) = analytic_derivative(family, &
                     bounds(corner_u), bounds(2+corner_tau), order_u, order_tau)
               end do
            end do
         end do
      end do
   end subroutine fill_nodes

   subroutine test_tensor_polynomials()
      real(dp), parameter :: bounds(4) = [-0.7_dp, 0.9_dp, 0.3_dp, 1.4_dp]
      real(dp) :: nodal(0:2,0:2,2,2), real_values(10), ad_values(10), expected, u, tau
      real(dp) :: max_real, max_ad, dense_real, dense_ad, scale
      integer :: power_u, power_tau, probe_u, probe_tau, component, polynomial

      max_real = 0.0_dp
      max_ad = 0.0_dp
      dense_real = 0.0_dp
      dense_ad = 0.0_dp
      do polynomial = 1, 37
         coefficients = 0.0_dp
         if (polynomial <= 36) then
            power_u = mod(polynomial-1, 6)
            power_tau = (polynomial-1)/6
            coefficients(power_u,power_tau) = 1.0_dp
         else
            do power_tau = 0, 5
               do power_u = 0, 5
                  coefficients(power_u,power_tau) = energy_scale*(-1.0_dp)**(power_u+power_tau)/ &
                     real((power_u+1)*(power_tau+1), dp)
               end do
            end do
         end if
         call fill_nodes(1, bounds, nodal)
         do probe_tau = 1, size(probes)
            tau = cell_coordinate(bounds(3), bounds(4), probes(probe_tau))
            do probe_u = 1, size(probes)
               u = cell_coordinate(bounds(1), bounds(2), probes(probe_u))
               call evaluate_pair(nodal, bounds, u, tau, real_values, ad_values)
               do component = 1, 10
                  expected = polynomial_derivative(u, tau, orders_u(component), orders_tau(component))
                  if (polynomial <= 36) then
                     scale = max(1.0_dp, abs(expected))
                     max_real = max(max_real, abs(real_values(component)-expected)/scale)
                     max_ad = max(max_ad, abs(ad_values(component)-expected)/scale)
                  else
                     scale = max(energy_scale, abs(expected))
                     dense_real = max(dense_real, abs(real_values(component)-expected)/scale)
                     dense_ad = max(dense_ad, abs(ad_values(component)-expected)/scale)
                  end if
               end do
            end do
         end do
      end do
      call report_error('36 tensor monomials, all ten real partials', max_real, 2.0e-11_dp)
      call report_error('36 tensor monomials, all ten AD partials', max_ad, 2.0e-11_dp)
      call report_error('dense degree-(5,5) polynomial real, erg/g scale', dense_real, 2.0e-11_dp)
      call report_error('dense degree-(5,5) polynomial AD, erg/g scale', dense_ad, 2.0e-11_dp)
   end subroutine test_tensor_polynomials

   subroutine arbitrary_nodes(nodal)
      real(dp), intent(out) :: nodal(0:2,0:2,2,2)
      integer :: order_u, order_tau, corner_u, corner_tau

      do corner_tau = 1, 2
         do corner_u = 1, 2
            do order_tau = 0, 2
               do order_u = 0, 2
                  nodal(order_u,order_tau,corner_u,corner_tau) = energy_scale* &
                     sin(real(1+order_u+3*order_tau+9*(corner_u-1)+18*(corner_tau-1), dp))
               end do
            end do
         end do
      end do
   end subroutine arbitrary_nodes

   subroutine test_endpoint_data()
      real(dp), parameter :: bounds(4) = [-2.0_dp, -0.4_dp, 1.0_dp, 1.7_dp]
      integer, parameter :: endpoint_components(8) = [1, 2, 3, 4, 5, 6, 8, 9]
      real(dp) :: nodal(0:2,0:2,2,2), real_values(10), ad_values(10), expected, max_error
      integer :: corner_u, corner_tau, component, endpoint_component

      call arbitrary_nodes(nodal)
      max_error = 0.0_dp
      do corner_tau = 1, 2
         do corner_u = 1, 2
            call evaluate_pair(nodal, bounds, bounds(corner_u), bounds(2+corner_tau), real_values, ad_values)
            do endpoint_component = 1, size(endpoint_components)
               component = endpoint_components(endpoint_component)
               expected = nodal(orders_u(component),orders_tau(component),corner_u,corner_tau)
               max_error = max(max_error, abs(real_values(component)-expected)/energy_scale, &
                  abs(ad_values(component)-expected)/energy_scale)
            end do
         end do
      end do
      call report_error('arbitrary endpoint jets through (2,1) and (1,2), real/AD', max_error, 5.0e-13_dp)
   end subroutine test_endpoint_data

   subroutine test_real_ad_agreement()
      real(dp), parameter :: bounds(4) = [17.1_dp, 17.6_dp, 12.2_dp, 13.05_dp]
      real(dp) :: nodal(0:2,0:2,2,2), real_values(10), ad_values(10), u, tau, scale, max_error
      integer :: probe_u, probe_tau, component

      call arbitrary_nodes(nodal)
      max_error = 0.0_dp
      do probe_tau = 1, size(probes)
         tau = cell_coordinate(bounds(3), bounds(4), probes(probe_tau))
         do probe_u = 1, size(probes)
            u = cell_coordinate(bounds(1), bounds(2), probes(probe_u))
            call evaluate_pair(nodal, bounds, u, tau, real_values, ad_values)
            do component = 1, 10
               scale = energy_scale/(bounds(2)-bounds(1))**orders_u(component)/ &
                  (bounds(4)-bounds(3))**orders_tau(component)
               max_error = max(max_error, abs(real_values(component)-ad_values(component))/scale)
            end do
         end do
      end do
      call report_error('arbitrary corner real/AD agreement incl. near/exact faces', max_error, 2.0e-12_dp)
   end subroutine test_real_ad_agreement

   subroutine test_exponential_convergence()
      real(dp), parameter :: fractions(5) = [0.11_dp, 0.29_dp, 0.47_dp, 0.71_dp, 0.89_dp]
      real(dp) :: bounds(4), nodal(0:2,0:2,2,2), real_values(10), ad_values(10)
      real(dp) :: errors(0:3,4), rates(0:3), h, u, tau, expected, ad_error
      integer :: level, probe_u, probe_tau, component, total_order

      errors = 0.0_dp
      ad_error = 0.0_dp
      do level = 1, 4
         h = 0.5_dp**(level-1)
         bounds = [0.3_dp, 0.3_dp+h, -0.2_dp, -0.2_dp+0.7_dp*h]
         call fill_nodes(2, bounds, nodal)
         do probe_tau = 1, size(fractions)
            tau = bounds(3) + fractions(probe_tau)*(bounds(4)-bounds(3))
            do probe_u = 1, size(fractions)
               u = bounds(1) + fractions(probe_u)*(bounds(2)-bounds(1))
               call evaluate_pair(nodal, bounds, u, tau, real_values, ad_values)
               do component = 1, 10
                  expected = exponential_derivative(u, tau, orders_u(component), orders_tau(component))
                  total_order = orders_u(component) + orders_tau(component)
                  errors(total_order,level) = max(errors(total_order,level), &
                     abs(real_values(component)-expected)/energy_scale)
               end do
               ad_error = max(ad_error, maxval(abs(real_values-ad_values))/energy_scale)
            end do
         end do
         write(*, '(a,f7.4,a,4(es12.4,1x))') 'exp h=', h, ' max errors / (erg/g scale), orders 0..3: ', errors(:,level)
      end do
      do total_order = 0, 3
         rates(total_order) = minval(log(errors(total_order,1:3)/errors(total_order,2:4))/log(2.0_dp))
      end do
      write(*, '(a,4(f8.4,1x))') 'exp minimum observed orders 0..3: ', rates
      call report_error('exp value convergence order >= 5.5', 5.5_dp-rates(0), 0.0_dp)
      call report_error('exp first partial convergence order >= 4.5', 4.5_dp-rates(1), 0.0_dp)
      call report_error('exp second partial convergence order >= 3.5', 3.5_dp-rates(2), 0.0_dp)
      call report_error('exp third partial convergence order >= 2.5', 2.5_dp-rates(3), 0.0_dp)
      call report_error('exp finest-grid value error', errors(0,4), 1.0e-10_dp)
      call report_error('exp real/AD agreement all partials', ad_error, 2.0e-10_dp)
   end subroutine test_exponential_convergence

   subroutine test_shared_faces()
      real(dp) :: left_bounds(4), right_bounds(4), left_nodes(0:2,0:2,2,2), right_nodes(0:2,0:2,2,2)
      real(dp) :: left_values(10), right_values(10), ad_values(10), u, tau, c2_error, transverse_error
      real(dp) :: jumps(2), face_fraction
      integer :: axis, probe, component, normal_third

      c2_error = 0.0_dp
      transverse_error = 0.0_dp
      jumps = 0.0_dp
      do axis = 1, 2
         if (axis == 1) then
            left_bounds = [-0.8_dp, 0.2_dp, -0.4_dp, 0.6_dp]
            right_bounds = [0.2_dp, 0.9_dp, -0.4_dp, 0.6_dp]
            normal_third = 7
         else
            left_bounds = [-0.4_dp, 0.6_dp, -0.8_dp, 0.2_dp]
            right_bounds = [-0.4_dp, 0.6_dp, 0.2_dp, 0.9_dp]
            normal_third = 10
         end if
         call fill_nodes(2, left_bounds, left_nodes)
         call fill_nodes(2, right_bounds, right_nodes)
         do probe = 0, 10
            face_fraction = real(probe, dp)/10.0_dp
            if (axis == 1) then
               u = left_bounds(2)
               tau = cell_coordinate(left_bounds(3), left_bounds(4), face_fraction)
            else
               u = cell_coordinate(left_bounds(1), left_bounds(2), face_fraction)
               tau = left_bounds(4)
            end if
            call evaluate_pair(left_nodes, left_bounds, u, tau, left_values, ad_values)
            call evaluate_pair(right_nodes, right_bounds, u, tau, right_values, ad_values)
            c2_error = max(c2_error, maxval(abs(left_values(1:6)-right_values(1:6)))/energy_scale)
            do component = 7, 10
               if (component == normal_third) cycle
               transverse_error = max(transverse_error, abs(left_values(component)-right_values(component))/energy_scale)
            end do
            jumps(axis) = max(jumps(axis), abs(left_values(normal_third)-right_values(normal_third))/energy_scale)
         end do
      end do
      call report_error('two-cell C2 faces in both orientations', c2_error, 2.0e-12_dp)
      call report_error('two-cell tangential/mixed third partials', transverse_error, 2.0e-11_dp)
      write(*, '(a,2(es12.4,1x))') 'Allowed normal C3 jumps / (erg/g scale), u and tau faces: ', jumps
      call report_error('nonzero normal C3 jumps are allowed and exercised', 1.0e-7_dp-minval(jumps), 0.0_dp)
   end subroutine test_shared_faces

   function physical_vector(native) result(physical)
      real(dp), intent(in) :: native(10)
      real(dp) :: physical(10)

      physical(1) = native(1)
      physical(2) = native(2)
      physical(3) = native(3) - 2.0_dp*native(2)
      physical(4) = native(4)
      physical(5) = native(5) - 2.0_dp*native(4)
      physical(6) = native(6) - 4.0_dp*native(5) + 4.0_dp*native(4)
      physical(7) = native(7)
      physical(8) = native(8) - 2.0_dp*native(7)
      physical(9) = native(9) - 4.0_dp*native(8) + 4.0_dp*native(7)
      physical(10) = native(10) - 6.0_dp*native(9) + 12.0_dp*native(8) - 8.0_dp*native(7)
   end function physical_vector

   subroutine test_sheared_coordinates()
      real(dp), parameter :: bounds(4) = [-2.3_dp, -1.1_dp, 13.7_dp, 14.3_dp]
      real(dp) :: nodal(0:2,0:2,2,2), native(10), ad_native(10), physical(10), ad_physical(10)
      real(dp) :: u, tau, x, expected, scale, native_error, physical_error, ad_error
      type(auto_diff_real_2var_order3) :: x_ad, tau_ad, u_ad, f_ad
      integer :: probe_u, probe_tau, component, status

      call fill_nodes(3, bounds, nodal)
      native_error = 0.0_dp
      physical_error = 0.0_dp
      ad_error = 0.0_dp
      do probe_tau = 1, size(probes)
         tau = cell_coordinate(bounds(3), bounds(4), probes(probe_tau))
         do probe_u = 1, size(probes)
            u = cell_coordinate(bounds(1), bounds(2), probes(probe_u))
            x = u + 2.0_dp*tau - 12.0_dp*ln10
            call evaluate_pair(nodal, bounds, u, tau, native, ad_native)
            physical = physical_vector(native)
            x_ad = x
            x_ad%d1val2 = 1.0_dp
            tau_ad = tau
            tau_ad%d1val1 = 1.0_dp
            u_ad = x_ad - 2.0_dp*tau_ad + 12.0_dp*ln10
            u_ad%val = u
            call evaluate_cell_ad(nodal, bounds(1), bounds(2), bounds(3), bounds(4), u_ad, tau_ad, f_ad, status)
            call require_ok(status)
            ad_physical = ad_vector(f_ad)
            do component = 1, 10
               expected = sheared_derivative(u, tau, orders_u(component), orders_tau(component))
               scale = max(energy_scale, abs(expected))
               native_error = max(native_error, abs(native(component)-expected)/scale, &
                  abs(ad_native(component)-expected)/scale)
               expected = physical_polynomial(x, tau, orders_u(component), orders_tau(component))
               scale = max(energy_scale, abs(expected))
               physical_error = max(physical_error, abs(physical(component)-expected)/scale)
               ad_error = max(ad_error, abs(ad_physical(component)-expected)/scale)
            end do
         end do
      end do
      call report_error('x=u+2*tau-12*ln10 native tensor derivatives', native_error, 2.0e-10_dp)
      call report_error('sheared physical partials through total degree 3', physical_error, 2.0e-10_dp)
      call report_error('sheared physical AD seeding through total degree 3', ad_error, 2.0e-10_dp)
   end subroutine test_sheared_coordinates

   subroutine test_nonlinear_ad_inputs()
      real(dp), parameter :: bounds(4) = [-0.7_dp, 0.9_dp, 0.3_dp, 1.4_dp]
      real(dp) :: nodal(0:2,0:2,2,2), max_error
      type(auto_diff_real_2var_order3) :: variable_one, variable_two, u_ad, tau_ad, f_ad, expected
      integer :: status

      coefficients = 0.0_dp
      coefficients(2,3) = 1.0_dp
      coefficients(5,0) = 1.0_dp
      coefficients(0,1) = 1.0_dp
      call fill_nodes(1, bounds, nodal)
      variable_one = 0.0_dp
      variable_one%d1val1 = 1.0_dp
      variable_two = 0.0_dp
      variable_two%d1val2 = 1.0_dp
      u_ad = 0.2_dp + variable_two + 0.3_dp*variable_one*variable_one
      tau_ad = 0.8_dp + variable_one - 0.2_dp*variable_two*variable_two + &
         0.15_dp*variable_one*variable_two + 0.01_dp*variable_two*variable_two*variable_two
      call evaluate_cell_ad(nodal, bounds(1), bounds(2), bounds(3), bounds(4), u_ad, tau_ad, f_ad, status)
      call require_ok(status)
      expected = u_ad*u_ad*tau_ad*tau_ad*tau_ad + u_ad*u_ad*u_ad*u_ad*u_ad + tau_ad
      max_error = maxval(abs(ad_vector(f_ad)-ad_vector(expected)))
      call report_error('AD preserves nonlinear input chain rules through order 3', max_error, 2.0e-11_dp)
   end subroutine test_nonlinear_ad_inputs

   subroutine expect_rejection(nodal, coordinates, expected_status)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), coordinates(6)
      integer, intent(in) :: expected_status
      type(freeeos_partials) :: partials
      type(auto_diff_real_2var_order3) :: u_ad, tau_ad, f_ad
      integer :: status

      call evaluate_cell_real(nodal, coordinates(1), coordinates(2), coordinates(3), coordinates(4), &
         coordinates(5), coordinates(6), partials, status)
      if (status /= expected_status .or. .not. all(ieee_is_nan(real_vector(partials)))) then
         failures = failures + 1
         write(*, '(a,2(i0,1x))') 'FAIL real rejection status, expected: ', status, expected_status
      end if
      u_ad = coordinates(5)
      tau_ad = coordinates(6)
      call evaluate_cell_ad(nodal, coordinates(1), coordinates(2), coordinates(3), coordinates(4), &
         u_ad, tau_ad, f_ad, status)
      if (status /= expected_status .or. .not. all(ieee_is_nan(ad_vector(f_ad)))) then
         failures = failures + 1
         write(*, '(a,2(i0,1x))') 'FAIL AD rejection status, expected: ', status, expected_status
      end if
      rejections = rejections + 2
   end subroutine expect_rejection

   subroutine test_invalid_inputs()
      real(dp), parameter :: valid_coordinates(6) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.4_dp, 0.8_dp]
      real(dp) :: nodal(0:2,0:2,2,2), coordinates(6), nonfinite(3), raw_ad(10)
      type(auto_diff_real_2var_order3) :: bad_ad, u_ad, tau_ad, f_ad
      integer :: coordinate, variant, order_u, order_tau, corner_u, corner_tau, axis, component, status
      integer :: failures_before, rejections_before

      failures_before = failures
      rejections_before = rejections
      nonfinite = [ieee_value(0.0_dp, ieee_quiet_nan), ieee_value(0.0_dp, ieee_positive_inf), &
         ieee_value(0.0_dp, ieee_negative_inf)]
      nodal = 0.0_dp
      do variant = 1, 3
         do coordinate = 1, 6
            coordinates = valid_coordinates
            coordinates(coordinate) = nonfinite(variant)
            call expect_rejection(nodal, coordinates, cell_invalid_input)
         end do
         do corner_tau = 1, 2
            do corner_u = 1, 2
               do order_tau = 0, 2
                  do order_u = 0, 2
                     nodal = 0.0_dp
                     nodal(order_u,order_tau,corner_u,corner_tau) = nonfinite(variant)
                     call expect_rejection(nodal, valid_coordinates, cell_invalid_input)
                  end do
               end do
            end do
         end do
      end do
      nodal = 0.0_dp
      do axis = 1, 2
         do variant = 1, 7
            coordinates = valid_coordinates
            select case (variant)
            case (1)
               coordinates(2*axis) = coordinates(2*axis-1)
            case (2)
               coordinates(2*axis) = -1.0_dp
            case (3)
               coordinates(2*axis-1) = -huge(1.0_dp)
               coordinates(2*axis) = huge(1.0_dp)
            case (4)
               coordinates(2*axis) = 1.0e-200_dp
            case (5)
               coordinates(2*axis) = 1.0e200_dp
            case (6)
               coordinates(2*axis) = tiny(1.0_dp)
            case (7)
               coordinates(2*axis) = 2.0_dp*huge(1.0_dp)**(-1.0_dp/3.0_dp)
            end select
            call expect_rejection(nodal, coordinates, cell_invalid_spacing)
         end do
         coordinates = valid_coordinates
         coordinates(4+axis) = nearest(coordinates(2*axis-1), -1.0_dp)
         call expect_rejection(nodal, coordinates, cell_out_of_bounds)
         coordinates = valid_coordinates
         coordinates(4+axis) = nearest(coordinates(2*axis), 1.0_dp)
         call expect_rejection(nodal, coordinates, cell_out_of_bounds)
      end do
      do variant = 1, 3
         do axis = 1, 2
            do component = 1, 10
               raw_ad = 0.0_dp
               raw_ad(1) = valid_coordinates(4+axis)
               raw_ad(component) = nonfinite(variant)
               bad_ad = auto_diff_real_2var_order3(raw_ad(1), raw_ad(2), raw_ad(3), raw_ad(4), raw_ad(5), &
                  raw_ad(6), raw_ad(7), raw_ad(8), raw_ad(9), raw_ad(10))
               u_ad = valid_coordinates(5)
               tau_ad = valid_coordinates(6)
               if (axis == 1) u_ad = bad_ad
               if (axis == 2) tau_ad = bad_ad
               call evaluate_cell_ad(nodal, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, u_ad, tau_ad, f_ad, status)
               if (status /= cell_invalid_input .or. .not. all(ieee_is_nan(ad_vector(f_ad)))) then
                  failures = failures + 1
                  write(*, '(a,i0)') 'FAIL nonfinite AD component rejection: ', component
               end if
               rejections = rejections + 1
            end do
         end do
      end do
      nodal = huge(1.0_dp)
      call expect_rejection(nodal, valid_coordinates, cell_arithmetic_failure)
      write(*, '(a,i0)') 'Invalid-input rejection calls exercised: ', rejections-rejections_before
      call report_error('invalid inputs rejected and all failed outputs poisoned', &
         real(failures-failures_before, dp), 0.0_dp)
   end subroutine test_invalid_inputs

end program test_biquintic
