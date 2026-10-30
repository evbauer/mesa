module arithmetic_weights
   use const_def, only: dp
   implicit none
   private
   public :: stable_native
   integer, parameter :: density_orders(10) = [0,1,0,2,1,0,3,2,1,0]
   integer, parameter :: thermal_orders(10) = [0,0,1,0,1,2,0,1,2,3]
   real(dp), parameter :: coefficients(0:5,6) = reshape([ &
      1.0_dp,0.0_dp,0.0_dp,-10.0_dp,15.0_dp,-6.0_dp, &
      0.0_dp,1.0_dp,0.0_dp,-6.0_dp,8.0_dp,-3.0_dp, &
      0.0_dp,0.0_dp,0.5_dp,-1.5_dp,1.5_dp,-0.5_dp, &
      0.0_dp,0.0_dp,0.0_dp,10.0_dp,-15.0_dp,6.0_dp, &
      0.0_dp,0.0_dp,0.0_dp,-4.0_dp,7.0_dp,-3.0_dp, &
      0.0_dp,0.0_dp,0.0_dp,0.5_dp,-1.0_dp,0.5_dp], [6,6])
contains
   subroutine stable_weights(coordinate, lower, upper, weights)
      real(dp), intent(in) :: coordinate, lower, upper
      real(dp), intent(out) :: weights(3,2,0:3)
      real(dp) :: spacing, fraction, polynomial, multiplier
      integer :: endpoint, node_order, order, degree, factor, basis
      spacing = upper-lower
      fraction = (coordinate-lower)/spacing
      do order = 0,3
         do endpoint = 1,2
            do node_order = 0,2
               basis = 3*(endpoint-1)+node_order+1
               polynomial = 0.0_dp
               do degree = 5,order,-1
                  multiplier = 1.0_dp
                  do factor = 0,order-1
                     multiplier = multiplier*real(degree-factor,dp)
                  end do
                  polynomial = polynomial*fraction+coefficients(degree,basis)*multiplier
               end do
               weights(node_order+1,endpoint,order) = polynomial*spacing**(node_order-order)
            end do
         end do
      end do
   end subroutine stable_weights

   subroutine stable_native(nodal, u0, u1, t0, t1, coordinate_u, coordinate_t, native)
      real(dp), intent(in) :: nodal(0:2,0:2,2,2), u0, u1, t0, t1, coordinate_u, coordinate_t
      real(dp), intent(out) :: native(10)
      real(dp) :: weights_u(3,2,0:3), weights_t(3,2,0:3)
      real(dp) :: anchor, datum, term, subtotal, correction, combined
      integer :: partial, endpoint_u, endpoint_t, derivative_u, derivative_t
      call stable_weights(coordinate_u,u0,u1,weights_u)
      call stable_weights(coordinate_t,t0,t1,weights_t)
      anchor = nodal(0,0,1,1)
      do partial = 1,10
         subtotal = 0.0_dp
         correction = 0.0_dp
         do endpoint_t = 1,2
            do endpoint_u = 1,2
               do derivative_t = 0,2
                  do derivative_u = 0,2
                     datum = nodal(derivative_u,derivative_t,endpoint_u,endpoint_t)
                     if (derivative_u == 0 .and. derivative_t == 0) datum = datum-anchor
                     term = datum*weights_u(derivative_u+1,endpoint_u,density_orders(partial))* &
                        weights_t(derivative_t+1,endpoint_t,thermal_orders(partial))
                     combined = subtotal+term
                     if (abs(subtotal) >= abs(term)) then
                        correction = correction+((subtotal-combined)+term)
                     else
                        correction = correction+((term-combined)+subtotal)
                     end if
                     subtotal = combined
                  end do
               end do
            end do
         end do
         native(partial) = subtotal+correction
      end do
      native(1) = native(1)+anchor
   end subroutine stable_native
end module arithmetic_weights
