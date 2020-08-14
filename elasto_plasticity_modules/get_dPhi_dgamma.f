c
c
c
      double precision function get_dPhi_dgamma( alpha_k, HillT_H, n_n1,
     & gamma_k, hardening_type, cm, hsv )
c
      use Tensor
      use cm_manager
      !implicit none
c
      type(Tensor4) :: A_inv, HillT_H
      type(Tensor2) :: Eye, n_n1
      double precision alpha_k, gamma_k
      double precision hardStress_R, yield_stress, shearMod_mu
      dimension cm(*), hsv(*)
      integer :: hardening_type
c Material parameters
      yield_stress = cm_get('yieldStress_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c      
      hardStress_R = get_hardeningStress_R(alpha_k, hardening_type, cm)
! @todo Here we invert a fourth order tensor
         ! A_inv = invert( identity4(Eye)
         !&                 + 2. * shearMod_mu * gamma_k
         !&                  / ( sqrt(2./3.) * (yield_stress-hardStress_R) )
         !&                  * HillT_H )
c
	! @todo The following only achieves superlinear behaviour close to the solution, before it's quadratic
	get_dPhi_dgamma = - 2. * shearMod_mu * (n_n1 .ddot. A_inv .ddot. n_n1)
     &                  + sqrt( 2./3. )
     &                    * get_d_R_d_gamma( alpha_k, gamma_k, 
     &                                       hardening_type, cm, hsv )
c      
      end function get_dPhi_dgamma
