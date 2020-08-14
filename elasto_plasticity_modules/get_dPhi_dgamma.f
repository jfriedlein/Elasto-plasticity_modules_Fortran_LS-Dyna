c
c
c
      double precision function get_dPhi_dgamma( alpha_k, HillT_H, n_n1,
     & gamma_k, hardening_type, cm, hsv )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      !implicit none
c
      type(Tensor4) :: HillT_H, A_inv
      type(Tensor2) :: Eye, n_n1
      double precision alpha_k, gamma_k
      double precision hardStress_R, yield_stress, shearMod_mu
      dimension cm(*), hsv(*)
      integer :: hardening_type
c Material parameters
      yield_stress = cm_get('yieldStress_____',cm)
      shearMod_mu =  cm_get('shearMod_mu_____',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c      
      hardStress_R = get_hardeningStress_R(alpha_k, hardening_type, cm)
      A_inv = inv( identity4(Eye)
     &                 + (2. * shearMod_mu * gamma_k
     &                  / ( sqrt(2./3.) * (yield_stress-hardStress_R) ))
     &                  * HillT_H )
c
c @todo Check why we get factors of 2 on the main diagonal for A_inv = inv(identity) instead of 0.5 as in DII
	! @todo The following only achieves superlinear behaviour close to the solution, before it's quadratic
	get_dPhi_dgamma = - 2. * shearMod_mu * (n_n1 .ddot. A_inv .ddot. n_n1)
     &                  + sqrt( 2./3. )
     &                    * get_d_R_d_gamma( alpha_k, gamma_k, 
     &                                       hardening_type, cm, hsv )
c      
      end function get_dPhi_dgamma
