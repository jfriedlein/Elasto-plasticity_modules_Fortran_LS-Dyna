c
c
c
      double precision function get_dPhi_dgamma( alpha_k, HillT_H, n_n1,
     & stress_t, gamma_k, Phi_k, cm_all, hsv )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      !implicit none
c
      type(Tensor4) :: HillT_H, A_inv, d_Ainv_dgamma
      type(Tensor2) :: Eye, n_n1, stress_t
      real(kind=8) alpha_k
      real gamma_k, d_R_d_gamma, Phi_k
      double precision hardStress_R, yield_stress, shearMod_mu
      dimension cm_all(2,*), hsv(*)
c Material parameters
      yield_stress = cm_get_pair('yieldStress_____',cm_all)
      shearMod_mu =  cm_get_pair('shearMod_mu_____',cm_all)
c Second order identity tensor
      Eye = identity2(Eye)
c      
      d_R_d_gamma = get_d_R_d_gamma( alpha_k, gamma_k, cm_all, hsv )
c
      hardStress_R = get_hardeningStress_R(alpha_k, cm_all)
c
      A_inv = inv( identity4(Eye)
     &                 + (2. * shearMod_mu * gamma_k
     &                  / ( sqrt(2./3.) * (yield_stress-hardStress_R) ))
     &                  * HillT_H )
c
      d_Ainv_dgamma
     &    = (A_inv .ddot. A_inv)
     &      .ddot. HillT_H
     &      * ( -2. * shearMod_mu )
     &           / ( sqrt(2./3.)*(yield_stress-hardStress_R) + Phi_k )
     &      * (1. + gamma_k/( sqrt(1.5)*Phi_k
     &                        + (yield_stress-hardStress_R) )
     &              * d_R_d_gamma)
c @note Be aware of the required parentheses here:
      get_dPhi_dgamma = ( (n_n1 .ddot. d_Ainv_dgamma) .ddot. stress_t )
     &                  + sqrt(2./3.) * d_R_d_gamma
c      
      end function get_dPhi_dgamma
