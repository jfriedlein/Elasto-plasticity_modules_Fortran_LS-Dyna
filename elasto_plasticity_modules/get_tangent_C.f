c
c
c
      type(Tensor4) function get_tangent_C( stress, HillT_H, alpha_k,
     &                          n_n1, gamma_k, hardening_type, cm, hsv )
c
      use Tensor
      use cm_manager
      use TensorXkinematics
      !implicit none
c
      type(Tensor2) :: stress, n_n1, Eye
      type(Tensor4) :: HillT_H, N_four, d_Tt_d_eps, E_e
      dimension cm(*), hsv(*)
      double precision gamma_k, alpha_k
      real yield_stress, shearMod_mu
      integer :: hardening_type
c Material parameters
      !lame_lambda = cm_get('lame_lambda_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
      bulkMod_kappa = cm_get('bulkMod_kappa___',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c @note
c Here we invert a fourth order tensor (lots of fun)
c
     	N_four = ( stress ** HillT_H ** stress)**(-1.5)
     &             * ( HillT_H * ( stress ** HillT_H ** stress )
     &                 - ((HillT_H**stress).dya.(stress**HillT_H)) )
          d_Tt_d_eps = bulkMod_kappa * ( Eye.dya.Eye )
     &				 + 2. * shearMod_mu *deviatoric_I4(Eye)
		E_e = inv( inv(d_Tt_d_eps) + gamma_k * N_four )
  
		get_tangent_C = E_e
     &                    - 1. / ( n_n1 ** E_e ** n_n1 - sqrt( 2./3. )
     &                      * get_d_R_d_gamma( alpha_k, gamma_k, 
     &                                       hardening_type, cm, hsv ) )
     &                      * ( (E_e**n_n1).dya.(n_n1**E_e) )
c      
      end function get_tangent_C
