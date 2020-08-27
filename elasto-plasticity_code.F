c ###############################################################################
c Small strain elasto-plasticity with elastoplastic modules hardening
c inputs:
c * material parameters (kappa, mu, hardMod_K, yieldStress, hardening_type, ...)
c * Hencky_strain_E
c * history (alpha, eps_p)
c outputs:
c * stress
c * history (with tangent 'es') 
      subroutine elasto_plasticity_code ( sstrain, hsv, cm,
     & stress, tangent_C )
        use Tensor
        use TensorXkinematics
        use hsv_manager
        use cm_manager
        use elasto_plasticity_modules
        !implicit none
c Declarations
 !      type(Tensor2), intent(in) :: sstrain
      type(Tensor2) :: sstrain
      type(Tensor2) :: stress, stress_t, stress_vol, stress_dev,
     &                 stress_dev_t, n, Eye, eps_p
      type(Tensor4) :: tangent_C, nxn, I_dev, IxI
      type(Tensor4) :: HillT_H
      double precision, dimension(*) :: hsv
      double precision Phi_t
      real lame_lambda, shearMod_mu, bulkMod_kappa
      real hardMod_K, yieldStress
      double precision alpha_n, R_t, d_lambda, norm_stress_dev_t
      double precision gamma_k, alpha_k, dPhi_dgamma
      real hardStress_R_inf, hardMod_K_exp
      integer k, hardening_type
c Material parameters
      real, dimension(6) :: cm
      lame_lambda = cm_get('lame_lambda_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
      bulkMod_kappa = lame_lambda + 2./3. * shearMod_mu
      hardening_type = cm_get('hardening_type__',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c History variables
      alpha_n = hsv_get_scalar('alpha', hsv)
      eps_p = hsv_get_symTen2('eps_p', hsv)
c
c Algorithm      
c Set up the Hill tensor for anisotropy or isotropy
      HillT_H = setup_Hill_tensor(cm)
c Trial stress tensor
      stress_t = get_stress_T_t( sstrain, hardening_type, cm, hsv)
c Trial yield function
      Phi_t = get_plastic_yield_fnc( stress_t, HillT_H, alpha_n,
     &                               hardening_type, cm )     
c Second order identity tensor
      Eye = identity2(Eye)
c Fourth order tensors
      IxI = Eye.dya.Eye
      I_dev = deviatoric_I4(Eye)
c Check the trial yield function
      if ( Phi_t < 1e-8 ) then
          stress = stress_t
          ! Keep the history unchanged for this elastic step
          tangent_C = bulkMod_kappa * IxI
     &                + 2. * shearMod_mu * I_dev
      else
          n = update_direction_n_n1 ( stress_t, HillT_H )
          gamma_k = 0.
          alpha_k = alpha_n
c           
          do k=1,10
              ! Compute update
              dPhi_dgamma = get_dPhi_dgamma( alpha_k, HillT_H, n, 
     &                            gamma_k, hardening_type, cm, hsv )
              gamma_k = gamma_k - Phi_t / dPhi_dgamma
              alpha_k = get_intVar_alpha( alpha_n, gamma_k,
     &                                    hardening_type, cm )
              stress = get_stress_k( stress_t, HillT_H, alpha_k,
     &                               gamma_k, hardening_type, cm )
              n = update_direction_n_n1 ( stress, HillT_H )
              Phi_t = get_plastic_yield_fnc( stress, HillT_H, alpha_k,
     &                                       hardening_type, cm)
              ! Check the new yield function
              if ( Phi_t < 1e-8 ) then
                  exit ! converged
              elseif ( k > 9 ) then
                  write(*,*) "elasto_plasticity_code<< Failed to 
     &converge in iterations on material point level with Phi=",Phi_t
                  pause
              endif
          enddo
c History         
          call hsv_set_scalar( alpha_k, 'alpha', hsv)
          call hsv_set_symTen2( (eps_p + gamma_k * n), 'eps_p', hsv)
c Tangent
          tangent_C = get_tangent_C( stress, HillT_H, alpha_k,
     &                            n, gamma_k, hardening_type, cm, hsv )
      endif
c
      return
      end