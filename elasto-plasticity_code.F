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
     & stress, tangent_C, failedIt )
        use Tensor
        use TensorXkinematics
        use hsv_manager
        use cm_manager
        use elasto_plasticity_modules
        !implicit none
c
c DECLARATIONS:
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
      double precision gamma_k, alpha_k, dPhi_dgamma, tol
      real hardStress_R_inf, hardMod_K_exp
      integer k, hardening_type, n_max_iterations
      logical failedIt
c Material parameters
      real, dimension(6) :: cm
      !lame_lambda = cm_get('lame_lambda_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
      !bulkMod_kappa = lame_lambda + 2./3. * shearMod_mu
      bulkMod_kappa = cm_get('bulkMod_kappa___',cm)
      hardening_type = cm_get('hardening_type__',cm)
c
c USER PARAMETERS:
c QP-level convergence tolerance
      tol = 1e-8
c Maximum nbr of iterations of qp level
      n_max_iterations = 15
c
c ALGORITHM: \n   
c Second order identity tensor
      Eye = identity2(Eye)
c Fourth order tensors
      IxI = Eye.dya.Eye
      I_dev = deviatoric_I4(Eye)
c History variables
      alpha_n = hsv_get_scalar('alpha', hsv)
      eps_p = hsv_get_symTen2('eps_p', hsv)
c Set up the Hill tensor for anisotropy or isotropy
      HillT_H = setup_Hill_tensor(cm)
c Trial stress tensor
      stress_t = get_stress_T_t( sstrain, hardening_type, cm, hsv)
c Trial yield function
      Phi_t = get_plastic_yield_fnc( stress_t, HillT_H, alpha_n,
     &                               hardening_type, cm )     
c Check the trial yield function
      if ( Phi_t < -tol ) then ! elastic
          ! The elastic trial stress assumption is correct,
          ! so we can accept the trial stress as the output stress
           stress = stress_t
          ! Keep the history unchanged for this elastic step
          ! Compute the elastic tangent
           tangent_C = bulkMod_kappa * IxI
     &                 + 2. * shearMod_mu * I_dev
      else ! plastic
          n = update_direction_n_n1 ( stress_t, HillT_H )
          gamma_k = 0.
          alpha_k = alpha_n
c           
          do k=1,n_max_iterations
              ! Compute subtangent and update the Lagrange multiplier
               dPhi_dgamma = get_dPhi_dgamma( alpha_k, HillT_H, n,  
     &                                        stress_t, gamma_k, Phi_t,
     &                                        hardening_type, cm, hsv )
               gamma_k = gamma_k - Phi_t / dPhi_dgamma
              ! Update the internal variable
               alpha_k = get_intVar_alpha( alpha_n, gamma_k,
     &                                    hardening_type, cm )
              ! Update the stress
               stress = get_stress_k( stress_t, HillT_H, alpha_k,
     &                                gamma_k, hardening_type, cm )
              ! Update the evolution direction (needed for anisotropy)
               n = update_direction_n_n1 ( stress, HillT_H )
              ! Compute the new yield function
              Phi_t = get_plastic_yield_fnc( stress, HillT_H, alpha_k,
     &                                       hardening_type, cm)
              ! Check the new yield function
               if ( abs(Phi_t) < tol ) then
                  exit ! converged
               elseif ( k > (n_max_iterations-1) ) then
         !             write(*,*) "elasto_plasticity_code<< Failed to 
         !&converge in iterations on material point level with Phi=",Phi_t
                  ! @todo Works basically, also the "further will be surpressed", but standard surrounding text is nonsense
                  call usermsg('elasto_plasticity_code<< Failed to 
     &converge in iterations on material point level ')
                  failedIt = .true.
               endif
          enddo
          ! Update the history         
           call hsv_set_scalar( alpha_k, 'alpha', hsv)
           call hsv_set_symTen2( (eps_p + gamma_k * n), 'eps_p', hsv)
          ! Tangent for plastic step
           tangent_C = get_tangent_C( stress, HillT_H, alpha_k, n, 
     &                                gamma_k, hardening_type, cm, hsv )
      endif
c
      return
      end subroutine