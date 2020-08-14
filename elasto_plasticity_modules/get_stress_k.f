c
c
c
      type(Tensor2) function get_stress_k( stress_t, HillT_H, alpha,
     &                                     gamma_k, hardening_type, cm )
c
      use Tensor
      use cm_manager
      !implicit none
c
      type(Tensor2) :: Eye, stress_t
      type(Tensor4) :: HillT_H
      dimension cm(*)
      double precision gamma_k, alpha
      real yield_stress, shearMod_mu
      integer :: hardening_type
c Material parameters
      yield_stress = cm_get('yieldStress_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c @note
c Here we invert a fourth order tensor (lots of fun)
c
         ! get_stress_k = invert ( identity4(Eye)
         !                        + 2. * shearMod_mu * gamma_k
         !                          / ( sqrt(2./3.)
         !                              * ( yield_stress - 
         !                                  get_hardeningStress_R( alpha, 
         !&                                             hardening_type, cm ))
         !                             ) * HillT_H )
         !                * stress_t
      get_stress_k = stress_t
      write(*,*) "Just testing, missing invert(fourth orer)"
c      
      end function get_stress_k
