c
c
c
      type(Tensor2) function get_stress_T_t( sstrain, hardening_type, 
     &                                       cm, hsv )
c
      use Tensor
      use hsv_manager
      use cm_manager
      !implicit none
c
      type(Tensor2) :: sstrain, Eye, eps_p
      dimension cm(*), hsv(*)
      integer :: hardening_type
      real(kind=8) bulkMod_kappa, shearMod_mu
c Material parameters
      !lame_lambda = cm_get('lame_lambda_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
      !bulkMod_kappa = lame_lambda + 2./3. * shearMod_mu
      bulkMod_kappa = cm_get('bulkMod_kappa___',cm)
c History variables
      eps_p = hsv_get_symTen2('eps_p', hsv)
c Second order identity tensor
      Eye = identity2(Eye)
c Volumetric part + deviatoric part
      get_stress_T_t = bulkMod_kappa * tr(sstrain) * Eye
     &                 + 2. * shearMod_mu * ( dev(sstrain) - eps_p )
c      
      end function get_stress_T_t
