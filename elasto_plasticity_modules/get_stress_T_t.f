c
c
c
      type(Tensor2) function get_stress_T_t( sstrain, cm_all, hsv )
c
      use Tensor
      use hsv_manager
      use cm_manager
      !implicit none
c
      type(Tensor2) :: sstrain, Eye, eps_p
      dimension cm_all(2,*), hsv(*)
      integer :: hardening_type
      real(kind=8) bulkMod_kappa, shearMod_mu
c Material parameters
      shearMod_mu = cm_get_pair('shearMod_mu_____',cm_all)
      bulkMod_kappa = cm_get_pair('bulkMod_kappa___',cm_all)
      hardening_type = cm_get_pair('hardening_type__',cm_all)
c History variables
      eps_p = hsv_get_symTen2('eps_p', hsv)
c Second order identity tensor
      Eye = identity2(Eye)
c Volumetric part + deviatoric part
      get_stress_T_t = bulkMod_kappa * tr(sstrain) * Eye
     &                 + 2. * shearMod_mu * ( dev(sstrain) - eps_p )
c      
      end function get_stress_T_t
