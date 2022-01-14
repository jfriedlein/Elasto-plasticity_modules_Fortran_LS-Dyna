c
c
c
      type(Tensor2) function get_stress_T_t( sstrain, cm_all, hsv,
     &                                       f1_in, f2_in  )
c
      use Tensor
      use hsv_manager
      use cm_manager
      !implicit none
c
      type(Tensor2) :: sstrain, Eye, eps_p_n
      dimension cm_all(2,*), hsv(*)
      real(kind=8) bulkMod_kappa, shearMod_mu
      real*8 :: f1, f2
      real(kind=8), optional :: f1_in, f2_in
c
      if ( present(f1_in) ) then
          f1=f1_in
      else
          f1=1.
      endif
c
      if ( present(f2_in) ) then
          f2=f2_in
      ! If only f1 is present but not f2,
      ! we use f1 as f for both f1 and f2 (f=f1=f2)
      elseif ( present(f1_in) ) then
          f2=f1
      else
          f2=1.
      endif
c
c Material parameters
      shearMod_mu = cm_get_pair('shearMod_mu_____',cm_all)
      bulkMod_kappa = cm_get_pair('bulkMod_kappa___',cm_all)
c History variables
      eps_p_n = hsv_get_symTen2('eps_p', hsv)
c Second order identity tensor
      Eye = identity2(Eye)
c Volumetric part + deviatoric part
      get_stress_T_t = f1 * bulkMod_kappa * tr(sstrain) * Eye
     &                 + 2.*f2*shearMod_mu * ( dev(sstrain) - eps_p_n )
c      
      end function get_stress_T_t
