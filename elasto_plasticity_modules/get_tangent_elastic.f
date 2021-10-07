c
c
c
      type(Tensor4) function get_tangent_elastic( cm, f1_in, f2_in )
c
      use Tensor
      use cm_manager
      use TensorXkinematics
      !implicit none
c
      type(Tensor2) :: Eye
      dimension cm(*)
      real(kind=8), optional :: f1_in, f2_in
      real(kind=8) :: f1, f2
c Material parameters
      shearMod_mu = cm_get('shearMod_mu_____',cm)
      bulkMod_kappa = cm_get('bulkMod_kappa___',cm)
c
      ! if "f1" and "f2" are input, then use them
       if ( present(f1_in) .and. present(f2_in) ) then
          f1=f1_in
          f2=f2_in
      ! in case only one of them is given, use them for both
      elseif ( present(f1_in) ) then
          f1=f1_in
          f2=f1
      ! without damage input, we use the standard elastic tangent
      else
          f1=1.
          f2=1.
      endif
c Second order identity tensor
      Eye = identity2(Eye)
c
      get_tangent_elastic = f1 * bulkMod_kappa * (Eye.dya.Eye)
     &                      + 2. * f2 * shearMod_mu * deviatoric_I4(Eye)
c      
      end function get_tangent_elastic
