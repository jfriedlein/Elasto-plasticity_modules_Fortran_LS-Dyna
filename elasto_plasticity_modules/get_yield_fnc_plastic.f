c
c
c
      real(kind=8) function get_yield_fnc_plastic( stress, alpha, cm,
     &                                             HillT_H_in, f4_in )
      ! @usage
      ! 1. For isotropic plasticity call
      !   Phi=get_yield_fnc_plastic( stress, alpha, cm )
      ! 2. For isotropic plasticity-damage call
      !   Phi=get_yield_fnc_plastic( stress, alpha, cm, f4_in=f4 )
      ! to skip the "HillT_H" as input and only provide "f4".
      ! 3. For anisotropic plasticity call
      !   Phi=get_yield_fnc_plastic( stress, alpha, cm, HillT_H )
      ! 4. For anisotropic plasticity-damage call
      !   Phi=get_yield_fnc_plastic( stress, alpha, cm, HillT_H, f4 )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      use enumerator_module
c
      implicit none
c
      type(Tensor2) :: stress, Eye
      type(Tensor4) :: HillT_H
      type(Tensor4), optional :: HillT_H_in
      real(kind=8) alpha
      real, dimension(*) :: cm
      integer :: anisotropy_type
      real,optional :: f4_in
      real(kind=8) :: f4, sigma_eff
c
      if ( present(HillT_H_in)) then
        HillT_H = HillT_H_in
      else
        HillT_H =  deviatoric_I4(Eye)
      endif
c
      if ( present(f4_in) ) then
          f4=f4_in
      else
          f4=1
      endif
c
      anisotropy_type = int(cm_get('anisotropy______',cm))
c
      if ( anisotropy_type == enum_P_iso
     &     .OR.
     &     anisotropy_type == enum_P_aniso_Hill48 ) then
        sigma_eff = get_yielding_norm( stress, HillT_H )
      elseif ( floor(anisotropy_type/10.) == enum_P_aniso_Yld91 ) then
        sigma_eff = get_stress_eff_Yld91(stress,cm)
      else
        write(*,*) 'get_stress_eff_Yld91<<
     &Provided anisotropy_type not defined'
        stop
      endif
c
      ! Plastic yield function
       get_yield_fnc_plastic =
     &      1./f4 * sigma_eff
     &      - sqrt(2./3.)
     &        * (
     &            cm_get('yieldStress_____',cm)
     &            - get_hardeningStress_R( alpha, cm )
     &          )
c      
      end function get_yield_fnc_plastic
