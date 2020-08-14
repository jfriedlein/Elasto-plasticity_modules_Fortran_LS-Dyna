c
c
c
      double precision function get_plastic_yield_fnc( stress, HillT_H,
     & alpha, hardening_type, cm )
c
      use Tensor
      use cm_manager
      !implicit none
c
      type(Tensor2) :: stress
      type(Tensor4) :: HillT_H
      dimension cm(*)
      integer :: hardening_type
c Plastic yield function
      get_plastic_yield_fnc = get_yielding_norm( stress, HillT_H )
     &            - sqrt(2./3.) * ( cm_get('yieldStress_____',cm)
     &          - get_hardeningStress_R( alpha, hardening_type, cm ) )
c      
      end function get_plastic_yield_fnc
