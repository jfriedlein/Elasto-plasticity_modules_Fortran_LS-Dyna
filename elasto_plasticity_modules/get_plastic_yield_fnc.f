c
c
c
      double precision function get_plastic_yield_fnc( stress, HillT_H,
     & alpha, hardening_type, cm, f4_in )
c
      use Tensor
      use cm_manager
      !implicit none
c
      type(Tensor2) :: stress
      type(Tensor4) :: HillT_H
      real(kind=8) alpha
      dimension cm(*)
      integer :: hardening_type
      real,optional :: f4_in
      real :: f4
c
      if ( present(f4_in) ) then
          f4=f4_in
      else
          f4=1
      endif
c Plastic yield function
      get_plastic_yield_fnc =
     &    1./f4 * get_yielding_norm( stress, HillT_H )
     &    - sqrt(2./3.)
     &      * ( 
     &            cm_get('yieldStress_____',cm)
     &            - get_hardeningStress_R( alpha, cm )
     &        )
c
      write(*,*) "get_plastic_yield_fnc<< 
     & OoO use get_yield_fnc_plastic(*) instead"
      stop
c
c      
      end function get_plastic_yield_fnc
