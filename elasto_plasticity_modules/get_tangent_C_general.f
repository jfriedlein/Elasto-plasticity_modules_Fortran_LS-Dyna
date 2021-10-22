c
c
c
      type(Tensor4) function get_tangent_C_general( stress, E_four,
     &                                              alpha_k, n_r_k,
     &                                              gamma_k, cm,
     &                                              hsv, n_s_k_in )
c
      use Tensor
      use cm_manager
      use TensorXkinematics
      !implicit none
c
      type(Tensor2) :: stress, n_r_k, n_s_k, Eye
      type(Tensor2), optional :: n_s_k_in
      type(Tensor4) :: E_four
      real, dimension(*) :: cm, hsv
      real(kind=8) alpha_k, gamma_k
      integer hardening_type
c
      if ( present(n_s_k_in) ) then
        n_s_k = n_s_k_in
      else
        n_s_k = n_r_k
      endif
c
      hardening_type = int(cm_get('hardening_type__',cm))
c Second order identity tensor
      Eye = identity2(Eye)
c
        get_tangent_C_general = E_four
     &              - 1. / ( n_s_k ** E_four ** n_r_k - sqrt( 2./3. )
     &                      * get_d_R_d_gamma( alpha_k, gamma_k,
     &                                       hardening_type, cm, hsv ) )
     &                      * ( (E_four**n_r_k).dya.(n_s_k**E_four) )
c      
      end function get_tangent_C_general
