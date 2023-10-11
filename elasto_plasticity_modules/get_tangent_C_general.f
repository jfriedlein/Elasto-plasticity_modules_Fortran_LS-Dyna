c
c
c
      type(Tensor4) function get_tangent_C_general( stress, E_four,
     &                                              alpha_k, n_r_k,
     &                                              gamma_k, cm_all,
     &                                              hsv, crv, nnpcrv,
     &                                              n_s_k_in )
c
      use Tensor
      use cm_manager
      use TensorXkinematics
      !implicit none
c
      include 'nlqparm' ! needed for parameter "lq1"
c
      type(Tensor2) :: stress, n_r_k, n_s_k, Eye, dB_gamma
      type(Tensor2), optional :: n_s_k_in
      type(Tensor4) :: E_four
      real*8, dimension(2,*), intent(in) :: cm_all
      real*8, intent(in), optional :: crv(lq1,2,*)
      integer, intent(in) :: nnpcrv(*)
      real*8, dimension(*), intent(in) ::  hsv
      real(kind=8) alpha_k, gamma_k
      integer hardening_type
      logical hardening_kinematic
c
      if ( present(n_s_k_in) ) then
        n_s_k = n_s_k_in
      else
        n_s_k = n_r_k
      endif
c
      hardening_type = int(cm_get_pair('hardening_type__',cm_all))
      hardening_kinematic = INT(cm_get_pair('kin_hard_type___',cm_all))
c Second order identity tensor
      Eye = identity2(Eye)
c kinematic hardening contribution
      if ( hardening_kinematic.NE.enum_kinHard_OFF ) then
        dB_gamma = get_dB_dgamma( n_r_k, gamma_k, cm_all, hsv )
      endif
c
        get_tangent_C_general = E_four
     &              - 1. / ( n_s_k ** E_four ** n_r_k
     &                       + n_s_k ** dB_gamma
     &                       - sqrt( 2./3. )
     &                         * get_d_R_d_gamma( alpha_k, gamma_k,
     &                                       cm_all, hsv, crv, nnpcrv )
     &                     ) * ( (E_four**n_r_k).dya.(n_s_k**E_four) )
c      
      end function get_tangent_C_general
