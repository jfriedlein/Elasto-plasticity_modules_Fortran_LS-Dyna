c
c
c
      type(Tensor2) function get_dB_dgamma( n_k, gamma_k, cm_all, hsv )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      implicit none
c
      type(Tensor2), intent(in) :: n_k
      type(Tensor2), dimension(3) :: B_i_n
      real*8, intent(in) :: gamma_k
      real*8, dimension(*), intent(in) :: hsv
      real*8, dimension(2,*), intent(in) :: cm_all
      real*8 :: a_omega
      real*8, dimension(3) :: k_1_i, k_2_i
      integer n_back_stresses
c Material parameters
      n_back_stresses = cm_get_pair('n_back_stresses_',cm_all)
      k_1_i(1) = cm_get_pair('backStr1_K1_____',cm_all) !2666.6667!5000.!77000
      k_1_i(2) = cm_get_pair('backStr2_K1_____',cm_all)
      k_2_i(1) = cm_get_pair('backStr1_K2_____',cm_all) !100.!409.75
      k_2_i(2) = cm_get_pair('backStr2_K2_____',cm_all)
      if (n_back_stresses==3) then
            k_1_i(3) = cm_get_pair('backStr3_K1_____',cm_all)
            k_2_i(3) = cm_get_pair('backStr3_K2_____',cm_all)
      else
            k_1_i(3) = 0.
            k_2_i(3) = 0.
      endif
c history
      B_i_n(1) = hsv_get_symTen2('B_1__', hsv)
      B_i_n(2) = hsv_get_symTen2('B_2__', hsv)
      if (n_back_stresses==3) then
         B_i_n(3) = hsv_get_symTen2('B_3__', hsv)
      else
         B_i_n(3) = 0.
      endif
c
      a_omega = k_1_i(1) / (1.+k_2_i(1)*gamma_k)
     &          + k_1_i(2) / (1.+k_2_i(2)*gamma_k)
     &          + k_1_i(3) / (1.+k_2_i(3)*gamma_k)
c
      get_dB_dgamma = (- k_2_i(1) / (1.+k_2_i(1)*gamma_k)**2) * B_i_n(1)
     &              + (- k_2_i(2) / (1.+k_2_i(2)*gamma_k)**2) * B_i_n(2)
     &              + (- k_2_i(3) / (1.+k_2_i(3)*gamma_k)**2) * B_i_n(3)
     &              + n_k * (
     &                          a_omega
     &                          + gamma_k * (
     &                                        k_1_i(1)*k_2_i(1)
     &                                        / (1.+k_2_i(1)*gamma_k)**2
     &                                        + k_1_i(2)*k_2_i(2)
     &                                        / (1.+k_2_i(2)*gamma_k)**2
     &                                        + k_1_i(3)*k_2_i(3)
     &                                        / (1.+k_2_i(3)*gamma_k)**2
     &                                       )
     &                       )
c      
c
      end function get_dB_dgamma
