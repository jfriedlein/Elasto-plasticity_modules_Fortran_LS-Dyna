c
c
c
      type(Tensor2) function get_dB_dgamma( n_k, gamma_k, cm, hsv )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      implicit none
c
      type(Tensor2), intent(in) :: n_k
      type(Tensor2), dimension(2) :: B_i_n
      real*8, intent(in) :: gamma_k
      real*8, dimension(*), intent(in) :: cm, hsv
      real*8 :: a_omega
      real*8, dimension(2) :: k_1_i, k_2_i
c Material parameters
      k_1_i(1)=0.!2666.6667!5000.!77000
      k_1_i(2)=0.
      k_2_i(1)=0.!100.!409.75
      k_2_i(2)=0.
c history
      B_i_n(1) = hsv_get_symTen2('B_1__', hsv)
      B_i_n(2) = hsv_get_symTen2('B_2__', hsv)
c
      a_omega = k_1_i(1) / (1.+k_2_i(1)*gamma_k)
     &          + k_1_i(2) / (1.+k_2_i(2)*gamma_k)
c
      get_dB_dgamma = (- k_2_i(1) / (1.+k_2_i(1)*gamma_k)**2) * B_i_n(1)
     &              + (- k_2_i(2) / (1.+k_2_i(2)*gamma_k)**2) * B_i_n(2)
     &              + n_k * (
     &                          a_omega
     &                          + gamma_k * (
     &                                        k_1_i(1)*k_2_i(1)
     &                                        / (1.+k_2_i(1)*gamma_k)**2
     &                                        + k_1_i(2)*k_2_i(2)
     &                                        / (1.+k_2_i(2)*gamma_k)**2
     &                                       )
     &                       )
c      
      end function get_dB_dgamma
