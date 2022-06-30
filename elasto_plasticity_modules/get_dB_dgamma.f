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
      type(Tensor2), dimension(:), allocatable :: B_i_n
      real*8, intent(in) :: gamma_k
      real*8, dimension(*), intent(in) :: hsv
      real*8, dimension(2,*), intent(in) :: cm_all
      real*8 :: a_omega
      real*8, dimension(:), allocatable :: k_1_i, k_2_i
      integer hardening_kinematic, n_back_stresses
      integer i
c Material parameters and history
       hardening_kinematic = INT(cm_get_pair('kin_hard_type___',cm_all))
       call init_kinematic_hardening(
     &                           !input->
     &                            hardening_kinematic, cm_all, hsv,
     &                           !output->
     &                            B_i_n, n_back_stresses, k_1_i, k_2_i )  
c
      a_omega=0.
      do i=1,n_back_stresses
            a_omega = a_omega + k_1_i(i) / (1.+k_2_i(i)*gamma_k)
      enddo
c
      get_dB_dgamma = 0.
      do i=1,n_back_stresses
            get_dB_dgamma = get_dB_dgamma
     &                      + (- k_2_i(i) / (1.+k_2_i(i)*gamma_k)**2)
     &                        * B_i_n(i)
      enddo
c
      get_dB_dgamma = get_dB_dgamma + n_k * a_omega
c
      do i=1,n_back_stresses
            get_dB_dgamma = get_dB_dgamma
     &                      + n_k * gamma_k
     &                        * k_1_i(i)*k_2_i(i)
     &                          / (1.+k_2_i(i)*gamma_k)**2
      enddo
c
      end function get_dB_dgamma
