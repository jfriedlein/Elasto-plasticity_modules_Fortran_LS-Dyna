c
      subroutine store_kinematic_hardening( gamma_k, n_k, cm_all, hsv )
c
      implicit none      
c         
      real*8, intent(in) :: gamma_k
      type(Tensor2), intent(in) :: n_k
      real, dimension(2,*), intent(in) :: cm_all
      real*8, dimension(*), intent(inout) :: hsv
c        
      integer hardening_kinematic
      integer i_k, n_back_stresses
      real*8, dimension(:), allocatable :: k_1_i, k_2_i
      type(Tensor2), dimension(:), allocatable :: B_i_n, B_i_k
c
      hardening_kinematic = INT(cm_get_pair('kin_hard_type___',cm_all))
c
      ! @todo currently we update the B-history after the tangent,
      ! because B_n is used for dB_gamma in the tangent, I am not
      ! sure this is correct/needed
      if ( hardening_kinematic==enum_kinHard_CR ) then
           call init_kinematic_hardening(
     &                           !input->
     &                            hardening_kinematic, cm_all, hsv,
     &                           !output->
     &                            B_i_n, n_back_stresses, k_1_i, k_2_i )
           allocate( B_i_k(n_back_stresses) )             

           do i_k=1,n_back_stresses
             B_i_k(i_k) = 1./(1.+k_2_i(i_k)*gamma_k)
     &                  * ( B_i_n(i_k) + gamma_k * n_k * k_1_i(i_k) )
           enddo
c           
           if ( n_back_stresses >= 1 ) then
              call hsv_set_symTen2(B_i_k(1),'B_1__', hsv)
           endif
           if ( n_back_stresses >= 2 ) then
              call hsv_set_symTen2(B_i_k(2),'B_2__', hsv)
           endif
           if ( n_back_stresses == 3) then              
              call hsv_set_symTen2(B_i_k(3),'B_3__', hsv)
           endif
c           
      endif ! kinHard

      end subroutine store_kinematic_hardening