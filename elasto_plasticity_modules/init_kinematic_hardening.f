c      
      subroutine init_kinematic_hardening(
     &                           !input->
     &                            hardening_kinematic, cm_all, hsv,
     &                           !output->
     &                            B_i_n, n_back_stresses,
     &                            k_1_i_out, k_2_i_out )
c
      use Tensor
      use cm_manager
      use hsv_manager
c      
      implicit none
c        
      integer, intent(in) :: hardening_kinematic
      real*8, dimension(2,*), intent(in) :: cm_all
      real*8, dimension(*), intent(in) ::  hsv
      integer, intent(out), optional :: n_back_stresses
      integer :: NBS
      type(Tensor2), dimension(:), allocatable, intent(out) :: B_i_n
      real*8, dimension(:), allocatable, optional, intent(out) 
     & :: k_1_i_out, k_2_i_out
      real*8, dimension(:), allocatable :: k_1_i, k_2_i
c
       select case( hardening_kinematic )
         case ( enum_kinHard_OFF ) ! no kinematic hardening
            allocate( B_i_n(0) )   
c             
            if ( present(n_back_stresses) ) then
                n_back_stresses = 0
            endif
c            
            if ( present(k_1_i_out) ) then
                k_1_i_out(:)=0.
            endif
c            
            if ( present(k_2_i_out) ) then
                k_2_i_out(:)=0.
            endif
c            
         case ( enum_kinHard_CR )
            NBS = cm_get_pair('n_back_stresses_',cm_all)
c
            if ( NBS < 1 .OR. NBS > 3 ) then
                write(*,*) "init_kinematic_hardening<<
     & If you choose kinematic hardening, please use 1, 2 or 3 back
     & stresses, not ",NBS
                call cstop ('E R R O R  T E R M I N A T I O N')
            endif
c            
            allocate( k_1_i(NBS) )
            allocate( k_2_i(NBS) )
            allocate( B_i_n(NBS) )    
c            
            if ( NBS >= 1 ) then
                k_1_i(1) = cm_get_pair('backStr1_K1_____',cm_all)
                k_2_i(1) = cm_get_pair('backStr1_K2_____',cm_all)
                B_i_n(1) = hsv_get_symTen2('B_1__', hsv)
            endif
            if ( NBS >= 2 ) then
                k_1_i(2) = cm_get_pair('backStr2_K1_____',cm_all)
                k_2_i(2) = cm_get_pair('backStr2_K2_____',cm_all)
                B_i_n(2) = hsv_get_symTen2('B_2__', hsv)
            endif
            if ( NBS == 3) then
              k_1_i(3) = cm_get_pair('backStr3_K1_____',cm_all)
              k_2_i(3) = cm_get_pair('backStr3_K2_____',cm_all)
              B_i_n(3) = hsv_get_symTen2('B_3__', hsv)
            endif
c            
            ! Only return NBS if requested. Also note that the variable
            ! "n_back_stresses" only exists if it was entered. That's why
            ! we save the data into "NBS", which always exists.
             if ( present(n_back_stresses) ) then
                n_back_stresses = NBS
             endif
c             
            if ( present(k_1_i_out) ) then
                allocate( k_1_i_out(NBS) )
                k_1_i_out(:)=k_1_i
            endif
c            
            if ( present(k_2_i_out) ) then
                allocate( k_2_i_out(NBS) )
                k_2_i_out(:)=k_2_i
            endif
c
        case default
            write( *, * ) 'init_kinematic_hardening<<
     & undefined kinematic hardening type'
            call cstop('E R R O R  T E R M I N A T I O N')
        end select ! hardening_kinematic
c
      return
      end subroutine init_kinematic_hardening       