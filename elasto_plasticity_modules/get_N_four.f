c
c
c
      type(Tensor4) function get_N_four( stress, HillT_H, cm_all )
c
        use Tensor
        use cm_manager
c
        implicit none
c
        type(Tensor2) :: stress
        type(Tensor4) :: HillT_H
        real, dimension(2,*), intent(in) :: cm_all
        integer :: anisotropy_type
c        
        anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
c
        if ( anisotropy_type == enum_P_iso
     &      .OR.
     &      floor(anisotropy_type/10.) == enum_P_aniso_Hill ) then
         get_N_four = ( stress ** HillT_H ** stress)**(-1.5)
     &                * ( HillT_H * ( stress ** HillT_H ** stress )
     &                    - ((HillT_H**stress).dya.(stress**HillT_H)) )
        else
         write(*,*) 'get_N_four<<
     &N_four not implemented for the given anisotropy type',
     &anisotropy_type
         stop
       endif
c
      end function get_N_four
      
