c
c
c
      subroutine setup_Hill_tensors( cm_all, HillT_H, HillT_H_s )
c
      use Tensor
      use cm_manager
      !implicit none
c
      real*8, dimension(2,*), intent(in) :: cm_all
      type(Tensor4), intent(out) :: HillT_H
      type(Tensor4), optional :: HillT_H_s
c     
      integer anisotropy_type
      logical Paniso_NonAssoc
      real*8, dimension(6) :: h_s_ij
c
      HillT_H = setup_Hill_tensor(cm_all)
c      
      anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
      Paniso_NonAssoc = ( anisotropy_type
     &                     == enum_P_aniso_Hill48_NonAssoc )
      if ( Paniso_NonAssoc ) then
         h_s_ij(1)=cm_get_pair('aniso_coeff_11_2',cm_all)
         h_s_ij(2)=cm_get_pair('aniso_coeff_22_2',cm_all)
         h_s_ij(3)=cm_get_pair('aniso_coeff_33_2',cm_all)
         h_s_ij(4)=cm_get_pair('aniso_coeff_12_2',cm_all)
         h_s_ij(5)=cm_get_pair('aniso_coeff_23_2',cm_all)
         h_s_ij(6)=cm_get_pair('aniso_coeff_31_2',cm_all)
         HillT_H_s = setup_Hill_tensor(cm_all, h_s_ij)
      else
         HillT_H_s = HillT_H
      endif
c
      end subroutine setup_Hill_tensors
