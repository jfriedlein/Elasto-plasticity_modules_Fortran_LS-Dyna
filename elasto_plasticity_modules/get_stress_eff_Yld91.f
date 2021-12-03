c
      real(kind=8) function get_stress_eff_Yld91( stress, cm_all )
      ! @note the effective stress is in accordance with the following yield function
      !  Phi = stress_eff - sqrt(2/3)*flow_stress
      ! meaning that the effective stress is scaled by sqrt(2/3)
      ! and not directly comparable with the von Mises stress 
      !  sigma_vM = sqrt(3/2)*||dev(stress)||
      ! but the deviatoric stress norm
      !  sigma_iso_eff = ||dev(stress)||
c
      use Tensor
      use cm_manager
      implicit none
c
      type(Tensor2) stress
      real, dimension(2,*) :: cm_all
      real(kind=8) :: J2tilde, J3tilde
      integer anisotropy_type
c
      anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
c
      ! Compute J2 and J3 for Yld91 according to [Cazacu 2019]
       J2tilde = get_Yld91_J2( stress, cm_all )
       J3tilde = get_Yld91_J3( stress, cm_all )
c
      ! Compute the effective stress according to [Cazacu 2019]
      ! Choose either FCC (m=8) or BCC (m=6)
      ! For instance if "enum_P_aniso_Yld91=2", a value of
      ! "anisotropy_type=28"=FCC and "anisotropy_type=26"=BCC
       if ( anisotropy_type==enum_P_aniso_Yld91FCC ) then
         ! FCC
          get_stress_eff_Yld91 = (
     &                             129.*J2tilde**4
     &                             - 324.*J2tilde*J3tilde**2
     &                           )**(1./8.)
       elseif ( anisotropy_type==enum_P_aniso_Yld91BCC ) then
         ! BCC
          get_stress_eff_Yld91 = (
     &                             33.*J2tilde**3 - 81/2.*J3tilde**2
     &                           )**(1./6.)
       else
         write(*,*) 'get_stress_eff_Yld91<<
     &Provided anisotropy_type not defined'
         stop
       endif
c
      ! Scale the effective stress by sqrt(2/3) for compatibility with
      ! our yield function (see above note)
       get_stress_eff_Yld91 = sqrt(2./3.) * get_stress_eff_Yld91
c
      end function get_stress_eff_Yld91
      
