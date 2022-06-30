c
c
c
      real(kind=8) function get_intVar_alpha( alpha_n, gamma_k,
     &  cm_all )
c
      use Tensor
      use cm_manager
      !implicit none
c
      double precision alpha_n
      real gamma_k
      dimension cm_all(2,*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( enum_hardening_linear ) ! linear hardening
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_intVar_alpha = (alpha_n + sqrt(2./3.) * gamma_k)
     &                     / ( 1. + sqrt(2./3.)
     &                         * cm_get_pair('hardMod_K_______',cm_all)
     &                         / cm_get_pair('hardStress_R_inf',cm_all)
     &                           * gamma_k)
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_potExp ) ! exponent and exponential
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_linExpExp )
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case ( enum_hardening_BM2013 ) !Bröcker and Matzenmiller 2013
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_loadCurve )
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case default
          write( *, * ) 'elpl-module-get_intVar_alpha<< 
     &Undefined hardening type'
      end select 
c      
      end function get_intVar_alpha
