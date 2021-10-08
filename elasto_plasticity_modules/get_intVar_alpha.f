c
c
c
      double precision function get_intVar_alpha( alpha_n, gamma_k,
     & hardening_type, cm )
c
      use Tensor
      use cm_manager
      !implicit none
c
      double precision alpha_n
      real gamma_k
      dimension cm(*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( enum_hardening_linear ) ! linear hardening
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_intVar_alpha = (alpha_n + sqrt(2./3.) * gamma_k)
     &                     / ( 1. + sqrt(2./3.)
     &                         * cm_get('hardMod_K_______',cm)
     &                         /cm_get('hardStress_R_inf',cm) * gamma_k)
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( enum_hardening_potExp ) ! exponent and exponential
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case default
          write( *, * ) 'elpl-module-get_intVar_alpha<< 
     &Undefined hardening type'
      end select 
c      
      end function get_intVar_alpha
