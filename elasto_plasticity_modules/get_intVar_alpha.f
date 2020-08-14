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
      double precision alpha_n, gamma_k
      dimension cm(*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( 0 ) ! linear hardening
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( 1 ) ! saturated alpha
         get_intVar_alpha = (alpha_n + sqrt(2./3.) * gamma_k)
     &                     / ( 1. + sqrt(2./3.)
     &                         * cm_get('hardMod_K_______',cm)
     &                         /cm_get('hardStress_R_inf',cm) * gamma_k)
        case( 2 ) ! saturated_Voce_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case( 3 ) ! saturated_Miehe_hard_stress
         get_intVar_alpha = alpha_n + sqrt(2./3.) * gamma_k
        case default
          write( *, * ) 'elpl-module-get_intVar_alpha<< 
     &Undefined hardening type'
      end select 
c      
      end function get_intVar_alpha
