c
c
c
      double precision function get_d_R_d_gamma( alpha_k, gamma_k,
     &  hardening_type, cm ,hsv )
c
      use Tensor
      use cm_manager
      use hsv_manager
      !implicit none
c
      double precision gamma_k, alpha_k
      dimension cm(*), hsv(*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( 0 ) ! linear hardening
         get_d_R_d_gamma = - cm_get('hardMod_K_______',cm)
     &                       * sqrt(2./3.)
        case( 1 ) ! saturated alpha
         get_d_R_d_gamma =  - cm_get('hardMod_K_______',cm)
     &                      * ( sqrt(2./3.)
     &						* (1. - cm_get('hardMod_K_______',cm)
     &                        /cm_get('hardStress_R_inf',cm)
     &                          * hsv_get_scalar('alpha', hsv)) ! to get alpha_n
     &					 	/ (1. + sqrt(2./3.) * cm_get('hardMod_K_______',cm)
     &                 /cm_get('hardStress_R_inf',cm) * gamma_k)**2 )
        case( 2 ) ! saturated_Voce_hard_stress
         get_d_R_d_gamma = - cm_get('hardMod_K_______',cm)
     &                        * exp( -cm_get('hardMod_K_______',cm)
     &                        /cm_get('hardStress_R_inf',cm) * alpha_k )
     &                        * sqrt(2./3.)
        case( 3 ) ! saturated_Miehe_hard_stress
         get_d_R_d_gamma = (- cm_get('hardMod_K_______',cm)
     &                    - cm_get('hardStress_R_inf',cm)
     &                      * cm_get('hardMod_K_exp___',cm)
     &                    * exp( -cm_get('hardMod_K_exp___',cm)
     &                           * alpha_k ))
     &                    * sqrt(2./3.)
        case default
          write( *, * ) 'elpl-module-get_d_R_d_gamma<< 
     &undefined hardening type'
      end select 
c      
      end function get_d_R_d_gamma
