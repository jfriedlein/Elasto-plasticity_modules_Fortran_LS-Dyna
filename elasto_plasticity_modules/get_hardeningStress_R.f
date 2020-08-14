c
c
c
      double precision function get_hardeningStress_R( alpha,
     & hardening_type, cm )
c
      use Tensor
      use cm_manager
      !implicit none
c
      double precision alpha
      dimension cm(*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( 0 ) ! linear hardening
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
        case( 1 ) ! saturated alpha
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
        case( 2 ) ! saturated_Voce_hard_stress
         get_hardeningStress_R = - cm_get('hardStress_R_inf',cm)
     &                       * ( 1
     &                            - exp(- cm_get('hardMod_K_______',cm)
     &                                    /cm_get('hardStress_R_inf',cm)
     &                                  * alpha) )
        case( 3 ) ! saturated_Miehe_hard_stress
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
     &                           - cm_get('hardStress_R_inf',cm)
     &                         * ( 1
     &                             - exp(-cm_get('hardMod_K_exp___',cm)
     &                               * alpha) )
        case default
          write( *, * ) 'elpl-module-get_hardeningStress_R<< 
     &undefined hardening type'
      end select 
c      
      end function get_hardeningStress_R
