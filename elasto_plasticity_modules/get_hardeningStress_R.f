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
      real(kind=8) alpha
      dimension cm(*)
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( enum_hardening_linear ) ! linear hardening
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_hardeningStress_R = - cm_get('hardStress_R_inf',cm)
     &                       * ( 1
     &                            - exp(- cm_get('hardMod_K_______',cm)
     &                                    /cm_get('hardStress_R_inf',cm)
     &                                  * alpha) )
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_hardeningStress_R = - cm_get('hardMod_K_______',cm) * alpha
     &                           - cm_get('hardStress_R_inf',cm)
     &                         * ( 1
     &                             - exp(-cm_get('hardMod_K_exp___',cm)
     &                               * alpha) )
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
          get_hardeningStress_R = cm_get('yieldStress_____',cm)
     &                            - cm_get('hardMod_K_______',cm)
     &                              * ( alpha+0.01 )
     &                                **cm_get('hardMod_K_exp___',cm)
        case( enum_hardening_potExp ) ! exponent and exponential
         get_hardeningStress_R = -cm_get('hardMod_K_______',cm)
     &                            * ( alpha+1e-3 )
     &                              **cm_get('hardExponent_n__',cm)
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
