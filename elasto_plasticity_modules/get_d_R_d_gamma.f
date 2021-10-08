c
c
c
      double precision function get_d_R_d_gamma( alpha_k, gamma_k,
     &  hardening_type, cm ,hsv )
      ! @todo Can we skip "hardening_type" as input and get it from "cm"
c
      use Tensor
      use cm_manager
      use hsv_manager
      !implicit none
c
      real gamma_k
      real(kind=8) alpha_k
      real, dimension (*) :: cm, hsv
      integer :: hardening_type
c Material parameters
      select case( hardening_type )
        case( enum_hardening_linear ) ! linear hardening
         get_d_R_d_gamma = - cm_get('hardMod_K_______',cm)
     &                       * sqrt(2./3.)
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_d_R_d_gamma =  - cm_get('hardMod_K_______',cm)
     &                      * ( sqrt(2./3.)
     &						* (1. - cm_get('hardMod_K_______',cm)
     &                        /cm_get('hardStress_R_inf',cm)
     &                          * hsv_get_scalar('alpha', hsv)) ! to get alpha_n
     &					 	/ (1. + sqrt(2./3.) * cm_get('hardMod_K_______',cm)
     &                 /cm_get('hardStress_R_inf',cm) * gamma_k)**2 )
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_d_R_d_gamma = - cm_get('hardMod_K_______',cm)
     &                        * exp( -cm_get('hardMod_K_______',cm)
     &                        /cm_get('hardStress_R_inf',cm) * alpha_k )
     &                        * sqrt(2./3.)
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_d_R_d_gamma = (- cm_get('hardMod_K_______',cm)
     &                    - cm_get('hardStress_R_inf',cm)
     &                      * cm_get('hardMod_K_exp___',cm)
     &                    * exp( -cm_get('hardMod_K_exp___',cm)
     &                           * alpha_k ))
     &                    * sqrt(2./3.)
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
         get_d_R_d_gamma = (-sqrt(2./3.)) *cm_get('hardMod_K_______',cm)
     &                     * cm_get('hardMod_K_exp___',cm)
     &                     * ( alpha_k+0.01 )
     &                     ** (cm_get('hardMod_K_exp___',cm)-1)
        case( enum_hardening_potExp ) ! exponent and exponential
         get_d_R_d_gamma = (
     &                        - cm_get('hardMod_K_______',cm)
     &                          * cm_get('hardExponent_n__',cm)
     &                          * ( alpha+1e-3 )
     &                            **( cm_get('hardExponent_n__',cm)-1 )
     &                        - cm_get('hardStress_R_inf',cm)
     &                          * cm_get('hardMod_K_exp___',cm)
     &                          * exp( -cm_get('hardMod_K_exp___',cm)
     &                                  * alpha_k )
     &                     ) * sqrt(2./3.)
        case default
          write( *, * ) 'elpl-module-get_d_R_d_gamma<< 
     &undefined hardening type'
      end select 
c      
      end function get_d_R_d_gamma
