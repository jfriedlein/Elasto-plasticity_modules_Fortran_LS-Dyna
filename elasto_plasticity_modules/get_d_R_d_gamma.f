c
c
c
      real(kind=8) function get_d_R_d_gamma( alpha_k, gamma_k,
     &  cm_all, hsv, crv, nnpcrv )
c
      use Tensor
      use cm_manager
      use hsv_manager
c
      include 'nlqparm'
c
      !implicit none
c
      real*8, intent(in) :: gamma_k
      real*8, intent(in) :: alpha_k
      real*8, dimension (2,*), intent(in) :: cm_all
      real*8, dimension (*), intent(in) :: hsv
      real*8, dimension(lq1,2,*), intent(in), optional :: crv
      integer, intent(in), optional :: nnpcrv(*)
      real*8 K, hardStress_R_inf, hardMod_K_exp,
     & hardExponent_n, expExponent_b, eid
c
c material parameters
      K = cm_get_pair('hardMod_K_______',cm_all)
      hardStress_R_inf = cm_get_pair('hardStress_R_inf',cm_all)
      hardMod_K_exp=cm_get_pair('hardMod_K_exp___',cm_all)
      hardExponent_n=cm_get_pair('hardExponent_n__',cm_all)
      expExponent_b=cm_get_pair('expExponent_b___',cm_all)
c
      select case( int(cm_get_pair('hardening_type__',cm_all)) )
        case( enum_hardening_linear ) ! linear hardening
         get_d_R_d_gamma = - K
     &                       * sqrt(2./3.)
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_d_R_d_gamma =  - K
     &                      * ( sqrt(2./3.)
     &						* (1. - K
     &                        /hardStress_R_inf
     &                          * hsv_get_scalar('alpha', hsv)) ! to get alpha_n
     &					 	/ (1. + sqrt(2./3.) * K
     &                 /hardStress_R_inf * gamma_k)**2 )
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_d_R_d_gamma = - K
     &                        * exp( -K / hardStress_R_inf * alpha_k )
     &                        * sqrt(2./3.)
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_d_R_d_gamma = (- K
     &                    - hardStress_R_inf
     &                      * hardMod_K_exp
     &                    * exp( -hardMod_K_exp
     &                           * alpha_k ))
     &                    * sqrt(2./3.)
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
         get_d_R_d_gamma = (-sqrt(2./3.)) *K
     &                     * hardMod_K_exp
     &                     * ( alpha_k+0.01 )
     &                     ** (hardMod_K_exp-1)
        case( enum_hardening_potExp ) ! exponent and exponential
         get_d_R_d_gamma = (
     &                        - K
     &                          * hardExponent_n
     &                          * ( alpha+1e-3 )
     &                            **( hardExponent_n-1 )
     &                        - hardStress_R_inf
     &                          * hardMod_K_exp
     &                          * exp( -hardMod_K_exp * alpha_k )
     &                     ) * sqrt(2./3.)
        case( enum_hardening_linExpExp )
         get_d_R_d_gamma = (
     &                    - K
     &                    - hardStress_R_inf * hardMod_K_exp
     &                      * exp( -hardMod_K_exp
     &                              * alpha_k**expExponent_b )
     &                      * expExponent_b*(alpha_k+1e-20)
     &                                       **(-1.+expExponent_b)
     &                     ) * sqrt(2./3.)
        case( enum_hardening_loadCurve )
         eid = cm_get_pair('loadCurve_ID____',cm_all)
         call crvval(crv,nnpcrv,eid,alpha_k,yval,get_d_R_d_gamma)
         get_d_R_d_gamma = (-1.) * get_d_R_d_gamma * sqrt(2./3.)
        case default
          write( *, * ) 'elpl-module-get_d_R_d_gamma<< 
     &undefined hardening type'
      end select 
c      
      end function get_d_R_d_gamma
