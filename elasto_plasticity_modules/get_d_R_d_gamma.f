c
c
c
      real(kind=8) function get_d_R_d_gamma( alpha_k, gamma_k,
     &  cm_all, hsv, crv, nnpcrv )
c
c @todo Add abs(alpha) to all hardening laws (what about tabular data?)
c
      use Tensor
      use cm_manager
      use hsv_manager
c
      include 'nlqparm'
c
      !implicit none
c
      real*8, intent(in) :: gamma_k ! only needed for a special hardening type
      real*8, intent(in) :: alpha_k
      real*8, dimension (2,*), intent(in) :: cm_all
      real*8, dimension (*), intent(in) :: hsv ! only needed for a special hardening type
      real*8, dimension(lq1,2,*), intent(in), optional :: crv
      integer, intent(in), optional :: nnpcrv(*)
      real*8 K, hardStress_R_inf, hardMod_K_exp,
     & hardExponent_n, expExponent_b, eid
      real*8, dimension(3) :: k1_i, k2_i
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
         write(*,*) "get_d_R_d_gamma<< saturatedAlpha hardening
     & is deactivated currently to avoid the need for gamma and hsv
     & as input. This is only relevant for AG-routines."
         call cstop('Not supported feature')
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
c     
        case( enum_hardening_VoceTriple )
            k1_i(1) = cm_get_pair('backStr1_K1_____',cm_all)
            k1_i(2) = cm_get_pair('backStr2_K1_____',cm_all)
            k1_i(3) = cm_get_pair('backStr3_K1_____',cm_all)
            k2_i(1) = cm_get_pair('backStr1_K2_____',cm_all)
            k2_i(2) = cm_get_pair('backStr2_K2_____',cm_all)
            k2_i(3) = cm_get_pair('backStr3_K2_____',cm_all)
c
      get_d_R_d_gamma = ( - 3./2. * k1_i(1)
     &                           * exp(- sqrt(3./2.) * k2_i(1)
     &                                  * alpha)
     &                    - 3./2. * k1_i(2)
     &                            * exp(- sqrt(3./2.) * k2_i(2)
     &                                   * alpha)
     &                    - 3./2. * k1_i(3)
     &                            * exp(- sqrt(3./2.) * k2_i(3)
     &                                   * alpha)      
     &                  ) * sqrt(2./3.)
c
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
         ! @note abs(alpha) is just used to protect the exponent from negative 
         ! values from memory errors, that would result in NaN            
         get_d_R_d_gamma = (
     &                    - K
     &                    - hardStress_R_inf * hardMod_K_exp
     &                      * exp( -hardMod_K_exp
     &                              * abs(alpha_k)**expExponent_b )
     &                      * expExponent_b*(abs(alpha_k)+1e-20)
     &                                       **(-1.+expExponent_b)
     &                     ) * sqrt(2./3.)
c     
        case ( enum_hardening_BM2013 ) !Bröcker and Matzenmiller 2013
         get_d_R_d_gamma = - (
     &                        exp(-K/hardStress_R_inf * alpha)
     &                        + hardMod_K_exp * ( K/hardStress_R_inf
     &                                            * alpha )
     &                                             **(expExponent_b-1.)
     &                                        * expExponent_b
     &                     ) * K * sqrt(2./3.)
      !write(*,*) "get_d_R_d_gamma",get_d_R_d_gamma
c
        case( enum_hardening_loadCurve )
         eid = cm_get_pair('loadCurve_ID____',cm_all)
         call crvval(crv,nnpcrv,eid,alpha_k,yval,get_d_R_d_gamma)
         get_d_R_d_gamma = (-1.) * get_d_R_d_gamma * sqrt(2./3.)
        case default
          write( *, * ) 'elpl-module-get_d_R_d_gamma<< 
     &undefined hardening type'
          call cstop('E R R O R  T E R M I N A T I O N!')
      end select 
c      
      end function get_d_R_d_gamma
