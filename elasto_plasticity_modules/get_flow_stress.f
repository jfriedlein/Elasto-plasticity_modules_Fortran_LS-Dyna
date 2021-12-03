c
c
c
      real(kind=8) function get_flow_stress( alpha, cm_all,
     &                                       crv, nnpcrv )
c
      use Tensor
      use cm_manager
c
      include 'nlqparm'
c
      !implicit none
c
      real(kind=8), intent(in) :: alpha
      real*8, intent(in) :: cm_all(2,*)
      real*8, intent(in), optional :: crv(lq1,2,*)
      integer, intent(in), optional :: nnpcrv(*)
      real*8 yieldStress, K, hardStress_R_inf, hardMod_K_exp,
     & hardExponent_n, expExponent_b, eid
c
      yieldStress=cm_get_pair('yieldStress_____',cm_all)
      K = cm_get_pair('hardMod_K_______',cm_all)
      hardStress_R_inf = cm_get_pair('hardStress_R_inf',cm_all)
      hardMod_K_exp=cm_get_pair('hardMod_K_exp___',cm_all)
      hardExponent_n=cm_get_pair('hardExponent_n__',cm_all)
      expExponent_b=cm_get_pair('expExponent_b___',cm_all)
c
c Material parameters
      select case( int(cm_get_pair('hardening_type__',cm_all)) )
        case( enum_hardening_linear ) ! linear hardening
         get_flow_stress = yieldStress + K * alpha
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_flow_stress = yieldStress + K * alpha
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_flow_stress = yieldStress + hardStress_R_inf
     &                       * ( 1
     &                            - exp(- K / hardStress_R_inf
     &                                  * alpha) )
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_flow_stress = yieldStress + K * alpha
     &                           + hardStress_R_inf
     &                         * ( 1.
     &                             - exp(-hardMod_K_exp
     &                               * alpha) )
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
          get_flow_stress =  K * ( alpha+0.01 )**hardMod_K_exp
        case( enum_hardening_potExp ) ! exponent and exponential
         get_flow_stress = yieldStress + K
     &                            * ( alpha+1e-3 )
     &                              **hardExponent_n
     &                           + hardStress_R_inf
     &                         * ( 1.
     &                             - exp(-hardMod_K_exp
     &                               * alpha) )
        case( enum_hardening_linExpExp )
         get_flow_stress = yieldStress + K * alpha
     &                           + hardStress_R_inf
     &                             * ( 1.
     &                                 - exp( -hardMod_K_exp
     &                                        * alpha**expExponent_b) )
        case( enum_hardening_loadCurve )
         eid = cm_get_pair('loadCurve_ID____',cm_all)
         call crvval(crv,nnpcrv,eid,alpha, get_flow_stress, slope)
        case default
          write( *, * ) 'elpl-module-get_flow_stress<<
     &undefined hardening type'
      end select 
c      
      end function get_flow_stress
