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
     &       hardExponent_n, expExponent_b, eid
      real*8, dimension(3) :: k1_i, k2_i
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
c     
        case( enum_hardening_VoceTriple )
            if (  cm_get_pair('kin_hard_type___',cm_all)
     &            == enum_kinHard_CR ) then
                  write(*,*) "get_flow_stress<< Trying to use 
     &the Voce-triple 
     &isotropic hardening with active kinematic hardening. The 
     &former option cannot be used combined with the latter, 
     &because it uses the parameters of the latter."
                  call cstop('E R R O R  T E R M I N A T I O N!')
            endif
c
            k1_i(1) = cm_get_pair('backStr1_K1_____',cm_all)
            k1_i(2) = cm_get_pair('backStr2_K1_____',cm_all)
            k1_i(3) = cm_get_pair('backStr3_K1_____',cm_all)
            k2_i(1) = cm_get_pair('backStr1_K2_____',cm_all)
            k2_i(2) = cm_get_pair('backStr2_K2_____',cm_all)
            k2_i(3) = cm_get_pair('backStr3_K2_____',cm_all)
c
            get_flow_stress = yieldStress
     &                     + sqrt(3./2.) * k1_i(1)/k2_i(1)
     &                       * ( 1
     &                            - exp(- sqrt(3./2.) * k2_i(1)
     &                                  * alpha) )
     &                     +  sqrt(3./2.)* k1_i(2)/k2_i(2)
     &                       * ( 1
     &                            - exp(- sqrt(3./2.) * k2_i(2)
     &                                  * alpha) )
     &                     +  sqrt(3./2.) * k1_i(3)/k2_i(3)
     &                       * ( 1
     &                            - exp(- sqrt(3./2.) * k2_i(3)
     &                                  * alpha) )    
c           
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
c     
        case ( enum_hardening_BM2013 ) !Bröcker and Matzenmiller 2013
        ! yieldStress = kappa_0
        ! hardStress_R_inf = kappa^inf
        ! K = E_kappa
        ! hardMod_K_exp = alpha_kappa
        ! expExponent_b = m_kappa                    
         get_flow_stress = yieldStress
     &                     + hardStress_R_inf
     &                       * ( 1. - exp(-K/hardStress_R_inf * alpha)
     &                           + hardMod_K_exp * ( K/hardStress_R_inf
     &                                               * alpha )
     &                                             **expExponent_b
     &                         )
      !write(*,*) "get_flow_stress",get_flow_stress
c     
        case( enum_hardening_loadCurve )
         eid = cm_get_pair('loadCurve_ID____',cm_all)
         call crvval(crv,nnpcrv,eid,alpha, get_flow_stress, slope)
        case default
          write( *, * ) 'elpl-module-get_flow_stress<<
     &undefined hardening type'
          call cstop('E R R O R  T E R M I N A T I O N')
      end select 
c      
      end function get_flow_stress
