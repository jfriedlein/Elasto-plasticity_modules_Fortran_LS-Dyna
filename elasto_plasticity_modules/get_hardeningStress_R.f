c
c
c
      real(kind=8) function get_hardeningStress_R( alpha, cm_all,
     &                                             crv, nnpcrv )
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
      real*8 K, hardStress_R_inf, hardMod_K_exp,
     & hardExponent_n, expExponent_b, eid
c
      K = cm_get_pair('hardMod_K_______',cm_all)
      hardStress_R_inf = cm_get_pair('hardStress_R_inf',cm_all)
      hardMod_K_exp=cm_get_pair('hardMod_K_exp___',cm_all)
      hardExponent_n=cm_get_pair('hardExponent_n__',cm_all)
      expExponent_b=cm_get_pair('expExponent_b___',cm_all)
c
      write(*,*) "get_hardeningStress_R<< deprecated"
      stop
c
c Material parameters
      select case( int(cm_get_pair('hardening_type__',cm_all)) )
        case( enum_hardening_linear ) ! linear hardening
         get_hardeningStress_R = - K * alpha
        case( enum_hardening_saturatedAlpha ) ! saturated alpha
         get_hardeningStress_R = - K * alpha
        case( enum_hardening_Voce ) ! saturated_Voce_hard_stress
         get_hardeningStress_R = - hardStress_R_inf
     &                       * ( 1
     &                            - exp(- K / hardStress_R_inf
     &                                  * alpha) )
        case( enum_hardening_linExp ) ! saturated_Miehe_hard_stress
         get_hardeningStress_R = - K * alpha
     &                           - hardStress_R_inf
     &                         * ( 1.
     &                             - exp(-hardMod_K_exp
     &                               * alpha) )
        case( enum_hardening_Swift ) ! Swift exponential hardening equations as in LS-Dyna (Manual Vol.II, *MAT_122_3D)
          get_hardeningStress_R = cm_get_pair('yieldStress_____',cm_all)
     &                            - K
     &                              * ( alpha+0.01 )
     &                                **hardMod_K_exp
        case( enum_hardening_potExp ) ! exponent and exponential
         get_hardeningStress_R = -K
     &                            * ( alpha+1e-3 )
     &                              **hardExponent_n
     &                           - hardStress_R_inf
     &                         * ( 1.
     &                             - exp(-hardMod_K_exp
     &                               * alpha) )
        case( enum_hardening_linExpExp )
         get_hardeningStress_R = - K * alpha
     &                           - hardStress_R_inf
     &                             * ( 1.
     &                                 - exp( -hardMod_K_exp
     &                                        * alpha**expExponent_b) )
        case default
          write( *, * ) 'elpl-module-get_hardeningStress_R<< 
     &undefined hardening type'
      end select 
c      
      end function get_hardeningStress_R
