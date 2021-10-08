!DEC$ IF .NOT. DEFINED (elasto_plasticity_modules_F)
!DEC$ DEFINE elasto_plasticity_modules_F
! -----------MODULE elasto-plasticity modules---------------------------
      module elasto_plasticity_modules

      use Tensor
      use TensorXLSDYNA
      use hsv_manager
c @todo add a docu
c
      integer, parameter :: enum_hardening_linear=0,
     &                      enum_hardening_saturatedAlpha=1,
     &                      enum_hardening_Voce=2,
     &                      enum_hardening_linExp=3,
     &                      enum_hardening_Swift=4,
     &                      enum_hardening_potExp=5
c
      contains
c New hardening laws need to be added in:
c * get_hardeningStress_R
c * get_d_R_d_gamma
c * get_intVar_alpha
c      
!      ------BEGIN FUNCTIONS-------------------------------------
        include './get_stress_T_t.f'
        include './get_d_R_d_gamma.f'
        include './get_dPhi_dgamma.f'
        include './get_hardeningStress_R.f'
        include './get_intVar_alpha.f'
        include './get_plastic_yield_fnc.f'
        include './get_stress_k.f'
        include './get_tangent_C.f'
        include './get_tangent_elastic.f'
        include './get_yielding_norm.f'
        include './update_direction_n_n1.f'
        include './setup_Hill_tensor.f'
c
      end module elasto_plasticity_modules
!DEC$ ENDIF
