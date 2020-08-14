      module elasto_plasticity_modules
      use Tensor
      use TensorXLSDYNA
      use hsv_manager
c @todo add a docu

      contains

!      ------BEGIN FUNCTIONS-------------------------------------
        include './get_stress_T_t.f'
        include './get_d_R_d_gamma.f'
        include './get_dPhi_dgamma.f'
        include './get_hardeningStress_R.f'
        include './get_intVar_alpha.f'
        include './get_plastic_yield_fnc.f'
        include './get_stress_k.f'
        include './get_tangent_C.f'
        include './get_yielding_norm.f'
        include './update_direction_n_n1.f'
c
      end module elasto_plasticity_modules
