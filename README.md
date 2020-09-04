# Elasto-plasticity_modules_LS-Dyna
A module containing elasto-plastic material models (Hill-Plasticity) with straightforward extension to various hardening laws. 

(also available for C++ library deal.II [here](https://github.com/jfriedlein/Elasto-plasticity_modules_dealii))

## What it does
We offer a framework to capture elasto-plastic material models up to anisotropic Hill-plasticity with various hardening laws. The framework can be found in the exemplary `MaterialModel.h` file and also contains subiterations on the qp level. Different hardenig laws defined by the hardening stress `R` and an evolution equation for the internal hardening variable `alpha` can be defined by only three equations (R, alpha, d_R_d_gamma in `elpl_equation_list.h`). The algorithm is general enough to produce quadratic convergence for any such defined hardening law. Currently, only isotropic hardening is supported, an extension for kinematic hardening might follow in the future.

@todo note on low efficiency, recommended for testing

@todo add the enumerator_list or a section of it

@todo ensure that it is standalone

@todo add support for Sacado (for arbitrary evolution equations, etc)

@todo add a verification example for anisotropy (either from external source for own verifi, or plate with a hole for internal verifi)

DOCU still missing!!!