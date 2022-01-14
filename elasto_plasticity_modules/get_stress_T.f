c
c
c
      type(Tensor2) function get_stress_T( sstrain, gamma_k, n_k,
     &                                     cm_all, hsv,
     &                                     f1, f2, f4  )
c
      use Tensor
      use hsv_manager
      use cm_manager
      !implicit none
c
      type(Tensor2) :: sstrain, Eye, eps_p_n, n_k
      dimension cm_all(2,*), hsv(*)
      integer :: hardening_type
      real(kind=8) bulkMod_kappa, shearMod_mu
      real*8 :: f1, f2, f4
      real*8 gamma_k
c
c Material parameters
      shearMod_mu = cm_get_pair('shearMod_mu_____',cm_all)
      bulkMod_kappa = cm_get_pair('bulkMod_kappa___',cm_all)
c History variables
      eps_p_n = hsv_get_symTen2('eps_p', hsv)
c Second order identity tensor
      Eye = identity2(Eye)
c Volumetric part + deviatoric part
      get_stress_T = f1 * bulkMod_kappa * tr(sstrain) * Eye
     &               + 2. * f2 * shearMod_mu
     &                    * ( dev(sstrain)
     &                        - (eps_p_n + gamma_k/f4 * n_k) )
c      
      end function get_stress_T
