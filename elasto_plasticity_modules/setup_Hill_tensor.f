c
c
c
      type(Tensor4) function setup_Hill_tensor( cm, h_ij_in )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      !implicit none
c
      real*8, dimension(*), intent(in) :: cm
      real*8, dimension(6), optional :: h_ij_in
      real*8, dimension(9) :: alpha_
      real*8, dimension(6,6) :: H_matrix
      type(Tensor1) :: a1, a2, a3
      type(Tensor1) :: dum1, dum2
      type(Tensor2), dimension(3,3) :: m_
      type(Tensor2) :: H_I, I_H, Eye
      real*8 :: sheetOrientation_theta_rad
      integer :: i,j,k,l
      integer :: anisotropy_type
      real*8, dimension(6) :: h_ij
c
c Check whether the Hill-coefficients have been provided as input.
      ! If they are input, use them ...
       if ( present(h_ij_in) ) then
          h_ij=h_ij_in
      ! ... else retrieve them from the material parameters "cm"
       else
          h_ij(1)=cm_get('HillCoeff_h11___',cm)
          h_ij(2)=cm_get('HillCoeff_h22___',cm)
          h_ij(3)=cm_get('HillCoeff_h33___',cm)
          h_ij(4)=cm_get('HillCoeff_h12___',cm)
          h_ij(5)=cm_get('HillCoeff_h23___',cm)
          h_ij(6)=cm_get('HillCoeff_h31___',cm)
       endif
c
      Eye = identity2(Eye)
c Retrieve the anisotropy type (iso, Hill, Yld91, ...)
      anisotropy_type = int(cm_get('anisotropy______',cm))
c Build the anisotropy tensor "setup_Hill_tensor"
      ! Isotropic
       if ( anisotropy_type == enum_P_iso  ) then
         setup_Hill_tensor = deviatoric_I4(Eye)
c
      ! Anisotropic Hill
       elseif ( anisotropy_type == enum_P_aniso_Hill48 ) then
	     ! set up the hill tensor based on the hill coefficients
	     ! Following the paper "Anisotropic additive plasticity in the logarithmic strain space"
	     ! by Miehe et al. eq. (3.40)-(3.46)
	      alpha_(1) = 2.d0/3.d0 * h_ij(1)**(-2)
	      alpha_(2) = 2.d0/3.d0 * h_ij(2)**(-2)
	      alpha_(3) = 2.d0/3.d0 * h_ij(3)**(-2)
	      alpha_(7) = 1.d0/3.d0 * h_ij(4)**(-2)
	      alpha_(8) = 1.d0/3.d0 * h_ij(5)**(-2)
	      alpha_(9) = 1.d0/3.d0 * h_ij(6)**(-2)
	      alpha_(4) = 0.5d0 * ( alpha_(3) - alpha_(1) - alpha_(2) )
	      alpha_(5) = 0.5d0 * ( alpha_(1) - alpha_(2) - alpha_(3) )
	      alpha_(6) = 0.5d0 * ( alpha_(2) - alpha_(1) - alpha_(3) )

	     ! The Hill tensor in matrix representation
          H_matrix(1,1) = alpha_(1)
          H_matrix(2,2) = alpha_(2)
          H_matrix(3,3) = alpha_(3)
          H_matrix(4,4) = 0.5d0 * alpha_(7)
          H_matrix(5,5) = 0.5d0 * alpha_(8)
          H_matrix(6,6) = 0.5d0 * alpha_(9)

          H_matrix(1,2) = alpha_(4)
          H_matrix(2,1) = H_matrix(1,2)

          H_matrix(2,3) = alpha_(5)
          H_matrix(3,2) = H_matrix(2,3)

          H_matrix(1,3) = alpha_(6)
          H_matrix(3,1) = H_matrix(1,3)
c
	     ! orthogonal basis by three orthogonal directions a_i
	     ! Retrieve the sheet orientation to the rolling direction
          sheetOrientation_theta_rad = cm_get('sheetOrientation',cm)
     &                                 / 180. * 4. * atan(1.)
		 ! first basis vector (for theta=0 equal to x-axis)
		  a1%a(1) = cos( sheetOrientation_theta_rad )
		  a1%a(2) = sin( sheetOrientation_theta_rad )
		  a1%a(3) = 0.d0
		 ! second basis vector (for theta=0 equal to y-axis)
		  a2%a(1) = - sin( sheetOrientation_theta_rad )
		  a2%a(2) =   cos( sheetOrientation_theta_rad )
		  a2%a(3) = 0.d0
		 ! third basis vector (for sheets always along z-axis)
          a3%a(1) = 0.d0
          a3%a(2) = 0.d0
		  a3%a(3) = 1.d0
c
         ! @todo Why does the a_(i) stuff not work?
         !      forall( i=1:3, j=1:3 )
         !&        m_(i,j) = 0.5 * ( (a_(i).dya.a_(j)) + (a_(j).dya.a_(i)) )
          m_(1,1) = 0.5d0 * ( (a1.dya.a1) + (a1.dya.a1) )
          m_(2,2) = 0.5d0 * ( (a2.dya.a2) + (a2.dya.a2) )
          m_(3,3) = 0.5d0 * ( (a3.dya.a3) + (a3.dya.a3) )
          m_(1,2) = 0.5d0 * ( (a1.dya.a2) + (a2.dya.a1) )
          m_(2,3) = 0.5d0 * ( (a2.dya.a3) + (a3.dya.a2) )
          m_(3,1) = 0.5d0 * ( (a3.dya.a1) + (a1.dya.a3) )
          m_(3,2) = m_(2,3)
          m_(2,1) = m_(1,2)
          m_(1,3) = m_(3,1)
c
	     ! @todo: What about the goofy factor of 2 for isotropy???
	       setup_Hill_tensor =
     &                 alpha_(1) * ( m_(1,1).dya.m_(1,1) )
     &		       + alpha_(2) * ( m_(2,2).dya.m_(2,2) )
     &		       + alpha_(3) * ( m_(3,3).dya.m_(3,3) )
     &		       + alpha_(4) * 0.5d0 * ( ( m_(1,1).dya.m_(2,2) )
     &                                + ( m_(2,2).dya.m_(1,1) ) ) *2.d0 !factor of 2?????
     &		       + alpha_(5) * 0.5d0 * ( ( m_(2,2).dya.m_(3,3) )
     &                                + ( m_(3,3).dya.m_(2,2) ) ) *2.d0 !factor of 2?????
     &		       + alpha_(6) * 0.5d0 * ( ( m_(1,1).dya.m_(3,3) )
     &                                + ( m_(3,3).dya.m_(1,1) ) ) *2.d0 !factor of 2?????
     &		       + alpha_(7) * 2.d0 * ( m_(1,2).dya.m_(2,1) )
     &		       + alpha_(8) * 2.d0 * ( m_(2,3).dya.m_(3,2) )
     &		       + alpha_(9) * 2.d0 * ( m_(1,3).dya.m_(3,1) )
c             
          ! Check whether the Hill tensor is purely deviatoric
           H_I = setup_Hill_tensor**Eye
           I_H = Eye**setup_Hill_tensor
           if ( ( norm(H_I) + norm(I_H)) > 1e-14 ) then
             write(*,*) "setup_Hill_tensor<< 
     &Hill Tensor not purely deviatoric,
     &wrong setup of the equations. Results in ",
     &norm(H_I)," instead of less than 1e-14 
     &(numercially zero). Hill Tensor not purely deviatoric."
              stop
	       endif
c
	  ! For Yld91-type anisotropy, we don't use the Hill tensor,
	  ! so we set it to zero to make sure that we receive an error,
	  ! whenever we use it at all
	   elseif ( floor(anisotropy_type/10.) == enum_P_aniso_Yld91 ) then
	     setup_Hill_tensor = 0.
	   else
         write(*,*) 'setup_Hill_tensor<<
     &Provided anisotropy_type not defined'
         stop
       endif
c      
      end function setup_Hill_tensor
