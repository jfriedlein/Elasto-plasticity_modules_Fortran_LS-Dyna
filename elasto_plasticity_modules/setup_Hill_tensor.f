c
c
c
      type(Tensor4) function setup_Hill_tensor( cm )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      !implicit none
c
      double precision, dimension(9) :: alpha_
      dimension cm(*)
      real, dimension(6,6) :: H_matrix
      type(Tensor1) :: a1, a2, a3
      type(Tensor1) :: dum1, dum2
      type(Tensor2), dimension(3,3) :: m_
      type(Tensor2) :: H_I, I_H, Eye
      real :: sheetOrientation_theta_rad
      integer :: i,j,k
      real anisotropy_active
c
      Eye = identity2(Eye)
c
      anisotropy_active = cm_get('anisotropy______',cm)
      ! Isotropic
      if ( anisotropy_active < 0.5 ) then
            setup_Hill_tensor = deviatoric_I4(Eye)
      ! Anisotropic
      elseif ( anisotropy_active > 0.5 ) then
	     ! set up the hill tensor based on the hill coefficients
	     ! Following the paper "Anisotropic additive plasticity in the logarithmic strain space"
	     ! by Miehe et al. eq. (3.40)-(3.46)
	     alpha_(1) = 2./3. * cm_get('HillCoeff_h11___',cm)**(-2)
	     alpha_(2) = 2./3. * cm_get('HillCoeff_h22___',cm)**(-2)
	     alpha_(3) = 2./3. * cm_get('HillCoeff_h33___',cm)**(-2)
	     alpha_(7) = 1./3. * cm_get('HillCoeff_h12___',cm)**(-2)
	     alpha_(8) = 1./3. * cm_get('HillCoeff_h23___',cm)**(-2)
	     alpha_(9) = 1./3. * cm_get('HillCoeff_h31___',cm)**(-2)
	     alpha_(4) = 0.5 * ( alpha_(3) - alpha_(1) - alpha_(2) )
	     alpha_(5) = 0.5 * ( alpha_(1) - alpha_(2) - alpha_(3) )
	     alpha_(6) = 0.5 * ( alpha_(2) - alpha_(1) - alpha_(3) )

	     ! The Hill tensor in matrix representation
            H_matrix(1,1) = alpha_(1)
            H_matrix(2,2) = alpha_(2)
            H_matrix(3,3) = alpha_(3)
            H_matrix(4,4) = 0.5 * alpha_(7)
            H_matrix(5,5) = 0.5 * alpha_(8)
            H_matrix(6,6) = 0.5 * alpha_(9)

            H_matrix(1,2) = alpha_(4)
            H_matrix(2,1) = H_matrix(1,2)

            H_matrix(2,3) = alpha_(5)
            H_matrix(3,2) = H_matrix(2,3)

            H_matrix(1,3) = alpha_(6)
            H_matrix(3,1) = H_matrix(1,3)

	     ! orthogonal basis by three orthogonal directions a_i
            sheetOrientation_theta_rad = cm_get('sheetOrientation',cm)
     &                                   /180. * 4. * atan(1.)
		     ! first basis vector (for theta=0° equal to x-axis)
		      a1%a(1) = cos( sheetOrientation_theta_rad )
		      a1%a(2) = sin( sheetOrientation_theta_rad )
		      a1%a(3) = 0.
		     ! second basis vector (for theta=0° equal to y-axis)
		      a2%a(1) = - sin( sheetOrientation_theta_rad )
		      a2%a(2) =   cos( sheetOrientation_theta_rad )
		      a2%a(3) = 0.
		     ! third basis vector (for sheets always along z-axis)
                a3%a(1) = 0.
                a3%a(2) = 0.
		      a3%a(3) = 1.

              ! @todo Why does the a_(i) stuff not work?
         !      forall( i=1:3, j=1:3 )
         !&        m_(i,j) = 0.5 * ( (a_(i).dya.a_(j)) + (a_(j).dya.a_(i)) )
                m_(1,1) = 0.5 * ( (a1.dya.a1) + (a1.dya.a1) )
                m_(2,2) = 0.5 * ( (a2.dya.a2) + (a2.dya.a2) )
                m_(3,3) = 0.5 * ( (a3.dya.a3) + (a3.dya.a3) )
                m_(1,2) = 0.5 * ( (a1.dya.a2) + (a2.dya.a1) )
                m_(2,1) = m_(1,2)
                m_(2,3) = 0.5 * ( (a2.dya.a3) + (a3.dya.a2) )
                m_(3,2) = m_(2,3)
                m_(3,1) = 0.5 * ( (a3.dya.a1) + (a1.dya.a3) )
                m_(1,3) = m_(3,1)

	     ! @todo: What about the goofy factor of 2 for isotropic???
	       setup_Hill_tensor = alpha_(1) * ( m_(1,1).dya.m_(1,1) )
     &		       + alpha_(2) * ( m_(2,2).dya.m_(2,2) )
     &		       + alpha_(3) * ( m_(3,3).dya.m_(3,3) )
     &		       + alpha_(4) * 0.5 * ( ( m_(1,1).dya.m_(2,2) )
     &                                + ( m_(2,2).dya.m_(1,1) ) ) * 2. !factor of 2?????
     &		       + alpha_(5) * 0.5 * ( ( m_(2,2).dya.m_(3,3) )
     &                                + ( m_(3,3).dya.m_(2,2) ) ) * 2. !factor of 2?????
     &		       + alpha_(6) * 0.5 * ( ( m_(1,1).dya.m_(3,3) )
     &                                + ( m_(3,3).dya.m_(1,1) ) ) * 2. !factor of 2?????
     &		       + alpha_(7) * 2. * ( m_(1,2).dya.m_(2,1) )
     &		       + alpha_(8) * 2. * ( m_(2,3).dya.m_(3,2) )
     &		       + alpha_(9) * 2. * ( m_(1,3).dya.m_(3,1) )

              H_I = setup_Hill_tensor**Eye
              I_H = Eye**setup_Hill_tensor
             
             if ( ( norm(H_I) + norm(I_H)) > 1e-14 ) then
                write(*,*) "HillT_H<< Hill Tensor not purely deviatoric,
     &wrong setup of the equations. Results in ",
     &norm(H_I)," instead of less than 1e-14 
     &(numercially zero). Hill Tensor not purely deviatoric"
                pause
	        endif
      endif
      
      end function setup_Hill_tensor