c
c
      real(kind=8) function get_Yld91_J3 ( stress, cm_all )
c
      use Tensor
      use cm_manager
c      
      implicit none      
c
      type(Tensor2) stress
      real, dimension(2,*) :: cm_all
      real(kind=8) :: a, b, c, f, g, h,
     &                sigmaXX, sigmaYY, sigmaZZ, 
     &                sigmaXY, sigmaYZ, sigmaXZ
c      
      ! Get the anisotropy coefficients
       a=cm_get_pair('aniso_coeff_11__',cm_all)
       b=cm_get_pair('aniso_coeff_22__',cm_all)
       c=cm_get_pair('aniso_coeff_33__',cm_all)
       h=cm_get_pair('aniso_coeff_12__',cm_all)
       f=cm_get_pair('aniso_coeff_23__',cm_all)
       g=cm_get_pair('aniso_coeff_31__',cm_all)
c      
      ! Save the stress components into separate variable for easier use 
       sigmaXX=stress%ab(1,1)
       sigmaYY=stress%ab(2,2)
       sigmaZZ=stress%ab(3,3)
       sigmaXY=stress%ab(1,2)
       sigmaYZ=stress%ab(2,3)
       sigmaXZ=stress%ab(1,3)
c
      ! Compute J3 for Yld91 according to [Cazacu 2019]
       get_Yld91_J3 = 2.*f*g*h*sigmaXY*sigmaXZ*sigmaYZ +
     &  1./27.*(((b + c)*sigmaXX - c*sigmaYY -
     & b*sigmaZZ)*(-c*sigmaXX + (c + a)*sigmaYY -
     &  a*sigmaZZ)*(-b*sigmaXX - a*sigmaYY + (a + b)*sigmaZZ)) - (f**2*
     & sigmaYZ**2)/3.
     & *((b + c)*sigmaXX - c*sigmaYY - b*sigmaZZ) - (g**2*sigmaXZ**2)/3.
     & *(-c*sigmaXX + (c + a)*sigmaYY - a*sigmaZZ) - (h**2*sigmaXY**2)/
     &    3.*(-b*sigmaXX - a*sigmaYY + (a + b)*sigmaZZ)
c     
      end function get_Yld91_J3
      
