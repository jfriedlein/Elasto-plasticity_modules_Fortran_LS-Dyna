c
c
      real(kind=8) function get_Yld91_J3 ( stress, cm )
c
      use Tensor
      use cm_manager
c      
      implicit none      
c
      type(Tensor2) stress
      real, dimension(*) :: cm
      real(kind=8) :: a, b, c, f, g, h,
     &                sigmaXX, sigmaYY, sigmaZZ, 
     &                sigmaXY, sigmaYZ, sigmaXZ
c      
      ! Get the anisotropy coefficients
       a=cm_get('HillCoeff_h11___',cm)
       b=cm_get('HillCoeff_h22___',cm)
       c=cm_get('HillCoeff_h33___',cm)
       h=cm_get('HillCoeff_h12___',cm)
       f=cm_get('HillCoeff_h23___',cm)
       g=cm_get('HillCoeff_h31___',cm)
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
      
