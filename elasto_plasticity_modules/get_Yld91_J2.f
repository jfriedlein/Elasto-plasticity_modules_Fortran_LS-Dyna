c
c
      real(kind=8) function get_Yld91_J2 ( stress, cm )
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
      ! Compute J2 for Yld91 according to [Cazacu 2019]
       get_Yld91_J2 = f**2*sigmaYZ**2 + g**2*sigmaXZ**2 +
     & h**2*sigmaXY**2 + ((b + c)*sigmaXX - c*sigmaYY - b*sigmaZZ)**2./
     & 18. + ((-c)*sigmaXX + (c + a)*sigmaYY - a*sigmaZZ)**2./
     & 18. + ((-b)*sigmaXX - a*sigmaYY + (a + b)*sigmaZZ)**2./18.
c     
      end function get_Yld91_J2
      