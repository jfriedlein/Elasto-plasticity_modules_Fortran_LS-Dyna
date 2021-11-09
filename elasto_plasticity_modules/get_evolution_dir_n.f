c
c
c
      type(Tensor2) function get_evolution_dir_n( stress, cm,
     &                                            HillT_H_in )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      use enumerator_module
c
      implicit none
c
      type(Tensor2), intent(in) :: stress
      type(Tensor2) :: Eye
      real, dimension(*), intent(in) :: cm
      type(Tensor4), optional :: HillT_H_in
      type(Tensor4) :: HillT_H
      real(kind=8) :: a, b, c, f, g, h,
     &                sigmaXX, sigmaYY, sigmaZZ,
     &                sigmaXY, sigmaYZ, sigmaXZ
      real(kind=8) :: dTbardJ2tilde, dTbardJ3tilde,
     &                sigmaTild11, sigmaTild22, sigmaTild33,
     &                J2tilde, J3tilde,
     &      dJ3tildedsigmaXXAY, dJ3tildedsigmaZZAY, dJ3tildedsigmaYYAY,
     &      dJ3tildedsigmaXYAY, dJ3tildedsigmaYZAY, dJ3tildedsigmaXZAY,
     &      dJ2tildedsigmaXXAY, dJ2tildedsigmaYYAY, dJ2tildedsigmaZZAY,
     &      dJ2tildedsigmaXYAY, dJ2tildedsigmaXZAY, dJ2tildedsigmaYZAY
      integer anisotropy_type
c
c
      if ( present(HillT_H_in)) then
        HillT_H = HillT_H_in
      else
        HillT_H =  deviatoric_I4(Eye)
      endif
c 
      anisotropy_type = int(cm_get('anisotropy______',cm))
c
       if ( anisotropy_type == enum_P_iso
     &      .OR.
     &      anisotropy_type == enum_P_aniso_Hill48 ) then
         get_evolution_dir_n = 1. / get_yielding_norm( stress, HillT_H )
     &                         * (HillT_H .ddot. stress)
c
       elseif ( floor(anisotropy_type/10.) == enum_P_aniso_Yld91 ) then
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

      ! Compute J2 and J3 for Yld91 according to [Cazacu 2019]
       J2tilde = get_Yld91_J2( stress, cm )
       J3tilde = get_Yld91_J3( stress, cm )
      !#################################################################
       sigmaTild11 = 1./3.*((b + c)*sigmaXX - c*sigmaYY - b*sigmaZZ)
       sigmaTild22 = 1./3.*(-c*sigmaXX + (c + a)*sigmaYY - a*sigmaZZ)
       sigmaTild33 = 1./3.*((-b)*sigmaXX - a*sigmaYY + (a+b)*sigmaZZ)

        dJ2tildedsigmaXXAY = 1./3.* sigmaTild11*(b + c)
     &   + 1./3.*(-c)*sigmaTild22
     &   + 1./3.*(-b)*sigmaTild33

        dJ2tildedsigmaYYAY = (1./3.)*sigmaTild11*(-c)
     &   + (1./3.)*(a + c)*sigmaTild22
     &   + (1./3.)*(-a)*sigmaTild33

        dJ2tildedsigmaZZAY = (1./3.)*sigmaTild11*(-b)
     &   + (1./3.)*(-a)*sigmaTild22
     &   + (1./3.)*(a + b)*sigmaTild33

        dJ2tildedsigmaXYAY = h**2*sigmaXY
        dJ2tildedsigmaXZAY = g**2*sigmaXZ
        dJ2tildedsigmaYZAY = f**2*sigmaYZ
      !#################################################################
       dJ3tildedsigmaXXAY = (1./27.)*(9.*b*h**2*sigmaXY**2
     &  + 9.*c*g**2*sigmaXZ**2 -
     &   9.*(b + c)*f**2*sigmaYZ**2 +
     &  c*(c*(sigmaXX - sigmaYY) + b*(sigmaXX - sigmaZZ))*(b*sigmaXX +
     & a*sigmaYY - (a + b)*sigmaZZ) +
     & (b + c)*((-c)*sigmaXX + (a + c)*sigmaYY -
     & a*sigmaZZ)*((-b)*sigmaXX - a*sigmaYY + (a + b)*sigmaZZ) +
     &   b*(c*(sigmaXX - sigmaYY) +
     & b*(sigmaXX - sigmaZZ))*(c*(sigmaXX - sigmaYY) +
     & a*(-sigmaYY + sigmaZZ)))

        dJ3tildedsigmaYYAY = (1./27.)*((-c)*(9.*(g*sigmaXZ
     &   - f*sigmaYZ)*(g*sigmaXZ + f*sigmaYZ) +
     & b*(2.*c*(sigmaXX - sigmaYY) + b*(sigmaXX - sigmaZZ))*(sigmaXX -
     &    sigmaZZ)) +
     & a**2.*(sigmaYY - sigmaZZ)*(2*b*(-sigmaXX + sigmaZZ) -
     & c*(2.*sigmaXX - 3.*sigmaYY + sigmaZZ)) +
     &   a*(9.*(h*sigmaXY - g*sigmaXZ)*(h*sigmaXY + g*sigmaXZ) -
     & b**2*(sigmaXX - sigmaZZ)**2 +
     &    c**2*(sigmaXX - sigmaYY)*(sigmaXX - 3.*sigmaYY + 2.*sigmaZZ)))

      dJ3tildedsigmaZZAY = -dJ3tildedsigmaXXAY - dJ3tildedsigmaYYAY

        dJ3tildedsigmaXYAY=f*g*h*sigmaXZ*sigmaYZ - (1./3.)*
     & sigmaXY*((-((-a)*h**2 + (a + b)*h**2))*sigmaXX -
     & a*h**2*sigmaYY + (a + b)*h**2*sigmaZZ)

      dJ3tildedsigmaXZAY=f*g*h*sigmaXY*sigmaYZ - (1./3.)*
     &sigmaXZ*((-((-a)*g**2 + (a + c)*g**2))*sigmaXX + (a + c)*g**2*
     &sigmaYY - a*g**2*sigmaZZ)

      dJ3tildedsigmaYZAY=f*g*h*sigmaXY*sigmaXZ - (1./3.)*sigmaYZ*
     & ((b*f**2 + c*f**2)*sigmaXX - c*f**2*sigmaYY - b*f**2*sigmaZZ)
      !#################################################################
      ! Choose either FCC (m=8) or BCC (m=6)
      ! For instance if "enum_P_aniso_Yld91=2", a value of
      ! "anisotropy_type=28"=FCC and "anisotropy_type=26"=BCC
       if ( anisotropy_type==enum_P_aniso_Yld91*10+8 ) then
         ! FCC
          dTbardJ2tilde = ((129./2.)*J2tilde**3 - (81./2.)*J3tilde**2)
     &            /(129.*J2tilde**4 - 324.*J2tilde*J3tilde**2)**(7./8.)

          dTbardJ3tilde = (-81.*J2tilde*J3tilde)
     &            /(129.*J2tilde**4 - 324.*J2tilde*J3tilde**2)**(7./8.)
       elseif ( anisotropy_type==enum_P_aniso_Yld91*10+6 ) then
         ! BCC
          dTbardJ2tilde = ((33./2.)*J2tilde**2)/(33.*J2tilde**3
     &                      - (81./2.)*J3tilde**2)**(5./6.)

          dTbardJ3tilde = ((-27./2.)*J3tilde)/(33.*J2tilde**3
     &                      - (81./2.)*J3tilde**2)**(5./6.)
       else
         write(*,*) "get_stress_eff_Yld91<< 
     &Provided anisotropy_type not defined"
         stop
      endif
      !#################################################################

        get_evolution_dir_n%ab(1,1) = dTbardJ2tilde*dJ2tildedsigmaXXAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaXXAY
        get_evolution_dir_n%ab(2,2) = dTbardJ2tilde*dJ2tildedsigmaYYAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaYYAY
        get_evolution_dir_n%ab(3,3) = dTbardJ2tilde*dJ2tildedsigmaZZAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaZZAY
        get_evolution_dir_n%ab(1,2) = dTbardJ2tilde*dJ2tildedsigmaXYAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaXYAY
        get_evolution_dir_n%ab(2,3) = dTbardJ2tilde*dJ2tildedsigmaYZAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaYZAY
        get_evolution_dir_n%ab(3,1) = dTbardJ2tilde*dJ2tildedsigmaXZAY
     &                                + dTbardJ3tilde*dJ3tildedsigmaXZAY
        get_evolution_dir_n%ab(2,1)=get_evolution_dir_n%ab(1,2)
        get_evolution_dir_n%ab(3,2)=get_evolution_dir_n%ab(2,3)
        get_evolution_dir_n%ab(1,3)=get_evolution_dir_n%ab(3,1)
c
       ! Scale the effective stress by sqrt(2/3) for compatibility with
       ! our yield function (see above note)
        get_evolution_dir_n=sqrt(2./3.)*get_evolution_dir_n
c
       else ! anisotropy_type
         write(*,*) 'get_stress_eff_Yld91<<
     &Provided anisotropy_type not defined'
         stop
       endif ! anisotropy_type
c      
      end function get_evolution_dir_n
c
