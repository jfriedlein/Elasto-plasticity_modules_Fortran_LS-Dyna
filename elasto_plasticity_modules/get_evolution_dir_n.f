c
c
c
      type(Tensor2) function get_evolution_dir_n( stress, cm_all,
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
      real*8, dimension(2,*), intent(in) :: cm_all
      type(Tensor4), optional :: HillT_H_in
      type(Tensor4) :: HillT_H
      real(kind=8) :: a, b, c, f, g, h,
     &                sigmaXX, sigmaYY, sigmaZZ,
     &                sigmaXY, sigmaYZ, sigmaXZ
      real(kind=8) :: dJ3tildedsigmaXXAY, dJ3tildedsigmaYYAY,
     &                dJ3tildedsigmaZZAY, dJ3tildedsigmaXYAY,
     &                dJ3tildedsigmaXZAY, dJ3tildedsigmaYZAY,
     &                J2tilde, J3tilde,
     &sigmaTild11,sigmaTild22,sigmaTild33,
     &dJ2tildedsigmaXXAY, dJ2tildedsigmaYYAY,
     &                dJ2tildedsigmaZZAY, dJ2tildedsigmaXYAY,
     &                dJ2tildedsigmaXZAY, dJ2tildedsigmaYZAY,
     & dTbardJ2tilde,dTbardJ3tilde
      integer anisotropy_type
c
c
      if ( present(HillT_H_in)) then
        HillT_H = HillT_H_in
      else
        HillT_H =  deviatoric_I4(Eye)
      endif
c 
      anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
c
       if ( anisotropy_type == enum_P_iso
     &      .OR.
     &      floor(anisotropy_type/10.) == enum_P_aniso_Hill ) then
         get_evolution_dir_n = 1. / get_yielding_norm( stress, HillT_H )
     &                         * (HillT_H .ddot. stress)
c
       elseif ( floor(anisotropy_type/10.) == enum_P_aniso_Yld91 ) then
c
      ! Get the anisotropy coefficients
       a=cm_get_pair('aniso_coeff_11__',cm_all)
       b=cm_get_pair('aniso_coeff_22__',cm_all)
       c=cm_get_pair('aniso_coeff_33__',cm_all)
       h=cm_get_pair('aniso_coeff_12__',cm_all)
       f=cm_get_pair('aniso_coeff_23__',cm_all)
       g=cm_get_pair('aniso_coeff_31__',cm_all)
c
      if ( a==1. .and. b==1. .and. c==1. .and.
     &     h==1. .and. f==1. .and. g==1. ) then
            !okay
      else
            write(*,*) "get_stress_eff_Yld91<< 
     &something wrong with anisotropic case. Do not use this.
     &Only usable for isotropy a=b=c=h=f=g=1"
         call cstop ('E R R O R  T E R M I N A T I O N')
      endif ! iso
       !write(*,*) "n=",get_evolution_dir_n
c
      ! Save the stress components into separate variable for easier use
       sigmaXX=stress%ab(1,1)
       sigmaYY=stress%ab(2,2)
       sigmaZZ=stress%ab(3,3)
       sigmaXY=stress%ab(1,2)
       sigmaYZ=stress%ab(2,3)
       sigmaXZ=stress%ab(1,3)

      ! Compute J2 and J3 for Yld91 according to [Cazacu 2019]
       J2tilde = get_Yld91_J2( stress, cm_all )
       J3tilde = get_Yld91_J3( stress, cm_all )
c
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

      ! #################################################################
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
       dJ3tildedsigmaXXAY = (1./27.)*(b*(9.*h**2*sigmaXY**2
     &  - 9.*f**2*sigmaYZ**2 + c**2*(sigmaXX - sigmaYY)
     &  *(3.*sigmaXX - sigmaYY - 2.*sigmaZZ)
     &  - a**2*(sigmaYY - sigmaZZ)**2) + 
     &    b**2*(sigmaXX - sigmaZZ)*(c*(3.*sigmaXX - 2*sigmaYY - sigmaZZ)
     &    + 2.*a*(-sigmaYY + sigmaZZ)) + c*(9.*g**2*sigmaXZ**2
     &    - 9.*f**2*sigmaYZ**2 - a*(sigmaYY - sigmaZZ)
     &    *(-2.*c*sigmaXX + a*sigmaYY + 2.*c*sigmaYY - a*sigmaZZ)))

        dJ3tildedsigmaYYAY = (1./27.)*((-c)*(9.*(g*sigmaXZ - f*sigmaYZ)
     &   *(g*sigmaXZ + f*sigmaYZ) + b*(2.*c*(sigmaXX - sigmaYY) 
     &   + b*(sigmaXX - sigmaZZ))*(sigmaXX - sigmaZZ)) + 
     & a**2*(sigmaYY - sigmaZZ)*(2.*b*(-sigmaXX + sigmaZZ) 
     & - c*(2.*sigmaXX - 3.*sigmaYY + sigmaZZ)) 
     & + a*(9.*(h*sigmaXY - g*sigmaXZ)*(h*sigmaXY + g*sigmaXZ)
     & - b**2*(sigmaXX - sigmaZZ)**2 + 
     & c**2*(sigmaXX - sigmaYY)*(sigmaXX - 3.*sigmaYY + 2.*sigmaZZ)))

      dJ3tildedsigmaZZAY = -dJ3tildedsigmaXXAY - dJ3tildedsigmaYYAY

        dJ3tildedsigmaXYAY=(1./3.)*h*(3.*f*g*sigmaXZ*sigmaYZ 
     &   + h*sigmaXY*(b*sigmaXX + a*sigmaYY - (a + b)*sigmaZZ))

      dJ3tildedsigmaXZAY=(1./3.)*g*(c*g*sigmaXZ*(sigmaXX - sigmaYY)
     & + 3.*f*h*sigmaXY*sigmaYZ + a*g*sigmaXZ*(-sigmaYY + sigmaZZ))

      dJ3tildedsigmaYZAY=f*g*h*sigmaXY*sigmaXZ
     & + (1./3.)*f**2*sigmaYZ*(c*(-sigmaXX + sigmaYY)
     & + b*(-sigmaXX + sigmaZZ))
      !#################################################################
      ! Choose either FCC (m=8) or BCC (m=6)
      ! For instance if "enum_P_aniso_Yld91=2", a value of
      ! "anisotropy_type=28"=FCC and "anisotropy_type=26"=BCC
       if ( anisotropy_type==enum_P_aniso_Yld91FCC ) then
         ! FCC
          dTbardJ2tilde = ((129./2.)*J2tilde**3 - (81./2.)*J3tilde**2)
     &            /(129.*J2tilde**4 - 324.*J2tilde*J3tilde**2)**(7./8.)

          dTbardJ3tilde = (-81.*J2tilde*J3tilde)
     &            /(129.*J2tilde**4 - 324.*J2tilde*J3tilde**2)**(7./8.)
       elseif ( anisotropy_type==enum_P_aniso_Yld91BCC ) then
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
         call cstop ('E R R O R  T E R M I N A T I O N')
       endif ! anisotropy_type
c      
      end function get_evolution_dir_n
c


      subroutine get_evolution_dirs_n( stress, cm_all,
     &                                 n_k, HillT_H_in, 
     &                                 HillT_H_s_in, n_s_k )
c
      implicit none
c
      type(Tensor2), intent(in) :: stress
      real*8, dimension(2,*), intent(in) :: cm_all
      type(Tensor2), intent(out) :: n_k
      type(Tensor4), optional, intent(in) :: HillT_H_in, HillT_H_s_in
      type(Tensor2), optional, intent(out) :: n_s_k
c
      integer anisotropy_type
      logical Paniso_NonAssoc
c
      if ( present(HillT_H_in)) then
        n_k = get_evolution_dir_n( stress, cm_all, HillT_H_in )
      else
        n_k = get_evolution_dir_n( stress, cm_all )
      endif
c
       anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
       Paniso_NonAssoc = ( anisotropy_type
     &                     == enum_P_aniso_Hill48_NonAssoc )    
c     
      if ( Paniso_NonAssoc
     &     .AND. present(HillT_H_s_in)
     &     .AND. present(n_s_k) ) then
        n_s_k = get_evolution_dir_n( stress, cm_all, HillT_H_s_in )
      elseif ( present(HillT_H_s_in) .XOR. present(n_s_k) ) then
         write(*,*) "get_evolution_dirs_n<<
     & If you provide HillT_H_s_in, you also need to provide n_s_k."
         call cstop ('E R R O R  T E R M I N A T I O N')
      else
         n_s_k = n_k
      endif
c
      return
      end subroutine get_evolution_dirs_n
