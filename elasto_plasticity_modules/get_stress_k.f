c
c
c
      type(Tensor2) function get_stress_k( stress_t, HillT_H, alpha,
     &                                     gamma_k, hardening_type, cm )
c
      use Tensor
      use TensorXkinematics
      use cm_manager
      !implicit none
c
      type(Tensor2) :: Eye, stress_t
      type(Tensor4) :: HillT_H, A, A_inv
      dimension cm(*)
      real gamma_k
      real(kind=8) alpha
      real yield_stress, shearMod_mu
      integer :: hardening_type
      integer :: i,j,k,l
c Material parameters
      yield_stress = cm_get('yieldStress_____',cm)
      shearMod_mu = cm_get('shearMod_mu_____',cm)
c Second order identity tensor
      Eye = identity2(Eye)
c @note
c Here we invert a fourth order tensor (lots of fun)
c
      A= identity4(Eye)
     &                     + (2. * shearMod_mu * gamma_k
     &                        / ( sqrt(2./3.)
     &                            * ( yield_stress - 
     &                                get_hardeningStress_R( alpha, 
     &                                             hardening_type, cm ))
     &                          )) * HillT_H
      get_stress_k = inv ( A ) ** stress_t
      
         ! A_inv = inv(A)
         ! write(*,*) "get_stress_k=", get_stress_k
         !   do i=1,3
         !       do j=1,3
         !           do k=1,3
         !               do l=1,3
         !                   if ( abs(A_inv%abcd(i,j,k,l))>1e-10 ) then
         !     write(*,*) "inv A[",i,j,k,l,"]=",
         !&     A_inv%abcd(i,j,k,l)
         !     endif
         !    enddo
         !    enddo
         !    enddo
         ! enddo
         ! 
         !         do i=1,3
         !       do j=1,3
         !           do k=1,3
         !               do l=1,3
         !                   if ( abs(A%abcd(i,j,k,l))>1e-10 ) then
         !                   write(*,*) "A[",i,j,k,l,"]=", A%abcd(i,j,k,l)
         !                   endif
         !    enddo
         !    enddo
         !    enddo
         ! enddo
         ! 
         !         do i=1,3
         !       do j=1,3
         !     write(*,*) "stress_t[",i,j,"]=",
         !&     stress_t%ab(i,j)
         !    enddo
         ! enddo
         !
         !               do i=1,3
         !       do j=1,3
         !     write(*,*) "stress_k[",i,j,"]=",
         !&     get_stress_k%ab(i,j)
         !    enddo
         !    enddo
      
c      
      end function get_stress_k
