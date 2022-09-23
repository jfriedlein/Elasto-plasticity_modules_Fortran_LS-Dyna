


      subroutine N5_to_N6( N5, N6)
            implicit none
            real*8, intent(in), dimension(5) :: N5
            real*8, intent(out), dimension(6) :: N6
            integer :: i
          N6(1:6) = 0.
          forall( i=1:5 ) N6(i+1) = N5(i)
          return
      end subroutine N5_to_N6

       subroutine N55_to_N66( N55, N66)
            implicit none
            real*8, intent(in), dimension(5,5) :: N55
            real*8, intent(out), dimension(6,6) :: N66
            integer :: i,j
          N66(1:6,1:6) = 0.
          forall( i=1:5, j=1:5 ) N66(i+1,j+1) = N55(i,j)
          return
      end subroutine N55_to_N66

      subroutine Yld2004(svec, cparams, a, mode, f, grad, hessian)
c**********************************************************************
c      Calculates Barlat's Yld2004 yield function and/or its gradient
c      and/or its Hessian expressed in the natural notation. Algorithm  
c      is based on the paper by Scherzinger, W. M.: A return mapping 
c      algorithm for isotropic and anisotropic plasticity models using
c      a line search method, Computer Methods in Applied Mechanics
c      and Engineering, 2017, DOI: 10.1016/j.cma.2016.11.026
c   [https://gitlab.com/ntnu-physmet/continuum-plasticity]
c
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (5)
c         It is stress tensor given in the natural notation
c (IN)    CPARAMS is REAL*8 array, dimension (:)
c         It is 18 anisotropic parameters of Yld2004 expected in the
c         natural notation, i.e. the first L anisotropic transformation 
c         matrix read
c
c         | 0  L12 L13  0   0   0 |         
c         | 0  L22 L23  0   0   0 |
c         | 0  L32 L33  0   0   0 |
c	    | 0   0   0  L44  0   0 |
c         | 0   0   0   0  L55  0 |
c         | 0   0   0   0   0  L66|
c
c         Coefficients in LPARAMS are ordered as
c         L'12, L'13, L'22, L'23, L'32, L'33, L'44, L'55, L'66, 
c         followed by coeffiients associated with the second anistropic 
c         transformation matrix L'' as
c         L''12, L''13, L''22, L''23, L''32, L''33, L''44, L''55, L''66
c (IN)    A is REAL*8
c         It is the exponent in Yld2004 (a must be >= 2)
c (IN)    MODE is INTEGER
c         Controls which of F, GRAD, HESSIAN is to be calculated
c         = 0  :  only F is calculated
c         = 1  :  only F and GRAD is calculated
c         else :  F, GRAD and HESSIAN are calculated
c (OUT)   F is REAL*8
c         It is scalar value of Yld2004 yield function
c (OUT)   GRAD is REAL*8 array, dimension (5)
c         It is gradient of Yld2004 yield function, computed in 
c         the natural notation
c (OUT)   HESSIAN is REAL*8 array, dimension (5,5)
c         It is the Hessian of Yld2004 yield function computed in 
c         the natural notation
c         
c**********************************************************************
      implicit none
c     ! in/out
      integer     mode
      real*8      svec(5), cparams(18), a, f, grad(5), hessian(5,5)
c     ! internal
      integer     i, j, emod
      real*8      s1vec(6), s2vec(6), scale, tol, asms(9), tmp(6),      
     +            tmpnew(6), tmp1, tmp2, tmp3
c
      real*8      L12, L13, L22, L23, L32, L33, L44, L55, L66,          
     +            M12, M13, M22, M23, M32, M33, M44, M55, M66,          
     +            s11, s12, s13, s21, s22, s23,                         
     +            eigvec11(3), eigvec12(3), eigvec13(3),                
     +            eigvec21(3), eigvec22(3), eigvec23(3),                
     +            s11s21, s11s22, s11s23,                               
     +            s12s21, s12s22, s12s23,                               
     +            s13s21, s13s22, s13s23,                               
     +            as11s21pa2, as11s22pa2, as11s23pa2,                   
     +            as12s21pa2, as12s22pa2, as12s23pa2,                   
     +            as13s21pa2, as13s22pa2, as13s23pa2,                   
     +            dfds11, dfds12, dfds13, dfds21, dfds22, dfds23,       
     +            ddfds11ds11, ddfds12ds12, ddfds13ds13,                
     +            ddfds21ds21, ddfds22ds22, ddfds23ds23,                
     +            ddfds11ds12, ddfds11ds13,                             
     +            ddfds12ds11, ddfds12ds13,                             
     +            ddfds13ds11, ddfds13ds12,                             
     +            ddfds21ds22, ddfds21ds23,                             
     +            ddfds22ds21, ddfds22ds23,                             
     +            ddfds23ds21, ddfds23ds22,                             
     +            ddfds11ds21, ddfds12ds22, ddfds13ds23,                
     +            ddfds11ds22, ddfds11ds23,                             
     +            ddfds12ds21, ddfds12ds23,                             
     +            ddfds13ds21, ddfds13ds22,                             
     +            ddfds21ds11, ddfds21ds12,                             
     +            ddfds21ds13, ddfds22ds11,                             
     +            ddfds22ds12, ddfds22ds13,                             
     +            ddfds23ds11, ddfds23ds12, ddfds23ds13,                
     +            e11xe11v(6), e11xe12v(6), e11xe13v(6),                
     +            e12xe11v(6), e12xe12v(6), e12xe13v(6),                
     +            e13xe11v(6), e13xe12v(6), e13xe13v(6),                
     +            e21xe21v(6), e21xe22v(6), e21xe23v(6),                
     +            e22xe21v(6), e22xe22v(6), e22xe23v(6),                
     +            e23xe21v(6), e23xe22v(6), e23xe23v(6),                
     +            ds11xds11(6,6), ds11xds12(6,6), ds11xds13(6,6),       
     +            ds12xds11(6,6), ds12xds12(6,6), ds12xds13(6,6),       
     +            ds13xds11(6,6), ds13xds12(6,6), ds13xds13(6,6),       
     +            ds21xds21(6,6), ds21xds22(6,6), ds21xds23(6,6),       
     +            ds22xds21(6,6), ds22xds22(6,6), ds22xds23(6,6),       
     +            ds23xds21(6,6), ds23xds22(6,6), ds23xds23(6,6),       
     +            ds11xds21(6,6), ds11xds22(6,6), ds11xds23(6,6),       
     +            ds12xds21(6,6), ds12xds22(6,6), ds12xds23(6,6),       
     +            ds13xds21(6,6), ds13xds22(6,6), ds13xds23(6,6),       
     +            ds21xds11(6,6), ds21xds12(6,6), ds21xds13(6,6),       
     +            ds22xds11(6,6), ds22xds12(6,6), ds22xds13(6,6),       
     +            ds23xds11(6,6), ds23xds12(6,6), ds23xds13(6,6),       
     +            Etmp1212(6,6), Etmp1221(6,6), Etmp2112(6,6),          
     +            Etmp2121(6,6), Etmp2323(6,6), Etmp2332(6,6),          
     +            Etmp3223(6,6), Etmp3232(6,6), Etmp3131(6,6),          
     +            Etmp3113(6,6), Etmp1331(6,6), Etmp1313(6,6),          
     +            E1212p(6,6), E2323p(6,6), E3131p(6,6),                
     +            E1212pp(6,6), E2323pp(6,6), E3131pp(6,6),             
     +            coefE1212p, coefE2323p, coefE3131p,                   
     +            coefE1212pp, coefE2323pp, coefE3131pp,                
     +            H1(6,6), H2(6,6), H3(6,6), H4(6,6), H5(6,6), H6(6,6), 
     +            H7(6,6), H8(6,6), LTxH7xL(5,5), MTxH8xM(5,5), 
     +            dsymLTxH3xM(5,5)             
c
      real*8      ov3, ov4, ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, 
     +            mov4, a2, a1, a1ovf, ova, ov2
      common /consts/ ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, ov2, ov3

      ! @note It is important to define all "common" variables
      ov2     = 0.5d0
	ov3     = 1.d0/3.d0
      ovsqrt3 = 1.d0/sqrt(3.d0)
      ovsqrt2 = 1.d0/sqrt(2.d0)
      ovsqrt6 = ovsqrt3*ovsqrt2
	sqrt2   = 1.d0/ovsqrt2
      sqrt3   = 1.d0/ovsqrt3

      ov4     = 0.25d0
      mov4    = -ov4
      a2      = a-2.d0
      a1      = a-1.d0 
	ova     = 1.d0/a
	tol     = 1.d-5
c      
      L12  = cparams(1)
      L13  = cparams(2)
      L22  = cparams(3)
      L23  = cparams(4)
      L32  = cparams(5)
      L33  = cparams(6)
      L44  = cparams(7)
      L55  = cparams(8)
      L66  = cparams(9)
      M12  = cparams(10)
      M13  = cparams(11)
      M22  = cparams(12)
      M23  = cparams(13)
      M32  = cparams(14)
      M33  = cparams(15)
      M44  = cparams(16)
      M55  = cparams(17)
      M66  = cparams(18)
c
      s1vec(1) = L12*svec(1) + L13*svec(2)
	s1vec(2) = L22*svec(1) + L23*svec(2)
	s1vec(3) = L32*svec(1) + L33*svec(2)
      s1vec(4) = L44*svec(3)
      s1vec(5) = L55*svec(4)
      s1vec(6) = L66*svec(5)
c
      s2vec(1) = M12*svec(1) + M13*svec(2)
	s2vec(2) = M22*svec(1) + M23*svec(2)
	s2vec(3) = M32*svec(1) + M33*svec(2)
      s2vec(4) = M44*svec(3)
      s2vec(5) = M55*svec(4)
      s2vec(6) = M66*svec(5)
c
      if (mode .EQ. 0) then
          emod = 0
      else
          emod = 1
      end if
c     Calculate eigenvalues and eigenvectors of transformed stresses
      call eig(s1vec, emod, s11, s12, s13, eigvec11, eigvec12, eigvec13)
      call eig(s2vec, emod, s21, s22, s23, eigvec21, eigvec22, eigvec23)
c      
c     
      asms(1) = abs(s11-s21)
      asms(2) = abs(s11-s22)
      asms(3) = abs(s11-s23)
c              
      asms(4) = abs(s12-s21)
      asms(5) = abs(s12-s22)
      asms(6) = abs(s12-s23)
c              
      asms(7) = abs(s13-s21)
      asms(8) = abs(s13-s22)
      asms(9) = abs(s13-s23)
c
      scale = asms(1)
      do i=2, 9
          if (asms(i) .GT. scale) scale = asms(i)
      enddo
c     
      f = 0.d0
      if (scale .GT. 1.d-16) then
          do i=1, 9
              f = f + (asms(i)/scale)**a
          enddo
c
c         Compute yield function F
          f = scale*(ov4*f)**ova
      end if
c
      if (mode .EQ. 0) return    
c
	a1ovf   = a1/f
c     rescaling eigenstresses by F
      s11 = s11/f
	s12 = s12/f
	s13 = s13/f
c
      s21 = s21/f
	s22 = s22/f
	s23 = s23/f
c
      s11s21 = s11-s21
      s11s22 = s11-s22
      s11s23 = s11-s23
c      
      s12s21 = s12-s21
      s12s22 = s12-s22
      s12s23 = s12-s23
c
      s13s21 = s13-s21
      s13s22 = s13-s22
      s13s23 = s13-s23
c
      as11s21pa2 = abs(s11s21)**a2
      as11s22pa2 = abs(s11s22)**a2
      as11s23pa2 = abs(s11s23)**a2
c      
      as12s21pa2 = abs(s12s21)**a2
      as12s22pa2 = abs(s12s22)**a2
      as12s23pa2 = abs(s12s23)**a2
c      
      as13s21pa2 = abs(s13s21)**a2
      as13s22pa2 = abs(s13s22)**a2
      as13s23pa2 = abs(s13s23)**a2
cc
      dfds11 = ov4*( s11s21*as11s21pa2 + s11s22*as11s22pa2            
     +        + s11s23*as11s23pa2 )
      dfds12 = ov4*( s12s21*as12s21pa2 + s12s22*as12s22pa2            
     +        + s12s23*as12s23pa2 )                                   
      dfds13 = ov4*( s13s21*as13s21pa2 + s13s22*as13s22pa2            
     +        + s13s23*as13s23pa2 )                                   
cc                                                                    
      dfds21 = ov4*(-s11s21*as11s21pa2 - s12s21*as12s21pa2            
     +        - s13s21*as13s21pa2 )                                   
      dfds22 = ov4*(-s11s22*as11s22pa2 - s12s22*as12s22pa2            
     +        - s13s22*as13s22pa2 )                                   
      dfds23 = ov4*(-s11s23*as11s23pa2 - s12s23*as12s23pa2            
     +        - s13s23*as13s23pa2 )
cc    
      ddfds11ds11 = a1ovf*( ov4*(as11s21pa2 + as11s22pa2 + as11s23pa2)
     +            - dfds11*dfds11 )                                   
      ddfds12ds12 = a1ovf*( ov4*(as12s21pa2 + as12s22pa2 + as12s23pa2)
     +            - dfds12*dfds12 )                                   
      ddfds13ds13 = a1ovf*( ov4*(as13s21pa2 + as13s22pa2 + as13s23pa2)
     +            - dfds13*dfds13 )                                   
      !                                                               
      ddfds21ds21 = a1ovf*( ov4*(as11s21pa2 + as12s21pa2 + as13s21pa2)
     +            - dfds21*dfds21 )                                   
      ddfds22ds22 = a1ovf*( ov4*(as11s22pa2 + as12s22pa2 + as13s22pa2)
     +            - dfds22*dfds22 )                                   
      ddfds23ds23 = a1ovf*( ov4*(as11s23pa2 + as12s23pa2 + as13s23pa2)
     +            - dfds23*dfds23 )
cc
      ddfds11ds12 = -a1ovf*dfds11*dfds12
      ddfds11ds13 = -a1ovf*dfds11*dfds13
      ddfds12ds11 = ddfds11ds12
      ddfds12ds13 = -a1ovf*dfds12*dfds13
      ddfds13ds11 = ddfds11ds13
      ddfds13ds12 = ddfds12ds13
c
      ddfds21ds22 = -a1ovf*dfds21*dfds22
      ddfds21ds23 = -a1ovf*dfds21*dfds23
      ddfds22ds21 = ddfds21ds22
      ddfds22ds23 = -a1ovf*dfds22*dfds23
      ddfds23ds21 = ddfds21ds23
      ddfds23ds22 = ddfds22ds23
cc
      ddfds11ds21 = a1ovf*(mov4*as11s21pa2 - dfds11*dfds21)
      ddfds12ds22 = a1ovf*(mov4*as12s22pa2 - dfds12*dfds22)
      ddfds13ds23 = a1ovf*(mov4*as13s23pa2 - dfds13*dfds23)
      ddfds11ds22 = a1ovf*(mov4*as11s22pa2 - dfds11*dfds22)
      ddfds11ds23 = a1ovf*(mov4*as11s23pa2 - dfds11*dfds23)
      ddfds12ds21 = a1ovf*(mov4*as12s21pa2 - dfds12*dfds21)
      ddfds12ds23 = a1ovf*(mov4*as12s23pa2 - dfds12*dfds23)
      ddfds13ds21 = a1ovf*(mov4*as13s21pa2 - dfds13*dfds21)
      ddfds13ds22 = a1ovf*(mov4*as13s22pa2 - dfds13*dfds22)
c
      ddfds21ds11 = ddfds11ds21
      ddfds21ds12 = ddfds12ds21
      ddfds21ds13 = ddfds13ds21
      ddfds22ds11 = ddfds11ds22
      ddfds22ds12 = ddfds12ds22
      ddfds22ds13 = ddfds13ds22
      ddfds23ds11 = ddfds11ds23
      ddfds23ds12 = ddfds12ds23
      ddfds23ds13 = ddfds13ds23
ccc
c
      tmp(1) = dfds11*eigvec11(1)*eigvec11(1)
     +       + dfds12*eigvec12(1)*eigvec12(1)
     +       + dfds13*eigvec13(1)*eigvec13(1)
      tmp(2) = dfds11*eigvec11(2)*eigvec11(2)
     +       + dfds12*eigvec12(2)*eigvec12(2)
     +       + dfds13*eigvec13(2)*eigvec13(2)
      tmp(3) = dfds11*eigvec11(3)*eigvec11(3)
     +       + dfds12*eigvec12(3)*eigvec12(3)
     +       + dfds13*eigvec13(3)*eigvec13(3)
	tmp(4) = dfds11*eigvec11(2)*eigvec11(3)
     +       + dfds12*eigvec12(2)*eigvec12(3)
     +       + dfds13*eigvec13(2)*eigvec13(3)
	tmp(5) = dfds11*eigvec11(1)*eigvec11(3)
     +       + dfds12*eigvec12(1)*eigvec12(3)
     +       + dfds13*eigvec13(1)*eigvec13(3)
	tmp(6) = dfds11*eigvec11(1)*eigvec11(2)
     +       + dfds12*eigvec12(1)*eigvec12(2)
     +       + dfds13*eigvec13(1)*eigvec13(2)
c
c     tmp into the natural notation
	tmpnew(1) = ovsqrt3*(tmp(1)+tmp(2)+tmp(3))
	tmpnew(2) = ovsqrt6*(2.d0*tmp(3)-tmp(1)-tmp(2))
	tmpnew(3) = ovsqrt2*(tmp(2)-tmp(1))
	tmpnew(4) = sqrt2*tmp(4)
	tmpnew(5) = sqrt2*tmp(5)
	tmpnew(6) = sqrt2*tmp(6)
c
c	Compute GRAD
c
	grad(1) = L12*tmpnew(1) + L22*tmpnew(2) + L32*tmpnew(3)
	grad(2) = L13*tmpnew(1) + L23*tmpnew(2) + L33*tmpnew(3)
      grad(3) = L44*tmpnew(4)
      grad(4) = L55*tmpnew(5)
      grad(5) = L66*tmpnew(6)
c	
      tmp(1) = dfds21*eigvec21(1)*eigvec21(1)
     +       + dfds22*eigvec22(1)*eigvec22(1)
     +       + dfds23*eigvec23(1)*eigvec23(1)
      tmp(2) = dfds21*eigvec21(2)*eigvec21(2)
     +       + dfds22*eigvec22(2)*eigvec22(2)
     +       + dfds23*eigvec23(2)*eigvec23(2)
	tmp(3) = dfds21*eigvec21(3)*eigvec21(3)
     +       + dfds22*eigvec22(3)*eigvec22(3)
     +       + dfds23*eigvec23(3)*eigvec23(3)
	tmp(4) = dfds21*eigvec21(2)*eigvec21(3)
     +       + dfds22*eigvec22(2)*eigvec22(3)
     +       + dfds23*eigvec23(2)*eigvec23(3)
	tmp(5) = dfds21*eigvec21(1)*eigvec21(3)
     +       + dfds22*eigvec22(1)*eigvec22(3)
     +       + dfds23*eigvec23(1)*eigvec23(3)
	tmp(6) = dfds21*eigvec21(1)*eigvec21(2)
     +       + dfds22*eigvec22(1)*eigvec22(2)
     +       + dfds23*eigvec23(1)*eigvec23(2)
c
c     tmp into the natural notation
	tmpnew(1) = ovsqrt3*(tmp(1)+tmp(2)+tmp(3))
	tmpnew(2) = ovsqrt6*(2.d0*tmp(3)-tmp(1)-tmp(2))
	tmpnew(3) = ovsqrt2*(tmp(2)-tmp(1))
	tmpnew(4) = sqrt2*tmp(4)
	tmpnew(5) = sqrt2*tmp(5)
	tmpnew(6) = sqrt2*tmp(6)
c
	grad(1) = grad(1) + M12*tmpnew(1) + M22*tmpnew(2) + M32*tmpnew(3)
	grad(2) = grad(2) + M13*tmpnew(1) + M23*tmpnew(2) + M33*tmpnew(3)
      grad(3) = grad(3) + M44*tmpnew(4)
      grad(4) = grad(4) + M55*tmpnew(5)
      grad(5) = grad(5) + M66*tmpnew(6)
c
c     Stop here if Hessian is not required 
c     and only F and GRAD are returned
      if (mode .EQ. 1) return
c
c
      call outer2vec(eigvec11, eigvec11, e11xe11v)
      call outer2vec(eigvec11, eigvec12, e11xe12v)
      call outer2vec(eigvec11, eigvec13, e11xe13v)
      call outer2vec(eigvec12, eigvec11, e12xe11v)
      call outer2vec(eigvec12, eigvec12, e12xe12v)
      call outer2vec(eigvec12, eigvec13, e12xe13v)
      call outer2vec(eigvec13, eigvec11, e13xe11v)
      call outer2vec(eigvec13, eigvec12, e13xe12v)
      call outer2vec(eigvec13, eigvec13, e13xe13v)
c
      call outer2vec(eigvec21, eigvec21, e21xe21v)
      call outer2vec(eigvec21, eigvec22, e21xe22v)
      call outer2vec(eigvec21, eigvec23, e21xe23v)
      call outer2vec(eigvec22, eigvec21, e22xe21v)
      call outer2vec(eigvec22, eigvec22, e22xe22v)
      call outer2vec(eigvec22, eigvec23, e22xe23v)
      call outer2vec(eigvec23, eigvec21, e23xe21v)
      call outer2vec(eigvec23, eigvec22, e23xe22v)
      call outer2vec(eigvec23, eigvec23, e23xe23v)
c    
c
	call symouter6(e11xe11v, e11xe11v, ds11xds11)
	call symouter6(e11xe11v, e12xe12v, ds11xds12)
	call symouter6(e11xe11v, e13xe13v, ds11xds13)
      call symouter6(e12xe12v, e11xe11v, ds12xds11)
	call symouter6(e12xe12v, e12xe12v, ds12xds12)
	call symouter6(e12xe12v, e13xe13v, ds12xds13)
	call symouter6(e13xe13v, e11xe11v, ds13xds11)
	call symouter6(e13xe13v, e12xe12v, ds13xds12)
	call symouter6(e13xe13v, e13xe13v, ds13xds13)
c             					
      call symouter6(e21xe21v, e21xe21v, ds21xds21)
	call symouter6(e21xe21v, e22xe22v, ds21xds22)
	call symouter6(e21xe21v, e23xe23v, ds21xds23)
      call symouter6(e22xe22v, e21xe21v, ds22xds21)
	call symouter6(e22xe22v, e22xe22v, ds22xds22)
	call symouter6(e22xe22v, e23xe23v, ds22xds23)
	call symouter6(e23xe23v, e21xe21v, ds23xds21)
	call symouter6(e23xe23v, e22xe22v, ds23xds22)
	call symouter6(e23xe23v, e23xe23v, ds23xds23)
c           					 		  		   
	call symouter6(e11xe11v, e21xe21v, ds11xds21)
	call symouter6(e11xe11v, e22xe22v, ds11xds22)
	call symouter6(e11xe11v, e23xe23v, ds11xds23)
      call symouter6(e12xe12v, e21xe21v, ds12xds21)
	call symouter6(e12xe12v, e22xe22v, ds12xds22)
	call symouter6(e12xe12v, e23xe23v, ds12xds23)
	call symouter6(e13xe13v, e21xe21v, ds13xds21)
	call symouter6(e13xe13v, e22xe22v, ds13xds22)
	call symouter6(e13xe13v, e23xe23v, ds13xds23)
c             					 		  		   
	call symouter6(e21xe21v, e11xe11v, ds21xds11)
	call symouter6(e21xe21v, e12xe12v, ds21xds12)
	call symouter6(e21xe21v, e13xe13v, ds21xds13)
      call symouter6(e22xe22v, e11xe11v, ds22xds11)
	call symouter6(e22xe22v, e12xe12v, ds22xds12)
	call symouter6(e22xe22v, e13xe13v, ds22xds13)
	call symouter6(e23xe23v, e11xe11v, ds23xds11)
	call symouter6(e23xe23v, e12xe12v, ds23xds12)
	call symouter6(e23xe23v, e13xe13v, ds23xds13)
c
c     note, only upper triangle of 6x6 sym matrices
c     H1 and H2 is calculated. H3 and H4 are not
c     symmetric, but H3 = H4T, so the upper
c     triangle is enough to calculate
	do i=1, 6
		do j=i, 6
              H1(i,j) = ddfds11ds11*ds11xds11(i,j)
     +				+ ddfds11ds12*ds11xds12(i,j)
     +				+ ddfds11ds13*ds11xds13(i,j)
     +				+ ddfds12ds11*ds12xds11(i,j)
     +				+ ddfds12ds12*ds12xds12(i,j)
     +				+ ddfds12ds13*ds12xds13(i,j)
     +				+ ddfds13ds11*ds13xds11(i,j)
     +				+ ddfds13ds12*ds13xds12(i,j)
     +				+ ddfds13ds13*ds13xds13(i,j)
c                                                 
	 		H2(i,j) = ddfds21ds21*ds21xds21(i,j)
     +				+ ddfds21ds22*ds21xds22(i,j)
     +				+ ddfds21ds23*ds21xds23(i,j)
     +				+ ddfds22ds21*ds22xds21(i,j)
     +				+ ddfds22ds22*ds22xds22(i,j)
     +				+ ddfds22ds23*ds22xds23(i,j)
     +				+ ddfds23ds21*ds23xds21(i,j)
     +				+ ddfds23ds22*ds23xds22(i,j)
     +				+ ddfds23ds23*ds23xds23(i,j)
c                                                 
	 		H3(i,j) = ddfds11ds21*ds11xds21(i,j)
     +				+ ddfds11ds22*ds11xds22(i,j)
     +				+ ddfds11ds23*ds11xds23(i,j)
     +				+ ddfds12ds21*ds12xds21(i,j)
     +				+ ddfds12ds22*ds12xds22(i,j)
     +				+ ddfds12ds23*ds12xds23(i,j)
     +				+ ddfds13ds21*ds13xds21(i,j)
     +				+ ddfds13ds22*ds13xds22(i,j)
     +				+ ddfds13ds23*ds13xds23(i,j)
c                                                 
	 		H4(i,j) = ddfds21ds11*ds21xds11(i,j)
     +				+ ddfds21ds12*ds21xds12(i,j)
     +				+ ddfds21ds13*ds21xds13(i,j)
     +				+ ddfds22ds11*ds22xds11(i,j)
     +				+ ddfds22ds12*ds22xds12(i,j)
     +				+ ddfds22ds13*ds22xds13(i,j)
     +				+ ddfds23ds11*ds23xds11(i,j)
     +				+ ddfds23ds12*ds23xds12(i,j)
     +				+ ddfds23ds13*ds23xds13(i,j)
          end do
      end do
c  
c
c     Eq (32) single prime
	call symouter6(e11xe12v, e11xe12v, Etmp1212)
	call symouter6(e11xe12v, e12xe11v, Etmp1221)
	call symouter6(e12xe11v, e11xe12v, Etmp2112)
	call symouter6(e12xe11v, e12xe11v, Etmp2121)
c
      do i=1, 6
          do j=i, 6
	        E1212p(i,j) = Etmp1212(i,j) + Etmp1221(i,j)
     +                    + Etmp2112(i,j) + Etmp2121(i,j)
          end do
      end do
c
      call symouter6(e12xe13v, e12xe13v, Etmp2323)
	call symouter6(e12xe13v, e13xe12v, Etmp2332)
	call symouter6(e13xe12v, e12xe13v, Etmp3223)
	call symouter6(e13xe12v, e13xe12v, Etmp3232)
c
      do i=1, 6
          do j=i, 6
      	    E2323p(i,j) = Etmp2323(i,j) + Etmp2332(i,j)
     +                    + Etmp3223(i,j) + Etmp3232(i,j)
          end do
      end do
c
	call symouter6(e13xe11v, e13xe11v, Etmp3131)
	call symouter6(e13xe11v, e11xe13v, Etmp3113)
	call symouter6(e11xe13v, e13xe11v, Etmp1331)
	call symouter6(e11xe13v, e11xe13v, Etmp1313)
c
      do i=1, 6
          do j=i, 6
	        E3131p(i,j) = Etmp3131(i,j) + Etmp3113(i,j)
     +                    + Etmp1331(i,j) + Etmp1313(i,j)
          end do
      end do
cc
c     Eq (32) double prime
	call symouter6(e21xe22v, e21xe22v, Etmp1212)
	call symouter6(e21xe22v, e22xe21v, Etmp1221)
	call symouter6(e22xe21v, e21xe22v, Etmp2112)
	call symouter6(e22xe21v, e22xe21v, Etmp2121)
c
      do i=1, 6
          do j=i, 6
	        E1212pp(i,j) = Etmp1212(i,j) + Etmp1221(i,j)
     +                     + Etmp2112(i,j) + Etmp2121(i,j)
          end do
      end do
c
      call symouter6(e22xe23v, e22xe23v, Etmp2323)
	call symouter6(e22xe23v, e23xe22v, Etmp2332)
	call symouter6(e23xe22v, e22xe23v, Etmp3223)
	call symouter6(e23xe22v, e23xe22v, Etmp3232)
c
      do i=1, 6
          do j=i, 6
	        E2323pp(i,j) = Etmp2323(i,j) + Etmp2332(i,j)
     +                     + Etmp3223(i,j) + Etmp3232(i,j)
          end do
      end do
c
	call symouter6(e23xe21v, e23xe21v, Etmp3131)
	call symouter6(e23xe21v, e21xe23v, Etmp3113)
	call symouter6(e21xe23v, e23xe21v, Etmp1331)
	call symouter6(e21xe23v, e21xe23v, Etmp1313)
c
      do i=1, 6
          do j=i, 6
	        E3131pp(i,j) = Etmp3131(i,j) + Etmp3113(i,j)
     +                     + Etmp1331(i,j) + Etmp1313(i,j)
          end do
      end do
c
c     single primed
	if (abs(s11-s12) .LT. tol) then
		coefE1212p = ddfds11ds11 - ddfds11ds12
	else
		coefE1212p = (dfds11 - dfds12)/(s11-s12)
	end if
c
	if (abs(s12-s13) .LT. tol) then
		coefE2323p = ddfds12ds12 - ddfds12ds13
	else
		coefE2323p = (dfds12 - dfds13)/(s12-s13)
	end if
c
	if (abs(s13-s11) .LT. tol) then
		coefE3131p = ddfds13ds13 - ddfds13ds11
	else
		coefE3131p = (dfds13 - dfds11)/(s13-s11)
      end if
c
cc    double primed
      if (abs(s21-s22) .LT. tol) then
		coefE1212pp = ddfds21ds21 - ddfds21ds22
	else
		coefE1212pp = (dfds21 - dfds22)/(s21-s22)
	end if
c
	if (abs(s22-s23) .LT. tol) then
		coefE2323pp = ddfds22ds22 - ddfds22ds23
	else
		coefE2323pp = (dfds22 - dfds23)/(s22-s23)
	end if
c
	if (abs(s23-s21) .LT. tol) then
		coefE3131pp = ddfds23ds23 - ddfds23ds21
	else
		coefE3131pp = (dfds23 - dfds21)/(s23-s21)
      end if
c	
c     note, only upper triangle of 6x6 sym matrices
c     H5 and H6 is calculated
	do i=1, 6
		do j=i, 6
			H5(i,j) = (coefE1212p*E1212p(i,j)
     +                +  coefE2323p*E2323p(i,j)         
     +                +  coefE3131p*E3131p(i,j))*ov2/f
c
			H6(i,j) = (coefE1212pp*E1212pp(i,j)
     +                +  coefE2323pp*E2323pp(i,j)       
     +                +  coefE3131pp*E3131pp(i,j))*ov2/f
c
          end do
          end do
c
c     summing up H1 with H5, and H2 with H6, before tarnsforming them
c     together by L and M, respectively
      do i=1, 6
          do j=i, 6
              H7(i,j) = H1(i,j) + H5(i,j)
              H8(i,j) = H2(i,j) + H6(i,j)
          end do
      end do
c
c
c     computing upper part of sym matrix (LT*H3*M + MT*H4*L) (Eq. 25)
c     note that H3 = H4T, then 
c     dsymLTxH3xM = LT*H3*M + (LT*H3*M)T = 2*sym(LT*H3*M)
      tmp1 = H3(1,1)*M12 + H3(1,2)*M22 + H3(1,3)*M32
      tmp2 = H4(1,2)*M12 + H3(2,2)*M22 + H3(2,3)*M32
      tmp3 = H4(1,3)*M12 + H4(2,3)*M22 + H3(3,3)*M32
      dsymLTxH3xM(1,1) =(L12*tmp1 + L22*tmp2 + L32*tmp3)*2.d0
      dsymLTxH3xM(1,2) = L13*tmp1 + L23*tmp2 + L33*tmp3
     +                 + M13*(H3(1,1)*L12 + H4(1,2)*L22 + H4(1,3)*L32)
     +                 + M23*(H3(1,2)*L12 + H3(2,2)*L22 + H4(2,3)*L32)
     +                 + M33*(H3(1,3)*L12 + H3(2,3)*L22 + H3(3,3)*L32)
      dsymLTxH3xM(1,3) = L44*(H4(1,4)*M12 + H4(2,4)*M22 + H4(3,4)*M32)
     +                 + M44*(H3(1,4)*L12 + H3(2,4)*L22 + H3(3,4)*L32)
      dsymLTxH3xM(1,4) = L55*(H4(1,5)*M12 + H4(2,5)*M22 + H4(3,5)*M32)
     +                 + M55*(H3(1,5)*L12 + H3(2,5)*L22 + H3(3,5)*L32)
      dsymLTxH3xM(1,5) = L66*(H4(1,6)*M12 + H4(2,6)*M22 + H4(3,6)*M32)
     +                 + M66*(H3(1,6)*L12 + H3(2,6)*L22 + H3(3,6)*L32)
      dsymLTxH3xM(2,2) =(L13*(H3(1,1)*M13 + H3(1,2)*M23 + H3(1,3)*M33)
     +                 + L23*(H4(1,2)*M13 + H3(2,2)*M23 + H3(2,3)*M33) 
     +                 + L33*(H4(1,3)*M13 + H4(2,3)*M23 + H3(3,3)*M33))
     +                 * 2.d0
      dsymLTxH3xM(2,3) = L44*(H4(1,4)*M13 + H4(2,4)*M23 + H4(3,4)*M33)
     +                 + M44*(H3(1,4)*L13 + H3(2,4)*L23 + H3(3,4)*L33)
      dsymLTxH3xM(2,4) = L55*(H4(1,5)*M13 + H4(2,5)*M23 + H4(3,5)*M33)
     +                 + M55*(H3(1,5)*L13 + H3(2,5)*L23 + H3(3,5)*L33)
      dsymLTxH3xM(2,5) = L66*(H4(1,6)*M13 + H4(2,6)*M23 + H4(3,6)*M33)
     +                 + M66*(H3(1,6)*L13 + H3(2,6)*L23 + H3(3,6)*L33)
      dsymLTxH3xM(3,3) = 2.d0*H3(4,4)*L44*M44
      dsymLTxH3xM(3,4) = H3(4,5)*L44*M55 + H4(4,5)*L55*M44
      dsymLTxH3xM(3,5) = H3(4,6)*L44*M66 + H4(4,6)*L66*M44
      dsymLTxH3xM(4,4) = 2.d0*H3(5,5)*L55*M55
      dsymLTxH3xM(4,5) = H3(5,6)*L55*M66 + H4(5,6)*L66*M55
      dsymLTxH3xM(5,5) = 2.d0*H3(6,6)*L66*M66
c
c     copmuting upper part of sym matrix LT*H7*L = LT*(H1+H5)*L (Eq. 25)
      LTxH7xL(1,1) = H7(1,1)*L12**2 + H7(2,2)*L22**2 + H7(3,3)*L32**2  
     +    + (H7(2,3)*L22*L32 + H7(1,3)*L12*L32 + H7(1,2)*L12*L22)*2.d0
      LTxH7xL(1,2) = L13*(H7(1,1)*L12 + H7(1,2)*L22 + H7(1,3)*L32)     
     +             + L23*(H7(1,2)*L12 + H7(2,2)*L22 + H7(2,3)*L32)     
     +             + L33*(H7(1,3)*L12 + H7(2,3)*L22 + H7(3,3)*L32)
      LTxH7xL(1,3) = L44*(H7(1,4)*L12 + H7(2,4)*L22 + H7(3,4)*L32)
      LTxH7xL(1,4) = L55*(H7(1,5)*L12 + H7(2,5)*L22 + H7(3,5)*L32)
      LTxH7xL(1,5) = L66*(H7(1,6)*L12 + H7(2,6)*L22 + H7(3,6)*L32)
      LTxH7xL(2,2) = H7(1,1)*L13**2 + H7(2,2)*L23**2 + H7(3,3)*L33**2  
     +    + (H7(2,3)*L23*L33 + H7(1,3)*L13*L33 + H7(1,2)*L13*L23)*2.d0
      LTxH7xL(2,3) = L44*(H7(1,4)*L13 + H7(2,4)*L23 + H7(3,4)*L33)
      LTxH7xL(2,4) = L55*(H7(1,5)*L13 + H7(2,5)*L23 + H7(3,5)*L33)
      LTxH7xL(2,5) = L66*(H7(1,6)*L13 + H7(2,6)*L23 + H7(3,6)*L33)
      LTxH7xL(3,3) = H7(4,4)*L44**2
      LTxH7xL(3,4) = H7(4,5)*L44*L55
      LTxH7xL(3,5) = H7(4,6)*L44*L66
      LTxH7xL(4,4) = H7(5,5)*L55**2
      LTxH7xL(4,5) = H7(5,6)*L55*L66
      LTxH7xL(5,5) = H7(6,6)*L66**2
c
c     copmuting upper part of sym matrix MT*H8*M = MT*(H2+H6)*M (Eq. 25)
      MTxH8xM(1,1) = H8(1,1)*M12**2 + H8(2,2)*M22**2 + H8(3,3)*M32**2  
     +    + (H8(2,3)*M22*M32 + H8(1,3)*M12*M32 + H8(1,2)*M12*M22)*2.d0
      MTxH8xM(1,2) = M13*(H8(1,1)*M12 + H8(1,2)*M22 + H8(1,3)*M32)     
     +             + M23*(H8(1,2)*M12 + H8(2,2)*M22 + H8(2,3)*M32)     
     +             + M33*(H8(1,3)*M12 + H8(2,3)*M22 + H8(3,3)*M32)
      MTxH8xM(1,3) = M44*(H8(1,4)*M12 + H8(2,4)*M22 + H8(3,4)*M32)
      MTxH8xM(1,4) = M55*(H8(1,5)*M12 + H8(2,5)*M22 + H8(3,5)*M32)
      MTxH8xM(1,5) = M66*(H8(1,6)*M12 + H8(2,6)*M22 + H8(3,6)*M32)
      MTxH8xM(2,2) = H8(1,1)*M13**2 + H8(2,2)*M23**2 + H8(3,3)*M33**2  
     +    + (H8(2,3)*M23*M33 + H8(1,3)*M13*M33 + H8(1,2)*M13*M23)*2.d0
      MTxH8xM(2,3) = M44*(H8(1,4)*M13 + H8(2,4)*M23 + H8(3,4)*M33)
      MTxH8xM(2,4) = M55*(H8(1,5)*M13 + H8(2,5)*M23 + H8(3,5)*M33)
      MTxH8xM(2,5) = M66*(H8(1,6)*M13 + H8(2,6)*M23 + H8(3,6)*M33)
      MTxH8xM(3,3) = H8(4,4)*M44**2
      MTxH8xM(3,4) = H8(4,5)*M44*M55
      MTxH8xM(3,5) = H8(4,6)*M44*M66
      MTxH8xM(4,4) = H8(5,5)*M55**2
      MTxH8xM(4,5) = H8(5,6)*M55*M66
      MTxH8xM(5,5) = H8(6,6)*M66**2
c
c     summing up the final hessian - due to the major symmetry of Hessian
c      only the upper triangle is calculated
      do i=1, 5
          do j=i, 5
              hessian(i,j) = LTxH7xL(i,j) + MTxH8xM(i,j) 
     +                     + dsymLTxH3xM(i,j)
          end do
      end do
c
      return
      end subroutine Yld2004

      subroutine outer2vec(e1, e2, v)
c**********************************************************************
c     Computes the outer product of E1 and E2 vector if dimension 3 and 
c     stores it in vector V in the natural notation
c**********************************************************************
	implicit none
	integer     i, j
	real*8      e1(3), e2(3), v(6), tmp(6)
      real*8      ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, ov2, ov3
      common /consts/ ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, ov2, ov3
c	
	tmp(1) = e1(1)*e2(1)
	tmp(2) = e1(2)*e2(2)
	tmp(3) = e1(3)*e2(3)
	tmp(4) = e1(2)*e2(3)
	tmp(5) = e1(1)*e2(3)
	tmp(6) = e1(1)*e2(2)
c	
	v(1) = ovsqrt3*(tmp(1)+tmp(2)+tmp(3))
	v(2) = ovsqrt6*(2.d0*tmp(3)-tmp(1)-tmp(2))
	v(3) = ovsqrt2*(tmp(2)-tmp(1))
	v(4) = sqrt2*tmp(4)
	v(5) = sqrt2*tmp(5)
	v(6) = sqrt2*tmp(6)
c
      return
      end subroutine outer2vec


      subroutine symouter6(e1, e2, M)
c**********************************************************************
c     Computes the upper triangle of the outer product of E1 and E2 
c     vectors of dimension 6 given in the natural notation and stores
c     it in 6x6 matrix M
c**********************************************************************
	implicit none
	integer     i, j
	real*8      e1(6), e2(6), M(6,6)
c		
      do i=1, 6
		do j=i, 6
			M(i,j) = e1(i)*e2(j)
		end do
      end do
c
      return
      end subroutine symouter6

      subroutine eig(svec, mode, e1, e2, e3, evec1, evec2, evec3)
c**********************************************************************
c      Calculates eigenvalues and eigenvectors of a 3x3 symmetric 
c      matrix written in the natural notation as 6x1 vector. Based on 
c      paper by Scherzinger, W. M. and Dohrmann, C. R., A robust  
c      algorithm for finding the eigenvalues and eigenvectors of 3x3 
c      symmetric matrices, Computer Methods in Applied Mechanics and 
c      Engineering, 2008, DOI: 10.1016/j.cma.2008.03.031
c
c**********************************************************************
c (IN)    SVEC is REAL*8 array, dimension (3,3)
c         It is a 6x1 vector representing a symmetric 3x3 tensor in 
c         the natural notation
c (IN)    MODE is integer
c         MODE = 0 returns only eigenvalues, eigenvectors will not be 
c         calculated
c         MODE = else returns both eigenvalues and eigenvectors
c (OUT)   E1, E2, E3 are REAL*8 
c         three eigenvalues
c (OUT)   EVEC1, EVEC2, EVEC3 are REAL*8 array, dimension (3)
c         three eigenvectors
c         
c**********************************************************************
      implicit none  
c     ! in/out
      integer     mode
      real*8      svec(6), e1, e2, e3, evec1(3), evec2(3), evec3(3)
c     ! internal
      integer     i, j, flag, imax, k, m
      real*8      J2, J3, PI, r(3), nr, sdott, Amat(3,3), u1(3), u2(3),
     +            t(3,2), As1(3), As2(3), AR11, AR22, AR23, AR32,   
     +            AR33, s1(3), s2(3), a1, a2, a3, cos3a, ns, w1(3)
      real*8      ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, ov2, ov3
      parameter (PI = 3.14159265359d0)
      common /consts/ ovsqrt3, ovsqrt2, ovsqrt6, sqrt2, sqrt3, ov2, ov3
c
c     build 3x3 deviatoric matrix from svec given in the natural notation
      Amat(1,1) = - ovsqrt6*svec(2) - ovsqrt2*svec(3)
      Amat(2,2) = - ovsqrt6*svec(2) + ovsqrt2*svec(3)
      Amat(3,3) = + 2.d0*ovsqrt6*svec(2)
      Amat(2,3) = ovsqrt2*svec(4)
      Amat(1,3) = ovsqrt2*svec(5)
      Amat(1,2) = ovsqrt2*svec(6)
      Amat(3,2) = Amat(2,3)
      Amat(3,1) = Amat(1,3)
      Amat(2,1) = Amat(1,2)
c
      J2 = 0.d0
      do i=2, 6
          J2 = J2 + svec(i)**2
      end do
      J2 = ov2*J2
c
      J3 = Amat(1,1)*(Amat(2,2)*Amat(3,3)-Amat(3,2)**2)
     +   + Amat(1,3)*(2.d0*Amat(2,1)*Amat(3,2) - Amat(2,2)*Amat(1,3))
     +   - Amat(3,3)*Amat(2,1)**2
c     
      if (J2 .LT. 1.d-30) then
          e1 = 0.d0
          e2 = 0.d0
          e3 = 0.d0
          evec1(1) = 1.d0
          evec1(2) = 0.d0
          evec1(3) = 0.d0
          evec2(1) = 0.d0
          evec2(2) = 1.d0
          evec2(3) = 0.d0
          evec3(1) = 0.d0
          evec3(2) = 0.d0
          evec3(3) = 1.d0
      goto 100
      end if
      cos3a = ov2*J3*(3.d0/J2)**(3.d0/2.d0)
c     to make cos3a within [-1, 1] interval
      cos3a = max(-1.d0, min(1.d0, cos3a))
c
      a1 = ov3*acos(cos3a)
      a3 = a1 + 2.d0*ov3*PI
      a2 = a1 + 4.d0*ov3*PI
c
      if (a1 .LT. PI/6.d0) then
          e1 = 2.d0*sqrt(J2*ov3)*cos(a1)
      else
          e1 = 2.d0*sqrt(J2*ov3)*cos(a3)
      end if
c      
      do i=1, 3
          Amat(i,i) = Amat(i,i) - e1
      end do
c
c     Find the largest column of Amat and store as s1
      ns = 0.d0
      do j=1, 3
          nr = Amat(1,j)**2+Amat(2,j)**2+Amat(3,j)**2
          if (nr .GT. ns) then
              ns = nr
              imax = j
              do i=1, 3
                  s1(i) = Amat(i,j)
              end do
          end if
      end do
c
      do i=1, 3
          s1(i) = s1(i)/sqrt(ns)
      end do
c
      m = 1
      do j=1, 3
          if (j .NE. imax) then
              sdott = s1(1)*Amat(1,j)+s1(2)*Amat(2,j)+s1(3)*Amat(3,j)
              do i=1, 3
                  t(i,m) = Amat(i,j) - sdott*s1(i)
              end do
              m = m+1
          end if
      end do
c
c     Find the largest t column and store as s2
      ns = 0.d0
      do j=1, 2
          nr = t(1,j)**2+t(2,j)**2+t(3,j)**2
          if (nr .GT. ns) then
              ns = nr
              do i=1, 3
                  s2(i) = t(i,j)
              end do
          end if
      end do
c
      do i=1, 3
          s2(i) = s2(i)/sqrt(ns)
      end do
c
c     First eigenvector v1
      evec1(1) = s1(2)*s2(3)-s2(2)*s1(3)
      evec1(2) = s1(3)*s2(1)-s2(3)*s1(1)
      evec1(3) = s1(1)*s2(2)-s2(1)*s1(2)
c
c     Build reduced form of A' matrix (Eq. 22)
      do i=1, 3
          Amat(i,i) = Amat(i,i) + e1
      end do
c
      AR11 = e1
      do i=1, 3
          As1(i) = Amat(i,1)*s1(1) + Amat(i,2)*s1(2) + Amat(i,3)*s1(3)
          As2(i) = Amat(i,1)*s2(1) + Amat(i,2)*s2(2) + Amat(i,3)*s2(3)
      end do
c
      AR22 = s1(1)*As1(1) + s1(2)*As1(2) + s1(3)*As1(3)
      AR23 = s1(1)*As2(1) + s1(2)*As2(2) + s1(3)*As2(3)
      AR32 = s2(1)*As1(1) + s2(2)*As1(2) + s2(3)*As1(3)
      AR33 = s2(1)*As2(1) + s2(2)*As2(2) + s2(3)*As2(3)
c
c     Find the remaining eigenvalues e2, e3 by the Wilkinsons shift
      e2 = ov2*(AR22+AR33) - ov2*sign(1.d0, AR22-AR33)
     +   * sqrt((AR22-AR33)**2 + 4.d0*AR23*AR32)
      e3 = AR22 + AR33 - e2
c
c     returns here if only eigenvalues are required
      if (mode .EQ. 0) goto 100
c
c     Find eigenvectors evec2 and evec3
      do i=1, 3
          Amat(i,i) = Amat(i,i) - e2
      end do
c
      do i=1, 3
          u1(i) = Amat(i,1)*s1(1) + Amat(i,2)*s1(2) + Amat(i,3)*s1(3)
          u2(i) = Amat(i,1)*s2(1) + Amat(i,2)*s2(2) + Amat(i,3)*s2(3)
      end do
c
      nr = u1(1)**2 + u1(2)**2 + u1(3)**2
      ns = u2(1)**2 + u2(2)**2 + u2(3)**2
c     if s1 and s2 are already second and third eigenvectors, then
c     both u1 and u2 and their norms equal zero (ELSE branch)
      if ((nr .GT. 1.d-30) .or. (ns .GT. 1.d-30)) then
          if (nr .GT. ns) then
              do i=1, 3
                  w1(i) = u1(i)/sqrt(nr)
              end do
          else
              do i=1, 3
                  w1(i) = u2(i)/sqrt(ns)
              end do
          end if
          evec2(1) = w1(2)*evec1(3)-evec1(2)*w1(3)
          evec2(2) = w1(3)*evec1(1)-evec1(3)*w1(1)
          evec2(3) = w1(1)*evec1(2)-evec1(1)*w1(2)
c
          evec3(1) = evec1(2)*evec2(3)-evec2(2)*evec1(3)
          evec3(2) = evec1(3)*evec2(1)-evec2(3)*evec1(1)
          evec3(3) = evec1(1)*evec2(2)-evec2(1)*evec1(2)
      else
          do i=1, 3
              evec2(i) = s1(i)
              evec3(i) = s2(i)
          end do
      end if
c
c     adding pressure to get final eigenvalues
100   e1 = e1 + ovsqrt3*svec(1)
      e2 = e2 + ovsqrt3*svec(1)
      e3 = e3 + ovsqrt3*svec(1)
c      
      return
      end subroutine eig
