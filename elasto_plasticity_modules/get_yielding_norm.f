c
c
c
      double precision function get_yielding_norm( stress, HillT_H )
c
      use Tensor
      !implicit none
c
      type(Tensor2) :: stress, tmp
      type(Tensor4) :: HillT_H
c
      get_yielding_norm = sqrt (stress .ddot. HillT_H .ddot. stress )
c      
      end function get_yielding_norm
