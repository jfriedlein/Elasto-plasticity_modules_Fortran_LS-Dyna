c
c
c
      real(kind=8) function get_yielding_norm( stress, HillT_H )
c
      use Tensor
      !implicit none
c
      type(Tensor2) :: stress
      type(Tensor4) :: HillT_H
      real(kind=8) tmp
c
      tmp = stress .ddot. HillT_H .ddot. stress
c In case the S:H:S expression gets too negative, we catch this to kick an 
c error message.
      if ( tmp < -1e-20 ) then
          write(*,*) "Yielding norm S:H:S got a bit too negative: ",tmp,
     &               ". Check this case and maybe increase the
     & tolerance for this error."
          call cstop ('E R R O R  T E R M I N A T I O N')
c Tiny negative numbers are okay, but we still need positive numbers for the
c square root, to avoid getting positive yielding norms even for nonsense
c S:H:S values, we limit the value to at least zero (max(...,0))
      else
          get_yielding_norm = sqrt( max(tmp,1e-20) )
      endif     
c      
      end function get_yielding_norm
