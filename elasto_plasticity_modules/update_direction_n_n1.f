c
c
c
      type(Tensor2) function update_direction_n_n1( stress, HillT_H )
c
      use Tensor
      !implicit none
c
      type(Tensor2) :: stress
      type(Tensor4) :: HillT_H
c 
      update_direction_n_n1 = 1. / get_yielding_norm( stress, HillT_H )
     &                        *  (HillT_H .ddot. stress)
c      
      write(*,*) "update_direction_n_n1<< 
     & OoO use get_evolution_dir_n(*) instead"
      stop
c
      end function update_direction_n_n1
