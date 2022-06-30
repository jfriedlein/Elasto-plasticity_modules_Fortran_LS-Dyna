c
      integer function return_solution_algo ( cm_all )
c
      use cm_manager
      use enumerator_module ! for enumerators
c
      implicit none
c
       real*8, intent(in) :: cm_all(2,*)
       integer hardening_kinematic, anisotropy_type
       integer algo, algo_best
c
       anisotropy_type = int(cm_get_pair('anisotropy______',cm_all))
       hardening_kinematic = INT(cm_get_pair('kin_hard_type___',cm_all))
c
      algo = INT( cm_get_pair('sol_algo________',cm_all) )
      ! Determine the most suitable "best" algorithm
         if ( hardening_kinematic==enum_kinHard_OFF ) then
           if ( anisotropy_type==enum_P_iso ) then ! isotropy
                  algo_best = enum_locIt_ClosestPointP
           elseif  ( floor(anisotropy_type/10.) == enum_P_aniso_Hill )
     &            then ! any kind of Hill48
                  algo_best = enum_locIt_ClosestPointP
           elseif  ( floor(anisotropy_type/10.) == enum_P_aniso_Yld91 )
     &            then ! any kind of Barlat91
                  algo_best = enum_locIt_CuttingPlane
           endif
      ! ... for all kinematic hardening models, we use CP
         else
           algo_best = enum_locIt_CuttingPlane
         endif
c
         if ( algo == enum_locIt_auto ) then ! choose algorithm automatically
            algo = algo_best
         elseif ( algo == enum_locIt_A4approach ) then
            write(*,*) "elasto_plasticity_addi_Px_Hx<< 
     & The A4-approach is currently not available."       
         else ! use the algorithm chosen by the user
            if (algo > algo_best) then
                write(*,*) "elasto_plasticity_addi_Px_Hx<< 
     & The chosen solution algorithm sol_algo is not suitable. Choose 
     & a different algorithm or the automatic (",enum_locIt_auto,
     & ") option."
                call cstop('E R R O R  T E R M I N A T I O N')
            else
                  ! keep the value of algo from cm "sol_algo"
            endif
         endif
c         
      return_solution_algo = algo
c
       end function return_solution_algo
