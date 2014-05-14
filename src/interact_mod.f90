! ---
! routines describing the interaction between photons and the dust
! ---
MODULE interact_mod
    USE datatype
    USE var_globalnew
    
    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE simu_type
    USE basic_type
    USE fluxes_type
    
    USE start_mod
    USE immediate_mod
    USE scatter_mod
    
    USE math_mod
  
    IMPLICIT NONE
  
    PRIVATE
  
    PUBLIC  :: interact
  
contains

  ! ################################################################################################
  ! interaction: steering routine
  ! ---
  SUBROUTINE interact(basics,grid, dust, rand_nr, fluxes, simu_var, t_dust, i_star_abs,dt_dust)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    !TYPE(Model_TYP),INTENT(IN)                        :: model
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Fluxes_TYP),INTENT(IN)                       :: fluxes
    
    TYPE(Simu_TYP),INTENT(INOUT)                     :: simu_var
    !--------------------------------------------------------------------------!

    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: t_dust
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: dt_dust
    REAL(kind=r2),DIMENSION(1:dust%n_dust, 1:dust%n_lam, 1:grid%n_cell),INTENT(INOUT)    :: i_star_abs
    !--------------------------------------------------------------------------!
    ! ---
!~     select case(simu_type)
!~ 
!~     case (1)
       CALL interact_temp(basics,grid,dust,rand_nr,fluxes,simu_var, t_dust, i_star_abs,dt_dust)
!~ 
!~     case (2,3,4,5)
!~        call interact_prim()
!~ 
!~     case default
!~        print *,"<!> subroutine interact: wrong input parameter. stopped."
!~        call stop_mc3d()
!~     end select
    
    ! update last point of interaction
    simu_var%pos_xyz_li(:) = simu_var%pos_xyz(:)
                    
  END SUBROUTINE interact
  

  ! ################################################################################################
  ! interaction type 1: pure MC radiative transfer:
  !                     scattering and absorption occur simultaneously
  ! ---
!~   subroutine interact_prim()
!~     use var_global
!~     use math_mod
!~     use start_mod
!~     use scatter_mod
!~     
!~     implicit none
!~     integer :: i_dust_action
!~     ! ---
!~     ! 1. select dust species for interaction
!~     call dust_select( i_dust_action )      
!~     
!~     ! 2. scattering & absorption
!~     call scatter( i_dust_action )
!~     call vecmat()
!~     call trafo( i_dust_action )
!~     
!~   end subroutine interact_prim
  

  ! ################################################################################################
  ! interaction type 2: temperature calculation:
  !                     interaction in immediate reemission scheme
  ! ---
  subroutine interact_temp(basics,grid,dust,rand_nr,fluxes,simu_var,t_dust, i_star_abs,dt_dust)

    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
!~     TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Fluxes_TYP),INTENT(IN)                      :: fluxes
    
    TYPE(Simu_TYP),INTENT(INOUT)                     :: simu_var
    !--------------------------------------------------------------------------!
        
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: t_dust
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: dt_dust
    REAL(kind=r2),DIMENSION(1:dust%n_dust, 1:dust%n_lam, 1:grid%n_cell),INTENT(INOUT)    :: i_star_abs
    !--------------------------------------------------------------------------!
    integer                                           :: i_dust_action 
    REAL(kind=r2)                                     :: rndx
    
    !--------------------------------------------------------------------------!
    ! ---
    ! 1. select dust species for interaction
    CALL dust_select( grid,dust,rand_nr,simu_var,i_dust_action)      
    
    ! 2. select type of interaction
    CALL RAN2(rand_nr, rndx)
    !print *, 'rand',rndx
    IF ( rndx < (dust%C_sca(i_dust_action,simu_var%nr_lam) / dust%C_ext(i_dust_action,simu_var%nr_lam)) ) then               

       ! 2.1 scattering
        CALL scatter( basics, rand_nr, dust, simu_var, i_dust_action )
        ! --- ---
        ! tbd: call trafo( i_dust_action )
        !      to be implemented to correctly consider the polarization state
        !      of the photon package during further radiative transfer
        ! --- ---
        CALL vecmat(simu_var)
    ELSE                
        CALL immediate( basics, rand_nr, grid,dust, simu_var, i_dust_action,t_dust, i_star_abs,dt_dust )
        CALL start_grain(basics, rand_nr, fluxes, simu_var)
    END IF
  end subroutine interact_temp
  

  ! ################################################################################################
  ! select dust species in current cell
  ! ---
  SUBROUTINE dust_select(grid,dust,rand_nr,simu_var,i_dust_action)

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(Simu_TYP),INTENT(INOUT)                    :: simu_var
    !--------------------------------------------------------------------------!
    integer,INTENT(OUT)                             :: i_dust_action 
    REAL(kind=r2)                                    :: rndx
    real(kind=r2) :: hd1, hd2
    !--------------------------------------------------------------------------!
    
    ! ---
    IF (dust%n_dust==1) THEN
       i_dust_action = 1

    ELSE
       ! approach: probability of interaction of the radiation with a dust grain species
       !           is in direct proportion to the
       !           a) the number density of these dust grains in the cell
       !      and  b) the extinction cross section of that species
       simu_var%prob_action(:) = grid%grd_dust_density(simu_var%nr_cell,:) * dust%C_ext(:,simu_var%nr_lam)
       
       CALL RAN2(rand_nr,rndx)    
       hd1 = rndx * sum( simu_var%prob_action(:) )

       i_dust_action = 1
       hd2           = simu_var%prob_action(i_dust_action)  +  hd1 * epsilon(1.0_r2)
       DO
          IF (hd2 >= hd1 ) THEN
             EXIT
          ELSE
             i_dust_action = i_dust_action + 1
             hd2           = hd2           + simu_var%prob_action(i_dust_action)
          END IF
       END DO
    END IF

  END SUBROUTINE dust_select

end module interact_mod
