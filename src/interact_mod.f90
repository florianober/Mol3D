! ---
! routines describing the interaction between photons and the dust
! ---
MODULE interact_mod
    USE datatype
    USE var_global
    
    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE photon_type
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

  ! ############################################################################
  ! interaction: steering routine
  ! ---
  SUBROUTINE interact(basics,grid, dust, rand_nr, fluxes, photon)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Fluxes_TYP),INTENT(IN)                      :: fluxes
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!

    CALL interact_temp(basics, grid, dust, rand_nr, fluxes, photon)

    
    ! update last point of interaction
    photon%pos_xyz_li(:) = photon%pos_xyz(:)
                    
  END SUBROUTINE interact

  ! ############################################################################
  ! interaction type 2: temperature calculation:
  !                     interaction in immediate reemission scheme
  ! ---
  subroutine interact_temp(basics, grid, dust, rand_nr, fluxes, photon)

    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Fluxes_TYP),INTENT(IN)                      :: fluxes
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!
        
    integer                                          :: i_dust
    REAL(kind=r2)                                    :: rndx
    
    !--------------------------------------------------------------------------!
    ! ---
    ! 1. select dust species for interaction
    i_dust = dust_select( grid, dust, rand_nr, photon)
    
    ! 2. select type of interaction
    CALL RAN2(rand_nr, rndx)

    IF ( rndx < dust%albedo(i_dust, photon%nr_lam) ) THEN
        photon%last_interaction_type = 'S'
        ! 2.1 scattering
        CALL scatter( basics, rand_nr, dust, photon, i_dust )
        ! apply rotation matrix -> get new direction of the photon package
        CALL vecmat(photon)
        ! update the stokes vektor to consider the correct polarization
        ! state (testing phase)
        CALL trafo(dust, photon, i_dust)
    ELSE
        
        photon%last_interaction_type = 'E'

        ! 2.1 immediate re-emission B&W
        CALL immediate(basics, rand_nr, grid, dust, photon, i_dust)
        CALL start_grain(basics, rand_nr, photon)
    END IF

  end subroutine interact_temp
  

    ! ##########################################################################
    ! select dust species in current cell
    ! ---
    FUNCTION dust_select(grid, dust, rand_nr, photon) RESULT(i_dust_action)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_dust_action
        REAL(kind=r2),DIMENSION(1:dust%n_dust)           :: prob_action
        REAL(kind=r2)                                    :: rndx
        real(kind=r2) :: hd1, hd2
        !----------------------------------------------------------------------!
        
        ! ---
        IF (dust%n_dust==1) THEN
           i_dust_action = 1

        ELSE
           ! approach: probability of interaction of the radiation with a dust 
           !           grain species is in direct proportion to the
           !           a) the number density of these dust grains in the cell
           !      and  b) the extinction cross section of that species
           prob_action(:) = grid%grd_dust_density(photon%nr_cell,:) *          &
                                 dust%C_ext(:,photon%nr_lam)
           
           CALL RAN2(rand_nr, rndx)    
           hd1 = rndx * sum( prob_action(:) )

           i_dust_action = 1
           hd2           = prob_action(i_dust_action)  +  hd1 * epsilon(1.0_r2)
           DO
              IF (hd2 >= hd1 ) THEN
                 EXIT
              ELSE
                 i_dust_action = i_dust_action + 1
                 hd2           = hd2           + prob_action(i_dust_action)
              END IF
           END DO
        END IF
    END FUNCTION dust_select

    
    ! ##########################################################################
    ! Scattering causes a change of stokes vector (I,Q,U,V).  This
    ! transformation occurs in two steps. 
    ! At first, stokes vector is transformed into the new r,l-plane
    ! after the PHI-rotation. 
    ! Therefore values of SIN2PH, COS2PH are necessary. 
    ! The second step contains the transformation of the stokes vector with the 
    ! according  (index is pointed by LOT) scattering matrix (S11,S12,S33,S34).
    ! By the help of  ALBEDO  the absorption is modeled, if this process takes 
    ! place at the same point ("dust grain")  as the scattering.
    ! ---
    SUBROUTINE trafo(dust, photon, i_dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN)                          :: dust
        TYPE(PHOTON_TYP), INTENT(INOUT)                     :: photon        
        INTEGER, INTENT(IN)                                 :: i_dust
        REAL(kind=r2)                                       :: QHELP, i_1
        !----------------------------------------------------------------------!
        ! ---
        ! 1. absorption
        photon%stokes(:) = photon%stokes(:) *                                  &
                           dust%albedo( i_dust, photon%nr_lam )
        i_1       = photon%stokes(1)

        ! 2. scattering
        ! Mathematical positive rotation of r,l-plane around
        !              p-axis by the angle PHI
        QHELP     =  photon%COS2PH * photon%stokes(2) -                        &
                     photon%SIN2PH * photon%stokes(3)
        photon%stokes(3) =  photon%SIN2PH * photon%stokes(2) +                 &
                            photon%COS2PH * photon%stokes(3)
        photon%stokes(2) =  QHELP
        
        ! Mathematical negative rotation around r-axis by THETA transformes
        ! the Stokes vector (I,Q,U,V) with the scattering matrix.
        ! LOT  points to the scattering matrix elements of photon angle which
        ! was determined by throwing dice before.
        photon%stokes(:) = matmul( real(dust%SME( :, :, i_dust, photon%nr_lam, &
                                photon%lot_th ),kind=r2), photon%stokes(:) )

        ! Calculation of an adaption factor NORMFA for the stokes vector so
        ! that intensity after transformation INORM = IBEFORE * ALBEDO.
        photon%stokes(:) = photon%stokes(:) * i_1/photon%stokes(1)

    END SUBROUTINE trafo

end module interact_mod
