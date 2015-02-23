MODULE observe_mod

    USE datatype
    USE var_global

    USE basic_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    USE source_type
    USE photon_type
    
    USE grd_mod
    USE math_mod
    USE transfer_mod
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: peel_off_photon, project_photon_on_map
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE project_photon_on_map(model, fluxes,  photon)
        USE photon_type

        IMPLICIT NONE

        !----------------------------------------------------------------------!
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
        
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_map, i_stokes
        INTEGER                                          :: x, y
        REAL(kind=r2)                                    :: angle
        !----------------------------------------------------------------------!
        DO i_map = 1, model%n_map  !tbd
            ! assume observer is at infinity
            angle = dot_product(photon%dir_xyz, model%dir_xyz(:, i_map))
            IF ( angle .ge. model%al_map(i_map) ) THEN
                CALL get_pixel_on_map(photon%pos_xyz_li,                       &
                                      model%e_x(:, i_map),                     &
                                      model%e_y(:, i_map),                     &
                                      model%px_model_length_x(i_map),          &
                                      model%px_model_length_y(i_map),          &
                                      model%n_bin_map+1, x, y)
                DO i_stokes = 1,4
                    !$omp atomic
                    fluxes%continuum_map(x, y, photon%nr_lam, i_stokes) =      &
                        fluxes%continuum_map(x, y, photon%nr_lam, i_stokes) +  &
                        photon%energy*photon%stokes(i_stokes)
                END DO
            END IF
        END DO

    END SUBROUTINE project_photon_on_map

    PURE SUBROUTINE get_pixel_on_map(p0_xyz, e_x, e_y, px_model_length_x, &
                                px_model_length_y, half_map_size, x, y)
        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        INTEGER,INTENT(OUT)                              :: x, y
        REAl(kind=r2), DIMENSION(1:3), INTENT(IN)        :: p0_xyz, e_x, e_y
        REAl(kind=r2), INTENT(IN)                        :: px_model_length_x
        REAl(kind=r2), INTENT(IN)                        :: px_model_length_y
        INTEGER, INTENT(IN)                              :: half_map_size
        !----------------------------------------------------------------------!
        
        
        x = nint(dot_product(p0_xyz, e_x) / &
                             px_model_length_x + half_map_size-1)
        y = nint(dot_product(p0_xyz, e_y) / &
                             px_model_length_y + half_map_size-1)

    END SUBROUTINE get_pixel_on_map

    SUBROUTINE peel_off_photon(model, grid, dust, fluxes, photon)
        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
        TYPE(PHOTON_TYP)                                 :: photon_peel
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_map
        REAL(kind=r2)                                    :: tau
        !----------------------------------------------------------------------!

        ! create a peel of photon on basis of the current photon and raytrace
        ! it in direction to the observer

        photon_peel = photon
        tau = 0.0_r2
        DO i_map = 1, model%n_map !tbd
            photon_peel%dir_xyz = model%dir_xyz(:, i_map)
            ! transport photon_peel package in direction to the observer

            CALL get_total_tau_direction(model, grid, dust, photon_peel, tau)
            ! When we arrive outside the model space, we know the total optical depth
            ! and can now reduce the amount of energy of this photon_peel package.
            ! Energy rescaling due to the observers distance will be done
            ! after the RT process.

            photon_peel%pos_xyz_li = photon_peel%pos_xyz
            IF (photon%kill) THEN
                photon_peel%energy = 0.0
            ELSE
                photon_peel%energy = photon%energy * exp(-tau)
            END IF
            CALL project_photon_on_map(model, fluxes, photon_peel)
        END DO

    END SUBROUTINE peel_off_photon

END MODULE observe_mod
