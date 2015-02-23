MODULE MCRT_mod

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
    USE observe_mod
    USE transfer_mod
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: MC_photon_transfer, monochromatic_RT
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE MC_photon_transfer(basics, grid,model,                          &
                                  dust, sources_in, fluxes,                    &
                                  deposit_energy, peel_off, i_lam_in)
        ! this is the main routine for MC photon transfer, used to calculate
        ! the temperature distribution and for monochromatic contiuum RT
        USE randgen_type
        USE start_mod
        USE interact_mod
        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(SOURCES),INTENT(IN)                         :: sources_in
        
        TYPE(Randgen_TYP)                                :: rand_nr
        
        !----------------------------------------------------------------------!
        LOGICAL, INTENT(IN)                              :: deposit_energy
        LOGICAL, INTENT(IN)                              :: peel_off
        INTEGER                                          :: seed
        INTEGER, INTENT(IN), OPTIONAL                    :: i_lam_in
        INTEGER(8)                                       :: i_phot, k_phot
        INTEGER(8)                                       :: kill_photon_count
        !$ INTEGER omp_get_thread_num 
        TYPE(PHOTON_TYP)                                 :: photon
        !----------------------------------------------------------------------!
        IF (.not.PRESENT(i_lam_in)) THEN
            print *,'| | | starting photon transfer ... [this may take a while]'
        END IF
        k_phot = model%no_photon/100
        kill_photon_count = 0
        ! initialize random number generator
        seed = -1
        !$omp parallel num_threads(basics%num_core) PRIVATE(rand_nr) FIRSTPRIVATE(seed)
        !$ seed = (omp_get_thread_num()+1)*(-1)

        CALL InitRandgen(rand_nr,seed,'RAN2')
        !----------------------------------------------------------------------!
        !$omp do schedule(dynamic) PRIVATE(i_phot) FIRSTPRIVATE(photon)
        DO i_phot=1,model%no_photon

            ! initiate photon
            IF (PRESENT(i_lam_in)) THEN
                CALL InitPhoton(photon, 1, 'Photon', i_lam_in)
                
            ELSE
                ! show progress
                IF (modulo(i_phot, k_phot) == 0 .or. i_phot==model%no_photon) THEN
                    write (*,'(A,I3,A)') " | | | | - progress : ", &
                    int(i_phot/real(model%no_photon)*100.0), ' % done...'//char(27)//'[A'
                END IF
                CALL InitPhoton(photon, 1, 'Photon')
            END IF
            
            ! 1. start photon (from included sources)
            CALL start_photon(basics, grid, model, rand_nr,                    &
                              dust, photon, sources_in)
            IF (photon%kill) THEN
                !$omp atomic
                kill_photon_count = kill_photon_count + 1
                CALL ClosePhoton(photon)  
                CYCLE
            END IF
            IF (PRESENT(i_lam_in) .or. peel_off) THEN
                CALL peel_off_photon(model, grid, dust, fluxes, photon)
            END IF

            ! 2. determine & go to next point of interaction
            CALL next_pos_const(model, rand_nr, grid,                          &
                                dust, photon, deposit_energy)

            ! 3. transfer through model space
            !    if still inside model space: interaction
            !                                 + go to next position
            !                           else: observe photon leaving
            !                                 the model space
            photon%n_interact = 0
            DO
                IF (photon%inside .and. (photon%n_interact < n_interact_max)) THEN
                    photon%n_interact = photon%n_interact +1

                    CALL interact(basics, grid, dust, rand_nr, fluxes, photon)

                    ! 3. peel of the photon if desired
                    !    in fact, we create a new photon and calculate 
                    !    the energy of this photon on the observers map
                    
                    IF (PRESENT(i_lam_in)) THEN
                        ! monochromatic RT without adding energy to the cells
                        IF (photon%last_interaction_type == 'S') THEN
                            CALL peel_off_photon(model, grid, dust, fluxes, photon)
                        ELSE IF (photon%last_interaction_type == 'E') THEN
                            photon%inside = .False.
                            CYCLE
                        END IF
                    ELSE IF (peel_off) THEN
                        CALL peel_off_photon(model, grid, dust, fluxes, photon)
                    END IF

                    CALL next_pos_const(model, rand_nr, grid,                  &
                                        dust, photon, deposit_energy)
                    CYCLE
                ELSE
                    IF (.not. photon%inside ) THEN
                        IF ( .not. photon%kill .and. deposit_energy .and.      &
                             .not. peel_off) THEN
                            ! observe photon
                            ! here we just count all photons, which are in
                            ! the observers direction.
                            ! This is very cheap, and gives correct but noisy
                            ! results, helpful for debuging
                            ! Thus, for high quality continuum images/seds
                            ! we highly recommend to use the
                            ! 'peel-off technique' (do_continuum_mc = .True.)
                            !
                            CALL project_photon_on_map(model, fluxes, photon)
                        ELSE IF (photon%kill) THEN
                            !$omp atomic
                            kill_photon_count = kill_photon_count + 1
                        END IF
                    END IF
                    EXIT
                END IF
            END DO
            CALL ClosePhoton(photon)
        END DO
        !$omp end do nowait
        !$omp end parallel
        IF (.not.PRESENT(i_lam_in)) THEN
            IF (show_error) PRINT '(A,F5.2,A)', ' | | | | ',                   &
                            REAL(kill_photon_count)/model%no_photon*100_r2,    &
                                   ' % photons killed               '
            print *, '| | | photon transfer finished                '
        END IF
    END SUBROUTINE MC_photon_transfer
        
    SUBROUTINE monochromatic_RT(basics, grid, model,                           &
                                 dust, sources_in, fluxes)
        ! RT at excat wavelengths (only star sources scattering yet)
        USE fileio, ONLY : save_continuum_map
        IMPLICIT NONE

        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(SOURCES),INTENT(IN)                         :: sources_in
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_dust
        INTEGER                                          :: i_tem
        INTEGER                                          :: i_lam
        !----------------------------------------------------------------------!
        ! reset maps
        fluxes%continuum_map = 0.0_r2
        DO i_dust = 1, dust%n_dust
            grid%cell_energy_sum(i_dust,:,1) = grid%cell_energy(i_dust,:) *    &
                                               (basics%PIx4 * grid%cell_vol(:))
        END DO
        !TbD: better adjustment of the wavelength table
        DO i_lam = 12, 55!dust%n_lam
            ! First we calculate the scattered light from the source
            WRITE (*,'(A,F6.2,A)') " | | | | - wavelength : ", &
                    dust%lam(i_lam), ' % done...'//char(27)//'[A'
            CALL MC_photon_transfer(basics, grid, model,                       &
                                    dust, sources_in, fluxes,                  &
                                    deposit_energy=.False.,                    &
                                    peel_off=.False.,                          &
                                    i_lam_in=i_lam)
            ! Second we calculate the re-emitted scattered light
            ! TbD (use the raytracing algorithm for dust re-emission for now)
            !
            
            ! We need to start photons from every cell in the disk.....
            !CALL MC_photon_transfer(basics, grid, model,                       &
            !                        dust, sources_in, fluxes,                  &
            !                        deposit_energy=.False.,                    &
            !                        i_lam_in=i_lam)
        END DO
        grid%cell_energy_sum(:,:,1) = 0.0_r2

        CALL save_continuum_map(model, basics, dust, fluxes, 2, peel_off=.False.)
    END SUBROUTINE monochromatic_RT

END MODULE MCRT_mod