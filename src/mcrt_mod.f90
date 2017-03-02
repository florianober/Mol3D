!------------------------------------------------------------------------------!
! This file is part of Mol3D.
!
!    Mol3D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Mol3D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Mol3D.  If not, see <http://www.gnu.org/licenses/>.
!
!    Diese Datei ist Teil von Mol3D.
!
!    Mol3D ist Freie Software: Sie können es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
!    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Mol3D wird in der Hoffnung, dass es nützlich sein wird, aber
!    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License für weitere Details.
!
!    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
!    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------!
MODULE mcrt_mod

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
    USE tools_mod
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
        USE interact_mod, ONLY : interact
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
        INTEGER(8)                                       :: i_phot
        INTEGER(8)                                       :: kill_photon_count
        !$ INTEGER omp_get_thread_num 
        TYPE(PHOTON_TYP)                                 :: photon
        !----------------------------------------------------------------------!
        IF (.not.PRESENT(i_lam_in)) THEN
            print *,'| | | starting photon transfer ... [this may take a while]'
        END IF

        kill_photon_count = 0
        ! initialize random number generator
        seed = -1
        !$omp parallel num_threads(basics%num_core) FIRSTPRIVATE(rand_nr, seed)
        !$ seed = (omp_get_thread_num()+1)*(-1)

        CALL InitRandgen(rand_nr, seed, basics%randgen_name)
        !----------------------------------------------------------------------!
        !$omp do schedule(dynamic) PRIVATE(i_phot) FIRSTPRIVATE(photon)
        DO i_phot=1,model%no_photon
            ! initiate photon
            IF (PRESENT(i_lam_in)) THEN
                CALL InitPhoton(photon, 1, 'Photon', i_lam_in)
            ELSE
                CALL InitPhoton(photon, 1, 'Photon')
            END IF

            ! show progress
            CALL show_progress(model%no_photon, i_phot, " | | | | - progress : ")

            ! 1. start photon (from included sources)
            CALL start_photon(basics, grid, model, dust, rand_nr, photon,      &
                              sources_in, fluxes)

            IF (photon%kill) THEN
                !$omp atomic
                kill_photon_count = kill_photon_count + 1
                CALL ClosePhoton(photon)  
                CYCLE
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
                    CALL interact(basics, model, grid, dust, rand_nr, fluxes, photon)

                    IF (PRESENT(i_lam_in) .and. (photon%last_interaction_type == 'E') ) THEN
                        ! monochromatic RT
                        ! photon is absorbed, -> next photon
                        photon%inside = .False.
                        EXIT
                    END IF

                    ! move on
                    CALL next_pos_const(model, rand_nr, grid,                  &
                                        dust, photon, deposit_energy)
                    CYCLE
                ELSE
                    IF (.not. photon%inside ) THEN
                        IF ( .not. photon%kill .and. .not. peel_off) THEN
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
            PRINT *, '| | | photon transfer finished                '
        END IF
    END SUBROUTINE MC_photon_transfer
        
    SUBROUTINE monochromatic_RT(basics, grid, model,                           &
                                dust, sources_in, fluxes)
        ! RT at exact wavelengths (only star sources scattering yet)
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
        INTEGER                                          :: i_lam
        !----------------------------------------------------------------------!
        ! reset maps
        fluxes%continuum_map = 0.0_r2

        !TbD: better adjustment of the wavelength table

        DO i_lam = 30, 30!dust%n_lam
            ! First we calculate the scattered light from the source
            WRITE (*,'(A,F7.3, A)') " | | | - wavelength : ",                  &
                                dust%lam(i_lam)*1e6, " micron"
            CALL MC_photon_transfer(basics, grid, model,                       &
                                    dust, sources_in, fluxes,                  &
                                    deposit_energy=.False.,                    &
                                    peel_off=basics%do_peel_off,               &
                                    i_lam_in=i_lam)
            ! Second we calculate the re-emitted scattered light
            ! TbD (use the raytracing algorithm for dust re-emission for now)
            ! THIS IS NOT FULLY IMPLEMENTED YET

            ! We need to start photons from every cell in the disk.....
            !CALL MC_photon_transfer(basics, grid, model,                       &
            !                        dust, sources_in, fluxes,                  &
            !                        deposit_energy=.False.,                    &
            !                        i_lam_in=i_lam)
        END DO

        CALL save_continuum_map(model, basics, dust, fluxes, 2,                &
                                peel_off=basics%do_peel_off)
    END SUBROUTINE monochromatic_RT

END MODULE mcrt_mod
