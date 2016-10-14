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
!  routines describing the interaction between photons and the dust
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
    USE observe_mod
    
    USE math_mod
  
    IMPLICIT NONE
  
    PRIVATE

    PUBLIC :: trafo, interact
CONTAINS

    ! ##########################################################################
    ! interaction: steering routine
    ! ---
    SUBROUTINE interact(basics, model, grid, dust, rand_nr, fluxes, photon)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Model_TYP),INTENT(IN)                       :: model
        
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        CALL interact_temp(basics, model, grid, dust, rand_nr, fluxes, photon)

        ! update last point of interaction
        photon%pos_xyz_li(:) = photon%pos_xyz(:)

    END SUBROUTINE interact
    ! ##########################################################################
    ! interaction type 1: pure MC radiative transport
    !                     scattering and absorption at the same time
    ! not used at the moment
    ! ---
    SUBROUTINE interact_mc(basics, grid, dust, rand_nr, photon)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        integer                                          :: i_dust
        !----------------------------------------------------------------------!
        ! ---
        ! 1. select dust species for interaction
        i_dust = dust_select( grid, dust, rand_nr, photon)

        IF (photon%stokes(1) .gt. 0.01_r2) THEN
            ! 2.1 absorption
            photon%stokes(:) = photon%stokes(:) *                              &
                               dust%albedo(i_dust, photon%nr_lam)
            ! 2.2 scattering
            CALL scatter( basics, rand_nr, dust, photon, i_dust )
            ! apply rotation matrix -> get new direction of the photon package

            CALL vecmat(photon)

            ! update the stokes vektor to consider the correct polarization
            ! state (testing phase)
            CALL trafo(dust, photon, i_dust)

            photon%last_interaction_type = 'S'
        ELSE
            photon%inside = .False.
        END IF
    END SUBROUTINE interact_mc

    ! ##########################################################################
    ! interaction type 2: temperature calculation:
    !                     interaction in immediate reemission scheme
    ! ---
    SUBROUTINE interact_temp(basics, model, grid, dust, rand_nr, fluxes, photon)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Model_TYP),INTENT(IN)                       :: model
        
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
            
        integer                                          :: i_dust
        REAL(kind=r2)                                    :: rndx
        
        !----------------------------------------------------------------------!
        ! ---
        ! 1. select dust species for interaction
        i_dust = dust_select( grid, dust, rand_nr, photon)
        
        ! 2. select type of interaction
        CALL GetNewRandomNumber(rand_nr, rndx)
        
        IF ( rndx < dust%albedo(i_dust, photon%nr_lam) ) THEN
            photon%last_interaction_type = 'S'
            ! 2.1 scattering
            ! 2.1.1 peel of the photon if desired
            !       in fact, we create a new photon and calculate 
            !       the energy of this photon on the observers map
            
            IF (basics%do_peel_off) THEN
                CALL peel_off_photon(model, grid, dust, fluxes, photon, i_dust)
            END IF
            CALL scatter( basics, rand_nr, dust, photon, i_dust )
            ! apply rotation matrix -> get new direction of the photon package

            CALL vecmat(photon)

            ! update the stokes vektor to consider the correct polarization
            ! state
            CALL trafo(dust, photon, i_dust)
        ELSE
            ! 2.2 immediate re-emission B&W
            photon%last_interaction_type = 'E'
            CALL immediate(basics, rand_nr, grid, dust, photon, i_dust)
            CALL start_grain(basics, rand_nr, photon)
        END IF

    END SUBROUTINE interact_temp

END MODULE interact_mod
