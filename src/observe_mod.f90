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

    USE scatter_mod
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
        REAL(kind=r2)                                    :: SINPHI, COSPHI
        REAL(kind=r2)                                    :: SIN2PH, COS2PH
        REAL(kind=r2)                                    :: QHELP
        REAL(kind=r2),DIMENSION(1:4)                     :: stokes
        REAL(kind=r2),DIMENSION(1:3,1:3)                 :: D_help
        !----------------------------------------------------------------------!
        DO i_map = 1, model%n_map  !tbd
            ! assume observer is at infinity
            !
            angle = dot_product(photon%dir_xyz, model%D_2obs(3, :, i_map))
            IF ( angle .ge. model%al_map(i_map) ) THEN
                !modify stokes vector according to observing direction
                stokes(:) = photon%stokes(:)
                ! we need to find the angle PHI
                !  -> we need to convert from the photon coordinate system 
                !     to the one of the observer
                !     this can be achived by rotate to the global
                !     coordinate system via the photon rotation matrix (photon%D_2global)
                !     and then rotate to the observer with the fixed rotation
                !     matrix of the corresponding observer map (model%D_2obs)
                !     The result is a matrix multiplication and the unknown
                !     angle phi can be read from the resulting matrix
                !
                !                                   / cos(phi)  sin(phi)    0 \
                !  D = D_2obs * photon%D_2global =  |-sin(phi)  cos(phi)    0 |
                !                                   \ 0         0           1 /
                !
                ! Note: This is an approximation and only valid if the
                !       propagation direction equals the direction to the
                !       observer.

                D_help = matmul(model%D_2obs(:, :, i_map), photon%D_2global)
                SINPHI = D_help(1,2)
                COSPHI = D_help(1,1)
                SIN2PH = 2.0_r2 * SINPHI * COSPHI
                COS2PH = 1.0_r2 - 2.0_r2 * SINPHI**2
                IF (D_help(3,3) < 0.92) THEN
                    print *, 'observe'
                    print *, model%D_2obs(:, :, i_map)
                    print *, ''
                    print *, photon%D_2global
                    print *, ''
                    print *, D_help
                    print *, '============================================================'
                END IF
                QHELP     =  COS2PH * stokes(2)  +  SIN2PH * stokes(3)
                stokes(3) = -SIN2PH * stokes(2)  +  COS2PH * stokes(3)
                stokes(2) =  QHELP

                CALL get_pixel_on_map(photon%pos_xyz_li,                       &
                                      model%D_2obs(1, :, i_map),               &
                                      model%D_2obs(2, :, i_map),               &
                                      model%px_model_length_x(i_map),          &
                                      model%px_model_length_y(i_map),          &
                                      model%n_bin_map+1, x, y)
                IF (x .le. model%n_bin_map*2 .and. x .ge. 0) THEN
                    IF (y .le. model%n_bin_map*2 .and. y .ge. 0) THEN
                        DO i_stokes = 1,4
                            !$omp atomic
                            fluxes%continuum_map(x, y, photon%nr_lam, i_stokes) =  &
                              fluxes%continuum_map(x, y, photon%nr_lam, i_stokes)  &
                              + photon%energy*stokes(i_stokes)
                        END DO
                    END IF
                END IF
            END IF
        END DO

    END SUBROUTINE project_photon_on_map

    PURE SUBROUTINE get_pixel_on_map(p0_xyz, e_x, e_y, px_model_length_x,      &
                                px_model_length_y, half_map_size, x, y)
        IMPLICIT NONE

        !----------------------------------------------------------------------!
        INTEGER,INTENT(OUT)                              :: x, y
        REAl(kind=r2), DIMENSION(1:3), INTENT(IN)        :: p0_xyz, e_x, e_y
        REAl(kind=r2), INTENT(IN)                        :: px_model_length_x
        REAl(kind=r2), INTENT(IN)                        :: px_model_length_y
        INTEGER, INTENT(IN)                              :: half_map_size
        !----------------------------------------------------------------------!

        x = nint(dot_product(p0_xyz, e_x) / px_model_length_x + half_map_size-1)
        y = nint(dot_product(p0_xyz, e_y) / px_model_length_y + half_map_size-1)

    END SUBROUTINE get_pixel_on_map

    SUBROUTINE peel_off_photon(model, grid, dust, fluxes, photon, i_dust)
        IMPLICIT NONE

        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
        TYPE(PHOTON_TYP)                                 :: photon_peel

        INTEGER, INTENT(IN)                              :: i_dust
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_map
        REAL(kind=r2)                                    :: tau
        REAL(kind=r2)                                    :: weighting
        REAL(kind=r2)                                    :: PHIPAR
        !----------------------------------------------------------------------!

        ! create a peel of photon on the basis of the current photon and raytrace
        ! raytrace it in direction to the observer

        photon_peel = photon
        tau = 0.0_r2
        DO i_map = 1, model%n_map !tbd
            ! direction of the new photon is the direction to the observer
            photon_peel%dir_xyz = model%D_2obs(3, :, i_map)
            ! transport photon_peel package in direction to the observer
            CALL get_total_tau_direction(model, grid, dust, photon_peel, tau)
            ! When we arrive outside the model space, we know the total optical depth
            ! and can now reduce the amount of energy of this photon_peel package.
            ! Energy rescaling due to the observers distance will be done
            ! after the RT process.
            !
            ! scattering in the observer's direction
            CALL vecmat(photon_peel, dust%D_ANG)
            !
            
            IF (photon%last_interaction_type == 'S') THEN !simulate a scatering event
                ! get the phi and theta angles for the given direction
                IF (dust%aniiso == 1 ) THEN
                    PHIPAR = (sqrt(photon_peel%stokes(2)**2 + photon_peel%stokes(3)**2) /     &
                              photon_peel%stokes(1)) * &
                              real(-dust%SME(1, 2, i_dust, photon_peel%nr_lam,photon_peel%lot_th) /  &
                              dust%SME(1, 1, i_dust, photon_peel%nr_lam, photon_peel%lot_th), kind=r2)

                    weighting =     1.0_r2/(2.0_r2 * PI) *                     &
                                    (1.0_r2 - PHIPAR * photon_peel%COS2PH) *   &
                                    dust%phase_pdf(photon_peel%lot_th, i_dust, &
                                    photon_peel%nr_lam)
                ELSEIF (dust%aniiso == 3 ) THEN
                    weighting = dust%phase_pdf(photon_peel%lot_th, i_dust, &
                                               photon_peel%nr_lam)*0.5_r2/PI
                ELSE
                    weighting = 0.25_r2/PI
                END IF
                
                ! Now, we need to convert the stokes vector
                CALL trafo(dust, photon_peel, i_dust)
                
            ELSE !peel off from a non scattering event (e.g. defind source)
                weighting = 0.25_r2/PI
            END IF
            
            photon_peel%energy = photon%energy * exp(-tau) * weighting
            
            photon_peel%pos_xyz_li = photon_peel%pos_xyz
            IF (photon%kill) THEN
                photon_peel%energy = 0.0
            END IF

            CALL project_photon_on_map(model, fluxes, photon_peel)
        END DO

    END SUBROUTINE peel_off_photon

END MODULE observe_mod
