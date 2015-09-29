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
MODULE simulation_mod

    USE datatype
    USE var_global

    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    USE source_type

    USE grd_mod
    USE math_mod
    USE fileio
    USE transfer_mod
    USE lvlpop_mod
    USE temp_mod
    USE linkedlist_mod
    USE MCRT_mod
    USE raytrace_mod

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: run_simu
    !--------------------------------------------------------------------------!

CONTAINS

    SUBROUTINE run_simu(basics, fluxes, grid, model, dust, gas, sources_in)

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Gas_TYP),INTENT(INOUT)                      :: gas
    TYPE(SOURCES),INTENT(IN)                         :: sources_in

    TYPE(l_list), POINTER                            :: pixel_list => null()
    TYPE(Pixel_TYP), POINTER                         :: pixel_p => null()
    TYPE(Pixel_TYP)                                  :: pixel_data
    TYPE(l_list), POINTER                            :: pixel_list_bak => null()
    !--------------------------------------------------------------------------!

    INTEGER                       :: i_map, i_r, i, j, no_pixel, k, l
         
    
    REAL(kind=r2)                 :: hd_stepwidth, dz_min, hd_rmax
         
    REAL(KIND=r2)                 :: pix_res_i, pix_res_j   ! pixelsize [arcsec]
    REAL(KIND=r2)                 :: unit_value             ! unit conversion 

    REAL(kind=r2), dimension(1:2)                :: coor_map
    REAL(kind=r2), dimension(:,:), allocatable :: inten_px


    !--------------------------------------------------------------------------!
    print *,''                                                    
    print *,"starting simulation"
    print *, "|"
    ! --- calculate or set temperature
    print *,"| temperature calculation"
    CALL set_temperature(basics, grid, model, dust, gas, sources_in,fluxes)
    print *, "| done!                 "
    print *, "|"
    ! --- save all results for later use
    ! first save the model itself
    print *, "| saving the model"
    CALL save_model(grid, basics)
    print *, "| done!                 "
    print *, "|"
    ! now, provide some extra visualisation output 
    ! 1: xz plane, 2: xy plane, 3: yz plane
    print *, "| generate visualisation output"
    CALL vis_plane_fits(grid, basics, model, 1, 256)
    CALL vis_plane_fits(grid, basics, model, 2, 256)
    CALL vis_plane_fits(grid, basics, model, 3, 256)
    print *, "| done!                 "
    print *, "|"

    ! Now we want to create images/sed's, either with an raytracing algorithm
    ! (fast and only option for line) or with pure MC RT
    ! including peel-off method (dust)

    ! MC RT
    IF (basics%do_continuum_mc) THEN
        PRINT *, "| monochromatic Monte-Carlo RT to create images/sed's"
        IF (basics%do_peel_off) PRINT *,"| | Peel-off technique enabled"
        CALL monochromatic_RT(basics, grid, model,                            &
                               dust, sources_in, fluxes)
        print *, "| done!                 "
        PRINT *, "|"
    END IF

    ! raytrace
    IF (basics%do_raytr) THEN
        ! --- first calculate level populations
        IF (basics%do_velo_ch_map ) THEN
            print *,"| calculate level populations"
            CALL calc_lvlpop(grid, model, gas)
            print *, "| done!                 "
            print *, "|"
        END IF
        hd_stepwidth = 0.2_r2

        ! define step width = f(radial grid cell index)
        dz_min = grid%co_mx_a(grid%n(1))
        do i_r=1,grid%n(1)
            if ( (grid%co_mx_a(i_r) - grid%co_mx_a(i_r-1)) < dz_min ) then
              dz_min = grid%co_mx_a(i_r) - grid%co_mx_a(i_r-1)
            end if
        end do
        dz_min = hd_stepwidth * dz_min

        ! maximum radial distance of a pixel to center [in pixel units]
        hd_rmax = real(model%n_bin_map,kind=r2) * sqrt(2.0_r2)

        ! ---
        ! raytracing
        print *, "| raytracing"

        do i_map=1, model%n_map               ! orientation of the map
!~             print *, "    - Map #", i_map, " of ", model%n_map

            IF (GetFluxesName(fluxes) == 'Jy_pix') THEN
                pix_res_i = model%px_model_length_x(i_map) /                   &
                            model%distance * PI /(3600.0_r2*180.0_r2)
                pix_res_j = model%px_model_length_y(i_map) /                   &
                            model%distance * PI /(3600.0_r2*180.0_r2)

                unit_value = 1.0e26_r2*pix_res_i*pix_res_j
            ELSE IF (GetFluxesName(fluxes) == 'T_mb') THEN
                unit_value = (con_c / gas%trans_freq(gas%tr_cat(1)))**2.0_r2 / &
                              2.0_r2/con_k 
            ELSE
                PRINT *, '  WARNING: requested unit not found'
                unit_value = 1.0_r2
            END IF

            PRINT *,"| | check for pixel"
            k = 1
            l = 0
            !$omp parallel num_threads(basics%num_core)
            !$omp do schedule(dynamic) private(i,j, coor_map)
            DO i = 0, 2*model%n_bin_map
!~             DO i = 126, 126
!~                 DO j = 198, 202! 2*model%n_bin_map
                DO j = 0, 2*model%n_bin_map
                    IF (l == int(k*(2.0*model%n_bin_map)**2*0.01)) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(l / real((2.0 *     &
                                 model%n_bin_map)**2)*100.0),                  &
                                 ' % done'//char(27)//'[A'
                        k = k + 1
                    END IF

                    coor_map(1) = model%px_model_length_x(i_map) *             &
                                  (REAL(i, KIND=r2) + 0.5) - model%r_ou /      &
                                  model%zoom_map(i_map)
                    coor_map(2) = model%px_model_length_y(i_map) *             &
                                  (REAL(j, KIND=r2) + 0.5) - model%r_ou /      &
                                  model%zoom_map(i_map)
                    
                    CALL get_individual_needed_px(basics, grid, model,         &
                                    i, j,                                      &
                                    model%px_model_length_x(i_map),            &
                                    model%px_model_length_y(i_map),            &
                                    coor_map(1), coor_map(2), i_map, pixel_list)
                    l = l + 1
                END DO ! coord(2), y
            END DO ! coord(1), x
            !$omp end do nowait
            !$omp end parallel

            no_pixel = GetSize(pixel_list)
            
            PRINT '(A,I7,A,I7,A)', ' | | do raytrace with ',                   &
                        no_pixel, ' (sub)-pixel '


            IF (basics%do_velo_ch_map ) THEN
                PRINT *, '| | calculating velocity channel maps'
                PRINT '(A,F8.2,A)',' | | | wavelength dust    :',              &
                                   dust%lam(dust%cont_map(1)) * 1e6, ' micron'
                PRINT '(A,F8.2,A)',' | | | central wavelength :', con_c /      &
                                   gas%trans_freq(gas%tr_cat(1))*1e6, ' micron'
            
                ALLOCATE (inten_px (-gas%i_vel_chan:gas%i_vel_chan,&
                                    1:gas%n_tr))

                inten_px = 0.0_r2
                k = no_pixel/100
                
                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i, pixel_data) firstprivate(inten_px)
                DO i = 1, no_pixel
                    !$omp critical
                    CALL GetLast(pixel_list, pixel_p)
                    pixel_data = pixel_p
                    CALL RemoveLast(pixel_list)
                    !$omp end critical

                    IF (modulo(i,k) == 0 .or. i == no_pixel) THEN
                        WRITE (*, '(A,I3,A)') ' | | | ',                       &
                                  int(i/real(no_pixel)*100.0),' % done' //     &
                                  char(27)//'[A'
                    END IF 
                    
                    inten_px(:,:) = raytrace_px(grid, model, dust, gas,        &
                                    pixel_data%pos_xy(1), pixel_data%pos_xy(2), i_map)
                
                    !$omp critical
                    fluxes%channel_map(pixel_data%pixel(1),pixel_data%pixel(2),:,:) = &
                        fluxes%channel_map(pixel_data%pixel(1),pixel_data%pixel(2),:,:) + &
                        inten_px(:,:) *                                        &
                        unit_value *                                           &
                        (pixel_data%size_xy(1) * pixel_data%size_xy(2)) /      &
                        (model%px_model_length_x(i_map)*                       &
                        model%px_model_length_y(i_map))
                    CALL AddElement(pixel_list_bak, pixel_data)
                    !$omp end critical
                END DO
                !$omp end do nowait
                !$omp end parallel
                DEALLOCATE(inten_px)
                print *, '| | | saving channel maps'

                CALL save_ch_map(model, basics, gas, fluxes, dust%n_dust)
                PRINT *, '| | done!'

                DO i = 1, no_pixel
                    CALL GetLast(pixel_list_bak, pixel_p)
                    pixel_data = pixel_p
                    CALL RemoveLast(pixel_list_bak)
                    CALL AddElement(pixel_list, pixel_data)
                END DO
            END IF

            IF (basics%do_continuum_raytrace ) THEN
                !reset flux map
                fluxes%continuum_map(:,:,:,:) = 0.0_r2
                
                ALLOCATE (inten_px (1:dust%n_lam,1:4))
                inten_px = 0.0_r2
                k = int(no_pixel/100.0)
                PRINT *, '| | calculating continuum maps'

                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i, pixel_data) firstprivate( inten_px)
                DO i = 1, no_pixel
                    IF (modulo(i,k) == 0 .or. i == no_pixel) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',                        &
                                  int(i/real(no_pixel)*100.0),' % done' //     &
                                  char(27)//'[A'
                    END IF
                    !$omp critical
                    CALL GetLast(pixel_list, pixel_p)
                    pixel_data = pixel_p
                    CALL RemoveLast(pixel_list)
                    !$omp end critical
                    inten_px(:,1) = raytrace_px(                               &
                                        grid, model, dust, gas,                &
                                        pixel_data%size_xy(1), pixel_data%size_xy(2), &
                                        pixel_data%pos_xy(1), pixel_data%pos_xy(2),  &
                                        i_map, 0)
                    !$omp critical
                    fluxes%continuum_map(pixel_data%pixel(1),pixel_data%pixel(2),:, :) =       &
                        fluxes%continuum_map(pixel_data%pixel(1),pixel_data%pixel(2),:,:) + &
                        inten_px(:,:) *                                   &
                        unit_value *                                      &
                        (pixel_data%size_xy(1) * pixel_data%size_xy(2)) /       &
                        (model%px_model_length_x(i_map)*                  &
                        model%px_model_length_y(i_map))
                    !$omp end critical
                END DO
                !$omp end do nowait
                !$omp end parallel
                print *, '| | | saving continuum maps'

                CALL save_continuum_map(model, basics, dust,                   &
                                        fluxes, 3, peel_off=.False.)
                DEALLOCATE(inten_px)
                PRINT *, '| | done!'
            END IF
        END DO ! orientation map.

        WRITE (*,"(A)") " | done! "
    END IF

    END SUBROUTINE run_simu

END MODULE simulation_mod
