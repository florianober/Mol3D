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
    use observe_mod
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: run_simu
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE run_simu(basics, fluxes ,grid , model, dust, gas, sources_in)
    
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
    
    !--------------------------------------------------------------------------!
    
    INTEGER                       :: i_map, i_r, i, j, no_pixel, k, l
         
    INTEGER,DIMENSION(0:2*model%n_bin_map,0:2*model%n_bin_map)  :: counter
    INTEGER,DIMENSION(:,:),ALLOCATABLE                          :: notopx
    
    REAL(kind=r2)                 :: hd_stepwidth, dz_min, hd_rmax
         
    REAL(KIND=r2)                 :: pix_res_i, pix_res_j   ! pixelsize [arcsec]
    REAL(KIND=r2)                 :: unit_value             ! unit conversion 
    
    real(kind=r2), dimension(1:2)                :: coor_map
    REAL(kind=r2), dimension(:,:), allocatable   :: calc_px
    REAL(kind=r2), dimension(:,:,:), allocatable :: inten_px
    REAL(kind=r2), dimension(:,:), allocatable   :: continuum_px
    
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
    CALL vis_plane(grid, basics, model, 1, 401)
    CALL vis_plane(grid, basics, model, 2, 401)
    CALL vis_plane(grid, basics, model, 3, 401)
    print *, "| done!                 "
    print *, "|"
    IF (basics%do_raytr) THEN
        ! --- first calculate level populations
        IF (basics%do_velo_ch_map ) THEN
            print *,"| calculate level populations"
            CALL calc_lvlpop(basics, grid , model, gas)
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
                unit_value =  (con_c/gas%trans_freq(gas%tr_cat(1)))**2.0_r2/2.0_r2/con_k 
            ELSE
                PRINT *, '  WARNING: requested unit not found'
                unit_value = 1.0_r2
            END IF
            
            PRINT *,"| | check for pixel"
            k = 1
            l = 0
            DO i = 0, 2*model%n_bin_map
                DO j = 0,2*model%n_bin_map
                    IF (l == int(k*(2*model%n_bin_map)**2*0.01)) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(l/real((2*model%n_bin_map)**2)*100),' % done'//char(27)//'[A'
                        k = k + 1
                    END IF

                    coor_map(1) = model%px_model_length_x(i_map) *             &
                                  (REAL(i, KIND=r2) + 0.5) - model%r_ou /      &
                                  model%zoom_map(1)
                    coor_map(2) = model%px_model_length_y(i_map) *             &
                                 (REAL(j, KIND=r2) + 0.5) - model%r_ou /       &
                                  model%zoom_map(1)
                    
                    CALL get_individual_needed_px(basics, grid, model,         &
                                    i, j,                                      &
                                    model%px_model_length_x(i_map)*0.5_r2,     &
                                    model%px_model_length_x(i_map)*0.5_r2,     &
                                    coor_map(1), coor_map(2), i_map, pixel_list)
                    l = l + 1
                END DO ! coord(2), y
            END DO ! coord(1), x   
                        
            no_pixel = GetSize(pixel_list)
            counter = 0
            ALLOCATE (calc_px (1:no_pixel,1:4))
            ALLOCATE (notopx (1:no_pixel,1:2))


            DO i = 1, no_pixel
                CALL GetLast(pixel_list,pixel_p)
                notopx(i,:) = pixel_p%pixel
                counter(notopx(i,1),notopx(i,2)) = counter(notopx(i,1),notopx(i,2)) +1
                calc_px(i,1:2) = pixel_p%pos_xy
                calc_px(i,3:4) = pixel_p%size_xy
                CALL RemoveLast(pixel_list)
            END DO
            
            PRINT '(A,I7,A,I7,A)', ' | | do raytrace with ',no_pixel,' pixel (', &
                          no_pixel-(2*model%n_bin_map +1)**2, ' subpixel)' 

            
            IF (basics%do_velo_ch_map ) THEN
                PRINT *, '| | calculating velocity channel maps'
                PRINT '(A,F8.2,A)',' | | | wavelength dust    :',dust%lam(dust%cont_map(1)) *1e6, ' micron'
                PRINT '(A,F8.2,A)',' | | | central wavelength :',con_c/gas%trans_freq(gas%tr_cat(1))*1e6, ' micron'
            
                ALLOCATE (inten_px (1:no_pixel,-gas%i_vel_chan:gas%i_vel_chan,1:gas%n_tr))

                inten_px = 0.0_r2
                k = no_pixel/100
                
                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i)
                DO i = 1, no_pixel
                    IF (modulo(i,k) == 0 .or. i == no_pixel) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(i/real(no_pixel)*100.0),' % done'//char(27)//'[A'
                    END IF 
                    
                    inten_px(i,:,:)     =   get_line_intensity_in_px(basics, grid,     &
                                            model, dust, gas,                  &
                                            calc_px(i,1), calc_px(i,2), i_map)
                END DO

                !$omp end do nowait
                !$omp end parallel

                DO i = 1, no_pixel
                    fluxes%channel_map(notopx(i,1),notopx(i,2),:,:) = &
                             fluxes%channel_map(notopx(i,1),notopx(i,2),:,:) + &
                             unit_value * inten_px(i,:,:) *                    &
                             (calc_px(i,3)*calc_px(i,4)*4.) /                  &
                             (model%px_model_length_x(i_map)*model%px_model_length_y(i_map))
                END DO
                
                DEALLOCATE(inten_px)
                print *, '| | | saving channel maps'

                CALL save_ch_map(model, basics, gas, fluxes)
                PRINT *, '| | done!'
                
            END IF
            
            IF (basics%do_continuum_map ) THEN
                ALLOCATE (continuum_px (1:no_pixel,1:dust%n_lam))
                continuum_px = 0.0_r2
                k = no_pixel/100
                PRINT *, '| | calculating continuum maps'
                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i) 
                DO i = 1, no_pixel
                    IF (modulo(i,k) == 0 .or. i == no_pixel) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(i/real(no_pixel)*100.0),' % done'//char(27)//'[A'
                    END IF
                    continuum_px(i,:)   =   get_continuum_intensity_in_px(grid,&
                                            model, dust, gas,                  &
                                            calc_px(i,1), calc_px(i,2), i_map)
    !~                 print *,i
                END DO

                !$omp end do nowait
                !$omp end parallel
                
                DO i = 1, no_pixel
                    fluxes%continuum_map(notopx(i,1),notopx(i,2),:) = &
                             fluxes%continuum_map(notopx(i,1),notopx(i,2),:) + &
                             unit_value * continuum_px(i,:) *                  &
                             (calc_px(i,3)*calc_px(i,4)*4.) /                  &
                             (model%px_model_length_x(i_map)*model%px_model_length_y(i_map))
                END DO
                DEALLOCATE(continuum_px)
                print *, '| | | saving continuum maps'
        
                CALL save_continuum_map(model, basics, dust, fluxes, 2)
                PRINT *, '| | done!'
            END IF
            DEALLOCATE( calc_px, notopx)
        END DO ! orientation map.
        
        WRITE (*,"(A)") " | done! "
    END IF
    
    
!~     print *, 'Simulation finished!'
    
    END SUBROUTINE run_simu
    
END MODULE simulation_mod
