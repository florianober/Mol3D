MODULE simulation_mod

    USE datatype
    USE var_globalnew

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
    USE tools_mod
    USE fileio
    USE transfer_mod
    USE string_mod
    USE error_mod
    USE lvlpop_mod
    USE temp_mod
    USE linkedlist_mod
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE :: get_intensity_px
    !--------------------------------------------------------------------------!
    PUBLIC :: run_simu!, primary_temp, primary_scatt, dust_RT
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
    !logical :: in_front, kill_photon
    
    INTEGER                       :: i_map, i_r, i, j, no_pixel, k, l
         
    INTEGER,DIMENSION(0:2*model%n_bin_map,0:2*model%n_bin_map)  :: counter
    INTEGER,DIMENSION(:,:),ALLOCATABLE                          :: notopx
    
    REAL(kind=r2)                 :: hd_stepwidth, dz_min, hd_rmax
         
    REAL(KIND=r2)                 :: rho_size_i, rho_size_j   ! [AU]
    REAL(KIND=r2)                 :: pix_res_i, pix_res_j     ! pixelsize [arcsec]
    REAL(KIND=r2)                 :: unit_value               ! unit conversion 
    
    real(kind=r2), dimension(1:2)                :: coor_map
    real(kind=r2), dimension(1:3)                :: ex, ey
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
    CALL vis_plane(grid, basics,model, 1,401)
    CALL vis_plane(grid, basics,model, 2,401)
    CALL vis_plane(grid, basics,model, 3,401)
    print *, "| done!                 "
    print *, "|"
    IF (basics%do_raytr) THEN
    
        ! --- first calculate level populations
        print *,"| calculate level populations"
        CALL calc_lvlpop(basics, grid , model, gas)
        print *, "| done!                 "
        print *, "|"
        hd_stepwidth = 0.2_r2
        
        ! 1. allocations
       
        ! 2. define step width = f(radial grid cell index)
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
           
            !unity vector pointing to the observer  
                
            grid%dir_xyz(1) = sin(model%th_map(i_map)) * sin(basics%PI2-model%ph_map(i_map))
            grid%dir_xyz(2) = sin(model%th_map(i_map)) * sin(model%ph_map(i_map))
            grid%dir_xyz(3) = sin(basics%PI2-model%th_map(i_map))
            !print '(3(ES15.6E3))',grid%dir_xyz
            
            ! vector marking the +x-direction in the map
            
            ex(1) = -sin(model%ph_map(i_map))
            ex(2) =  sin(basics%PI2-model%ph_map(i_map))
            ex(3) =  0.0_r2
            
            ! vector marking the +y-direction in the map
            
            ey(1) = sin(basics%PI2-model%th_map(i_map)) * (-sin(basics%PI2-model%ph_map(i_map)))
            ey(2) = sin(basics%PI2-model%th_map(i_map)) * (-sin(model%ph_map(i_map)))
            ey(3) = sin(model%th_map(i_map))
            
            
            ! calculate the size of each px  
            !user unit, but should be fixed to [AU]
            rho_size_i = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)   
            rho_size_j = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)
            
            pix_res_i  = rho_size_i/model%distance *PI /(3600.0_r2*180.0_r2)
            pix_res_j  = rho_size_j/model%distance *PI /(3600.0_r2*180.0_r2)
            
            IF (GetFluxesName(fluxes) == 'Jy_pix') THEN
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

                    coor_map(1) = rho_size_i * (REAL(i, KIND=r2) + 0.5) - model%r_ou/model%zoom_map(1)
                    coor_map(2) = rho_size_j * (REAL(j, KIND=r2) + 0.5) - model%r_ou/model%zoom_map(1)
                    
                    CALL get_map_px(basics, grid, model, &
                                  i, j, rho_size_i*0.5_r2, rho_size_j*0.5_r2,     &
                                  coor_map(1), coor_map(2), ex, ey, pixel_list)
                    l = l + 1
                END DO ! coord(2), j
            END DO ! coord(1), i   
                        
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
                k = 1
                
                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i) 
                DO i = 1, no_pixel
                    IF (i == int(k*no_pixel*0.01)) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(i/real(no_pixel)*100),' % done'//char(27)//'[A'
                        k = k + 1
                    END IF
                    
                    inten_px(i,:,:)     =   get_intensity_px(basics, grid, &
                                            model, dust, gas, &
                                            calc_px(i,1), calc_px(i,2), ex, ey)
                END DO

                !$omp end do nowait
                !$omp end parallel

                DO i = 1, no_pixel
                    fluxes%channel_map(notopx(i,1),notopx(i,2),:,:) = &
                             fluxes%channel_map(notopx(i,1),notopx(i,2),:,:) +  unit_value* &
                             inten_px(i,:,:)*(calc_px(i,3)*calc_px(i,4)*4.)/(rho_size_i*rho_size_j)
                END DO
                
                
                
                DEALLOCATE(inten_px)
                print *, '| | | saving channel maps'

                CALL save_ch_map(model, basics, gas, fluxes)
                PRINT *, '| | done!'
                
            END IF
            
            IF (basics%do_continuum_map ) THEN
                ALLOCATE (continuum_px (1:no_pixel,1:dust%n_lam))
                continuum_px = 0.0_r2
                k = 1
                PRINT *, '| | calculating continuum maps'
                !$omp parallel num_threads(basics%num_core)
                !$omp do schedule(dynamic) private(i) 
                DO i = 1, no_pixel
                    IF (i == int(k*no_pixel*0.01)) THEN
                        WRITE (*,'(A,I3,A)') ' | | | ',int(i/real(no_pixel)*100),' % done...'//char(27)//'[A'
                        k = k + 1
                    END IF
                    continuum_px(i,:)   =   get_continuum_px(grid, &
                                            model, dust, gas, &
                                            calc_px(i,1), calc_px(i,2), ex, ey)
    !~                 print *,i
                END DO

                
                !$omp end do nowait
                !$omp end parallel
                
                DO i = 1, no_pixel
                    fluxes%continuum_map(notopx(i,1),notopx(i,2),:) = &
                             fluxes%continuum_map(notopx(i,1),notopx(i,2),:) +  unit_value* &
                             continuum_px(i,:)*(calc_px(i,3)*calc_px(i,4)*4.)/(rho_size_i*rho_size_j)
                END DO
                DEALLOCATE(continuum_px)
                print *, '| | | saving continuum maps'
        
                CALL save_continuum_map(model, basics, dust, fluxes)
                PRINT *, '| | done!'
            END IF
            DEALLOCATE( calc_px, notopx)
        END DO ! orientation map.
        
        WRITE (*,"(A)") " | done! "
    END IF
    
    
!~     print *, 'Simulation finished!'
    
    END SUBROUTINE run_simu
    
 

    FUNCTION get_intensity_px(basics, grid, model, dust, gas, coor_map1, coor_map2, ex, ey)     &
                                        RESULT(px_intensity)
        
    IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Gas_TYP),INTENT(IN)                         :: gas
        
        !--------------------------------------------------------------------------!  
        REAL(KIND=r2), DIMENSION(1:3),INTENT(IN)         :: ex
        REAL(KIND=r2), DIMENSION(1:3),INTENT(IN)         :: ey
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        
        
        REAL(KIND=r2), DIMENSION(-gas%i_vel_chan:gas%i_vel_chan,1:gas%n_tr) :: px_intensity
        
        
        REAL(KIND=r2)                                    :: ray_len
        REAL(KIND=r2)                                    :: dz
        REAL(KIND=r2)                                    :: dz_new
        REAL(KIND=r2)                                    :: dz_sum
        REAL(KIND=r1)                                    :: j_dust
        REAL(KIND=r1)                                    :: alpha_dust
        REAL(KIND=r1)                                    :: velo_dir_xyz
        
        REAL(KIND=r2)                                    :: d_l
        REAL(KIND=r2)                                    :: cell_d_l
        REAL(KIND=r2)                                    :: cell_sum
        REAL(KIND=r1)                                    :: ray_minA
        REAL(KIND=r1)                                    :: abs_err
        REAL(KIND=r1)                                    :: rel_err
        REAL(KIND=r1)                                    :: epsi
        REAL(KIND=r2)                                    :: epsr2
        
        
        REAL(KIND=r1)                                              :: gauss_val
        REAL(KIND=r1), DIMENSION(-gas%i_vel_chan:gas%i_vel_chan)   :: intensity
        REAL(KIND=r1)                                              :: intensity_new
        REAL(KIND=r1)                                              :: intensity_new2
        REAL(KIND=r1)                                              :: j_ges
        REAL(KIND=r1)                                              :: j_ul
        REAL(KIND=r1)                                              :: expo

        REAL(KIND=r1)                                              :: alpha_ges
        REAL(KIND=r1)                                              :: alpha_ul
        
        
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz_new
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz_cell
        
        REAL(KIND=r2), DIMENSION(1:6)                      :: RK_k
        
         
        INTEGER                                          :: vch
        INTEGER                                          :: tr, k
        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
    
        LOGICAL                                          :: kill_photon, log_size
        !--------------------------------------------------------------------------!
        !save photon way, just for debugging
!~         open(unit=1, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'rays.dat', &
!~         action="write", status="unknown", form="formatted")
        

        rel_err = 1.0e-8
        abs_err = 1.0e-20
        
        ! reset intensity in current ray
        j_dust     = 0.0
        j_ges      = 0.0
        j_ul       = 0.0
        alpha_dust = 0.0
        alpha_ges  = 0.0
        alpha_ul   = 0.0
        
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2 ! need (length of ray)**2
        CASE('cylindrical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2
        CASE('cartesian')
            ray_len = 1.
            print *,'TbD, not implemented yet'
            stop
        CASE DEFAULT
            print *, 'selected coordinate system not found, simulation'
            stop
        END SELECT
        
        pos_xyz(:)  = -grid%dir_xyz(:) * sqrt(ray_len) + coor_map1*ex + coor_map2*ey
        px_intensity(:,:) = 0.0  
        intensity(:)      = 0.0  
        ray_minA          = 20.0_r2
        epsr2             = EPSILON(dz)
        dz                = epsr2
        dz_sum            =  0.0_r2
        
        log_size = .False.
        
        IF ( ray_len .gt. 0.0_r2 ) THEN
                       

            nr_cell      = get_cell_nr(grid,pos_xyz)
            DO WHILE ( dz_sum*(1.0_r2+epsr2*1.0e3) .lt. 2.0_r2* sqrt(ray_len) )
                IF ( nr_cell==0 ) THEN

                    CALL path_skip( grid, pos_xyz,grid%dir_xyz, &
                                   pos_xyz_new,nr_cell_new,d_l)
                    dz_sum = dz_sum + d_l
                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new
                    
                END IF

                CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, grid%dir_xyz)

                ! At this point, we have an entrance point(pos_xyz) and exit point (pos_xyz_new) and 
                ! the length of the cell path (d_l) 
                IF ( (nr_cell /= 0 ) .and. (grid%grd_mol_density(nr_cell) .gt. 1.0e-200_r2) ) THEN
                    DO tr = 1, gas%n_tr   ! be careful with more than one transition at once
                        alpha_dust =    grid%grd_dust_density(nr_cell,1) * &
                                        dust%C_ext(1,dust%num_lam_map(dust%cont_map(tr)))
                        
                        j_dust       =  grid%grd_dust_density(nr_cell,1) * &
                                        dust%C_abs(1,dust%num_lam_map(dust%cont_map(tr))) * &
                                        planckhz(grid%t_dust(nr_cell, 1 ),&
                                        con_c/dust%lam(dust%num_lam_map(dust%cont_map(tr))))
!~                             j_dust = 0.0_r2
!~                         alpha_dust = 0.0_r2
                        j_ul =      grid%grd_mol_density(nr_cell)                     *   &
                                    grid%lvl_pop(gas%trans_upper(gas%tr_cat(tr)),nr_cell)  *   &
                                    gas%trans_einstA(gas%tr_cat(tr)) * &
                                    basics%linescale*grid%cell_gauss_a(nr_cell)
!~                             j_ul = 0.0_r2

                        alpha_ul =      grid%grd_mol_density(nr_cell)                     *   &
                                        (grid%lvl_pop(gas%trans_lower(gas%tr_cat(tr)),nr_cell) *   &
                                        gas%trans_einstB_l(gas%tr_cat(tr))                 -   &
                                        grid%lvl_pop(gas%trans_upper(gas%tr_cat(tr)),nr_cell)  *   &
                                        gas%trans_einstB_u(gas%tr_cat(tr))) * &
                                        basics%linescale*grid%cell_gauss_a(nr_cell)
                                        
!~                             alpha_ul = 0.0_r2
                        DO vch = -gas%i_vel_chan, gas%i_vel_chan
                            cell_d_l = d_l
                            dz = cell_d_l * model%ref_unit
                            cell_sum = 0.0_r2
                            pos_xyz_cell = pos_xyz
                            intensity_new   = 0.0
                            intensity_new2  = 0.0
                            DO WHILE (cell_sum .lt. d_l)
                                RK_k(:) = 0.0
                                
                                DO k = 1,6
                                    IF (velo_type == 1 ) THEN
                                    ! use the analytical velocity distribution
                                        velo_dir_xyz   =  dot_product(Set_velo(pos_xyz_cell+ &
                                                      cell_d_l*grid%dir_xyz*RK_c(k),model%kep_const),grid%dir_xyz) 

                                    ELSEIF (velo_type == 2 ) THEN
                                    ! linear interpolation of the velocity 
                                        velo_dir_xyz   = Get_velo(dot_product(grid%velo(nr_cell,:),grid%dir_xyz), &
                                                              dot_product(grid%velo(nr_cell_new,:),grid%dir_xyz), &
                                                              cell_d_l*RK_c(k)/d_l)
                                    END IF
                                    expo = -((gas%velo_channel(vch)-velo_dir_xyz)**2*      &
                                                     grid%cell_gauss_a2(nr_cell))
                                    gauss_val  =  exp(expo)
                                    !gauss_val    = get_expo(expo)
                                    j_ges          = j_ul  * gauss_val  + j_dust
                                    alpha_ges      = alpha_ul * gauss_val  + alpha_dust
                                    !j_ges         = get_opa(j_ul,gauss_val,j_dust)
                                    !alpha_ges     = get_opa(alpha_ul,gauss_val,alpha_dust) 
                                    
                                    RK_k(k) = ( -alpha_ges*(intensity(vch) + &
                                                dz*dot_product(RK_a(:,k),RK_k(:)) ) + j_ges)
                                END DO ! for all k
                        
                                intensity_new  = intensity(vch) + dz*(dot_product(RK_b1(:),RK_k(:)))
                                intensity_new2 = intensity(vch) + dz*(dot_product(RK_b2(:),RK_k(:)))
                                
                                epsi= abs(intensity_new2-intensity_new)/&
                                         ( rel_err*abs(intensity_new) + abs_err)
                                         
    !~                             dz_new = 0.9_r2*dz*exp(-log(epsi)*0.2_r2)
                                dz_new = 0.9*dz*epsi**(-0.2)
                                IF ( epsi .le. 1 ) THEN
                                    intensity(vch) = intensity_new
                                    pos_xyz_cell = pos_xyz_cell+cell_d_l*grid%dir_xyz
                                    cell_sum = cell_sum + cell_d_l
    !~                                 print '(4(ES15.6E3))', norm(pos_xyz), intensity(-1:1)
                                    dz = MIN(dz_new,4*dz)
                                    cell_d_l = dz*model%ref_unitn
                                    IF ( cell_sum + cell_d_l .gt. d_l) THEN
                                        cell_d_l = d_l - cell_sum
                                        dz = cell_d_l*model%ref_unit
                                    END IF
                                
                                ELSE
                                    dz = MAX(dz_new,0.25*dz)
                                    cell_d_l = dz*model%ref_unitn
                                END IF
                            END DO  ! end walk inside one cell
                        END DO !vch
!~                     print *, nr_cell
                    END DO !transitions , please do only one transition per run
                END IF
                
                pos_xyz = pos_xyz_new
!~                 print '(161(ES15.6E3))', norm(pos_xyz), j_ul(-gas%i_vel_chan:gas%i_vel_chan)
                nr_cell = nr_cell_new
                dz_sum = dz_sum + d_l
            END DO !walk in z direction towards observer
            
        px_intensity(:,1) = intensity
        END IF
!~         close(unit=1)
    END FUNCTION get_intensity_px
    
    FUNCTION get_continuum_px(grid, model, dust, gas, coor_map1, coor_map2, ex, ey)     &
                                        RESULT(px_intensity)
        
    IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Gas_TYP),INTENT(IN)                         :: gas
        
        !--------------------------------------------------------------------------!  
        REAL(KIND=r2), DIMENSION(1:3),INTENT(IN)         :: ex
        REAL(KIND=r2), DIMENSION(1:3),INTENT(IN)         :: ey
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        
        
        REAL(KIND=r2), DIMENSION(1:dust%n_lam)           :: px_intensity
        
        
        REAL(KIND=r2)                                    :: ray_len
        REAL(KIND=r2)                                    :: dz
        REAL(KIND=r2)                                    :: dz_new
        REAL(KIND=r2)                                    :: dz_sum
        REAL(KIND=r2)                                    :: j_dust
        REAL(KIND=r2)                                    :: alpha_dust
        
        REAL(KIND=r2)                                    :: d_l
        REAL(KIND=r2)                                    :: cell_d_l
        REAL(KIND=r2)                                    :: cell_sum
        REAL(KIND=r1)                                    :: ray_minA
        REAL(KIND=r1)                                    :: abs_err
        REAL(KIND=r1)                                    :: rel_err
        REAL(KIND=r1)                                    :: epsi
        REAL(KIND=r2)                                    :: epsr2
        
        
        REAL(KIND=r2), DIMENSION(1:dust%n_lam)                     :: intensity
        REAL(KIND=r2)                                              :: intensity_new
        REAL(KIND=r2)                                              :: intensity_new2
        REAL(KIND=r2)                                              :: j_ges

        REAL(KIND=r2)                                              :: alpha_ges

        
        
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz_new
        REAL(KIND=r2), DIMENSION(1:3)                      :: pos_xyz_cell
        
        REAL(KIND=r2), DIMENSION(1:6)                      :: RK_k
        
         
        INTEGER                                          :: i_lam
        INTEGER                                          :: tr, k
        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
    
        LOGICAL                                          :: kill_photon
        !--------------------------------------------------------------------------!

        rel_err = 1.0e-8
        abs_err = 1.0e-20
        
        ! reset intensity in current ray
        j_dust     = 0.0
        j_ges      = 0.0
        alpha_dust = 0.0
        alpha_ges  = 0.0

        
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2 ! need (length of ray)**2
        CASE('cylindrical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2
        CASE('cartesian')
            ray_len = 1.
            print *,'TbD, not implemented yet'
            stop
        CASE DEFAULT
            print *, 'selected coordinate system not found, simulation'
            stop
        END SELECT
        
        pos_xyz(:)  = -grid%dir_xyz(:) * sqrt(ray_len) + coor_map1*ex + coor_map2*ey
        px_intensity(:)   = 0.0  
        intensity(:)      = 0.0  
        ray_minA          = 20.0_r2
        epsr2             = EPSILON(dz)
        dz                = epsr2
        dz_sum            =  0.0_r2
        
        IF ( ray_len .gt. 0.0_r2 ) THEN

            nr_cell      = get_cell_nr(grid,pos_xyz)
            DO WHILE ( dz_sum*(1.0_r2+epsr2*1.0e3) .lt. 2.0_r2* sqrt(ray_len) )
                IF ( nr_cell==0 ) THEN

                    CALL path_skip( grid, pos_xyz,grid%dir_xyz, &
                                   pos_xyz_new,nr_cell_new,d_l)
                    dz_sum = dz_sum + d_l
                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new
                    
                END IF

                CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, grid%dir_xyz)


                ! At this point, we have an entrance point(pos_xyz) and exit point (pos_xyz_new) and 
                ! the length of the cell path (d_l) 
                IF ( (nr_cell /= 0 ) .and. (grid%grd_mol_density(nr_cell) .gt. 1.0e-200_r2) ) THEN
                    DO tr = 1, gas%n_tr   ! be careful with more than one transition at once
                        DO i_lam = 1, dust%n_lam
                            alpha_dust =    grid%grd_dust_density(nr_cell,1) * &
                                        dust%C_ext(1,i_lam)
                            
                            j_dust   =  grid%grd_dust_density(nr_cell,1) * &
                                        dust%C_abs(1,i_lam) * &
                                        planckhz(grid%t_dust(nr_cell, 1 ),&
                                        con_c/dust%lam(i_lam))
                            cell_d_l = d_l
                            dz = cell_d_l * model%ref_unit
                            cell_sum = 0.0_r2
                            pos_xyz_cell = pos_xyz
                            intensity_new   = 0.0_r2
                            intensity_new2  = 0.0_r2
                            DO WHILE (cell_sum .lt. d_l)
                                RK_k(:) = 0.0
                                
                                
                                DO k = 1,6

                                    j_ges          = j_dust
                                    alpha_ges      = alpha_dust
!~                                     IF (intensity(i_lam) <  1.e-200_r2 .and. j_ges < 1e-200_r2) THEN
!~                                         RK_k(k) = 0.0_r2
!~                                     ELSE
                                        RK_k(k) = ( -alpha_ges*(intensity(i_lam) + &
                                                dz*dot_product(RK_a(:,k),RK_k(:)) ) + j_ges)
!~                                     END IF
                                END DO ! for all k
                        
                                intensity_new  = intensity(i_lam) + dz*(dot_product(RK_b1(:),RK_k(:)))
                                intensity_new2 = intensity(i_lam) + dz*(dot_product(RK_b2(:),RK_k(:)))
                                IF (intensity_new < 0.0_r2 .or.intensity_new2 < 0.0_r2 ) THEN
                                    intensity_new  = 0.0_r2
                                    intensity_new2 = 0.0_r2
                                
                                END IF
!~                                 IF (intensity_new < 0.0_r2) THEN
!~                                     print *, 'here', i_lam, j_dust, alpha_dust
!~                                     print *, grid%t_dust(nr_cell, 1 )
!~                                     print *, intensity_new, intensity_new2
!~                                     print *,planckhz(grid%t_dust(nr_cell, 1 ),&
!~                                         con_c/dust%lam(i_lam))
!~                                     print *, intensity(i_lam)
!~                                     print *, intensity(:)
!~                                     print *, RK_k
!~                                     stop
!~                                 END IF
                                epsi= abs(intensity_new2-intensity_new)/&
                                         ( rel_err*abs(intensity_new) + abs_err)
                                         
                                dz_new = 0.9*dz*epsi**(-0.2)
                                IF ( epsi .le. 1 ) THEN
                                    intensity(i_lam) = intensity_new
                                    pos_xyz_cell = pos_xyz_cell+cell_d_l*grid%dir_xyz
                                    cell_sum = cell_sum + cell_d_l
                                    dz = MIN(dz_new,4*dz)
                                    cell_d_l = dz*model%ref_unitn
                                    IF ( cell_sum + cell_d_l .gt. d_l) THEN
                                        cell_d_l = d_l - cell_sum
                                        dz = cell_d_l*model%ref_unit
                                    END IF
                                
                                ELSE

                                    dz = MAX(dz_new,0.25*dz)
                                    cell_d_l = dz*model%ref_unitn
                                END IF

                            END DO  ! end walk inside one cell
                        END DO ! i_lam
!~                     print *, nr_cell
                    END DO !transitions , please do only one transition per run
                END IF
                
                pos_xyz = pos_xyz_new
!~                 print '(161(ES15.6E3))', norm(pos_xyz), j_ul(-gas%i_vel_chan:gas%i_vel_chan)
                nr_cell = nr_cell_new
                dz_sum = dz_sum + d_l
            END DO !walk in z direction towards observer
            
        px_intensity(:) = intensity
        END IF
    END FUNCTION get_continuum_px
    
    RECURSIVE SUBROUTINE get_map_px(basics, grid, model, i,j,  &
                                        xxres,yyres, coor_map1, coor_map2, ex, ey, pixel_list)

        
    IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model

        TYPE(l_list),POINTER,INTENT(INOUT)               :: pixel_list
        TYPE(Pixel_TYP)                                  :: pixel_data
        
        !--------------------------------------------------------------------------!  
        REAL(KIND=r2), DIMENSION(3),INTENT(IN)           :: ex
        REAL(KIND=r2), DIMENSION(3),INTENT(IN)           :: ey
        REAL(KIND=r2),INTENT(IN)                         :: xxres
        REAL(KIND=r2),INTENT(IN)                         :: yyres
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        

        
        REAL(KIND=r2)                                    :: ray_len
        REAL(KIND=r2)                                    :: dz
        REAL(KIND=r2)                                    :: dz_sum
        
        REAL(KIND=r2)                                    :: d_l
        REAL(KIND=r2)                                    :: ray_minA

        REAL(KIND=r2), DIMENSION(3)                      :: pos_xyz
        REAL(KIND=r2), DIMENSION(3)                      :: pos_xyz_new

        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
        INTEGER,INTENT(IN)                               :: i,j

    
        LOGICAL                                          :: kill_photon, log_size
        !--------------------------------------------------------------------------!

        
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2 ! need (length of ray)**2
        CASE('cylindrical')
            ray_len = model%r_ou**2-coor_map1**2-coor_map2**2
        CASE('cartesian')
            ray_len = model%r_ou**2
            print *,'TbD, not implemented yet'
            stop
        CASE DEFAULT
            print *, 'selected coordinate system not found, simulation'
            stop
        END SELECT
        pos_xyz(:)  = -grid%dir_xyz(:) * sqrt(ray_len) + coor_map1*ex + coor_map2*ey

        ray_minA          = 20.0
        dz                = EPSILON(dz)
        dz_sum            =  0.0_r2
        log_size = .False.   
        
        IF ( ray_len .gt. 0.0_r2 ) THEN

            nr_cell      = get_cell_nr(grid,pos_xyz)
!~             print *, 2.0_r2* sqrt(ray_len)
            DO WHILE ( dz_sum*(1.0_r2+epsilon(dz_sum)*1.0e3) .lt. 2.0_r2* sqrt(ray_len) )
                IF ( nr_cell==0 ) THEN

                    CALL path_skip( grid, pos_xyz,grid%dir_xyz, &
                                   pos_xyz_new,nr_cell_new,d_l)
                    dz_sum = dz_sum + d_l
                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new
                    

                END IF
                
                
                
                CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, grid%dir_xyz)
                
                ray_minA = MIN(ray_minA, grid%cell_minA(nr_cell))
                IF ( xxres*yyres .gt. ray_minA ) THEN
                    log_size = .True.
                    EXIT
                END IF

                pos_xyz = pos_xyz_new
                nr_cell = nr_cell_new

                dz_sum = dz_sum + d_l

            END DO !walk in z direction towards observer
            
        END IF
        IF ( log_size) THEN      
            

            CALL        get_map_px(basics, grid, model,                                &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                      &
                        coor_map1+xxres*0.5_r2, coor_map2+yyres*0.5_r2,                    &
                        ex, ey, pixel_list)

            CALL        get_map_px(basics, grid,model,                             &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                      &
                        coor_map1+xxres*0.5_r2, coor_map2-yyres*0.5_r2,                    &
                        ex, ey, pixel_list)
                        
            CALL        get_map_px(basics, grid,model,                             &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                      &
                        coor_map1-xxres*0.5_r2, coor_map2+yyres*0.5_r2,                    &
                        ex, ey,pixel_list )
                        
            CALL        get_map_px(basics, grid, model,                             &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                      &
                        coor_map1-xxres*0.5_r2, coor_map2-yyres*0.5_r2,                    &
                        ex, ey, pixel_list)                                                         
        ELSE
            pixel_data%pixel  =  (/i,j/)
            pixel_data%pos_xy  = (/coor_map1,coor_map2/)
            pixel_data%size_xy = (/xxres,yyres/)
            CALL AddElement(pixel_list,pixel_data)
            CONTINUE
        END IF

    END SUBROUTINE get_map_px
    
END MODULE simulation_mod
