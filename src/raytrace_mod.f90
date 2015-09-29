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
! Module containing all subroutines for the raytracing (of lines and continuum)
!        this module needs to be cleaned and reviewed!
! ---
MODULE raytrace_mod
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
    USE model_mod
    USE transfer_mod

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: raytrace_px, get_individual_needed_px
    !--------------------------------------------------------------------------!
    INTERFACE raytrace_px
        MODULE PROCEDURE get_line_intensity_in_px,                             &
                         get_continuum_intensity_in_px
    END INTERFACE

CONTAINS
    FUNCTION get_line_intensity_in_px(grid, model, dust, gas,                  &
                              coor_map1, coor_map2, i_map)                     &
                              RESULT(px_intensity)

    IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Gas_TYP),INTENT(IN)                         :: gas
        !----------------------------------------------------------------------!
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        INTEGER,INTENT(IN)                               :: i_map

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

        REAL(KIND=r2), DIMENSION(1:3)                   :: pos_xyz
        REAL(KIND=r2), DIMENSION(1:3)                   :: pos_xyz_new
        REAL(KIND=r2), DIMENSION(1:3)                   :: pos_xyz_cell

        REAL(KIND=r2), DIMENSION(1:6)                   :: RK_k

        INTEGER                                         :: vch
        INTEGER                                         :: tr, k
        INTEGER                                         :: nr_cell
        INTEGER                                         :: nr_cell_new

        LOGICAL                                         :: kill_photon, log_size
        !----------------------------------------------------------------------!

        ! reset intensity in current ray
        j_dust     = 0.0
        j_ges      = 0.0
        j_ul       = 0.0
        alpha_dust = 0.0
        alpha_ges  = 0.0
        alpha_ul   = 0.0
        velo_dir_xyz = 0.0
        
        ray_len = get_ray_length(GetGridName(grid), model%r_ou,                &
                                 coor_map1, coor_map2)

        ! model%D_2obs(1, :, i_map) gives the x-axis in the observes plane
        ! model%D_2obs(2, :, i_map) gives the y-axis in the observes plane
        ! model%D_2obs(3, :, i_map) gives the direction to the observer (z-axis)

        pos_xyz(:) = -model%D_2obs(3, :, i_map) * sqrt(ray_len) +              &
                      coor_map1*model%D_2obs(1, :, i_map) +                    &
                      coor_map2*model%D_2obs(2, :, i_map)
        px_intensity(:,:) = 0.0  
        intensity(:)      = 0.0  
        ray_minA          = 20.0_r2
        epsr2             = EPSILON(dz)
        dz                = epsr2
        dz_sum            =  0.0_r2

        log_size = .False.

        IF ( ray_len .gt. 0.0_r2 ) THEN
            nr_cell = get_cell_nr(grid,pos_xyz)
            DO WHILE ( dz_sum*(1.0_r2+epsr2*1.0e3) .lt. 2.0_r2* sqrt(ray_len) )
                IF ( nr_cell==0 ) THEN
                    CALL path_skip( grid, pos_xyz,model%D_2obs(3, :, i_map), &
                                    pos_xyz_new, nr_cell_new, d_l)
                    dz_sum = dz_sum + d_l
                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new
                    IF (.not. check_inside(nr_cell, grid) ) EXIT
                END IF

                CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new,   &
                           d_l, kill_photon, model%D_2obs(3, :, i_map))

                ! At this point, we have an entrance point(pos_xyz) and
                ! exit point (pos_xyz_new) and 
                ! the length of the cell path (d_l) 
                IF ( (nr_cell /= 0 ) .and. (grid%grd_mol_density(nr_cell) .gt. 1.0e-200_r2) ) THEN
                    DO tr = 1, gas%n_tr   ! more than one transition at once is
                                          ! not yet fully implemented
                        alpha_dust =    sum(grid%grd_dust_density(nr_cell,:) * &
                                        dust%C_ext(:,dust%num_lam_map(dust%cont_map(tr))))

                        j_dust       =  sum(grid%grd_dust_density(nr_cell,:) * &
                                        dust%C_abs(:,dust%num_lam_map(dust%cont_map(tr))) * &
                                        planckhz(grid%t_dust(nr_cell, :), &
                                        dust%nu(dust%num_lam_map(dust%cont_map(tr)))))

                        j_ul =      grid%grd_mol_density(nr_cell)                     *   &
                                    grid%lvl_pop(gas%trans_upper(gas%tr_cat(tr)),nr_cell)  *   &
                                    gas%trans_einstA(gas%tr_cat(tr)) * &
                                    linescale*grid%cell_gauss_a(nr_cell)

                        alpha_ul =  grid%grd_mol_density(nr_cell) *            &
                                    (grid%lvl_pop(gas%trans_lower(gas%tr_cat(tr)),nr_cell) *   &
                                    gas%trans_einstB_l(gas%tr_cat(tr))                 -   &
                                    grid%lvl_pop(gas%trans_upper(gas%tr_cat(tr)),nr_cell)  *   &
                                    gas%trans_einstB_u(gas%tr_cat(tr))) *      &
                                    linescale*grid%cell_gauss_a(nr_cell)

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
                                        ! analytical velocity distribution
                                        velo_dir_xyz = dot_product(get_velocity(model, &
                                                        pos_xyz_cell+          &
                                                        cell_d_l *             &
                                                        model%D_2obs(3, :, i_map) *      &
                                                        RK_c(k)),              &
                                                        model%D_2obs(3, :, i_map))

                                    ELSEIF (velo_type == 2 ) THEN
                                        ! linear interpolation of the velocity 
                                        velo_dir_xyz = interpolate_velo(       &
                                                            dot_product(grid%velo(nr_cell,:), &
                                                            model%D_2obs(3, :, i_map)),       &
                                                            dot_product(grid%velo(nr_cell_new,:), &
                                                            model%D_2obs(3, :, i_map)),       &
                                                            cell_d_l*RK_c(k)/d_l)
                                    END IF
                                    expo = -((gas%velo_channel(vch)-velo_dir_xyz)**2*      &
                                                     grid%cell_gauss_a2(nr_cell))
                                    gauss_val  = exp(expo)

                                    j_ges = j_ul  * gauss_val  + j_dust
                                    alpha_ges = alpha_ul * gauss_val  + alpha_dust

                                    RK_k(k) = ( -alpha_ges*(intensity(vch) + &
                                                dz*dot_product(RK_a(:,k),RK_k(:)) ) + j_ges)
                                END DO ! for all k

                                intensity_new  = intensity(vch) + dz*(dot_product(RK_b1(:),RK_k(:)))
                                intensity_new2 = intensity(vch) + dz*(dot_product(RK_b2(:),RK_k(:)))
                                
                                epsi= abs(intensity_new2-intensity_new)/&
                                         ( rel_err*abs(intensity_new) + abs_err)
                                ! dz_new = 0.9_r2*dz*exp(-log(epsi)*0.2_r2)
                                dz_new = 0.9*dz*epsi**(-0.2)
                                IF ( epsi .le. 1 ) THEN
                                    intensity(vch) = intensity_new
                                    pos_xyz_cell = pos_xyz_cell + cell_d_l*model%D_2obs(3, :, i_map)
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
                        END DO !vch
                    END DO !transitions , please do only one transition per run
                END IF

                pos_xyz = pos_xyz_new
                nr_cell = nr_cell_new
                dz_sum = dz_sum + d_l
            END DO !walk in z direction towards observer

        px_intensity(:,1) = intensity
        END IF
    END FUNCTION get_line_intensity_in_px
    
    FUNCTION get_continuum_intensity_along_one_path(grid, model, dust,        &
                              coor_map1, coor_map2, i_map)                    &
                              RESULT(intensity)
        
    IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        !------------------------------------------------------------- --------!
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        INTEGER, INTENT(IN)                              :: i_map

        REAL(KIND=r2)                                    :: ray_len
        REAL(KIND=r2)                                    :: dz
        REAL(KIND=r2)                                    :: dz_new
        REAL(KIND=r2)                                    :: j_dust
        REAL(KIND=r2)                                    :: alpha_dust

        REAL(KIND=r2)                                    :: d_l
        REAL(KIND=r2)                                    :: cell_d_l
        REAL(KIND=r1)                                    :: epsi
        REAL(KIND=r2)                                    :: cell_sum

        REAL(KIND=r2), DIMENSION(1:dust%n_lam)           :: intensity
        REAL(KIND=r2)                                    :: intensity_new
        REAL(KIND=r2)                                    :: intensity_new2
        REAL(KIND=r2)                                    :: j_ges

        REAL(KIND=r2)                                    :: alpha_ges

        REAL(KIND=r2), DIMENSION(1:3)                    :: pos_xyz
        REAL(KIND=r2), DIMENSION(1:3)                    :: pos_xyz_new
        REAL(KIND=r2), DIMENSION(1:3)                    :: pos_xyz_cell

        REAL(KIND=r2), DIMENSION(1:6)                    :: RK_k
        
         
        INTEGER                                          :: i_lam
        INTEGER                                          :: k
        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
    
        LOGICAL                                          :: kill_photon
        !----------------------------------------------------------------------!

        ! reset intensity in current ray
        j_dust     = 0.0
        j_ges      = 0.0
        alpha_dust = 0.0
        alpha_ges  = 0.0

        ray_len = get_ray_length(GetGridName(grid), model%r_ou,                &
                                 coor_map1, coor_map2)
        
        pos_xyz(:) = -model%D_2obs(3, :, i_map) * sqrt(ray_len) +              &
                      coor_map1*model%D_2obs(1, :, i_map) +                    &
                      coor_map2*model%D_2obs(2, :, i_map)

        intensity(:)      = 0.0  
        dz                = EPSILON(dz)

        IF ( ray_len .gt. 0.0_r2 ) THEN

            nr_cell      = get_cell_nr(grid, pos_xyz)
            DO WHILE (check_inside(nr_cell, grid))
                IF ( nr_cell==0 ) THEN
                    CALL path_skip( grid, pos_xyz,model%D_2obs(3, :, i_map), &
                                   pos_xyz_new,nr_cell_new,d_l)

                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new
                    IF (.not. check_inside(nr_cell, grid) ) EXIT
                END IF

                CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new,   &
                           d_l, kill_photon, model%D_2obs(3, :, i_map))

                ! At this point, we have an entrance point(pos_xyz) and exit point (pos_xyz_new) and 
                ! the length of the cell path (d_l) 
                IF ( sum(grid%grd_dust_density(nr_cell,:)) .gt. 1.0e-200_r2 ) THEN
                    DO i_lam = 1, dust%n_lam
                        alpha_dust = sum(grid%grd_dust_density(nr_cell,:) * &
                                     dust%C_ext(:,i_lam))
                        
                        j_dust   =  sum(grid%grd_dust_density(nr_cell,:) * &
                                    dust%C_abs(:,i_lam) * &
                                    planckhz(grid%t_dust(nr_cell, : ),&
                                    dust%nu(i_lam)))
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
                                RK_k(k) = ( -alpha_ges*(intensity(i_lam) + &
                                            dz*dot_product(RK_a(:,k),RK_k(:)) ) + j_ges)
                            END DO ! for all k
                    
                            intensity_new  = intensity(i_lam) + dz*(dot_product(RK_b1(:),RK_k(:)))
                            intensity_new2 = intensity(i_lam) + dz*(dot_product(RK_b2(:),RK_k(:)))
                            IF (intensity_new < 0.0_r2 .or. intensity_new2 < 0.0_r2 ) THEN
                                intensity_new  = 0.0_r2
                                intensity_new2 = 0.0_r2
                            END IF
                            epsi= abs(intensity_new2-intensity_new) /      &
                                     ( rel_err*abs(intensity_new) + abs_err)
                                     
                            dz_new = 0.9*dz*epsi**(-0.2)
                            IF ( epsi .le. 1 ) THEN
                                intensity(i_lam) = intensity_new
                                pos_xyz_cell = pos_xyz_cell+cell_d_l*model%D_2obs(3, :, i_map)
                                cell_sum = cell_sum + cell_d_l
                                dz = MIN(dz_new, 4*dz)
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
                END IF
                
                pos_xyz = pos_xyz_new
                nr_cell = nr_cell_new

            END DO !walk in z direction towards observer
        END IF
    END FUNCTION get_continuum_intensity_along_one_path

    RECURSIVE FUNCTION get_continuum_intensity_in_px(grid, model, dust, gas,  &
                              xxres, yyres, coor_map1, coor_map2, i_map,      &
                              depth, intensity_mid)                           &
                              RESULT(px_intensity)
        
    IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(Gas_TYP),INTENT(IN)                         :: gas
        
        !----------------------------------------------------------------------!
        REAL(KIND=r2),INTENT(IN)                         :: xxres
        REAL(KIND=r2),INTENT(IN)                         :: yyres
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        INTEGER, INTENT(IN)                              :: i_map
        INTEGER, INTENT(IN)                              :: depth
        REAL(KIND=r2), DIMENSION(1:dust%n_lam),OPTIONAL  :: intensity_mid
        REAL(KIND=r2), DIMENSION(1:dust%n_lam)           :: intensity_central
        
        REAL(KIND=r2), DIMENSION(4)                      :: coor_x_arr
        REAL(KIND=r2), DIMENSION(4)                      :: coor_y_arr

        REAL(KIND=r2), DIMENSION(1:dust%n_lam)           :: px_intensity
        REAL(KIND=r2), DIMENSION(1:dust%n_lam, 1:4)      :: intensity
        REAL(KIND=r2)                                    :: diff
        INTEGER                                          :: k

        !----------------------------------------------------------------------!
        intensity = 0.0_r2
        intensity_central = 0.0_r2
        px_intensity(:) = 0.0_r2
        !central Pixel
        
        IF (.not. present(intensity_mid)) THEN
            
            intensity_central = get_continuum_intensity_along_one_path(grid,   &
                            model, dust, coor_map1, coor_map2, i_map)
        ELSE
            intensity_central = intensity_mid
        END IF
        !Subpixel
        !Pixel 1
        coor_x_arr(1) = coor_map1+xxres*0.25_r2
        coor_y_arr(1) = coor_map2+yyres*0.25_r2
        
        !Pixel 2
        coor_x_arr(2) = coor_map1+xxres*0.25_r2
        coor_y_arr(2) = coor_map2-yyres*0.25_r2
        
        !Pixel 3
        coor_x_arr(3) = coor_map1-xxres*0.25_r2
        coor_y_arr(3) = coor_map2+yyres*0.25_r2
        
        !Pixel 4
        coor_x_arr(4) = coor_map1-xxres*0.25_r2
        coor_y_arr(4) = coor_map2-yyres*0.25_r2

        DO k = 1, 4
            intensity(:, k) = get_continuum_intensity_along_one_path(grid,     &
                               model, dust, coor_x_arr(k), coor_y_arr(k), i_map)
            diff = maxval(abs(intensity(:, k)                    &
                     - intensity_central)/(intensity_central+epsilon(1.0_r2)))
            IF (diff .gt. 0.1 .and. depth .lt. 4) EXIT
        END DO
        

        IF (diff .gt. 0.1 .and. depth .lt. 4) THEN
            DO k = 1, 4
                intensity(:, k) = get_continuum_intensity_in_px(grid, model, dust, gas,  &
                              xxres*0.5_r2, yyres*0.5_r2, coor_x_arr(k), coor_y_arr(k),    &
                              i_map, depth+1)
            END DO
        ELSE
!~         !$omp critical
!~         open(unit=1, file=TRIM('pixel_pos.dat'), &
!~                     action="write", status="unknown", form="formatted",position="append")
!~         DO k = 1,4
!~             WRITE(unit=1, fmt=*) coor_x_arr(k), coor_y_arr(k), depth
!~         END DO
!~         CLOSE(unit=1)
!~         !$omp end critical
        END IF
        ! only the first stokes component
        px_intensity(:) = sum(intensity, dim=2) * 0.25_r2

    END FUNCTION get_continuum_intensity_in_px
    
    RECURSIVE SUBROUTINE get_individual_needed_px(basics, grid, model, i,j,    &
                                    xxres,yyres, coor_map1, coor_map2,         &
                                    i_map, pixel_list)
    USE linkedlist_mod
    IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model

        TYPE(l_list),POINTER,INTENT(INOUT)               :: pixel_list
        TYPE(Pixel_TYP)                                  :: pixel_data
        
        !----------------------------------------------------------------------!
        REAL(KIND=r2),INTENT(IN)                         :: xxres
        REAL(KIND=r2),INTENT(IN)                         :: yyres
        REAL(KIND=r2),INTENT(IN)                         :: coor_map1
        REAL(KIND=r2),INTENT(IN)                         :: coor_map2
        INTEGER,INTENT(IN)                               :: i_map
        
        REAL(KIND=r2)                                    :: coor_x
        REAL(KIND=r2), DIMENSION(4)                      :: coor_x_arr
        REAL(KIND=r2)                                    :: coor_y
        REAL(KIND=r2), DIMENSION(4)                      :: coor_y_arr
        
        REAL(KIND=r2)                                    :: ray_len
        
        REAL(KIND=r2)                                    :: d_l
        REAL(KIND=r2)                                    :: ray_minA

        REAL(KIND=r2), DIMENSION(3)                      :: pos_xyz
        REAL(KIND=r2), DIMENSION(3)                      :: pos_xyz_new

        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
        INTEGER                                          :: k
        INTEGER,INTENT(IN)                               :: i,j

    
        LOGICAL                                         :: kill_photon, log_size
        !----------------------------------------------------------------------!
        ray_minA = 20.0
        log_size = .False.
        
        !Pixel 1
        coor_x_arr(1) = coor_map1+xxres*0.25_r2
        coor_y_arr(1) = coor_map2+yyres*0.25_r2
        
        !Pixel 2
        coor_x_arr(2) = coor_map1+xxres*0.25_r2
        coor_y_arr(2) = coor_map2-yyres*0.25_r2
        
        !Pixel 3
        coor_x_arr(3) = coor_map1-xxres*0.25_r2
        coor_y_arr(3) = coor_map2+yyres*0.25_r2
        
        !Pixel 4
        coor_x_arr(4) = coor_map1-xxres*0.25_r2
        coor_y_arr(4) = coor_map2-yyres*0.25_r2

        DO k = 1, 4
            coor_x = coor_x_arr(k)
            coor_y = coor_y_arr(k)
            ray_len = get_ray_length(GetGridName(grid), model%r_ou,            &
                                     coor_x, coor_y)
            pos_xyz(:) = -model%D_2obs(3, :, i_map) * sqrt(ray_len) +          &
                          coor_x * model%D_2obs(1, :, i_map) +  &
                          coor_y * model%D_2obs(2, :, i_map)

            IF ( ray_len .gt. 0.0_r2 .and. .not. log_size) THEN
                nr_cell  = get_cell_nr(grid, pos_xyz)
                DO WHILE (check_inside(nr_cell, grid))
                    IF ( nr_cell==0 ) THEN
                        ray_minA = MIN(ray_minA, grid%cell_minA(0))
                        CALL path_skip( grid, pos_xyz, model%D_2obs(3, :, i_map), &
                                       pos_xyz_new,nr_cell_new,d_l)
                        
                        pos_xyz = pos_xyz_new
                        nr_cell = nr_cell_new
                        
                        IF (.not. check_inside(nr_cell, grid) ) EXIT
                    END IF
                    ray_minA = MIN(ray_minA, grid%cell_minA(nr_cell))
                    IF ( xxres*yyres .gt. ray_minA) THEN
                        log_size = .True.
                    END IF
                    CALL path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, &
                               d_l, kill_photon, model%D_2obs(3, :, i_map) )

                    pos_xyz = pos_xyz_new
                    nr_cell = nr_cell_new

                END DO !walk in z direction towards observer
            END IF
        END DO

        IF ( log_size) THEN
            CALL        get_individual_needed_px(basics, grid, model,          &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                       &
                        coor_map1+xxres*0.25_r2, coor_map2+yyres*0.25_r2,      &
                        i_map, pixel_list)

            CALL        get_individual_needed_px(basics, grid, model,          &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                       &
                        coor_map1+xxres*0.25_r2, coor_map2-yyres*0.25_r2,      &
                        i_map, pixel_list)
                        
            CALL        get_individual_needed_px(basics, grid, model,          &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                       &
                        coor_map1-xxres*0.25_r2, coor_map2+yyres*0.25_r2,      &
                        i_map,pixel_list )
                        
            CALL        get_individual_needed_px(basics, grid, model,          &
                        i,j, xxres*0.5_r2, yyres*0.5_r2,                       &
                        coor_map1-xxres*0.25_r2, coor_map2-yyres*0.25_r2,      &
                        i_map, pixel_list)                                     
        ELSE
            pixel_data%pixel  =  (/i, j/)
            pixel_data%pos_xy  = (/coor_map1, coor_map2/)
            pixel_data%size_xy = (/xxres, yyres/)
            !$omp critical
            CALL AddElement(pixel_list, pixel_data)
            !$omp end critical
            CONTINUE
        END IF

    END SUBROUTINE get_individual_needed_px

    FUNCTION get_ray_length(grid_name, r_ou, coor_map1, coor_map2)        &
                                 RESULT(ray_len)
        ! calculates (length of required ray)**2
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(KIND=r2)                                     :: ray_len
        REAL(KIND=r2), INTENT(IN)                         :: r_ou
        REAL(KIND=r2), INTENT(IN)                         :: coor_map1
        REAL(KIND=r2), INTENT(IN)                         :: coor_map2
        CHARACTER(len=*), INTENT(IN)                      :: grid_name
        !----------------------------------------------------------------------!
        SELECT CASE(TRIM(grid_name))
        
        CASE('spherical')
            ray_len = r_ou**2-coor_map1**2-coor_map2**2 
        CASE('cylindrical')
            ray_len = r_ou**2-coor_map1**2-coor_map2**2
        CASE('cartesian')
            ! not sure if this ray is sufficient...
            ray_len = r_ou**2-coor_map1**2-coor_map2**2
        CASE DEFAULT
            PRINT *, 'selected coordinate system not found, simulation'
            STOP
        END SELECT
        ray_len = ray_len *(1.0_r2 - 1.0e6_r2*EPSILON(1.0_r2))
    END FUNCTION

    ELEMENTAL FUNCTION interpolate_velo(velo1, velo2, x) RESULT(velo_x)

        ! interpolate velocity 
        !
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2),INTENT(IN)                           :: velo1, velo2
        REAL(kind=r2)                                      :: velo_x
        REAL(kind=r2),INTENT(IN)                           :: x
        !----------------------------------------------------------------------!

        velo_x = (velo2-velo1)*x + velo1

    END FUNCTION interpolate_velo

END MODULE raytrace_mod
