MODULE fileio

    USE var_globalnew
    USE datatype
    USE grid_type
    USE basic_type
    USE dust_type
    USE model_type
    USE fluxes_type
    USE gas_type
!~     USE math_mod
    
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: vis_plane, save_ch_map, sv_temp, load_temp_dist, save_input, save_model, &
              save_boundaries, save_continuum_map
    !--------------------------------------------------------------------------!
CONTAINS    

    SUBROUTINE save_boundaries(grid,model,basics)
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Model_TYP), INTENT(IN)                  :: model
        !--------------------------------------------------------------------------!
        INTEGER                                      :: i
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=16)                            :: a
        CHARACTER(len=16)                            :: b
        CHARACTER(len=16)                            :: c
        !--------------------------------------------------------------------------!
        filename = TRIM(basics%path_results)//Getproname(basics)//'_cell_boundaries.dat'
        open(unit=1, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
            
        ! write header:
        write(unit=1,fmt='(A)') '# cell boundaries: '
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            write(unit=1,fmt='(A)') '# id            r            theta            phi'
            
        CASE('cylindrical')
            write(unit=1,fmt='(A)') '# id            r            phi              z'
            
        CASE('cartesian')
            write(unit=1,fmt='(A)') '# id            x            y                z'
        CASE DEFAULT   
            print *, 'selected coordinate system not found, save_boundaries'
            stop
        END SELECT
        write(unit=1,fmt='(A)') ''
        DO i = 0,maxval(grid%n)
        
            IF ( i .gt. grid%n(1) ) THEN
            
                WRITE(a,fmt='(A)') ''
            ELSE
                WRITE(a,fmt='(F16.8)') grid%co_mx_a(i)
            END IF
            
            IF ( i .gt. grid%n(2) ) THEN
                WRITE(b,fmt='(A)') ''
            ELSE
                WRITE(b,fmt='(F16.8)') grid%co_mx_b(i)
            END IF
            
            IF ( i .gt. grid%n(3) ) THEN
                WRITE(c,fmt='(A)') ''
            ELSE
                WRITE(c,fmt='(F16.8)') grid%co_mx_c(i)
            END IF
            
!~             write(unit=1,fmt='(I5,3(F16.8))') i,a,b,c
            write(unit=1,fmt='(I5,3(A))') i,a,b,c
        
        END DO
        CLOSE(unit=1)
!~         stop
    END SUBROUTINE  save_boundaries
        
    SUBROUTINE save_model(grid, basics)
    
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        !--------------------------------------------------------------------------!
        INTEGER                                      :: plane, pix
        
        REAL(kind=r2),DIMENSION(1:3)                 :: caco
        REAL(kind=r2)                                :: dxy
                
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=256)                           :: outname
        
        INTEGER                                      :: i_cell
        INTEGER                                      :: i,j,k
        !--------------------------------------------------------------------------!
        filename = TRIM(basics%path_results)//Getproname(basics)//'_model.dat'
        open(unit=1, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
        
        ! write header:
        write(unit=1,fmt='(12A)') '#cell  ','mid_x  ','mid_y  ','mid_z  ',  'no_dens_dust  ','no_dens_mol  ', &
                                    'no_dens_H  ', 'no_dens_Hp  ', 'no_dens_Ho  ','T_dust  ', 'T_gas  ', &
                                    'abs_velocity'
        
        DO i_cell=1,grid%n_cell
            
            write(unit=1,fmt='(I8,11(ES15.6E3))') i_cell, grid%cellmidcaco(i_cell,:), grid%grd_dust_density(i_cell,1), &
                                                   grid%grd_mol_density(i_cell) , &
                                                   grid%grd_col_density(i_cell,1:3), &
                                                   grid%t_dust(i_cell,1), grid%t_gas(i_cell), &
                                                   grid%absvelo(i_cell)
            
        END DO
        CLOSE(unit=1)
        
!~         stop
    END SUBROUTINE save_model
    
    

    SUBROUTINE vis_plane(grid, basics, model, plane, pix) 
    
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Model_TYP), INTENT(IN)                  :: model
        !--------------------------------------------------------------------------!
        INTEGER                                      :: plane
        INTEGER                                      :: pix
        
        REAL(kind=r2),DIMENSION(1:3)                 :: caco
        REAL(kind=r2)                                :: dxy
                
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        
        INTEGER                                      :: i_cell
        INTEGER                                      :: i_x, i_y
        
        !--------------------------------------------------------------------------!
        fileext = '.dat'

        IF ( plane ==  1) THEN                                       !view xz-plane
            filename = TRIM(basics%path_results)//Getproname(basics)//'_visual_xz'
            
        ELSE IF ( plane == 2) THEN                                   !view xy-plane
            filename = TRIM(basics%path_results)//Getproname(basics)//'_visual_xy'
            
        ELSE IF ( plane == 3) THEN                                   !view yz-plane
            filename = TRIM(basics%path_results)//Getproname(basics)//'_visual_yz'
        ELSE 
            print *, 'plane not specified'
            RETURN
        END IF
        
        outname = TRIM(filename)//fileext
        open(unit=1, file=TRIM(outname), &
            action="write", status="unknown", form="formatted")

        IF ( plane == 3) THEN
            print *,'visualize  -> yz - plane'
            dxy = 2.0*model%r_ou/pix
            write(unit=1,fmt='(I6.4,A)') pix ,'    # no of pixel of the map'
            write(unit=1,fmt=*) ''
        
            DO i_x=0,pix-1
                DO i_y=0,pix-1
                    caco(1) = 0.0_r2
                    caco(2) = -model%r_ou+(0.5 + i_x) * dxy
                    caco(3) = -model%r_ou+(0.5 + i_y) * dxy
                    IF ( check_inside(caco,grid,model) ) THEN
                        i_cell = get_cell_nr(grid,caco)
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(2), &
                            caco(3), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%absvelo(i_cell), &
                            grid%delta_t_dust(i_cell,1)
                    ELSE
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(2), &
                            caco(3), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2
                    
                    END IF
                END DO
            END DO 
            
        ELSEIF ( plane == 2) THEN
            print *,'visualize  -> xy - plane'
            dxy = 2.0*model%r_ou/pix
            write(unit=1,fmt='(I6.4,A)') pix ,'    # no of pixel of the map'
            write(unit=1,fmt=*) ''
        
            DO i_x=0,pix-1
                DO i_y=0,pix-1
                    caco(1) = -model%r_ou+(0.5 + i_x) * dxy
                    caco(2) = -model%r_ou+(0.5 + i_y) * dxy
                    caco(3) = 0.0_r2
                    IF ( check_inside(caco,grid,model) ) THEN
                        i_cell = get_cell_nr(grid,caco)
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(1), &
                            caco(2), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%absvelo(i_cell), &
                            grid%delta_t_dust(i_cell,1)
                    ELSE
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(1), &
                            caco(2), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2
                    
                    END IF
                END DO
            END DO 

        ELSE IF (plane == 1) THEN
            print *,'visualize  -> xz - plane'
            dxy = 2.0*model%r_ou/pix
            write(unit=1,fmt='(I6.4,A)') pix ,'    # no of pixel of the map'
            write(unit=1,fmt=*) ''
            
            DO i_x=0,pix-1
                DO i_y=0,pix-1
                    caco(1) = -model%r_ou+(0.5 + i_x) * dxy
                    caco(2) = 0.0_r2
                    caco(3) = -model%r_ou+(0.5 + i_y) * dxy
                    IF ( check_inside(caco,grid,model) ) THEN
                        i_cell = get_cell_nr(grid,caco)
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(1), &
                            caco(3), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%absvelo(i_cell), &
                            grid%delta_t_dust(i_cell,1)
                    ELSE
                        write(unit=1,fmt='(11(ES15.6E3))') &
                            caco(1), &
                            caco(3), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2
                    END IF
                END DO
            END DO     
        END IF

        close(unit=1)
        
        
    END SUBROUTINE vis_plane
    
    SUBROUTINE sv_temp(basics, grid)
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        !--------------------------------------------------------------------------!
    
        CALL sv_temp_x(basics, grid)
!~         Call sv_temp_dist(basics, grid)
    END SUBROUTINE sv_temp
    
     SUBROUTINE load_temp_dist(basics, grid,tmp_found)
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(INOUT)                :: grid
        !--------------------------------------------------------------------------!
                
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        CHARACTER(len=256)                           :: waste, coord
        
        INTEGER                                      :: i_cell, i
        INTEGER                                      :: a,b,c
        LOGICAL, INTENT(INOUT)                       :: tmp_found
        !--------------------------------------------------------------------------!
        
        
        print '(2A)', '  loading temperature distribution from: ', TRIM(basics%pronam_old)
        fileext = '.dat'
        filename = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_temp_dist'
        outname = TRIM(filename)//fileext
        open(unit=1, file=TRIM(outname), &
            action="read", status="unknown", form="formatted")
        
        READ(unit=1,fmt=*) coord
        READ(unit=1,fmt=*) a, waste
        READ(unit=1,fmt=*) b, waste
        READ(unit=1,fmt=*) c, waste
        
        tmp_found = .TRUE.
        IF (.not. (TRIM(coord) .eq. TRIM(GetGridName(grid)))) THEN
            tmp_found = .FALSE.
        ELSE IF ( .not. a .eq. grid%n(1)) THEN
            tmp_found = .FALSE.
        ELSE IF ( .not. b .eq. grid%n(2)) THEN
            tmp_found = .FALSE.
        ELSE IF ( .not. c .eq. grid%n(3)) THEN
            tmp_found = .FALSE.
        END IF
        
        IF (tmp_found) THEN
            DO i_cell=0,grid%n_cell
                READ(unit=1,fmt=*) i,grid%t_dust(i_cell,1)
            END DO
            print * ,'  temperature successfully loaded '
        END IF
        CLOSE(unit=1)
        
    END SUBROUTINE load_temp_dist
       
    SUBROUTINE sv_temp_dist(basics, grid)
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        !--------------------------------------------------------------------------!
                
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        
        INTEGER                                      :: i_cell
        !--------------------------------------------------------------------------!
        
        fileext = '.dat'
        filename = TRIM(basics%path_results)//Getproname(basics)//'_temp_dist'
        outname = TRIM(filename)//fileext
        open(unit=1, file=TRIM(outname), &
            action="write", status="unknown", form="formatted")
        
        WRITE(unit=1,fmt=*) TRIM(GetGridName(grid))//'   #coordinate system name'
        WRITE(unit=1,fmt='(I6.4,A)') grid%n(1) ,'    #number of first component '
        WRITE(unit=1,fmt='(I6.4,A)') grid%n(2) ,'    #number of second component '
        WRITE(unit=1,fmt='(I6.4,A)') grid%n(3) ,'    #number of third component '
        WRITE(unit=1,fmt=*) ''
        
        
        DO i_cell=0,grid%n_cell
            WRITE(unit=1,fmt=*) i_cell,grid%t_dust(i_cell,1)
        END DO
        CLOSE(unit=1)
    END SUBROUTINE sv_temp_dist
    
    
    SUBROUTINE sv_temp_x(basics, grid)
    
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        !--------------------------------------------------------------------------!
        INTEGER                                      ::  i_a, i_b, i_c, i_dust, i_cell
        REAL(kind=r2)                                ::  hd_r, hd_tem
        
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        !--------------------------------------------------------------------------!
        fileext = '.dat'
        filename = TRIM(basics%path_results)//Getproname(basics)//'_temp_x'
        outname = TRIM(filename)//fileext
        open(unit=1, file=TRIM(outname), &
            action="write", status="unknown", form="formatted")
        SELECT CASE(GetGridName(grid))
            
            CASE('spherical')
                !print *, 'sp temp_x'
                i_b = ceiling(real(grid%n(2))/2.0)
                i_c = 1
                
                ! tbd: for all dust species (here: for first only)
                i_dust = 1
                DO i_a = 1, grid%n(1)
                    hd_r   = (grid%co_mx_a(i_a-1) + grid%co_mx_a(i_a))/2.0_r2
                    i_cell = grid%cell_idx2nr(i_a,i_b,i_c)

                    IF (grid%grd_dust_density(i_cell,i_dust) .eq. 0.0) THEN
                        hd_tem = 0.0_r2
                    ELSE
                        hd_tem = grid%t_dust(i_cell,i_dust)
                    END IF
                   write(unit=1,fmt=*) hd_r, hd_tem
                END DO
    
            CASE('cylindrical')
                !print *, 'cy temp_x'
                i_b = 1
                i_c = ceiling(real(grid%n(3))/2.0)
                
                ! tbd: for all dust species (here: for first only)
                i_dust = 1
                DO i_a = 1, grid%n(1)
                    hd_r   = (grid%co_mx_a(i_a-1) + grid%co_mx_a(i_a))/2.0_r2
                    i_cell = grid%cell_idx2nr(i_a,i_b,i_c)

                    IF (grid%grd_dust_density(i_cell,i_dust) .eq. 0.0) THEN
                        hd_tem = 0.0_r2
                    ELSE
                        hd_tem = grid%t_dust(i_cell,i_dust)
                    END IF
                   write(unit=1,fmt=*) hd_r, hd_tem
                END DO

            CASE('cartesian')
                print *, 'TbD, not finished yet, save temp_x'

                
            CASE DEFAULT
                print *, 'selected coordinate system not found, save temp_x '
        END SELECT
        close(unit=1)
    
    
    END SUBROUTINE sv_temp_x
    
    SUBROUTINE save_input(basics,input_file)
    ! not used at the moment
    IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        CHARACTER(len=*)                             :: input_file
        CHARACTER(len=256)                           :: command

        !--------------------------------------------------------------------------!
      
        command = 'cp '//TRIM(input_file)//' '//TRIM(basics%path_results)//Getproname(basics)//'_input_file.dat'
        ! I'm not sure if the SYSTEM command is a good choice.
        CALL SYSTEM(command)
        
    END SUBROUTINE save_input
    
    SUBROUTINE save_continuum_map(model, basics, dust, fluxes)
        
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                  :: model
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Dust_TYP), INTENT(IN)                   :: dust
        TYPE(Fluxes_TYP), INTENT(IN)                 :: fluxes
        !--------------------------------------------------------------------------!
        INTEGER                                      :: i_lam, i, j
        !--------------------------------------------------------------------------!
        
        
        OPEN(unit=1, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'continuum_map.dat', &
            action="write", status="unknown", form="formatted")
        
        OPEN(unit=2, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'continuum_sed.dat', &
            action="write", status="unknown", form="formatted")
        
        !write file header
        WRITE(unit=1,fmt='(I6.4,A)') dust%n_lam, '     #number of wavelengths'
        WRITE(unit=1,fmt='(I6.4,A)') 2*model%n_bin_map+1 , '   #size of map'
        WRITE(unit=1,fmt=*) ''
        
        WRITE(unit=2,fmt='(A)') '#wavelength [m]    flux [Jy]'
        WRITE(unit=2,fmt=*) ''
            
        
        !write results
            
        DO i_lam =1, dust%n_lam
!~             print *, i_lam
            WRITE (unit=2,fmt='(2(ES15.6E3))') dust%lam(i_lam), sum(fluxes%continuum_map(:,:,i_lam))

            WRITE(unit=1,fmt=*) ''
            WRITE (unit=1,fmt='(ES15.6E3,A)') dust%lam(i_lam), '   #wavelength'
            WRITE(unit=1,fmt=*) ''
            DO i = 0, 2*model%n_bin_map
                DO j = 0, 2*model%n_bin_map   
                    WRITE(unit=1,fmt='(I5,I5,1(ES15.6E3))')   &
                            i, &
                            j, &
                            fluxes%continuum_map(i,j,i_lam)
                END DO
            END DO
        END DO
            
        CLOSE(unit=1)
        CLOSE(unit=2)

    END SUBROUTINE save_continuum_map
    
    SUBROUTINE save_ch_map(model, basics, gas, fluxes)
        
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                  :: model
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Gas_TYP), INTENT(IN)                    :: gas
        TYPE(Fluxes_TYP), INTENT(IN)                 :: fluxes
        !--------------------------------------------------------------------------!
        INTEGER                                      :: tr, vch, i, j
        !--------------------------------------------------------------------------!
        
        
        OPEN(unit=1, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'velo_ch_map.dat', &
            action="write", status="unknown", form="formatted")
            
        OPEN(unit=2, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'velo_ch_mapsum.dat', &
            action="write", status="unknown", form="formatted")
            
        OPEN(unit=3, file=TRIM(basics%path_results)//Getproname(basics)//'_'//'velo_ch_mapint.dat', &
            action="write", status="unknown", form="formatted")
        
        !write file header
        
        WRITE(unit=1,fmt='(I6.4,A)') gas%i_vel_chan*2+1, '     #number of velo-channels'
        WRITE(unit=1,fmt='(I6.4,A)') 2*model%n_bin_map+1 , '   #size of map'
        WRITE(unit=1,fmt=*) ''
        
        WRITE(unit=3,fmt='(I6.4,A)') 2*model%n_bin_map+1 , '   #size of map'
        WRITE(unit=3,fmt=*) ''
        
        
        !write results
            
        DO tr=1,gas%n_tr
            DO vch=-gas%i_vel_chan, gas%i_vel_chan
            
                WRITE(unit=2,fmt='(2(ES15.6E3))')   gas%velo_channel(vch), &
                                                    sum(fluxes%channel_map(:,:,vch,tr))
                                                    
                WRITE(unit=1,fmt=*) ''
                WRITE (unit=1,fmt='(ES15.6E3,A)') gas%velo_channel(vch), '   #velo-channel'
                WRITE(unit=1,fmt=*) ''
                DO i = 0, 2*model%n_bin_map
                    DO j = 0, 2*model%n_bin_map   
                        WRITE(unit=1,fmt='(I5,I5,1(ES15.6E3))')   &
                                i, &
                                j, &
                                fluxes%channel_map(i,j,vch,tr)
                        
                        IF ( vch == 0 ) THEN 
                            WRITE(unit=3,fmt='(I5,I5,1(ES15.6E3))')   &
                                        i, &
                                        j, &
                                        sum(fluxes%channel_map(i,j,:,tr))/ &
                                        real(gas%i_vel_chan, kind=r2)*gas%vel_max*1e-3
                        END IF
                    END DO
                    !WRITE(unit=1,fmt=*) ''
                    !WRITE(unit=3,fmt=*) ''
                END DO
                !WRITE(unit=1,fmt=*) ''
                !WRITE(unit=3,fmt=*) ''
            END DO
            WRITE(unit=2,fmt=*) ' '
            WRITE(unit=2,fmt=*) ' '
            
        END DO
        CLOSE(unit=1)
        CLOSE(unit=2)
        CLOSE(unit=3)
    END SUBROUTINE save_ch_map
    
!~     ! save map (entire stokes vector)
!~     SUBROUTINE sv_stokes(basics,grid,model,fluxes,dust)
!~ 
!~         IMPLICIT NONE
!~         !--------------------------------------------------------------------------!
!~         TYPE(Basic_TYP), INTENT(IN)                  :: basics
!~         TYPE(Grid_TYP), INTENT(IN)                   :: grid
!~         TYPE(Model_TYP), INTENT(IN)                  :: model
!~         TYPE(Fluxes_TYP), INTENT(IN)                 :: fluxes
!~         TYPE(Dust_TYP), INTENT(IN)                   :: dust
!~         !--------------------------------------------------------------------------!
!~         INTEGER                                      :: i1
!~         INTEGER                                      :: i2
!~         INTEGER                                      :: i_stokes
!~         INTEGER                                      :: i_map
!~         INTEGER                                      :: i_lam_map
!~         !--------------------------------------------------------------------------!
!~         ! ---
!~         DO i_stokes=1,1!4
!~             print *, "saving stokes parameter map(s) to ", &
!~             TRIM(basics%path_results)//Getproname(basics)//".stokes_map."//fluxes%stokes_ext(i_stokes)
!~ 
!~             open(unit=1, file=TRIM(basics%path_results)//Getproname(basics)// &
!~                 ".stokes_map."//fluxes%stokes_ext(i_stokes), &
!~                 action="write", status="unknown", form="formatted")
!~        
!~                write(unit=1,fmt=*) "# Stokes ", fluxes%stokes_ext(i_stokes)
!~                write(unit=1,fmt=*) "# --------"
!~                write(unit=1,fmt=*) model%n_map,     "# n_map"
!~                write(unit=1,fmt=*) dust%n_lam_map, "# n_lam_map"
!~                write(unit=1,fmt=*) model%n_bin_map, "# n_bin_map"
!~                write(unit=1,fmt=*) "# --------"
!~                write(unit=1,fmt=*) "# flux [Jy], intensity [W/m/sr]"
!~                write(unit=1,fmt=*) "# --------"
!~ 
!~             do i_map=1, model%n_map
!~                 write(unit=1,fmt=*) rad2grad(model%th_map(i_map)), "# theta (map) [deg]"
!~                 write(unit=1,fmt=*) rad2grad(model%ph_map(i_map)), "# phi   (map) [deg]"
!~              sel molecule mass   :  4.4827E-07 M_sun

!~                 do i_lam_map=1, 1 !TBD!~!
!~                     write(unit=1,fmt=*) "# --------"
!~                     write(unit=1,fmt=*) dust%lam(dust%num_lam_map(i_lam_map))*1.0e+6_r2, "# wavelength [micron]"
!~                     write(unit=1,fmt=*) "# --------"
!~              
!~                     do i1=-model%n_bin_map, model%n_bin_map
!~                         do i2=-model%n_bin_map, model%n_bin_map
!~                             write(unit=1,fmt=*) &
!~                                 cnv_lum2Jy( fluxes%stokes_map( i_stokes,i_map,i_lam_map,i1,i2 ), &
!~                                 dust%lam(dust%num_lam_map(i_lam_map)), i_map,model%distance ),&
!~                                 cnv_Wmsr(   fluxes%stokes_map( i_stokes,i_map,i_lam_map,i1,i2 ), i_map )
!~                         end do
!~                     end do
!~           
!~                 end do
!~             end do
!~ 
!~             close(unit=1)
!~         end do
!~ 
!~     end subroutine sv_stokes
  
    

END MODULE fileio
