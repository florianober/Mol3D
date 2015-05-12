MODULE fileio

    USE var_global
    USE datatype
    USE grid_type
    USE basic_type
    USE dust_type
    USE model_type
    USE fluxes_type
    USE gas_type
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
        CHARACTER(len=30), PARAMETER :: file_a  = "input/grid/a_boundaries.dat"
        CHARACTER(len=30), PARAMETER :: file_b  = "input/grid/b_boundaries.dat"
        CHARACTER(len=30), PARAMETER :: file_c  = "input/grid/c_boundaries.dat"
    !--------------------------------------------------------------------------!
    PUBLIC :: vis_plane, vis_plane_fits, save_ch_map, sv_temp, save_input,     &
              save_model, save_boundaries, read_boundaries, read_no_cells,     &
              save_continuum_map, read_model
    !--------------------------------------------------------------------------!

CONTAINS

     SUBROUTINE read_no_cells(n_a, n_b, n_c, b_file_a, b_file_b, b_file_c)
        !----------------------------------------------------------------------!
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_a
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_b
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_c
        
        INTEGER, INTENT(OUT)                :: n_a
        INTEGER, INTENT(OUT)                :: n_b
        INTEGER, INTENT(OUT)                :: n_c
        !----------------------------------------------------------------------!
        
        IF (PRESENT(b_file_a)) THEN
            open(unit=1, file=b_file_a, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_a, &
                action="read", status="unknown", form="formatted")
        END IF

        read(unit=1, fmt=*) n_a
        close(unit=1)
        
        IF (PRESENT(b_file_b)) THEN
            open(unit=1, file=b_file_b, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_b, &
                action="read", status="unknown", form="formatted")
        END IF

        read(unit=1, fmt=*) n_b
        close(unit=1)

        IF (PRESENT(b_file_c)) THEN
            open(unit=1, file=b_file_c, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_c, &
                action="read", status="unknown", form="formatted")
        END IF

        read(unit=1, fmt=*) n_c
        close(unit=1)

    END SUBROUTINE read_no_cells

    SUBROUTINE read_boundaries(grid, b_file_a, b_file_b, b_file_c)
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(INOUT)                :: grid
        !----------------------------------------------------------------------!
        INTEGER                                      :: i, i_abc
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_a
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_b
        CHARACTER(len=*), INTENT(IN), OPTIONAL       :: b_file_c
        CHARACTER(len=252)                           :: fmt_kind
        REAL(kind=r2)                                :: value_in

        !----------------------------------------------------------------------!

        ! 1 coordinate
        IF (PRESENT(b_file_a)) THEN
            open(unit=1, file=b_file_a, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_a, &
                action="read", status="unknown", form="formatted")
        END IF

        ! read header lines
        read(unit=1,fmt=*) i_abc ! number of r cells
        IF ( i_abc /= grid%n(1)  )  THEN
            PRINT *, "ERROR: No of 1 coordinate is not consistent"
            STOP
        END IF
        read(unit=1,fmt=*) fmt_kind
        DO i = 0, grid%n(1)
            read(unit=1,fmt=*) value_in
            IF ( TRIM(fmt_kind) == 'log10' ) THEN
                grid%co_mx_a(i) = 10.0**(value_in)
            ELSEIF (TRIM(fmt_kind) == 'noscale' ) THEN
                grid%co_mx_a(i) = value_in
            ELSE
                PRINT *, "ERROR: 1 coordinate has no known kind"
                STOP
            END IF
        END DO
        close(unit=1)

        ! ---
        ! 2 coordinate

        IF (PRESENT(b_file_b)) THEN
            open(unit=1, file=b_file_b, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_b, &
                action="read", status="unknown", form="formatted")
        END IF

        ! read header lines
        read(unit=1,fmt=*) i_abc ! number of th cells
        IF ( i_abc /= grid%n(2)  ) THEN
            PRINT *, "ERROR: No of 2 coordinate is not consistent"
            STOP
        END IF
        read(unit=1,fmt=*) fmt_kind

        DO i = 0, grid%n(2)
            read(unit=1,fmt=*) value_in
            IF ( TRIM(fmt_kind) == 'log10' ) THEN
                grid%co_mx_b(i) = 10.0**(value_in)
            ELSEIF (TRIM(fmt_kind) == 'noscale' ) THEN
                grid%co_mx_b(i) = value_in
            ELSE
                PRINT *, "ERROR: 2 coordinate has no known kind"
                STOP
            END IF
        END DO
        close(unit=1)

        ! ---
        ! 3 coordinate
        IF (PRESENT(b_file_c)) THEN
            open(unit=1, file=b_file_c, &
                action="read", status="unknown", form="formatted")
        ELSE
            open(unit=1, file=file_c, &
                action="read", status="unknown", form="formatted")
        END IF

        ! read header lines
        read(unit=1,fmt=*) i_abc ! number of th cells
        IF ( i_abc /= grid%n(3)  ) THEN
            PRINT *, "ERROR: No of 3 coordinate is not consistent"
            STOP
        END IF
        read(unit=1,fmt=*) fmt_kind

        DO i = 0, grid%n(3)
            read(unit=1,fmt=*) value_in
            IF ( TRIM(fmt_kind) == 'log10' ) THEN
                grid%co_mx_c(i) = 10.0**(value_in)
            ELSEIF (TRIM(fmt_kind) == 'noscale' ) THEN
                grid%co_mx_c(i) = value_in
            ELSE
                PRINT *, "ERROR: 3 coordinate has no known kind"
                STOP
            END IF
        END DO
        close(unit=1)

    END SUBROUTINE read_boundaries

    SUBROUTINE save_boundaries(grid, basics)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        !----------------------------------------------------------------------!
        INTEGER                                      :: i
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=16)                            :: a
        CHARACTER(len=16)                            :: b
        CHARACTER(len=16)                            :: c
        !----------------------------------------------------------------------!
        filename = TRIM(basics%path_results)//Getproname(basics)//'_cell_boundaries.dat'
        open(unit=1, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
            
        filename = TRIM(basics%path_results)//Getproname(basics)//'_a_boundaries.dat'
        open(unit=2, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
            
        filename = TRIM(basics%path_results)//Getproname(basics)//'_b_boundaries.dat'
        open(unit=3, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
            
        filename = TRIM(basics%path_results)//Getproname(basics)//'_c_boundaries.dat'
        open(unit=4, file=TRIM(filename), &
            action="write", status="unknown", form="formatted")
            
        ! write header:
        ! for the general file
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
        ! for individual files
        write(unit=2,fmt=*) grid%n(1)
        write(unit=2,fmt='(A)') 'noscale'
        write(unit=3,fmt=*) grid%n(2)
        write(unit=3,fmt='(A)') 'noscale'
        write(unit=4,fmt=*) grid%n(3)
        write(unit=4,fmt='(A)') 'noscale'
        
        ! write boundaries
        DO i = 0,maxval(grid%n)
        
            IF ( i .gt. grid%n(1) ) THEN
            
                WRITE(a,fmt='(A)') ''
            ELSE
                WRITE(a,fmt='(F16.8)') grid%co_mx_a(i)
                WRITE(unit=2,fmt=*) grid%co_mx_a(i)
            END IF
            
            IF ( i .gt. grid%n(2) ) THEN
                WRITE(b,fmt='(A)') ''
            ELSE
                WRITE(b,fmt='(F16.12)') grid%co_mx_b(i)
                WRITE(unit=3,fmt=*) grid%co_mx_b(i)
            END IF
            
            IF ( i .gt. grid%n(3) ) THEN
                WRITE(c,fmt='(A)') ''
            ELSE
                WRITE(c,fmt='(F16.12)') grid%co_mx_c(i)
                WRITE(unit=4,fmt=*) grid%co_mx_c(i)
            END IF
            
            write(unit=1,fmt='(I5,3(A))') i,a,b,c
        
        END DO
        CLOSE(unit=1)
        CLOSE(unit=2)
        CLOSE(unit=3)
        CLOSE(unit=4)

    END SUBROUTINE  save_boundaries

    SUBROUTINE read_model(grid, filename)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(INOUT)                :: grid
        !----------------------------------------------------------------------!
        INTEGER                                      :: i_cell
        INTEGER                                      :: n_dust
        INTEGER                                      :: k
        INTEGER                                      :: entries
        INTEGER                                      :: sta, u, bs, rw, nfound
        INTEGER,DIMENSION(2)                         :: naxes
        REAL(kind=r2), DIMENSION(:), ALLOCATABLE     :: line
        CHARACTER(len=*), INTENT(IN)                 :: filename
        CHARACTER(len=256)                           :: comment
        LOGICAL                                      :: anyf
        !----------------------------------------------------------------------!
        sta = 0
        ! get a new u(nit) number
        call ftgiou(u,sta)
        ! init fits file
        call ftopen(u,TRIM(filename),0,bs,sta)
        ! check axis
        call ftgknj(u, 'NAXIS', 1, 2, naxes, nfound, sta)
        IF (nfound /= 2) THEN
            print *, "ERROR, NAXIS keyword not found in model fits file"
            STOP
        END IF
        ! check number of dust species
        call ftgkyj(u, 'N_DUST', n_dust, comment, sta)

        IF (sta == 1) THEN
            print *, "ERROR, N_DUST keyword not found in model fits file"
            STOP
        END IF
        entries = n_dust*2 + 11
        
        IF (n_dust /= grid%nh_n_dust .or. naxes(1) /= entries) THEN
            print *, "ERROR, the number of dust species is not consistent in model fits file"
            STOP
        END IF
        ALLOCATE(line(1:entries))
        IF (naxes(2) /= grid%n_cell) THEN 
            print *, "ERROR, the number of grid cells in model fits file   &
                     &is not correct"
            STOP
        END IF
        k = grid%n_cell/100
        DO i_cell = 1, grid%n_cell
            IF (modulo(i_cell,k) == 0 .or. i_cell == grid%n_cell) THEN
                WRITE (*,'(A,I3,A)') ' | | ',int(i_cell/                   &
                       real(grid%n_cell)*100.0),' % done'//char(27)//'[A'
            END IF
        
            call ftgpvd(u,1,(i_cell-1)*entries + 1,entries,1e-200_r2,line,anyf,sta)
            ! set dust density and H2, the observed molecules can be set by
            ! the model module
            IF (sta /= 0) THEN
                print *,'ERROR'
                print *,line
                print *,i_cell
                stop
            END IF
            grid%grd_dust_density(i_cell, 1:n_dust) = line(4:3+n_dust)

            grid%grd_col_density(i_cell,1:3) = line(n_dust+5:n_dust+7)

            ! set temperature
            grid%t_dust(i_cell,:) = line(n_dust+8:2*n_dust+7)
            grid%t_gas(i_cell)    = line(2*n_dust + 9)
            ! set velocity
            grid%velo(i_cell,:) = line(2*n_dust+9:2*n_dust+11)
        END DO
        ! close the fits file
        call ftclos(u, sta)
        ! free the (u)nit number
        call ftfiou(u, sta)
        DEALLOCATE(line)

    END SUBROUTINE read_model
        
    SUBROUTINE save_model(grid, basics)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        !----------------------------------------------------------------------!
        INTEGER                                      :: i_cell
        INTEGER                                      :: k
        INTEGER                                      :: entries
        INTEGER                                      :: sta, u
        !----------------------------------------------------------------------!
        ! save the hole model (after temperature calculation)
        k = grid%n_cell/100
        sta = 0
        entries = 11 + 2*grid%nh_n_dust

        ! get a new u(nit) number
        call ftgiou(u,sta)
        ! init fits file
        call ftinit(u,'!'//TRIM(basics%path_results)//Getproname(basics)//'_model.fits',1,sta)
        ! write header
        call ftphpr(u,.true.,-64,2,(/entries, grid%n_cell/),0,1,.true.,sta)
        ! write array to fits file
        DO i_cell=1,grid%n_cell
            IF (modulo(i_cell,k) == 0 .or. i_cell == grid%n_cell) THEN
                WRITE (*,'(A,I3,A)') ' | | ',                                  &
                         int(i_cell/real(grid%n_cell)*100.0),                  &
                         ' % done'//char(27)//'[A'
            END IF
            call ftpprd(u,1,(i_cell-1)*entries+ 1,entries,                     &
                        (/ grid%cellmidcaco(i_cell,:),                         &
                        grid%grd_dust_density(i_cell,:),                       &
                        grid%grd_mol_density(i_cell) ,                         &
                        grid%grd_col_density(i_cell,1:3),                      &
                        REAL(grid%t_dust(i_cell,:),kind=r2),                   &
                        REAL(grid%t_gas(i_cell),kind=r2),                      &
                        REAL(grid%velo(i_cell,:),kind=r2)/), sta)
        END DO
        ! add keywords
        call ftpkyj(u,'N_DUST', grid%nh_n_dust, 'Number of dust species', sta)
        ! close the fits file
        call ftclos(u, sta)
        ! free the (u)nit number
        call ftfiou(u, sta)
    END SUBROUTINE save_model

    SUBROUTINE vis_plane(grid, basics, model, plane, pix)
        !old one, has been replaced by the fits version

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Model_TYP), INTENT(IN)                  :: model
        !----------------------------------------------------------------------!
        INTEGER                                      :: plane
        INTEGER                                      :: pix
        
        REAL(kind=r2),DIMENSION(1:3)                 :: caco
        REAL(kind=r2)                                :: dxy
                
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        
        INTEGER                                      :: i_cell
        INTEGER                                      :: i_x, i_y
        
        !----------------------------------------------------------------------!
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
            print *,'| |  -> yz - plane'
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
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(2), &
                            caco(3), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%velo(i_cell,:)
                    ELSE
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(2), &
                            caco(3), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2
                    
                    END IF
                END DO
            END DO 
            
        ELSEIF ( plane == 2) THEN
            print *,'| |  -> xy - plane'
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
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(1), &
                            caco(2), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%velo(i_cell,:)
                    ELSE
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(1), &
                            caco(2), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2
                    
                    END IF
                END DO
            END DO 

        ELSE IF (plane == 1) THEN
            print *,'| |  -> xz - plane'
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
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(1), &
                            caco(3), &
                            grid%grd_dust_density(i_cell,1), &
                            grid%grd_mol_density(i_cell), &
                            grid%grd_col_density(i_cell,1:3), &
                            grid%t_dust(i_cell,1), &
                            grid%t_gas(i_cell), &
                            grid%velo(i_cell,:)
                    ELSE
                        write(unit=1,fmt='(12(ES15.6E3))') &
                            caco(1), &
                            caco(3), &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2, &
                            0.0_r2, &
                            0.0_r2, &
                            0.0_r2,0.0_r2,0.0_r2
                    END IF
                END DO
            END DO     
        END IF
        close(unit=1)

    END SUBROUTINE vis_plane
    
    SUBROUTINE vis_plane_fits(grid, basics, model, plane, pix) 

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Model_TYP), INTENT(IN)                  :: model
        !----------------------------------------------------------------------!
        INTEGER                                      :: plane
        INTEGER                                      :: pix
        
        REAL(kind=r2),DIMENSION(1:3)                 :: caco
        REAL(kind=r2)                                :: dxyz
        REAL(kind=r2), DIMENSION(:,:,:), ALLOCATABLE :: map_out
                
        CHARACTER(len=256)                           :: filename
        CHARACTER(len=5)                             :: fileext
        CHARACTER(len=256)                           :: outname
        
        INTEGER                                      :: i_cell, i_entry
        INTEGER                                      :: start
        INTEGER                                      :: i_x, i_y

        INTEGER                                      :: k
        INTEGER                                      :: entries
        INTEGER                                      :: sta, u
        !----------------------------------------------------------------------!
        ! save the hole model (after temperature calculation)
        k = grid%n_cell/100
        sta = 0
        entries = 8 + 2*grid%nh_n_dust

        ALLOCATE( map_out(0:pix-1,0:pix-1,1:entries))
        ! generate the output map
        dxyz = 2.0*model%r_ou/pix
        DO i_y=0, pix-1
            DO i_x=0,pix-1
                IF (plane == 1) THEN
                    caco(1) = -model%r_ou+(0.5 + i_x) * dxyz
                    caco(2) =  0.0_r2
                    caco(3) = -model%r_ou+(0.5 + i_y) * dxyz
                ELSE IF (plane == 2) THEN
                    caco(1) = -model%r_ou+(0.5 + i_x) * dxyz
                    caco(2) = -model%r_ou+(0.5 + i_y) * dxyz
                    caco(3) = 0.0_r2
                ELSE IF (plane == 3) THEN
                    caco(1) = 0.0_r2
                    caco(2) = -model%r_ou+(0.5 + i_x) * dxyz
                    caco(3) = -model%r_ou+(0.5 + i_y) * dxyz
                END IF
                IF ( check_inside(caco, grid, model) ) THEN
                    i_cell = get_cell_nr(grid, caco)
                    map_out(i_x,i_y,:) = (/grid%grd_dust_density(i_cell,:),&
                                grid%grd_mol_density(i_cell),              &
                                grid%grd_col_density(i_cell,1:3),          &
                                REAL(grid%t_dust(i_cell,:), kind=r2),      &
                                REAL(grid%t_gas(i_cell), kind=r2),         &
                                REAL(grid%velo(i_cell,:), kind=r2)/)
                ELSE
                    map_out(i_x,i_y,:) = 0.0_r2
                END IF
            END DO
        END DO

        ! now save the map
        fileext = '.fits'

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

        outname = TRIM(filename)//fileext//'.gz'
        ! get a new u(nit) number
        call ftgiou(u, sta)
        
        ! init fits file
        call ftinit(u,'!'//outname,1,sta)

        ! write header
        call ftphpr(u,.true.,-64,3,(/pix, pix, entries/),0,1,.true.,sta)
        ! write array to fits file
        DO i_entry = 1, entries
            start = 1 + (i_entry-1) * pix**2
            call ftpprd(u,1, start, pix**2,   &
                    (/map_out(:, :, i_entry)/),sta)
        END DO
        ! add essential keywords
        CALL add_essential_fits_keys(u, grid%nh_n_dust, pix, model%r_ou,&
                                     GetModelName(model))
        ! close the fits file
        call ftclos(u, sta)
        ! free the (u)nit number
        call ftfiou(u, sta)

        DEALLOCATE(map_out)
    END SUBROUTINE vis_plane_fits

    SUBROUTINE add_essential_fits_keys(u, n_dust, pix, r_max, ref_u_str)
        IMPLICIT NONE
        !----------------------------------------------------------------------!

        INTEGER, INTENT(IN)                          :: u
        INTEGER, INTENT(IN)                          :: n_dust
        INTEGER, INTENT(IN)                          :: pix

        REAL(kind=r2), INTENT(IN)                    :: r_max
        
        CHARACTER(len=*), INTENT(IN)                 :: ref_u_str

        INTEGER                                      :: sta
        !----------------------------------------------------------------------!
        sta = 0

        ! add the keys
        CALL ftpkyj(u,'N_DUST', n_dust, 'Number of dust species', sta)
        CALL ftpkys(u,'CTYPE1', 'LINEAR', 'Unit 1 type', sta)
        CALL ftpkyd(u,'CRVAL1', 0.0, 1, 'Value 1', sta)
        CALL ftpkyj(u,'CRPIX1', int((pix-1)/2)+1, 'Pixel where CRVAL1 is defined', sta)
        CALL ftpkyd(u,'CDELT1', r_max*2/pix, 1, 'Delta 1', sta)
        CALL ftpkys(u,'CUNIT1', TRIM(ref_u_str), 'Unit 1', sta)

        CALL ftpkys(u,'CTYPE2', 'LINEAR', 'Unit 2 type', sta)
        CALL ftpkyd(u,'CRVAL2', 0.0, 1, 'Value 1', sta)
        CALL ftpkyj(u,'CRPIX2', int((pix-1)/2)+1, 'Pixel where CRVAL2 is defined', sta)
        CALL ftpkyd(u,'CDELT2', r_max*2/pix, 1, 'Delta 2', sta)
        CALL ftpkys(u,'CUNIT2', TRIM(ref_u_str), 'Unit 2', sta)

    END SUBROUTINE add_essential_fits_keys
    
    SUBROUTINE sv_temp(basics, grid)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        !----------------------------------------------------------------------!
    
        CALL sv_temp_x(basics, grid)

    END SUBROUTINE sv_temp
    
    SUBROUTINE sv_temp_x(basics, grid)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Grid_TYP), INTENT(IN)                   :: grid
        !----------------------------------------------------------------------!
        INTEGER                                      ::  i_a, i_b, i_c
        INTEGER                                      ::  i_dust, i_cell
        REAL(kind=r2)                                ::  hd_r, hd_tem
        
        CHARACTER(len=252)                           :: filename
        CHARACTER(len=4)                             :: fileext
        CHARACTER(len=256)                           :: outname
        CHARACTER(len=6)                             :: num2string
        !----------------------------------------------------------------------!

        fileext = '.dat'
        DO i_dust = 1, grid%nh_n_dust
            write(num2string,'(I2.2)') i_dust
            fileext = '.dat'
            filename = TRIM(basics%path_results)//Getproname(basics)//'_temp_x_'
            outname = TRIM(filename)//TRIM(num2string)//fileext
            open(unit=1, file=TRIM(outname), &
                action="write", status="unknown", form="formatted")
            SELECT CASE(GetGridName(grid))
                CASE('spherical')

                    i_b = ceiling(real(grid%n(2))/2.0)
                    i_c = 1
                    
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
                    i_b = ceiling(real(grid%n(2))/2.0)
                    i_c = ceiling(real(grid%n(3))/2.0)

                    DO i_a = ceiling(real(grid%n(1))/2.0), grid%n(1)
                        hd_r   = (grid%co_mx_a(i_a-1) + grid%co_mx_a(i_a))/2.0_r2
                        i_cell = grid%cell_idx2nr(i_a, i_b, i_c)

                        IF (grid%grd_dust_density(i_cell, i_dust) .eq. 0.0) THEN
                            hd_tem = 0.0_r2
                        ELSE
                            hd_tem = grid%t_dust(i_cell, i_dust)
                        END IF
                       write(unit=1,fmt=*) hd_r, hd_tem
                    END DO

                CASE DEFAULT
                    print *, 'selected coordinate system not found, save temp_x'
            END SELECT
            close(unit=1)
        END DO
    
    END SUBROUTINE sv_temp_x
    
    SUBROUTINE save_input(basics,input_file)
    ! not used at the moment
    IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        CHARACTER(len=*)                             :: input_file
        CHARACTER(len=256)                           :: command

        !----------------------------------------------------------------------!
      
        command = 'cp '//TRIM(input_file)//' '//TRIM(basics%path_results)//    &
                        Getproname(basics)//'_input_file.dat'
        ! I'm not sure if the SYSTEM command is a good choice.
        CALL SYSTEM(command)
        
    END SUBROUTINE save_input
    
    SUBROUTINE save_continuum_map(model, basics, dust, fluxes, mode, peel_off)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                  :: model
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Dust_TYP), INTENT(IN)                   :: dust
        TYPE(Fluxes_TYP), INTENT(IN)                 :: fluxes
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)                            :: mode
        INTEGER                                        :: i_lam, sta, u
        INTEGER                                        :: i_stokes
        INTEGER                                        :: i_map
        REAL(kind=r2)                                  :: unit_value
        CHARACTER(len=100)                             :: outname
        CHARACTER(len=256)                             :: fits_name
        REAL(kind=r2),DIMENSION(:,:,:,:), ALLOCATABLE  :: map_out
        REAL(kind=r2),DIMENSION(1:dust%n_lam, 1:4)     :: sed
        LOGICAL, INTENT(IN)                            :: peel_off
        !----------------------------------------------------------------------!
        sed(:, :) = 0.0_r2
        ALLOCATE(map_out(0:2*model%n_bin_map,                                  &
                         0:2*model%n_bin_map, 1:dust%n_lam, 1:4))
        map_out = fluxes%continuum_map
        IF (mode .eq. 1) THEN
            ! finite wavelength/frequency bin
            outname = '_bin'

            DO i_map = 1, model%n_map
                map_out(:,:,1,:) = 0.0_r2
                DO i_lam = 2, dust%n_lam
                    IF (peel_off) THEN
                        unit_value =                                           &
                                 1.0_r2/dust%d_nu(i_lam) *                     &
                                 1.0e26_r2 /                                   &
                                 (4.0_r2 * PI *(model%distance*con_pc)**2)
                    ELSE
                        unit_value =                                           &
                                 1.0_r2/dust%d_nu(i_lam)/                      &
                                 ((1.0_r2-model%al_map(i_map))*(PI*2.0_r2)) *  &
                                 1.0e26_r2 /                                   &
                                 (model%distance*con_pc)**2
                    END IF
                    map_out(:, :, i_lam, :) =                                  &
                             map_out(:, :, i_lam, :) * unit_value
                    DO i_stokes = 1, 4
                        sed(i_lam, i_stokes) = sum(map_out(:,:,i_lam, i_stokes))
                    END DO
                END DO
            END DO
        ELSEIF  (mode .eq. 2) THEN
            ! monochromatic RT results
            outname = '_mono'
            DO i_map = 1, model%n_map
                DO i_lam = 1, dust%n_lam
                    IF (peel_off) THEN
                        unit_value =  dust%lam(i_lam)**2/con_c/                &
                                 (4.0_r2 * PI *(model%distance*con_pc)**2)     &
                                 *1.0e26_r2
                    ELSE
                        unit_value =  dust%lam(i_lam)**2/con_c/                &
                                 ((1.0_r2-model%al_map(i_map))*(PI*2.0_r2)) /  &
                                 ((model%distance*con_pc)**2) *                &
                                 1.0e26_r2
                    END IF
                    map_out(:, :, i_lam, :) =                                  &
                             map_out(:, :, i_lam, :) * unit_value
                    DO i_stokes = 1, 4
                        sed(i_lam, i_stokes) = sum(map_out(:,:,i_lam, i_stokes))
                    END DO
                END DO
            END DO
        ELSEIF  (mode .eq. 3) THEN
            outname = '_raytrace'
            DO i_lam = 1, dust%n_lam
                DO i_stokes = 1, 4
                    sed(i_lam, i_stokes) = sum(map_out(:, :, i_lam, i_stokes))
                END DO
            END DO
            ! no unit conversion here
        ELSE
            PRINT *, 'ERROR, wrong mode to write continuum map'
            
            DEALLOCATE(map_out)
            RETURN
        END IF
        u = 0
        sta = 0
        ! get a new u(nit) number
        CALL ftgiou(u,sta)
        ! init fits file
        fits_name = TRIM(basics%path_results)//Getproname(basics)//          &
                    '_continuum_map'//TRIM(outname)//'.fits.gz'
        CALL ftinit(u,'!'//fits_name,1,sta)
        ! write header
        CALL ftphpr(u,.true.,-64,4,(/2*model%n_bin_map+1,                  &
                                    2*model%n_bin_map+1,                   &
                                    dust%n_lam, 4/),                       &
                    0,1,.true.,sta)
        ! write array to fits file
        CALL ftpprd(u,1,1,(2*model%n_bin_map+1)**2*dust%n_lam*4, map_out(:,:,:,:), sta)
        
        ! add other header keywords 
        CALL add_essential_fits_keys(u, dust%n_dust,                           &
                                     2*model%n_bin_map+1, model%r_ou,          &
                                     GetModelName(model))
        ! close the fits file
        CALL ftclos(u, sta)
        ! free the (u)nit number
        CALL ftfiou(u, sta)
        DEALLOCATE(map_out)
        !write sed in ascii format
        DO i_stokes = 1, 4
            OPEN(unit=2, file=TRIM(basics%path_results)// Getproname(basics)// &
                              '_continuum_sed' // TRIM(outname) // '_' //      &
                              TRIM(fluxes%stokes_ext(i_stokes)) // '.dat',     &
                action="write", status="unknown", form="formatted")
            WRITE(unit=2,fmt='(A)') '#wavelength [m]    flux [Jy]'
            WRITE(unit=2,fmt=*) ''
            DO i_lam =1, dust%n_lam
                WRITE (unit=2,fmt='(2(ES15.6E3))') dust%lam(i_lam),            &
                                                   sed(i_lam, i_stokes)
            END DO
            CLOSE(unit=2)
        END DO

    END SUBROUTINE save_continuum_map

    SUBROUTINE save_ch_map(model, basics, gas, fluxes, n_dust)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                  :: model
        TYPE(Basic_TYP), INTENT(IN)                  :: basics
        TYPE(Gas_TYP), INTENT(IN)                    :: gas
        TYPE(Fluxes_TYP), INTENT(IN)                 :: fluxes
        INTEGER, INTENT(IN)                          :: n_dust
        !----------------------------------------------------------------------!
        INTEGER                                      :: vch, u, sta
        !----------------------------------------------------------------------!
        
        ! write velocity channel map to fits file
        sta = 0
        
        ! get a new u(nit) number
        call ftgiou(u,sta)
        ! init fits file
        call ftinit(u,'!'//TRIM(basics%path_results)//Getproname(basics)//'_velo_ch_map.fits.gz',1,sta)
        ! write header
        call ftphpr(u,.true.,-32,3,(/2*model%n_bin_map+1,2*model%n_bin_map+1,gas%i_vel_chan*2+1/),0,1,.true.,sta)
        ! write array to fits file
        call ftpprd(u,1,1,(2*model%n_bin_map+1)**2*(gas%i_vel_chan*2+1),fluxes%channel_map(:,:,:,1),sta)
        ! add other header keywords 
        CALL add_essential_fits_keys(u, n_dust,                                &
                                     2*model%n_bin_map+1, model%r_ou,          &
                                     GetModelName(model))

        ! close the fits file
        call ftclos(u, sta)
        ! free the (u)nit number
        call ftfiou(u, sta)
        
        ! write velocity integrated map to fits file
        sta = 0
        
        ! get a new u(nit) number
        call ftgiou(u,sta)
        ! init fits file
        call ftinit(u,'!'//TRIM(basics%path_results)//Getproname(basics)//'_velo_int_map.fits.gz',1,sta)
        ! write header
        call ftphpr(u,.true.,-64,2,(/2*model%n_bin_map+1,2*model%n_bin_map+1/),0,1,.true.,sta)
        ! write array to fits file
        call ftpprd(u,1,1,(2*model%n_bin_map+1)**2,sum(fluxes%channel_map(:,:,:,1),DIM=3)/ &
                    real(gas%i_vel_chan, kind=r2)*gas%vel_max*1e-3,sta)
        ! add other header keywords
!~         call ftpkyj(u,'EXPOSURE',1500,'Total Exposure Time',sta)
        ! close the fits file
        call ftclos(u, sta)
        ! free the (u)nit number
        call ftfiou(u, sta)
        
        
        !write full spectrum in ascii format 
        OPEN(unit=2, file=TRIM(basics%path_results)//Getproname(basics)//'_spectrum.dat', &
            action="write", status="unknown", form="formatted")
        
        ! write file header


        DO vch=-gas%i_vel_chan, gas%i_vel_chan
        
            WRITE(unit=2,fmt='(2(ES15.6E3))')   gas%velo_channel(vch), &
                                                sum(fluxes%channel_map(:,:,vch,1))
        END DO
        
        CLOSE(unit=2)
        
    END SUBROUTINE save_ch_map

END MODULE fileio
