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
!  grid specific routines
! 
! ---
MODULE grd_mod
    
    USE datatype
    USE var_global

    USE basic_type
    USE dust_type
    USE grid_type
    USE gas_type
    USE model_type

    USE math_mod
    USE fileio
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC ::   make_grid
    !--------------------------------------------------------------------------!

CONTAINS

  ! ############################################################################
  ! make grid
  ! ---
    SUBROUTINE make_grid(basics, grid, model, dust, gas)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        TYPE(Dust_TYP), INTENT(IN)                  :: dust
        TYPE(Gas_TYP), INTENT(INOUT)                :: gas
        
        !----------------------------------------------------------------------!

        INTEGER                                     :: i_a
        INTEGER                                     :: i_b
        INTEGER                                     :: i_c
        INTEGER                                     :: i_cell
        INTEGER                                     :: i_dust
        
        REAL(kind=r2)                               :: hd_totalmass
        REAL(kind=r2)                               :: hd_dusttotalmass
        REAL(kind=r1)                               :: hd1
        REAL(kind=r2)                               :: hd2

        REAL(kind=r2), DIMENSION(:), ALLOCATABLE    :: hd_arr1
        REAL(kind=r2), DIMENSION(1:3)               :: moco
        !----------------------------------------------------------------------!

        ! ---
        ! 1. set boundaries of each cell for the selected coordinate system

        Call set_boundaries(grid, model,basics)

        ! ---
        ! 2. - set cell numbering
        !    - calculate cell volumes
        !    - estimate number of grain / cell
        print *, "preparing grid [this may take a while]"

        ! dust-free region around center (r < r_in):
        i_cell      = 0
    
        ! inside model: grid (r >= r_in):
        do i_a=1, grid%n(1)
            do i_b=1, grid%n(2)
                do i_c=1, grid%n(3)
                    i_cell = i_cell + 1

                    ! a/b/c index => cell number
                    grid%cell_idx2nr(i_a, i_b, i_c) = i_cell 

                    ! cell number => a/b/c index
                    grid%cell_nr2idx(1, i_cell) = i_a
                    grid%cell_nr2idx(2, i_cell) = i_b
                    grid%cell_nr2idx(3, i_cell) = i_c

                    ! find cell volume and min A cell-Wall
                    CALL calc_cell_properties(grid,model,i_a,i_b,i_c)
                    ! calculate midpoint coordinate
                    moco(1) = (grid%co_mx_a(i_a) + grid%co_mx_a(i_a - 1)) / 2.0_r2
                    moco(2) = (grid%co_mx_b(i_b) + grid%co_mx_b(i_b - 1)) / 2.0_r2
                    moco(3) = (grid%co_mx_c(i_c) + grid%co_mx_c(i_c - 1)) / 2.0_r2

                    ! get cartesian coordinate from moco
                    grid%cellmidcaco(i_cell,:) = mo2ca(grid, moco)
                end do
            end do
        end do

        ! define ghost cells:
        ! id = 0 : inner null (spherical, cylindrical)
        ! id = n_cell + 1 : outer shell of the model
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            grid%cell_idx2nr(grid%n(1)+1, 1:, 1:) = grid%n_cell + 1
            grid%cell_idx2nr(1:, grid%n(2)+1, 1:) = grid%n_cell + 1 ! should not appear
            grid%cell_idx2nr(1:, 1:, grid%n(3)+1) = grid%n_cell + 1 ! should not appear
        CASE('cylindrical')
            grid%cell_idx2nr(grid%n(1)+1, 1:, 1:) = grid%n_cell + 1
            grid%cell_idx2nr(1:, grid%n(2)+1, 1:) = grid%n_cell + 1 ! should not appear
            grid%cell_idx2nr(:, :, grid%n(3)+1) = grid%n_cell + 1
            grid%cell_idx2nr(:, :, 0) = grid%n_cell + 1
            
        CASE('cartesian')
            grid%cell_idx2nr(grid%n(1)+1, :, :) = grid%n_cell + 1
            grid%cell_idx2nr(0, :, :) = grid%n_cell + 1
            grid%cell_idx2nr(:, grid%n(2)+1, :) = grid%n_cell + 1
            grid%cell_idx2nr(:, 0, :) = grid%n_cell + 1
            grid%cell_idx2nr(:, :, grid%n(3)+1) = grid%n_cell + 1
            grid%cell_idx2nr(:, :, 0) = grid%n_cell + 1
            
        CASE DEFAULT
            print *, 'selected coordinate system not found, setting ghost cells'
            stop
        END SELECT
        !TbD volume and area of zero cell

        grid%cell_minA(0) = minval(grid%cell_minA(1:grid%n_cell))


        ! fill grid with disk properties,e.g. temp, density, velocity...
        CALL set_grid_properties(basics, grid, dust, gas, model)


        ! verification: maximum cell number (i.e., total number of cells) = total number of cells
        if ( i_cell /= grid%n_cell ) then
            print *,"subroutine mk_grd(): wrong number of ESCs: ", i_cell
            print *,"should be                                : ", grid%n_cell 
            stop
        end if

        ! verification
        hd1 = sum(grid%cell_vol(:) / model%ref_unit**3)
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            hd2 = (4.0_r2/3.0_r2) * PI * (model%r_ou**3 - model%r_in**3)
        CASE('cylindrical')
            hd2 = 2.0_r2 * model%r_ou * PI *(model%r_ou**2 - model%r_in**2)
        CASE('cartesian')
            hd2 = (2.0_r2 * model%r_ou)**3
        CASE DEFAULT
            print *, 'selected coordinate system not found, verification of volume space'
            stop
        END SELECT

        if ( abs(hd1-hd2)/hd1 > 0.01_r2 ) then
            print *,"ERROR"
            print *,"    Difference between the VOLUME OF THE MODEL SPACE"
            print *,"    and the TOTAL VOLUME OF THE SINGLE ESCs is larger than 1% :"
            print *,"    - volume( model space)        : ", hd2, " AU^3"
            print *,"    - total volume of single ESCs : ", hd1, " AU^3"
            print *
            print *,"Please verify that you used a correct (odd/even) number of grid cells"
            stop
        end if

        print *, '                                                             '
        !----------------------------------------------------------------------!
        !

        IF (.not. basics%old_model) THEN

            hd_totalmass = 0.0_r2
            DO i_dust=1, dust%n_dust
                
                hd_totalmass =  hd_totalmass + sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 * &
                                dust%r_dust(i_dust)**3 * dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
            END DO
            ! rescale mass to user-defined value
            grid%grd_dust_density(:,:) = grid%grd_dust_density(:,:) * model%mass/hd_totalmass
            grid%Nv(:,:)               = grid%Nv(:,:)               * model%mass/hd_totalmass

        END IF

        ! dust mass verification
        hd_totalmass = 0.0_r2
        hd_dusttotalmass = 0.0_r2
        do i_dust=1, dust%n_dust
            hd_dusttotalmass =  hd_dusttotalmass + sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 * &
                            dust%r_dust(i_dust)**3 * dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
        end do
        
        ! estimate gas mass and rescale 
        ! Nv_mol,  defined molecule
        ! Nv_col,  collision partners 1 = H , 2 = H para ...,3 = H ortho
        !---- H2
        ! total gas mass
        hd_totalmass =  sum(grid%Nv_col(:,1)) * col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun + &
                        sum(grid%Nv_mol(:)) * gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
                 
        
        ! rescale total mass
        grid%grd_col_density(:,1) = grid%grd_col_density(:,1) * hd_dusttotalmass*gas%dust_mol_ratio/hd_totalmass
        
        grid%Nv_col(:,1) = grid%grd_col_density(:,1) * REAL(grid%cell_vol(:),kind=r2)

        !---- H2 para

        grid%Nv_col(:,2) =  0.25*grid%Nv_col(:,1) 
                            
        grid%grd_col_density(:,2) = grid%Nv_col(:,2)/grid%cell_vol(:)
        
        !---- H2 ortho

        grid%Nv_col(:,3) =  0.75*grid%Nv_col(:,1) 

        grid%grd_col_density(:,3) = grid%Nv_col(:,3)/grid%cell_vol(:)
        

        !---- selected molecule
        
        grid%grd_mol_density(:) = grid%grd_mol_density(:) * hd_dusttotalmass*gas%dust_mol_ratio/hd_totalmass
        grid%Nv_mol(:) = grid%grd_mol_density(:) * REAL(grid%cell_vol(:),kind=r2)

        !TbD: all other partners:
        !grid%Nv_col(:,4:6)            = 0.0_r2
        !grid%grd_col_density(:,4:6)   = 0.0_r2
        !grid%grd_col_density(0,:)     = 0.0_r2
        
        hd_totalmass = 0.0_r2
        hd_dusttotalmass = 0.0_r2
        ! calculate total mass

        ! dust mass
        DO i_dust=1, dust%n_dust
            hd_dusttotalmass =  hd_dusttotalmass +                             &
                sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 *                    &
                dust%r_dust(i_dust)**3 *                                       &
                dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
        END DO

        print '(A,ES11.4,A)',' dust mass           : '                         &
               ,hd_dusttotalmass, ' M_sun'
        IF (dust%n_dust > 1) THEN
            DO i_dust=1, dust%n_dust
            print '(A,I2,A,ES11.4,A)','          species ', i_dust, ' : ',     &
                sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 *                    &
                dust%r_dust(i_dust)**3 *                                       &
                dust%den_dust(i_dust)*1.0e+3_r2 / M_sun, ' M_sun'
            END DO
        END IF

        hd_totalmass = hd_dusttotalmass +                                      &
            ! H mass
            sum(grid%Nv_col(:,1)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun +   &
            ! selected molecule mass
            sum(grid%Nv_mol(:))     *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
        print *, ''
        print '(A,ES11.4,A)',' H2 mass             : '                               &
              ,sum(grid%Nv_col(:,1)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
        print '(A,ES11.4,A)','          orthogonal : '&
              ,sum(grid%Nv_col(:,3)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
              
        print '(A,ES11.4,A)','            parallel : '&
              ,sum(grid%Nv_col(:,2)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
              
        print *, ''
        print '(A,ES11.4,A)',' molecule mass       : '&
              ,sum(grid%Nv_mol(:)) *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun, ' M_sun'

        print '(A)', ' ---------------------------------------'
        print '(A,ES11.4,A)', ' total disk mass     : ', hd_totalmass, ' M_sun'
        print *,''
        !----------------------------------------------------------------------!
        ! ---
        ! set smallest step width for photon transfer
        ! - goal: avoid step widths smaller than the amount that will change
        !         a r2 type floating point
        !         variable (in addition operation)
        ! - factor: 10.0_r2: margin (for safety)
        grid%d_l_min = model%r_ou * epsilon(1.0_r2) * 10.0_r2

        ! check if d_l_min is significantly larger (factor 10) than smallest grid cell
        allocate( hd_arr1(1:grid%n(1)) )
        do i_a=1, grid%n(1)
            hd_arr1(i_a) = grid%co_mx_a(i_a) - grid%co_mx_a(i_a-1)
        end do

        if (grid%d_l_min > 10.0_r2*minval(hd_arr1(:))) then
            print *, "error   : mk_grd: radial extent of smallest cell is too small."
            print *, "solution: (a) decrease step width factor [sf] and/or"
            print *, "          (b) reduce number of grid cells"
        else
            print '(A,ES11.4,A,A)'," smallest step width : ", grid%d_l_min,    &
                  ' ', GetModelName(model)
        end if

        deallocate( hd_arr1 )

    END SUBROUTINE make_grid

    ! #########################################################################

    SUBROUTINE calc_cell_properties(grid, model, i_a, i_b, i_c)
        !generic routine
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)              :: grid
        TYPE(Model_TYP), INTENT(IN)                :: model
        !----------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_a, i_b, i_c
        !----------------------------------------------------------------------!
        SELECT CASE(GetGridName(grid))
        CASE('spherical')
            CALL calc_cell_properties_sp(grid, model, i_a, i_b, i_c)

        CASE('cylindrical')
            CALL calc_cell_properties_cy(grid, model, i_a, i_b, i_c)

        CASE('cartesian')
            CALL calc_cell_properties_ca(grid, model, i_a, i_b, i_c)

        CASE DEFAULT
            print *, 'selected coordinate system not found, calc_cell_properties'
            stop
        END SELECT


    END SUBROUTINE calc_cell_properties

    SUBROUTINE calc_cell_properties_ca(grid, model, i_x, i_y, i_z)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)              :: grid
        TYPE(Model_TYP), INTENT(IN)                :: model
        !----------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_x, i_y, i_z
        INTEGER                                    :: i_cell
        REAL(kind=r2)                              :: dx
        REAL(kind=r2)                              :: dy
        REAL(kind=r2)                              :: dz
        !----------------------------------------------------------------------!
        
        i_cell = grid%cell_idx2nr(i_x, i_y, i_z)
        dx = abs(grid%co_mx_a(i_x-1) - grid%co_mx_a(i_x))
        dy = abs(grid%co_mx_b(i_y-1) - grid%co_mx_b(i_y))
        dz = abs(grid%co_mx_c(i_z-1) - grid%co_mx_c(i_z))

        grid%cell_minA(i_cell) = MIN( dx * dy, dx * dz, dz * dy)

        ! volume of individual cells [m^3]
        
        grid%cell_vol(i_cell) = dx * dy * dz * model%ref_unit**3

    END SUBROUTINE calc_cell_properties_ca
    
    
    SUBROUTINE calc_cell_properties_cy(grid, model, i_r, i_ph, i_z)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)              :: grid
        TYPE(Model_TYP), INTENT(IN)                :: model
        !----------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_r, i_ph, i_z
        INTEGER                                    :: i_cell
        REAL(kind=r2)                              :: dth
        REAL(kind=r2)                              :: dz
        !----------------------------------------------------------------------!
        
        i_cell = grid%cell_idx2nr(i_r,i_ph,i_z)
        dz = abs(grid%co_mx_c(i_z-1) - grid%co_mx_c(i_z))
        dth = 2.0_r2*PI / real(grid%n(2), kind=r2)
        
        grid%cell_minA(i_cell) = MIN( 0.5* dth *  (grid%co_mx_a(i_r)**2-grid%co_mx_a(i_r-1)**2), &
                                    ( dth * grid%co_mx_a(i_r-1) * dz), &
                                    ( (grid%co_mx_a(i_r)-grid%co_mx_a(i_r-1))*dz))
                                    
        ! volume of individual cells [m^3]
        
        grid%cell_vol(i_cell) = 0.5_r2 * dz * dth * &
                               (grid%co_mx_a(i_r)**2-grid%co_mx_a(i_r-1)**2) * &
                                model%ref_unit**3

    END SUBROUTINE calc_cell_properties_cy

    SUBROUTINE calc_cell_properties_sp(grid,model, i_r, i_th, i_ph)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !----------------------------------------------------------------------!
        INTEGER,INTENT(IN)                          :: i_r, i_th, i_ph
        INTEGER                                     :: i_cell
        !----------------------------------------------------------------------!
        
        i_cell = grid%cell_idx2nr(i_r, i_th , i_ph)

        grid%cell_minA(i_cell) = MIN( &
                                grid%co_mx_a(i_r-1)**2.0_r2 * &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1) ) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1)) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 )* &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) )

        ! volume of individual cells [m^3]
        grid%cell_minA(i_cell) = grid%cell_minA(i_cell)
        grid%cell_vol(i_cell) = &
                                abs( &
                                (2.0_r2*PI/3.0_r2) * &
                                (grid%co_mx_a(i_r)**3 - grid%co_mx_a(i_r-1)**3)   * &
                                (sin(grid%co_mx_b(i_th)) - sin(grid%co_mx_b(i_th-1))) * &
                                ((grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1))/(2.0_r2*PI)) &
                                ) * model%ref_unit**3

    END SUBROUTINE calc_cell_properties_sp


    !  ! #######################################################################
    SUBROUTINE set_boundaries(grid,model,basics)
        !generic routine
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        CHARACTER(len=252)                          :: file_a, file_b, file_c
        !----------------------------------------------------------------------!
        IF (basics%old_model) THEN
            ! if we use an old model, we simply read the boundaries from the input files
            file_a = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_a_boundaries.dat'
            file_b = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_b_boundaries.dat'
            file_c = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_c_boundaries.dat'
            CALL read_boundaries(grid, file_a, file_b, file_c)
        
        ELSE
            SELECT CASE(GetGridName(grid))
            
            CASE('spherical')
                CALL set_boundaries_sp(grid, model)

            CASE('cylindrical')
                CALL set_boundaries_cy(grid, model)

            CASE('cartesian')
                CALL set_boundaries_ca(grid, model)

            CASE DEFAULT
                print *, 'selected coordinate system not found, set_boundaries'
                stop
            END SELECT
        END IF
        
        CALL save_boundaries(grid,basics)
    
    END SUBROUTINE set_boundaries
    
    SUBROUTINE set_boundaries_cy(grid,model)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !----------------------------------------------------------------------!
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_z
        
        REAL(kind=r2)                               :: dr1
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dz
        REAL(kind=r2)                               :: sf
        !----------------------------------------------------------------------!
        
        SELECT CASE(GetGridType(grid))
        CASE(1)
            ! normal grid
            
            ! ---
            ! 1.1. rho
            ! co_mx_a( 0 ): r_in
            grid%co_mx_a(0) = model%r_in

            ! dr1: radial extension of first cell
            dr1 = (model%r_ou - model%r_in) * (grid%sf-1.0_r2)/ (grid%sf**grid%n(1) - 1.0_r2)

            ! set outer ring radii:
            ! co_mx_a( 1 ): r_in + dr1
            ! co_mx_a(n(1)): r_ou
            do i_r=1, grid%n(1)
                grid%co_mx_a(i_r) = grid%co_mx_a(0) + dr1 * (grid%sf**i_r - 1.0_r2) / (grid%sf-1.0_r2)
            end do
            
            ! ---
            ! 1.2 theta
            ! co_mx_b(  0 ) = 0    (  0°)
            ! co_mx_b(n(3)) = 2xPI (360°)
            
            grid%co_mx_b(0) = 0.0_r2
            !grid%co_mx_b(0) = -PI
        
            dth = 2.0_r2*PI / real(grid%n(2), kind=r2)
            do i_th=1, grid%n(2)
                grid%co_mx_b(i_th) =  grid%co_mx_b(0) + dth * real(i_th, kind=r2)
            end do
            
            ! ---
            ! 1.3. z
            ! co_mx_c( 0 )  =- r_out
            ! co_mx_c(n(3)) =  r_out
            ! dz should not be constant, because we want to have small cells in the inner region
            !     -> we use sf too, but not the user given value but a fixed one
            
    !~         grid%co_mx_c(grid%n(3)*0.5) = 0.0_r2
    !~         sf = 1.08_r2
            sf = grid%sf+0.01
            !dz = 2.0_r2*model%r_ou / real(grid%n(3), kind=r2)
            dz = model%r_ou * (sf-1.0_r2)/ (sf**((grid%n(3)-1)*0.5) - 1.0_r2)
            do i_z=1, int(grid%n(3)*0.5)
                grid%co_mx_c(int((grid%n(3)-1)*0.5+i_z+1)) =            dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
                grid%co_mx_c(int((grid%n(3)-1)*0.5-i_z)) =  -1.0_r2 * dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
            end do
            grid%co_mx_c(int((grid%n(3)-1)*0.5)) = -1.0_r2/3.0_r2 * dz 
            grid%co_mx_c(int((grid%n(3)-1)*0.5+1)) = 1.0_r2/3.0_r2 * dz 

        CASE(2)
            !grid based on Tobias self similar disk simulations
            ! ---
            ! 1.1. rho
            ! co_mx_a( 0 ): r_in
            grid%co_mx_a(0) = model%r_in
            DO i_r = 1,grid%n(1)
                grid%co_mx_a(i_r) = model%r_in*(model%r_ou/model%r_in)**(float(i_r)/grid%n(1))
            END DO
            ! ---
            ! 1.2 theta
            ! co_mx_b(  0 ) = 0    (  0°)
            ! co_mx_b(n(3)) = 2xPI (360°)
            
            grid%co_mx_b(0) = 0.0_r2
            !grid%co_mx_b(0) = -PI
        
            dth = 2.0_r2*PI / real(grid%n(2), kind=r2)
            DO i_th=1, grid%n(2)
                grid%co_mx_b(i_th) =  grid%co_mx_b(0) + dth * real(i_th, kind=r2)
            END DO
            
            ! ---
            ! 1.3. z
            ! co_mx_c( 0 )  =- r_out
            ! co_mx_c(n(3)) =  r_out
            grid%co_mx_c(0)         = -model%r_ou
            DO i_th = 0, grid%n(3)-2
                grid%co_mx_c(i_th+1) = 520.0*sinh(5.0*(2.0*(float(i_th)/(grid%n(3)-2))-1))/sinh(5.0)
            END DO
            grid%co_mx_c(grid%n(3)) = model%r_ou

        CASE(9)
           ! ---
            ! This case is for a very general input of a spherical grid,
            ! coordinates defined in external files
            ! The style of the input file has to be consistent
            
            CALL read_boundaries(grid)

            ! ---
            ! 1 r
            ! test for consistency
            ! co_mx_a(  0 ) = r_in
            ! co_mx_a(n(1)) = r_ou
            
            IF ( abs(model%r_in - grid%co_mx_a(0)) .gt. 1.0e3_r2*epsilon(model%r_in) ) THEN
                PRINT *, "ERROR: Inner rim given in input file is not consistent with the grid"
                PRINT *, "r_in               : ", model%r_in
                PRINT *, "innerst cell bound : ", grid%co_mx_a(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_a(grid%n(1))) .gt. 1.0e3_r2*epsilon(model%r_ou) ) THEN
                PRINT *, "ERROR: Outer rim given in input file is not consistent with the grid"
                PRINT *, "r_ou             : ", model%r_ou
                PRINT *, "outer cell bound : ", grid%co_mx_a(grid%n(1))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_a(0) = model%r_in
            grid%co_mx_a(grid%n(1)) = model%r_ou

            ! ---
            ! 2 phi
            ! test for consistency
            ! co_mx_b(  0 ) = 0         (   0°)
            ! co_mx_b(n(2)) = 2 * PI    ( 360°)

            IF ( abs(grid%co_mx_b(0)) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: lower phi coordinate is not 0.0"
                PRINT *, "lowest cell boundary : ", grid%co_mx_b(0)
                STOP
            ELSEIF ( abs(2.0 * PI - grid%co_mx_b(grid%n(2))) .gt. 1.0e4_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: upper phi coordinate is not 2 * PI"
                PRINT *, "upper cell boundary : ", grid%co_mx_b(grid%n(2))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_b(0) = 0.0_r2
            grid%co_mx_b(grid%n(2)) = 2_r2 * PI

            ! ---
            ! 3 z
            ! test for consistency
            ! co_mx_c(  0 ) = -r_ou
            ! co_mx_c(n(3)) = r_ou

            IF ( abs(model%r_ou + grid%co_mx_c(0)) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: lower z coordinate is not -r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_c(grid%n(3))) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: upper z coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(grid%n(3))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_c(0) = -model%r_ou
            grid%co_mx_c(grid%n(3)) = model%r_ou
        END SELECT

    END SUBROUTINE set_boundaries_cy
    
    SUBROUTINE set_boundaries_ca(grid, model)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !----------------------------------------------------------------------!
        INTEGER                                     :: i_x
        INTEGER                                     :: i_y
        INTEGER                                     :: i_z
        
        REAL(kind=r2)                               :: dx
        REAL(kind=r2)                               :: dy
        REAL(kind=r2)                               :: dz
        REAL(kind=r2)                               :: sf
        !----------------------------------------------------------------------!
        
        SELECT CASE(GetGridType(grid)) 
        CASE(1)
            ! normal grid (TbD, use a sinh spacing)
            
            ! ---
            ! 1.1. x
            ! co_mx_a( 0 )  =- r_out
            ! co_mx_a(n(1)) =  r_out
            ! dz should not be constant, because we want to have small cells in
            !        the inner region 
            !     -> we use sf too, but not the user given value but a fixed one

            sf = grid%sf+0.01
            dx = model%r_ou * (sf-1.0_r2)/ (sf**((grid%n(1)-1)*0.5) - 1.0_r2)
            do i_x=1, int(grid%n(1)*0.5)
                grid%co_mx_a(int((grid%n(1)-1)*0.5+i_x+1)) = dx * (sf**i_x - 1.0_r2) / (sf-1.0_r2)
                grid%co_mx_a(int((grid%n(1)-1)*0.5-i_x)) = -1.0_r2 * dx * (sf**i_x - 1.0_r2) / (sf-1.0_r2)
            end do
            grid%co_mx_a(int((grid%n(1)-1)*0.5)) = -1.0_r2/3.0_r2 * dx 
            grid%co_mx_a(int((grid%n(1)-1)*0.5+1)) = 1.0_r2/3.0_r2 * dx
            
            ! ---
            ! 1.2 y
            ! co_mx_b( 0 )  =- r_out
            ! co_mx_b(n(2)) =  r_out
            ! dz should not be constant, because we want to have small cells in
            !        the inner region (TbD, use a sinh spacing)
            !     -> we use sf too, but not the user given value but a fixed one

            sf = grid%sf+0.01
            dy = model%r_ou * (sf-1.0_r2)/ (sf**((grid%n(2)-1)*0.5) - 1.0_r2)
            do i_y=1, int(grid%n(2)*0.5)
                grid%co_mx_b(int((grid%n(2)-1)*0.5+i_y+1)) = dy * (sf**i_y - 1.0_r2) / (sf-1.0_r2)
                grid%co_mx_b(int((grid%n(2)-1)*0.5-i_y)) = -1.0_r2 * dy * (sf**i_y - 1.0_r2) / (sf-1.0_r2)
            end do
            grid%co_mx_b(int((grid%n(2)-1)*0.5)) = -1.0_r2/3.0_r2 * dy 
            grid%co_mx_b(int((grid%n(2)-1)*0.5+1)) = 1.0_r2/3.0_r2 * dy 


            ! ---
            ! 1.3. z
            ! co_mx_c( 0 )  =- r_out
            ! co_mx_c(n(3)) =  r_out
            ! dz should not be constant, because we want to have small cells in
            !        the inner region (TbD, use a sinh spacing)
            !     -> we use sf too, but not the user given value but a fixed one

            sf = grid%sf+0.01
            dz = model%r_ou * (sf-1.0_r2)/ (sf**((grid%n(3)-1)*0.5) - 1.0_r2)
            do i_z=1, int(grid%n(3)*0.5)
                grid%co_mx_c(int((grid%n(3)-1)*0.5+i_z+1)) = dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
                grid%co_mx_c(int((grid%n(3)-1)*0.5-i_z)) = -1.0_r2 * dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
            end do
            grid%co_mx_c(int((grid%n(3)-1)*0.5)) = -1.0_r2/3.0_r2 * dz 
            grid%co_mx_c(int((grid%n(3)-1)*0.5+1)) = 1.0_r2/3.0_r2 * dz 

        CASE(9)
           ! ---
            ! This case is for a very general input of a spherical grid,
            ! coordinates defined in external files
            ! The style of the input file has to be consistent
            
            CALL read_boundaries(grid)

            ! ---
            ! 1 x
            ! test for consistency
            ! co_mx_a(  0 ) = -r_ou
            ! co_mx_a(n(1)) = r_ou

            IF ( abs(model%r_ou + grid%co_mx_a(0)) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: lower x coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_a(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_a(grid%n(1))) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: upper x coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_a(grid%n(1))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_a(0) = -model%r_ou
            grid%co_mx_a(grid%n(1)) = model%r_ou

            ! ---
            ! 2 y
            ! test for consistency
            ! co_mx_b(  0 ) = -r_ou
            ! co_mx_b(n(2)) = r_ou
            
            IF ( abs(model%r_ou + grid%co_mx_b(0)) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: lower y coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_b(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_b(grid%n(2))) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: upper y coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(grid%n(2))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_b(0) = -model%r_ou
            grid%co_mx_b(grid%n(2)) = model%r_ou

            ! ---
            ! 3 z
            ! test for consistency
            ! co_mx_c(  0 ) = -r_ou
            ! co_mx_c(n(3)) = r_ou

            IF ( abs(model%r_ou + grid%co_mx_c(0)) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: lower z coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_c(grid%n(3))) .gt. 1.0e3_r2*epsilon(model%r_ou)) THEN
                PRINT *, "ERROR: upper z coordinate is not r_ou"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(grid%n(3))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_c(0) = -model%r_ou
            grid%co_mx_c(grid%n(3)) = model%r_ou
            
        END SELECT

    END SUBROUTINE set_boundaries_ca

    SUBROUTINE set_boundaries_sp(grid,model)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        !----------------------------------------------------------------------!
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_ph
        
        REAL(kind=r2)                               :: dr1
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dph
        REAL(kind=r2)                               :: value_in
        !----------------------------------------------------------------------!
        SELECT CASE(GetGridType(grid))
        
        CASE(1)
            ! ---
            ! 1.1. rho
            ! co_mx_a( 0 ): r_in
            grid%co_mx_a(0) = model%r_in

            ! dr1: radial extension of first cell
            dr1 = (model%r_ou - model%r_in) * (grid%sf-1.0_r2)/                &
                  (grid%sf**grid%n(1) - 1.0_r2)

            ! set outer ring radii:
            ! co_mx_a( 1 ): r_in + dr1
            ! co_mx_a(n(1)): r_ou
            do i_r=1, grid%n(1)
                grid%co_mx_a(i_r) = grid%co_mx_a(0) + dr1 *                    &
                                    (grid%sf**i_r - 1.0_r2) / (grid%sf-1.0_r2)
            end do

            ! ---
            ! 1.2. theta
            ! co_mx_b(  0 ) = -PI/2 (-90°)
            ! co_mx_b(n(2)) = +PI/2 (+90°)
            grid%co_mx_b(0) = -PI/2.0_r2
            
            dth = PI / real(grid%n(2), kind=r2)
            do i_th=1, grid%n(2)
                grid%co_mx_b(i_th) = grid%co_mx_b(0) +  dth * real(i_th, kind=r2)
            end do
            
            ! this line prevents to create a cone with an open angle with exact
            ! 90 deg, which in fact is no cone anymore,
            ! maybe we ca find a better solution, but it is working
            IF (MOD(grid%n(2),2) == 0 ) THEN
                grid%co_mx_b(grid%n(2)/2) = grid%co_mx_b(grid%n(2)/2+1) * 1.0e-5_r2
            END IF
            ! ---
            ! 1.3 phi
            ! co_mx_c(  0 ) = 0    (  0°)
            ! co_mx_c(n(3)) = 2xPI (360°)
            grid%co_mx_c(0) = 0.0_r2
        
            dph = 2.0_r2*PI / real(grid%n(3), kind=r2)
            do i_ph=1, grid%n(3)
                grid%co_mx_c(i_ph) =  grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do
            
        CASE(2)
            ! here we have a linear spherical grid
            ! 
            !
            ! ---
            ! 1.1. rho
            ! co_mx_a( 0 ): r_in
            grid%co_mx_a(0) = model%r_in

            ! dr1: radial extension
            dr1 = (model%r_ou - model%r_in) / real(grid%n(1), kind=r2)

            ! set outer ring radii:
            ! co_mx_a( 1 ): r_in + dr1
            ! co_mx_a(n(1)): r_ou
            do i_r=1, grid%n(1)
                grid%co_mx_a(i_r) = grid%co_mx_a(0) + dr1 * real(i_r, kind=r2)
            end do
            
            ! ---
            ! 1.2. theta
            ! co_mx_b(  0 ) = -PI/2 (-90°)
            ! co_mx_b(n(2)) = +PI/2 (+90°)
            grid%co_mx_b(0) = -PI/2.0_r2

            dth = PI / real(grid%n(2), kind=r2)
            do i_th=1, grid%n(2)
                grid%co_mx_b(i_th) = grid%co_mx_b(0) + dth * real(i_th, kind=r2)
            end do
            ! this line prevents to create a cone with an open angle with exact
            ! 90 deg, which in fact is no cone anymore,
            ! maybe we ca find a better solution, but it is working
            IF (MOD(grid%n(2),2) == 0 ) THEN
                grid%co_mx_b(grid%n(2)/2) = grid%co_mx_b(grid%n(2)/2+1) * 1.0e-5_r2
            END IF
            ! ---
            ! 1.3 phi
            ! co_mx_c(  0 ) = 0    (  0°)
            ! co_mx_c(n(3)) = 2xPI (360°)
            grid%co_mx_c(0) = 0.0_r2
        
            dph = 2.0_r2*PI / real(grid%n(3), kind=r2)
            do i_ph=1, grid%n(3)
                grid%co_mx_c(i_ph) = grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do
            
        CASE(3)
            ! ---
            ! 1.1 r
            open(unit=1, file="input/grid/logr.dat", &
                    action="read", status="unknown", form="formatted")
                    
            DO i_r = 0, grid%n(1)
                read(unit=1,fmt=*) value_in
                grid%co_mx_a(i_r) = 10.0**(value_in)
            END DO
            close(unit=1)
            ! ---
            ! 1.2 th
            open(unit=1, file="input/grid/theta.dat", &
                    action="read", status="unknown", form="formatted")
                    
            DO i_th = 0, grid%n(2)
                read(unit=1,fmt=*) value_in
                grid%co_mx_b(i_th) = value_in
            END DO
            close(unit=1)
        
            ! ---
            ! 1.3 phi
            ! co_mx_c(  0 ) = 0    (  0°)
            ! co_mx_c(n(3)) = 2xPI (360°)
            grid%co_mx_c(0) = 0.0_r2
        
            dph = 2.0_r2*PI / real(grid%n(3), kind=r2)
            do i_ph=1, grid%n(3)
                grid%co_mx_c(i_ph) = grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do

        CASE(9)
            ! ---
            ! This case is for a very general input of a spherical grid,
            ! coordinates defined in external files
            ! The style of the input file has to be consistent
            
            CALL read_boundaries(grid)

            ! ---
            ! 1 r
            ! test for consistency
            ! co_mx_a(  0 ) = r_in
            ! co_mx_a(n(1)) = r_ou
            
            IF ( abs(model%r_in - grid%co_mx_a(0)) .gt. 1.0e3_r2*epsilon(model%r_in) ) THEN
                PRINT *, "ERROR: Inner rim given in input file is not consistent with the grid"
                PRINT *, "r_in               : ", model%r_in
                PRINT *, "innerst cell bound : ", grid%co_mx_a(0)
                STOP
            ELSEIF ( abs(model%r_ou - grid%co_mx_a(grid%n(1))) .gt. 1.0e3_r2*epsilon(model%r_ou) ) THEN
                PRINT *, "ERROR: Outer rim given in input file is not consistent with the grid"
                PRINT *, "r_ou             : ", model%r_ou
                PRINT *, "outer cell bound : ", grid%co_mx_a(grid%n(1))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_a(0) = model%r_in
            grid%co_mx_a(grid%n(1)) = model%r_ou
            ! ---
            ! 2 theta
            ! test for consistency
            ! co_mx_b(  0 ) = -PI/2    (- 90°)
            ! co_mx_b(n(2)) =  PI/2    (  90°)

            IF ( abs(PI/2_r2 + grid%co_mx_b(0)) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: lower theta coordinate is not -PI/2"
                PRINT *, "lowest cell boundary : ", grid%co_mx_b(0)
                STOP
            ELSEIF ( abs(PI/2_r2 - grid%co_mx_b(grid%n(2))) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: upper theta coordinate is not PI/2"
                PRINT *, "upper cell boundary : ", grid%co_mx_b(grid%n(2))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_b(0) = -PI/2_r2
            grid%co_mx_b(grid%n(2)) = PI/2_r2
            
            ! ---
            ! 3 phi
            ! test for consistency
            ! co_mx_c(  0 ) = 0         (   0°)
            ! co_mx_c(n(3)) = 2 * PI    ( 360°)

            IF ( abs(grid%co_mx_c(0)) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: lower phi coordinate is not 0.0"
                PRINT *, "lowest cell boundary : ", grid%co_mx_c(0)
                STOP
            ELSEIF ( abs(2.0 * PI - grid%co_mx_c(grid%n(3))) .gt. 1.0e4_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: upper phi coordinate is not 2 * PI"
                PRINT *, "upper cell boundary : ", grid%co_mx_c(grid%n(3))
                STOP
            END IF
            ! set the boundaries to the Mol3D values to ensure consistency
            grid%co_mx_c(0) = 0.0_r2
            grid%co_mx_c(grid%n(3)) = 2_r2 * PI
        END SELECT
    
    END SUBROUTINE set_boundaries_sp
    

    !  ! #######################################################################


    SUBROUTINE set_grid_properties(basics, grid, dust, gas, model)
    ! set grid density and temperature! (and velocity) 
    ! prepared for general coordinates
    
        USE model_mod
        USE parser_mod
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Dust_TYP), INTENT(IN)                  :: dust
        TYPE(GAS_TYP), INTENT(INOUT)                :: gas
        TYPE(Model_TYP), INTENT(IN)                 :: model
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        TYPE(Vector3d)                              :: p_vec
        !----------------------------------------------------------------------!
        INTEGER                                     :: i_cell, i, k
        INTEGER                                     :: i_r, i_th, i_ph
        INTEGER                                     :: io, N_dummy

        REAL(kind=r2), DIMENSION(3)                 :: moco, velo_hlp
        REAL(kind=r2), DIMENSION(3,3)               :: S
        REAL(kind=r2)                               :: value_in
        REAL(kind=r1)                               :: P_xy, P_z
        CHARACTER(len=256)                          :: waste
        !----------------------------------------------------------------------!
        ! ---

        !CALL parse('R_gap_in', R_gap_in,'input/additional.dat')
        !CALL parse('R_gap_ou', R_gap_ou,'input/additional.dat')

        IF (basics%old_model) THEN
            print '(2A)', ' | loading model parameter from: ', TRIM(basics%pronam_old)
            CALL read_model(grid, dust, TRIM(basics%path_results) //           &
                                  TRIM(basics%pronam_old)//'_model.fits')

            PRINT *, "| done!                        "
            CLOSE(unit=1)

        ELSE IF (GetGridType(grid) == 3 .and. GetGridName(grid) == 'spherical' ) THEN
            print *, 'Tobias density'
            OPEN(unit=1, file="input/grid/tobias_disk.dat", &
                action="read", status="unknown", form="formatted")
            DO i = 1,11
                !read header
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            moco(3) = PI
            DO i = 1, (grid%n(1))*(grid%n(2))*grid%n(3)
                READ(unit=1,fmt=*,iostat=io) moco(1),moco(2),value_in

                moco(1) = 10.0**(moco(1))
                i_cell = get_cell_nr(grid,mo2ca(grid,moco))
                grid%grd_dust_density(i_cell,:) = 10.0**(value_in)
                grid%grd_col_density(i_cell,:)  = 10.0**(value_in)
            END DO
            CLOSE(unit=1)
        ELSE IF (GetGridType(grid) == 6 .and. GetGridName(grid) == 'spherical' ) THEN
            ! note: 
            !       this one is just to read in Marios disk models
            !       (should be removed in near future)
            print *,'  read model from text_file'
            OPEN(unit=1, file="input/grid/HRFIX_410.dat", &
                         action="read", status="old", form="formatted")
            DO i = 1,4
            !read header
               READ(unit=1,fmt=*,iostat=io) waste
            END DO
            k = 1
            p_vec%c_sys = "spherical"
            DO i_r = 1, grid%n(1)
                DO i_th = 2, grid%n(2)-1
                    DO i_ph = 1, grid%n(3)
                        i_cell  = grid%cell_idx2nr(i_r,i_th,i_ph)
                        IF (i_cell >= int(k*grid%n_cell*0.01)) THEN
                            WRITE (*,'(A,I3,A)') ' | | | ',                    &
                                    int(i_cell/real(grid%n_cell)*100.0),       &
                                    '% done'//char(27)//'[A'
                            k = k +1
                        END IF
                        p_vec%comp = ca2sp(grid%cellmidcaco(i_cell,:))
                        READ(unit=1,fmt=*,iostat=io) N_dummy, velo_hlp, value_in
                        grid%grd_dust_density(i_cell,:) = value_in
                        grid%grd_col_density(i_cell,:)  = value_in
                        
                        CALL generate_rotation_matrix(p_vec,S)
                        grid%velo(i_cell,:)  =  matmul(S,velo_hlp)
                    END DO
                END DO
            END DO
            CLOSE(unit=1)
        ELSE IF (GetGridType(grid) == 9) THEN
            ! this is the most general way to read in the different 
            ! distribution (density, temperature, velocity...).
            ! It reads the "input/grid/model.fits" file that 
            ! should be a symbolic link to the actual model file.

            PRINT *, "Read model distribution from file"
            CALL read_model(grid, dust, "input/grid/model.fits")

        ELSE
            ! analytical density distribution, defined in the model_mod.f90 file
            DO i_cell = 1, grid%n_cell
                ! number of particles / cell / species;
                ! assumption:  density in each cell is constant and equal to the
                !              value at the cell center
                P_xy = sqrt(grid%cellmidcaco(i_cell, 1)**2 +                   &
                            grid%cellmidcaco(i_cell, 2)**2)
                P_z  = abs(grid%cellmidcaco(i_cell, 3))

                ! density at midpoint coordinate (number of particles / m^3)
                ! here, we set the density distribution for dust and H2 
                ! tbd: This is the same distribution for all elements, this
                !      should be generalized in future
                !      
                !
                ! set the density for the dust component
                
                !###############################################################
                ! add your custom density distribution here
                !

                grid%grd_dust_density(i_cell,:) = get_density(model,       &
                                    grid%cellmidcaco(i_cell, :)) * dust%sidi(:)
                    
                ! set the density for all other collision partner (H, He,...)
                grid%grd_col_density(i_cell,:)  = get_density(model,       &
                                                  grid%cellmidcaco(i_cell, :))

                !###############################################################
            END DO
        END IF

        ! now do loop over all cells and calculate remaining things
        ! here we set the molecule density with respect to the dust density 
        DO i_cell = 1, grid%n_cell

            ! number of particles / cell / species;
            ! assumption:  density in each cell is constant and equal to
            !              the value at the cell center
            P_xy = sqrt(grid%cellmidcaco(i_cell, 1)**2 +                       &
                        grid%cellmidcaco(i_cell, 2)**2)
            P_z  = abs(grid%cellmidcaco(i_cell, 3))
            
            ! set the density for the selected molecule
            ! don't load this from the old model, so we are able to distribute 
            ! the molecules in the way we want

            grid%grd_mol_density(i_cell)   = grid%grd_col_density(i_cell, 1) * &
                                             gas%mol_abund

            ! resulting number of particles/molecules

            grid%Nv(i_cell,:)     = grid%grd_dust_density(i_cell,:) *          &
                                        REAL(grid%cell_vol(i_cell),kind=r2)
            grid%Nv_mol(i_cell)   = grid%grd_mol_density(i_cell)    *          &
                                        REAL(grid%cell_vol(i_cell),kind=r2)
            grid%Nv_col(i_cell,:) = grid%grd_col_density(i_cell,:)  *          &
                                        REAL(grid%cell_vol(i_cell),kind=r2)
            ! set velocity, in a future release we should generalize this
            !
            IF ( .not. basics%old_model .and. .not. GetGridType(grid) == 9) THEN
                grid%velo(i_cell,:)  = get_velocity(model,                     &
                                                    grid%cellmidcaco(i_cell,:),&
                                                    'KEPLERIAN')
            END IF
            grid%absvelo(i_cell) = norm(REAL(grid%velo(i_cell,:),kind=r2))
        END DO
    END SUBROUTINE set_grid_properties

end module grd_mod
