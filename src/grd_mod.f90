! ---
! grid specific routines
! (optimized for specific geometry)
! ---
MODULE grd_mod
    
    USE datatype
    USE var_globalnew
    USE basic_type
    USE dust_type
    USE grid_type
    USE gas_type
    USE model_type
    USE math_mod
    USE tools_mod
    USE fileio
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC ::   make_grid, &
                Set_velo, &
                Get_velo
    !--------------------------------------------------------------------------!
    
CONTAINS




  ! ################################################################################################
  ! make grid
  ! ---
    SUBROUTINE make_grid(basics, grid, model, dust, gas)
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        TYPE(Dust_TYP), INTENT(IN)                  :: dust
        TYPE(Gas_TYP), INTENT(INOUT)                :: gas
        
        !------------------------------------------------------------------------!

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
        
        LOGICAL                                     :: mass_dens
        !------------------------------------------------------------------------!
        
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
             
                    ! r/th/ph index => cell number
                    grid%cell_idx2nr(i_a,i_b,i_c) = i_cell 
             
                    ! cell number => r/th/ph index
                    grid%cell_nr2idx(1,i_cell) = i_a
                    grid%cell_nr2idx(2,i_cell) = i_b
                    grid%cell_nr2idx(3,i_cell) = i_c
                    
                    
                    ! find cell volume and min A cell-Wall
                    CALL calc_cell_properties(grid,model,i_a,i_b,i_c)
                    ! calculate midpoint coordinate
                    moco(1) = ( grid%co_mx_a( i_a ) + grid%co_mx_a( i_a  -1) ) / 2.0_r2
                    moco(2) = ( grid%co_mx_b(i_b) + grid%co_mx_b(i_b -1) ) / 2.0_r2
                    moco(3) = ( grid%co_mx_c(i_c) + grid%co_mx_c(i_c -1) ) / 2.0_r2
                        
                    ! get cartesian coordinate from moco
                    grid%cellmidcaco(i_cell,:) = mo2ca(grid, moco)
                    
                end do
            end do
        end do
        
        grid%cell_minA(0) = model%r_in**2*PI
!~         print *, minval(grid%cell_minA), grid%cell_minA(0)
!    ! fill grid with disk properties,e.g. temp, density, velocity...
        CALL set_grid_properties(basics,grid,gas,model)   
    ! verification: maximum cell number (i.e., total number of cells) = total number of cells
        if ( i_cell /= grid%n_cell ) then
            print *,"subroutine mk_grd(): wrong number of ESCs: ", i_cell
            print *,"should be                                : ", grid%n_cell 
            !call stop_mc3d()!TBD!
        end if
       
!   ! verification
        hd1 = sum(grid%cell_vol(:) / model%ref_unit**3)
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            hd2 = (4.0_r2/3.0_r2) * PI * (model%r_ou**3 - model%r_in**3)
        CASE('cylindrical')
            hd2 = 2.0_r2 * model%r_ou * PI *(model%r_ou**2 - model%r_in**2)
        CASE('cartesian')
            print *, 'TbD,..verification of volume space'
            stop
        CASE DEFAULT
            print *, 'selected coordinate system not found, verification of volume space'
            stop
        END SELECT
!~         print *, abs(hd1-hd2)/hd1
        if ( abs(hd1-hd2)/hd1 > 0.01_r2 ) then
            print *,"!!! Warning:in subroutine make_grid() "
            print *,"    Difference between the VOLUME OF THE MODEL SPACE"
            print *,"    and the TOTAL VOLUME OF THE SINGLE ESCs is larger than 1% :"
            print *,"    - volume( model space)        : ", hd2, " AU^3"
            print *,"    - total volume of single ESCs : ", hd1, " AU^3"
            print *
            print *,"mc3d stopped."
            stop
        end if
        
        print *, '                                                               '
        !------------------------------------------------------------------------!
        !    ! estimate mass
        !  scale to given mass, or to given density (at some point here 100 AU)
        mass_dens = .True.
!~         mass_dens = .False.
        IF (.not. basics%old_model) THEN
            IF (mass_dens) THEN
                hd_totalmass = 0.0_r2
                do i_dust=1, dust%n_dust
!~                     stop
                    hd_totalmass =  hd_totalmass + sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 * &
                                    dust%r_dust(i_dust)**3 * dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
                end do
                ! rescale mass to user-defined value
                grid%grd_dust_density(:,:) = grid%grd_dust_density(:,:) * model%mass/hd_totalmass
                grid%Nv(:,:)               = grid%Nv(:,:)               * model%mass/hd_totalmass
            ELSE
                i_cell = get_cell_nr(grid,(/100.0_r2,0.0_r2,0.0_r2/))
                hd_totalmass = 22895.6132878962_r2/grid%grd_dust_density(i_cell,1)
                grid%Nv(:,:)               = grid%Nv(:,:) * hd_totalmass
                grid%grd_dust_density(:,:) = grid%grd_dust_density(:,:) * hd_totalmass

            END IF
        END IF
        ! for immediate_temp
!~         DO i_dust = 1, dust%n_dust
!~             grid%Nv_r(:,i_dust) = grid%Nv(:,i_dust) * dust%r_dust(i_dust)**2 * basics%PI2x4
!~             DO i_cell = 1, grid%n_cell
!~                 IF (grid%Nv_r(i_cell,i_dust) <= 1.0_r2 ) THEN
!~                     grid%t_dust(i_cell,i_dust)      = basics%t_dust_min
!~                 END IF
!~             END DO
!~         END DO
        

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
        grid%Nv_col(:,4:6)            = 0.0_r2
        grid%grd_col_density(:,4:6)   = 0.0_r2
        grid%grd_col_density(0,:)     = 0.0_r2
        
        ! calculate total mass
        hd_totalmass =  &
                        ! dust mass
                        sum(grid%Nv(:,1)) * 4.0_r2*PI/3.0_r2 * &
                            dust%r_dust(1)**3 * dust%den_dust(1)*1.0e+3_r2 / M_sun + &
                        ! H mass
                        sum(grid%Nv_col(:,1)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun + &
                        ! selected molecule mass
                        sum(grid%Nv_mol(:))     *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
        
!~         print *, grid%Nv_col(46,1)
        print '(A,ES11.4,A)',' dust mass           : '&
               ,sum(grid%Nv(:,1)) * 4.0_r2*PI/3.0_r2 * &
                            dust%r_dust(1)**3 * dust%den_dust(1)*1.0e+3_r2 / M_sun, ' M_sun'

        print '(A,ES11.4,A)',' H2 ortho mass       : '&
              ,sum(grid%Nv_col(:,3)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
              
        print '(A,ES11.4,A)',' H2 para mass        : '&
              ,sum(grid%Nv_col(:,2)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
              
        print '(A,ES11.4,A)',' H2 total mass       : '&
              ,sum(grid%Nv_col(:,1)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun, ' M_sun'
                         
        print '(A,ES11.4,A)',' sel molecule mass   : '&
              ,sum(grid%Nv_mol(:)) *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun, ' M_sun'  
        
        print '(A)', ' ---------------------------------------'
        print '(A,ES11.4,A)',' total disk mass     : ',hd_totalmass, ' M_sun'
        print *,''
        !------------------------------------------------------------------------!
        ! ---
        ! set smallest step width for photon transfer
        ! - goal: avoid step widths smaller than the amount that will change a r2 type floating point
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
            print '(A,ES11.4,A,A)'," smallest step width : ", grid%d_l_min,' ', GetModelName(model)
        end if

        deallocate( hd_arr1 )
        
        
    END SUBROUTINE make_grid
    
    !  ! ################################################################################################
    
    SUBROUTINE calc_cell_properties(grid,model,i_a,i_b,i_c)
        !generic routine
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !------------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_a,i_b,i_c
        !------------------------------------------------------------------------!
        SELECT CASE(GetGridName(grid))
        CASE('spherical')
            CALL calc_cell_properties_sp(grid,model,i_a,i_b,i_c)
            
        CASE('cylindrical')
            CALL calc_cell_properties_cy(grid,model,i_a,i_b,i_c)
            !print *, 'TbD,..calc_cell_properties'
            !stop
            
        CASE('cartesian')
            print *, 'TbD,..calc_cell_properties'
            stop      
             
        CASE DEFAULT
            print *, 'selected coordinate system not found, calc_cell_properties'
            stop
        END SELECT
        
        
    END SUBROUTINE calc_cell_properties
    
    
    SUBROUTINE calc_cell_properties_cy(grid,model,i_r,i_ph,i_z)

        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)              :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !------------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_r,i_ph,i_z
        INTEGER                                     :: i_cell
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dz
        !------------------------------------------------------------------------!
        
        i_cell = grid%cell_idx2nr(i_r,i_ph,i_z)
        dz = abs(grid%co_mx_c(i_z-1) - grid%co_mx_c(i_z))
        dth = 2.0_r2*PI / real(grid%n(2), kind=r2)
        
        grid%cell_minA(i_cell) = MIN( 0.5* dth *  (grid%co_mx_a(i_r)**2-grid%co_mx_a(i_r-1)**2), &
                                    ( dth * grid%co_mx_a(i_r-1) * dz), &
                                    ( (grid%co_mx_a(i_r)-grid%co_mx_a(i_r-1))*dz))
                                    
!~         grid%cell_minA(i_cell) = grid%cell_minA(i_cell)*1.0e3
        ! volume of individual cells [m^3]
        
        grid%cell_vol(i_cell) = 0.5_r2 * dz * dth * &
                               (grid%co_mx_a(i_r)**2-grid%co_mx_a(i_r-1)**2) * &
                                model%ref_unit**3
!~         print *,grid%cell_vol(i_cell)
    END SUBROUTINE calc_cell_properties_cy
    
    
    
    
    
    SUBROUTINE calc_cell_properties_sp(grid,model,i_r,i_th,i_ph)

        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !------------------------------------------------------------------------!
        INTEGER,INTENT(IN)                         :: i_r,i_th,i_ph
        INTEGER                                     :: i_cell
        !------------------------------------------------------------------------!
        
        i_cell = grid%cell_idx2nr(i_r,i_th,i_ph)
!~         grid%cell_minA(i_cell) = ( &
!~                                 grid%co_mx_a(i_r-1)**2.0_r2 * &
!~                                 (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) * &
!~                                 abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1) ) + &
!~                                 (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 ) * &
!~                                 abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1)) + &
!~                                 (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 )* &
!~                                 (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) )/3.0_r2
        grid%cell_minA(i_cell) = MIN( &
                                grid%co_mx_a(i_r-1)**2.0_r2 * &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1) ) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1)) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 )* &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) )
!~         print *, grid%co_mx_a(i_r-1)
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
    
    
    
    
    !  ! ################################################################################################    
    SUBROUTINE set_boundaries(grid,model,basics)
        !generic routine
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        CHARACTER(len=252)                          :: file_a,file_b,file_c
        !------------------------------------------------------------------------!
        IF (basics%old_model) THEN
            ! if we use an old model, we simply read the boundaries from the input files
            file_a = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_a_boundaries.dat'
            file_b = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_b_boundaries.dat'
            file_c = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_c_boundaries.dat'
            CALL read_boundaries(grid,file_a,file_b,file_c)
        
        ELSE
            SELECT CASE(GetGridName(grid))
            
            CASE('spherical')
                CALL set_boundaries_sp(grid,model)
                
            CASE('cylindrical')
                CALL set_boundaries_cy(grid,model)
                
            CASE('cartesian')
                print *, 'TbD,..set_boundaries'
                stop
            CASE DEFAULT
                print *, 'selected coordinate system not found, set_boundaries'
                stop
            END SELECT
        END IF
        
        CALL save_boundaries(grid,basics)
    
    END SUBROUTINE set_boundaries
    
    SUBROUTINE set_boundaries_cy(grid,model)
    
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !------------------------------------------------------------------------!
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_z
        
        REAL(kind=r2)                               :: dr1
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dz
        REAL(kind=r2)                               :: sf
        !------------------------------------------------------------------------!        
        
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
            !      -> we use sf too, but not the user given value but a fixed one
            
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
            
            
        END SELECT

    END SUBROUTINE set_boundaries_cy
    
    SUBROUTINE set_boundaries_sp(grid,model)
    
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(INOUT)              :: model
        !------------------------------------------------------------------------!
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_ph
        INTEGER                                     :: i_abc
        INTEGER                                     :: i
        
        REAL(kind=r2)                               :: dr1
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dph
        REAL(kind=r2)                               :: trash
        REAL(kind=r2)                               :: value_in
        
        Character(256)                              :: waste
        !------------------------------------------------------------------------!
        SELECT CASE(GetGridType(grid))
        
        CASE(1)
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
            ! 1.2. theta
            ! co_mx_b(  0 ) = -PI/2 (-90°)
            ! co_mx_b(n(2)) = +PI/2 (+90°)
            grid%co_mx_b(0) = -PI/2.0_r2

            dth = PI / real(grid%n(2), kind=r2)
            do i_th=1, grid%n(2)
                grid%co_mx_b(i_th) = grid%co_mx_b(0) +  dth * real(i_th, kind=r2)
            end do
        
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
                grid%co_mx_b(i_th) = grid%co_mx_b(0) +  dth * real(i_th, kind=r2)
            end do
            
            ! ---
            ! 1.3 phi
            ! co_mx_c(  0 ) = 0    (  0°)
            ! co_mx_c(n(3)) = 2xPI (360°)
            grid%co_mx_c(0) = 0.0_r2
        
            dph = 2.0_r2*PI / real(grid%n(3), kind=r2)
            do i_ph=1, grid%n(3)
                grid%co_mx_c(i_ph) =  grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do
            
        CASE(3)
            print *, 'pluto interface'
            ! here we have a linear spherical grid
            ! this is motivated by a given density distribution from the pluto code
            ! to save endless memory, the grid is equally spaced, but only between 
            ! in theta direction, r and ph are unchanged
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
            IF (grid%n(1) /= 100) print *, 'WARNING: no of r angles should be 100'
            do i_r=1, grid%n(1)
                grid%co_mx_a(i_r) = grid%co_mx_a(0) + dr1 * real(i_r, kind=r2)
            end do
            
            ! ---
            ! 1.2. theta
            ! co_mx_b(  0 ) = -PI/2 (-90°)
            ! co_mx_b(n(2)) = +PI/2 (+90°)
            
            grid%co_mx_b(0)         = -PI/2.0_r2
            grid%co_mx_b(grid%n(2)) =  PI/2.0_r2
            
            ! give boundaries, derived from input file, check it!!!
            grid%co_mx_b(1)            = -0.30099091
            grid%co_mx_b(grid%n(2)-1)  =  0.30099091
            IF (grid%n(2) /= 65) print *, 'WARNING: no of th angles should be 65'

            dth = (grid%co_mx_b(grid%n(2)-1) - grid%co_mx_b(1)) / real(grid%n(2)-2, kind=r2)
            do i_th=2, grid%n(2)-1
                grid%co_mx_b(i_th) = grid%co_mx_b(1) +  dth * real(i_th-1, kind=r2)
            end do
            
            ! ---
            ! 1.3 phi
            ! co_mx_c(  0 ) = 0    (  0°)
            ! co_mx_c(n(3)) = 2xPI (360°)
            grid%co_mx_c(0) = 0.0_r2
            IF (grid%n(3) /= 124) print *, 'WARNING: no of ph angles should be 124'
            dph = 2.0_r2*PI / real(grid%n(3), kind=r2)
            do i_ph=1, grid%n(3)
                grid%co_mx_c(i_ph) =  grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do
            
        CASE(4)
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
                grid%co_mx_c(i_ph) =  grid%co_mx_c(0) + dph * real(i_ph, kind=r2)
            end do
        
        
        CASE(5)
            ! This is an old and very rubbish implemented version. it should work, but...;)
            !
            ! user given spherical grid, here interface to PLUTO code
            ! read grid from file, as ph and th are linear spaced, the no of cells are
            ! allready calculated und we have just du space them evenly
            ! the r cells are taken from the input file
            ! here we have to reset the given modelparameter
            ! let's discuss this later
            open(unit=1, file="input/grid/pluto_disk.dat", &
                action="read", status="unknown", form="formatted")

            DO i = 1,23
                read(unit=1,fmt=*) waste
            END DO
            
            ! make r boundaries
            DO i_abc = 0,grid%n(1)-1
                read(unit=1,fmt=*) i, trash, grid%co_mx_a(i_abc)
            END DO
            ! transform cm to AU
            grid%co_mx_a(:) = grid%co_mx_a(:)/con_au/100.0_r2
            
            ! we assume an exponential increasing r coordinate, therefore we can calculate 
            ! the effective outer boundary == r_ou
            ! 
            grid%co_mx_a(grid%n(1)) = exp((log(grid%co_mx_a(1))-log(grid%co_mx_a(0)))*(grid%n(1)+1)+ &
                                       log(grid%co_mx_a(0)))
            
            ! now adjust the innercoordinate to the user given inner rim
            ! please be aware, that this includes the assumption, that the given PLUTO disk 
            ! model is somehow scalable (self similar)
            grid%co_mx_a = grid%co_mx_a * model%r_in/grid%co_mx_a(0)
            
            ! now adjust the outer radius the given one
            
            model%r_ou = grid%co_mx_a(grid%n(1))
            print *, model%r_ou
            
            ! make th boundaries
            grid%co_mx_b(0) = -PI/2.0_r2
            
            DO i_abc=1, grid%n(2)
                grid%co_mx_b(i_abc) = grid%co_mx_b(0) +  PI / real(grid%n(2), kind=r2) * real(i_abc, kind=r2)
            end do

            ! make ph boundaries
            grid%co_mx_c(0) = 0.0
            DO i_abc=1, grid%n(3)
                grid%co_mx_c(i_abc) = grid%co_mx_c(0) +  2.0*PI / real(grid%n(3), kind=r2) * real(i_abc, kind=r2)
            end do

            close(unit=1)

        CASE(9)
            ! ---
            ! This case is for a very general input of a spherical grid, coordinates defined in external files
            ! The style of the input file has to be consistent
            
            CALL read_boundaries(grid, "input/grid/spherical_r.dat", &
                                 "input/grid/spherical_theta.dat","input/grid/spherical_phi.dat")

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
                PRINT *, "r_ou             : ", model%r_in
                PRINT *, "outer cell bound : ", grid%co_mx_a(grid%n(1))
                STOP
            END IF
 

            ! ---
            ! 2 theta
            ! test for consistency
            ! co_mx_b(  0 ) = -PI/2    (- 90°)
            ! co_mx_b(n(2)) =  PI/2    (  90°)

            IF ( abs(PI/2 + grid%co_mx_b(0)) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: lower theta coordinate is not -PI/2"
                PRINT *, "lowest cell boundary : ", grid%co_mx_b(0)
                STOP
            ELSEIF ( abs(PI/2 - grid%co_mx_b(grid%n(2))) .gt. 1.0e3_r2*epsilon(PI) ) THEN
                PRINT *, "ERROR: upper theta coordinate is not PI/2"
                PRINT *, "upper cell boundary : ", grid%co_mx_b(grid%n(2))
                STOP
            END IF
            
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


        END SELECT
    
    END SUBROUTINE set_boundaries_sp
    

    !  ! ################################################################################################


    SUBROUTINE set_grid_properties(basics,grid,gas,model)
    ! set grid density and temperature! (and velocity) 
    ! prepared for gerneral coordinates
    
        USE model_mod
        USE parser_mod
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(GAS_TYP), INTENT(INOUT)                :: gas
        TYPE(Model_TYP), INTENT(IN)                 :: model
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        TYPE(Vector3d)                              :: p_vec
        !------------------------------------------------------------------------!
        
        INTEGER                                     :: i_cell, i_cell_in, i, k
        INTEGER                                     :: i_in, i_r, i_th, i_ph
        INTEGER                                     :: io
        INTEGER,DIMENSION(1:3)                      :: pluto_n
        
        
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE      :: pluto_r
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE      :: pluto_ph
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE      :: pluto_th
        
        REAL(kind=r2), DIMENSION(10)                :: line
        REAL(kind=r2), DIMENSION(3)                 :: moco, velo_hlp
        REAL(kind=r2), DIMENSION(3,3)               :: S
        REAL(kind=r2)                               :: R_gap_in, R_gap_ou
        REAL(kind=r2)                               :: value_in
        REAL(kind=r1)                               :: P_xy, P_z
        CHARACTER(len=256)                          :: filename
        CHARACTER(len=256)                          :: waste
        
        !------------------------------------------------------------------------!
        ! ---

        CALL parse('R_gap_in',R_gap_in,'input/additional.dat')
        CALL parse('R_gap_ou',R_gap_ou,'input/additional.dat')
            
        IF (basics%old_model) THEN
        
            print '(2A)', ' | loading model parameter from: ', TRIM(basics%pronam_old)
            ! load the density stored in a model file of the old project
            filename = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_model.dat'
            OPEN(unit=1, file=TRIM(filename), &
                action="read", status="unknown", form="formatted")
            READ(unit=1,fmt=*) waste
            k = 1
            DO i_cell = 1, grid%n_cell
                IF (i_cell >= int(k*grid%n_cell*0.01)) THEN
                    WRITE (*,'(A,I3,A)') ' | | ',int(i_cell/real(grid%n_cell)*100.0),' % done'//char(27)//'[A'
                    k = k + 1
                END IF
            
                READ(unit=1,fmt=*,iostat=io) i_cell_in, line
                IF (io < 0 .or. i_cell_in /= i_cell) THEN
                    PRINT *,'ERROR in model file ('//TRIM(filename)//') cell not found',i_cell
                    STOP
                END IF
                
                ! set dust density and H2, the observed molecules can be set by the model module
                grid%grd_dust_density(i_cell,1) = line(4)
                
                grid%grd_col_density(i_cell,1:3) = line(6:8)
                
                ! set temperature
                grid%t_dust(i_cell,1) = line(9)
                grid%t_gas(i_cell)    = line(10)

            END DO
            print *,"                          "
            CLOSE(unit=1)
        ELSE IF (basics%pluto_data) THEN
            OPEN(unit=1, file="input/grid/sp1_reduced.dat", &
                action="read", status="unknown", form="formatted")
            ! read density from provided file
            ! first read header
            DO i = 1,4
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            READ(unit=1,fmt=*,iostat=io) pluto_n(1)
            READ(unit=1,fmt=*,iostat=io) waste
            READ(unit=1,fmt=*,iostat=io) pluto_n(2)
            READ(unit=1,fmt=*,iostat=io) pluto_n(3)
            
            ALLOCATE( pluto_r(1:pluto_n(1)), pluto_th(1:pluto_n(2)),pluto_ph(1:pluto_n(3)))
            
            
            DO i = 1,3
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            
            ! read r coordinates
            DO i = 1,pluto_n(1)
                READ(unit=1,fmt=*,iostat=io) i_in, pluto_r(i)
            END DO
            pluto_r = pluto_r*grid%sf
            DO i = 1,2
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            
            ! read th coordinates and correct for mc3d grid with pi/2
            DO i = 1,pluto_n(2)
                READ(unit=1,fmt=*,iostat=io) i_in, pluto_th(i)
            END DO
            pluto_th = pluto_th - PI/2.0
            DO i = 1,2
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            
            !read ph coordinates
            DO i = 1,pluto_n(3)
                READ(unit=1,fmt=*,iostat=io) i_in, pluto_ph(i)
            END DO
            ! now read density distribution
            
            DO i = 1,3
                READ(unit=1,fmt=*,iostat=io) waste
            END DO
            
            DO i = 1,pluto_n(1)*pluto_n(2)*pluto_n(3)
                
                READ(unit=1,fmt=*,iostat=io) i_r, i_th, i_ph, value_in

                moco = (/pluto_r(i_r),pluto_th(i_th),pluto_ph(i_ph)/)
                i_cell = get_cell_nr(grid,mo2ca(grid,moco))
                    
!~                 IF (i_ph == 1) THEN    
                
                grid%grd_dust_density(i_cell,:) = value_in
                grid%grd_col_density(i_cell,:) = value_in
!~                 END IF
                
            END DO
            
            DEALLOCATE(pluto_r,pluto_th,pluto_ph)
            CLOSE(unit=1)
        ELSE IF (GetGridType(grid) == 4 .and. GetGridName(grid) == 'spherical' ) THEN
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
!~                 print *, i
                moco(1) = 10.0**(moco(1))
                i_cell = get_cell_nr(grid,mo2ca(grid,moco))
                grid%grd_dust_density(i_cell,:) = 10.0**(value_in)
                grid%grd_col_density(i_cell,:)  = 10.0**(value_in)
            END DO
            CLOSE(unit=1)
        ELSE IF (GetGridType(grid) == 9 .and. GetGridName(grid) == 'spherical' ) THEN            
            print *,'  read model from text_file'
            OPEN(unit=1, file="input/grid/model.dat", &
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
                            WRITE (*,'(A,I3,A)') ' | | | ',int(i_cell/real(grid%n_cell)*100.0),' &
                                                    &% done'//char(27)//'[A'
                            k = k +1
                        END IF
                        p_vec%comp = ca2sp(grid%cellmidcaco(i_cell,:))
                        READ(unit=1,fmt=*,iostat=io) pluto_n, velo_hlp, value_in
                        grid%grd_dust_density(i_cell,:) = value_in
                        grid%grd_col_density(i_cell,:)  = value_in
                        
                        CALL generate_rotation_matrix(p_vec,S)
                        grid%velo(i_cell,:)  =  matmul(S,velo_hlp)
                    END DO
                END DO
            END DO
        
            CLOSE(unit=1)
        
        
        ELSE 
            DO i_cell = 1, grid%n_cell
                ! number of particles / cell / species;
                ! assumption:  density in each cell is constant and equal to the value in the cell center
                
                P_xy = sqrt(grid%cellmidcaco(i_cell,1)**2+grid%cellmidcaco(i_cell,2)**2)
                P_z  = abs(grid%cellmidcaco(i_cell,3))
                ! density at midpoint coordinate (number of particles / m^3)
                ! here, we set the density distribution for dust and H2 
                ! tbd: This is the same distribution for all elements, this should be generalized
                !      in future
                !
                ! set the density for the dust component
                !###########################################################################
                ! add your custom density distribution here
                !
                IF (P_xy .lt. R_gap_in .or. P_xy .gt. R_gap_ou) THEN
!~                 IF (abs(atan(P_z/P_xy)) < 0.30 ) THEN
                    grid%grd_dust_density(i_cell,:)    = get_den(model,grid%cellmidcaco(i_cell,:))
                    
                    ! set the density for all other elements (H, He,...)
                    grid%grd_col_density(i_cell,:)     = get_den(model,grid%cellmidcaco(i_cell,:))
                ELSE
                    grid%grd_dust_density(i_cell,:)    = 0.0_r2
                    grid%grd_col_density(i_cell,:)     = 0.0_r2
!~                 END IF
                
                END IF
                !###########################################################################
            END DO
        END IF

        
        ! now do loop over all cells and calculate remaining things
        ! here we set the molecule density with respect to the input dust density 
        DO i_cell = 1, grid%n_cell

            ! number of particles / cell / species;
            ! assumption:  density in each cell is constant and equal to the value in the cell center
            P_xy = sqrt(grid%cellmidcaco(i_cell,1)**2+grid%cellmidcaco(i_cell,2)**2)
            P_z  = abs(grid%cellmidcaco(i_cell,3))
            
            ! set the density for the selected molecule
            ! don't load this from the old model, so we are able to distribute the molecules
            ! in the way we want
!~                     IF (P_z .gt. P_xy**1.8/1000. .and. P_z .lt. P_xy**1.2/5.) THEN
            grid%grd_mol_density(i_cell)   = grid%grd_col_density(i_cell,1)*gas%mol_abund
!~                     ELSE
!~                         grid%grd_mol_density(i_cell)   = 0.0
!~                     END IF


            ! resulting number of particles/molecules

            grid%Nv(i_cell,:)     = grid%grd_dust_density(i_cell,:) * REAL(grid%cell_vol(i_cell),kind=r2)
            grid%Nv_mol(i_cell)   = grid%grd_mol_density(i_cell)    * REAL(grid%cell_vol(i_cell),kind=r2)
            grid%Nv_col(i_cell,:) = grid%grd_col_density(i_cell,:)  * REAL(grid%cell_vol(i_cell),kind=r2)
            ! set velocity, in a future release we should generalize this
            !
            IF (GetGridType(grid) /= 9) THEN
                grid%velo(i_cell,:)  = Set_velo(grid%cellmidcaco(i_cell,:),model%kep_const)
            END IF
            grid%absvelo(i_cell) = norm(REAL(grid%velo(i_cell,:),kind=r2))
        END DO
    END SUBROUTINE set_grid_properties
    
    
    PURE FUNCTION Set_velo(caco,kep_const) RESULT(velo)
    
        ! define your velocity distribution here! velo in cartesian coordinates
        !
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)                 :: caco
        REAL(kind=r1), DIMENSION(1:3)                            :: velo
        REAL(kind=r2)                                            :: konst, expp, r
        REAL(kind=r1),INTENT(IN)                                 :: kep_const
        !------------------------------------------------------------------------!
        ! we assume pure keplerian rotation, therefore v(r) ~ r**-0.5
!~         konst = 26000.0_r2 
        konst = kep_const  !
        
        expp  = -1.5  ! is a result of keplerian rotation and coordinate system conversion
        r     = (caco(1)**2+caco(2)**2)**(expp*0.5)*konst
        
        velo(1) = (-1.0) * caco(2) * r
        velo(2) =          caco(1) * r
        velo(3) = 0.0
        
    END FUNCTION Set_velo
    

    ELEMENTAL FUNCTION Get_velo(velo1,velo2,x) RESULT(velo_x)
    
        ! define your velocity distribution here! velo in cartesian coordinates
        !
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL(kind=r2),INTENT(IN)                           :: velo1,velo2
        REAL(kind=r2)                                      :: velo_x
        REAL(kind=r2),INTENT(IN)                           :: x
        !------------------------------------------------------------------------!

        velo_x = (velo2-velo1)*x + velo1
        
        
    END FUNCTION Get_velo
    


end module grd_mod
