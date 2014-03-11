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
                Get_velo, &
                check_inside
    !--------------------------------------------------------------------------!
    
CONTAINS

  ! ################################################################################################
  ! check: is current position of photon stil inside the model space?
!  ! ---
    FUNCTION check_inside(caco,model) result(check_inside_result)
    
    IMPLICIT NONE
    !------------------------------------------------------------------------!
        
    TYPE(Model_TYP), INTENT(IN)                 :: model
    
    REAL(kind=r2), dimension(:), intent(in)   :: caco
    
    LOGICAL                                      :: check_inside_result   
    !------------------------------------------------------------------------!
    check_inside_result = (norm(caco) <= model%r_ou)
    
    END FUNCTION check_inside


  ! ################################################################################################
  ! make grid
  ! ---
    SUBROUTINE make_grid(basics, grid, model, dust, gas)
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        TYPE(Model_TYP), INTENT(IN)                 :: model
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
        REAL(kind=r1)                               :: hd2

        REAL(kind=r2), DIMENSION(:), ALLOCATABLE    :: hd_arr1
        
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
                end do
            end do
        end do
        grid%cell_minA(0) = model%r_in**2*PI
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
        
        print *, ''
        !------------------------------------------------------------------------!
        !    ! estimate mass
        !  scale to given mass, or to given density (at some point here 100 AU)
        mass_dens = .True.
!~         mass_dens = .False.
        IF (.not. basics%old_model) THEN
            IF (mass_dens) THEN
                hd_totalmass = 0.0_r2
                do i_dust=1, dust%n_dust
                    if (dust%sidi_par(i_dust,3)==0.0_r2) then
                        hd2 = dust%sidi_par(i_dust,1)
                    else
                        hd2 = (1.0_r2+dust%sidi_par(i_dust,3)) / (4.0_r2+dust%sidi_par(i_dust,3)) &
                        * (dust%sidi_par(i_dust,2)**(4.0_r2+dust%sidi_par(i_dust,3)) - &
                        dust%sidi_par(i_dust,1)**(4.0_r2+dust%sidi_par(i_dust,3))) &
                        / (dust%sidi_par(i_dust,2)**(1.0_r2+dust%sidi_par(i_dust,3)) - &
                        dust%sidi_par(i_dust,1)**(1.0_r2+dust%sidi_par(i_dust,3)))
                        hd2 = hd2**(1.0_r2/3.0_r2)
                    end if
                    hd_totalmass =  hd_totalmass + sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 * &
                                    hd2**3 * dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
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
        DO i_dust = 1, dust%n_dust
            grid%Nv_r(:,i_dust) = grid%Nv(:,i_dust) * dust%r_dust(i_dust)**2 * basics%PI2x4
            DO i_cell = 1, grid%n_cell
                IF (grid%Nv_r(i_cell,i_dust) <= 1.0_r2 ) THEN
                    grid%t_dust(i_cell,i_dust)      = basics%t_dust_min
                END IF
            END DO
        END DO
        
!~         print *,grid%grd_dust_density(get_cell_nr(grid,(/100.0_r2,0.0_r2,0.0_r2/)),1)
!~         print *,grid%Nv(get_cell_nr(grid,(/100.0_r2,0.0_r2,0.0_r2/)),1)
        ! dust mass verification
        hd_totalmass = 0.0_r2
        hd_dusttotalmass = 0.0_r2
        do i_dust=1, dust%n_dust
            if (dust%sidi_par(i_dust,3)==0.0_r2) then
                hd2 = dust%sidi_par(i_dust,1)
            else
                hd2 = (1.0_r2+dust%sidi_par(i_dust,3)) / (4.0_r2+dust%sidi_par(i_dust,3)) &
                * (dust%sidi_par(i_dust,2)**(4.0_r2+dust%sidi_par(i_dust,3)) - &
                dust%sidi_par(i_dust,1)**(4.0_r2+dust%sidi_par(i_dust,3))) &
                / (dust%sidi_par(i_dust,2)**(1.0_r2+dust%sidi_par(i_dust,3)) - &
                dust%sidi_par(i_dust,1)**(1.0_r2+dust%sidi_par(i_dust,3)))
                hd2 = hd2**(1.0_r2/3.0_r2)
            end if
            hd_dusttotalmass =  hd_dusttotalmass + sum(grid%Nv(:,i_dust)) * 4.0_r2*PI/3.0_r2 * &
                            hd2**3 * dust%den_dust(i_dust)*1.0e+3_r2 / M_sun
        end do
        
        ! estimate gas mass and rescale 
        ! Nv_mol,  defined molecule
        ! Nv_col,  collision partners 1 = H , 2 = H para ...,3 = H ortho
        !---- H2
        ! total gas mass
        hd_totalmass =  sum(grid%Nv_col(:,1)) * col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun + &
                        sum(grid%Nv_mol(:)) * gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
                 
        
        ! rescale total mass
!~         grid%grd_col_density(:,1) = grid%grd_col_density(:,1) * model%mass*gas%dust_mol_ratio/hd_totalmass
        grid%grd_col_density(:,1) = grid%grd_col_density(:,1) * hd_dusttotalmass*gas%dust_mol_ratio/hd_totalmass
        
        grid%Nv_col(:,1) = grid%grd_col_density(:,1) * REAL(grid%cell_vol(:),kind=r2)
!~         grid%Nv_col(:,1) =  grid%Nv_col(:,1) * &
!~                             model%mass*gas%dust_mol_ratio/(1.0_r2+gas%mol_abund)/hd_totalmass                      
        
!~         grid%grd_col_density(:,1) = grid%Nv_col(:,1)/grid%cell_vol(:)
        
        !---- H2 para

        grid%Nv_col(:,2) =  0.25*grid%Nv_col(:,1) 
                            
        grid%grd_col_density(:,2) = grid%Nv_col(:,2)/grid%cell_vol(:)
        
        !---- H2 ortho

        grid%Nv_col(:,3) =  0.75*grid%Nv_col(:,1) 
                            
        grid%grd_col_density(:,3) = grid%Nv_col(:,3)/grid%cell_vol(:)
        

        !---- selected molecule
        
        grid%grd_mol_density(:) = grid%grd_mol_density(:) * hd_dusttotalmass*gas%dust_mol_ratio/hd_totalmass
        grid%Nv_mol(:) = grid%grd_mol_density(:) * REAL(grid%cell_vol(:),kind=r2)
!~         hd_totalmass =  sum(grid%Nv_mol(:)) *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
        
!~         grid%Nv_mol(:) =    grid%Nv_mol(:) * &
!~                             model%mass*gas%dust_mol_ratio/(1.0_r2+1.0_r2/gas%mol_abund)/hd_totalmass
!~                             
!~         grid%grd_mol_density(:) =   grid%grd_mol_density(:) * &
!~                                     model%mass*gas%dust_mol_ratio/(1.0_r2+1.0_r2/gas%mol_abund)/hd_totalmass

        !TbD: all other partners:
        grid%Nv_col(:,4:6)            = 0.0_r2
        grid%grd_col_density(:,4:6)   = 0.0_r2
        grid%grd_col_density(0,:)     = 0.0_r2
        
        ! calculate total mass
        hd_totalmass =  &
                        ! dust mass
                        sum(grid%Nv(:,1)) * 4.0_r2*PI/3.0_r2 * &
                            hd2**3 * dust%den_dust(1)*1.0e+3_r2 / M_sun + &
                        ! H mass
                        sum(grid%Nv_col(:,1)) *col_p_weight(1)/con_Na*1.0e-3_r2/ M_sun + &
                        ! selected molecule mass
                        sum(grid%Nv_mol(:))     *gas%mol_weight/con_Na*1.0e-3_r2/ M_sun
        
!~         print *, grid%Nv_col(46,1)
        print '(A,ES11.4,A)',' dust mass           : '&
               ,sum(grid%Nv(:,1)) * 4.0_r2*PI/3.0_r2 * &
                            hd2**3 * dust%den_dust(1)*1.0e+3_r2 / M_sun, ' M_sun'

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
!~         stop
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
        !dz  = 2.0_r2*model%r_ou / real(grid%n(3), kind=r2)
        dth = 2.0_r2*PI / real(grid%n(2), kind=r2)
        
        grid%cell_minA(i_cell) = MIN( 0.5* dth *  (grid%co_mx_a(i_r)**2-grid%co_mx_a(i_r-1)**2), &
                                    ( dth * grid%co_mx_a(i_r-1) * dz), &
                                    ( (grid%co_mx_a(i_r)-grid%co_mx_a(i_r-1))*dz))
        grid%cell_minA(i_cell) = grid%cell_minA(i_cell)
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
        
        grid%cell_minA(i_cell) = MIN( &
                                grid%co_mx_a(i_r-1)**2.0_r2 * &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1) ) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 ) * &
                                abs(grid%co_mx_b(i_th) - grid%co_mx_b(i_th-1)) , &
                                (grid%co_mx_a(i_r)**2.0_r2-grid%co_mx_a(i_r-1)**2.0_r2 )* &
                                (grid%co_mx_c(i_ph) - grid%co_mx_c(i_ph-1) ) )
                               
        ! volume of individual cells [m^3]
        
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
        TYPE(Model_TYP), INTENT(IN)                 :: model
        TYPE(Basic_TYP), INTENT(IN)                 :: basics
        !------------------------------------------------------------------------!
        
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
        
        CALL save_boundaries(grid,model,basics)
    
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
        sf = 1.08_r2
!~         sf = grid%sf
        !dz = 2.0_r2*model%r_ou / real(grid%n(3), kind=r2)
        dz = model%r_ou * (sf-1.0_r2)/ (sf**((grid%n(3)-1)*0.5) - 1.0_r2)
        do i_z=1, int(grid%n(3)*0.5)
            grid%co_mx_c(int((grid%n(3)-1)*0.5+i_z+1)) =            dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
            grid%co_mx_c(int((grid%n(3)-1)*0.5-i_z)) =  -1.0_r2 * dz * (sf**i_z - 1.0_r2) / (sf-1.0_r2)
        end do
        grid%co_mx_c(int((grid%n(3)-1)*0.5)) = -1.0_r2/3.0_r2 * dz 
        grid%co_mx_c(int((grid%n(3)-1)*0.5+1)) = 1.0_r2/3.0_r2 * dz 

    END SUBROUTINE set_boundaries_cy
    
    SUBROUTINE set_boundaries_sp(grid,model)
    
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        
        TYPE(Grid_TYP), INTENT(INOUT)               :: grid
        TYPE(Model_TYP), INTENT(IN)                 :: model
        !------------------------------------------------------------------------!
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_ph
        
        REAL(kind=r2)                               :: dr1
        REAL(kind=r2)                               :: dth
        REAL(kind=r2)                               :: dph
        !------------------------------------------------------------------------!        
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
        !------------------------------------------------------------------------!
        
        INTEGER                                     :: i_a
        INTEGER                                     :: i_b 
        INTEGER                                     :: i_c 
        INTEGER                                     :: i_cell, i_cell_in
        INTEGER                                     :: io
        

        
        REAL(kind=r2), DIMENSION(1:3)               :: caco
        REAL(kind=r2), DIMENSION(1:3)               :: moco
        REAL(kind=r2), DIMENSION(10)                :: line
        REAL(kind=r2)                               :: R_gap_in, R_gap_ou
        REAL(kind=r1)                               :: P_xy, P_z
        CHARACTER(len=256)                          :: filename
        CHARACTER(len=256)                          :: waste
        !------------------------------------------------------------------------!
        ! ---

        CALL parse('R_gap_in',R_gap_in,'input/additional.dat')
        CALL parse('R_gap_ou',R_gap_ou,'input/additional.dat')

        IF (basics%old_model) THEN
        
            print '(2A)', '  loading model parameter from: ', TRIM(basics%pronam_old)
            ! load the density stored in a model file of the old project
            filename = TRIM(basics%path_results)//TRIM(basics%pronam_old)//'_model.dat'
            OPEN(unit=1, file=TRIM(filename), &
                action="read", status="unknown", form="formatted")
            READ(unit=1,fmt=*) waste
        END IF
        


        i_cell = 0
        do i_a=1, grid%n(1)
            do i_b=1, grid%n(2)
                do i_c=1, grid%n(3)
                    i_cell = i_cell + 1

                    ! number of particles / cell / species;
                    ! assumption:  density in each cell is constant and equal to the value in the cell center

                    ! midpoint coordinate
                    moco(1) = ( grid%co_mx_a( i_a ) + grid%co_mx_a( i_a  -1) ) / 2.0_r2
                    moco(2) = ( grid%co_mx_b(i_b) + grid%co_mx_b(i_b -1) ) / 2.0_r2
                    moco(3) = ( grid%co_mx_c(i_c) + grid%co_mx_c(i_c -1) ) / 2.0_r2
                    
                    caco = mo2ca(grid, moco)
                    P_xy = sqrt(caco(1)**2+caco(2)**2)
                    P_z  = abs(caco(3))
                    grid%cellmidcaco(i_cell,:) = caco
                    
                    IF (basics%old_model) THEN

                        READ(unit=1,fmt=*,iostat=io) i_cell_in, line
!~                         print *, i_cell_in, i_cell
                        IF (io < 0 .or. i_cell_in /= i_cell) THEN
                            PRINT *,'ERROR in model file ('//TRIM(filename)//') cell not found',i_cell
                            STOP
                        END IF
                        
                        ! set dust density, all other molecules can be set by the model_mod module
                        grid%grd_dust_density(i_cell,1) = line(4)
                        
                        grid%grd_col_density(i_cell,1:3) = line(6:8)
                        
                        ! set temperature
                        grid%t_dust(i_cell,1) = line(9)
                        grid%t_gas(i_cell)    = line(10)
                        

                    ELSE
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
                            grid%grd_dust_density(i_cell,:)    = get_den(model,caco(:))
                            
                            ! set the density for all other elements (H, He,...)
                            grid%grd_col_density(i_cell,:)     = get_den(model,caco(:))
                        END IF
                        !###########################################################################
                    END IF
                    
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
                    ! set velocity, in a future release we should generalize this as well somehow
                    ! 
                    grid%velo(i_cell,:)  = Get_velo(caco,model%kep_const)
                    
                    
                    grid%absvelo(i_cell) = norm(REAL(grid%velo(i_cell,:),kind=r2))
                    
                end do
            end do
        end do
        
        IF (basics%old_model) THEN
            CLOSE(unit=1)
        END IF
        
    END SUBROUTINE set_grid_properties
    
    
    PURE FUNCTION Get_velo(caco,kep_const) RESULT(velo)
    
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
        
    END FUNCTION Get_velo
    


end module grd_mod
