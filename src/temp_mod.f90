MODULE temp_mod
      
    USE datatype
    USE var_global
    
    USE basic_type
    USE randgen_type
    USE dust_type
    USE grid_type
    USE model_type
    USE gas_type
    USE fluxes_type
    USE photon_type
    USE source_type

    USE model_mod, ONLY : get_temperature
    
    USE MCRT_mod, ONLY : MC_photon_transfer

    USE fileio


    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: set_temperature
    
CONTAINS

  ! ############################################################################
  ! simulation of primary source radiation:  
  ! ---------------------------------------
  !
  ! ---
  
    SUBROUTINE set_temperature(basics,grid,model,dust,gas,sources_in,fluxes)


    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Gas_TYP),INTENT(INOUT)                      :: gas
    TYPE(Fluxes_TYP) ,INTENT(INOUT)                  :: fluxes
    TYPE(SOURCES),INTENT(IN)                         :: sources_in
    !--------------------------------------------------------------------------!
    REAL(kind=r2), DIMENSION(1:3)                    :: caco
    REAL(kind=r2), DIMENSION(1:3)                    :: moco
    REAL(kind=r2)                                    :: min_cell_energy
    REAL(kind=r2)                                    :: v_turb
    
    INTEGER                                          :: i_cell
    INTEGER                                          :: i_dust
    INTEGER                                          :: i_tem
    INTEGER                                          :: i_a
    INTEGER                                          :: i_b
    INTEGER                                          :: i_c
    
    !--------------------------------------------------------------------------!
    IF (.not. basics%old_model) THEN
        IF ( basics%do_MC_temperature) THEN
            IF (.not. SourcesInitialized(sources_in)) THEN
                print *, 'ERROR: There is no source defined'
                stop
            END IF
            CALL MC_temp(basics,grid,model,dust,sources_in,fluxes)
        ELSE IF (.not. basics%do_MC_temperature) THEN
            print *, '| | using analytical temperature distribution'
            i_cell = 0
            DO i_a=1, grid%n(1)
                DO i_b=1, grid%n(2)
                    DO i_c=1, grid%n(3)
                        i_cell = i_cell + 1
                        ! set temperature analytical 
                        moco(1) = ( grid%co_mx_a( i_a ) + grid%co_mx_a( i_a  -1) ) / 2.0_r2
                        moco(2) = ( grid%co_mx_b(i_b) + grid%co_mx_b(i_b -1) ) / 2.0_r2
                        moco(3) = ( grid%co_mx_c(i_c) + grid%co_mx_c(i_c -1) ) / 2.0_r2
                        
                        caco = mo2ca(grid, moco)
                        DO i_dust = 1, dust%n_dust
                            grid%t_dust(i_cell, i_dust) = get_temperature(caco)
                        END DO
                    END DO
                END DO
            END DO
        END IF
    END IF
    
    DO i_cell = 1,grid%n_cell
        ! set correct internal cell energy
        DO i_dust = 1, dust%n_dust
            i_tem = min(int((grid%t_dust(i_cell, i_dust) -         &
                         basics%t_dust_min)/basics%d_tem), basics%n_tem)
            min_cell_energy = dust%QB(i_tem, i_dust)
            IF (grid%cell_energy(i_dust,i_cell) .lt. min_cell_energy) THEN
                grid%cell_energy(i_dust,i_cell) = min_cell_energy
            END IF
        END DO

        ! setting the gas temperature equal to the dust temperature
        ! we need to gerneralize this
!~         grid%t_gas(i_cell) = grid%t_dust(i_cell, 1)
        
    END DO
    
    ! CALCULATE PROBERTIES, WHICH NEED TEMPERATURE INFORMATIONS
    
    IF (basics%do_velo_ch_map) THEN
        !
        DO i_cell = 1,grid%n_cell

            ! set line width in each cell (in fact the inverse value)
            IF (grid%t_gas(i_cell) > 1.0e-2) THEN
                v_turb = 100.0_r2 ! 100.0 eq turb. line width TbD
                !v_turb = grid%velo(i_cell,3)! use i_z (as an indicator for turbulence)
                
                grid%cell_gauss_a(i_cell) = 1.0_r2 /                           &
                                (sqrt(2.0_r2 * con_k * grid%t_gas(i_cell) /    &
                                     (gas%mol_weight * 1.0e-3_r2 / con_Na) +   &
                                      v_turb))
            ELSE
                grid%cell_gauss_a(i_cell) = 0.0_r2
            END IF                
            ! save quadratic line width value for faster calculations
            grid%cell_gauss_a2(i_cell)  =   grid%cell_gauss_a(i_cell)**2
        END DO
    END IF
    END SUBROUTINE set_temperature

    SUBROUTINE MC_temp(basics, grid, model,                                    &
                       dust, sources_in, fluxes)

        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(SOURCES),INTENT(IN)                         :: sources_in
        
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_dust
        INTEGER                                          :: i_tem
        REAL(kind=r2)                                    :: min_cell_energy
                
        !----------------------------------------------------------------------!
        ! ---
        ! some preprations
        ! reset flux map
        fluxes%continuum_map(:,:,:,:) = 0.0_r2
        ! set internal cell energy
        DO i_dust = 1, dust%n_dust
            ! add a min temperature/energy (here 3 K)
            i_tem = int((3.0_r2 - basics%t_dust_min)/basics%d_tem)
            min_cell_energy = dust%QB(i_tem, i_dust)
            WHERE (grid%cell_energy(i_dust,:) .lt. min_cell_energy)
                grid%cell_energy(i_dust,:) = min_cell_energy
            END WHERE
            grid%cell_energy_sum(i_dust,:,1) = grid%cell_energy(i_dust,:) *    &
                                               (basics%PIx4 * grid%cell_vol(:))
        END DO
        ! photon transfer
        print *,'| | calculate temperature with Monte Carlo method'
        IF (basics%do_peel_off) PRINT *,"| | Peel-off technique enabled"
        CALL MC_photon_transfer(basics, grid,model,                            &
                                dust, sources_in, fluxes,                      &
                                deposit_energy=.True.,                         &
                                peel_off=basics%do_peel_off)
        ! prepare & save final results
        DO i_dust=1,dust%n_dust
            ! This is done, because in future we may want go back to pure lucy
            ! iterations, thus we have a total energy array
            ! 
            grid%cell_energy(i_dust,:) = grid%cell_energy_sum(i_dust,:,1)/ &
                                         (basics%PIx4 * grid%cell_vol(:))
            grid%cell_energy_sum(i_dust,:,:) = 0.0_r2
        END DO

        CALL save_continuum_map(model, basics, dust, fluxes, 1, basics%do_peel_off)

        CALL temp_final(basics, grid, dust)

    END SUBROUTINE MC_temp
      
  SUBROUTINE temp_final(basics, grid, dust)
    USE math_mod, ONLY : binary_search
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust

    !--------------------------------------------------------------------------!
    integer                                          :: i_dust, i_tem, nr_cell
    real(kind=r2)                                    :: hd2
    !--------------------------------------------------------------------------!
    
    ! ---
    print *, "| | final temperature calculation"
    hd2 = 0.0_r2
    do i_dust=1, dust%n_dust
        do nr_cell=1, grid%n_cell
            hd2    =  (grid%cell_energy(i_dust,nr_cell) )
         
            ! [2] find corresponding temperature from QB integral
            i_tem = MIN(binary_search(hd2,dust%QB(:,i_dust))-1,basics%n_tem)

            grid%t_dust(nr_cell,i_dust) = &
                ( ( (hd2 - dust%QB(i_tem-1,i_dust)) / (dust%QB(i_tem,i_dust) - &
                     dust%QB(i_tem-1,i_dust)) ) +                              &
                     real(i_tem-1,kind=r2) ) *                                 &
                     basics%d_tem 
        end do
    end do

    END SUBROUTINE temp_final
END MODULE temp_mod
