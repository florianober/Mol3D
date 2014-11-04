MODULE temp_mod
      
    USE datatype
    USE var_globalnew
    
    USE basic_type
    USE randgen_type
    USE dust_type
    USE grid_type
    USE model_type
    USE gas_type
    USE fluxes_type
    USE photon_type
    USE source_type
    
    USE transfer_mod
    USE start_mod
    USE interact_mod
    USE immediate_mod
    
    USE math_mod
    USE tools_mod
    USE fileio


    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: set_temperature!, primary_scatt, dust_RT
    
CONTAINS

  ! ################################################################################################
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
    
    INTEGER                                          :: i_cell
    INTEGER                                          :: i_a
    INTEGER                                          :: i_b
    INTEGER                                          :: i_c
    
    !--------------------------------------------------------------------------! 
    IF (.not. basics%old_model) THEN
        IF ( basics%calc_tmp) THEN
            IF (.not. SourcesInitialized(sources_in)) THEN
                print *, 'ERROR: There is no source defined'
                stop
            END IF
            CALL primary_temp(basics,grid,model,dust,sources_in,fluxes)
        ELSE IF (.not. basics%calc_tmp) THEN
            print *, '| | using analytical temperature distribution'
            i_cell = 0
            DO i_a=1, grid%n(1)
                DO i_b=1, grid%n(2)
                    DO i_c=1, grid%n(3)
                        i_cell = i_cell + 1
                        ! set temperature ( analytical(done) or import from mc3d (TbD!~!) )
                        moco(1) = ( grid%co_mx_a( i_a ) + grid%co_mx_a( i_a  -1) ) / 2.0_r2
                        moco(2) = ( grid%co_mx_b(i_b) + grid%co_mx_b(i_b -1) ) / 2.0_r2
                        moco(3) = ( grid%co_mx_c(i_c) + grid%co_mx_c(i_c -1) ) / 2.0_r2
                        
                        caco = mo2ca(grid, moco)
                        grid%t_dust(i_cell,1)  = Get_temp(caco)

                    END DO
                END DO
            END DO
        END IF
        
    END IF
    
    DO i_cell = 1,grid%n_cell
        IF (grid%grd_dust_density(i_cell,1) <= 1.0e-34 ) THEN
            grid%t_dust(i_cell,1) = 0.0
        END IF
        grid%t_gas(i_cell) = grid%t_dust(i_cell,1)! setting the gas temperature equal to the dust temperature
        
        ! include some kind of freeze out Temperature (TbD: define this more generally)

!~         IF ( grid%t_gas(i_cell) <= 35.0 .or. grid%t_gas(i_cell) >= 60.0 ) THEN
!~             grid%grd_mol_density(i_cell)   = 0.0
!~         END IF
    END DO
    
    CALL sv_temp(basics, grid)      ! The propose of this routine has changed, nowadays it
                                    ! is just for saving the x midplane temperature

    
    ! CALCULATE PROBERTIES, WHICH NEED TEMPERATURE INFORMATIONS
    
    IF (basics%do_raytr) THEN
        !
        DO i_cell = 1,grid%n_cell

            ! set line width in each cell (in fact the inverse value)
            grid%cell_gauss_a(i_cell) = 1.0_r2/(sqrt(2.0_r2*con_k*grid%t_dust(i_cell,1)/ &  
                                                (gas%mol_weight*1.0e-3_r2/con_Na) &
                                                +100.0_r2**2)) ! 100.0 eq turb. line width TbD
                                                !+grid%velo(i_cell,3)**2))   ! use i_z (as an indicator for turbolence)
                                                !+grid%v_turb(i_cell)**2))   ! use turb array (TbD)
                                                
            ! save quadratic line width value for faster calculations
            
            grid%cell_gauss_a2(i_cell)  =   grid%cell_gauss_a(i_cell)**2
            
        END DO
    END IF
    END SUBROUTINE set_temperature
  
    PURE FUNCTION Get_temp(caco) RESULT(temps)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)                  :: caco
        REAL(kind=r2)                                             :: temps
        REAL(kind=r2)                                             :: konst, P_xy, P_h
        
        !------------------------------------------------------------------------!

            !define your temperature distribution here!
            !print *, 'using an analytical temperature distribution, isotherm in z direction'
            konst = 100.0_r2
            
            P_xy  = sqrt(caco(1)**2+caco(2)**2)
            P_h     = 5.0_r2 * (P_xy/100.0_r2)**1.125_r2
            
           ! temps   = 500.0_r2*P_xy**(-0.5) * (2.0-exp(-0.5*(abs(caco(3))/P_h)**2))
            temps   = 200.0_r2*P_xy**(-0.5)

    END FUNCTION Get_temp
  
    SUBROUTINE primary_temp(basics,grid,model,dust,sources_in,fluxes)


        IMPLICIT NONE
        
        !--------------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                   :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(SOURCES),INTENT(IN)                         :: sources_in
    
        TYPE(Randgen_TYP)                                :: rand_nr
        
        !--------------------------------------------------------------------------!  
        INTEGER                                          :: i_phot, i_dust
        INTEGER                                          :: seed, k_phot
        !$ INTEGER omp_get_thread_num 
                
        TYPE(PHOTON_TYP)                                 :: photon
        
        !--------------------------------------------------------------------------! 
        ! ---
        ! some preprations


        ! photon transfer
        print *,'| | calculate temperature with Monte Carlo method'
        print *,'| | | starting photon transfer ... [this may take a while]'
        !--------------------------------------------------------------------------!
        grid%t_dust(:,:)  = basics%t_dust_min
        grid%t_dust(0,:)  = 0.0_r2
        k_phot = model%n_star_emi/100

        ! initialize random number generator
        seed = -1
        !$omp parallel num_threads(basics%num_core) PRIVATE(rand_nr) FIRSTPRIVATE(seed)
        !$ seed = (omp_get_thread_num()+1)*(-1)
!~         !$ print *,seed
        CALL InitRandgen(rand_nr,seed,'RAN2')
        !--------------------------------------------------------------------------! 
!~         DO lucy = 1,5 ! TbD, in future we will use lucy iterations again, 
                         !      -> parallelisation should be more efficient
        !$omp do schedule(dynamic) PRIVATE(i_phot) FIRSTPRIVATE(photon)
        DO i_phot=1,model%n_star_emi
            ! show progress
            IF (modulo(i_phot, k_phot) == 0 .or. i_phot==model%n_star_emi) THEN
                    write (*,'(A,I3,A)') " | | | - progress : ", &
                    int(i_phot/real(model%n_star_emi)*100.0), ' % done...'//char(27)//'[A'
            END IF


            ! initiate photon
            CALL InitPhoton(photon, 1, 'Photon',dust%n_dust)
            photon%observe = .True.
            ! 1. start photon (from primary source only)
            CALL start_photon(basics,grid, model, rand_nr, fluxes, dust, photon, sources_in)

            ! 2. determine & go to next point of interaction
            CALL next_pos_const(model, rand_nr, grid, dust, photon)

            ! 3. transfer through model space
            !    if still inside model space: interaction
            !                                 + go to next position
            !                           else: observe photon leaving the model space
            photon%n_interact = 0
            DO
                IF (photon%inside .and. (photon%n_interact < n_interact_max)) THEN
                    photon%n_interact = photon%n_interact +1

                    CALL interact(basics, grid, dust, rand_nr, fluxes, photon)
                    CALL next_pos_const(model, rand_nr, grid, dust, photon)
                    cycle
                ELSE 
                    IF (photon%inside) THEN
                        kill_photon_count = kill_photon_count +1 
                    ELSEIF (.not. photon%inside .and. photon%observe) THEN
                            ! observe photon
!~                         print *, photon%inside
!~                         print *, kill_photon_count    
                        CALL observe_photon(basics, fluxes, model, grid, photon)
                    END IF
                    EXIT
                END IF
            END DO
            CALL ClosePhoton(photon)
        END DO
        !$omp end do nowait
        !$omp end parallel
        
        print *, '| | | photon transfer finished                '
        ! prepare & save final results
        DO i_dust=1,dust%n_dust
            grid%cell_energy(i_dust,:) = grid%cell_energy_sum(i_dust,:,1)/ &
                                         (basics%PIx4 * grid%cell_vol(:))
            grid%cell_energy_sum(i_dust,:,:) = 0.0_r2
        END DO
        CALL temp_final(basics, grid, dust)

    END SUBROUTINE primary_temp
    
    SUBROUTINE observe_photon(basics, fluxes, model, grid, photon)
        ! This routine should be moved to an own module (e.g. observe_mod)
        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                      :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
!~         TYPE(Dust_TYP),INTENT(IN)                        :: dust
!~         TYPE(SOURCES),INTENT(IN)                         :: sources_in
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
!~         TYPE(Randgen_TYP)                                :: rand_nr
        
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_map
        REAL(kind=r2), DIMENSION(3)                      :: ex, ey
        REAL(kind=r2)                                    :: angle
!~         INTEGER                                          :: x, y
        ! TbD: this should be done in advance and we should write observation
        !      variables in an own type
        i_map = 1
        grid%dir_xyz(1) = sin(model%th_map(i_map)) * sin(basics%PI2-model%ph_map(i_map))
        grid%dir_xyz(2) = sin(model%th_map(i_map)) * sin(model%ph_map(i_map))
        grid%dir_xyz(3) = sin(basics%PI2-model%th_map(i_map))
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
!~         rho_size_i = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)   
!~         rho_size_j = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)
!~         
!~         pix_res_i  = rho_size_i/model%distance *PI /(3600.0_r2*180.0_r2)
!~         pix_res_j  = rho_size_j/model%distance *PI /(3600.0_r2*180.0_r2)
        angle = acos(dot_product(photon%dir_xyz,grid%dir_xyz))/PI*180
        IF ( angle .lt. 5) THEN
            print *, angle
            print *, dot_product(photon%pos_xyz_li,ex)
            print *, dot_product(photon%pos_xyz_li,ey)
            
        END IF
        
    END SUBROUTINE observe_photon
  
      
END MODULE temp_mod

