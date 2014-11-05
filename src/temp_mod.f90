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
        INTEGER                                          :: i_dust, i_lam
        INTEGER                                          :: seed
        INTEGER(8)                                       :: i_phot, k_phot
        REAL(kind=r2),dimension(1:dust%n_lam)            :: sed
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
        k_phot = model%no_photon/100

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
        DO i_phot=1,model%no_photon
            ! show progress
            IF (modulo(i_phot, k_phot) == 0 .or. i_phot==model%no_photon) THEN
                    write (*,'(A,I3,A)') " | | | - progress : ", &
                    int(i_phot/real(model%no_photon)*100.0), ' % done...'//char(27)//'[A'
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
                        CALL observe_photon(basics, fluxes, model, dust,grid, photon)
!~                         print *, i_phot
                    END IF
                    EXIT
                END IF
            END DO
            CALL ClosePhoton(photon)
        END DO
        !$omp end do nowait
        !$omp end parallel
!~         fluxes%continuum_map_temp = fluxes%continuum_map_temp*1.0e26_r2/(model%distance*con_pc)**2
!~         print *, sum(fluxes%continuum_map_temp(:,:,:))
        print * , sum(fluxes%continuum_map_temp(:,:,:))
!~         DO i_lam = 1, dust%n_lam
!~             fluxes%continuum_map_temp(:,:,i_lam) = fluxes%continuum_map_temp(:,:,i_lam)/dust%d_lam(i_lam)
!~             sed(i_lam) = sum(fluxes%continuum_map_temp(:,:,i_lam))
!~         END DO
!~         print *, integ1(dust%lam(:), sed(:), 1, dust%n_lam )
        
        CALL save_continuum_map(model, basics, dust, fluxes, 1)
        
        print *, '| | | photon transfer finished                '
        stop
        ! prepare & save final results
        DO i_dust=1,dust%n_dust
            grid%cell_energy(i_dust,:) = grid%cell_energy_sum(i_dust,:,1)/ &
                                         (basics%PIx4 * grid%cell_vol(:))
            grid%cell_energy_sum(i_dust,:,:) = 0.0_r2
        END DO
        CALL temp_final(basics, grid, dust)

    END SUBROUTINE primary_temp
    
    SUBROUTINE observe_photon(basics, fluxes, model, dust, grid, photon)
        ! This routine should be moved to an own module (e.g. observe_mod)
        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Fluxes_TYP),INTENT(INOUT)                      :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
!~         TYPE(SOURCES),INTENT(IN)                         :: sources_in
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
!~         TYPE(Randgen_TYP)                                :: rand_nr
        
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_map
        REAL(kind=r2), DIMENSION(3)                      :: ex, ey
        INTEGER                                          :: x, y
        REAL(kind=r2)                                    :: angle, al
        REAL(kind=r2)                                    :: v_product
        REAL(kind=r2)                                    :: unit_value
        REAL(kind=r2)                                    :: rho_size_x, rho_size_y
        REAL(kind=r2)                                    :: pix_res_x, pix_res_y
        


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
        
        !
        !user unit, but should be fixed to [AU]
!~         rho_size_x = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)   
!~         rho_size_y = model%r_ou/(REAL(model%n_bin_map,KIND=r2)+0.5_r2)/model%zoom_map(1)
        rho_size_x = model%r_ou/(REAL(model%n_bin_map,KIND=r2))/model%zoom_map(1)   
        rho_size_y = model%r_ou/(REAL(model%n_bin_map,KIND=r2))/model%zoom_map(1)
!~         
        pix_res_x  = rho_size_x/model%distance *PI /(3600.0_r2*180.0_r2)
        pix_res_y  = rho_size_y/model%distance *PI /(3600.0_r2*180.0_r2)
!~         unit_value = 1.0e26_r2*pix_res_x*pix_res_y
        ! calculate the angle
!~         angle = dot_product(photon%dir_xyz,grid%dir_xyz)
        v_product = dot_product(photon%dir_xyz,grid%dir_xyz)
        if (v_product .gt. 0.0_r2 ) THEN
        angle = acos(v_product)*180.0_r2/PI
!~         ELSE
!~             angle = acos(v_product)*180.0_r2/PI+180.0_r2
!~         END IF
        ! cos(5.0_r2/180.0_r2*PI should be calculated before
!~         IF ( angle .lt. cos(5.0_r2/180.0_r2*PI)) THEN
!~         print *, angle
        al = 45.0_r2
!~         al = 90.0_r2/180.0_r2*PI
!~         if (angle -abs( arcdis( photon%dir_xyz(3),sqrt( photon%dir_xyz(1)**2 + photon%dir_xyz(2)**2), &
!~                            atanx(photon%dir_xyz(2), photon%dir_xyz(1) ), &
!~                            cos(model%th_map(i_map)),sin(model%th_map(i_map)), model%ph_map(i_map) ) )*180.0/PI &
!~                             .gt. 1.0) THEN

!~         END IF
        IF ( angle .le. al ) THEN
!~         IF ( angle .lt. 0.0873 ) THEN
!~             print *, angle
!~             print *, acos(dot_product(photon%dir_xyz,grid%dir_xyz))*180/PI
!~             
            ! now find position on map
!~             print *, dot_product(photon%pos_xyz_li,ex)
!~             print *, dot_product(photon%pos_xyz_li,ey)
!~             print *, photon%pos_xyz_li
            x = dot_product(photon%pos_xyz_li,ex)/rho_size_x + REAL(model%n_bin_map+1,KIND=r2)
            y = dot_product(photon%pos_xyz_li,ey)/rho_size_y + REAL(model%n_bin_map+1,KIND=r2)
!~             unit_value = dust%lam(photon%nr_lam)**2/con_c
!~             print *,angle
!~             unit_value = (1.0_r2-cos(al*PI/360.0))/(PI*4.0_r2)
            unit_value = 1.0_r2/((1.0_r2-cos(al*PI/180.0))*(PI*2.0_r2))
!~             print *, unit_value
!~             unit_value = 1/(2*PI*(1.0-angle))
!~             print *, unit_value
!~             stop
            fluxes%continuum_map_temp(x,y,photon%nr_lam) = fluxes%continuum_map_temp(x,y,photon%nr_lam) + &
                                                           photon%energy*unit_value
!~             print *, photon%energy
!~             print *, unit_value
!~             print *, dust%lam(photon%nr_lam)
!~             print *, unit_value * photon%energy
!~ 
!~             stop
        END IF
        
    END SUBROUTINE observe_photon
  
      
END MODULE temp_mod

