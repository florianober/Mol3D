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
    TYPE(Fluxes_TYP) ,INTENT(IN)                     :: fluxes
    TYPE(SOURCES),INTENT(IN)                         :: sources_in
    !--------------------------------------------------------------------------!  
    REAL(kind=r2), DIMENSION(1:3)                  :: caco
    REAL(kind=r2), DIMENSION(1:3)                  :: moco
    REAL(kind=r2), DIMENSION(1:gas%col_trans)      :: col_mtr_tmp_ul
    REAL(kind=r2), DIMENSION(1:gas%col_trans)      :: col_mtr_tmp_lu
    
    INTEGER                                          :: i_cell, hi_i
    INTEGER                                          :: i_a
    INTEGER                                          :: i_b
    INTEGER                                          :: i_c
    INTEGER                                          :: j,k
    
    !--------------------------------------------------------------------------! 
    IF (.not. basics%old_model) THEN
        IF ( basics%calc_tmp) THEN
            IF (.not. SourcesInitialized(sources_in)) THEN
                print *, 'ERROR: There is no source defined'
                stop
            END IF
            CALL primary_temp(basics,grid,model,dust,sources_in,fluxes)
        ELSE IF (.not. basics%calc_tmp) THEN
            print *, '  using analytical temperature distribution'
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

    
    ! CALCULATE ALL PROBERTIES, WHICH NEED TEMPERATURE INFORMATIONS
    
    IF (basics%do_raytr) THEN
        PRINT *, '  temperature included, now calculate some additional parameters'
        ALLOCATE(grid%col_finalcolmatrixup(1:gas%col_trans,1:grid%n_cell), &
                 grid%col_finalcolmatrixlow(1:gas%col_trans,1:grid%n_cell) )
        grid%col_finalcolmatrixup  = 0.0_r2
        grid%col_finalcolmatrixlow = 0.0_r2
        DO i_cell = 1,grid%n_cell

            ! set line width in each cell (in fact the inverse value)
            grid%cell_gauss_a(i_cell) = 1.0_r2/(sqrt(2.0_r2*con_k*grid%t_dust(i_cell,1)/ &  
                                                (gas%mol_weight*1.0e-3_r2/con_Na) &
                                                +100.0_r2**2)) ! 100.0 eq trub. line width TbD!~!
                                                
            ! save quadratic line width value for faster calculations
            
            grid%cell_gauss_a2(i_cell)  =   grid%cell_gauss_a(i_cell)**2
            ! calculate (interpolate) collision parameters c = f(T)
            DO j = 1, gas%col_partner
                k = gas%col_id(j)
                IF( any(gas%col_upper(:,k) == 0) ) THEN
                    CONTINUE
                ELSE
                    IF (grid%t_gas(i_cell) .lt. MINVAL(gas%col_alltemps(:,k)) ) THEN
                    hi_i = 1
                    
                    ELSE IF (grid%t_gas(i_cell) .gt. MAXVAL(gas%col_alltemps(:,k)) ) THEN
                    
                    hi_i = gas%col_temps - 1
                    ELSE

                    hi_i = binary_search(REAL(grid%t_gas(i_cell),kind=r2), REAL(gas%col_alltemps(:,k),kind=r2))
                    END IF
                            
                    col_mtr_tmp_ul(:) = &
                                        ipol2(REAL(gas%col_alltemps(hi_i,k),kind=r2), REAL(gas%col_alltemps(hi_i+1,k),kind=r2),   &
                                        REAL(gas%col_colmatrix(k,:,hi_i),kind=r2),REAL(gas%col_colmatrix(k,:,hi_i+1),kind=r2), &
                                        REAL(grid%t_gas(i_cell),kind=r2))* grid%grd_col_density(i_cell,k)
                                        
                    col_mtr_tmp_lu(:) = &
                            col_mtr_tmp_ul(:)* &
                            gas%g_level(gas%col_upper(:,k))/gas%g_level(gas%col_lower(:,k)) * &
                            exp(-con_h * &
                            (gas%energylevel(gas%col_upper(:,k))*con_c*100.0_r2 - &
                            gas%energylevel(gas%col_lower(:,k))*con_c*100.0_r2)  &
                            /(con_k*grid%t_gas(i_cell)))
                            
                    grid%col_finalcolmatrixup(:,i_cell) =   grid%col_finalcolmatrixup(:,i_cell) + &
                                                            col_mtr_tmp_ul(:)
                    grid%col_finalcolmatrixlow(:,i_cell) =  grid%col_finalcolmatrixlow(:,i_cell) + &
                                                            col_mtr_tmp_lu(:)
                END IF 
            END DO
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
        TYPE(Fluxes_TYP),INTENT(IN)                      :: fluxes
        TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
        TYPE(Model_TYP),INTENT(IN)                       :: model
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        TYPE(SOURCES),INTENT(IN)                         :: sources_in
    !~     TYPE(Gas_TYP),INTENT(INOUT)                      :: gas
    
        TYPE(Randgen_TYP)                                 :: rand_nr
        
        !--------------------------------------------------------------------------!  

        INTEGER                                          :: i_lam, i_phot, i_dust
        INTEGER                                          :: seed, k_phot
!~         !$ INTEGER omp_get_thread_num 
                
        !# shared variables
        TYPE(PHOTON_TYP)                                 :: photon
        
        !--------------------------------------------------------------------------! 
        ! ---
        ! some preprations


        ! photon transfer
        print *,' calculate temperature via immediate reemission concept'
        print *,' starting photon transfer ... [this may take a while]'
        !--------------------------------------------------------------------------!
        grid%t_dust(:,:)  = basics%t_dust_min
        grid%t_dust(0,:)  = 0.0_r2


!~         !$ call omp_set_num_threads( basics%num_core )
        
!~         !$omp parallel PRIVATE(i_phot,seed,photon,rand_nr,t_dust,grd_d_l,i_star_abs)
        
        
        ! initialize random number generator
        seed = -1
!~         !$ seed = omp_get_thread_num()+1
        !print *, 'here'
        CALL InitRandgen(rand_nr,seed,'RAN2')
        !--------------------------------------------------------------------------! 
!~         DO lucy = 1,5 ! TbD, in future we will use lucy iterations again, 
                         !      -> parallelisation should be more efficient
        k_phot = model%n_star_emi/100
        DO i_phot=1,model%n_star_emi
            ! show progress
!~                 !$OMP MASTER
            IF (modulo(i_phot, k_phot) == 0 .or. i_phot==model%n_star_emi) THEN
                    write (*,'(A,I3,A)') "   - photon : ", int(i_phot/real(model%n_star_emi)*100.0), ' % done...'//char(27)//'[A'
            END IF
!~             !$OMP END MASTER
!~             !$omp do schedule(static)

            ! initiate photon
            CALL InitPhoton(photon, 1, 'Photon',dust%n_dust)
            
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
                    END IF
                    EXIT
                END IF
            END DO
            CALL ClosePhoton(photon)
!~                 !$omp end do nowait
        END DO
        print *, ' photon transfer finished     '
        ! prepare & save final results
           
        ! [solution 2] apply mean intensity method
        !              => higher accuracy: from now on this temperature distribution will be used
!~         !$OMP CRITICAL
!~             grid%t_dust = t_dust
!~         !$OMP END CRITICAL
!~         DEALLOCATE (t_dust)
!~         !$omp end parallel
        DO i_dust=1,dust%n_dust
            grid%cell_energy(i_dust,:) = grid%cell_energy_sum(i_dust,:,1)/ (basics%PIx4 * grid%cell_vol(:))
            grid%cell_energy_sum(i_dust,:,:) = 0.0_r2
        END DO
        CALL temp_final2(basics, grid, dust)
        !stop
        ! save SED (quick&dirty SED)
        !  call sv_stokes_sed()     
    END SUBROUTINE primary_temp
  

!~   ! ################################################################################################
!~   ! simulation of primary source radiation:  
!~   ! ---------------------------------------
!~   !       photon_type = 4: calculation of primary source scattered light map
!~   !                   5: calculation of primary source SED
!~   ! ---
!~   subroutine primary_scatt()
!~     use datatype
!~     use var_global
!~     use start_mod
!~     use interact_mod
!~     use observe_mod
!~     use sv_results_mod
!~     use transfer_mod
!~     use immediate_mod
!~     use math_mod
!~     use tools_mod
!~ 
!~     implicit none
!~ 
!~     integer :: i_lam_map, i_phot, i_phot_show 
!~     ! ---
!~     ! preparation
!~     i_phot_show = nint( n_star_emi/10.0_r2 )
!~ 
!~     ! photon transfer
!~     print *, "    photon transfer ... [this may take a while]"
!~ 
!~     do i_lam_map=1, n_lam_map
!~        print *,"   - wavelength           : ", i_lam_map, " / ", n_lam_map
!~ 
!~        do i_phot=1, nint(n_star_emi)
!~           if (photon_type==4) then
!~              if ( modulo(i_phot,i_phot_show)==0 ) then
!~                 print *,"     - photon counter [%] : ", &
!~                      nint( 100.0_r2 * real(i_phot,kind=r2) / n_star_emi )
!~              end if
!~           end if
!~           
!~           ! 1. start photon (from primary source only)
!~           call start_prim( num_lam_map(i_lam_map) )
!~ 
!~           ! 2. determine & go to next point of interaction
!~           call next_pos_const()
!~ 
!~           ! 3. transfer through model space
!~           !    if still inside model space: interaction
!~           !                                 + go to next position
!~           !                           else: observe photon leaving the model space
!~           do
!~              if (inside .and. (stokes(1)>i_min)) then
!~                 call interact()
!~                 call next_pos_const()
!~                 cycle
!~              else
!~                 call observe_MC(i_lam_map)
!~                 exit
!~              end if
!~           end do
!~        end do
!~     end do
!~ 
!~     ! ---
!~     ! prepare & save final results
!~     select case(photon_type)
!~    
!~     case(4)
!~        ! save scattered light maps
!~        call sv_stokes_map()
!~        
!~     case(5)
!~        ! save SED (scattered light SED only)
!~        call sv_stokes_sed()
!~ 
!~     case default
!~        print *,"<!> subroutine primary_scatt: wrong input parameter. stopped."
!~        call stop_mc3d()
!~     end select
!~        
!~   end subroutine primary_scatt
!~   
!~ 
!~   ! ################################################################################################
!~   ! dust reemission map / SED
!~   subroutine dust_RT()
!~     use var_global
!~     use reemission_mod
!~     use sv_results_mod
!~     use tools_mod
!~     
!~     implicit none
!~     ! ---
!~     ! 1. run simulation, depending on type of radiative transfer chosen
!~     select case(ree_type)
!~     case(1)
!~        ! mc radiative transfer (including scattering)
!~        call reemission_mc()
!~     case(2)
!~        ! raytracing (no scattering); argument: get_map = .true.
!~        call reemission_raytrace()
!~     case default
!~        print *,"<!> subroutine dust_RT [a]: wrong input parameter. stopped."
!~        call stop_mc3d()
!~     end select
!~     
!~     ! ---
!~     ! 2. save results
!~     select case(photon_type)
!~     case(2)
!~        call sv_stokes_map()
!~     case(3)
!~        call sv_stokes_sed()
!~     case default
!~        print *,"<!> subroutine dust_RT [b]: wrong input parameter. stopped."
!~        call stop_mc3d()
!~     end select
!~        
!~   end subroutine dust_RT
      
END MODULE temp_mod

