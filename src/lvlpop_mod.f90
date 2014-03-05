MODULE lvlpop_mod

    USE datatype
    USE var_globalnew

    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    
    USE grd_mod
    USE math_mod
    USE tools_mod
    USE transfer_mod
    USE string_mod
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE :: pop_LTE, pop_FEP, pop_LVG, elem_LVG, create_matrix!, elem_LVG
    !--------------------------------------------------------------------------!
    PUBLIC :: calc_lvlpop
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE calc_lvlpop(basics,grid , model, gas)
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Gas_TYP),INTENT(IN)                         :: gas
    
    INTEGER                                          :: i
    
    CHARACTER(len=252)                               :: filename
    CHARACTER(len=252)                               :: fmtstr
    CHARACTER(len=4)                                 :: fileext
    CHARACTER(len=256)                               :: outname
    
    !--------------------------------------------------------------------------!  
!~     CALL cpu_time(t1)
    SELECT CASE(GetGasType(gas))
    CASE(1)
        print *, ' using LTE method'
        CALL pop_LTE(grid, gas)
    CASE(2)
        print *, ' using FEP method'
        CALL pop_FEP(grid, gas)
    CASE(3)
        print *, ' using LVG method'
        CALL pop_LVG(basics, grid, model, gas)   
            
    CASE DEFAULT
        print *, 'ERROR: in calc_lvlpop!'
        print *, 'ERROR: selected method is not implemended yet'
        stop
    END SELECT
    DEALLOCATE(grid%col_finalcolmatrixup,grid%col_finalcolmatrixlow)
    
    print *, 'populations calculated'
    fileext = '.dat'
    filename = TRIM(basics%path_results)//Getproname(basics)//'_lvl_pop'
    outname = TRIM(filename)//fileext
    OPEN(unit=1, file=TRIM(outname), &
        action="write", status="unknown", form="formatted")
    WRITE(fmtstr,'(A,I0,A)') '(I8,',gas%egy_lvl,'(ES15.6E3))'

    DO i = 1, grid%n_cell

        WRITE(unit=1,fmt=fmtstr) i, grid%lvl_pop(i,:)
    
    END DO
    CLOSE(unit=1)
!~     CALL cpu_time(t2)
!~     write (*,'(a,1pg12.4)')    'cpu_time:     ', t2-t1
    END SUBROUTINE calc_lvlpop

    
    !--------------------------------------------------------------------------!
    ! Routines to calculate level populations:
    !--------------------------------------------------------------------------!
    
    SUBROUTINE pop_LTE(grid, gas)
    
    ! ---    
    ! assumtion: T_dust = T_kin = T_gas 
    ! good for optical thick molecules like CO
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Gas_TYP),INTENT(IN)                         :: gas
    
    INTEGER                                          :: i_cell
    !--------------------------------------------------------------------------!        
        
    DO i_cell = 1, grid%n_cell
        IF ( grid%t_dust(i_cell,1) .lt. 1.0e-20 ) THEN
            grid%lvl_pop(i_cell,:) = 0.0
        ELSE
            grid%lvl_pop(i_cell,:) =  gas%g_level(:)/gas%g_level(1) *exp(-con_h/(con_k*grid%t_dust(i_cell,1)) &
                                    *(gas%energylevel(:))*con_c*100.0_r2)
                                    
            grid%lvl_pop(i_cell,:) = grid%lvl_pop(i_cell,:)/sum(grid%lvl_pop(i_cell,:)) !normalize array
        END IF
    END DO
    
    END SUBROUTINE pop_LTE
    
    
    SUBROUTINE pop_FEP(grid,gas)
    
    ! ---    free escape probability
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(INOUT)                                   :: grid
    TYPE(Gas_TYP),INTENT(IN)                                       :: gas
    
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl)           :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                         :: c
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl)                       :: J_mid, J_ext
    REAL(kind=r1)                                                  :: sum_p
    INTEGER                                                        :: i_cell, i
    !--------------------------------------------------------------------------!        
    
    DO i=1, gas%trans_lvl
!~             J_ext(i) = 1e-14_r2*planckhz(1000.0_r2,gas%trans_freq(i)) + planckhz(2.72_r2,gas%trans_freq(i))
            J_ext(i) = planckhz(2.75_r1,gas%trans_freq(i))
    END DO
    J_mid = J_ext
    
    DO i_cell = 1, grid%n_cell
        IF (grid%grd_mol_density(i_cell) .lt. 1.0e-200_r2) CYCLE
        CALL create_matrix(grid,gas,J_mid,i_cell,A,c)

        CALL solvestateq(A,gas%egy_lvl,c,grid%lvl_pop(i_cell,:))
        
        
        ! verify the result (sum over all populations should be 1)

        sum_p = sum(grid%lvl_pop(i_cell,:))

        IF ( abs(1.0-sum_p) > 1.0e4*epsilon(sum_p) ) THEN
            print *, 'Warning, level polulation do not sum to 1'
            print *, abs(1.0-sum_p)
        END IF
    END DO

    END SUBROUTINE pop_FEP
    
    
    SUBROUTINE pop_LVG(basics, grid, model, gas)
    !
    ! In the LVG (for a keplerian disk) approach we estimate the
    ! coherence length L at each cells mitpoint. Here, we assume that the disk is rotating in the
    ! equatorial disk plane only!
    ! Maybe we can generalize this in future
    !
    ! ---    
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Gas_TYP),INTENT(IN)                         :: gas
    
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl)           :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                         :: c
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                         :: old_pop
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl)                       :: J_mid, J_ext
    REAL(kind=r1)                                                  :: sum_p
    REAL(kind=r2)                                                  :: R_mid, L
    
    INTEGER                                           :: i_cell, i, max_iteration
    !--------------------------------------------------------------------------!        

!~     print *, 'ERROR: method under development'
    
    max_iteration = 200
    !calculate external radiation field

    J_ext(:) = planckhz(2.75,gas%trans_freq(:))

   
    DO i_cell = 1,grid%n_cell
        
        ! Now iterate
    

        IF (grid%grd_mol_density(i_cell) .lt. 1.0e-40_r2) CYCLE
        
                
!~         DO i = 1, gas%egy_lvl
!~             grid%lvl_pop(i_cell,i) =  gas%g_level(i)/gas%g_level(1) *exp(-con_h/(con_k*grid%t_dust(i_cell,1)) &
!~                                 *(gas%energylevel(i))*con_c*100.0_r2)
!~                                 
!~         END DO
!~         grid%lvl_pop(i_cell,:) = grid%lvl_pop(i_cell,:)/sum(grid%lvl_pop(i_cell,:)) !normalize array
        
        
        R_mid = sqrt(grid%cellmidcaco(i_cell,1)**2+grid%cellmidcaco(i_cell,2)**2)
        L = R_mid * sqrt(2.0/3.0 /( grid%cell_gauss_a(i_cell) * grid%absvelo(i_cell))) * model%ref_unit
        
        DO i = 1, max_iteration  ! more iteration needed if cells are far from LTE solution
            old_pop = grid%lvl_pop(i_cell,:)

            !---------------------------------------------------------------------------------------
            J_mid = elem_LVG(grid%grd_mol_density(i_cell),grid%lvl_pop(i_cell,gas%trans_upper(:)),  & 
                             grid%lvl_pop(i_cell,gas%trans_lower(:)),grid%cell_gauss_a(i_cell),     &
                             basics%linescale,gas%trans_einstA(:),gas%trans_einstB_u(:),         &
                             gas%trans_einstB_l(:), L, J_ext(:))
            !---------------------------------------------------------------------------------------
            
            CALL create_matrix(grid,gas,J_mid,i_cell,A,c)
            
            CALL solvestateq(A,gas%egy_lvl,c,grid%lvl_pop(i_cell,:))
            
            ! verify the result (sum over all populations should be 1)

            sum_p = sum(grid%lvl_pop(i_cell,:))
            IF ( abs(1.0-sum_p) > 1.0e4*epsilon(sum_p) ) THEN
                print *, 'Warning, level polulation do not sum to 1'
                print *, abs(1.0-sum_p)
            ELSE IF (isnan(sum_p)) THEN
                print *, 'lvl_populations are nan',i_cell, i
                print *, old_pop
                print *, 'J_mid'
                print *, J_mid
                stop
            END IF
            IF ( MAXVAL(abs(grid%lvl_pop(i_cell,:)-old_pop)/old_pop) .lt. 1.0e-2) THEN
                EXIT

            END IF
        END DO
        IF (i == max_iteration +1 ) print *, 'Warning: maximum iteration needed in cell:', i_cell
!~         print *, 'iter:',  i, i_cell
!~             print *,MAXVAL(abs(grid%lvl_pop(i_cell,:)-old_pop)/old_pop)
    END DO
    END SUBROUTINE pop_LVG
    
    ELEMENTAL FUNCTION elem_LVG(mol_dens,lvl_pop_u,lvl_pop_l,doppler,line_konst,A,B_u,B_l, L,J_ext) &
                               RESULT(J_mid)
        ! ---    
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        REAL(kind=r2),INTENT(IN)                                 :: mol_dens
        REAL(kind=r1),INTENT(IN)                                 :: lvl_pop_u
        REAL(kind=r1),INTENT(IN)                                 :: lvl_pop_l
        REAL(kind=r1),INTENT(IN)                                 :: doppler
        REAL(kind=r1),INTENT(IN)                                 :: line_konst
        REAL(kind=r1),INTENT(IN)                                 :: A
        REAL(kind=r1),INTENT(IN)                                 :: B_u
        REAL(kind=r1),INTENT(IN)                                 :: B_l
        REAL(kind=r2),INTENT(IN)                                 :: L
        REAL(kind=r2),INTENT(IN)                                 :: J_ext
        
        REAL(kind=r1)                                            :: j
        REAL(kind=r1)                                            :: alpha
        REAL(kind=r1)                                            :: beta
        REAL(kind=r1)                                            :: S
        REAL(kind=r1)                                            :: tau
        REAL(kind=r2)                                            :: J_mid

        
        !--------------------------------------------------------------------------!
        j       =  mol_dens * lvl_pop_u * A * line_konst*doppler
        alpha   =  mol_dens *(lvl_pop_l*B_l-lvl_pop_u*B_u)*line_konst*doppler
                
!~         IF ( alpha .lt. 1e-20 .or. isnan(alpha) ) THEN 
        IF ( alpha .lt. 1e-20  ) THEN 
            S     = 0.0_r2
            alpha = 0.0
        ELSE
            S = j/alpha
        END IF
                    
        tau = alpha * L
                    
        IF (tau .lt. 1.0e-6) THEN
            beta = 1.0 - 0.5*tau
        ELSE
            beta = (1.0 - exp(-tau))/tau
        END IF

        J_mid = (1.0 - beta)*S+beta*J_ext

    END FUNCTION elem_LVG
    
    
    PURE SUBROUTINE create_matrix(grid,gas,J_mid,i_cell,A,c)
    ! ---    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(IN)                                          :: grid
    TYPE(Gas_TYP),INTENT(IN)                                           :: gas
    
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl),INTENT(OUT)   :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl),INTENT(OUT)                 :: c
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl),INTENT(IN)                :: J_mid

    INTEGER,INTENT(IN)                                                 :: i_cell
    INTEGER                                                            :: col_tr
    INTEGER                                                            :: i, j, k
    !--------------------------------------------------------------------------!    
    A(:,:) = 0.0
    c(:)   = 0.0
    DO i=1, gas%egy_lvl
        DO j=1, gas%trans_lvl 
            ! find line and fill matrix
            !
            IF (gas%trans_upper(j) == i) THEN
                k = gas%trans_lower(j)
                A(i,k) =  A(i,k) + gas%trans_einstB_l(j)*J_mid(j)
                A(i,i) =  A(i,i) - gas%trans_einstA(j)-gas%trans_einstB_u(j)*J_mid(j)
                !PRINT *, - gas%trans_einstA(j)-gas%trans_einstB_u(j)*J_mid(j)
            ELSE IF (gas%trans_lower(j) == i) THEN
                k = gas%trans_upper(j)
                A(i,k) = A(i,k) + gas%trans_einstA(j) + gas%trans_einstB_u(j)*J_mid(j)
                A(i,i) = A(i,i) - gas%trans_einstB_l(j)*J_mid(j)
            END IF
        END DO
!~         DO j = 1,gas%col_partner
            DO col_tr=1, gas%col_trans
                IF (gas%col_lower(col_tr,gas%col_id(1)) == i) THEN
                    k = gas%col_upper(col_tr,gas%col_id(1))
                    A(i,k) = A(i,k) + grid%col_finalcolmatrixup(col_tr,i_cell)
                    A(i,i) = A(i,i) - grid%col_finalcolmatrixlow(col_tr,i_cell)
                ELSE IF ( gas%col_upper(col_tr,gas%col_id(1)) == i) THEN
                    k = gas%col_lower(col_tr,gas%col_id(1))
                    A(i,k) = A(i,k) + grid%col_finalcolmatrixlow(col_tr,i_cell)
                    A(i,i) = A(i,i) - grid%col_finalcolmatrixup(col_tr,i_cell)
                END IF
            END DO
!~         END DO
    END DO
    

    ! The first(last?, first gives better/cleaner results) entry should be one to solve linear equation 
    c(1) = 1.0_r2
!~     c(gas%egy_lvl) = 1.0_r2

!~     A(gas%egy_lvl,:) = 1.0_r2
    A(1,:) = 1.0_r2

    END SUBROUTINE create_matrix
    
    
END MODULE lvlpop_mod
