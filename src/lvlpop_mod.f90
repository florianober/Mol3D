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

MODULE lvlpop_mod

    USE datatype
    USE var_global

    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    
    USE grd_mod
    USE math_mod
    USE fileio, ONLY : vis_T_exc_fits
    USE model_mod, ONLY : get_J_ext

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE :: pop_LTE, pop_FEP, pop_LVG, elem_LVG,                            &
               create_matrix, calc_collision_parameter
    !--------------------------------------------------------------------------!
    PUBLIC :: calc_lvlpop
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE calc_lvlpop(basics, grid, model, gas)

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Gas_TYP),INTENT(IN)                         :: gas
    !--------------------------------------------------------------------------!

    SELECT CASE(GetGasType(gas))
    CASE(1)
        print *, '| | using LTE method'
        CALL pop_LTE(grid, gas)
    CASE(2)
        print *, '| | using FEP method'
        CALL pop_FEP(basics, grid, gas)
    CASE(3)
        print *, '| | using LVG method'
        CALL pop_LVG(basics, grid, model, gas)
    CASE DEFAULT
        print *, 'ERROR: in calc_lvlpop routine (lvlpop_mod.f90)!'
        print *, 'ERROR: the selected method is not implemended'
        stop
    END SELECT
    IF (basics%do_save_T_exc) THEN
        print *, '| | saving excitation temperatures'
        CALL vis_T_exc_fits(basics, grid, model, gas)
    END IF
    END SUBROUTINE calc_lvlpop

    !--------------------------------------------------------------------------!
    ! Routines to calculate level populations:
    !--------------------------------------------------------------------------!

    SUBROUTINE pop_LTE(grid, gas)

    ! ---    
    ! assumtion: T_dust = T_kin = T_gas 
    ! usable und easy assumption for optical thick molecules like CO

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Gas_TYP),INTENT(IN)                         :: gas

    INTEGER                                          :: i_cell
    !--------------------------------------------------------------------------!

    DO i_cell = 1, grid%n_cell
        IF ( grid%t_dust(i_cell,1) .lt. 1.0e-20 ) THEN
            grid%lvl_pop(:,i_cell) = 0.0
        ELSE
            grid%lvl_pop(:,i_cell) =  gas%g_level(:)/gas%g_level(1) *          &
                                      exp(-con_h/(con_k*grid%t_gas(i_cell))*&
                                      (gas%energylevel(:))*con_c*100.0_r2)
                                    
            grid%lvl_pop(:,i_cell) = grid%lvl_pop(:,i_cell) /                  &
                                     sum(grid%lvl_pop(:,i_cell))!normalize array
        END IF
    END DO

    END SUBROUTINE pop_LTE


    SUBROUTINE pop_FEP(basics, grid, gas)

    ! ---    free escape probability
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                                 :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                               :: grid
    TYPE(Gas_TYP),INTENT(IN)                                   :: gas
    !--------------------------------------------------------------------------!
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl)       :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                     :: c
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                     :: new_pop
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl)                   :: J_mid
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl)                   :: J_ext
    REAL(kind=r1)                                              :: sum_p
    REAL(kind=r2),DIMENSION(1:gas%col_trans,1:2)               :: final_col_para
    INTEGER                                                    :: i_cell, k
    !--------------------------------------------------------------------------!
    new_pop    = 1.0e-30
    new_pop(1) = 1.0
    ! get external radiation field -> definition in model_mod.f90
    J_ext(:) = get_J_ext(gas%trans_freq(:))

    J_mid = J_ext
    k = grid%n_cell/100
    !$omp parallel num_threads(basics%num_core)
    !$omp do schedule(dynamic) private(i_cell,new_pop,final_col_para,sum_p, A,c)
    DO i_cell = 1, grid%n_cell
        IF (modulo(i_cell, k) == 0 .or. i_cell == grid%n_cell) THEN
            WRITE (*,'(A,I3,A)') ' | | | ',int(i_cell/real(grid%n_cell)*100.0),&
                                 ' % done'//char(27)//'[A'
        END IF

        !IF (grid%grd_mol_density(i_cell) .lt. 1.0e-200_r2) CYCLE
        final_col_para = calc_collision_parameter(grid,gas,i_cell)
        CALL create_matrix(gas, J_mid, A, c, final_col_para)

        CALL solve_eq(A,gas%egy_lvl,c,new_pop)

        ! verify the result (sum over all populations should be 1)

        sum_p = sum(new_pop(:))

        IF ( abs(1.0-sum_p) > 1.0e4*epsilon(sum_p) ) THEN
            print *, 'Warning, level polulation do not sum to 1'
            print *, abs(1.0-sum_p)
        END IF

        grid%lvl_pop(:,i_cell) = new_pop
    END DO
    !$omp end do nowait
    !$omp end parallel
    
    END SUBROUTINE pop_FEP
    
    
    SUBROUTINE pop_LVG(basics, grid, model, gas)
    !
    ! In the LVG (for a keplerian disk) approach we estimate the
    ! coherence length L at each cells mitpoint. Here, we assume that the disk 
    ! is rotating in the equatorial disk plane only!
    ! Maybe we can generalize this in future
    !
    ! ---    
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Model_TYP),INTENT(IN)                       :: model
    TYPE(Gas_TYP),INTENT(IN)                         :: gas
    !--------------------------------------------------------------------------!
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl)       :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                     :: c
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                     :: old_pop
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl)                     :: new_pop
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl)                   :: J_mid, J_ext
    REAL(kind=r1)                                              :: sum_p
    REAL(kind=r2)                                              :: R_mid, L
    REAL(kind=r2),DIMENSION(1:gas%col_trans,1:2)               :: final_col_para

    INTEGER                                       :: i_cell, max_iteration
    INTEGER                                       :: i, j, k
    !--------------------------------------------------------------------------!

    max_iteration = 100
    ! get external radiation field -> definition in model_mod.f90
    J_ext(:) = get_J_ext(gas%trans_freq(:))

    k = grid%n_cell/100
    !$omp parallel num_threads(basics%num_core)
    !$omp do schedule(dynamic) private(i_cell, new_pop, old_pop, J_mid,        &
    !$omp                              final_col_para, R_mid, L, sum_p, i, A, c)
    
    DO i_cell=1, grid%n_cell
        IF (modulo(i_cell,k) == 0 .or. i_cell == grid%n_cell) THEN
            WRITE (*,'(A,I3,A)') ' | | | ',                                    &
                      int(i_cell/real(grid%n_cell)*100.0),                     &
                      ' % done'//char(27)//'[A'
        END IF
        !starting values
        new_pop(:) = 1e-30
        old_pop(:) = 1e-30
        old_pop(1) = 1.0
        new_pop(1) = 1.0
        A(:,:) = 0.0
        c = 0.0

        IF (grid%cell_gauss_a(i_cell) * grid%absvelo(i_cell) < 1.0e-16_r2) CYCLE
        R_mid = sqrt(grid%cellmidcaco(i_cell, 1)**2 +                          &
                     grid%cellmidcaco(i_cell, 2)**2)

        L = R_mid * sqrt(2.0_r2/3.0_r2 /( grid%cell_gauss_a(i_cell) *          &
                            grid%absvelo(i_cell))) * model%ref_unit
        final_col_para = calc_collision_parameter(grid, gas, i_cell)

        ! Now iterate
        DO i = 1, max_iteration
            old_pop = new_pop
            !------------------------------------------------------------------!
            J_mid = elem_LVG(grid%grd_mol_density(i_cell),                     &
                             new_pop(gas%trans_upper(:)),                      &
                             new_pop(gas%trans_lower(:)),                      &
                             grid%cell_gauss_a(i_cell),                        &
                             gas%trans_einstA(:),                              &
                             gas%trans_einstB_u(:),                            &
                             gas%trans_einstB_l(:),                            &
                             L,                                                &
                             J_ext(:))
            !------------------------------------------------------------------!

            CALL create_matrix(gas, J_mid, A, c, final_col_para)
            
            CALL solve_eq(A, gas%egy_lvl, c, new_pop)
            
            ! verify the result (sum over all populations should be 1)

            sum_p = sum(new_pop(:))
            IF ( abs(1.0-sum_p) > 1.0e4 * epsilon(sum_p) ) THEN
                print *, 'Warning, level polulation do not sum to 1'
                print *, abs(1.0-sum_p)
            ELSE IF (sum_p /= sum_p) THEN
                print *, 'ERROR: lvl_populations are nan', i_cell, i
                print *, old_pop
                print *, 'news pop'
                print *, new_pop
                print *, 'J_mid'
                print *, J_mid
                print *, ''
                print *, L, R_mid
                print *, ''
                print *, grid%cell_gauss_a(i_cell) * grid%absvelo(i_cell)
                print *, grid%cell_gauss_a(i_cell) , grid%absvelo(i_cell)
                stop
            END IF
            j = MAXLOC(new_pop, 1)
            IF ( abs(new_pop(j)-old_pop(j)) / (old_pop(j)+EPSILON(old_pop(j))) &
                        .lt. 1.0e-2 .and. i > 1) EXIT
        END DO
        IF (i == max_iteration +1 .and. show_error) THEN
            print *, 'Warning: maximum iteration needed in cell:', i_cell
            print *, i_cell, i
            print *, ''
            print *, abs(new_pop(MAXLOC(new_pop, 1))-old_pop(MAXLOC(new_pop, 1)))/ &
                        (old_pop(MAXLOC(new_pop, 1))+EPSILON(old_pop))
            print *, ''
        END IF
        grid%lvl_pop(:,i_cell) = new_pop
    END DO

    !$omp end do nowait
    !$omp end parallel

    END SUBROUTINE pop_LVG

    ELEMENTAL FUNCTION elem_LVG(mol_dens, lvl_pop_u, lvl_pop_l, doppler,       &
                                A, B_u, B_l, L, J_ext)                         &
                                RESULT(J_mid)
        ! ---    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2),INTENT(IN)                                 :: mol_dens
        REAL(kind=r1),INTENT(IN)                                 :: lvl_pop_u
        REAL(kind=r1),INTENT(IN)                                 :: lvl_pop_l
        REAL(kind=r1),INTENT(IN)                                 :: doppler

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

        !----------------------------------------------------------------------!
        j       =  mol_dens * lvl_pop_u * A * linescale*doppler
        alpha   =  mol_dens *(lvl_pop_l*B_l-lvl_pop_u*B_u)*linescale*doppler

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

        J_mid = (1.0 - beta) * S + beta * J_ext
    END FUNCTION elem_LVG


    PURE SUBROUTINE create_matrix(gas, J_mid, A, c, final_col_para)
    ! ---    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Gas_TYP),INTENT(IN)                                         :: gas
    !--------------------------------------------------------------------------!
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl,1:gas%egy_lvl),INTENT(OUT) :: A
    REAL(kind=r1),DIMENSION(1:gas%egy_lvl),INTENT(OUT)               :: c
    REAL(kind=r2),DIMENSION(1:gas%trans_lvl),INTENT(IN)              :: J_mid
    REAL(kind=r2),DIMENSION(1:gas%col_trans,1:2),INTENT(IN)   :: final_col_para

    INTEGER                                                          :: col_tr
    INTEGER                                                          :: i, j, k
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
            ELSE IF (gas%trans_lower(j) == i) THEN
                k = gas%trans_upper(j)
                A(i,k) = A(i,k) + gas%trans_einstA(j) + gas%trans_einstB_u(j)*J_mid(j)
                A(i,i) = A(i,i) - gas%trans_einstB_l(j)*J_mid(j)
            END IF
        END DO
        DO col_tr=1, gas%col_trans
            IF (gas%col_lower(col_tr,gas%col_id(1)) == i) THEN
                k = gas%col_upper(col_tr,gas%col_id(1))
                A(i,k) = A(i,k) + final_col_para(col_tr,1)
                A(i,i) = A(i,i) - final_col_para(col_tr,2)
            ELSE IF ( gas%col_upper(col_tr,gas%col_id(1)) == i) THEN
                k = gas%col_lower(col_tr,gas%col_id(1))
                A(i,k) = A(i,k) + final_col_para(col_tr,2)
                A(i,i) = A(i,i) - final_col_para(col_tr,1)
            END IF
        END DO
    END DO

    ! The first entry should be one
    ! to solve linear equation 
    c(1) = 1.0_r2
    A(1,:) = 1.0_r2

    END SUBROUTINE create_matrix

    FUNCTION calc_collision_parameter(grid,gas,i_cell) RESULT(final_col_para)
        ! ---
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Grid_TYP),INTENT(IN)                              :: grid
        TYPE(Gas_TYP),INTENT(IN)                               :: gas

        REAL(kind=r2),DIMENSION(1:gas%col_trans,1:2)           :: final_col_para
        REAL(kind=r2), DIMENSION(1:gas%col_trans)              :: col_mtr_tmp_ul
        REAL(kind=r2), DIMENSION(1:gas%col_trans)              :: col_mtr_tmp_lu


        INTEGER,INTENT(IN)                                     :: i_cell
        INTEGER                                                :: j, k, hi_i
        !----------------------------------------------------------------------!
        final_col_para(:, :)  = 0.0_r2
        DO j = 1, gas%col_partner
            k = gas%col_id(j)
            IF( any(gas%col_upper(:,k) == 0) ) THEN
                CONTINUE
            ELSE
                IF (grid%t_gas(i_cell) .lt. MINVAL(gas%col_alltemps(:, k)) ) THEN
                    hi_i = 1
                ELSE IF (grid%t_gas(i_cell) .gt. MAXVAL(gas%col_alltemps(:, k)) ) THEN
                    hi_i = gas%col_temps - 1
                ELSE

                    hi_i = binary_search(grid%t_gas(i_cell), gas%col_alltemps(:,k))
                END IF

                col_mtr_tmp_ul(:) = &
                                    ipol2(gas%col_alltemps(hi_i,k),            &
                                          gas%col_alltemps(hi_i+1,k),          &
                                          gas%col_colmatrix(k,:,hi_i),         &
                                          gas%col_colmatrix(k,:,hi_i+1),       &
                                          grid%t_gas(i_cell)) *                &
                                    grid%grd_col_density(i_cell, k)

                col_mtr_tmp_lu(:) =                                            &
                        col_mtr_tmp_ul(:) *                                    &
                        gas%g_level(gas%col_upper(:,k)) /                      &
                        gas%g_level(gas%col_lower(:,k)) *                      &
                        exp(-con_h *                                           &
                        (gas%energylevel(gas%col_upper(:,k))*con_c*100.0_r2 -  &
                        gas%energylevel(gas%col_lower(:,k))*con_c*100.0_r2) /  &
                        (con_k*grid%t_gas(i_cell)))

                final_col_para(:,1) = final_col_para(:,1) + col_mtr_tmp_ul(:)
                final_col_para(:,2) = final_col_para(:,2) + col_mtr_tmp_lu(:)
            END IF
        END DO

    END FUNCTION calc_collision_parameter

END MODULE lvlpop_mod
