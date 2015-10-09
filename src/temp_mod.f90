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
    USE math_mod, ONLY : binary_search, ipol2

    USE fileio


    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: set_temperature
    
CONTAINS

    SUBROUTINE set_temperature(basics, grid, model, dust, gas, sources_in,fluxes)

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
    REAL(kind=r2)                                    :: QB
    
    INTEGER                                          :: i_cell
    INTEGER                                          :: i_dust
    INTEGER                                          :: i_tem
    INTEGER                                          :: i_a
    INTEGER                                          :: i_b
    INTEGER                                          :: i_c
    
    !--------------------------------------------------------------------------!

    IF (.not. basics%old_model) THEN
        ! First set minimal internal cell energy
        DO i_dust = 1, dust%n_dust
            grid%cell_energy(i_dust, :) = dust%QB(0, i_dust)*                  &
                                          basics%PIx4 * grid%cell_vol(:)
        END DO
        IF ( basics%do_MC_temperature) THEN
            IF (.not. SourcesInitialized(sources_in)) THEN
                print *, 'ERROR: There is no source of radiation defined'
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
                        moco(1) = (grid%co_mx_a(i_a) + grid%co_mx_a(i_a - 1)) / 2.0_r2
                        moco(2) = (grid%co_mx_b(i_b) + grid%co_mx_b(i_b - 1)) / 2.0_r2
                        moco(3) = (grid%co_mx_c(i_c) + grid%co_mx_c(i_c - 1)) / 2.0_r2
                        
                        caco = mo2ca(grid, moco)
                        DO i_dust = 1, dust%n_dust
                            grid%t_dust(i_cell, i_dust) = get_temperature(caco)
                            
                            ! calculate corresponding internal cell energy
                            i_tem = binary_search(grid%t_dust(i_cell, i_dust), &
                                                  dust%tem_tab(:))-1
                            i_tem = MIN(i_tem, basics%n_tem-1)
                            QB = ipol2(REAL(dust%tem_tab(i_tem), kind=r2),     &
                                       REAL(dust%tem_tab(i_tem + 1), kind=r2), &
                                       dust%QB(i_tem, i_dust),                 &
                                       dust%QB(i_tem + 1, i_dust),             &
                                       REAL(grid%t_dust(i_cell, i_dust), kind=r2) )
                                                    
                            grid%cell_energy(i_dust, i_cell) = QB *            &
                                            PI * 4.0_r2 * grid%cell_vol(i_cell)
                        END DO
                    END DO
                END DO
            END DO
        END IF
        ! the gas temperature equals the dust temperature (of dust species 1)
        ! tbd.: this should be generalized
        grid%t_gas(:) = grid%t_dust(:, 1)
    END IF

    ! CALCULATE PROBERTIES, WHICH NEED TEMPERATURE INFORMATIONS

    IF (basics%do_velo_ch_map) THEN
        !
        DO i_cell = 1,grid%n_cell

            ! set line width in each cell (in fact the inverse value)
            IF (grid%t_gas(i_cell) > 1.0e-2) THEN
                
                grid%cell_gauss_a(i_cell) = 1.0_r2 /                           &
                                (sqrt(2.0_r2 * con_k * grid%t_gas(i_cell) /    &
                                     (gas%mol_weight * 1.0e-3_r2 / con_Na) +   &
                                      gas%v_turb**2))
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
        ! reset flux map
        fluxes%continuum_map(:,:,:,:) = 0.0_r2

        ! photon transfer
        print *,'| | calculate temperature with Monte Carlo method'
        IF (basics%do_peel_off) PRINT *,"| | Peel-off technique enabled"
        CALL MC_photon_transfer(basics, grid,model,                            &
                                dust, sources_in, fluxes,                      &
                                deposit_energy=.True.,                         &
                                peel_off=basics%do_peel_off)
        ! prepare & save final results
        print *, "| | save contiunuum maps and sed's (bined mode)"
        CALL save_continuum_map(model, basics, dust, fluxes, 1, basics%do_peel_off)

        CALL temp_final(basics, grid, dust)

    END SUBROUTINE MC_temp
      
  SUBROUTINE temp_final(basics, grid, dust)

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
            hd2 = (grid%cell_energy(i_dust, nr_cell) )/                        &
                  (basics%PIx4 * grid%cell_vol(nr_cell))
         
            ! [2] find corresponding temperature from QB integral
            i_tem = MIN(binary_search(hd2, dust%QB(:,i_dust))- 1, basics%n_tem)
            IF (i_tem == 0) THEN
                grid%t_dust(nr_cell, i_dust) = dust%tem_tab(i_tem)
            ELSE
                grid%t_dust(nr_cell, i_dust) =                                 &
                                    ipol2( dust%QB(i_tem-1, i_dust),           &
                                    dust%QB(i_tem, i_dust),                    &
                                    REAL(dust%tem_tab(i_tem-1), kind=r2),      &
                                    REAL(dust%tem_tab(i_tem), kind=r2), hd2 )
            END IF
        end do
    end do

    END SUBROUTINE temp_final
END MODULE temp_mod
