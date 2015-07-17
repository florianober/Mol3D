! ---
! routines required for the concept of immediate temperature correction
!         [by Bjorkman & Wood, 2001]
! ---
MODULE immediate_mod
    USE datatype
    USE var_global

    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE photon_type
    USE basic_type
    USE fluxes_type

    USE math_mod

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: immediate

CONTAINS

    ! ##########################################################################
    ! immediate reemission RT concept
    ! see Bjorkman & Wood 2001, ApJ
    ! ---
        SUBROUTINE immediate(basics, rand_nr, grid, dust,photon, i_dust_action)
      
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!

        INTEGER, INTENT(IN)                              :: i_dust_action

        INTEGER                                          :: i_lam, i_tem
        REAL(kind=r2)                                    :: rndx
        !----------------------------------------------------------------------!
        ! 1. find energy in cell
        !
        CALL immediate_temp( basics, grid, dust, photon%nr_cell,               &
                             i_dust_action, i_tem)

        ! 2. chose a wavelength in the difference
        !    SED randomly == new wavelength, the photon energy is not
        !    changed
        CALL GetNewRandomNumber(rand_nr,rndx)
        i_lam = MIN(binary_search(rndx,                                    &
                    dust%QdB_dT_l_cdf(:,i_tem,i_dust_action))+1,           &
                    dust%n_lam)
        photon%nr_lam       = i_lam
    END SUBROUTINE immediate


    ! ##########################################################################
    ! estimation of the  dust temperature in the cell
    ! ---
    SUBROUTINE immediate_temp( basics,grid,dust, i_cell, i_dust, i_tem)
  
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        
        !----------------------------------------------------------------------!
            
        integer, intent(in)                              :: i_dust
        integer, intent(in)                              :: i_cell

        integer,INTENT(OUT)                              :: i_tem
        real(kind=r2)                                    :: hd1
        !----------------------------------------------------------------------!
        ! ---
        ! 
        hd1    = (grid%cell_energy_sum(i_dust, i_cell,1)/                      &
                 (basics%PIx4 * grid%cell_vol(i_cell)))

        i_tem = MIN(binary_search(hd1, dust%QB(:,i_dust))-1, basics%n_tem)
        IF ( (hd1 - dust%QB(i_tem-1,i_dust))/                                  &
             (dust%QB(i_tem,i_dust) - dust%QB(i_tem-1,i_dust))                 &
                     .lt. 0.5_r2 ) i_tem = i_tem-1

    END SUBROUTINE immediate_temp

END MODULE immediate_mod
