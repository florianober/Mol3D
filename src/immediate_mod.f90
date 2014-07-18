! ---
! routines required for the concept of immediate reemission [by Bjorkman & Wood, 2001]
! ---
MODULE immediate_mod
    USE datatype
    USE var_globalnew
    
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
    PUBLIC :: immediate, immediate_temp, temp_final
contains

  ! ################################################################################################
  ! immediate reemission RT concept (routine from mc3d but modified)
  ! see Bjorkman & Wood 2001, ApJ
  ! ---
  SUBROUTINE immediate(basics,rand_nr,grid,dust,photon, i_dust_action)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!
        
    integer, intent(in)                              :: i_dust_action

    integer       :: i_tem_2
    integer       :: i_lam, i_tem
    real(kind=r2) :: t_dust_new, hd1, rndx
    real(kind=r2) :: i_tem_hd2
    !--------------------------------------------------------------------------!
    
    !1.+ 2. estimate new temperature -> t_dust(nr_cell,i_dust_action) & store old temp
    CALL immediate_temp( basics, grid,dust,photon, i_dust_action, t_dust_new, i_tem)
                         
    ! ----------------------------------------------------------------------------------------------
    ! 3. estimate difference in SED

    i_tem_hd2 = (t_dust_new - basics%t_dust_min) / basics%d_tem
    i_tem_2   = nint( i_tem_hd2 )
    

    ! chose a wavelength in the difference SED randomly == new wavelength
    CALL RAN2(rand_nr,rndx)

    i_lam = MIN(binary_search(rndx,dust%QdB_dT_l_cdf(:,i_tem_2,i_dust_action))+1, dust%n_lam)

    photon%nr_lam         = i_lam
    photon%current_albedo = dust%albedo(:,i_lam)
  end subroutine immediate

  
  
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
    print *, "    final temperature calculation ..."
    hd2 = 0.0_r2
    do i_dust=1, dust%n_dust
        do nr_cell=1, grid%n_cell
            hd2    =  (grid%cell_energy(i_dust,nr_cell) )
         
            ! [2] find corresponding temperature from QB integral
            i_tem = MIN(binary_search(hd2,dust%QB(:,i_dust))-1,basics%n_tem)

            grid%t_dust(nr_cell,i_dust) = &
                ( ( (hd2 - dust%QB(i_tem-1,i_dust)) / (dust%QB(i_tem,i_dust) - dust%QB(i_tem-1,i_dust)) ) &
                + real(i_tem-1,kind=r2) ) &
                * basics%d_tem
        end do
    end do


  end subroutine temp_final
  

  ! ################################################################################################
  ! estimation of the  dust temperature distribution
  ! ---
  subroutine immediate_temp( basics,grid,dust,photon, i_dust, temp_new, i_tem)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!
        
    integer, intent(in)                             :: i_dust

    integer,INTENT(OUT)                             :: i_tem
    real(kind=r2)                                   :: hd1
    real(kind=r2),INTENT(OUT)                       :: temp_new
    !--------------------------------------------------------------------------!
    ! ---
    !TbD, check the if case .... sounds a bit arbitrary (from mc3d .....)
    if ( grid%Nv(photon%nr_cell,i_dust) > 1.0e-16_r2) then
        ! absorb photon
            hd1    = (grid%cell_energy_sum(i_dust,photon%nr_cell,1)/ &
                     (basics%PIx4 * grid%cell_vol(photon%nr_cell)))

            i_tem = MIN(binary_search(hd1,dust%QB(:,i_dust))-1,basics%n_tem)

            temp_new = &
                ( ( (hd1 - dust%QB(i_tem-1,i_dust)) / (dust%QB(i_tem,i_dust) - dust%QB(i_tem-1,i_dust)) ) &
                + real(i_tem-1,kind=r2) ) &
                * basics%d_tem
    ELSE
        temp_new = basics%t_dust_min
        i_tem    = 1
    END IF
  end subroutine immediate_temp
  
end module immediate_mod

