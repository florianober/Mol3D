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
    PUBLIC :: immediate, immediate_temp, temp_final2!,  temp_final1
contains

  ! ################################################################################################
  ! immediate reemission RT concept (routine from mc3d but modified)
  ! see Bjorkman & Wood 2001, ApJ
  ! ---
  SUBROUTINE immediate(basics,rand_nr,grid,dust,photon, i_dust_action,t_dust, i_star_abs,dt_dust )
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                    :: photon
    !--------------------------------------------------------------------------!
        
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: t_dust
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)      :: dt_dust
    REAL(kind=r2),DIMENSION(1:dust%n_dust, 1:dust%n_lam, 1:grid%n_cell),INTENT(INOUT)    :: i_star_abs
    !--------------------------------------------------------------------------!
    integer, intent(in)                             :: i_dust_action

    integer       :: nr_lam_old, i_tem_1, i_tem_2
    integer       :: i_lam
    real(kind=r2) :: t_dust_old,t_dust_new, hd1, rndx
    real(kind=r2) :: i_tem_hd1, i_tem_hd2
    REAL(KIND=r2),DIMENSION(1:dust%n_lam)  :: pt_1, pt_2
    !--------------------------------------------------------------------------!
    
    !1.+ 2. estimate new temperature -> t_dust(nr_cell,i_dust_action) & store old temp
    CALL immediate_temp( basics, grid,dust,photon, i_dust_action, t_dust_old, t_dust_new,&
                         t_dust, i_star_abs)
    ! 3. estimate difference in SED
    i_tem_hd1 = (t_dust_old - basics%t_dust_min) / basics%d_tem
    !print *,t_dust(photon%nr_cell,1)
    dt_dust(photon%nr_cell,1) = abs(t_dust_old - t_dust_new) 
    i_tem_1   = nint( i_tem_hd1 ) 
    
    nr_lam_old  = photon%nr_lam

    ! interpolate: tabulated planck function @ t_dust_old
!~     pt_1 = ipol2( &
!~                 real(      i_tem_1,kind=r2),     real(i_tem_1+1,kind=r2), &
!~             dust%planck_tab(i_tem_1,:), dust%planck_tab(i_tem_1+1,:), &
!~             i_tem_hd1 )
        
       ! interpolate: tabulated planck function @ t_dust(i_dust_action,nr_cell)
    i_tem_hd2 = (t_dust_new - basics%t_dust_min) / basics%d_tem
    i_tem_2   = nint( i_tem_hd2 )
!~     if (i_tem_2 > basics%n_tem) then
!~         print *,"Failure : in subroutine immediate(): i_tem > n_tem"
!~         print *,"Solution: Chose more photons or fewer grid cells"
!~         stop
!~     else
!~         pt_2   = ipol2( &
!~                     real(      i_tem_2,kind=r2),     real(i_tem_2+1,kind=r2), &
!~                     dust%planck_tab(i_tem_2,:), dust%planck_tab(i_tem_2+1,:), &
!~                     i_tem_hd2 )
!~        end if

!~     photon%diff_planck(:) = dust%C_abs(i_dust_action,:) * ( pt_2 - pt_1 )
!~     print *, i_tem_2, t_dust_new
!~     photon%diff_planck(:) = dust%QdB_dT_l(i_dust_action,i_tem_2, :)!*(t_dust_new - t_dust_old)

    ! ----------------------------------------------------------------------------------------------
    ! difference > 0 ?
!~     if ( sum(abs(photon%diff_planck(:))) /= 0.0_r2 ) then
    
       ! normalize the difference-SED
!~        photon%diff_planck(:) = photon%diff_planck(:) / integ1(dust%lam(:), &
!~                                  photon%diff_planck(:), 1, dust%n_lam)

       ! chose a wavelength in the difference SED randomly == new wavelength
        CALL RAN2(rand_nr,rndx)
        i_lam = 0
!~         print *,'before seach'
!~         i_lam = MIN(binary_search(rndx,ipol2(real(i_tem_2,kind=r2),real(i_tem_2+1,kind=r2), &
!~                                    dust%QdB_dT_l_norm(i_dust_action,i_tem_2,:),         &
!~                                    dust%QdB_dT_l_norm(i_dust_action,i_tem_2+1,:),i_tem_hd2))+1, &
!~                                    dust%n_lam)
        i_lam = MIN(binary_search(rndx,dust%QdB_dT_l_norm(i_dust_action,i_tem_2,:))+1, dust%n_lam)
!~         if (i_lam > dust%n_lam) THEN
!~             print *,rndx,dust%QdB_dT_l_norm(i_dust_action,i_tem_2,i_lam),i_lam
!~         end if
!~         print *,'after seach'
!~         print *,i_lam, dust%n_lam
!~     else
!~         i_lam = dust%n_lam-1
!~        
!~     endif
    photon%nr_lam         = i_lam
!~     print *, i_lam
!~     photon%energy         = photon%energy * dust%d_lam(nr_lam_old)/dust%d_lam(photon%nr_lam)
!~     photon%energy         = photon%energy
    photon%current_albedo = dust%albedo(:,i_lam)
  end subroutine immediate


  ! ################################################################################################
  ! - estimation of the final dust temperature distribution
  !   solution 1: via immediate_temp
  ! ---
!~   subroutine temp_final1()
!~     use var_global
!~ 
!~     implicit none
!~     integer :: i_dust
!~     ! ---
!~     print *, "    final temperature calculation [solution 1] ..."
!~     n_wrong_temp = 0
!~     do i_dust=1, n_dust
!~        do nr_cell=1, n_cell
!~           call immediate_temp( i_dust )
!~           if (t_dust(nr_cell,i_dust) > t_dust_max) then
!~              n_wrong_temp = n_wrong_temp +1
!~           end if
!~        end do
!~     end do
!~   end subroutine temp_final1
  

!~     FUNCTION temp_solver(basics, model, grid, dust,i_cell,N) RESULT (temp)
!~         ! calculate the temperatur for a given cell number and number of interactions (N)
!~         IMPLICIT NONE
!~         !--------------------------------------------------------------------------!
!~         TYPE(Basic_TYP),INTENT(IN)                       :: basics
!~         TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
!~         TYPE(Dust_TYP),INTENT(IN)                        :: dust
!~         TYPE(Model_TYP),INTENT(IN)                       :: model
!~ 
!~         !--------------------------------------------------------------------------!
!~         integer                                          :: i_dust, i_tem, i_cell
!~         real(kind=r2)                                    :: hd2
!~         real(kind=r2)                                    :: temp
!~         !--------------------------------------------------------------------------!
!~         
!~         temp = 0.0_r2
!~         
!~         hd1(:) = dust%C_abs(i_dust,:)  *  model%ref_unit*grid%grd_d_l(nr_cell,:)
!~         hd2    = (1.0_r2 / (basics%PIx4 * grid%cell_vol(nr_cell))) * integ1( dust%lam(:), hd1(:), 1, dust%n_lam)
!~              
!~         ! [2] find corresponding temperature from QB integral
!~         ! tbd: apply faster search algorithm
!~         ! tbd: interpolate: to get smooth temperature distribution instead of "steps"
!~         i_tem = 0
!~         do 
!~            if ( dust%QB(i_dust,i_tem) > hd2 ) then
!~               exit
!~           else
!~             i_tem = i_tem + 1
!~             if (i_tem == basics%n_tem) then
!~                exit
!~             end if
!~             cycle
!~          end if
!~         end do
!~               grid%t_dust(nr_cell,i_dust) =  dust%tem_tab(i_tem)
!~         
!~ 
!~     END FUNCTION temp_solver
  
  
  SUBROUTINE temp_final2(basics, model, grid, dust)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(Model_TYP),INTENT(IN)                       :: model

    !--------------------------------------------------------------------------!
    integer                                          :: i_dust, i_tem, nr_cell
    real(kind=r2)                                    :: hd2
    !--------------------------------------------------------------------------!
    
    ! ---
    print *, "    final temperature calculation [solution 2] ..."
!~     allocate( hd1(1:dust%n_lam) )
    hd2 = 0.0_r2
    do i_dust=1, dust%n_dust
        do nr_cell=1, grid%n_cell
            ! [1] absorption
            ! tbd: if Nv=0 => t=0
            hd2    =  (grid%cell_energy_sum(i_dust,nr_cell) )&
                    / (basics%PIx4 * grid%cell_vol(nr_cell))
         
            ! [2] find corresponding temperature from QB integral
            i_tem = MIN(binary_search(hd2,dust%QB(i_dust,:))-1,basics%n_tem)

            grid%t_dust(nr_cell,i_dust) = &
                ( ( (hd2 - dust%QB(i_dust,i_tem-1)) / (dust%QB(i_dust,i_tem) - dust%QB(i_dust,i_tem-1)) ) &
                + real(i_tem-1,kind=r2) ) &
                * basics%d_tem
        end do
    end do


  end subroutine temp_final2
  

  ! ################################################################################################
  ! estimation of the  dust temperature distribution
  ! ---
  subroutine immediate_temp( basics,grid,dust,photon, i_dust_action, temp_old, temp_new,& 
                             t_dust, i_star_abs)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                    :: photon
    !--------------------------------------------------------------------------!
        
    REAL(kind=r1),DIMENSION(0:grid%n_cell,1:dust%n_dust),INTENT(INOUT)                  :: t_dust
    REAL(kind=r2),DIMENSION(1:dust%n_dust, 1:dust%n_lam, 1:grid%n_cell),INTENT(INOUT)    :: i_star_abs
    !--------------------------------------------------------------------------!

    integer, intent(in)                             :: i_dust_action

    integer                                          :: i_tem
    real(kind=r2)                                    :: hd1
    real(kind=r2),INTENT(OUT)                       :: temp_new
    real(kind=r2),INTENT(OUT)                       :: temp_old
    !--------------------------------------------------------------------------!
    ! ---
    !print *, 'here1', grid%Nv_r(photon%nr_cell,i_dust_action)
    !TbD, check the if case .... sounds arbitrary
    if ( grid%Nv(photon%nr_cell,i_dust_action) > 1.0e-16_r2) then
        ! absorb photon
            hd1    = (grid%cell_energy_sum(i_dust_action,photon%nr_cell)/ &
                     (basics%PIx4 * grid%cell_vol(photon%nr_cell)))

            i_tem = MIN(binary_search(hd1,dust%QB(i_dust_action,:))-1,basics%n_tem)

            temp_new = &
                ( ( (hd1 - dust%QB(i_dust_action,i_tem-1)) / (dust%QB(i_dust_action,i_tem) - dust%QB(i_dust_action,i_tem-1)) ) &
                + real(i_tem-1,kind=r2) ) &
                * basics%d_tem
            ! in some cases the temperature may be > maximum temperature
            ! mark these cells
            ! tbd: solve this bug
            if ( temp_new > basics%t_dust_max) then
                temp_new = basics%t_dust_min
            end if
        
            temp_old = t_dust(photon%nr_cell,i_dust_action)
            t_dust(photon%nr_cell,i_dust_action)      = temp_new
    ELSE
        temp_old = basics%t_dust_min
        temp_new = basics%t_dust_min
    END IF
  end subroutine immediate_temp
  
end module immediate_mod

