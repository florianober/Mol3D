! ---
! photon transfer through model space
! (optimized for specific geometry)
! ---
module transfer_mod

    USE datatype
    USE var_globalnew
    
    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE simu_type
    
    USE math_mod
    USE tools_mod
    USE grd_mod

    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC  ::  path, &
                path_skip,next_pos_const
    !--------------------------------------------------------------------------!
contains

    ! ################################################################################################
    ! photon transfer through model space
    ! here: constant density in each ESC is assumed
    ! ---
    subroutine next_pos_const(model, rand_nr, grid, dust, simu_var,grd_d_l)

        IMPLICIT NONE
        !--------------------------------------------------------------------------!    
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Model_TYP),INTENT(IN)                        :: model
        TYPE(Dust_TYP),INTENT(IN)                         :: dust
        
        TYPE(Simu_TYP),INTENT(INOUT)                      :: simu_var
        !--------------------------------------------------------------------------!
        
        REAL(kind=r2),DIMENSION(:,:),INTENT(INOUT)      :: grd_d_l
        !--------------------------------------------------------------------------
            
        logical                                            :: kill_photon
        !integer                                            :: nr_cell_new
        real(kind=r2)                                      :: tau_end, d_l, d_tau, rndx
 
        !real(kind=r2), dimension(1:3)                      :: pos_xyz_new
        
        !--------------------------------------------------------------------------!
        ! ---
        ! 1. determine optical depth (distance) to next point of interaction
        CALL RAN2(rand_nr,rndx)
        tau_end = -log(1.0_r2 - rndx)
        !print *, 'tau', tau_end,'rndx', rndx
        ! 2. go to next point of interaction
!~         print *,'now'
        DO
            ! 2.1. inside the dust sublimation radius: skip
            IF (simu_var%nr_cell==0) THEN
!~                 print *, 'here?'
                CALL path_skip( grid, simu_var,rand_nr)
!~                 print *, 'here? huch'
                simu_var%pos_xyz = simu_var%pos_xyz_new
                simu_var%nr_cell = simu_var%nr_cell_new
                IF (check_inside(simu_var%pos_xyz_new(:),grid, model) ) THEN
                    kill_photon = .True.
                    exit
                END IF

        
            END IF
            ! 2.2. determine: - geometric path length in current cell (d_l [ref_unit]);
            !                 - new entry point in neighbouring cell
            !                 - new cell number
            CALL path( grid, simu_var%pos_xyz, simu_var%pos_xyz_new, simu_var%nr_cell, &
                        simu_var%nr_cell_new, d_l, kill_photon, simu_var%dir_xyz)
            !print *, 'pos',simu_var%pos_xyz
            IF (kill_photon) THEN
                ! stop transfer of this photon package
                d_l = 0.0_r2
                ! set condition for leaving the loop
                d_tau = tau_end * 1.00001_r2

                ! move photon package outside model space ! maybe we should generalize this,
                ! but however it should work this way
                simu_var%pos_xyz_new(1)   = model%r_ou * 1.00001_r2
                simu_var%pos_xyz_new(2:3) = 0.0_r2
              
            ELSE
                ! 2.3. determine corresponding optical path length        (d_tau)
                d_tau = (d_l * model%ref_unit) * sum(dust%C_ext(:,simu_var%nr_lam) *&
                         grid%grd_dust_density(simu_var%nr_cell,:))
                !print *, 'd_tau', sum(dust%C_ext(:,simu_var%nr_lam)),simu_var%nr_lam

            END IF

            ! 2.4. check status
            IF (d_tau > tau_end) THEN
                !print *,'higher',d_tau, simu_var%nr_cell
                ! optical depth and thus final cell reached, 
                ! but exact position of new point of interaction has yet to be determined
                ! ---
                ! a) adjust geometrical fractional path
                
                d_l = d_l * tau_end / d_tau
              
                ! b) consider only fraction of path in current cell (i.e., nr_cell_new = nr_cell)
                simu_var%pos_xyz_new(:) = simu_var%pos_xyz(:) + simu_var%dir_xyz(:) * d_l

                ! c) store information about (fractional) path through last cell
                grd_d_l(simu_var%nr_cell,simu_var%nr_lam) = grd_d_l(simu_var%nr_cell,simu_var%nr_lam) + &
                                                                              d_l * simu_var%c_in_akt
                ! d) new point of interaction
                simu_var%pos_xyz(:) = simu_var%pos_xyz_new(:)

                ! e) nr_cell: remains unchanged (photon does not leave current cell)
                ! simu_var%nr_cell = simu_var%nr_cell_new

                exit
      
            ELSE
                !print *, 'lower',d_tau, simu_var%nr_cell
                ! optical depth has not yet been reached
                IF ( check_inside(simu_var%pos_xyz_new(:),grid, model) ) THEN
                    ! photon still inside model space
                    ! ---
                    ! a) store information about (fractional) path through last cell
                    
                    !#######################################################
                    ! be careful here
                    grd_d_l(simu_var%nr_cell,simu_var%nr_lam) = grd_d_l(simu_var%nr_cell,simu_var%nr_lam) +& 
                                                                                  d_l * simu_var%c_in_akt
                    !#######################################################
                    ! b) set new starting point; adjust optical depth
                    simu_var%pos_xyz(:) = simu_var%pos_xyz_new(:)             
                    simu_var%nr_cell    = simu_var%nr_cell_new
                    tau_end    = tau_end - d_tau
                    cycle

                ELSE
                    ! photon already outside model space
                    ! ---
                    ! a) store information about (fractional) path through last cell
                    grd_d_l(simu_var%nr_cell,simu_var%nr_lam) = grd_d_l(simu_var%nr_cell,simu_var%nr_lam) + & 
                                                                                 d_l * simu_var%c_in_akt
                 
                    ! b) set flag
                    simu_var%inside = .false.
                    exit
                END IF
            END IF
        END DO
!~         print *,'now_let_go'
    end subroutine next_pos_const
  
  
  ! ################################################################################################
  ! determine geometrical distance between given point and next cell boundary
  ! ---
  SUBROUTINE path_skip( grid,simu_var,rand_nr)
    ! steering routine
    IMPLICIT NONE
    !--------------------------------------------------------------------------!

    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Simu_TYP),INTENT(INOUT)                     :: simu_var
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    !--------------------------------------------------------------------------!
    
    SELECT CASE(GetGridName(grid))
        
    CASE('spherical')
        CALL path_skip_sp( grid,simu_var)
    CASE('cylindrical')
        CALL path_skip_cy( grid,simu_var,rand_nr)

    CASE('cartesian')
        print *, 'TbD, not finished yet, path'
        stop
    CASE DEFAULT
        print *, 'selected coordinate system not found, path'
        stop
    END SELECT

  END SUBROUTINE path_skip
  
  SUBROUTINE path_skip_cy( grid, simu_var,rand_nr)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!

    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Simu_TYP),INTENT(INOUT)                     :: simu_var
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    !--------------------------------------------------------------------------!
!~     integer :: i_r, i_th, i_ph
    
    real(kind=r2) :: d_l1,d_l2,d_l
    real(kind=r2) :: a,b,c, rndx
!~     real(kind=r2) :: hd_r1, hd1, hd2, hd3, hd_th, hd_ph
    real(kind=r2), dimension(1:3) :: pos_xyz,pos_xyz_new, dir_xyz
    !--------------------------------------------------------------------------!
    ! ---
    ! starting point: pos_xyz
    ! direction     : dir_xyz 
    ! cell number   : nr_cell
    ! ---
    dir_xyz = simu_var%dir_xyz
    pos_xyz = simu_var%pos_xyz
    
    a = sum(dir_xyz(1:2)**2)
    
    IF (a .lt. 1.0e-9) THEN
        d_l1 = grid%co_mx_a(grid%n(1))+1.1_r2
    ELSE
        b = sum(dir_xyz(1:2)*pos_xyz(1:2))
        c = sum(pos_xyz(1:2)**2)-grid%co_mx_a(0)**2
        
        d_l1 =  -b/a + ((b/a)**2-c/a)**0.5
    END IF

    d_l = d_l1
    d_l = d_l +epsilon(1.0_r2)*1.0e6_r2*d_l
    
    pos_xyz_new    = pos_xyz + d_l * dir_xyz
    simu_var%nr_cell_new = get_cell_nr( grid, pos_xyz_new)
!~     IF ( ((pos_xyz_new(1)**2+pos_xyz_new(2)**2)**0.5 .lt. 2.0_r2)) THEN
!~         print *, 'i_cell is still zero in skip path cy'
!~         print *, ((pos_xyz_new(1)**2+pos_xyz_new(2)**2)**0.5 .lt. 2.0_r2)
!~         print *, (pos_xyz_new(1)**2+pos_xyz_new(2)**2)**0.5
!~         print *, d_l
!~         print *, pos_xyz
!~         print *, pos_xyz_new
!~     END IF
    
    simu_var%pos_xyz_new = pos_xyz_new
  END SUBROUTINE path_skip_cy
  
  subroutine path_skip_sp( grid, simu_var )
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!

    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Simu_TYP),INTENT(INOUT)                     :: simu_var
    !--------------------------------------------------------------------------!
    integer :: i_r, i_th, i_ph
    
    real(kind=r2) :: d_l
    real(kind=r2) :: hd_r1, hd1, hd2, hd3, hd_th, hd_ph
    real(kind=r2), dimension(1:3) :: p0_vec, d_vec
    !--------------------------------------------------------------------------!
    
    ! ---
    ! starting point: pos_xyz
    ! direction     : dir_xyz 
    ! cell number   : nr_cell
    ! ---
    p0_vec(:) = simu_var%pos_xyz(:)
    d_vec(:)  = simu_var%dir_xyz(:)

    ! ---
    ! 1. determine distance to inner shell radius (r_in)
    hd_r1 = grid%co_mx_a(0)
    
    hd1 = dot_product(p0_vec, d_vec)
    hd2 = dot_product( d_vec, d_vec)
    hd3 = dot_product(p0_vec,p0_vec)
    d_l = -(hd1/hd2) + sqrt((hd1/hd2)**2 - (hd3-(hd_r1**2))/hd2)

    ! -
    ! 2. new coordinate
    !    rem.: entering the inner cell by 1e-4
    d_l = d_l  +  ((grid%co_mx_a(1) - grid%co_mx_a(0)) / 1.0e+4_r2)

    simu_var%pos_xyz_new(:) = simu_var%pos_xyz(:) +   d_l * simu_var%dir_xyz(:)
    
    ! -
    ! 3. determine new cell number
    i_r = 1
    
    hd_th = atan2( simu_var%pos_xyz_new(3), sqrt(simu_var%pos_xyz_new(1)**2 + simu_var%pos_xyz_new(2)**2))
    i_th = floor( real(grid%n(2),kind=r2) * (hd_th + PI*0.5_r2)/PI ) +1
    
    hd_ph = atan3( simu_var%pos_xyz_new(2), simu_var%pos_xyz_new(1) )
    i_ph = floor( real(grid%n(3),kind=r2) * hd_ph/PI*2.0_r2       ) +1

    ! check if number of <th> cell is in range [1,n_th]
    if      (i_th < 1   ) then
       i_th = 1
       i_ph = i_ph + grid%n(3)/2
    else if (i_th > grid%n(2)) then
       i_th = 1
       i_ph = i_ph + grid%n(3)/2             
    end if
    
    ! check if number of <ph> cell is in range [1,n_ph]
    if      (i_ph < 1   ) then
       i_ph = grid%n(3)
    else if (i_ph > grid%n(3)) then
       i_ph = 1             
    end if

    simu_var%nr_cell_new = grid%cell_idx2nr( i_r, i_th, i_ph )
    !print *, simu_var%pos_xyz
  end subroutine path_skip_sp


  ! ################################################################################################
  ! determine geometrical distance between given point and next cell boundary
  ! ---
  SUBROUTINE path( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, dir_xyz)

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(IN)                         :: grid
    !--------------------------------------------------------------------------!
    
    real(kind=r2),               intent(out)          :: d_l
    real(kind=r2), dimension(1:3), intent(in)        :: pos_xyz
    real(kind=r2), dimension(1:3), intent(out)       :: pos_xyz_new
    real(kind=r2), dimension(1:3), intent(in)        :: dir_xyz
    integer,                     intent(in)           :: nr_cell
    integer,                     intent(out)          :: nr_cell_new

    logical,                     intent(out)          :: kill_photon
    !logical                                           :: show_error
    !--------------------------------------------------------------------------!
    
    
!~     grid%counter = grid%counter +1
    SELECT CASE(GetGridName(grid))
        
    CASE('spherical')
        CALL path_sp( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, dir_xyz)
    CASE('cylindrical')

        CALL path_cy( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, dir_xyz)

    CASE('cartesian')
        print *, 'TbD, not finished yet, path'
        stop
    CASE DEFAULT
        print *, 'selected coordinate system not found, path'
        stop
    END SELECT

  END SUBROUTINE path
  
  SUBROUTINE path_cy( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, dir_xyz)
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(IN)                         :: grid
    !--------------------------------------------------------------------------!
    
    real(kind=r2),               intent(out)         :: d_l
    real(kind=r2), dimension(1:3), intent(in)        :: pos_xyz
    real(kind=r2), dimension(1:3), intent(out)       :: pos_xyz_new
    real(kind=r2), dimension(1:3), intent(in)        :: dir_xyz
    
    integer,                     intent(in)           :: nr_cell
    integer,                     intent(out)          :: nr_cell_new

    logical,                     intent(out)          :: kill_photon
    !logical                                           :: show_error
    !--------------------------------------------------------------------------!
    real(kind=r2),dimension(6) :: d_l_selc
    real(kind=r2) :: d_r, p_r, d_ph, p_ph, a, p, q, sq
    integer       :: i_r, i_ph, i_z, i, loca
    
    !--------------------------------------------------------------------------!
    ! ---
    ! default:
    kill_photon = .false.
    
    i_r  = grid%cell_nr2idx(1,nr_cell)
    i_ph = grid%cell_nr2idx(2,nr_cell)
    i_z  = grid%cell_nr2idx(3,nr_cell)
    
!~     print *, i_r, i_ph, i_z
    !pos_xyz(1) = -pos_xyz(1)
    
    d_r     = sqrt(dir_xyz(1)**2+dir_xyz(2)**2)
    p_r     = sqrt(pos_xyz(1)**2+pos_xyz(2)**2)
    
    d_ph    = atan3(dir_xyz(2),dir_xyz(1))
    p_ph    = atan3(pos_xyz(2),pos_xyz(1))
    
    d_l_selc(:) = 0.0_r2
    
    ! find z component
    
    IF ( abs(dir_xyz(3)) .gt. epsilon(d_r) ) THEN
    
        d_l_selc(1) = (grid%co_mx_c(i_z)-pos_xyz(3))/dir_xyz(3)
        d_l_selc(2) = (grid%co_mx_c(i_z-1)-pos_xyz(3))/dir_xyz(3)
    ELSE 
        d_l_selc(1) = -1.0_r2
        d_l_selc(2) = -1.0_r2
    END IF
    
    ! find rho component

    p = 2.0_r2*(dir_xyz(1)*pos_xyz(1)+dir_xyz(2)*pos_xyz(2)) / &
            (dir_xyz(1)**2+dir_xyz(2)**2)
    q = (pos_xyz(1)**2+pos_xyz(2)**2-grid%co_mx_a(i_r)**2)/ &
        (dir_xyz(1)**2+dir_xyz(2)**2)
    sq = p**2/4.0-q
    !print *,'1',sq
    
    !print *, -p/2.0_r2 + sqrt(p**2/4.0_r2-q)
    
    IF ( sq .gt. epsilon(sq) ) THEN
    
        d_l_selc(3) = -p/2.0_r2 + sqrt(p**2/4.0_r2-q)
    
    ELSE
        d_l_selc(3) = -1.0_r2
    END IF
    
    q = (pos_xyz(1)**2+pos_xyz(2)**2-grid%co_mx_a(i_r-1)**2)/ &
        (dir_xyz(1)**2+dir_xyz(2)**2)  
        
    sq = p**2/4.0-q       
    
    IF ( sq .gt. epsilon(sq) ) THEN
    
        d_l_selc(4) = -p/2.0_r2 - sqrt(p**2/4.0_r2-q)
    
    ELSE
        d_l_selc(4) = -1.0_r2
    END IF

    
    ! find phi component   ! test this for more than 1 grid cell in phi direction!
!~     IF (grid%n(2) .gt. 1 ) THEN
    IF (d_ph .gt. epsilon(d_ph) .and. grid%n(2) .gt. 1  ) THEN
        a = 1.0_r2/TAN(grid%co_mx_b(i_ph))
        d_l_selc(5) = (pos_xyz(1)-pos_xyz(2)*a)/(dir_xyz(2)*a-dir_xyz(1))
        
        a = 1.0_r2/TAN(grid%co_mx_b(i_ph-1))
        d_l_selc(6) = (pos_xyz(1)-pos_xyz(2)*a)/(dir_xyz(2)*a-dir_xyz(1))
    ELSE
        d_l_selc(5) = -1.0_r2
        d_l_selc(6) = -1.0_r2
    END IF
    !print *, d_r
    !print *, d_l_selc(:)
    
    ! find shortest path and leaving boundary
    
    d_l = 1.0e15_r2
    DO i=1,6
!~         print *, d_l_selc(i)
        IF ( d_l_selc(i) .gt. epsilon(d_l) .and. d_l_selc(i) .lt. d_l ) THEN
            d_l = d_l_selc(i)
            loca = i 
        END IF
    END DO
    d_l = d_l + 1.0e6_r2*epsilon(d_l)
    nr_cell_new = get_cell_nr(grid,d_l*dir_xyz + pos_xyz)
    pos_xyz_new = d_l*dir_xyz + pos_xyz
    
    
    IF (d_l < grid%d_l_min) then
       ! step width too small
        kill_photon       = .true.
        kill_photon_count = kill_photon_count +1
        d_l = d_l+1.0e3_r2*epsilon(d_l) 
    end if

    !print *, pos_xyz
    !print *, pos_xyz_new
    !print *, nr_cell
    !print *, nr_cell_new
  
  END SUBROUTINE path_cy
  
  
  SUBROUTINE path_sp( grid, pos_xyz, pos_xyz_new, nr_cell, nr_cell_new, d_l, kill_photon, dir_xyz)
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Grid_TYP),INTENT(IN)                         :: grid
    !--------------------------------------------------------------------------!
    
    real(kind=r2),               intent(out)         :: d_l
    real(kind=r2), dimension(1:3), intent(in)        :: pos_xyz
    real(kind=r2), dimension(1:3), intent(out)       :: pos_xyz_new
    real(kind=r2), dimension(1:3), intent(in)        :: dir_xyz
    
    
    integer,                     intent(in)           :: nr_cell
    integer,                     intent(out)          :: nr_cell_new

    logical,                     intent(out)          :: kill_photon
    !logical                                           :: show_error
    !--------------------------------------------------------------------------!
    integer                                           :: hi1, i1, p_i1
    integer, dimension(1)                             :: hi_arr1
    real(kind=r2) :: hd_r1, hd_ph1, hd_th1, hd1, hd2, hd3, hd4, hd5, hd6, hd7
    real(kind=r2), dimension(1:6)                     :: d_lx, d_lx_sel
    real(kind=r2), dimension(1:3)                     :: p0_vec, d_vec
    !real(kind=r2), dimension(1:3)                     :: spco
    real(kind=r2), dimension(1:2)                     :: hd_arr1
    
    !--------------------------------------------------------------------------!
    ! ---
    ! default:
    kill_photon = .false.

    ! starting point: pos_xyz (grid type variable)
    ! direction     : dir_xyz (grid type variable)
    ! cell number   : nr_cell (grid type variable)
    p0_vec(:) = pos_xyz(:)
    d_vec(:)  = dir_xyz(:)

    
    ! ---
    ! 1. determine distance to the boundaries of the cell
    ! 1.1. r_min, r_max
    ! inner boundary:
    hd_r1 = grid%co_mx_a( grid%cell_nr2idx(1,nr_cell) -1 )

    hd1 = dot_product(p0_vec, d_vec)
    hd2 = dot_product(p0_vec,p0_vec)

    hd3 = hd1**2 - (hd2-(hd_r1**2))
    if (hd3 < 0.0_r2) then
       d_lx(1) = -1.0_r2
    else
       hd_arr1(1) = -hd1 + sqrt(hd3)
       hd_arr1(2) = -hd1 - sqrt(hd3)
       if (hd_arr1(1)<0.0_r2 .or. hd_arr1(2)<0.0_r2) then
          d_lx(1) = maxval(hd_arr1(:))
       else
          d_lx(1) = minval(hd_arr1(:))
       end if
    end if

    ! outer boundary:
    hd_r1 = grid%co_mx_a( grid%cell_nr2idx(1,nr_cell)    )
    hd3 = hd1**2 - (hd2-(hd_r1**2))
    if (hd3 < 0.0_r2) then
       d_lx(2) = -1.0_r2
    else
       hd_arr1(1) = -hd1 + sqrt(hd3)
       hd_arr1(2) = -hd1 - sqrt(hd3)
       if (hd_arr1(1)<0.0_r2 .or. hd_arr1(2)<0.0_r2) then
          d_lx(2) = maxval(hd_arr1(:))
       else
          d_lx(2) = minval(hd_arr1(:))
       end if
    end if
    
    IF  (d_lx(2) .eq. 0.0_r2 ) THEN
        d_lx(2) = grid%d_l_min *1.1_r2
    ELSE IF (d_lx(1) .eq. 0.0_r2 ) THEN
        d_lx(1) = grid%d_l_min *1.1_r2
    END IF
    !print *,d_lx(1), d_lx(2)
    ! -
    ! 1.2. theta_min, theta_max
    if (grid%n(2) > 1) then
       ! inner
       hd_th1 = grid%co_mx_b( grid%cell_nr2idx(2,nr_cell) -1 )
       hd1 = tan(hd_th1)**2
       hd2 = hd1*( d_vec(1)**2        +  d_vec(2)**2       ) -  d_vec(3)**2
       hd3 = hd1*( d_vec(1)*p0_vec(1) +  d_vec(2)*p0_vec(2)) -  d_vec(3)*p0_vec(3)
       hd4 = hd1*(p0_vec(1)**2        + p0_vec(2)**2       ) - p0_vec(3)**2
       
       hd5 = hd3/hd2
       hd6 = hd4/hd2
       
       hd7 = hd5**2 - hd6
       if (hd7 < 0.0_r2) then
          d_lx(3) = -1.0_r2
       else
          hd_arr1(1) = -hd5 + sqrt(hd7)
          hd_arr1(2) = -hd5 - sqrt(hd7)
          if (hd_arr1(1)<0.0_r2 .or. hd_arr1(2)<0.0_r2) then
             d_lx(3) = maxval(hd_arr1(:))
          else
             d_lx(3) = minval(hd_arr1(:))
          end if
       end if
       
       ! outer
       hd_th1 = grid%co_mx_b( grid%cell_nr2idx(2,nr_cell)    )
       hd1 = tan(hd_th1)**2
       hd2 = hd1*( d_vec(1)**2        +  d_vec(2)**2       ) -  d_vec(3)**2
       hd3 = hd1*( d_vec(1)*p0_vec(1) +  d_vec(2)*p0_vec(2)) -  d_vec(3)*p0_vec(3)
       hd4 = hd1*(p0_vec(1)**2        + p0_vec(2)**2       ) - p0_vec(3)**2
       
       hd5 = hd3/hd2
       hd6 = hd4/hd2
       
       hd7 = hd5**2 - hd6
       if (hd7 < 0.0_r2) then
          d_lx(4) = -1.0_r2
       else
          hd_arr1(1) = -hd5 + sqrt(hd7)
          hd_arr1(2) = -hd5 - sqrt(hd7)
          if (hd_arr1(1)<0.0_r2 .or. hd_arr1(2)<0.0_r2) then
             d_lx(4) = maxval(hd_arr1(:))
          else
             d_lx(4) = minval(hd_arr1(:))
          end if
       end if
       
       ! define sign of step width in <th> direction
       if      (dir_xyz(3)>0.0_r2) then
          d_lx(3) = -abs(d_lx(3))
       else if (dir_xyz(3)<0.0_r2) then
          d_lx(4) = -abs(d_lx(4))
       end if
    else
       d_lx(3:4) = -1.0_r2
    end if


    ! -
    ! 1.3. phi_min, phi_max
    if (grid%n(3) > 1) then
       ! inner
       hd_ph1  = grid%co_mx_c( grid%cell_nr2idx(3,nr_cell) -1 )
       hd1     = tan(hd_ph1)
       d_lx(5) = (p0_vec(2) - p0_vec(1)*hd1) / (d_vec(1)*hd1 - d_vec(2))
       
       ! outer
       hd_ph1  = grid%co_mx_c( grid%cell_nr2idx(3,nr_cell)    )
       hd1     = tan(hd_ph1)
       !print *,hd1
       d_lx(6) = (p0_vec(2) - p0_vec(1)*hd1) / (d_vec(1)*hd1 - d_vec(2))
    else
       d_lx(5:6) = -1.0_r2
    end if
    ! -
    ! 2. determine closest boundary & new new cell number
    ! check: is there any boundary in forward direction of the photon (distance > 0)?
    if (maxval(d_lx(:)) <= 0.0_r2) then
       ! closest boundary NOT found :-(
       kill_photon       = .true.
       kill_photon_count = kill_photon_count +1
   
       if (show_error) then
          print *, "*** problem in sr path():"
          print *, "*** no cell boundary in forward direction."
          
          print *, "pos_xyz"
          print *, real(pos_xyz_new)
          ! -
          print *, "rtp"
          print *, &
               real(norm(pos_xyz_new)), &
               !real(rad2grad(atan2(pos_xyz_new(3),norm(pos_xyz_new(1:2))))), &
               real(rad2grad(atan3(pos_xyz_new(2),pos_xyz_new(1))))
          ! -
          print *, "cell #", grid%cell_nr2idx(:,nr_cell)
          
          print *, "---"
          print *, "check line: 'd_l = d_lx_sel(hi1) * 1.0000001_r2'"
       end if
    end if

    ! determine d_l
    d_lx_sel = 0.0_r2
    p_i1 = 0
    do i1=1, 6
       if (d_lx(i1) > 0.0_r2) then
          p_i1 = p_i1 + 1
          d_lx_sel(p_i1) = d_lx(i1)
       end if
    end do
    hi_arr1 = minloc( d_lx_sel(1:p_i1) )
    hi1     = hi_arr1(1)
    
    
    d_l = d_lx_sel(hi1) +1.0e6_r2*epsilon(d_l)
    
    ! check: is the calculated step width large than the minimum allowed step width?
    if (d_l < grid%d_l_min) then
       ! step width too small
       kill_photon       = .true.
       kill_photon_count = kill_photon_count +1
       d_l = d_l+1.0e5_r2*epsilon(d_l) 
       !print *, 'here'
!~        if (show_error) then
!~           print *, "*** problem in sr path():"
!~           print *, "*** calculated current step width too small"
!~           
!~           print *, "d_l, d_l_min"
!~           print *, d_l, grid%d_l_min
!~        end if
    end if

    ! -
    ! 3.new position
    pos_xyz_new(:) = pos_xyz(:)  +  d_l * dir_xyz(:)
    ! 
    ! 4. new cell number
    !spco        = ca2sp(pos_xyz_new)
    !IF (norm(pos_xyz_new(:)) .gt. grid%co_mx_a(grid%n(1)) ) THEN  
    nr_cell_new = get_cell_nr( grid,pos_xyz_new )
    !ELSE
    !    nr_cell_new = nr_cell
    !    kill_photon = .true.
    !    kill_photon_count = kill_photon_count +1
    !END IF
  end subroutine path_sp

end module transfer_mod

