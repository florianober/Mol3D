!----------------------------------------------------------------------------!
! def Grid_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE Grid_type
  
    USE datatype
    USE var_globalnew
    USE model_type
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE 
    !
    ! multiple coordinates possible:
    !
    ! spherical :  r = / 1 = rho  \
    !                  | 2 = theta|
    !                  \ 3 = phi  /
    ! 
    ! in work:
    ! cylindrical :  r = / 1 = rho  \
    !                    | 2 = theta|
    !                    \ 3 = z    /
    !
    ! not yet implemented:
    ! cartesian : r = / 1 = x \
    !                 | 2 = y |
    !                 \ 3 = z /
    !
    !--------------------------------------------------------------------------!
    TYPE Grid_TYP
        TYPE(Common_TYP) :: grdtype                     ! -----------------    !
        !-----------------------------------------------------------------------!
        REAL(kind=r2)                             :: sf
        REAL(kind=r2)                             :: d_l_min
        
        REAL(kind=r2),DIMENSION(:),POINTER        :: co_mx_a
        REAL(kind=r2),DIMENSION(:),POINTER        :: co_mx_b
        REAL(kind=r2),DIMENSION(:),POINTER        :: co_mx_c
        
        REAL(kind=r1),DIMENSION(:),POINTER        :: cell_minA
        REAL(kind=r1),DIMENSION(:),POINTER        :: cell_gauss_a
        REAL(kind=r1),DIMENSION(:),POINTER        :: cell_gauss_a2
        ! cell_vol(i_cell)  : volume of cell #i_cell [m^3]
        REAL(kind=r2),DIMENSION(:),POINTER        :: cell_vol
        !  Nv(i_cell, i_dust): number of grains of species #i_dust in cell #i_cell
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: Nv       
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: Nv_col        
        REAL(kind=r2),DIMENSION(:),POINTER        :: Nv_mol        
        REAL(kind=r1),DIMENSION(:,:),POINTER      :: lvl_pop        
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: Nv_r
        !  grd_density: this is the number density of the dust component (from mc3d)
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: grd_dust_density
        !  grd_coldensity: this is the number density of all(!) possible collision partner for the
        !                  selected molecule
        !                  first dimension: cell identifier
        !                  second dimension: collision partner
        !                                    1) H2
        !                                    2) H2 para
        !                                    3) H2 ortho
        !                                    4) He
        !                                    5) Elektron
        !                                    6) ...
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: grd_col_density
        !  grd_mol_density: this is the number density of the selected molecule
        !
        REAL(kind=r2),DIMENSION(:),POINTER         :: grd_mol_density
        
        
        REAL(kind=r1), DIMENSION(:,:), POINTER   :: col_finalcolmatrixup
        REAL(kind=r1), DIMENSION(:,:), POINTER   :: col_finalcolmatrixlow
        
        REAL(kind=r1),DIMENSION(:,:),POINTER      :: velo
        REAL(kind=r1),DIMENSION(:),POINTER        :: absvelo
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: cellmidcaco
        REAL(kind=r2),DIMENSION(:,:,:),POINTER    :: i_star_abs
        REAL(kind=r1),DIMENSION(:,:),POINTER      :: t_dust
        REAL(kind=r1),DIMENSION(:,:),POINTER      :: delta_t_dust
        REAL(kind=r1),DIMENSION(:),POINTER        :: t_gas
        REAL(kind=r2), DIMENSION(:),POINTER       :: ddust
        REAL(kind=r2),DIMENSION(:,:),POINTER      :: grd_d_l
        REAL(kind=r2),DIMENSION(1:3)               :: dir_xyz
        !REAL(kind=r2),DIMENSION(1:3)              :: pos_xyz
        
        INTEGER,DIMENSION(1:3)                    :: n
        INTEGER                                   :: n_cell
        INTEGER                                   :: i_cell
        INTEGER                                   :: nh_n_dust
        INTEGER,DIMENSION(:,:,:),POINTER          :: cell_idx2nr
        INTEGER,DIMENSION(:,:),POINTER            :: cell_nr2idx
                
        INTEGER                                   :: counter
        
    END TYPE Grid_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Grid_TYP, &
        ! methods
        InitGrid, &
        CloseGrid, &
        GetGridType, &
        GetGridName, &
        GridInitialized, &
        get_cell_nr, &
        mo2ca, &
        check_inside
    
!~     PRIVATE :: &
!~         get_cell_nr_sp, &
!~         get_cell_nr_cy
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitGrid(this,ut,un, model, n_a, sf, n_b, n_c, n_dust, n_lam, egy_lvl)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP)        :: this
        TYPE(Model_TYP)       :: model  
        
        REAL(kind=r2)         :: sf
        REAL(kind=r2)         :: trash
        REAL(kind=r2),DIMENSION(1:2)         :: d_angle
            
        INTEGER               :: ut
        INTEGER,DIMENSION(1:3)  :: n
        INTEGER               :: i
        INTEGER               :: i_abc
        INTEGER               :: n_a
        INTEGER               :: n_b
        INTEGER               :: n_c
        INTEGER               :: n_dust
        INTEGER               :: n_lam
        INTEGER               :: egy_lvl       !from gas type
        
        CHARACTER(LEN=*)      :: un
        CHARACTER(LEN=256)    :: waste
        
        !------------------------------------------------------------------------!
        INTENT(IN)            :: ut, un, model, n_a, sf, n_b, n_c, n_dust, n_lam, egy_lvl
        INTENT(INOUT)         :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%grdtype,ut,un)
        SELECT CASE(GetGridname(this))
        
        CASE('spherical')
            print *, 'using a spherical grid'
        
            SELECT CASE(GetGridType(this))
            
            CASE(1,2,3)
                ! case 1 logarithm r spaced grid
                ! case 2 linear r spaced grid
                ! we can use the same definitions
                this%n(1)   = n_a
                this%sf     = sf
                this%n(2)   = n_b
                this%n(3)   = n_c
            
            CASE(4)
                ! case 1 logarithm r spaced grid
                ! case 2 linear r spaced grid
                ! we can use the same definitions
                this%n(1)   = 200
                this%sf     = 1.000
                this%n(2)   = 201
                this%n(3)   = 1
            
            CASE(5)
                ! interface to Pluto data
                ! 
                ! was a test and is working, but not in use
                open(unit=1, file="input/grid/pluto_disk.dat", &
                action="read", status="unknown", form="formatted")
                
                DO i = 1,6
                    read(unit=1,fmt=*) waste
                END DO
                read(unit=1,fmt=*) n(1)
                this%n(1) = n(1)
                read(unit=1,fmt=*) n(2)
                read(unit=1,fmt=*) n(3)
                DO i = 1,14
                    read(unit=1,fmt=*) waste
                END DO
                DO i_abc = 1,this%n(1)+2
                    read(unit=1,fmt=*) waste
                END DO
                ! calculate no of th angles
                !
                
                IF (n(2) >1 ) THEN
                    DO i_abc = 1, 2
                        read(unit=1,fmt=*) i, trash, d_angle(i_abc)
                    END DO
                    
                    DO i_abc = 1, n(2)-2
                        read(unit=1,fmt=*) waste
                    END DO
                    read(unit=1,fmt=*) waste
                    read(unit=1,fmt=*) waste
                    
                ELSE
                    read(unit=1,fmt=*) i, trash, d_angle(1)
                    d_angle(2) = PI
                    read(unit=1,fmt=*) waste
                    read(unit=1,fmt=*) waste
                END IF
                ! make sure, this%n(2) is odd!
                this%n(2) = anint(PI/abs(d_angle(2)-d_angle(1)))
                
                
                ! calculate no of ph angles
                IF (n(3) >1 ) THEN
                    DO i_abc = 1, 2
                        read(unit=1,fmt=*) i, trash, d_angle(i_abc)
                    END DO
                ELSE
                    read(unit=1,fmt=*) i, trash, d_angle(1)
                    d_angle(2) = 2.0*PI
                END IF
                ! make sure, this%n(3) is even or 1
                this%n(3) = anint(2.0*PI/abs(d_angle(2)-d_angle(1)))
                
                close(unit=1)
                this%sf     = 1.0
                
            CASE DEFAULT
                print *, 'selected coordinate type not found. try grid_type = [1,2]'
                stop
            END SELECT
            
        CASE('cylindrical')
            print *, 'using a cylindrical grid'
            SELECT CASE(GetGridType(this))
            CASE(1)
                this%n(1)   = n_a
                this%sf     = sf
                this%n(2)   = n_b
                this%n(3)   = n_c
            CASE(2)
                !here we fix grid values for torbias self similar results
                this%n(1)   = 84
                this%sf     = sf
                this%n(2)   = 1
                this%n(3)   = 103
            CASE DEFAULT
                print *, 'selected coordinate type not found. try grid_type = [1]'
                stop
            END SELECT 
            
        CASE('cartesian')
            print *, 'TbD,..initiate grid, cartesian grid is not implemented yet'
            stop
        CASE DEFAULT
            print *, 'selected coordinate system not found.'
            stop
        END SELECT
        
        ! calculate and show final no of cell
        this%n_cell = this%n(1) * this%n(2) * this%n(3)
        print *, "Number of ESCs: ", this%n_cell
        
        ! now allocate all grid dependend properties and intitiate them to zero
        
        ALLOCATE( &
            this%co_mx_a(  0:this%n(1) ), &
            this%co_mx_b( 0:this%n(2) ), &
            this%co_mx_c( 0:this%n(3) ), &
            this%ddust( 1:n_dust), &
            this%cell_vol( 0:this%n_cell ), &
            this%cell_minA( 0:this%n_cell ), &
            this%cell_gauss_a( 0:this%n_cell ), &
            this%cell_gauss_a2( 0:this%n_cell ), &
            this%Nv(       0:this%n_cell, 1:n_dust ), &
            this%Nv_col(       0:this%n_cell, 1:6 ), &
            this%Nv_mol(       0:this%n_cell), &
            this%Nv_r(     0:this%n_cell, 1:n_dust ), &
            this%grd_dust_density( 0:this%n_cell, 1:n_dust ), &
            this%grd_mol_density( 0:this%n_cell), &
            this%grd_col_density( 0:this%n_cell,1:6), &
            this%velo( 0:this%n_cell,1:3 ), &
            this%cellmidcaco(0:this%n_cell,1:3), &
            this%absvelo( 0:this%n_cell), &
            this%cell_idx2nr( 0:this%n(1), 0:this%n(2), 0:this%n(3) ), &
            this%cell_nr2idx( 1:3, 0:this%n_cell ), &
            this%t_dust(  0:this%n_cell, 1:n_dust ), &
            this%delta_t_dust(  0:this%n_cell, 1:n_dust ), &
            this%t_gas(  0:this%n_cell), &
            this%grd_d_l( 1:this%n_cell, 1:n_lam ) , &
            this%i_star_abs( 1:n_dust, 1:n_lam, 1:this%n_cell), &
            this%lvl_pop( 0:this%n_cell, 1:egy_lvl ) )
            
        this%nh_n_dust = n_dust
        this%co_mx_a(:)           = 0.0_r2
        this%co_mx_b(:)           = 0.0_r2
        this%co_mx_c(:)           = 0.0_r2
        this%cell_vol(:)          = 0.0_r2
        this%cell_minA(:)         = 0.0_r2
        this%cell_gauss_a(:)      = 0.0
        this%cell_gauss_a2(:)     = 0.0
        this%cellmidcaco(:,:)     = 0.0_r2
        this%ddust(:)             = 0.0_r2
        this%Nv(:,:)              = 0.0_r2
        this%Nv_mol(:)            = 0.0_r2
        this%Nv_col(:,:)          = 0.0_r2
        this%Nv_r(:,:)            = 0.0_r2
        this%grd_dust_density(:,:)= 0.0_r2
        this%grd_col_density(:,:) = 0.0_r2
        this%grd_mol_density(:)   = 0.0_r2
        this%velo(:,:)            = 0.0_r2
        this%absvelo(:)           = 0.0_r2
        this%cell_idx2nr(:,:,:)   = 0
        this%cell_nr2idx(:,:)     = 0
        this%t_dust(:,:)          = 0.0
        this%delta_t_dust(:,:)    = 0.0
        this%t_gas(:)             = 0.0
        this%grd_d_l(:,:)         = 0.0_r2
        this%i_star_abs(:,:,:)    = 0.0_r2
        this%d_l_min              = 0.0_r2
        this%dir_xyz(:)           = 0.0_r2
        this%lvl_pop(:,:)         = 0.0
        
        this%counter            = 0
    END SUBROUTINE InitGrid

    
    SUBROUTINE CloseGrid(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%grdtype)
        DEALLOCATE( &
            this%co_mx_a, &
            this%co_mx_b, &
            this%co_mx_c, &
            this%cell_vol, &
            this%cell_minA, &
            this%cell_gauss_a, &
            this%cell_gauss_a2, &
            this%ddust, &
            this%Nv, &
            this%Nv_col, &
            this%Nv_mol, &
            this%Nv_r, &
            this%grd_dust_density, &
            this%grd_col_density, &
            this%grd_mol_density, &
            this%cellmidcaco, &
            this%velo, &
            this%absvelo, &
            this%cell_idx2nr, &
            this%cell_nr2idx, &
            this%t_dust, &
            this%delta_t_dust, &
            this%t_gas, &
            this%grd_d_l, &
            this%lvl_pop, &
            this%i_star_abs)
        
        
    END SUBROUTINE CloseGrid


    PURE FUNCTION GetGridType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%grdtype)
    END FUNCTION GetGridType


    PURE FUNCTION GetGridName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%grdtype)
    END FUNCTION GetGridName

    PURE FUNCTION GridInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Grid_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%grdtype)
    END FUNCTION GridInitialized

    FUNCTION mo2ca(grid,moco) RESULT(caco)
       ! generic function for all geometrics
        USE datatype
        USE var_globalnew
        USE math_mod, ONLY : sp2ca, cy2ca
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                  :: grid
        !------------------------------------------------------------------------!
        
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)   :: moco
        REAL(kind=r2), DIMENSION(1:3)               :: caco
        !------------------------------------------------------------------------!
        SELECT CASE(GetGridName(grid))
        
        CASE('spherical')
            caco = sp2ca(moco)
        CASE('cylindrical')
            caco = cy2ca(moco)
        CASE('cartesian')
            caco = moco
 
        CASE DEFAULT
            print *, 'selected coordinate system not found, mo2ca'
            stop
        END SELECT
        
    
    END FUNCTION mo2ca
    
    ! ################################################################################################
    ! check: is current position of photon stil inside the model space?
    ! ---
    FUNCTION check_inside(caco,this,model) result(check_inside_result)
    USE math_mod, ONLY : norm
    
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Grid_TYP), INTENT(IN)                  :: this
    TYPE(Model_TYP), INTENT(IN)                 :: model
    
    REAL(kind=r2), dimension(3), intent(in)     :: caco
    
    LOGICAL                                     :: check_inside_result   
    !------------------------------------------------------------------------!
    SELECT CASE(GetGridName(this))
        
        CASE('spherical')
            check_inside_result = (norm(caco) <= model%r_ou)
            
        CASE('cylindrical')
            check_inside_result = (sqrt( dot_product( caco(1:2), caco(1:2) ) ) <= model%r_ou &
                                  .and. abs(caco(3)) <= model%r_ou )
        
        CASE('cartesian')
            check_inside_result = (all(abs(caco) <= model%r_ou))
        CASE DEFAULT
            print *, 'selected coordinate system not found, set_boundaries'
            stop
    END SELECT
    
    END FUNCTION check_inside
    
    
    FUNCTION get_cell_nr(this,caco) RESULT(get_cell_nr_result)
        ! generic function for all geometrics
        USE datatype
        USE var_globalnew
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                  :: this
        !------------------------------------------------------------------------!
        
        INTEGER                                     :: get_cell_nr_result 
        
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)  :: caco
        !------------------------------------------------------------------------!
        
        SELECT CASE(GetGridName(this))
        
        CASE('spherical')
            get_cell_nr_result = get_cell_nr_sp(this,caco)
            
        CASE('cylindrical')
            get_cell_nr_result = get_cell_nr_cy(this,caco)
        CASE('cartesian')
            print *, 'TbD, not finished yet, cell number'
            stop
            
        CASE DEFAULT
            print *, 'selected coordinate system not found'
            stop
        END SELECT
    
    
    END FUNCTION get_cell_nr
    
    FUNCTION get_cell_nr_cy(this,caco) RESULT(get_cell_nr_result)
        USE datatype
        USE var_globalnew
        USE math_mod, ONLY : atan3, binary_search
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                  :: this
        !------------------------------------------------------------------------!
        
        INTEGER                                     :: i_r
        INTEGER                                     :: i_ph
        INTEGER                                     :: i_z
        INTEGER                                     :: get_cell_nr_result 
        
        REAL(kind=r2), DIMENSION(1:3),INTENT(in)  :: caco
        REAL(kind=r2)  :: h_r, h_ph, d_ph
        !------------------------------------------------------------------------!
        
!~         print *,'starting'
        h_r  = sqrt(caco(1)**2+caco(2)**2)
        h_ph = atan3(caco(2),caco(1))

        ! 1.1 get i_r
        
        ! --
        !dr1 = (this%co_mx_a(this%n(1)) - this%co_mx_a(0) ) * (this%sf-1.0_r2)/ (this%sf**this%n(1) - 1.0_r2)
        !i_r = floor( log(1.0_r2 + (this%sf-1.0_r2)*(h_r-this%co_mx_a(0))/dr1) / log(this%sf) )  + 1
        
        
        IF (this%n(1) == 1) THEN
            i_r = 1
        ELSE IF (h_r < this%co_mx_a(0) )THEN
            i_r = 0
        ELSE IF (h_r >= this%co_mx_a(this%n(1))) THEN
            i_r = this%n(1)  !this should be outside the model space (check this)
        ELSE
            !print *, this%co_mx_a(100), spco(1)
            i_r = binary_search(h_r,this%co_mx_a)
            !stop
        END IF


        ! 1.2 get i_ph
        ! --
        d_ph = PI*2.0_r2/this%n(2)
        i_ph = int(h_ph/d_ph)+1      
        
        IF  ( h_ph .ge. PI*2.0_r2 ) THEN
            i_ph = this%n(2)
        END IF


        ! 1.3 get i_z
        ! --

        IF ( abs(caco(3)) >=  this%co_mx_c(this%n(3)) ) THEN
            i_z = this%n(3)
        ELSE
            i_z = binary_search(caco(3),this%co_mx_c )
        END IF

        
        get_cell_nr_result = this%cell_idx2nr(i_r,i_ph,i_z)
    END FUNCTION get_cell_nr_cy
    
  ! ################################################################################################
  ! determine cell number from sperical coordinates. still used from mc3d, should be rewritten
  ! ---
    FUNCTION get_cell_nr_sp(this,caco) RESULT(get_cell_nr_result)
        USE datatype
        USE var_globalnew
        USE math_mod, ONLY : ca2sp, binary_search
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Grid_TYP), INTENT(IN)                  :: this
        !------------------------------------------------------------------------!
        
        INTEGER                                     :: i_r
        INTEGER                                     :: i_th
        INTEGER                                     :: i_ph
        INTEGER                                     :: get_cell_nr_result 
        
        REAL(kind=r2), DIMENSION(1:3)               :: spco
        REAL(kind=r2), DIMENSION(1:3),INTENT(in)    :: caco
        !------------------------------------------------------------------------!
        
        spco = ca2sp(caco)
        ! ---
        ! theta

        IF (this%n(2) == 1) THEN
            i_th = 1
        ELSE
            i_th = binary_search(spco(2),this%co_mx_b)
        END IF
        IF ( i_th .gt. this%n(2)) THEN
			i_th = this%n(2)			
        END IF
        
        ! ---
        ! phi
        
        IF (this%n(3) == 1) THEN
            i_ph = 1
        ELSE
            i_ph = binary_search(spco(3),this%co_mx_c)
        END IF
        
        
        ! ---
        ! r
        IF (this%n(1) == 1) THEN
            i_r = 1
        ELSE IF (spco(1) .gt. this%co_mx_a(this%n(1))) THEN
            i_r = this%n(1)  !this should be outside the model space (check this)
        ELSE
            !print *, this%co_mx_a(100), spco(1)
            i_r = binary_search(spco(1),this%co_mx_a)
            !stop
        END IF
        
        get_cell_nr_result = this%cell_idx2nr( i_r, i_th, i_ph )

    END FUNCTION get_cell_nr_sp
    
!  ! ################################################################################################
!  ! determine cell number from cartesian coordinate
!  ! using an analytical expression;
!  ! only applicable in grids - without local subgrids
!  !                          - where cell boundaries in r direction are described by step factor sf
!  ! ---
!  function get_cell_nr2(caco) result(get_cell_nr2_result)
!    use datatype
!    use var_global
!    use math_mod

!    real(kind=r2), dimension(:), intent(in) :: caco
!    integer :: get_cell_nr2_result
!    integer :: i_r, i_th, i_ph
!    real(kind=r2) :: hd_r, hd_th, hd_ph
!    ! ---
!    ! r
!    hd_r = norm(caco)
!    if (hd_r <= r_ou) then
!       i_r = floor( log(1.0_r2 + (sf-1.0_r2)*(hd_r-r_in)/dr1) / log(sf) )  + 1
!    else
!       i_r = n(1)
!    end if

!    ! th
!    if (n(2) > 1) then
!       hd_th = atan2( caco(3), sqrt(caco(1)**2 + caco(2)**2))
!       i_th = floor( real(n(2),kind=r2) * (hd_th + PI2)/PI ) +1
!       if (i_th > n(2)) then
!          i_th = n(2)
!       end if
!    else
!       i_th = 1
!    end if

!    ! ph
!    if (n(3) > 1) then
!       hd_ph = atan3( caco(2), caco(1) )
!       i_ph = floor( real(n(3),kind=r2) * hd_ph/PIx2       ) +1
!       if (i_ph > n(3)) then
!          i_ph = n(3)
!       end if
!    else
!       i_ph = 1
!    end if

!    get_cell_nr2_result = cell_idx2nr( i_r, i_th, i_ph )

!  end function get_cell_nr2

End Module Grid_type
