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
        mo2ca
    
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
        REAL(kind=r2)         :: r_in
            
        INTEGER               :: ut
        INTEGER               :: n_a
        INTEGER               :: n_b
        INTEGER               :: n_c
        INTEGER               :: n_dust
        INTEGER               :: n_lam
        INTEGER               :: egy_lvl       !from gas type
        
        CHARACTER(LEN=*)      :: un
        
        !------------------------------------------------------------------------!
        INTENT(IN)            :: ut, un, model, n_a, sf, n_b, n_c, n_dust, n_lam, egy_lvl
        INTENT(INOUT)         :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%grdtype,ut,un)

        r_in       = model%r_in  !no use yet! ?
        this%n(1)   = n_a
        this%sf     = sf
        this%n(2)   = n_b
        this%n(3)   = n_c
        
        this%n_cell = this%n(1) * this%n(2) * this%n(3)
        print *, "Number of ESCs: ", this%n_cell
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
        
        
        h_r  = sqrt(caco(1)**2+caco(2)**2)
        h_ph = atan3(caco(2),caco(1))
!~         print *, h_r, h_ph
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
        
        !d_z = this%co_mx_c(this%n(3))*2.0_r2/this%n(3)

        !print *,this%co_mx_c(37)
        
!~ 		IF (caco(3) == 0) THEN
!~ 			print *, caco
!~ 			print *, h_r
!~ 		END IF	
        !print *,i_z
        !i_z = int((caco(3)+ this%co_mx_c(this%n(3)))/d_z)+1 
        !stop
        IF ( abs(caco(3)) >=  this%co_mx_c(this%n(3)) ) THEN
            i_z = this%n(3)
        ELSE
            i_z = binary_search(caco(3),this%co_mx_c )
        END IF
!~         print *, i_r, i_ph, i_z

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
!         p_th = 0
!         
!         do
!             if (this%co_mx_b(p_th) > spco(2)) then
!                 i_th = p_th
!                 !print *, this%co_mx_th(p_th)
!                 !print *, 'here'
!                 exit
!             else if (p_th >= this%n(2)) then
!                 ! in the case of rounding errors
!                 ! tbd: check, if this can be avoided
!                 i_th = this%n(2)
!                 exit
!             else
!                 p_th = p_th +1
!                 cycle
!             end if
!         end do
        IF (this%n(2) == 1) THEN
            i_th = 1
        ELSE
            i_th = binary_search(spco(2),this%co_mx_b)
        END IF
        
        !IF (i_th /= test_ph) THEN
        !    print *, i_th , test_ph
        !END IF
        
        
        ! phi
!~         p_ph = 0
!~         do
!~             if (this%co_mx_c(p_ph) >= spco(3) ) then
!~                 i_ph = p_ph
!~                 exit
!~             else if (p_ph >= this%n(3)) then
!~                 ! in the case of rounding errors
!~                 ! tbd: check, if this can be avoided
!~                 i_ph = this%n(3)
!~                 exit
!~             else
!~                 p_ph = p_ph +1
!~                 cycle
!~             end if
!~         end do
        
        IF (this%n(3) == 1) THEN
            i_ph = 1
        ELSE
            i_ph = binary_search(spco(3),this%co_mx_c)
        END IF

    !    ! r

!        p_r = 0
!        do
!            if (p_r == this%n(1)) then
!                i_r = p_r
!                exit
!            end if
!            
!            if (this%co_mx_a(p_r) >= spco(1)) then
!                i_r = p_r
!                exit
!            else
!                p_r = p_r +1
!                cycle
!            end if
!        end do
 
        IF (this%n(1) == 1) THEN
            i_r = 1
        ELSE IF (spco(1) .gt. this%co_mx_a(this%n(1))) THEN
            i_r = this%n(1)  !this should be outside the model space (check this)
        ELSE
            !print *, this%co_mx_a(100), spco(1)
            i_r = binary_search(spco(1),this%co_mx_a)
            !stop
        END IF
        
        !IF (i_r /= test_ph) THEN
        !    print *, i_r , test_ph
        !END IF
        
!~         IF ( spco(1) > this%co_mx_a(this%n(1)) )   THEN
!~             print *, i_r, 'i_r bigger'
!~         END IF
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
