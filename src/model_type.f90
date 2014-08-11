!----------------------------------------------------------------------------!
! def Model_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE model_type
  
    USE datatype
    USE var_globalnew
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    USE math_mod, ONLY : grad2rad
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Model_TYP
        TYPE(Common_TYP) :: modltype                    ! -----------------    !
        !-----------------------------------------------------------------------!
        REAL(kind=r2)        :: ref_unit
        REAL(kind=r2)        :: ref_unitn
        REAL(kind=r2)        :: r_in
        REAL(kind=r2)        :: r_ou
        REAL(kind=r2)        :: mass
        REAL(kind=r1)        :: r_star
        REAL(kind=r1)        :: M_star
        REAL(kind=r1)        :: t_star
        REAL(kind=r1)        :: l_star
        REAL(kind=r1)        :: distance
        REAL(kind=r1)        :: kep_const
        
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: th_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: ph_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: al_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: zoom_map
        
        INTEGER               :: n_map
        INTEGER               :: n_bin_map
        INTEGER               :: n_r_sub
        INTEGER               :: n_pix_sub
        INTEGER               :: n_star_emi
        
        
    END TYPE Model_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Model_TYP, &
        ! methods
        InitModel, &
        CloseModel, &
        GetModelType, &
        GetModelName, &
        ModelInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitModel(this, ut , un, ref_u, r_in, r_ou, &
                    mass_dust, t_eff, r_star,M_star, n_map, distance, &
                    n_star_emi,th_map,ph_map,zoom_map,al_map,n_bin_map)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Model_TYP)       :: this
        
        INTEGER               :: ut
        INTEGER               :: n_map
        INTEGER               :: n_bin_map
        INTEGER               :: i_map
        INTEGER               :: n_star_emi
        
        CHARACTER(LEN=*)      :: un
        
        REAL(kind=r2)         :: ref_u
        REAL(kind=r2)         :: distance
        REAL(kind=r2)         :: r_in
        REAL(kind=r2)         :: r_ou
        REAL(kind=r2)         :: mass_dust
        REAL(kind=r2)         :: t_eff
        REAL(kind=r2)         :: r_star
        REAL(kind=r2)         :: M_star
        
        REAL(kind=r2),DIMENSION(1:n_map)         :: th_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: ph_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: al_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: zoom_map

        
        !------------------------------------------------------------------------!
        INTENT(IN)          ::  ut,un,ref_u, r_in, r_ou, &
                                mass_dust, t_eff, r_star, M_star, n_map, distance, &
                                th_map,ph_map,zoom_map,al_map,n_bin_map,n_star_emi
        !tbd: add l_star
        INTENT(INOUT)       :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%modltype,ut,un)
        
        this%ref_unit = ref_u
        this%ref_unitn = 1.0_r2/ref_u
        
        this%t_star = t_eff
        this%r_star = r_star*R_sun
        this%M_star = M_star*M_sun
        
        this%kep_const = (con_gamma * this%M_star/con_AU)**0.5  ! in units of m/s
!~         print *, this%kep_const
        this%l_star = 4.0_r2 * PI * this%r_star**2 * con_sigma * this%t_star**4
!~         print *, " => Stellar luminosity [L_sun]: ", this%l_star/L_sun
        ALLOCATE(   this%th_map(1:n_map), &
                    this%ph_map(1:n_map), &
                    this%zoom_map(1:n_map), &
                    this%al_map(1:n_map)   )
                    
        this%th_map     = th_map
        this%ph_map     = ph_map
        this%zoom_map   = zoom_map
        this%al_map     = al_map
        this%distance   = distance
        this%mass       = mass_dust
        this%n_bin_map  = n_bin_map
        this%r_in       = r_in
        this%r_ou       = r_ou
        this%n_map      = n_map
        this%n_star_emi = n_star_emi
        
        this%n_r_sub   = 2          !raytracing: nr of subdiv. in radial direction
        this%n_pix_sub = 11         !raytracing: ... in r/phi dir. in project. on pix.
        
        
        DO i_map=1, this%n_map
            this%th_map(i_map) = grad2rad(th_map(i_map))
            this%ph_map(i_map) = grad2rad(ph_map(i_map))
        END DO
    END SUBROUTINE InitModel


    SUBROUTINE CloseModel(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%modltype)
        DEALLOCATE( this%th_map, &
                    this%ph_map, &
                    this%al_map, &
                    this%zoom_map )
        
        
    END SUBROUTINE CloseModel


    PURE FUNCTION GetModelType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%modltype)
    END FUNCTION GetModelType


    PURE FUNCTION GetModelName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%modltype)
    END FUNCTION GetModelName

    PURE FUNCTION ModelInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Model_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%modltype)
    END FUNCTION ModelInitialized
    


End Module model_type
