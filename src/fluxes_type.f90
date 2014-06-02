!----------------------------------------------------------------------------!
! def fluxes_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE fluxes_type
  
    USE datatype
    USE var_globalnew
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Fluxes_TYP
        TYPE(Common_TYP) :: fltype                       ! ------------------- !
        !----------------------------------------------------------------------!
        CHARACTER(len=1), dimension(1:4)               :: stokes_ext
        
        REAL(kind=r2), dimension(1:4)                  :: stokes, stokes_ini
        !REAL(kind=r2), dimension(:,:,:,:,:), POINTER   :: stokes_map
        REAL(kind=r2), DIMENSION(:,:,:,:), POINTER     :: channel_map
        REAL(kind=r2), DIMENSION(:,:,:), POINTER       :: continuum_map
        INTEGER                :: n_dust_emi
        
    END TYPE Fluxes_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Fluxes_TYP, &
        ! methods
        InitFluxes, &
        CloseFluxes, &
        GetFluxesType, &
        GetFluxesName, &
        FluxesInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitFluxes(this,ut,un, emi_dust, n_bin_map, n_tr, vch,n_lam)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Fluxes_TYP)       :: this
        INTEGER                :: emi_dust
        INTEGER                :: n_bin_map
        INTEGER                :: n_tr
        INTEGER                :: ut
        INTEGER                :: vch
        INTEGER                :: n_lam
        CHARACTER(LEN=*)       :: un
        !------------------------------------------------------------------------!
        INTENT(IN)             :: ut,un, emi_dust, n_bin_map, n_tr, vch,n_lam
        INTENT(INOUT)          :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%fltype,ut,un)

        ! stokes parameters
        this%stokes_ext(1) = "I"
        this%stokes_ext(2) = "Q"
        this%stokes_ext(3) = "U"
        this%stokes_ext(4) = "V"
        
        this%stokes_ini(1)  = 1.0_r2
        this%stokes_ini(2)  = 0.0_r2
        this%stokes_ini(3)  = 0.0_r2
        this%stokes_ini(4)  = 0.0_r2
        
        this%stokes(1)  = 1.0_r2
        this%stokes(2)  = 0.0_r2
        this%stokes(3)  = 0.0_r2
        this%stokes(4)  = 0.0_r2
        
        this%n_dust_emi = emi_dust
        
        
        
        allocate(  &
                  this%channel_map(0:2*n_bin_map, 0:2*n_bin_map,-vch:vch,1:n_tr), &
                  this%continuum_map(0:2*n_bin_map, 0:2*n_bin_map,1:n_lam))
                  
        !this%stokes_map(:,:,:,:,:) = 0.0_r2
        this%channel_map(:,:,:,:)    = 0.0_r2
        this%continuum_map(:,:,:)    = 0.0_r2

    END SUBROUTINE InitFluxes


    SUBROUTINE CloseFluxes(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Fluxes_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%fltype)
        DEALLOCATE(this%channel_map)
        
        
    END SUBROUTINE CloseFluxes

    PURE FUNCTION GetFluxesType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Fluxes_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%fltype)
    END FUNCTION GetFluxesType


    PURE FUNCTION GetFluxesName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Fluxes_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%fltype)
    END FUNCTION GetFluxesName

    PURE FUNCTION FluxesInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Fluxes_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%fltype)
    END FUNCTION FluxesInitialized
    


End Module fluxes_type
