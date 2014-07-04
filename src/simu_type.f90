!----------------------------------------------------------------------------!
! def Simu_TYP: This type contains all variables needed for the actual simulations
!               These variables may and should change throught simulation
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE simu_type
  
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
    TYPE Simu_TYP
        TYPE(Common_TYP) :: mtype                      ! -----------------     !
        !-----------------------------------------------------------------------!
        LOGICAL                                          :: inside
        
        INTEGER                                          :: n_interact
        INTEGER                                          :: nr_cell
        INTEGER                                          :: nr_cell_new
        INTEGER                                          :: nr_lam
        
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_li
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_new
        REAL(kind=r2), DIMENSION(1:3)                   :: dir_xyz
        REAL(kind=r2), DIMENSION(1:3,1:3)               :: D
        REAL(kind=r2), DIMENSION(1:4)                   :: stokes
        REAL(kind=r2)                                   :: SINPHI, COSPHI
        REAL(kind=r2)                                   :: SIN2PH, COS2PH
        REAL(kind=r2)                                   :: SINTHE, COSTHE
        
        REAL(kind=r2)                                  :: energy
        REAL(kind=r2), DIMENSION(:),POINTER            :: prob_action
        REAL(kind=r2), DIMENSION(:),POINTER            :: current_albedo
        

    END TYPE Simu_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Simu_TYP, &
        ! methods
        InitSimu, &
        CloseSimu, &
        GetSimuType, &
        GetSimuName, &
        SimuInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitSimu(this, ut, un, n_dust, n_lam)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Simu_TYP)        :: this
        
        INTEGER               :: ut, n_dust, n_lam
        
        CHARACTER(LEN=*)      :: un
        
        !------------------------------------------------------------------------!
        INTENT(IN)            :: ut,un ,n_dust , n_lam
        INTENT(INOUT)         :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        
        this%inside      = .true.
        this%n_interact  = 0
        this%nr_cell     = 0
        this%nr_cell_new = 0
        this%nr_lam      = 0
        
        this%pos_xyz     = 0.0_r2
        this%pos_xyz_new = 0.0_r2
        this%pos_xyz_li  = 0.0_r2
        this%dir_xyz     = 0.0_r2
        
        this%D           = 0.0_r2
        this%stokes      = 0.0_r2
        
        this%SINPHI      = 0.0_r2 
        this%COSPHI      = 0.0_r2 
        this%SINTHE      = 0.0_r2 
        this%COSTHE      = 0.0_r2
        this%SIN2PH      = 0.0_r2
        this%COS2PH      = 0.0_r2
        
        this%energy    = 0.0_r2

        ALLOCATE ( this%prob_action( 1:n_dust), &
                   this%current_albedo( 1:n_dust))
        this%prob_action(:)   = 0.0_r2
        this%current_albedo   = 0.0_r2

    END SUBROUTINE InitSimu


    SUBROUTINE CloseSimu(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Simu_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        
        DEALLOCATE(this%prob_action)
        
    END SUBROUTINE CloseSimu


    PURE FUNCTION GetSimuType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Simu_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSimuType


    PURE FUNCTION GetSimuName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Simu_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetSimuName

    PURE FUNCTION SimuInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Simu_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
    END FUNCTION SimuInitialized
    

End Module simu_type
