!------------------------------------------------------------------------------!
! def PHOTON_TYP: This type includes all informations of the active photon
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE photon_type
  
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE


    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE PHOTON_TYP
        TYPE(Common_TYP) :: mtype                      ! -----------------     !
        !----------------------------------------------------------------------!
        LOGICAL                                         :: inside
        LOGICAL                                         :: peel_off
        LOGICAL                                         :: observe
        
        INTEGER                                         :: n_interact
        INTEGER                                         :: nr_cell
        INTEGER                                         :: nr_cell_new
        INTEGER                                         :: nr_lam
        
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_li
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_new
        REAL(kind=r2), DIMENSION(1:3)                   :: dir_xyz
        REAL(kind=r2), DIMENSION(1:3,1:3)               :: D
        REAL(kind=r2), DIMENSION(1:4)                   :: stokes
        REAL(kind=r2)                                   :: SINPHI, COSPHI
        REAL(kind=r2)                                   :: SIN2PH, COS2PH
        REAL(kind=r2)                                   :: SINTHE, COSTHE
        
        REAL(kind=r2)                                   :: energy
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE         :: prob_action
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE         :: current_albedo
        

    END TYPE PHOTON_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        PHOTON_TYP, &
        ! methods
        InitPhoton, &
        ClosePhoton, &
        GetPhotonType, &
        GetPhotonName, &
        PhotonInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitPhoton(this, ut, un, n_dust, observe)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP)        :: this
        LOGICAL                 :: observe
        
        INTEGER                 :: ut, n_dust
        
        CHARACTER(LEN=*)        :: un

!~         LOGICAL,OPTIONAL        :: observe
        !----------------------------------------------------------------------!
        INTENT(IN)              :: ut,un ,n_dust, observe
        INTENT(INOUT)           :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        
        this%inside      = .True.
        this%observe     =  observe
        this%peel_off    = .False.

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

    END SUBROUTINE InitPhoton


    SUBROUTINE ClosePhoton(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        
        DEALLOCATE(this%prob_action,this%current_albedo)
        
    END SUBROUTINE ClosePhoton


    PURE FUNCTION GetPhotonType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetPhotonType


    PURE FUNCTION GetPhotonName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetPhotonName

    PURE FUNCTION PhotonInitialized(this) RESULT(i)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(PHOTON_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !--------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
    END FUNCTION PhotonInitialized
    

End Module photon_type
