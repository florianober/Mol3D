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
        LOGICAL                                         :: kill
        LOGICAL                                         :: fixed_lam
        
        CHARACTER(1)                                    :: last_interaction_type
        
        INTEGER                                         :: n_interact
        INTEGER                                         :: nr_cell
        INTEGER                                         :: nr_cell_new
        INTEGER                                         :: nr_lam
        INTEGER                                         :: lot_th
        
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
        REAL(kind=r2)                                   :: energy_scale


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
        PhotonInitialized, &
        vecmat
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitPhoton(this, ut, un, i_lam_in)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP)                            :: this
        
        INTEGER                                     :: ut
        INTEGER, INTENT(IN), OPTIONAL               :: i_lam_in
        
        CHARACTER(LEN=*)        :: un
        !----------------------------------------------------------------------!
        INTENT(IN)              :: ut, un
        INTENT(INOUT)           :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        
        this%inside      = .True.
        this%kill        = .False.

        this%n_interact  = 0
        this%nr_cell     = 0
        this%nr_cell_new = 0
        IF (PRESENT(i_lam_in)) THEN
            this%nr_lam = i_lam_in
            this%fixed_lam = .True.
        ELSE
            this%nr_lam  = 0
            this%fixed_lam = .False.
        END IF
        this%lot_th      = 0
        
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
        this%last_interaction_type = 'N'

    END SUBROUTINE InitPhoton


    SUBROUTINE ClosePhoton(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        
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
    
    SUBROUTINE vecmat(this)
        IMPLICIT NONE
        !--------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT)                     :: this
        
        real(kind=r2), dimension(1:3,1:3)                   :: D_help
        REAL(kind=r2)                                       :: R,L,P
        !--------------------------------------------------------------------!

        R = -this%SINTHE * this%SINPHI
        L =  this%SINTHE * this%COSPHI
        P =  this%COSTHE

        this%dir_xyz(2) = -( this%D(1,1) * R  +  this%D(1,2) * L  +  this%D(1,3) * P )
        this%dir_xyz(1) =    this%D(2,1) * R  +  this%D(2,2) * L  +  this%D(2,3) * P
        this%dir_xyz(3) =    this%D(3,1) * R  +  this%D(3,2) * L  +  this%D(3,3) * P

        D_help(1,1) =   this%COSPHI
        D_help(2,1) =   this%SINPHI
        D_help(3,1) =   0.0_r2
        D_help(1,2) = - this%SINPHI*this%COSTHE
        D_help(2,2) =   this%COSPHI*this%COSTHE
        D_help(3,2) = - this%SINTHE
        D_help(1,3) = - this%SINPHI*this%SINTHE
        D_help(2,3) =   this%COSPHI*this%SINTHE
        D_help(3,3) =   this%COSTHE

        this%D = matmul(this%D,D_help)
    END SUBROUTINE vecmat

End Module photon_type
