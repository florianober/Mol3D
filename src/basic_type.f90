!----------------------------------------------------------------------------!
! def basic_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE basic_type
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE
    ! ---
    ! special variables
    ! ---
    !              ( 1 0 0 )                      ( 1 0 0 0 )
    ! mat_ident3 = ( 0 1 0 )         mat_ident4 = ( 0 1 0 0 )
    !              ( 0 0 1 )                      ( 0 0 1 0 )
    !                                             ( 0 0 0 1 )
    ! PIx2  = PI * 2
    ! PIx4  = PI * 4
    ! PIx34 = PI * 3/4
    ! PI2   = PI / 2
    ! PI4   = PI / 4
    ! PI2x4 = PI^2 * 4

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Basic_TYP
        TYPE(Common_TYP) :: mtype                       ! -----------------    !
        !----------------------------------------------------------------------!
        INTEGER          :: pnamelen
        INTEGER          :: concept_ps
        INTEGER          :: num_core
        
        REAL(kind=r2)    :: PIx2
        REAL(kind=r2)    :: PIx4
        REAL(kind=r2)    :: PIx34
        REAL(kind=r2)    :: PI2
        REAL(kind=r2)    :: PI4
        REAL(kind=r2)    :: PI2x4
        
        REAL(kind=r1)    :: linescale
        
        ! store all variables from file
        INTEGER          :: n_tem
        
        REAL(kind=r1)    :: d_tem
        REAL(kind=r1)    :: t_dust_max
        REAL(kind=r1)    :: t_dust_min
        

        REAL(kind=r1), dimension(1:3,1:3) :: mat_ident3
        REAL(kind=r1), dimension(1:4,1:4) :: mat_ident4
        
        
        CHARACTER(len=256)                :: pronam
        CHARACTER(len=256)                :: pronam_old
        CHARACTER(len=256)                :: path_results
        LOGICAL                           :: MCRT
        LOGICAL                           :: calc_tmp
        LOGICAL                           :: old_model
        LOGICAL                           :: pluto_data
        LOGICAL                           :: do_raytr
        LOGICAL                           :: do_continuum_map
        LOGICAL                           :: do_velo_ch_map

    END TYPE Basic_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Basic_TYP, &
        ! methods
        InitBasic, &
        CloseBasic, &
        GetSimuType, &
        GetName, &
        BasicInitialized, &
        Getproname
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitBasic(this, ut, un, pname, presult, concept, calc_tmp,old_pname, old_model, &
                            do_raytr, do_continuum_map,do_velo_ch_map,n_tem, t_dust_min, t_dust_max, num_core,pluto_data)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP)       :: this
        
        INTEGER               :: ut
        INTEGER               :: concept
        INTEGER               :: n_tem
        INTEGER               :: num_core
        
        CHARACTER(LEN=*)      :: un
        CHARACTER(LEN=*)      :: pname
        CHARACTER(LEN=*)      :: old_pname
        CHARACTER(LEN=*)      :: presult
        
        REAL(kind=r1)         :: t_dust_max
        REAL(kind=r1)         :: t_dust_min
        
        LOGICAL               :: calc_tmp
        LOGICAL               :: old_model
        LOGICAL               :: pluto_data
        LOGICAL               :: do_raytr
        LOGICAL               :: do_continuum_map
        LOGICAL               :: do_velo_ch_map

        !----------------------------------------------------------------------!
        INTENT(IN)          :: ut,un,presult, pname, concept, num_core, &
                                n_tem, t_dust_max, t_dust_min, calc_tmp,old_model,old_pname, &
                                pluto_data,do_continuum_map,do_velo_ch_map
        INTENT(INOUT)       :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
    
        this%PIx2  = PI * 2.0_r2
        this%PIx4  = PI * 4.0_r2
        this%PIx34 = PI * 3.0_r2/4.0_r2
        this%PI2   = PI / 2.0_r2
        this%PI4   = PI / 4.0_r2
        this%PI2x4 = PI**2 * 4.0_r2
        
        this%num_core = num_core
        this%linescale = con_h*con_c/(4.0*PI*sqrt(PI))

        ! identity matrix [3x3]
        this%mat_ident3(1,1) = 1.0
        this%mat_ident3(1,2) = 0.0
        this%mat_ident3(1,3) = 0.0
           
        this%mat_ident3(2,1) = 0.0
        this%mat_ident3(2,2) = 1.0
        this%mat_ident3(2,3) = 0.0
           
        this%mat_ident3(3,1) = 0.0
        this%mat_ident3(3,2) = 0.0
        this%mat_ident3(3,3) = 1.0
        
        ! identity matrix [4x4]      
        this%mat_ident4(1,1) = 1.0_r2
        this%mat_ident4(1,2) = 0.0_r2
        this%mat_ident4(1,3) = 0.0_r2
        this%mat_ident4(1,4) = 0.0_r2
           
        this%mat_ident4(2,1) = 0.0
        this%mat_ident4(2,2) = 1.0
        this%mat_ident4(2,3) = 0.0
        this%mat_ident4(2,4) = 0.0
           
        this%mat_ident4(3,1) = 0.0
        this%mat_ident4(3,2) = 0.0
        this%mat_ident4(3,3) = 1.0
        this%mat_ident4(3,4) = 0.0
           
        this%mat_ident4(4,1) = 0.0
        this%mat_ident4(4,2) = 0.0
        this%mat_ident4(4,3) = 0.0
        this%mat_ident4(4,4) = 1.0


        this%MCRT       = .false.
        this%n_tem = n_tem
        this%t_dust_max = t_dust_max
        this%t_dust_min = t_dust_min
        this%calc_tmp   = calc_tmp
        this%old_model    = old_model
        this%pluto_data    = pluto_data
        this%do_raytr   = do_raytr
        this%do_continuum_map   = do_continuum_map
        this%do_velo_ch_map = do_velo_ch_map
        
        this%d_tem = (this%t_dust_max - this%t_dust_min) / (real(this%n_tem, kind=r2)-1)
        
        this%pronam = pname    
        this%pronam_old = old_pname    
        this%path_results = presult
        this%pnamelen = LEN_TRIM(pname)
        this%concept_ps = concept
        
    END SUBROUTINE InitBasic


    SUBROUTINE CloseBasic(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
    END SUBROUTINE CloseBasic


    PURE FUNCTION GetSimuType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSimuType


    PURE FUNCTION GetName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetName

    PURE FUNCTION BasicInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Basic_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
    END FUNCTION BasicInitialized
    
    FUNCTION Getproname(this) RESULT(pname)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Basic_TYP), INTENT(IN) :: this
          CHARACTER(len=this%pnamelen) :: pname
          !------------------------------------------------------------------------!
          pname = this%pronam
    END FUNCTION Getproname
    

End Module basic_type
