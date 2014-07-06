!----------------------------------------------------------------------------!
! def SOURCE_TYP: This type includes all informations of all sources added
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE source_type
  
    USE datatype
    USE var_globalnew
    USE dust_type
    USE math_mod
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE


    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
        TYPE SOURCE_TYP
        !-----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: Luminosity
        REAL(kind=r2),DIMENSION(1:3)                         :: pos_xyz
        REAL(kind=r2),DIMENSION(:),POINTER                   :: wave_cpf
        INTEGER                                              :: s_type

    END TYPE SOURCE_TYP
    
    
    TYPE SOURCES
        TYPE(Common_TYP) :: mtype                       ! -----------------     !
        !-----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: L_total
        TYPE(SOURCE_TYP),DIMENSION(:),POINTER                :: source
        REAL(kind=r2),DIMENSION(:),POINTER                   :: source_cpf
        REAL(kind=r2),DIMENSION(:),POINTER                   :: lam
        
        INTEGER                                              :: n_sources
        INTEGER                                              :: n_lam

    END TYPE SOURCES
    
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        SOURCES, &
        ! methods
        InitSources, &
        CloseSources, &
        GetSourcesType, &
        GetSourcesName, &
        SourcesInitialized, &
        AddSources, &
        GetNewSource, &
        GetNewLam
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitSources(this, ut, un,dust)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES)          :: this
        TYPE(Dust_TYP)         :: dust
        INTEGER                :: ut
        
        CHARACTER(LEN=*)       :: un
        
        !------------------------------------------------------------------------!
        INTENT(IN)             :: ut,un
        INTENT(INOUT)          :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        this%n_sources = 0
        ALLOCATE(this%lam(1:dust%n_lam))
        this%n_lam = dust%n_lam
        
        this%lam   = dust%lam
        
        
    END SUBROUTINE InitSources

    SUBROUTINE AddSources(this,T_star,R_star,pos_xyz,s_type)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES)         :: this
        INTEGER               :: i_source, s_type, i
        REAL(kind=r2)         :: L
        REAL(kind=r1)         :: R_star,T_star
        REAL(kind=r2),DIMENSION(1:3)               :: pos_xyz
        TYPE(SOURCE_TYP),DIMENSION(:), ALLOCATABLE :: source_tmp
        
        !------------------------------------------------------------------------!
        INTENT(INOUT)         :: this
        INTENT(IN)            :: pos_xyz, s_type, R_star,T_star
        !------------------------------------------------------------------------!
        i_source = this%n_sources + 1
        
        IF (this%n_sources > 0) THEN
            ALLOCATE(source_tmp(1:i_source))
            source_tmp(1:this%n_sources) = this%source
            DEALLOCATE(this%source)
            DEALLOCATE(this%source_cpf)
            ALLOCATE(this%source(1:i_source), &
                     this%source_cpf(1:i_source))
            this%source_cpf = 0.0_r2
            this%source = source_tmp
            DEALLOCATE(source_tmp)
            
        ELSE
            ALLOCATE(this%source(1:i_source), &
                     this%source_cpf(1:i_source))
        END IF
        
        IF (s_type == 1) THEN 
            ! source is a point star
            !
!~             L = con_sigma*T_star**4*PI*4.0_r2*R_star**2
            ALLOCATE(this%source(i_source)%wave_cpf(1:this%n_lam))
            this%source(i_source)%wave_cpf(:) = 0.0_r2
            L = 0.0_r2
            L = integ1(this%lam(:),PI * 4.0_r2 * PI * R_star**2 * planck(T_star,this%lam(:)) , 1, this%n_lam)
            DO i = 1,this%n_lam
                
                IF ( i == 1) THEN
                    this%source(i_source)%wave_cpf(i) = 0.0_r2 
                ELSE 
                    this%source(i_source)%wave_cpf(i) = this%source(i_source)%wave_cpf(i-1) + &
                                integ1(this%lam(:),PI * 4.0_r2 * PI * R_star**2 * planck(T_star,this%lam(:)) , i-1, i) /&
                                L
                END IF
            END DO
            this%source(i_source)%pos_xyz     = pos_xyz
            this%source(i_source)%s_type      = s_type
            this%n_sources                    = i_source
            
            this%source(i_source)%Luminosity  = L
        ELSE
            print *, 'ERROR: source type not implemented'
            stop
        END IF
!~         print *, L
       
        
        ! update source cpf
        this%L_total = sum(this%source(:)%Luminosity)
         
        DO i = 1,i_source
            IF (i == 1) THEN
                this%source_cpf(i) = this%source(i)%Luminosity /this%L_total
            ELSE
                this%source_cpf(i) = this%source_cpf(i-1) + this%source(i)%Luminosity /this%L_total
            END IF
        END DO
        
        

    END SUBROUTINE AddSources    
    

    SUBROUTINE CloseSources(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        DEALLOCATE(this%source)
        
    END SUBROUTINE CloseSources


    PURE FUNCTION GetSourcesType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)      :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSourcesType


    PURE FUNCTION GetSourcesName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)       :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetSourcesName
    
    FUNCTION GetNewSource(this,rndx) RESULT(i)
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER                   :: i
        !------------------------------------------------------------------------!
        IF (this%n_sources == 1) THEN
            i = 1
        ELSE
            i = binary_search(rndx, this%source_cpf) +1
        END IF
        
    
    END FUNCTION GetNewSource
    
    FUNCTION GetNewLam(this,i_source,rndx) RESULT(i)
        
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER, INTENT(IN)       :: i_source
        INTEGER                   :: i
        !------------------------------------------------------------------------!

        i = MIN(binary_search(rndx, this%source(i_source)%wave_cpf) +1,this%n_lam)
        
    
    END FUNCTION GetNewLam


    PURE FUNCTION SourcesInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(SOURCES), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
          
          IF (i .and. this%n_sources == 0) THEN
            i = .False.
          END IF
    END FUNCTION SourcesInitialized
    

End Module source_type
