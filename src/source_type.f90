!----------------------------------------------------------------------------!
! def SOURCE_TYP: This type includes all informations of all sources added
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE source_type
    USE datatype
    USE var_global
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
        !----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: Luminosity
        REAL(kind=r2),DIMENSION(1:3)                         :: pos_xyz
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: wave_pdf
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: wave_cdf
        INTEGER                                              :: s_type

    END TYPE SOURCE_TYP

    TYPE SOURCES
        TYPE(Common_TYP) :: mtype
        !----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: L_total
        TYPE(SOURCE_TYP),DIMENSION(:),ALLOCATABLE            :: source
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: source_cdf
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: lam
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: d_lam
        
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

    SUBROUTINE InitSources(this, ut, un, dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES)          :: this
        TYPE(Dust_TYP)         :: dust
        INTEGER                :: ut
        
        CHARACTER(LEN=*)       :: un
        
        !----------------------------------------------------------------------!
        INTENT(IN)             :: ut,un
        INTENT(INOUT)          :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        this%n_sources = 0
        ALLOCATE(this%lam(1:dust%n_lam))
        ALLOCATE(this%d_lam(1:dust%n_lam))
        this%n_lam = dust%n_lam
        
        this%lam   = dust%lam
        this%d_lam   = dust%d_lam

    END SUBROUTINE InitSources

    SUBROUTINE AddSources(this, s_type, pos_xyz, T_star, R_star, L_star)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES)                              :: this
        INTEGER                                    :: i_source, s_type, i
        REAL(kind=r2)                              :: L, B
        REAL(kind=r1),OPTIONAL                     :: R_star, L_star
        REAL(kind=r1)                              :: T_star
        REAL(kind=r2),DIMENSION(1:3)               :: pos_xyz
        TYPE(SOURCE_TYP),DIMENSION(:), ALLOCATABLE :: source_tmp

        !----------------------------------------------------------------------!
        INTENT(INOUT)         :: this
        INTENT(IN)            :: pos_xyz, s_type, R_star,T_star
        !----------------------------------------------------------------------!
        i_source = this%n_sources + 1
        B = 0.0_r2
        L = 0.0_r2

        IF (this%n_sources > 0) THEN
            ALLOCATE(source_tmp(1:i_source))
            source_tmp(1:this%n_sources) = this%source
            DEALLOCATE(this%source)
            DEALLOCATE(this%source_cdf)
            ALLOCATE(this%source(1:i_source), &
                     this%source_cdf(1:i_source))
            this%source_cdf = 0.0_r2
            this%source(1:this%n_sources) = source_tmp
            DEALLOCATE(source_tmp)

        ELSE
            ALLOCATE(this%source(1:i_source), &
                     this%source_cdf(1:i_source))
        END IF

        IF (s_type == 1) THEN 
            ! source is a point star
            !
            ALLOCATE(this%source(i_source)%wave_cdf(1:this%n_lam))
            ALLOCATE(this%source(i_source)%wave_pdf(1:this%n_lam))
            this%source(i_source)%wave_cdf(:) = 0.0_r2
            this%source(i_source)%wave_pdf(:) = 0.0_r2
            B = integ1(this%lam(:),planck(T_star,this%lam(:)) , 1, this%n_lam)
            this%source(i_source)%wave_pdf(:) = planck(T_star,this%lam(:))/B
            DO i = 1, this%n_lam
                IF ( i == 1) THEN
                    this%source(i_source)%wave_cdf(i) =  this%source(i_source)%wave_pdf(i)
                ELSE
                    this%source(i_source)%wave_cdf(i) = this%source(i_source)%wave_cdf(i-1) + &
                                integ1(this%lam(:),this%source(i_source)%wave_pdf(:) , i-1, i)
                END IF
            END DO

            IF (present(L_star)) THEN
                L = L_star
            ELSE
                IF (present(R_star)) THEN
                    L = integ1(this%lam(:),PI * 4.0_r2 * PI * R_star**2 *      &
                        planck(T_star,this%lam(:)) , 1, this%n_lam)
                ELSE
                    print *,'ERROR: star radius not given'
                    return
                END IF
            END IF

        ELSE IF (s_type == 2) THEN
            ! source is a point source, but given is the luminosity and Temperature
            !
            ALLOCATE(this%source(i_source)%wave_cdf(1:this%n_lam))
            ALLOCATE(this%source(i_source)%wave_pdf(1:this%n_lam))
            this%source(i_source)%wave_cdf(:) = 0.0_r2
            this%source(i_source)%wave_pdf(:) = 0.0_r2
            
            B = integ1(this%lam(:),planck(T_star,this%lam(:)) , 1, this%n_lam)
            this%source(i_source)%wave_pdf(:) = planck(T_star,this%lam(:))/B
            DO i = 1,this%n_lam
                IF ( i == 1) THEN
                    this%source(i_source)%wave_cdf(i) = 0.0_r2 
                ELSE 
                    this%source(i_source)%wave_cdf(i) = this%source(i_source)%wave_cdf(i-1) + &
                                integ1(this%lam(:),this%source(i_source)%wave_pdf(:) , i-1, i)
                END IF
            END DO
            IF (present(L_star)) THEN
                L = L_star
            ELSE
                print *,'ERROR: star luminosity not given'
                return
            END IF
        ELSE
            print *, 'ERROR: source type not implemented'
            stop
        END IF
        this%source(i_source)%pos_xyz     = pos_xyz
        this%source(i_source)%s_type      = s_type
        this%n_sources                    = i_source
        this%source(i_source)%Luminosity  = L
        ! update source cdf
        this%L_total = sum(this%source(:)%Luminosity)
        DO i = 1,i_source
            IF (i == 1) THEN
                this%source_cdf(i) = this%source(i)%Luminosity /this%L_total
            ELSE
                this%source_cdf(i) = this%source_cdf(i-1) +                    &
                                     this%source(i)%Luminosity /this%L_total
            END IF
        END DO

    END SUBROUTINE AddSources


    SUBROUTINE CloseSources(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        DEALLOCATE(this%source)
        
    END SUBROUTINE CloseSources


    PURE FUNCTION GetSourcesType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)      :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSourcesType


    PURE FUNCTION GetSourcesName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)       :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetSourcesName
    
    FUNCTION GetNewSource(this,rndx) RESULT(i)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER                   :: i
        !----------------------------------------------------------------------!
        IF (this%n_sources == 1) THEN
            i = 1
        ELSE
            i = binary_search(rndx, this%source_cdf)+1
        END IF

    END FUNCTION GetNewSource
    
    FUNCTION GetNewLam(this,i_source,rndx) RESULT(i)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER, INTENT(IN)       :: i_source
        INTEGER                   :: i
        !----------------------------------------------------------------------!
        i = binary_search(rndx, this%source(i_source)%wave_cdf) + 1 
        
    END FUNCTION GetNewLam

    PURE FUNCTION SourcesInitialized(this) RESULT(i)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        LOGICAL :: i
        !----------------------------------------------------------------------!
        i = Initialized_common(this%mtype)
          
        IF (i .and. this%n_sources == 0) THEN
            i = .False.
        END IF
    END FUNCTION SourcesInitialized

END MODULE source_type
