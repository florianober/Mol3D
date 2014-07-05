!----------------------------------------------------------------------------!
! def SOURCE_TYP: This type includes all informations of all sources added
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE source_type
  
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
        TYPE SOURCE_TYP
        !-----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: Luminosity
        REAL(kind=r2),DIMENSION(1:3)                         :: pos_xyz
        
        INTEGER                                              :: i_source

    END TYPE SOURCE_TYP
    
    
    TYPE SOURCES
        TYPE(Common_TYP) :: mtype                       ! -----------------     !
        !-----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: Luminosity
        TYPE(SOURCE_TYP),DIMENSION(:),POINTER                :: source
        
        INTEGER                                              :: n_sources
        INTEGER                                              :: source_init

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
        AddSources
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitSources(this, ut, un, n_sources)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES)          :: this
        INTEGER               :: ut, n_sources
        
        CHARACTER(LEN=*)      :: un
        
        !------------------------------------------------------------------------!
        INTENT(IN)            :: ut,un, n_sources
        INTENT(INOUT)         :: this
        !------------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        this%n_sources = n_sources
        this%source_init  = 0
        
        ALLOCATE(this%source(1:n_sources))

    END SUBROUTINE InitSources

    SUBROUTINE AddSources(this,L,pos_xyz)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(SOURCES)         :: this
        INTEGER               :: i_source
        REAL(kind=r1)         :: L
        REAL(kind=r2),DIMENSION(1:3)         :: pos_xyz
        
        !------------------------------------------------------------------------!
        INTENT(INOUT)         :: this
        INTENT(IN)            :: L, pos_xyz
        !------------------------------------------------------------------------!
        i_source = this%source_init + 1
        IF (i_source > this%n_sources) THEN
			print *, 'ERROR, adding more sources than defined'
			STOP
        END IF
        this%source(i_source)%Luminosity = L
        this%source(i_source)%pos_xyz    = pos_xyz
        this%source(i_source)%i_source   = i_source
        this%source_init  = this%source_init +1
        

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

    PURE FUNCTION SourcesInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(SOURCES), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
          
          IF (i .and. this%n_sources /= this%source_init) THEN
			i = .False.
          END IF
    END FUNCTION SourcesInitialized
    

End Module source_type
