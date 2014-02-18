!----------------------------------------------------------------------------!
! basic data and methods common to all objects 
! inspired by fosite by T. Illenseer
!----------------------------------------------------------------------------!
MODULE common_type

  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! common data structure
  !--------------------------------------------------------------------------!
  TYPE Common_TYP
     PRIVATE
     INTEGER           :: type
     CHARACTER(LEN=32) :: name
     INTEGER           :: error               ! error code                   !
     LOGICAL           :: init = .FALSE.      ! init status                  !
  END TYPE Common_TYP
  SAVE
  ! these variables should be the same for all objects
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Common_TYP, &
       ! methods
       InitCommon, &
       CloseCommon, &
       GetType, &
       GetName, &
       Initialized
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitCommon(this,t,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP)  :: this
    INTEGER           :: t
    CHARACTER(LEN=*)  :: n
    !------------------------------------------------------------------------!
    INTENT(IN)        :: t,n
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    this%type = t
    this%name = n
    this%error = 0
    this%init  = .TRUE.
  END SUBROUTINE InitCommon

  SUBROUTINE CloseCommon(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    this%init = .FALSE.
  END SUBROUTINE CloseCommon


  PURE FUNCTION GetType(this) RESULT(t)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    INTEGER :: t
    !------------------------------------------------------------------------!
    t = this%type
  END FUNCTION GetType


  PURE FUNCTION GetName(this) RESULT(n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: n
    !------------------------------------------------------------------------!
    n = this%name
  END FUNCTION GetName

  PURE FUNCTION Initialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = this%init
  END FUNCTION Initialized


END MODULE common_type
