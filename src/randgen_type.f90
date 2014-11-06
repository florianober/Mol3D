!----------------------------------------------------------------------------!
! def randgen_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE randgen_type
  
  USE datatype
  USE common_type, &
       GetType_common => GetType, GetName_common => GetName, &
       Initialized_common => Initialized
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! 
  !--------------------------------------------------------------------------!

  INTEGER, PARAMETER, PRIVATE       :: M=714025, IA=1366, IC=150889

  TYPE Randgen_TYP
     TYPE(Common_TYP)               :: gentype            !   RAN2 ...       !
     REAL(kind=r2)                  :: RM, rndx
     INTEGER                        :: IDUM, IY, IFF
     INTEGER, dimension(97)        :: IR
  END TYPE Randgen_TYP
  SAVE

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Randgen_TYP, &
       ! methods
       InitRandgen, &
       CloseRandgen, &
       GetGenSeed, &
       GetGenName, &
       RandgenInitialized, &
       RAN2
  !--------------------------------------------------------------------------!
CONTAINS

  SUBROUTINE InitRandgen(this,seed,un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP)   :: this
    INTEGER             :: seed
    INTEGER             :: J
    CHARACTER(LEN=*)    :: un
    !------------------------------------------------------------------------!
    INTENT(IN)          :: seed,un
    INTENT(INOUT)       :: this
    !------------------------------------------------------------------------!
    CALL InitCommon(this%gentype,seed,un)

    this%IDUM = seed
    this%IFF  = 0
    this%IR   = 0
    this%IY   = 0

    this%RM=1.0/M
    if  (this%IDUM<0 .or. this%IFF==0)  then
       this%IFF  = 1
       this%IDUM = modulo(IC-this%IDUM,M)
       do J=1,97
          this%IDUM = modulo(IA*this%IDUM+IC,M)
          this%IR(J)= this%IDUM
       end do
       this%IDUM = modulo(IA*this%IDUM+IC,M)
       this%IY   = this%IDUM
    endif

  END SUBROUTINE InitRandgen


  SUBROUTINE CloseRandgen(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP), INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    CALL CloseCommon(this%gentype)
  END SUBROUTINE CloseRandgen

  SUBROUTINE RAN2(this,randno)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP), INTENT(INOUT) :: this
    INTEGER                          :: J
    REAL(kind=r2),INTENT(OUT)        :: randno
    !------------------------------------------------------------------------!
      
    J = 1+(97*this%IY)/M
    if (J > 97 .or. J < 1) then
        print *, "subroutine RAN2(): failed."
        stop
    end if
    this%IY   = this%IR( J)
    this%rndx = this%IY*this%RM
    this%IDUM = modulo(IA*this%IDUM+IC,M)
    this%IR(J)= this%IDUM
    
    randno = this%rndx
  END SUBROUTINE RAN2

  PURE FUNCTION GetGenSeed(this) RESULT(ut)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP), INTENT(IN) :: this
    INTEGER :: ut
    !------------------------------------------------------------------------!
    ut = GetType_common(this%gentype)
  END FUNCTION GetGenSeed


  PURE FUNCTION GetGenName(this) RESULT(un)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: un
    !------------------------------------------------------------------------!
    un = GetName_common(this%gentype)
  END FUNCTION GetGenName

  PURE FUNCTION RandgenInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Randgen_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%gentype)
  END FUNCTION RandgenInitialized

End Module randgen_type
