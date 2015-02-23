MODULE error_mod

    USE datatype

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: criticals  

CONTAINS

    SUBROUTINE criticals(depth, msg)

        IMPLICIT NONE
    
        !----------------------------------------------------------------------!
        INTEGER            :: depth
        CHARACTER(len=*)   :: msg
        !----------------------------------------------------------------------!

        IF ( depth == 1) THEN

        ELSE IF (depth == 2) THEN
             print '(2A)', 'WARNING: ', msg       
        
        ELSE IF (depth == 3) THEN
            print '(2A)', 'CRITICAL ERROR: ', msg
            print '(A)', 'stoping Mol3D!'
            stop
        ELSE
            print *, 'wrong input in subroutine crititals (var: depth), ->     &
                      & stoping Mol3D!'
            stop
        END IF

    END SUBROUTINE criticals

END MODULE error_mod
