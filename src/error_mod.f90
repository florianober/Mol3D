MODULE error_mod

    USE datatype

    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: criticals  
  

CONTAINS

    SUBROUTINE criticals(deep, msg)

        IMPLICIT NONE
    
        !--------------------------------------------------------------------------!
        INTEGER            :: deep
        CHARACTER(len=*)  :: msg
        !--------------------------------------------------------------------------!

        IF ( deep == 1) THEN
        
        ELSE IF (deep == 2) THEN
             print '(2A)', 'WARNING: ', msg       
        
        ELSE IF (deep == 3) THEN
            print '(2A)', 'CRITICAL ERROR: ', msg
            print '(A)', 'stoping mol3d!'
            stop
        ELSE
            print *, 'wrong input in subroutine crititals (var: deep), -> stoping mol3d!'
            stop
        END IF

    END SUBROUTINE criticals


END MODULE error_mod
