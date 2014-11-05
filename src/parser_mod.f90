!
! This module has been developed to handle input files in a more dynamic way
!
!
MODULE parser_mod

    USE datatype
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC :: parse
    !--------------------------------------------------------------------------!
    INTERFACE parse
        MODULE PROCEDURE parse_int_i1, parse_int_i2, parse_real, parse_str, parse_log
    END INTERFACE
    !--------------------------------------------------------------------------!
    
    CONTAINS
    
    SUBROUTINE parse_int_i1(keyword,var,data_file)

        IMPLICIT NONE
        ! 
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN) :: keyword 
        CHARACTER(len=*),INTENT(IN) :: data_file
        INTEGER,INTENT(INOUT)       :: var
        
        CHARACTER(len=256)          :: val_st
        INTEGER                     :: iost
        !------------------------------------------------------------------------!
        CALL get_value(keyword,data_file,val_st,iost)
        IF ( iost == 0 ) THEN
            READ(val_st,fmt=*,IOSTAT=iost) var
        END IF
        IF ( iost /= 0 ) THEN
            PRINT '(2A)','ERROR in input file ('//TRIM(data_file)//') with keyword: ', keyword
            STOP
        END IF
    END SUBROUTINE parse_int_i1
    
    SUBROUTINE parse_int_i2(keyword,var,data_file)

        IMPLICIT NONE
        ! 
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN) :: keyword 
        CHARACTER(len=*),INTENT(IN) :: data_file
        INTEGER(8),INTENT(INOUT)    :: var
        
        CHARACTER(len=256)          :: val_st
        INTEGER                     :: iost
        !------------------------------------------------------------------------!
        CALL get_value(keyword,data_file,val_st,iost)
        IF ( iost == 0 ) THEN
            READ(val_st,fmt=*,IOSTAT=iost) var
        END IF
        IF ( iost /= 0 ) THEN
            PRINT '(2A)','ERROR in input file ('//TRIM(data_file)//') with keyword: ', keyword
            STOP
        END IF
    END SUBROUTINE parse_int_i2
    
    
    SUBROUTINE parse_real(keyword,var,data_file)

        IMPLICIT NONE
        ! 
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN) :: keyword 
        CHARACTER(len=*),INTENT(IN) :: data_file
        REAL(kind=r2),INTENT(INOUT) :: var
        
        CHARACTER(len=256)          :: val_st
        INTEGER                     :: iost
        !------------------------------------------------------------------------!
        CALL get_value(keyword,data_file,val_st,iost)
        IF ( iost == 0 ) THEN
            READ(val_st,fmt=*,IOSTAT=iost) var
        END IF
        IF ( iost /= 0 ) THEN
            PRINT '(2A)','ERROR in input file ('//TRIM(data_file)//') with keyword: ', keyword
            STOP
        END IF


    END SUBROUTINE parse_real
    
    
    SUBROUTINE parse_str(keyword,var,data_file)

        IMPLICIT NONE
        ! 
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN)    :: keyword 
        CHARACTER(len=*),INTENT(IN)    :: data_file
        CHARACTER(len=*),INTENT(INOUT) :: var
        
        CHARACTER(len=256)          :: val_st
        INTEGER                     :: iost
        !------------------------------------------------------------------------!
        CALL get_value(keyword,data_file,val_st,iost)
        IF ( iost == 0 ) THEN
            READ(val_st,fmt='(A)',IOSTAT=iost) var
        END IF
        IF ( iost /= 0 ) THEN
            PRINT '(2A)','ERROR in input file ('//TRIM(data_file)//') with keyword: ', keyword
            STOP
        END IF


    END SUBROUTINE parse_str
    
    SUBROUTINE parse_log(keyword,var,data_file)

        IMPLICIT NONE
        ! 
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN)    :: keyword 
        CHARACTER(len=*),INTENT(IN)    :: data_file
        LOGICAL,INTENT(INOUT)          :: var
        
        CHARACTER(len=256)          :: val_st
        INTEGER                     :: iost
        !------------------------------------------------------------------------!
        CALL get_value(keyword,data_file,val_st,iost)
        IF ( iost == 0 ) THEN
            READ(val_st,fmt='(L1)',IOSTAT=iost) var
        END IF
        IF ( iost /= 0 ) THEN
            PRINT '(2A)','ERROR in input file ('//TRIM(data_file)//') with keyword: ', keyword
            STOP
        END IF


    END SUBROUTINE parse_log
    
    
    !This is the main work
    !
    SUBROUTINE get_value(keyword,data_file,val_st,iost)

        IMPLICIT NONE
        !This routine seaches each line for the keyword and returns the value as a string
        !
        !------------------------------------------------------------------------!
        CHARACTER(len=*),INTENT(IN)      :: keyword 
        CHARACTER(len=*),INTENT(IN)      :: data_file
        CHARACTER(len=256),INTENT(OUT)   :: val_st
        INTEGER,INTENT(OUT)              :: iost
        CHARACTER(len=256)               :: line
        
        INTEGER                          :: i,j
        INTEGER                          :: br_l, br_r
        INTEGER                          :: io
        INTEGER                          :: key_found
        
        
        !------------------------------------------------------------------------!
        iost = 0
        
        OPEN(unit=1, file=data_file, &
             action="read", status="unknown", form="formatted")
        key_found = 0
        line = ''
        DO WHILE (key_found == 0)
            ! read the next line of the input file
            !
            READ(unit=1,fmt='(A)',iostat=io) line
            
            IF (io < 0) THEN
                !found end of file
                !PRINT *, 'ERROR: Keyword not found in input file'
                iost = -1
                EXIT 
            END IF
            j = len_trim(line)
            
            ! search for the '=' symbol in this line
            !
            DO i=1,j
                IF (line(i:i) == '=' ) THEN
                    EXIT
                END IF
            END DO
            
            ! test if the keyword belong to this entry
            !
            IF (line(1:i-1) == keyword) THEN
                key_found = 1
                !now, find the brackets {}
                DO br_l = i+1, j-1
                    IF (line(br_l:br_l) == '{' ) THEN
                        EXIT
                    END IF
                END DO
                
                DO br_r = br_l+1, j
                    IF (line(br_r:br_r) == '}' ) THEN
                        EXIT
                    END IF
                END DO
                
                !check if brackets have been found.
                IF (br_l == j .or. br_r == j+1 ) THEN
                    iost = 1
                    !PRINT *, 'ERROR: no value found (brackets not correct?)'
                    EXIT
                END IF

                !now set the file value (inside the {} brackets) to the given var 
                !
                val_st = line(br_l+1:br_r-1)
            END IF
            
        END DO
        CLOSE(unit=1)
    
    END SUBROUTINE get_value
    
END MODULE parser_mod
