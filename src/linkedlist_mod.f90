MODULE linkedlist_mod

    USE datatype

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    
    TYPE Pixel_TYP
        !----------------------------------------------------------------------!
        
        INTEGER, DIMENSION(1:2)                          :: pixel
        
        REAL(kind=r2), DIMENSION(1:2)                    :: pos_xy
        REAL(kind=r2), DIMENSION(1:2)                    :: size_xy

    END TYPE Pixel_TYP

    TYPE l_list
        PRIVATE
        TYPE (l_list), POINTER        :: next => null()
        
        !----------------------------------------------------------------------!
        
        TYPE(Pixel_TYP)               :: p_data
        
        !----------------------------------------------------------------------!
    END TYPE l_list
    SAVE
    !--------------------------------------------------------------------------!
    PUBLIC :: &
        ! types
        l_list, &
        Pixel_TYP, &
        ! methods
        AddElement, &
        ShowAll, &
        RemoveAll, &
        RemoveLast, &
        GetSize, &
        GetLast
!~         GetLastElement
    CONTAINS

        SUBROUTINE AddElement(this,p_data)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(l_list), POINTER              :: new
        TYPE(l_list), POINTER              :: temp
        
        Type(Pixel_TYP)                    :: p_data
        
        !----------------------------------------------------------------------!
        INTENT(INOUT)                      :: this
        INTENT(IN)                         :: p_data
        !----------------------------------------------------------------------!
        
        ALLOCATE(new)
        
        new%p_data = p_data

        
        
        IF (ASSOCIATED(this) .eqv. .FALSE.) THEN
            this => new
        ELSE
            
            temp      => this
            this      => new
            this%next => temp
        
        END IF
        
        
        END SUBROUTINE AddElement
        
        SUBROUTINE ShowAll(this)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(l_list), POINTER              :: temp
        !----------------------------------------------------------------------!
        INTENT(INOUT)                      :: this
        !----------------------------------------------------------------------!
        
        temp => this
        
        DO WHILE (associated(temp))
            print *,  temp%p_data%pixel(1)
            temp => temp%next
        END DO
        
        END SUBROUTINE ShowAll
        
        
        SUBROUTINE RemoveAll(this)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(l_list), POINTER              :: temp

        !----------------------------------------------------------------------!
        INTENT(INOUT)                      :: this
        !----------------------------------------------------------------------!

        DO WHILE (associated(this))
            temp => this
            this => this%next
            
            DEALLOCATE(temp)
         
        END DO
        
        END SUBROUTINE RemoveAll
        
        
        SUBROUTINE RemoveLast(this)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(l_list), POINTER              :: temp
        
        !----------------------------------------------------------------------!
        INTENT(INOUT)                      :: this
        !----------------------------------------------------------------------!

        
        temp => this
        
        IF ( ASSOCIATED(temp) ) THEN
            
            this   => temp%next
            DEALLOCATE(temp)
            
        END IF
        
        END SUBROUTINE RemoveLast
        
        SUBROUTINE GetLast(this,last)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(Pixel_TYP), POINTER           :: last
        
        !----------------------------------------------------------------------!
        INTENT(IN)                         :: this
        INTENT(OUT)                        :: last
        !----------------------------------------------------------------------!
        
        IF ( ASSOCIATED(this) ) THEN
            
            last => this%p_data
        ELSE 
            last => null()
        
        END IF
        
        END SUBROUTINE GetLast
        
        FUNCTION GetSize(this) RESULT(size_int)
        !----------------------------------------------------------------------!
        TYPE(l_list), POINTER              :: this
        TYPE(l_list), POINTER              :: temp
        INTEGER                            :: size_int
        
        !----------------------------------------------------------------------!
        INTENT(IN)                         :: this
        !----------------------------------------------------------------------!
        size_int = 0
        
        temp => this
        
        DO WHILE  ( ASSOCIATED(temp) )
            
            size_int = size_int +1
            temp => temp%next
            
        END DO
        
        END FUNCTION GetSize
        

END MODULE linkedlist_mod

!~ PROGRAM tests
!~     
!~     USE datatype
!~     USE linked_list
!~     
!~     IMPLICIT NONE 
!~     !-----------------------------------------------------------------------!
!~     TYPE(Pixel_TYP) ,POINTER                       :: pixel
!~     TYPE(l_list),POINTER                           :: Liste => null()
!~     
!~     INTEGER                                        :: i,j
!~     !-----------------------------------------------------------------------!
!~ 
!~     
!~     ALLOCATE(pixel)
!~     DO i = 1, 1
!~         DO j = 1, 2
!~             pixel%pixel(1) = i
!~             pixel%pixel(2) = j
!~             pixel%pos_xy   = 0.0_r2
!~             pixel%size_xy   = 1.1_r2
!~             CALL AddElement(Liste,pixel)
!~         END DO
!~     END DO
!~     
!~  
!~     
!~     
!~ 
!~ END PROGRAM tests
