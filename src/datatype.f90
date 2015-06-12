!
! definion of data types
! ---
! goal:    equip Mol3D with the flexibility to change the accuracy of floating
!          point variables,
!          in particular in the case of "double precision" variables (here: r2)
!
! how2use: e.g. let a be a variable of type double precision:
!
!                 use datatype
!                 real(kind=r2) :: a
!                 a = 2.0_r2
!

MODULE datatype
    IMPLICIT NONE

    ! real
    INTEGER, PARAMETER, PUBLIC :: r1=selected_real_kind(1)

    ! double precision
    INTEGER, PARAMETER, PUBLIC :: r2=selected_real_kind(p=12)
    
    ! quad precision
    ! INTEGER, PARAMETER, PUBLIC :: r3=selected_real_kind(p=21)

END MODULE datatype
