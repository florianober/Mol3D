! ===
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
! ===
module datatype
 implicit none

 ! real
 integer, parameter, public :: r1=selected_real_kind(1)

 ! double precision
 integer, parameter, public :: r2=selected_real_kind(p=12)
 
 ! quad precision
 ! integer, parameter, public :: r3=selected_real_kind(p=21)
 
 
end module datatype
