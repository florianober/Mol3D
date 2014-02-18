! ---
! definition of density distribution
! --- 
! **** THIS IS THE PART WHICH IS OF INTEREST FOR MOST USERS OF THE MC3D ****
! ************ EDIT THIS PART YOURSELF AND RE-COMPILE THE CODE *************
! ---
module model_mod
    use datatype
    use var_globalnew
    implicit none
    public :: dustden!, get_n_den_par 
contains

  ! ################################################################################################
  ! return: number of user-defined parameters required to describe the density distribution
  !         given in 'subroutine dustden' below
  ! ---
  !subroutine get_n_den_par
  !  use var_global
    ! ---

    ! **********************************************************************************************
    ! BEGIN: Number of parameters
    ! ---
    !n_den_par = 0
    ! ---
    ! END  : Number of parameters
    ! **********************************************************************************************

  !end subroutine get_n_den_par


  ! ################################################################################################
  ! return: particle number density [ 1 / m^3 ]
  ! ---
  ! RULES: HOW TO DEFINE YOUR DENSITY DISRTRIBUTION
  !
  ! a) which are the variables I can make use of?
  !
  !    [1] r_in
  !        inner radius of your dust configuration (e.g. inner disk / shell radius)
  !        unit: your chosen reference unit (e.g., AU)
  !
  !    [2] r_ou
  !        outer radius of your dust configuration (e.g. outer disk / shell radius)
  !        unit: your chosen reference unit (e.g., AU)
  !
  !    [3] density
  !        this variable should contain the dust number density
  !        (in arbitrary units, because the dust mass will be normalized later)
  !    
  !    [4] den_par(1:n_den_par)
  !        - additional variables you may use to define parameters of your density distribution
  !          (type 'r2': equivalent to double precision; see datatype.f90)
  !
  !    [5] poscax(1:3)
  !        - contains cartesian coordinate of the point at which the density distribution
  !          shall be defined:
  !            caco(1) = x coordinate
  !            caco(2) = y coordinate
  !            caco(3) = z coordinate
  !        - unit: your chosen reference unit (e.g., AU)
  !
  !    [6] your own variables:
  !        should start with "P_" (e.g., P_start, P_radius, P_test, ...)
  !        since this prefix is reserverd for variables to be used by the public user
  !
  !    General:
  !    --------
  !    VARIABLES [1] - [4] : defined in the terminal dialog with the user
  !    VARIABLE  [5]       : defined during the runtime of the code
  !    VARIABLE  [6]       : either fixed value (in the code) or defined during runtime
  ! 
  !    ---
  !    less important (but available) variables:
  !
  !    [7] r_star
  !        - stellar radius [m]
  !
  !    [8] ref_unit
  !        - conversion factor of your reference unit to [m]
  !
  !
  ! b) example: see below
  ! ---
    FUNCTION dustden(model,caco) RESULT(density)
        USE model_type, ONLY : Model_TYP
        USE parser_mod
        
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                      :: model
        REAL(kind=r2), DIMENSION(1:3), INTENT(IN)        :: caco

        REAL(kind=r2)                                    :: density
        REAL(kind=r2)                                    :: P_xy, P_z, P_h
        REAL(kind=r2)                                    :: R_ou2, R_in2
        !--------------------------------------------------------------------------!
        
        ! ---

        ! ******************************************************************************************
        ! DEFINE YOUR DENSITY DISTRIBUTION HERE
        ! ---
        ! BEGIN: DEFINITION DENSITY DISTRIBUTION
        P_xy = sqrt(caco(1)**2+caco(2)**2)
        P_z  = abs(caco(3))
        R_ou2 = 20.0_r2
        R_in2 = 10.0_r2
!~         CALL parse('R_in2',R_in2,'input/input.dat')
!~         CALL parse('R_ou2',R_ou2,'input/input.dat')
        
!~         IF ( (P_xy >= 20.0_r2) .and. (P_xy < model%r_ou) ) THEN
!~            P_h     = 20.0_r2 * (P_xy/100.0_r2)**1.125_r2
!~            density = (P_xy/100.0_r2)**(-2.625_r2)  *  exp( -0.5_r2 * (P_z/P_h)**2 )
!~            !density = 100.0_r2
!~         ELSE
!~            density = 0.0_r2
!~         END IF
!~         print *, caco(2)
!~         IF  (abs(caco(2)) .lt. 5) print *, 'here'
!~         IF ( (P_z .lt. 5) .and. (abs(caco(1)-95.0) .lt. 5) .and. (abs(caco(2)) .lt. 5) ) THEN
!~             density = 100.0_r2
!~             print *, 'here'
!~         ELSE
!~             density = 0.0_r2
!~         END IF
!~         print *, (P_xy > R_ou2 .or. P_xy < R_in2), (P_xy >= model%r_in) .and. (P_xy < model%r_ou)
        IF ( (P_xy >= model%r_in) .and. (P_xy < model%r_ou) .and. (P_xy > R_ou2 .or. P_xy < R_in2) ) THEN
            P_h     = 20.0_r2 * (P_xy/100.0_r2)**1.125_r2
            density = (P_xy/100.0_r2)**(-2.625_r2)  *  exp( -0.5_r2 * (P_z/P_h)**2 )
        ELSE
            density = 0.0_r2
        END IF

        ! END: DEFINITION DENSITY DISTRIBUTION
        ! ******************************************************************************************
        
        ! DO NEVER REMOVE THIS LINE
        
        density = density *model%mass !old relict
        !grid%ddust(:) = model%mass * density * sidi(:)
    END FUNCTION dustden
  
end module model_mod

