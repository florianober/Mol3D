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
    public :: get_den!, get_n_den_par 
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
    FUNCTION get_den(model,caco) RESULT(density)
        USE model_type, ONLY : Model_TYP
        
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                      :: model
        REAL(kind=r2), DIMENSION(1:3), INTENT(IN)        :: caco

        REAL(kind=r2)                                    :: density
        REAL(kind=r2)                                    :: P_xy, P_z, P_h
        !--------------------------------------------------------------------------!
        
        ! ---

        ! ******************************************************************************************
        ! DEFINE YOUR DENSITY DISTRIBUTION HERE
        ! ---
        ! BEGIN: DEFINITION DENSITY DISTRIBUTION
        P_xy = sqrt(caco(1)**2+caco(2)**2)
        P_z  = abs(caco(3))

        IF ( (P_xy >= model%r_in) .and. (P_xy < model%r_ou) ) THEN
            P_h     = 20.0_r2 * (P_xy/100.0_r2)**1.125_r2
            density = (P_xy/100.0_r2)**(-2.625_r2)  *  exp( -0.5_r2 * (P_z/P_h)**2 )
        ELSE
            density = 0.0_r2
        END IF

        ! END: DEFINITION DENSITY DISTRIBUTION
        ! ******************************************************************************************
        
        ! DO NEVER REMOVE THIS LINE
        
        density = density *model%mass !old relict from mc3d, tbd: can be removed
        
    END FUNCTION get_den
  
end module model_mod

