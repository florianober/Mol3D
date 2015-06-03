! ---
! In this module the analytical (disk) model is defined
! -> please adjust for your needs and don't forget to recompile the Mol3D !
! --- 

MODULE model_mod
    USE datatype
    USE var_global

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: get_density, get_velocity, get_temperature

CONTAINS
    ! -------------------------------------------------------------------------!
    ! density
    ! -------------------------------------------------------------------------!
    FUNCTION get_density(model, caco, model_name_in) RESULT(density)
        ! generic routine to set the density distribution
        USE model_type, ONLY : Model_TYP
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                      :: model
        REAL(kind=r2), DIMENSION(1:3), INTENT(IN)        :: caco

        REAL(kind=r2)                                    :: density

        CHARACTER(len=*), OPTIONAL                       :: model_name_in
        CHARACTER(len=25)                                :: model_name
        !----------------------------------------------------------------------!
        IF (.not. present(model_name_in)) THEN
            model_name = 'FLARED_DISK'
        ELSE
            model_name = model_name_in
        END IF
        
        SELECT CASE(TRIM(model_name))

        CASE('FLARED_DISK')
            density = get_density_disk(caco, r_in=model%r_in, r_ou=model%r_ou, &
                                   h_0=10.0_r2, alpha=2.625_r2, beta=1.125_r2)
        CASE DEFAULT
            PRINT *, "ERROR, density model '"//TRIM(model_name)//              &
                     "' not defined (see model_mod.f90)"
            STOP
        END SELECT
    END FUNCTION get_density

    PURE FUNCTION get_density_disk(caco, r_in, r_ou, h_0, alpha, beta)         &
                                   RESULT(density)
        !
        ! A simple Shakura & Sunyaev (1973) accretion disk
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(1:3), INTENT(IN)        :: caco
        REAL(kind=r2), INTENT(IN)                        :: r_in, r_ou
        REAL(kind=r2), INTENT(IN)                        :: alpha, beta, h_0
        
        REAL(kind=r2)                                    :: density
        REAL(kind=r2)                                    :: P_xy, P_z, P_h
        !----------------------------------------------------------------------!
        
        P_xy = sqrt(caco(1)**2 + caco(2)**2)
        P_z  = abs(caco(3))
        P_h     = h_0 * (P_xy/100.0_r2)**beta
        
        IF ( (P_xy >= r_in) .and. (P_xy < r_ou) ) THEN
            density = (P_xy/100.0_r2)**(-alpha) * exp(-0.5_r2 * (P_z/P_h)**2)
        ELSE
            density = 0.0_r2
        END IF
    END FUNCTION get_density_disk

    ! -------------------------------------------------------------------------!
    ! velocity
    ! -------------------------------------------------------------------------!
    FUNCTION get_velocity(model, caco, model_name_in) RESULT(velo)
        ! generic routine to set the velocity
        USE model_type, ONLY : Model_TYP
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN)                      :: model
        REAL(kind=r2), DIMENSION(1:3), INTENT(IN)        :: caco

        REAL(kind=r1), DIMENSION(1:3)                    :: velo
        CHARACTER(len=*), OPTIONAL                       :: model_name_in
        CHARACTER(len=25)                                :: model_name
        !----------------------------------------------------------------------!
        IF (.not. present(model_name_in)) THEN
            model_name = 'KEPLERIAN'
        ELSE
            model_name = model_name_in
        END IF
        
        SELECT CASE(TRIM(model_name))

        CASE('KEPLERIAN')
            velo = get_velocity_keplerian(caco, model%kep_const)
        CASE DEFAULT
            PRINT *, "ERROR, velocity model '"//TRIM(model_name)//             &
                     "' not defined (see model_mod.f90)"
            STOP
        END SELECT
    END FUNCTION get_velocity

    PURE FUNCTION get_velocity_keplerian(caco, kep_const) RESULT(velo)
    
        ! simple Keplerian rotation
        !
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)              :: caco
        REAL(kind=r1), DIMENSION(1:3)                         :: velo
        REAL(kind=r2)                                         :: konst, expp, r
        REAL(kind=r1),INTENT(IN)                              :: kep_const
        !----------------------------------------------------------------------!
        ! we assume pure keplerian rotation, therefore v(r) ~ r**-0.5

        konst = kep_const  !
        
        expp = -1.5  ! is a result of keplerian rotation and
                     ! coordinate system conversion
        r = (caco(1)**2+caco(2)**2 + EPSILON(1.0_r2))**(expp*0.5) * konst
        
        velo(1) = (-1.0) * caco(2) * r
        velo(2) =          caco(1) * r
        velo(3) = 0.0
        
    END FUNCTION get_velocity_keplerian
    ! -------------------------------------------------------------------------!
    ! temperature
    ! -------------------------------------------------------------------------!
    PURE FUNCTION get_temperature(caco) RESULT(temps)
        ! analytical temperature distribution
        !
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(1:3),INTENT(IN)         :: caco
        REAL(kind=r2)                                    :: temps
        REAL(kind=r2)                                    :: konst, P_xy, P_h
        
        !----------------------------------------------------------------------!

        !define your own temperature distribution here!
        konst = 100.0_r2
        
        P_xy = sqrt(caco(1)**2 + caco(2)**2 + EPSILON(1.0))
        P_h  = 5.0_r2 * (P_xy/100.0_r2)**1.125_r2

        !temps   = 500.0_r2*P_xy**(-0.5) * (2.0-exp(-0.5*(abs(caco(3))/P_h)**2))
        temps   = 200.0_r2*P_xy**(-0.5)

    END FUNCTION get_temperature
END MODULE model_mod

