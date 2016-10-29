!------------------------------------------------------------------------------!
! This file is part of Mol3D.
!
!    Mol3D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Mol3D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Mol3D.  If not, see <http://www.gnu.org/licenses/>.
!
!    Diese Datei ist Teil von Mol3D.
!
!    Mol3D ist Freie Software: Sie können es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
!    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Mol3D wird in der Hoffnung, dass es nützlich sein wird, aber
!    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License für weitere Details.
!
!    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
!    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------!
! In this module the analytical (disk) model is defined
! -> please adjust for your needs and don't forget to recompile Mol3D !
! --- 

MODULE model_mod
    USE datatype
    USE var_global

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: get_density, get_velocity, get_temperature, get_J_ext

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

    ! -------------------------------------------------------------------------!
    ! external radiation field (required to calculate lvl populations)
    ! -------------------------------------------------------------------------!
    FUNCTION get_J_ext(freq) RESULT(J_ext)
        USE math_mod, ONLY : planckhz
        ! analytical temperature distribution
        !
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(:), INTENT(IN)  :: freq
        REAL(kind=r2), DIMENSION(1:size(freq))              :: J_ext
        !----------------------------------------------------------------------!

        !define your own external radiation field here
!~         J_ext(:) = 1e-14_r2*planckhz(10000.0, freq(:)) +                       &
!~                             planckhz(2.73, freq(:))
        J_ext(:) = planckhz(2.73, freq(:))
    END FUNCTION get_J_ext
END MODULE model_mod
