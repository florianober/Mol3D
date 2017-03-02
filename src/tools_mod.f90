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
! Module tools_mod;
!        contains several tools, which are often used everywhere
! ---
MODULE tools_mod

    USE datatype
    USE var_global

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE 
    !--------------------------------------------------------------------------!
    PUBLIC  :: Mol3D_Error
    !--------------------------------------------------------------------------!

    CONTAINS
        SUBROUTINE Mol3D_Error(Error_msg, Error_type)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN), OPTIONAL                         :: Error_type
        INTEGER                                               :: Error_code
        CHARACTER(len=*), INTENT(IN), OPTIONAL                :: Error_msg
        CHARACTER(len=256)                                    :: Error_msg_show
        !----------------------------------------------------------------------!
        
        IF (PRESENT(ERROR_type)) THEN
            Error_code = Error_type
        ELSE
            Error_code = 1
        END IF

        IF (PRESENT(Error_msg)) THEN
            Error_msg_show = Error_msg
        ELSE
            Error_msg_show = "no reason specified"
        END IF
        
        !----------------------------------------------------------------------!
        PRINT *, ""
        PRINT *, "ERROR: "  // TRIM(Error_msg_show)

        IF (Error_code == 1) THEN
            PRINT *, "Mol3D has stopped"
            PRINT *, "If you don't know why this has happend,"
            PRINT *, "please ask fober@astrophysik.uni-kiel.de"
            STOP
        END IF
        END SUBROUTINE Mol3D_Error

END MODULE tools_mod

