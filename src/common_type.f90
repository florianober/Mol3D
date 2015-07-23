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
! basic data and methods common to all objects 
! inspired by fosite by T. Illenseer
!------------------------------------------------------------------------------!
MODULE common_type

  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  PRIVATE
  ! common data structure
  !----------------------------------------------------------------------------!
  TYPE Common_TYP
     PRIVATE
     INTEGER           :: type
     CHARACTER(LEN=32) :: name
     INTEGER           :: error               ! error code                     !
     LOGICAL           :: init = .FALSE.      ! init status                    !
  END TYPE Common_TYP
  SAVE
  ! these variables should be the same for all objects
  !----------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Common_TYP, &
       ! methods
       InitCommon, &
       CloseCommon, &
       GetType, &
       GetName, &
       Initialized
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitCommon(this,t,n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP)  :: this
    INTEGER           :: t
    CHARACTER(LEN=*)  :: n
    !------------------------------------------------------------------------!
    INTENT(IN)        :: t,n
    INTENT(OUT)       :: this
    !------------------------------------------------------------------------!
    this%type = t
    this%name = n
    this%error = 0
    this%init  = .TRUE.
  END SUBROUTINE InitCommon

  SUBROUTINE CloseCommon(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)    :: this
    !------------------------------------------------------------------------!
    this%init = .FALSE.
  END SUBROUTINE CloseCommon


  PURE FUNCTION GetType(this) RESULT(t)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    INTEGER :: t
    !------------------------------------------------------------------------!
    t = this%type
  END FUNCTION GetType


  PURE FUNCTION GetName(this) RESULT(n)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: n
    !------------------------------------------------------------------------!
    n = this%name
  END FUNCTION GetName

  PURE FUNCTION Initialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Common_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = this%init
  END FUNCTION Initialized


END MODULE common_type
