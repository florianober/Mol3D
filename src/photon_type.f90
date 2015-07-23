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
! def PHOTON_TYP: This type includes all informations of the active photon
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE photon_type
  
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE


    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE PHOTON_TYP
        TYPE(Common_TYP) :: mtype                      ! -----------------     !
        !----------------------------------------------------------------------!
        LOGICAL                                         :: inside
        LOGICAL                                         :: kill
        LOGICAL                                         :: fixed_lam
        
        CHARACTER(1)                                    :: last_interaction_type
        
        INTEGER                                         :: n_interact
        INTEGER                                         :: nr_cell
        INTEGER                                         :: nr_cell_new
        INTEGER                                         :: nr_lam
        INTEGER                                         :: lot_th
        
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_li
        REAL(kind=r2), DIMENSION(1:3)                   :: pos_xyz_new
        REAL(kind=r2), DIMENSION(1:3)                   :: dir_xyz
        REAL(kind=r2), DIMENSION(1:3,1:3)               :: D_2global
        REAL(kind=r2), DIMENSION(1:4)                   :: stokes
        REAL(kind=r2)                                   :: SINPHI, COSPHI
        REAL(kind=r2)                                   :: SIN2PH, COS2PH
        REAL(kind=r2)                                   :: SINTHE, COSTHE
        
        REAL(kind=r2)                                   :: energy
        REAL(kind=r2)                                   :: energy_scale


    END TYPE PHOTON_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        PHOTON_TYP, &
        ! methods
        InitPhoton, &
        ClosePhoton, &
        GetPhotonType, &
        GetPhotonName, &
        PhotonInitialized, &
        vecmat, &
        update_angle
    INTERFACE vecmat
        MODULE PROCEDURE vecmat_angle, vecmat_dir
    END INTERFACE
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitPhoton(this, ut, un, i_lam_in)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP)                            :: this
        
        INTEGER                                     :: ut
        INTEGER, INTENT(IN), OPTIONAL               :: i_lam_in
        
        CHARACTER(LEN=*)        :: un
        !----------------------------------------------------------------------!
        INTENT(IN)              :: ut, un
        INTENT(INOUT)           :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        
        this%inside      = .True.
        this%kill        = .False.

        this%n_interact  = 0
        this%nr_cell     = 0
        this%nr_cell_new = 0
        IF (PRESENT(i_lam_in)) THEN
            this%nr_lam = i_lam_in
            this%fixed_lam = .True.
        ELSE
            this%nr_lam  = 0
            this%fixed_lam = .False.
        END IF
        this%lot_th      = 0
        
        this%pos_xyz     = 0.0_r2
        this%pos_xyz_new = 0.0_r2
        this%pos_xyz_li  = 0.0_r2 ! position of last interaction
        this%dir_xyz     = 0.0_r2

        this%D_2global    = 0.0_r2
        this%stokes      = 0.0_r2
        
        this%SINPHI      = 0.0_r2 
        this%COSPHI      = 0.0_r2 
        this%SINTHE      = 0.0_r2 
        this%COSTHE      = 0.0_r2
        this%SIN2PH      = 0.0_r2
        this%COS2PH      = 0.0_r2
        
        this%energy    = 0.0_r2
        this%last_interaction_type = 'N'

    END SUBROUTINE InitPhoton


    SUBROUTINE ClosePhoton(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        
    END SUBROUTINE ClosePhoton


    PURE FUNCTION GetPhotonType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetPhotonType


    PURE FUNCTION GetPhotonName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetPhotonName

    PURE FUNCTION PhotonInitialized(this) RESULT(i)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN) :: this
        LOGICAL :: i
        !----------------------------------------------------------------------!
        i = Initialized_common(this%mtype)
    END FUNCTION PhotonInitialized

    SUBROUTINE vecmat_angle(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT)                     :: this

        REAL(kind=r2), DIMENSION(1:3,1:3)                   :: D_help
        REAL(kind=r2), DIMENSION(1:3)                       :: dir_rlp
        !----------------------------------------------------------------------!
        
        ! old MC3D version -> not used anymore
        ! D_test = this%D_2global
        ! R = -this%SINTHE * this%SINPHI
        ! L =  this%SINTHE * this%COSPHI
        ! P =  this%COSTHE

        ! this%dir_xyz(2) = -(this%D_2global(1,1) * R  +  this%D_2global(1,2) * L  +  this%D_2global(1,3) * P)
        ! this%dir_xyz(1) =  this%D_2global(2,1) * R  +  this%D_2global(2,2) * L  +  this%D_2global(2,3) * P
        ! this%dir_xyz(3) =  this%D_2global(3,1) * R  +  this%D_2global(3,2) * L  +  this%D_2global(3,3) * P
        
        ! D_help = set_matrix(this)
        ! D_test = matmul(D_test, D_help)

        ! D_help(1,1) =   this%COSPHI
        ! D_help(2,1) =   this%SINPHI
        ! D_help(3,1) =   0.0_r2
        ! D_help(1,2) = - this%SINPHI*this%COSTHE
        ! D_help(2,2) =   this%COSPHI*this%COSTHE
        ! D_help(3,2) = - this%SINTHE
        ! D_help(1,3) = - this%SINPHI*this%SINTHE
        ! D_help(2,3) =   this%COSPHI*this%SINTHE
        ! D_help(3,3) =   this%COSTHE

        dir_rlp(1) = - this%SINTHE * this%SINPHI
        dir_rlp(2) =   this%SINTHE * this%COSPHI
        dir_rlp(3) =   this%COSTHE
        
        this%dir_xyz = matmul(this%D_2global(:,:), dir_rlp)
        D_help = set_matrix(this)

        this%D_2global = matmul(this%D_2global, D_help)
        
    END SUBROUTINE vecmat_angle

    SUBROUTINE vecmat_dir(this, d_ang)
        USE math_mod, ONLY : atan3, solve_eq

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT)                     :: this
        
        real(kind=r2), dimension(1:3,1:3)                   :: D_help
        REAL(kind=r2), dimension(1:3)                       :: dir_help
        REAL(kind=r2), dimension(1:3)                       :: dir_rlp
        REAL(kind=r2)                                       :: theta, phi
        REAL(kind=r2)                                       :: d_ang
        !----------------------------------------------------------------------!

        ! the direction in the global coordinate system is given,
        ! now we need to calculate the theta and phi angle in the photon frame

        ! get the direction in the photon frame
        dir_rlp = matmul(transpose(this%D_2global), this%dir_xyz)

        phi = atan3(- dir_rlp(1), dir_rlp(2))
        this%SINPHI = sin(phi)
        this%COSPHI = cos(phi)
        
        theta = acos(dir_rlp(3))
        this%COSTHE = dir_rlp(3)
        this%SINTHE = sqrt(1-this%COSTHE**2)
        ! now we can calculate the remaining values

        this%SIN2PH = 2.0_r2 * this%SINPHI * this%COSPHI
        this%COS2PH = 1.0_r2 - 2.0_r2 * this%SINPHI**2
        this%lot_th = int(theta / d_ang) + 1

        ! fill the rotation matrix
        ! new definition
        D_help = set_matrix(this)

        this%D_2global = matmul(this%D_2global, D_help)

    END SUBROUTINE vecmat_dir

    PURE FUNCTION set_matrix(this) RESULT(D_rot)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(IN)                        :: this

        REAL(kind=r2), DIMENSION(1:3,1:3)                   :: D_rot
        !----------------------------------------------------------------------!
        ! old definition
        D_rot(1,1) =   this%COSPHI
        D_rot(2,1) =   this%SINPHI
        D_rot(3,1) =   0.0_r2
        D_rot(1,2) = - this%SINPHI*this%COSTHE
        D_rot(2,2) =   this%COSPHI*this%COSTHE
        D_rot(3,2) = - this%SINTHE
        D_rot(1,3) = - this%SINPHI*this%SINTHE
        D_rot(2,3) =   this%COSPHI*this%SINTHE
        D_rot(3,3) =   this%COSTHE
    END FUNCTION set_matrix

    PURE SUBROUTINE update_angle(this)
        !calculate SIN2PH/COS2PH for the stokes rotation
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(PHOTON_TYP), INTENT(INOUT)                        :: this
        !----------------------------------------------------------------------!
        this%SIN2PH = 2.0_r2 * this%SINPHI * this%COSPHI
        this%COS2PH = 1.0_r2 - 2.0_r2 * this%SINPHI**2
    END SUBROUTINE update_angle

End Module photon_type
