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
! def Model_TYP
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE model_type
  
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    USE math_mod, ONLY : grad2rad
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Model_TYP
        TYPE(Common_TYP) :: modltype                    ! -----------------    !
        !----------------------------------------------------------------------!
        REAL(kind=r2)        :: ref_unit
        REAL(kind=r2)        :: ref_unitn
        REAL(kind=r2)        :: r_in
        REAL(kind=r2)        :: r_ou
        REAL(kind=r2)        :: mass
        REAL(kind=r1)        :: r_star
        REAL(kind=r1)        :: M_star
        REAL(kind=r1)        :: t_star
        REAL(kind=r1)        :: l_star
        REAL(kind=r1)        :: distance
        REAL(kind=r1)        :: kep_const
        
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: th_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: ph_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: al_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: zoom_map
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: px_model_length_x
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE         :: px_model_length_y
        
        REAL(kind=r2),DIMENSION(:,:),ALLOCATABLE       :: dir_xyz
        REAL(kind=r2),DIMENSION(:,:),ALLOCATABLE       :: e_x
        REAL(kind=r2),DIMENSION(:,:),ALLOCATABLE       :: e_y
        REAL(kind=r2),DIMENSION(:,:,:),ALLOCATABLE     :: D_2obs 
        
        
        INTEGER               :: n_map
        INTEGER               :: n_bin_map
        INTEGER(8)            :: no_photon
        
        
    END TYPE Model_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Model_TYP, &
        ! methods
        InitModel, &
        CloseModel, &
        GetModelType, &
        GetModelName, &
        ModelInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitModel(this, ut , un, ref_u, r_in, r_ou, &
                    mass_dust, t_eff, r_star,M_star, n_map, distance, &
                    no_photon,th_map,ph_map,zoom_map,al_map,n_bin_map)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP)       :: this
        
        INTEGER               :: ut
        INTEGER               :: n_map
        INTEGER               :: n_bin_map
        INTEGER               :: i_map
        INTEGER(8)            :: no_photon
        
        CHARACTER(LEN=*)      :: un
        
        REAL(kind=r2)         :: ref_u
        REAL(kind=r2)         :: distance
        REAL(kind=r2)         :: r_in
        REAL(kind=r2)         :: r_ou
        REAL(kind=r2)         :: mass_dust
        REAL(kind=r2)         :: t_eff
        REAL(kind=r2)         :: r_star
        REAL(kind=r2)         :: M_star
        
        REAL(kind=r2),DIMENSION(1:n_map)         :: th_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: ph_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: al_map
        REAL(kind=r2),DIMENSION(1:n_map)         :: zoom_map

        
        !----------------------------------------------------------------------!
        INTENT(IN)          ::  ut,un,ref_u, r_in, r_ou, &
                                mass_dust, t_eff, r_star, M_star, n_map, distance, &
                                th_map,ph_map,zoom_map,al_map,n_bin_map,no_photon
        INTENT(INOUT)       :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%modltype,ut,un)
        
        this%ref_unit = ref_u
        this%ref_unitn = 1.0_r2/ref_u
        
        this%t_star = t_eff
        this%r_star = r_star*R_sun
        this%M_star = M_star*M_sun
        
        this%kep_const = (con_gamma * this%M_star/con_AU)**0.5  ! in units of m/s
        this%l_star = 4.0_r2 * PI * this%r_star**2 * con_sigma * this%t_star**4

        ALLOCATE(   this%th_map(1:n_map), &
                    this%ph_map(1:n_map), &
                    this%zoom_map(1:n_map), &
                    this%al_map(1:n_map), &
                    this%px_model_length_x(1:n_map), &
                    this%px_model_length_y(1:n_map), &
                    this%dir_xyz(1:3,1:n_map), &
                    this%e_x(1:3,1:n_map), &
                    this%e_y(1:3,1:n_map), &
                    this%D_2obs(1:3,1:3, 1:n_map) )
                    
        this%th_map     = th_map
        this%ph_map     = ph_map
        this%zoom_map   = zoom_map
        this%al_map     = al_map
        this%distance   = distance
        this%mass       = mass_dust
        this%n_bin_map  = n_bin_map
        this%r_in       = r_in
        this%r_ou       = r_ou
        this%n_map      = n_map
        this%no_photon  = no_photon

        DO i_map=1, this%n_map
            this%th_map(i_map) = grad2rad(th_map(i_map))
            this%ph_map(i_map) = grad2rad(ph_map(i_map))
            this%al_map(i_map) = cos(grad2rad(al_map(i_map)))
            
!~             this%dir_xyz(1, i_map) = sin(this%th_map(i_map)) * cos(this%ph_map(i_map))
!~             this%dir_xyz(2, i_map) = sin(this%th_map(i_map)) * sin(this%ph_map(i_map))
!~             this%dir_xyz(3, i_map) = cos(this%th_map(i_map))
!~             ! vector marking the +x-direction in the observers map
!~                 
!~             this%e_x(1, i_map) = -sin(this%ph_map(i_map))
!~             this%e_x(2, i_map) =  cos(this%ph_map(i_map))
!~             this%e_x(3, i_map) =  0.0_r2
!~             
!~             ! vector marking the +y-direction in the observers map
!~             
!~             this%e_y(1, i_map) = cos(this%th_map(i_map)) * (-cos(this%ph_map(i_map)))
!~             this%e_y(2, i_map) = cos(this%th_map(i_map)) * (-sin(this%ph_map(i_map)))
!~             this%e_y(3, i_map) = sin(this%th_map(i_map))

            ! Rotationsmatrix to convert a vector from the global coordinate 
            ! system into the observers one

            ! version of S. Wolf -> consistent to MC3D
            
            ! +x-direction in the observers map
            this%D_2obs(1, 1, i_map) = -sin(this%ph_map(i_map))
            this%D_2obs(1, 2, i_map) =  cos(this%ph_map(i_map))
            this%D_2obs(1, 3, i_map) =  0.0_r2
            ! +y-direction in the observers map
            this%D_2obs(2, 1, i_map) = cos(this%th_map(i_map)) * (-cos(this%ph_map(i_map)))
            this%D_2obs(2, 2, i_map) = cos(this%th_map(i_map)) * (-sin(this%ph_map(i_map)))
            this%D_2obs(2, 3, i_map) = sin(this%th_map(i_map))
            ! +z-direction in the observers map == torwards the observer
            this%D_2obs(3, 1, i_map) = sin(this%th_map(i_map)) * cos(this%ph_map(i_map))
            this%D_2obs(3, 2, i_map) = sin(this%th_map(i_map)) * sin(this%ph_map(i_map))
            this%D_2obs(3, 3, i_map) = cos(this%th_map(i_map))

            ! version of O. Fischer (maybe more intuitive and comparable to other codes)
            !
            ! this%D_2obs(1, 1, i_map) =  cos(this%ph_map(i_map))
            ! this%D_2obs(2, 1, i_map) = -cos(this%th_map(i_map)) * sin(this%ph_map(i_map))
            ! this%D_2obs(3, 1, i_map) = -sin(this%th_map(i_map)) * sin(this%ph_map(i_map))

            ! this%D_2obs(1, 2, i_map) = sin(this%ph_map(i_map))
            ! this%D_2obs(2, 2, i_map) = cos(this%th_map(i_map)) * cos(this%ph_map(i_map))
            ! this%D_2obs(3, 2, i_map) = sin(this%th_map(i_map)) * cos(this%ph_map(i_map))

            ! this%D_2obs(1, 3, i_map) =  0.0_r2
            ! this%D_2obs(2, 3, i_map) = -sin(this%th_map(i_map))
            ! this%D_2obs(3, 3, i_map) =  cos(this%th_map(i_map))

            ! user unit, but should be fixed to [AU]
            this%px_model_length_x(i_map) = this%r_ou/(REAL(this%n_bin_map,KIND=r2)+0.5_r2)/this%zoom_map(i_map)   
            this%px_model_length_y(i_map) = this%r_ou/(REAL(this%n_bin_map,KIND=r2)+0.5_r2)/this%zoom_map(i_map)
        END DO
    END SUBROUTINE InitModel

    SUBROUTINE CloseModel(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%modltype)
        DEALLOCATE( this%th_map, &
                    this%ph_map, &
                    this%al_map, &
                    this%zoom_map, &
                    this%px_model_length_x, &
                    this%px_model_length_y, &
                    this%dir_xyz, &
                    this%e_x, &
                    this%e_y, &
                    this%D_2obs)

    END SUBROUTINE CloseModel

    PURE FUNCTION GetModelType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%modltype)
    END FUNCTION GetModelType


    PURE FUNCTION GetModelName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Model_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%modltype)
    END FUNCTION GetModelName

    PURE FUNCTION ModelInitialized(this) RESULT(i)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(Model_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !--------------------------------------------------------------------!
          i = Initialized_common(this%modltype)
    END FUNCTION ModelInitialized

END MODULE model_type
