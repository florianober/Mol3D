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
! def basic_TYP
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE basic_type
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE
    ! ---
    ! special variables
    ! ---
    !              ( 1 0 0 )                      ( 1 0 0 0 )
    ! mat_ident3 = ( 0 1 0 )         mat_ident4 = ( 0 1 0 0 )
    !              ( 0 0 1 )                      ( 0 0 1 0 )
    !                                             ( 0 0 0 1 )
    ! PIx2  = PI * 2
    ! PIx4  = PI * 4
    ! PIx34 = PI * 3/4
    ! PI2   = PI / 2
    ! PI4   = PI / 4
    ! PI2x4 = PI^2 * 4

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Basic_TYP
        TYPE(Common_TYP) :: mtype                                              !
        !----------------------------------------------------------------------!
        INTEGER          :: pnamelen
        INTEGER          :: num_core
        
        REAL(kind=r2)    :: PIx2
        REAL(kind=r2)    :: PIx4
        REAL(kind=r2)    :: PIx34
        REAL(kind=r2)    :: PI2
        REAL(kind=r2)    :: PI4
        REAL(kind=r2)    :: PI2x4
                
        ! store all variables from file
        INTEGER          :: n_tem
        
        REAL(kind=r1)    :: d_tem
        REAL(kind=r1)    :: t_dust_max
        REAL(kind=r1)    :: t_dust_min
        

        REAL(kind=r1), dimension(1:3,1:3) :: mat_ident3
        REAL(kind=r1), dimension(1:4,1:4) :: mat_ident4
        
        
        CHARACTER(len=256)                :: pronam
        CHARACTER(len=256)                :: pronam_old
        CHARACTER(len=256)                :: path_results
        CHARACTER(len=256)                :: input_file
        CHARACTER(len=256)                :: new_input_file
        CHARACTER(len=20)                 :: randgen_name
        CHARACTER(len=20)                 :: velocity_type
        CHARACTER(len=20)                 :: temperature_type
        
        LOGICAL                           :: do_MC_temperature
        LOGICAL                           :: old_model
        LOGICAL                           :: do_raytr
        LOGICAL                           :: do_continuum_raytrace
        LOGICAL                           :: do_continuum_mc
        LOGICAL                           :: do_peel_off
        LOGICAL                           :: do_velo_ch_map

    END TYPE Basic_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Basic_TYP, &
        ! methods
        InitBasic, &
        CloseBasic, &
        GetSimuType, &
        GetName, &
        BasicInitialized, &
        Getproname
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitBasic(this, ut, un, pname, presult,                         &
                         do_MC_temperature, old_pname, old_model,              &
                         do_continuum_raytrace, do_continuum_mc,               &
                         do_velo_ch_map, do_peel_off,                          &
                         n_tem, t_dust_min, t_dust_max, num_core,              &
                         input_file, new_input_file, randgen_name)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP)       :: this
        
        INTEGER               :: ut
        INTEGER               :: n_tem
        INTEGER               :: num_core
        
        CHARACTER(LEN=*)      :: un
        CHARACTER(LEN=*)      :: pname
        CHARACTER(LEN=*)      :: old_pname
        CHARACTER(LEN=*)      :: presult
        CHARACTER(len=*)      :: input_file
        CHARACTER(len=*)      :: new_input_file
        CHARACTER(len=*)      :: randgen_name
        
        REAL(kind=r1)         :: t_dust_max
        REAL(kind=r1)         :: t_dust_min
        LOGICAL               :: old_model
        
        LOGICAL               :: do_MC_temperature
        LOGICAL               :: do_continuum_raytrace
        LOGICAL               :: do_continuum_mc
        LOGICAL               :: do_velo_ch_map
        LOGICAL               :: do_peel_off

        !----------------------------------------------------------------------!
        INTENT(IN)          :: ut, un, presult, pname, num_core,               &
                               n_tem, t_dust_max, t_dust_min,                  &
                               do_MC_temperature, old_model, old_pname,        &
                               do_continuum_raytrace, do_continuum_mc,         &
                               do_velo_ch_map, do_peel_off, randgen_name
        INTENT(INOUT)       :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
    
        this%PIx2  = PI * 2.0_r2
        this%PIx4  = PI * 4.0_r2
        this%PIx34 = PI * 3.0_r2/4.0_r2
        this%PI2   = PI / 2.0_r2
        this%PI4   = PI / 4.0_r2
        this%PI2x4 = PI**2 * 4.0_r2
        
        this%num_core = num_core

        ! identity matrix [3x3]
        this%mat_ident3(1,1) = 1.0
        this%mat_ident3(1,2) = 0.0
        this%mat_ident3(1,3) = 0.0
           
        this%mat_ident3(2,1) = 0.0
        this%mat_ident3(2,2) = 1.0
        this%mat_ident3(2,3) = 0.0
           
        this%mat_ident3(3,1) = 0.0
        this%mat_ident3(3,2) = 0.0
        this%mat_ident3(3,3) = 1.0
        
        ! identity matrix [4x4]      
        this%mat_ident4(1,1) = 1.0_r2
        this%mat_ident4(1,2) = 0.0_r2
        this%mat_ident4(1,3) = 0.0_r2
        this%mat_ident4(1,4) = 0.0_r2
           
        this%mat_ident4(2,1) = 0.0
        this%mat_ident4(2,2) = 1.0
        this%mat_ident4(2,3) = 0.0
        this%mat_ident4(2,4) = 0.0
           
        this%mat_ident4(3,1) = 0.0
        this%mat_ident4(3,2) = 0.0
        this%mat_ident4(3,3) = 1.0
        this%mat_ident4(3,4) = 0.0
           
        this%mat_ident4(4,1) = 0.0
        this%mat_ident4(4,2) = 0.0
        this%mat_ident4(4,3) = 0.0
        this%mat_ident4(4,4) = 1.0

        this%n_tem = n_tem
        this%t_dust_max = t_dust_max
        this%t_dust_min = t_dust_min
        this%do_MC_temperature   = do_MC_temperature
        this%old_model  = old_model

        this%do_continuum_raytrace  = do_continuum_raytrace
        this%do_continuum_mc  = do_continuum_mc
        this%do_peel_off = do_peel_off
        this%do_velo_ch_map = do_velo_ch_map

        this%randgen_name = randgen_name
        this%velocity_type = 'ANALYTICAL_KEPLERIAN'   ! possible values: ANALYTICAL_KEPLERIAN, INTERPOLATE
        this%temperature_type = 'MC'       ! possible values: MC, ANALYTICAL

        IF (do_continuum_raytrace .or. do_velo_ch_map) THEN
            this%do_raytr = .True.
        ELSE
            this%do_raytr = .False.
        END IF

        this%d_tem = (this%t_dust_max - this%t_dust_min) /                     &
                     (real(this%n_tem, kind=r2)-1)
        
        this%pronam = pname    
        this%pronam_old = old_pname    
        this%path_results = presult
        this%pnamelen = LEN_TRIM(pname)
        this%input_file = input_file
        this%new_input_file = new_input_file
        
        
    END SUBROUTINE InitBasic


    SUBROUTINE CloseBasic(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
    END SUBROUTINE CloseBasic


    PURE FUNCTION GetSimuType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSimuType


    PURE FUNCTION GetName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetName

    PURE FUNCTION BasicInitialized(this) RESULT(i)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(Basic_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !--------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
    END FUNCTION BasicInitialized
    
    FUNCTION Getproname(this) RESULT(pname)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(Basic_TYP), INTENT(IN) :: this
          CHARACTER(len=this%pnamelen) :: pname
          !--------------------------------------------------------------------!
          pname = this%pronam
    END FUNCTION Getproname
    

End Module basic_type
