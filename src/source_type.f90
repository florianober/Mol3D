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
! def SOURCE_TYP: This type includes all informations of all sources added
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE source_type
    USE datatype
    USE var_global
    USE dust_type
    USE math_mod
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE


    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
        TYPE SOURCE_TYP
        !----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: Luminosity
        REAL(kind=r2),DIMENSION(1:3)                         :: pos_xyz
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: wave_pdf
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: wave_cdf
        INTEGER                                              :: s_type
        INTEGER                                              :: emission_concept

    END TYPE SOURCE_TYP

    TYPE SOURCES
        TYPE(Common_TYP) :: mtype
        !----------------------------------------------------------------------!
        REAL(kind=r2)                                        :: L_total
        TYPE(SOURCE_TYP),DIMENSION(:),ALLOCATABLE            :: source
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: source_cdf
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: lam
        REAL(kind=r2),DIMENSION(:),ALLOCATABLE               :: d_lam
        
        INTEGER                                              :: n_sources
        INTEGER                                              :: n_lam

    END TYPE SOURCES

    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        SOURCES, &
        ! methods
        InitSources, &
        CloseSources, &
        GetSourcesType, &
        GetSourcesName, &
        SourcesInitialized, &
        AddSources, &
        GetNewSource, &
        GetNewLam, &
        Emission_from_source
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitSources(this, ut, un, dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES)          :: this
        TYPE(Dust_TYP)         :: dust
        INTEGER                :: ut
        
        CHARACTER(LEN=*)       :: un
        
        !----------------------------------------------------------------------!
        INTENT(IN)             :: ut,un
        INTENT(INOUT)          :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
        this%n_sources = 0

        ! we need to know the wavelength range
        ! at the moment we use the range coming along with the dust proberties
        ! we could thing about a more clever/general solution ....

        ALLOCATE(this%lam(1:dust%n_lam))
        ALLOCATE(this%d_lam(1:dust%n_lam))
        this%n_lam = dust%n_lam
        this%lam     = dust%lam
        this%d_lam   = dust%d_lam

    END SUBROUTINE InitSources

    SUBROUTINE AddSources(this, s_type, pos_xyz, T_star, R_star, L_star)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES)                              :: this
        INTEGER                                    :: i_source, s_type, i
        REAL(kind=r2)                              :: L, B
        REAL(kind=r1),OPTIONAL                     :: R_star, L_star
        REAL(kind=r1)                              :: T_star
        REAL(kind=r2),DIMENSION(1:3)               :: pos_xyz
        TYPE(SOURCE_TYP),DIMENSION(:), ALLOCATABLE :: source_tmp

        !----------------------------------------------------------------------!
        INTENT(INOUT)         :: this
        INTENT(IN)            :: pos_xyz, s_type, R_star,T_star
        !----------------------------------------------------------------------!
        i_source = this%n_sources + 1
        B = 0.0_r2
        L = 0.0_r2

        IF (this%n_sources > 0) THEN
            ALLOCATE(source_tmp(1:i_source))
            source_tmp(1:this%n_sources) = this%source
            DEALLOCATE(this%source)
            DEALLOCATE(this%source_cdf)
            ALLOCATE(this%source(1:i_source), &
                     this%source_cdf(1:i_source))
            this%source_cdf = 0.0_r2
            this%source(1:this%n_sources) = source_tmp
            DEALLOCATE(source_tmp)

        ELSE
            ALLOCATE(this%source(1:i_source), &
                     this%source_cdf(1:i_source))
        END IF

        IF (s_type == 1) THEN 
            ! source is a point star
            !
            ALLOCATE(this%source(i_source)%wave_cdf(1:this%n_lam))
            ALLOCATE(this%source(i_source)%wave_pdf(1:this%n_lam))
            this%source(i_source)%wave_cdf(:) = 0.0_r2
            this%source(i_source)%wave_pdf(:) = 0.0_r2
            B = integ1(this%lam(:),planck(T_star,this%lam(:)) , 1, this%n_lam)
            this%source(i_source)%wave_pdf(:) = planck(T_star,this%lam(:))/B
            DO i = 1, this%n_lam
                IF ( i == 1) THEN
                    this%source(i_source)%wave_cdf(i) =  this%source(i_source)%wave_pdf(i)
                ELSE
                    this%source(i_source)%wave_cdf(i) = this%source(i_source)%wave_cdf(i-1) + &
                                integ1(this%lam(:),this%source(i_source)%wave_pdf(:) , i-1, i)
                END IF
            END DO

            IF (present(L_star)) THEN
                L = L_star
            ELSE
                IF (present(R_star)) THEN
                    L = integ1(this%lam(:),PI * 4.0_r2 * PI * R_star**2 *      &
                        planck(T_star,this%lam(:)) , 1, this%n_lam)
                ELSE
                    print *,'ERROR: star radius not given'
                    return
                END IF
            END IF
            this%source(i_source)%emission_concept = 1 ! (isotropic)
        ELSE IF (s_type == 2) THEN
            ! source is a point source, but given is the Luminosity and
            ! effective Temperature
            !
            ALLOCATE(this%source(i_source)%wave_cdf(1:this%n_lam))
            ALLOCATE(this%source(i_source)%wave_pdf(1:this%n_lam))
            this%source(i_source)%wave_cdf(:) = 0.0_r2
            this%source(i_source)%wave_pdf(:) = 0.0_r2
            
            B = integ1(this%lam(:),planck(T_star,this%lam(:)) , 1, this%n_lam)
            this%source(i_source)%wave_pdf(:) = planck(T_star,this%lam(:))/B
            DO i = 1,this%n_lam
                IF ( i == 1) THEN
                    this%source(i_source)%wave_cdf(i) = 0.0_r2 
                ELSE 
                    this%source(i_source)%wave_cdf(i) = this%source(i_source)%wave_cdf(i-1) + &
                                integ1(this%lam(:),this%source(i_source)%wave_pdf(:) , i-1, i)
                END IF
            END DO
            IF (present(L_star)) THEN
                L = L_star
            ELSE
                print *,'ERROR: star luminosity not given'
                return
            END IF
            this%source(i_source)%emission_concept = 1 ! (isotropic)
        ELSE
            print *, 'ERROR: source type not implemented'
            stop
        END IF
        this%source(i_source)%pos_xyz     = pos_xyz
        this%source(i_source)%s_type      = s_type
        this%n_sources                    = i_source
        this%source(i_source)%Luminosity  = L
        ! update source cdf
        this%L_total = sum(this%source(:)%Luminosity)
        DO i = 1,i_source
            IF (i == 1) THEN
                this%source_cdf(i) = this%source(i)%Luminosity /this%L_total
            ELSE
                this%source_cdf(i) = this%source_cdf(i-1) +                    &
                                     this%source(i)%Luminosity /this%L_total
            END IF
        END DO

    END SUBROUTINE AddSources

    SUBROUTINE Emission_from_source(source, rand_nr,                           &
                                    SINPHI, COSPHI, SINTHE, COSTHE)
        USE randgen_type
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCE_TYP), INTENT(IN)                   :: source
        TYPE(Randgen_TYP),INTENT(INOUT)                :: rand_nr

        REAL(kind=r2), INTENT(OUT)                     :: SINPHI
        REAL(kind=r2), INTENT(OUT)                     :: COSPHI
        REAL(kind=r2), INTENT(OUT)                     :: SINTHE
        REAL(kind=r2), INTENT(OUT)                     :: COSTHE
        
        REAL(kind=r2)                                  :: rndx1, rndx2
        !----------------------------------------------------------------------!
        ! We need to calculate the
        IF ( source%emission_concept == 1) THEN
            ! (isotropic)
            CALL GetNewRandomNumber(rand_nr, rndx1)
            CALL GetNewRandomNumber(rand_nr, rndx2)
            CALL isotropic_sphere(rndx1, rndx2, SINPHI, COSPHI, SINTHE, COSTHE)
        ! TbD: more emission concepts, expanded sources, dark limb darkening,
        !      sunspots, ....
        ELSE
            print *, "Type of emission concept for the selected source not known"
            stop
        END IF
    END SUBROUTINE Emission_from_source

    SUBROUTINE CloseSources(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        DEALLOCATE(this%source)
        
    END SUBROUTINE CloseSources


    PURE FUNCTION GetSourcesType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)      :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetSourcesType


    PURE FUNCTION GetSourcesName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN)       :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetSourcesName
    
    FUNCTION GetNewSource(this,rndx) RESULT(i)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER                   :: i
        !----------------------------------------------------------------------!
        IF (this%n_sources == 1) THEN
            i = 1
        ELSE
            i = binary_search(rndx, this%source_cdf)+1
        END IF

    END FUNCTION GetNewSource
    
    FUNCTION GetNewLam(this,i_source,rndx) RESULT(i)
        
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        REAL(kind=r2), INTENT(IN) :: rndx
        INTEGER, INTENT(IN)       :: i_source
        INTEGER                   :: i
        !----------------------------------------------------------------------!
        i = binary_search(rndx, this%source(i_source)%wave_cdf) + 1 
        
    END FUNCTION GetNewLam

    PURE FUNCTION SourcesInitialized(this) RESULT(i)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(SOURCES), INTENT(IN) :: this
        LOGICAL :: i
        !----------------------------------------------------------------------!
        i = Initialized_common(this%mtype)
          
        IF (i .and. this%n_sources == 0) THEN
            i = .False.
        END IF
    END FUNCTION SourcesInitialized

END MODULE source_type
