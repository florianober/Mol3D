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
! scattering routines (this module is mainly imported from MC3D)
! ---
MODULE scatter_mod
    USE datatype
    USE var_global
    
    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE photon_type
    USE basic_type
    USE fluxes_type

    USE math_mod

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!

    PUBLIC  :: scatter

CONTAINS

    ! ##########################################################################
    ! Scattering [1]: Determination of the scattering direction
    ! ---
    SUBROUTINE scatter(basics, rand_nr, dust, photon, i_dust_sca)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP), INTENT(IN)                       :: basics
        TYPE(Randgen_TYP), INTENT(INOUT)                  :: rand_nr
        TYPE(Dust_TYP), INTENT(IN)                        :: dust

        TYPE(PHOTON_TYP), INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)                               :: i_dust_sca
        !----------------------------------------------------------------------!
        ! ---
        IF      (dust%aniiso == 1) then         ! anisotropic mie scattering
            CALL miesca(basics, rand_nr, dust, photon, i_dust_sca)

        ELSE IF (dust%aniiso == 2) then         ! isotropic scattering
            CALL isosca(rand_nr, photon, i_dust_sca)

        ELSE IF (dust%aniiso == 3) then         ! henyey-greenstein scattering
            CALL hgsca(rand_nr, dust, photon, i_dust_sca)

        END IF
        ! update the SIN2PH and COS2PH values
        CALL update_angle(photon)

    END SUBROUTINE scatter

    ! ##########################################################################
    ! SCATTERING BY SPHERICAL PARTICLES OR ELECTRONS
    !
    ! The scattering direction (SINPHI,COSPHI,SINTHE,COSTHE,SIN2PH,COS2PH)
    ! is determined with the help of the Monte Carlo Method.
    ! In the vector  SCAANG  indexes pointing to the scattering directions
    ! are distributed according to the characteristics of scattering
    ! (SCAANG was calculated in SUBROUTINE BEGIN).
    ! By the help of SCAANG for a discret number of scattering angles 
    ! theta between 0 and 180 degree with a distance of  DANG  degree are 
    ! drawing by lots.  LOT adresses the number of direction drawing by 
    ! lot. Azimuthal angle phi is drawing by lot basing on linear 
    ! polarization of the incident "photon", expressed by values of  I,Q,U  
    ! and the scattering angle theta.
    ! PHI = 0...360 degr.,  THETA = 0...180 degr.
    ! ---
  SUBROUTINE miesca(basics, rand_nr, dust, photon, i_dust)

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP),INTENT(IN)                       :: basics
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Dust_TYP),INTENT(IN)                        :: dust
        
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        integer, intent(in)                              :: i_dust

        logical       :: hl1
        real(kind=r2) :: HELP,PHI,PHI1,PHIPAR,ANGLE,GAMMA,hd1,hd2, rndx
        !----------------------------------------------------------------------!
        ! ---
        !                     *** THETA-Determination ***
        ! Determination of the scattering angle  ANGLE  and the index  LOT 
        ! which shows to the accompanying elements of the scattering matrix.
        ! LOT=1 shows to THETA=0 degree and LOT=x shows to THETA=(x-1) * DANG
        ! degree.
        CALL GetNewRandomNumber(rand_nr, rndx)
        photon%lot_th = dust%SCAANG(i_dust, photon%nr_lam,                     &
                                    nint(real(dust%nrndpt,kind=r2)*rndx) )
        ANGLE  = (photon%lot_th-1) * dust%D_ANG
        photon%SINTHE = sin(ANGLE)
        photon%COSTHE = cos(ANGLE)

        !                      *** PHI'-Determination ***
        ! Parameter  PHIPAR  declares deviation from an uniform distribution
        ! of the azimuthal angles PHI'. 
        PHIPAR = (sqrt(photon%stokes(2)**2 + photon%stokes(3)**2) /            &
                  photon%stokes(1)) * &
                  real(-dust%SME(1, 2, i_dust, photon%nr_lam,photon%lot_th) /  &
                  dust%SME(1, 1, i_dust, photon%nr_lam, photon%lot_th), kind=r2)

        CALL GetNewRandomNumber(rand_nr, rndx)

        if (abs(PHIPAR) < 0.1_r2) then
           PHI = rndx * basics%PIx2                ! PHI = rndx * 2.0_r2 * PI
        else
           do
              HELP = rndx * basics%PIx4            ! HELP = rndx * 4.0_r2 * PI
              PHI  = HELP  +  PHIPAR * sin(HELP)

              ! Calculation of PHI by an iterative solution of Kepler's
              ! equation. The iteration will end if DABS(PHI-PHI1)<0.0175
              ! (about 1 degree). 
              do
                 PHI1 = PHI
                 hd1  = abs(PHI - PHI1)
                 PHI  = HELP  +  PHIPAR * sin(PHI1)
                 hd2  = abs(PHI - PHI1)

                 if(abs(PHI - PHI1) <= 0.0175_r2) then 
                    hl1 = .false.
                    exit
                 else
                    if ((hd1-hd2) < 1.0e-15_r2) then
                         hl1 = .true.   ! endlos-schleife -> neuen wert auslosen
                         CALL GetNewRandomNumber(rand_nr,rndx)
                         exit
                    endif
                 endif
              enddo
              if (.not. hl1) then
                  exit
              endif
           enddo
           PHI = PHI / 2.0_r2
           ! Between the direction angle phi and the angle phi'
           ! calculated by inverse of distribution function exists the   
           ! following connection:  phi = phi' + 180 deg - gamma ; 
           ! gamma=Uein/Qein.
           GAMMA = 0.0_r2 ! not sure if needed
           if (photon%stokes(2) /= 0.0_r2) then
              GAMMA  = 0.5_r2 * atan2(photon%stokes(3),photon%stokes(2))
              if  (photon%stokes(3) < 0.0_r2) then
                 GAMMA = PI + GAMMA
              endif
           else
              if     (photon%stokes(3) < 0.0_r2) then
                 GAMMA = basics%PIx34    ! GAMMA = PI * 3.0_r2/4.0_r2
              elseif (photon%stokes(3) > 0.0_r2) then
                 GAMMA = basics%PI4      ! GAMMA = PI / 4.0_r2
              endif
           endif
           PHI = PI - GAMMA + PHI
        endif

        photon%SINPHI = sin(PHI)
        photon%COSPHI = cos(PHI)

    END SUBROUTINE miesca


    ! ##########################################################################
    ! isotropic scattering
    ! ---
    SUBROUTINE isosca(rand_nr, photon, i_dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        integer, intent(in)                              :: i_dust
        REAL(kind=r2)                                    :: rndx1, rndx2
        !----------------------------------------------------------------------!
        ! ---
        CALL GetNewRandomNumber(rand_nr, rndx1)
        CALL GetNewRandomNumber(rand_nr, rndx2)
        CALL isotropic_sphere(rndx1, rndx2, photon%SINPHI, photon%COSPHI,      &
                                            photon%SINTHE, photon%COSTHE)
        photon%lot_th = nint( rad2grad(dasico( photon%SINTHE, photon%COSTHE )) )
        IF (photon%lot_th==0) THEN
           photon%lot_th = 1
        END IF

    END SUBROUTINE isosca

    ! ##########################################################################
    ! henyey-greenstein scattering
    ! NOT tested
    ! ---
    SUBROUTINE hgsca(rand_nr, dust, photon, i_dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Dust_TYP),INTENT(IN)                        :: dust

        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
        integer, intent(in)                              :: i_dust
        REAL(kind=r2)                                    :: rndx
        REAL(kind=r2)                                    :: ANGLE
        !----------------------------------------------------------------------!

        ! ---
        ! theta
        CALL GetNewRandomNumber(rand_nr, rndx)
        photon%lot_th = dust%SCAANG(i_dust, photon%nr_lam,                     &
                                    nint(real(dust%nrndpt, kind=r2)*rndx))
        ANGLE  = photon%lot_th * dust%D_ANG
        photon%SINTHE = sin(ANGLE)
        photon%COSTHE = cos(ANGLE)

        ! phi
        CALL GetNewRandomNumber(rand_nr, rndx)
        photon%SINPHI = sin( rndx * PI * 2.0_r2)
        photon%COSPHI = cos( rndx * PI * 2.0_r2)

    END SUBROUTINE hgsca

END MODULE scatter_mod

