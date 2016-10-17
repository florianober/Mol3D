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
    USE roots

    IMPLICIT NONE
!~     INTERFACE
!~         PURE SUBROUTINE func(x,fx,plist)
!~             USE datatype
!~             IMPLICIT NONE
!~             REAL(kind=r2), INTENT(IN)  :: x
!~             REAL(kind=r2), INTENT(IN), DIMENSION(:), OPTIONAL :: plist
!~             REAL(kind=r2), INTENT(OUT) :: fx
!~         END SUBROUTINE func
!~     END INTERFACE

    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!

    PUBLIC  :: scatter, trafo

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
            CALL isosca(rand_nr, photon)

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
        real(kind=r2) :: HELP,PHI,PHI1,PHIPAR,ANGLE,GAMMA_ANGLE,hd1,hd2, rndx


        REAL(kind=r2) :: root, xm, dx_rel, dx_acc, plist(2), HELP1
        INTEGER       :: iter, error, iter2
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
        HELP = rndx * basics%PIx4            ! HELP = rndx * 4.0_r2 * PI
        plist(1) = PHIPAR
        plist(2) = HELP
        
        CALL GetRoot_BrentDekker(func,-0.0_r2, 4.0_r2*PI,root,error,plist,HELP,iter)
        HELP1 = root * 0.5_r2
        PHI = HELP1

!        IF (abs(PHIPAR) < 0.1_r2) THEN
!            PHI = rndx * basics%PIx2                ! PHI = rndx * 2.0_r2 * PI
!        ELSE
!            hl1 = .true.
!            DO WHILE (hl1)
!                HELP = rndx * basics%PIx4            ! HELP = rndx * 4.0_r2 * PI
!                PHI  = HELP  +  PHIPAR * sin(HELP)
!
!                ! Calculation of PHI by an iterative solution of Kepler's
!                ! equation. The iteration will end if DABS(PHI-PHI1)<0.0175
!                ! (about 1 degree).
!                PHI1 = 0.0_r2
!                DO
!                    hd1  = abs(PHI - PHI1)
!                    PHI1 = PHI
!                    PHI  = HELP  +  PHIPAR * sin(PHI1)
!                    hd2  = abs(PHI - PHI1)
!                    iter2 = iter2 + 1
!
!                    if(abs(PHI - PHI1) <= 0.0175_r2) then 
!                        hl1 = .false.
!                        exit
!                    else
!                        if (abs(hd1-hd2) < con_eps) then
!                            CALL GetNewRandomNumber(rand_nr,rndx)
!                            exit
!                        endif
!                    endif
!                enddo
!            enddo
!            PHI = PHI / 2.0_r2
!        END IF
        !IF (PHIPAR > 0.14) THEN
        !    print *, iter, iter2
        !    print *, abs(HELP1-PHI)/PHI
        !    print *, ''
        !END IF
        !IF (PHIPAR > 0.2) STOP
        ! Between the direction angle phi and the angle phi'
        ! calculated by inverse of distribution function exists the   
        ! following connection:  phi = phi' + 180 deg - GAMMA_ANGLE ; 
        ! tan(2GAMMA_ANGLE) = Uein/Qein.
        
        GAMMA_ANGLE = 0.0_r2 ! not sure if needed
        if (photon%stokes(2) /= 0.0_r2) then
            GAMMA_ANGLE  = 0.5_r2 * atan2(photon%stokes(3),photon%stokes(2))
            if  (photon%stokes(3) < 0.0_r2) then
                GAMMA_ANGLE = PI + GAMMA_ANGLE
            endif
        else
            if     (photon%stokes(3) < 0.0_r2) then
                GAMMA_ANGLE = basics%PIx34    ! GAMMA_ANGLE = PI * 3.0_r2/4.0_r2
            elseif (photon%stokes(3) > 0.0_r2) then
                GAMMA_ANGLE = basics%PI4      ! GAMMA_ANGLE = PI / 4.0_r2
            endif
        endif
        PHI = PI - GAMMA_ANGLE + PHI
        photon%SINPHI = sin(PHI)
        photon%COSPHI = cos(PHI)

    END SUBROUTINE miesca


    ! ##########################################################################
    ! isotropic scattering
    ! ---
    SUBROUTINE isosca(rand_nr, photon)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
        !----------------------------------------------------------------------!
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
        ANGLE  = (photon%lot_th -1)* dust%D_ANG
        photon%SINTHE = sin(ANGLE)
        photon%COSTHE = cos(ANGLE)

        ! phi
        CALL GetNewRandomNumber(rand_nr, rndx)
        photon%SINPHI = sin( rndx * PI * 2.0_r2)
        photon%COSPHI = cos( rndx * PI * 2.0_r2)

    END SUBROUTINE hgsca

    ! ##########################################################################
    ! Scattering causes a change of stokes vector (I,Q,U,V).  This
    ! transformation occurs in two steps. 
    ! At first, stokes vector is transformed into the new r,l-plane
    ! after the PHI-rotation. 
    ! Therefore values of SIN2PH, COS2PH are necessary. 
    ! The second step contains the transformation of the stokes vector with the 
    ! according  (index is pointed by LOT) scattering matrix (S11,S12,S33,S34).
    ! ---
    SUBROUTINE trafo(dust, photon, i_dust)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN)                          :: dust
        TYPE(PHOTON_TYP), INTENT(INOUT)                     :: photon
        INTEGER, INTENT(IN)                                 :: i_dust
        REAL(kind=r2)                                       :: QHELP, i_1
        !----------------------------------------------------------------------!
        ! ---
        i_1       = photon%stokes(1)
        ! scattering
        ! Mathematical positive rotation of r,l-plane around
        !              p-axis by the angle PHI
        QHELP     =  photon%COS2PH * photon%stokes(2) +                        &
                     photon%SIN2PH * photon%stokes(3)
        photon%stokes(3) =  - photon%SIN2PH * photon%stokes(2) +               &
                              photon%COS2PH * photon%stokes(3)
        photon%stokes(2) =  QHELP
        ! Mathematical negative rotation around r-axis by THETA transformes
        ! the Stokes vector (I,Q,U,V) with the scattering matrix.
        ! LOT  points to the scattering matrix elements of photon angle which
        ! was determined by throwing dice before.
        photon%stokes(:) = matmul(dust%SME( :, :, i_dust, photon%nr_lam,       &
                                  photon%lot_th ), photon%stokes(:) )

        ! normalize the stokes vector
        photon%stokes(:) = photon%stokes(:) * i_1/photon%stokes(1)

    END SUBROUTINE trafo

    PURE SUBROUTINE func(x,fx,plist)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        REAL(kind=r2), INTENT(IN)  :: x
        REAL(kind=r2), INTENT(IN), DIMENSION(:), OPTIONAL :: plist
        REAL(kind=r2), INTENT(OUT) :: fx
        !----------------------------------------------------------------------!
        fx = x - plist(1) * SIN(x) - plist(2)
    END SUBROUTINE func
    
END MODULE scatter_mod

