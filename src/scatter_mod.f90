! ---
! scattering routines
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

  ! ################################################################################################
  ! Scattering [1]: Determination of the scattering direction
  ! ---
  SUBROUTINE scatter(basics,rand_nr,dust,photon, i_dust_sca )
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                    :: photon
    !--------------------------------------------------------------------------!
    integer, intent(in)                             :: i_dust_sca
    
    !--------------------------------------------------------------------------!
    ! ---
    IF      (dust%aniiso == 1) then         ! * anisotropic scattering *
        CALL miesca(basics,rand_nr,dust,photon, i_dust_sca )
    ELSE IF (dust%aniiso == 2) then         ! * isotropic scattering   *
        !call isosca()
        PRINT *, 'TbD, not implemented yet'
        STOP
    ELSE IF (dust%aniiso == 3) then         ! * henyey-greenstein scattering *
        !call hgsca(i_dust_sca)
        PRINT *, 'TbD, not implemented yet'
        STOP
    end if
  END SUBROUTINE scatter


  ! ################################################################################################
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
  ! lot. Azimuthal angle phi is drawing by lot  basing on linear 
  ! polarization of the incident "photon", expressed by values of  I,Q,U  
  ! and the scattering angle theta.
  ! PHI = 0...360 degr.,  THETA = 0...180 degr.
  ! ---
  SUBROUTINE miesca(basics,rand_nr,dust,photon, i_dust)
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    !TYPE(Grid_TYP),INTENT(INOUT)                     :: grid
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                    :: photon
    !--------------------------------------------------------------------------!
    integer, intent(in)                             :: i_dust
    integer                                          :: lot_th 

    logical       :: hl1
    real(kind=r2) :: HELP,PHI,PHI1,PHIPAR,ANGLE,GAMMA,hd1,hd2, rndx
    !--------------------------------------------------------------------------!
    ! ---
    !                     *** THETA-Determination ***
    ! Determination of the scattering angle  ANGLE  and the index  LOT 
    ! which shows to the accompanying elements of the scattering matrix.
    ! LOT=1 shows to THETA=0 degree and LOT=x shows to THETA=(x-1) * DANG
    ! degree.
    CALL RAN2(rand_nr,rndx)
    lot_th = dust%SCAANG(i_dust, photon%nr_lam, nint(real(dust%nrndpt,kind=r2)*rndx) )
    ANGLE  = (lot_th-1) * dust%D_ANG
    photon%SINTHE = sin(ANGLE)
    photon%COSTHE = cos(ANGLE)

    !                      *** PHI'-Determination ***
    ! Parameter  PHIPAR  declares deviation from an uniform distribution
    ! of the azimuthal angles PHI'. 
    PHIPAR = ( sqrt(photon%stokes(2)**2 + photon%stokes(3)**2) / photon%stokes(1) ) * &
        real(-dust%SME(1,2,i_dust,photon%nr_lam,lot_th) / dust%SME(1,1,i_dust,photon%nr_lam,lot_th), kind=r2)

    CALL RAN2(rand_nr,rndx)

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
                     hl1 = .true.   ! endlos-schleife :-( -> neuen wert auslosen
                     CALL RAN2(rand_nr,rndx)
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
    photon%SIN2PH = 2.0_r2 * photon%SINPHI * photon%COSPHI
    photon%COS2PH = 1.0_r2 - 2.0_r2 * photon%SINPHI**2
  END SUBROUTINE miesca


  ! ################################################################################################
  ! isotropic scattering
  ! ---
!~   subroutine isosca()
!~     use datatype
!~     use var_global
!~     use math_mod
!~ 
!~     implicit none
!~     ! ---
!~     ! theta
!~     call RAN2()
!~     SINTHE = sqrt(4.0_r2 * (rndx - rndx**2))
!~     COSTHE = 1.0_r2 - 2.0_r2*rndx
!~     lot_th = nint( rad2grad(dasico( SINTHE, COSTHE )) )
!~     if (lot_th==0) then
!~        lot_th = 1
!~     end if
!~ 
!~     ! phi
!~     call RAN2()
!~     SINPHI = sin( rndx * PIx2 )
!~     COSPHI = cos( rndx * PIx2 )
!~     SIN2PH = 2.0_r2 * SINPHI * COSPHI
!~     COS2PH = 1.0_r2 - 2.0_r2 * SINPHI**2
!~ 
!~   end subroutine isosca


  ! ################################################################################################
  ! henyey-greenstein scattering
  ! ---
!~   subroutine hgsca(i_dust)
!~     use datatype
!~     use var_global
!~     use math_mod
!~ 
!~     implicit none
!~     integer, intent(in) :: i_dust
!~ 
!~     real(kind=r2) :: ANGLE
!~     ! ---
!~     ! theta
!~     call RAN2()
!~     lot_th = SCAANG( i_dust, nr_lam, nint(real(nrndpt,kind=r2)*rndx) )
!~     ANGLE  = lot_th * D_ANG
!~     SINTHE = sin(ANGLE)
!~     COSTHE = cos(ANGLE)
!~ 
!~     ! phi
!~     call RAN2()
!~     SINPHI = sin( rndx * PIx2 )
!~     COSPHI = cos( rndx * PIx2 )
!~     SIN2PH = 2.0_r2 * SINPHI * COSPHI
!~     COS2PH = 1.0_r2 - 2.0_r2 * SINPHI**2
!~   end subroutine hgsca
END MODULE scatter_mod

