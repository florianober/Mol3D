! ---
! start photon
! ---
MODULE start_mod
    USE datatype
    USE var_global
    
    USE grid_type
    USE dust_type
    USE model_type
    USE randgen_type
    USE photon_type
    USE basic_type
    USE fluxes_type
    USE source_type
    
    USE math_mod
  
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC  :: start_photon, start_grain

contains

  ! ############################################################################
  ! start photon from random source
  ! ---
  SUBROUTINE start_photon(basics, grid, model, rand_nr, dust,photon, sources_in)
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(MODEL_TYP),INTENT(IN)                       :: model
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(SOURCES),INTENT(IN)                         :: sources_in
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!

    integer                                          :: i_lam
    integer                                          :: i_source
    REAL(kind=r2)                                    :: rndx
    !--------------------------------------------------------------------------!
    i_lam = 1
    ! 1. get new source to start from
    CALL RAN2(rand_nr,rndx)
    i_source = GetNewSource(sources_in, rndx)

    ! 2. direction
    CALL emission_from_source(sources_in%source(i_source), rand_nr,            &
                              photon%SINPHI,                                   &
                              photon%COSPHI,                                   &
                              photon%SINTHE,                                   &
                              photon%COSTHE)

    ! calculate SIN2PH/COS2PH for the stokes rotation
    CALL update_angle(photon)

    ! 3. get new wavelength
    IF (.not. photon%fixed_lam) THEN
        DO
            CALL RAN2(rand_nr,rndx)
            IF (rndx .le. 1.0e-300_r2) CYCLE
            IF (rndx .gt. maxval(sources_in%source(i_source)%wave_cdf) ) THEN
                i_lam = -1
                photon%kill = .True.
                photon%inside  = .False.
            ELSE
                i_lam = GetNewLam(sources_in, i_source, rndx)
                photon%kill = .False.
                photon%inside  = .True.
            END IF
            EXIT
        END DO
    END IF
    ! 4. starting point = location of source
    photon%pos_xyz(:) = sources_in%source(i_source)%pos_xyz

    ! 5. set corresponding cell number
    photon%nr_cell = get_cell_nr(grid, photon%pos_xyz(:))

    ! ---
    ! [B] initialize new photon parameters
    ! stokes vector
    ! starting with unpolarized photon
    photon%stokes(:) = (/1.0_r2, 0.0_r2, 0.0_r2, 0.0_r2/)

    IF (photon%fixed_lam) THEN
        ! when using a fixed wavelength, we need to weightend the photons
        ! energy with the source pdf
        photon%energy = sources_in%source(i_source)%wave_pdf(photon%nr_lam) *  &
                        sources_in%L_total /                                   &
                        model%no_photon
    ELSE
        photon%nr_lam         = i_lam          ! initial wavelength of photon

        ! set energy of photon package   
        photon%energy         = sources_in%L_total /                           &
                                model%no_photon
    END IF
    ! photon inside model space

    photon%D_2global(:,:) = basics%mat_ident3(:,:)  ! initialize rotation matrix

    ! ---
    ! [C] apply rotation matrix
    CALL vecmat(photon)
    ! ---
    ! last point of interaction = start position

    photon%pos_xyz_li(:) = photon%pos_xyz(:)
    photon%last_interaction_type = 'N'
  END SUBROUTINE start_photon


  ! ############################################################################
  ! emission from dust grain at current position
  ! ---
  SUBROUTINE start_grain(basics, rand_nr, photon)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!
    REAL(kind=r2)                                    :: rndx1, rndx2
    !--------------------------------------------------------------------------!
    ! ---
    ! [A] determine starting point & direction
    ! 1. point of emission:
    ! current position

    ! 2. direction of emission
    ! (isotropic)
    CALL RAN2(rand_nr, rndx1)
    CALL RAN2(rand_nr, rndx2)
    CALL isotropic_sphere(rndx1, rndx2,                                        &
                              photon%SINPHI,                                   &
                              photon%COSPHI,                                   &
                              photon%SINTHE,                                   &
                              photon%COSTHE)

    ! calculate SIN2PH/COS2PH for the stokes rotation
    CALL update_angle(photon)
    ! ---
    ! [B] initialize new photon parameters
    ! unpolarized stokes vector
    photon%stokes(:) = (/1.0_r2, 0.0_r2, 0.0_r2, 0.0_r2/)

    ! photon inside model space
    photon%inside    = .true.

    ! initialize rotation matrix
    !photon%D(:,:)    = basics%mat_ident3(:,:)
    photon%D_2global(:,:)    = basics%mat_ident3(:,:)

    ! ---
    ! [C] apply rotation matrix
    CALL vecmat(photon)
    
  END SUBROUTINE start_grain

END MODULE start_mod

