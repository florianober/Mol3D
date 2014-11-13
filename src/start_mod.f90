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
    PUBLIC  :: start_photon, start_grain!, start_cell, trafo

contains

  ! ################################################################################################
  ! start photon from random source
  ! ---
  SUBROUTINE start_photon(basics,grid,model,rand_nr,fluxes,dust,photon, sources_in)
    
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(MODEL_TYP),INTENT(IN)                       :: model
    TYPE(Grid_TYP),INTENT(IN)                        :: grid
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(Dust_TYP),INTENT(IN)                        :: dust
    TYPE(FLUXES_TYP),INTENT(IN)                      :: fluxes
    TYPE(SOURCES),INTENT(IN)                         :: sources_in
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                   :: photon
    !--------------------------------------------------------------------------!

    integer                                          :: i_lam
    integer                                          :: i_source
    REAL(kind=r2)                                    :: rndx
    !--------------------------------------------------------------------------!

    IF (basics%concept_ps == 1) THEN
        ! concept_ps==1: point source, isotropic emission
        ! 1. direction (isotropic)
        CALL RAN2(rand_nr, rndx)
        
        photon%SINPHI = sin( rndx * basics%PIx2 )
        photon%COSPHI = cos( rndx * basics%PIx2 )
        photon%SIN2PH = 2.0_r2 * photon%SINPHI * photon%COSPHI
        photon%COS2PH = 1.0_r2 - 2.0_r2 * photon%SINPHI**2
        
        CALL RAN2(rand_nr, rndx)
        
        photon%SINTHE = sqrt(4.0_r2 * (rndx - rndx**2))
        photon%COSTHE = 1.0_r2 - 2.0_r2*rndx
        
        ! 2. get new source
        CALL RAN2(rand_nr,rndx)
        i_source = GetNewSource(sources_in,rndx)
        ! 3. get new wavelength
        CALL RAN2(rand_nr,rndx)
        i_lam = GetNewLam(sources_in,i_source,rndx)
        
        ! 4. starting point = location of source
        photon%pos_xyz(:) = sources_in%source(i_source)%pos_xyz

        ! 5. star is located in inner dust-free region
        photon%nr_cell = get_cell_nr(grid,photon%pos_xyz(:))
        
    ELSE
           print *, "Type of emission concept for primary source not known [concept_ps]"
           print *, "in subroutine start_prim in primary_temp."
           stop
    END IF

    ! ---
    ! [B] initialize new photon parameters
    ! stokes vector
    photon%stokes(:) = fluxes%stokes_ini(:)                    !starting with unpolarized photon
        
    photon%nr_lam         = i_lam                              ! initial wavelength of photon    
!~     photon%current_albedo = dust%albedo(:,i_lam)
    
!~     photon%energy         = dust%c_in_star(i_lam)
!~     photon%energy         = sources_in%source(i_source)%Luminosity/&
!~                            (model%no_photon*sources_in%source_cdf(i_source))    ! set energy of photon package
    photon%energy         = sources_in%L_total / &
                            model%no_photon   ! set energy of photon package
    
    
    
    photon%inside         = .true.                             ! photon inside model space
    photon%D(:,:)         = basics%mat_ident3(:,:)             ! initialize rotation matrix
        
    ! ---
    ! [C] apply rotation matrix
    CALL vecmat(photon)

    ! ---
    ! last point of interaction = start position
    
    photon%pos_xyz_li(:) = photon%pos_xyz(:)

  END SUBROUTINE start_photon
  

  ! ################################################################################################
  ! prepare photon launch from start from cell
  ! : for dust MCRT
  ! ---
!~   subroutine start_cell(i_lam, i_cell)
!~     use datatype
!~     use var_global
!~     use math_mod
!~ 
!~     implicit none
!~     integer, intent(in) :: i_lam, i_cell
!~     
!~     real(kind=r2) :: hd_r1, hd_r2, hd_th1, hd_th2, hd_thm, hd_ph1, hd_ph2
!~     real(kind=r2), dimension(1:3) :: hd_sph
!~     ! ---
!~     ! 1. randomly select position in current cell
!~     ! 1.a: <r>
!~     hd_r1  = co_mx_r( cell_nr2idx(1,i_cell) -1 )     ! r.min
!~     hd_r2  = co_mx_r( cell_nr2idx(1,i_cell)    )     ! r.max
!~     
!~     call RAN2()
!~     hd_sph(1) = (rndx*(hd_r2**3 - hd_r1**3) + hd_r1**3)**(1.0_r2/3.0_r2)
!~ 
!~     ! ---
!~     ! 1.b: <th>
!~     ! tbd: check (was copied only from mc3d.v3)
!~     hd_th1 = co_mx_th( cell_nr2idx(2,i_cell) -1 )     ! th.min
!~     hd_th2 = co_mx_th( cell_nr2idx(2,i_cell)    )     ! th.max
!~     hd_thm = (hd_th1 + hd_th2) / 2.0_r2               ! mean theta
!~     
!~     call RAN2()
!~     hd_sph(2) = acos( cos(hd_th2) + rndx*(cos(hd_th1)-cos(hd_th2)) )
!~     if (hd_thm < 0.0_r2) then
!~        hd_sph(2) = - hd_sph(2)
!~     end if
!~     
!~     ! ---
!~     ! 1.c: <ph>
!~     hd_ph1 = co_mx_ph( cell_nr2idx(3,i_cell) -1 )     ! ph.min
!~     hd_ph2 = co_mx_ph( cell_nr2idx(3,i_cell)    )     ! ph.max
!~     
!~     call RAN2()
!~     hd_sph(3) = hd_ph1 + rndx*abs(hd_ph2-hd_ph1)
!~     
!~     ! ---
!~     ! 1.d: convert (r,th,ph) => pos_xyz
!~     call sp2ca( hd_sph(:), pos_xyz(:) )
!~     
!~     ! 1.e: last point of interaction := start position
!~     pos_xyz_li(:) = pos_xyz(:)
!~ 
!~     ! ---
!~     ! 2. set stokes parameters for unpolarized emission from dust grain
!~     call start_grain()
!~     
!~     ! ---
!~     ! 3. set remaining photon parameters 
!~     photon%nr_lam   = i_lam              ! initial wavelength of photon    
!~     photon%inside   = .true.             ! photon inside model space
!~     photon%energy   = c_in_dust(nr_lam)  ! set energy of photon package
!~     photon%nr_cell  = i_cell             ! set cell number
!~ 
!~   end subroutine start_cell


  ! ################################################################################################
  ! emission from dust grain at current position
  ! ---
  SUBROUTINE start_grain(basics,rand_nr,fluxes,photon)
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP),INTENT(IN)                       :: basics
    TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
    TYPE(FLUXES_TYP),INTENT(IN)                      :: fluxes
    
    TYPE(PHOTON_TYP),INTENT(INOUT)                    :: photon
    !--------------------------------------------------------------------------!
    REAL(kind=r2)                                    :: rndx
    !--------------------------------------------------------------------------!
    ! ---
    ! [A] determine starting point & direction
    ! 1. point of emission:
    ! current position

    ! 2. direction of emission
    CALL RAN2(rand_nr,rndx)
    photon%SINPHI = sin( rndx * basics%PIx2 )
    photon%COSPHI = cos( rndx * basics%PIx2 )
    photon%SIN2PH = 2.0_r2 * photon%SINPHI * photon%COSPHI
    photon%COS2PH = 1.0_r2 - 2.0_r2 * photon%SINPHI**2
    
    CALL RAN2(rand_nr,rndx)
    photon%SINTHE = sqrt(4.0_r2 * (rndx - rndx**2))
    photon%COSTHE = 1.0_r2 - 2.0_r2*rndx

    ! ---
    ! [B] initialize new photon parameters
    ! stokes vector
    photon%stokes(:) = fluxes%stokes_ini(:)

    ! photon inside model space
    photon%inside    = .true.

    ! initialize rotation matrix
    photon%D(:,:)    = basics%mat_ident3(:,:)
    
    ! ---
    ! [C] apply rotation matrix
    CALL vecmat(photon)

  END SUBROUTINE start_grain



  ! ################################################################################################
  ! Scattering causes a change of stokes vector (I,Q,U,V).  This transformation occurs in two steps. 
  ! At first, stokes vector is transformed into the new r,l-plane after the PHI-rotation. 
  ! Therefore values of SIN2PH, COS2PH are necessary. 
  ! The second step contains the transformation of the stokes vector with the according 
  ! (index is pointed by LOT) scattering matrix (S11,S12,S33,S34).
  ! By the help of  ALBEDO  the absorption is modeled, if this process takes place 
  ! at the same point ("dust grain")  as the scattering.
  ! ---
!~   subroutine trafo( i_dust_action )
!~     use datatype
!~     use var_global
!~ 
!~     integer, intent(in) :: i_dust_action
!~     real(kind=r2)       :: QHELP, i_1
!~     ! ---
!~     ! 1. absorption
!~     stokes(:) = stokes(:) * albedo( i_dust_action, nr_lam )
!~     i_1       = stokes(1)
!~ 
!~     ! 2. scattering
!~     ! Mathematical positive rotation of r,l-plane around p-axis by the angle  PHI
!~     QHELP     =  COS2PH * stokes(2)  -  SIN2PH * stokes(3)
!~     stokes(3) =  SIN2PH * stokes(2)  +  COS2PH * stokes(3)
!~     stokes(2) =  QHELP
!~     
!~     ! Mathematical negative rotation around r-axis by THETA transformes
!~     ! the Stokes vector (I,Q,U,V) with the scattering matrix.
!~     ! LOT  points to the scattering matrix elements of this angle which
!~     ! was determined by throwing dice before.
!~     stokes(:) = matmul( real(SME( :,:, i_dust_action, nr_lam, lot_th ),kind=r2), stokes(:) )
!~ 
!~     ! Calculation of an adaption factor NORMFA for the stokes vector so
!~     ! that intensity after transformation INORM = IBEFORE * ALBEDO.
!~     stokes(:) = stokes(:) * i_1/stokes(1)
!~ 
!~   end subroutine trafo

end module start_mod

