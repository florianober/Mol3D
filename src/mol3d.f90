! ---
! Mol3D: (SELF-CONSISTENT) N-LTE LINE AND DUST CONTINUUM RADIATIVE TRANSFER CODE
!
! MC3D: RADIATIVE TRANSFER CODE (Mol3D consists of some parts)
!       [by S. Wolf - wolf@astrophysik.uni-kiel.de]
!
! Feedback to:    Florian Ober - fober@astrophysik.uni-kiel.de 
!
!
! ---
! Publication policy:
! If you make use of this code or of parts of it,
! please refer to the following paper:
! Wolf, S., 2003, Comp. Phys. Comm., 150, 99
!         (see http://adsabs.harvard.edu/abs/2003CoPhC.150...99W)
! ---
!
! ---
PROGRAM Mol3D
    USE var_global
    USE datatype

    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    USE source_type

    USE initiate, ONLY       : inimol, cleanup
    USE fileio, ONLY         : vis_plane, save_ch_map
    USE grd_mod, ONLY        : make_grid
    USE simulation_mod, ONLY : run_simu
    USE math_mod

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP)        :: basics
    TYPE(Fluxes_TYP)       :: fluxes
    TYPE(Grid_TYP)         :: grid
    TYPE(Model_TYP)        :: model
    TYPE(Dust_TYP)         :: dust
    TYPE(Gas_TYP)          :: gas
    TYPE(SOURCES)          :: sources_in
    REAL                   :: t0,t1, second
    INTEGER                :: hour, minute

    !$ double precision omp_get_wtime
    !$ double precision wt0,wt1   

    !--------------------------------------------------------------------------!
    ! ---
    CALL cpu_time(t0)
    !$ wt0=omp_get_wtime()

    PRINT *, "======================================================"
    PRINT *, "Mol3D (version "//mol3d_version//")                      Welcome!"
    PRINT *, "======================================================"

    !--------------------------------------------------------------------------!
    ! 1. initiate code 
    !--------------------------------------------------------------------------!
        CALL inimol(basics,fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 2. all initialisation done, now set up the grid!
    !--------------------------------------------------------------------------!
        CALL make_grid(basics, grid, model, dust, gas)
    !--------------------------------------------------------------------------!
    ! 3. run simulation
    !--------------------------------------------------------------------------!
        CALL run_simu(basics, fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 4. clean up
    !--------------------------------------------------------------------------!
        CALL cleanup(basics, fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 5. time measurement
    !--------------------------------------------------------------------------!
        CALL cpu_time(t1)
        hour   = int((t1-t0)/REAL(3600))
        minute = int((t1-t0)/REAL(60)) - hour*60
        second = t1-t0 - hour*3600 - minute*60
        !$ wt1=omp_get_wtime()
        PRINT *,'|'
        PRINT '(3a,i3,a,i2,a,f5.2,a,f10.2,a)', " simulation ",                 &
                    Getproname(basics), " took",                               &
                    hour, 'h ', minute, 'm ', second, 's    (', t1-t0 ,'sec)'
        !$ write (*,'(a,1pg12.4)') 'omp_get_wtime:', wt1-wt0
    !--------------------------------------------------------------------------!

    PRINT *,''
    PRINT *, 'Bye bye'
    PRINT *

END PROGRAM Mol3D
