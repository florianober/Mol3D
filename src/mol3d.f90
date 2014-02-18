! ---
! MOL3D: (SELF-CONSISTENT) RADIATIVE LINE TRANSFER CODE BASED ON THE MONTE CARLO METHOD 
!
! MC3D: RADIATIVE TRANSFER CODE mol3d is based on [by S. Wolf - wolf@astrophysik.uni-kiel.de]
!
! Feedback to:    Florian Ober - fober@astrophysik.uni-kiel.de 
!
!
! ---
! Publication policy:
! If you make use of this code or of parts of it, please refer to the following paper:
! Wolf, S., 2003, Comp. Phys. Comm., 150, 99
!                               (see http://adsabs.harvard.edu/abs/2003CoPhC.150...99W)
! ---
!
! ---
program mol3d
    USE var_globalnew
    USE datatype
    
    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    
    USE initiate, ONLY       : inimol, cleanup
    USE fileio, ONLY         : vis_plane, save_ch_map!, sv_stokes
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
    REAL                   :: t0,t1
    
    !$ double precision omp_get_wtime
    !$ double precision wt0,wt1   
    
!~     INTEGER                :: u,l, i
!~     REAL(kind=r2), DIMENSION(5,5)          :: A
!~     REAL(kind=r2), DIMENSION(5)            :: c,x
    
    !--------------------------------------------------------------------------! 
    ! ---
    call cpu_time(t0)
    !$ wt0=omp_get_wtime()
    
    print *, "==========================================================================="
    print *, "mol3d (version "//mol3d_version//")                                           Welcome!"
    print *, "==========================================================================="
    
    !--------------------------------------------------------------------------!
    ! 1. initiate code [level 1: model setup]
    !--------------------------------------------------------------------------!    
        CALL inimol(basics,fluxes, grid, model, dust , gas)
    !--------------------------------------------------------------------------!
    ! 2. all initialisation done, now set up the grid!
    !--------------------------------------------------------------------------! 
        CALL make_grid(basics, grid, model, dust, gas)
    !--------------------------------------------------------------------------!
    ! 3. run simulation
    !--------------------------------------------------------------------------!
        CALL run_simu(basics, fluxes, grid, model, dust, gas)
    
    !--------------------------------------------------------------------------!
    ! 4. provide some extra visualisations
    !--------------------------------------------------------------------------!
        CALL vis_plane(grid, basics, 1,'densH2',801)
!~         CALL vis_plane(grid, basics, 2,'densH2',801)
        CALL vis_plane(grid, basics, 1,'temp',801)
        CALL vis_plane(grid, basics, 1,'velo',801)
!~         CALL vis_plane(grid, basics, 2,'temp',801)
!~         CALL vis_plane(grid, basics, 2,'velo',801)
    
    !--------------------------------------------------------------------------!
    ! 5. clean up
    !--------------------------------------------------------------------------!
        CALL cleanup(basics, fluxes, grid, model, dust, gas)
    !--------------------------------------------------------------------------!
    ! 6. time measurement
    !--------------------------------------------------------------------------!
    
        CALL cpu_time(t1)
        !$ wt1=omp_get_wtime()
        write (*,'(a,1pg12.4)')    'cpu_time:     ', t1-t0
        !$ write (*,'(a,1pg12.4)') 'omp_get_wtime:', wt1-wt0
        
    !--------------------------------------------------------------------------!
    print *,''
    print *,"Simulation <", Getproname(basics), "> finished."
    print *
end program mol3d
