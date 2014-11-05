MODULE initiate

    USE datatype
    USE var_globalnew

    USE basic_type, ONLY   : Basic_TYP, InitBasic, CloseBasic
    USE randgen_type, ONLY : Randgen_TYP, InitRandgen, CloseRandgen
    USE fluxes_type, ONLY  : Fluxes_TYP, InitFluxes, CloseFluxes
    USE model_type, ONLY   : Model_TYP, InitModel, CloseModel
    USE grid_type, ONLY    : Grid_TYP, InitGrid, CloseGrid
    USE dust_type, ONLY    : Dust_TYP, InitDust, CloseDust
    USE gas_type, ONLY     : Gas_TYP, InitGas, CloseGas
    USE source_type
    
    USE parser_mod
    USE fileio, ONLY       : save_input
    
    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    PRIVATE
    !--------------------------------------------------------------------------!
    PUBLIC ::   inimol, &
                cleanup
    !--------------------------------------------------------------------------!
  
CONTAINS

SUBROUTINE inimol(basics, fluxes, grid, model, dust, gas, sources_in)

    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP)                                  :: basics
    TYPE(Fluxes_TYP)                                 :: fluxes
    TYPE(Grid_TYP)                                   :: grid
    TYPE(Model_TYP)                                  :: model
    TYPE(Dust_TYP)                                   :: dust
    TYPE(Gas_TYP)                                    :: gas
    TYPE(SOURCES)                                    :: sources_in
    !--------------------------------------------------------------------------!
    
    CHARACTER(len=256)                               :: proname
    CHARACTER(len=256)                               :: old_proname
    CHARACTER(len=256)                               :: new_input_file
    CHARACTER(len=256)                               :: input_file
    CHARACTER(len=256)                               :: gas_cat_name
    CHARACTER(len=256)                               :: outname
    CHARACTER(len=256)                               :: help
    CHARACTER(len=256)                               :: r_path
    CHARACTER(len=256)                               :: flux_unit
    CHARACTER(len=4)                                 :: ref_u_str
    CHARACTER(len=32)                                :: grid_name
    CHARACTER(len=32)                                :: in_arg
    CHARACTER(len=8), ALLOCATABLE, DIMENSION(:)      :: dust_cat
    REAL(kind=r2)                                    :: ref_u
    REAL(kind=r2)                                    :: r_in
    REAL(kind=r2)                                    :: r_ou
    REAL(kind=r2)                                    :: sf
    REAL(kind=r2)                                    :: mass_dust
    REAL(kind=r2)                                    :: molratio
    REAL(kind=r2)                                    :: abundance
    REAL(kind=r2)                                    :: vel_max
    REAL(kind=r2)                                    :: t_eff
    REAL(kind=r2)                                    :: R_star
    REAL(kind=r2)                                    :: M_star
    REAL(kind=r2)                                    :: distance
    REAL(kind=r2)                                    :: sizexp
    REAL(kind=r1)                                    :: t_dust_max
    REAL(kind=r1)                                    :: t_dust_min
    
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: dens_dust
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: th_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: ph_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: zoom_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: al_map
    
    INTEGER, ALLOCATABLE, DIMENSION(:)               :: tr_cat
    
    INTEGER                                          :: photon_type
    INTEGER                                          :: n_a
    INTEGER                                          :: n_b
    INTEGER                                          :: n_c
    INTEGER                                          :: concept_ps
    INTEGER                                          :: n_map
    INTEGER                                          :: n_dust_emi
    INTEGER(8)                                       :: no_photon
    INTEGER                                          :: n_bin_map
    INTEGER                                          :: n_dust
    INTEGER                                          :: n_tr
    INTEGER                                          :: i_vel_chan
    INTEGER                                          :: mode
    INTEGER                                          :: aniso
    INTEGER                                          :: n_scatt_th
    INTEGER                                          :: n_tem
    INTEGER                                          :: grid_type
    INTEGER                                          :: nrndpt  !! resolution of scattering distribution"
    
    INTEGER                                          :: num_core !number of cores used in parallel part
    INTEGER,DIMENSION(8)                             :: dtime
    
    LOGICAL                                          :: calc_tmp
    LOGICAL                                          :: old_model
    LOGICAL                                          :: pluto_data  !remark: this should be removed
    LOGICAL                                          :: do_raytr
    LOGICAL                                          :: do_velo_ch_map
    LOGICAL                                          :: do_continuum_map
    !--------------------------------------------------------------------------!
    INTENT(INOUT)                                    :: basics, &
                                                        fluxes, &
                                                        grid, & 
                                                        model, &
                                                        dust, &
                                                        gas, &
                                                        sources_in

    print *, "Code initialization"
    IF ( iargc() == 1 ) THEN
        ! use command line given input name
        CALL getarg(1, in_arg)
        new_input_file = in_arg
    ELSE
        new_input_file = 'input/input.dat'
    END IF

    
    CALL parse('r_path',r_path,new_input_file)       !define results path
    CALL parse('proname',proname,new_input_file)     !project name
    

    ! Tell the project name
    PRINT '(2A)'," Simulation name: ", TRIM(proname)
    
    ref_u      = con_AU                                 !AU at the moment for prop. disks!
    ref_u_str  = "[AU]"                                 !ref unit
    concept_ps = 1                                      !Emission concept (1) Point-like, isotropic
    
    ! set further parameters 
    ! these variables are written to global variables, as mc3d used to do this, this should be changed!
    ! mostly done!
    n_tem             = 40000
    n_scatt_th        = 361
    nrndpt            = 10000
!~     acc_select_level  = 1.0e-6 ! old and not in use anymore
!~     i_lam_show        = 10 ! old and not in use anymore
    t_dust_min        = 0
    t_dust_max        = 2000.0
    n_interact_max    = 100000
    
    show_error = .False.                             ! show some minor warnings
    velo_type  = 1                                   ! analytical velocity distribution
!~     velo_type  = 2                                ! lin. interpolated velocity distribution

    ! use results of earlier calculations
    !
    CALL parse('old_model',old_model,new_input_file)
    
    IF (old_model) THEN
        ! use the input file of the old calculation to ensure that there is no change inside the model itself
        CALL parse('old_proname',old_proname,new_input_file)
        IF (old_proname == proname) THEN
            print *, 'ERRORL: The given project name equals the old project name. Please adjust!'
            STOP
        END IF
        
        input_file = TRIM(r_path)//TRIM(old_proname)//'_input_file.dat'
    ELSE
        ! use the parameters given in the old input_file as model parameters
        input_file = new_input_file
        old_proname = 'notdefined'
    END IF
        
    ! save input to file
    outname = TRIM(r_path)//TRIM(proname)//'_input_file.dat'
    
    open(unit=3, file=TRIM(outname), &
            action="write", status="unknown", form="formatted")
    
    call date_and_time(values=dtime)
    WRITE(help,fmt='(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') dtime(1),'-', dtime(2),'-',dtime(3),' / ', dtime(5),':', dtime(6)
    
    WRITE(unit=3,fmt='(A)')  '##################################################################################'
    WRITE(unit=3,fmt='(A)')  '#                            mol3d input file                                    #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# All parameters that can be changed without recompiling mol3d are listed here.  #'
    WRITE(unit=3,fmt='(A)')  '# Please be aware that some parameters not intendet to be changeable yet.        #'
    WRITE(unit=3,fmt='(A)')  '# Therefore, they are note listed here, but accessable via the                   #'
    WRITE(unit=3,fmt='(A)')  '# "src/initiate.f90" source file.                                                    #'
    WRITE(unit=3,fmt='(A)')  '# If you have questions, please feel free to ask me.                             #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# author         : Florian Ober                                                  #'
    WRITE(unit=3,fmt='(A)')  '# email          : fober@astrophysik.uni-kiel.de                                 #'
    WRITE(unit=3,fmt='(A)')  '# mol3d version  : '//mol3d_version//'                                                      #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# execution date : '//TRIM(help)//'                                            #'
    WRITE(unit=3,fmt='(A)')  '##################################################################################'
    WRITE(unit=3,fmt='(A)')  ''
    WRITE(unit=3,fmt='(A)')  ''
    

    WRITE(unit=3,fmt='(A)') 'proname = {'//TRIM(proname)//'}                 project name'
    WRITE(unit=3,fmt='(A)') 'r_path = {'//TRIM(r_path)//'}                  define results path'
    
    WRITE(help,fmt='(L1)') old_model
    WRITE(unit=3,fmt='(A)') 'old_model = {'//TRIM(help)//'}'
    
    WRITE(unit=3,fmt='(A)') 'old_proname = {'//TRIM(old_proname)//&
        '}           the project name of the older calculation'
    
    CALL parse('num_core',num_core,new_input_file)  ! number of cores used (raytracing, temperature calculation (TbD))
    
    WRITE(help,fmt='(I3.3)') num_core
    !$ PRINT '(A,I3,A)', " Parallel mode, using: ", num_core, " cores"
    WRITE(unit=3,fmt='(A)') 'num_core = {'//TRIM(help)// &
    '}                     number of cores used (raytracing, temperature calculation (TbD))'
    WRITE(unit=3,fmt='(A)') '' 
    
    ! calculate dust temperature (with monte carlo method), or use an analytical expression
    ! only valid if old_model = False
    CALL parse('calc_temp',calc_tmp,new_input_file) 
    
    WRITE(help,fmt='(L1)') calc_tmp
    WRITE(unit=3,fmt='(A)') 'calc_temp = {'//TRIM(help)//&
    '}                      calculate temperature (with monte carlo method)'

!~     CALL parse('do_raytr',do_raytr,new_input_file)
!~     WRITE(help,fmt='(L1)') do_raytr
!~     WRITE(unit=3,fmt='(A)') 'do_raytr = {'//TRIM(help)//'}'
!~     WRITE(unit=3,fmt='(A)') '' 
    
!~     pluto_data = .True.
    pluto_data = .False.
    do_continuum_map = .True.
    do_velo_ch_map   = .False.
    IF (do_continuum_map .or. do_velo_ch_map) THEN
        do_raytr = .True.
    ELSE
        do_raytr = .False.
    END IF
    
    CALL InitBasic(basics,photon_type,'Reemission map',proname,r_path, concept_ps, calc_tmp, old_proname, &
                    old_model, &
                    do_raytr,do_continuum_map,do_velo_ch_map, n_tem, t_dust_min, t_dust_max, &
                    num_core, pluto_data)
    !--------------------------------------------------------------------------! 
     
    
    CALL parse('r_in',r_in,input_file)
    WRITE(help,fmt='(ES15.6)') r_in
    WRITE(unit=3,fmt='(A)') 'r_in = {'//TRIM(help)//&
    '}             inner radius'

    CALL parse('r_ou',r_ou,input_file)
    WRITE(help,fmt='(ES15.6)') r_ou    
    WRITE(unit=3,fmt='(A)') 'r_ou = {'//TRIM(help)//&
    '}             outer radius'
    
    
    CALL parse('mass_dust',mass_dust,input_file)
    WRITE(help,fmt='(ES15.6)') mass_dust  
    WRITE(unit=3,fmt='(A)') 'mass_dust = {'//TRIM(help)//&
    '}        dust mass'
    
    CALL parse('t_eff',t_eff,input_file)
    WRITE(help,fmt='(ES15.6)') t_eff
    WRITE(unit=3,fmt='(A)') 't_eff = {'//TRIM(help)//&
    '}            effective Temperature'
    
    CALL parse('R_star',R_star,input_file)
    WRITE(help,fmt='(ES15.6)') R_star
    WRITE(unit=3,fmt='(A)') 'R_star = {'//TRIM(help)//&
    '}           radius of star'
    
    CALL parse('M_star',M_star,new_input_file)
    WRITE(help,fmt='(ES15.6)') M_star
    WRITE(unit=3,fmt='(A)') 'M_star = {'//TRIM(help)//&
    '}           stellar mass'
    
    WRITE(unit=3,fmt='(A)') ''
    n_map     = 1                          !numper of maps to calculate
    CALL parse('distance',distance,new_input_file)
    WRITE(help,fmt='(ES15.6)') distance
    WRITE(unit=3,fmt='(A)') 'distance = {'//TRIM(help)//&
    '}         object distance [pc]'
    
    CALL parse('n_bin_map',n_bin_map,new_input_file)
    WRITE(help,fmt='(I4.4)') n_bin_map
    WRITE(unit=3,fmt='(A)') 'n_bin_map = {'//TRIM(help)// &
    '}                   half number of pixel, total pixel = 2 * n_bin_map1 + 1'
    
    
    CALL parse('no_photon',no_photon,input_file)
    WRITE(help,fmt='(ES9.2)') real(no_photon,kind=r2)
    WRITE(unit=3,fmt='(A)') 'no_photon = {'//TRIM(help)// &
    '}            number of photons per wavelength'
    
    ALLOCATE(th_map(1:n_map))
    ALLOCATE(ph_map(1:n_map))
    ALLOCATE(zoom_map(1:n_map))
    ALLOCATE(al_map(1:n_map))
    
    CALL parse('th_map',th_map(1),new_input_file)
    WRITE(help,fmt='(F6.2)') th_map
    WRITE(unit=3,fmt='(A)') 'th_map = {'//TRIM(help)//&
    '}                    inclination to observer'
    
    CALL parse('ph_map',ph_map(1),new_input_file) 
    WRITE(help,fmt='(F6.2)') ph_map
    WRITE(unit=3,fmt='(A)') 'ph_map = {'//TRIM(help)//&
    '}                   '
    WRITE(unit=3,fmt='(A)') ''
    
!~    zoom_map(1) = 1.0_r2
    
    CALL parse('zoom_map',zoom_map(1),new_input_file) 
    WRITE(help,fmt='(F6.2)') zoom_map(1)
    WRITE(unit=3,fmt='(A)') 'zoom_map = {'//TRIM(help)//&
    '}                  factor to zoom inside the model >= 1 !'
    WRITE(unit=3,fmt='(A)') ''
    
    
    al_map(1)   = 1
    
    CALL InitModel(model,1, ref_u_str, ref_u, r_in, r_ou, mass_dust, t_eff, R_star, M_star, n_map, distance, &
                   no_photon,th_map,ph_map,zoom_map,al_map,n_bin_map)
                    
    DEALLOCATE( th_map, &
                ph_map, &
                zoom_map, &
                al_map)
    
    !--------------------------------------------------------------------------!
    
    
    n_tr = 1                 !Number of transitions, please do not change this value. Only calculate
                             ! one transition at one time
    ALLOCATE( tr_cat(1:n_tr) )
    
    CALL parse('molratio',molratio,new_input_file)
    WRITE(help,fmt='(ES15.6)') molratio
    WRITE(unit=3,fmt='(A)') 'molratio = {'//TRIM(help)//&
    '}         gas-dust ratio'
    
    CALL parse('abundance',abundance,new_input_file)
    WRITE(help,fmt='(ES15.6)') abundance
    WRITE(unit=3,fmt='(A)') 'abundance = {'//TRIM(help)//&
    '}        selected molecule abundance to H'

    !set the method to calculate the level populations
    ! 1 := LTE
    ! 2 := FEP
    ! 3 := LVG
    ! .
    ! 9 := read user file (TbD)
    CALL parse('line_mode',mode,new_input_file)
    WRITE(help,fmt='(I1)') mode
    WRITE(unit=3,fmt='(A)') 'line_mode = {'//TRIM(help)// &
    '}                      set the method to calculate the level populations,1 := LTE  2 := FEP  3 := LVG  9 := file (TbD)'
    
    !tr_cat(1) = 4            !transition (tr) in molecular cat   !tbd
    CALL parse('line',tr_cat(1),new_input_file)
    WRITE(help,fmt='(I1)') tr_cat(1)
    WRITE(unit=3,fmt='(A)') 'line = {'//TRIM(help)// &
    '}                           line transition no (see gas_input file)  '
    WRITE(unit=3,fmt='(A)') ''
    CALL parse('gas_cat_name',gas_cat_name,new_input_file)
    WRITE(unit=3,fmt='(A)') 'gas_cat_name = {'//TRIM(gas_cat_name)// &
    '}              gas catalog name (in input/mol) (Leiden database style)'
    
    CALL parse('i_vel_chan',i_vel_chan,new_input_file)
    WRITE(help,fmt='(I3)') i_vel_chan
    WRITE(unit=3,fmt='(A)') 'i_vel_chan = {'//TRIM(help)// &
    '}                   (half) Number of velocity channels'
    
    
    CALL parse('vel_max',vel_max,new_input_file)    
    WRITE(help,fmt='(F8.2)') vel_max
    WRITE(unit=3,fmt='(A)') 'vel_max = {'//TRIM(help)//&
    '}                  max velocity in spectrum (in m/s)'
    
    CALL InitGas(gas, mode, gas_cat_name, molratio, abundance, n_tr, tr_cat, i_vel_chan, vel_max)
    
    DEALLOCATE( tr_cat )
    
    !--------------------------------------------------------------------------! 
    
    n_dust = 1                 !Number of dust species"
    ALLOCATE( dens_dust(1:n_dust) )
    
    dens_dust(1) = 2.5         !Density [g/cm^3] of dust grain species 1"
    sizexp = 0.0_r2            !Exponent of the dust grain size distribution  mrn  [-3.5]
    
    aniso  = 1                 !Scattering 1) Anisotropic 2) Isotropic
    ALLOCATE(dust_cat(1:n_dust))
    
    CALL parse('dust_cat_name',dust_cat(1),input_file)
    WRITE(unit=3,fmt='(A)') 'dust_cat_name = {'//TRIM(dust_cat(1))// &
    '}           dust catalog name (in input/dust) (in mc3d style)'
    WRITE(unit=3,fmt='(A)') ''
    !call get_n_den_par()
    !if (n_den_par>0) then
    !    allocate(den_par(1:n_den_par)) !"Number of free parameters in the current model: ", n_den_par
    !    den_par(1) = 2.34
    !end if
    CALL InitDust(dust, basics, gas, model, 1, 'mrn', n_dust, dens_dust, sizexp, aniso, dust_cat,n_scatt_th, nrndpt)
    DEALLOCATE( dens_dust, &
                dust_cat)

    
    !--------------------------------------------------------------------------!
    
    IF (pluto_data .and. .not. old_model ) THEN
        ! fix the cell properties
        grid_name = 'spherical'
        grid_type = 3
        n_a = 100
        n_b = 65
        n_c = 124
        CALL parse('sf',sf,input_file)  ! in fact, it is used to scale the grid
        model%r_in = sf*model%r_in
        model%r_ou = sf*model%r_ou
    ELSE
        ! this is the normal way, read the grid properties from the input file
        CALL parse('grid_name',grid_name,input_file)  
          ! possible values: spherical (mostly mc3d style)
          !                  cylindrical (in progress, working but there are some bugs, test for yourself)
          !                  cartesian   (planed)
          !
          ! EXAMPLE:
          ! 
          ! spherical :  r = / 1 = rho  \
          !                  | 2 = theta|
          !                  \ 3 = phi  /
          ! n(1) = n_r :  number of point in r direction
          ! n(2) = n_th : number of point in th direction, must be odd
          ! n(3) = n_ph : number of point in oh direction, must be 1 or even
          ! sf   = step factor for logarithmic r scale
          ! ..
        CALL parse('grid_type',grid_type,input_file)
        CALL parse('n_a',n_a,input_file)
        CALL parse('n_b',n_b,input_file) 
        CALL parse('n_c',n_c,input_file)
        CALL parse('sf',sf,input_file) 
        
        IF (pluto_data) THEN
            model%r_in = sf*model%r_in
            model%r_ou = sf*model%r_ou
        END IF

    END IF
    
    ! write informations about the grid in the log/input file
    
    WRITE(unit=3,fmt='(A)') 'grid_name = {'//TRIM(grid_name)// &
        '}              possible values: spherical (OK), cylindrical (testing), cartesian (planned)'
        
    WRITE(help,fmt='(I1)') grid_type
    WRITE(unit=3,fmt='(A)') 'grid_type = {'//TRIM(help)// &
    '}                      default = 1. 2 and higher values for user defined grids'
    
    WRITE(help,fmt='(I4)') n_a
    WRITE(unit=3,fmt='(A)') 'n_a = {'//TRIM(help)// &
    '}                     '

    WRITE(help,fmt='(I4)') n_b
    WRITE(unit=3,fmt='(A)') 'n_b = {'//TRIM(help)// &
    '}                     '
    
    WRITE(help,fmt='(I4)') n_c
    WRITE(unit=3,fmt='(A)') 'n_c = {'//TRIM(help)// &
    '}                     ' 

    WRITE(help,fmt='(F12.9)') sf
    WRITE(unit=3,fmt='(A)') 'sf = {'//TRIM(help)//&
    '}                  stepfactor for logarithm grid'
    WRITE(unit=3,fmt='(A)') ''
    
    ! give the type of the specified grid 
    ! 1 == analytical version (prefered)
    ! 2 == user input version, under construction
    !
    
!   Some examples here:
!~     SELECT CASE(grid_type)
    
!~     CASE('cylindrical')
!~         grid_type = 1
!~         n_a  = 100
!~         sf   = 1.00001
!~         n_b  = 70          
!~         n_c  = 101   
!~         
!~     CASE('spherical')
!~         grid_type = 1
!~         n_a  = 300
!~         sf   = 1.01
!~         n_b  = 257
!~         n_c  = 1
    
!~     END SELECT
    CALL InitGrid(grid,grid_type,grid_name, n_a, sf, n_b, n_c,dust%n_dust, & 
                  gas%egy_lvl)
    !--------------------------------------------------------------------------! 
    
    CALL InitSources(sources_in,1,'sources',dust)
    ! -------------------------------------
    ! Mode 1: give R_star and T_star
    !   CALL AddSources(sources_in,1,(/0.0_r2,0.0_r2,0.0_r2/), R_star=model%r_star, T_star=model%t_star)
    ! Mode 2: give T_star and Luminosity
    !   CALL AddSources(sources_in,2,(/0.0_r2,0.0_r2,0.0_r2/), T_star=model%T_star, L_star=REAL(L_sun*1.907343,kind=r1))
    ! -------------------------------------
    
    ! 1 source (used as the primary source, defined in the input file)
    CALL AddSources(sources_in,1,(/0.0_r2,0.0_r2,0.0_r2/), R_star=model%r_star, T_star=model%t_star)
    ! 2 source (added by hand, could be a embedded planet or whatever)
!~     CALL AddSources(sources_in,2,(/sf*4.95_r2,-sf*0.72_r2,0.0_r2/), T_star=1000.0, L_star=REAL(1e-4*L_sun,kind=r1))
    print *, sources_in%L_total
    print *, sources_in%L_total/(4*PI)
    ! .. source..
    
    print '(A,I3)', " Sources found: ", sources_in%n_sources
    print '(A,F5.2,A)'," Total Luminosity included: ", sources_in%L_total/L_sun, " L_sun"
    !--------------------------------------------------------------------------! 
    ! Stokes vector: Unpolarized radiation (assumed initial state)
    n_dust_emi = 0
    CALL parse('flux_unit',flux_unit,new_input_file)
    WRITE(unit=3,fmt='(A)') 'flux_unit = {'//TRIM(flux_unit)// &
    '}                 possible values: Jy_pix, T_mb'
    
    CALL InitFluxes(fluxes,1,flux_unit, n_dust_emi, &
                    model%n_bin_map,gas%n_tr, gas%i_vel_chan,dust%n_lam)
    
    
    !--------------------------------------------------------------------------! 
    !  Save input file to project results
    ! 

!~     CALL save_input(basics,input_file)
    print *, "Code initialization done !"
    Close(unit=3)
    
end subroutine inimol

SUBROUTINE cleanup(basics ,fluxes ,grid , model, dust, gas, sources_in)

    IMPLICIT NONE
    
    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP)                                  :: basics
    TYPE(Fluxes_TYP)                                 :: fluxes
    TYPE(Grid_TYP)                                   :: grid
    TYPE(Model_TYP)                                  :: model
    TYPE(Dust_TYP)                                   :: dust
    TYPE(gas_TYP)                                    :: gas
    TYPE(SOURCES)                                    :: sources_in
    !--------------------------------------------------------------------------!
    INTENT(INOUT)                                    :: basics, &
                                                        fluxes, &
                                                        grid, & 
                                                        model, &
                                                        dust, &
                                                        gas, &
                                                        sources_in
                                                    
    CALL CloseBasic(basics)
    CALL CloseDust(dust)
    CALL CloseFluxes(fluxes)
    CALL CloseGrid(grid)
    CALL CloseModel(model)
    CALL CloseGas(gas)
    CALL CloseSources(sources_in)
    

END SUBROUTINE cleanup

END MODULE initiate
