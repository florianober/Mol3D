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
MODULE initiate

    USE datatype
    USE var_global

    USE basic_type, ONLY   : Basic_TYP, InitBasic, CloseBasic
    USE fluxes_type, ONLY  : Fluxes_TYP, InitFluxes, CloseFluxes
    USE model_type, ONLY   : Model_TYP, InitModel, CloseModel
    USE grid_type, ONLY    : Grid_TYP, InitGrid, CloseGrid
    USE dust_type, ONLY    : Dust_TYP, InitDust, CloseDust
    USE gas_type, ONLY     : Gas_TYP, InitGas, CloseGas
    USE source_type

    USE parser_mod
    USE fileio, ONLY       : read_no_cells

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
    CHARACTER(len=252)                               :: file_a, file_b, file_c
    CHARACTER(len=4)                                 :: ref_u_str
    CHARACTER(len=20)                                :: randgen_name
    CHARACTER(len=32)                                :: grid_name
    CHARACTER(len=256)                               :: in_arg
    CHARACTER(len=8), ALLOCATABLE, DIMENSION(:)      :: dust_cat
    REAL(kind=r2)                                    :: ref_u
    REAL(kind=r2)                                    :: r_in
    REAL(kind=r2)                                    :: r_ou
    REAL(kind=r2)                                    :: sf
    REAL(kind=r2)                                    :: mass_dust
    REAL(kind=r2)                                    :: molratio
    REAL(kind=r2)                                    :: abundance
    REAL(kind=r2)                                    :: vel_max
    REAL(kind=r2)                                    :: T_star
    REAL(kind=r2)                                    :: R_star
    REAL(kind=r2)                                    :: M_star
    REAL(kind=r2)                                    :: distance
    REAL(kind=r2)                                    :: sizexp
    REAL(kind=r1)                                    :: t_dust_max
    REAL(kind=r1)                                    :: t_dust_min
    REAL(kind=r2)                                    :: v_turb

    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: dens_dust
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: plist
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: th_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: ph_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: zoom_map
    REAL(kind=r2), ALLOCATABLE, DIMENSION(:)         :: al_map

    INTEGER, ALLOCATABLE, DIMENSION(:)               :: tr_cat

    INTEGER                                          :: photon_type
    INTEGER                                          :: n_a
    INTEGER                                          :: n_b
    INTEGER                                          :: n_c
    INTEGER                                          :: n_map
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
    INTEGER                                          :: io_error

    ! resolution of scattering distribution"
    INTEGER                                          :: nrndpt
    ! number of cores used in parallel mode
    INTEGER                                          :: num_core
    INTEGER,DIMENSION(8)                             :: dtime

    LOGICAL                                          :: do_MC_temperature
    LOGICAL                                          :: old_model
    LOGICAL                                          :: do_velo_ch_map
    LOGICAL                                          :: do_continuum_raytrace
    LOGICAL                                          :: do_continuum_mc
    LOGICAL                                          :: peel_off
    LOGICAL                                          :: do_save_T_exc
    !--------------------------------------------------------------------------!
    INTENT(INOUT)                                    :: basics, &
                                                        fluxes, &
                                                        grid, &
                                                        model, &
                                                        dust, &
                                                        gas, &
                                                        sources_in
    !--------------------------------------------------------------------------!
    print *, "Code initialization"

    !----------------------------  Init Basics  -------------------------------!

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

    ref_u      = con_AU           ! AU at the moment for  protoplanetary disks

    ref_u_str  = "AU"             ! reference unit

    ! set further parameters


    n_scatt_th        = 361
    nrndpt            = 10000
    n_tem             = 500
    t_dust_min        = 2.73
    t_dust_max        = 3000.0
    n_interact_max    = 100000

    show_error = .False.               ! show some minor warnings

    ! use results of earlier calculations
    !
    CALL parse('old_model', old_model, new_input_file)

    IF (old_model) THEN
        ! use the input file of the old calculation to ensure model consistency
        CALL parse('old_proname',old_proname,new_input_file)
        IF (old_proname == proname) THEN
            print *, 'ERROR: The given project name equals the old project name. Please adjust!'
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

    CALL date_and_time(values=dtime)
    WRITE(help,fmt='(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') dtime(1),'-', dtime(2),'-',dtime(3),' / ', dtime(5),':', dtime(6)

    WRITE(unit=3,fmt='(A)')  '##################################################################################'
    WRITE(unit=3,fmt='(A)')  '#                            Mol3D input file                                    #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# All parameters that can be changed without recompiling Mol3D are listed here.  #'
    WRITE(unit=3,fmt='(A)')  '# Please be aware that some parameters are not intended to be changeable yet.    #'
    WRITE(unit=3,fmt='(A)')  '# They are not listed here, but accessable via the                               #'
    WRITE(unit=3,fmt='(A)')  '# "src/initiate.f90" source file.                                                #'
    WRITE(unit=3,fmt='(A)')  '# If you have questions, please feel free to ask me.                             #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# author         : Florian Ober                                                  #'
    WRITE(unit=3,fmt='(A)')  '# email          : fober@astrophysik.uni-kiel.de                                 #'
    WRITE(unit=3,fmt='(A)')  '# git commit     : '//TRIM(mol3d_version)//'                      #'
    WRITE(unit=3,fmt='(A)')  '#                                                                                #'
    WRITE(unit=3,fmt='(A)')  '# execution date : '//TRIM(help)//'                                            #'
    WRITE(unit=3,fmt='(A)')  '##################################################################################'
    WRITE(unit=3,fmt='(A)')  ''
    WRITE(unit=3,fmt='(A)')  ''


    WRITE(unit=3,fmt='(A)') 'proname = {'//TRIM(proname)//'}                  project name'
    WRITE(unit=3,fmt='(A)') 'r_path = {'//TRIM(r_path)//'}                  define results path'

    WRITE(help,fmt='(L1)') old_model
    WRITE(unit=3,fmt='(A)') 'old_model = {'//TRIM(help)//'}'

    WRITE(unit=3,fmt='(A)') 'old_proname = {'//TRIM(old_proname)//&
        '}           the project name of the older calculation'


    CALL GETENV('OMP_NUM_THREADS', help)

    READ(help, '(i10)',iostat=io_error) num_core
    IF (io_error /= 0 .or. LEN(TRIM(help)) == 0) THEN
        ! read from input file
        CALL parse('num_core', num_core, new_input_file)
    END IF

    IF (num_core < 1) num_core = 1 ! just to prevent negativ or zero cores

    WRITE(help,fmt='(I3.3)') num_core
    !$ PRINT '(A,I3,A)', " Parallel mode, using: ", num_core, " cores"
    WRITE(unit=3,fmt='(A)') 'num_core = {'//TRIM(help)// &
    '}                     number of cores used '
    WRITE(unit=3,fmt='(A)') ''

    ! calculate dust temperature (with monte carlo method), or use an analytical
    ! expression only valid if old_model = False
    IF (.not. old_model) THEN
        CALL parse('do_MC_temperature',do_MC_temperature,new_input_file)

        WRITE(help,fmt='(L1)') do_MC_temperature
        WRITE(unit=3,fmt='(A)') 'do_MC_temperature = {'//TRIM(help)//&
        '}              calculate temperature (with monte carlo method)'
        CALL parse('do_MC_temperature',do_MC_temperature,new_input_file)

    ELSE
        do_MC_temperature = .False.
    END IF

    CALL parse('do_peel_off',peel_off,new_input_file)
    WRITE(help,fmt='(L1)') peel_off
    WRITE(unit=3,fmt='(A)') 'do_peel_off = {'//TRIM(help)//&
    '}                    use peel-off technique'

    CALL parse('do_continuum_mono', do_continuum_mc, new_input_file)
    WRITE(help,fmt='(L1)') do_continuum_mc
    WRITE(unit=3,fmt='(A)') 'do_continuum_mono = {'//TRIM(help)//&
    '}              make full Monte Carlo continuum maps and seds'


    CALL parse('do_continuum_raytrace', do_continuum_raytrace, new_input_file)
    WRITE(help,fmt='(L1)') do_continuum_raytrace
    WRITE(unit=3,fmt='(A)') 'do_continuum_raytrace = {'//TRIM(help)//'}'
    CALL parse('do_velo_ch_map', do_velo_ch_map, new_input_file)
    WRITE(help,fmt='(L1)') do_velo_ch_map
    WRITE(unit=3,fmt='(A)') 'do_velo_ch_map = {'//TRIM(help)//'}'
    !
    !  randgen_name = name of the desired random number generator,
    !                 possible values at the moment:
    !            XX   'MT'       : Mersenne Twister (not available at the moment)
    !                              Copyright (C) 1997 Makoto Matsumoto and
    !                                                 Takuji Nishimura
    !                 'TAUS88'   : L'Ecuyer, P. (1996)
    !                 'LFSR113'  : L'Ecuyer, P. (1999)
    !                 'KISS99'   : George Marsaglia 1999 (prefered)
    !                 'CONG'     : very basic linear congruential generator
    !                 'COMPILER' : The compiler own
    randgen_name = 'KISS99'
    ! save the calculated excitation temperature?
    do_save_T_exc = .False.
    CALL InitBasic(basics,photon_type,'Reemission map', proname, r_path,       &
                   do_MC_temperature, old_proname,                             &
                   old_model, do_continuum_raytrace,                           &
                   do_continuum_mc, do_velo_ch_map, peel_off, do_save_T_exc,   &
                   n_tem, t_dust_min, t_dust_max, num_core,                    &
                   input_file, new_input_file, randgen_name)

    !------------------------------  Init Model  ------------------------------!

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

    CALL parse('T_star', T_star,input_file)
    WRITE(help,fmt='(ES15.6)') T_star
    WRITE(unit=3,fmt='(A)') 'T_star = {'//TRIM(help)//&
    '}             effective star Temperature'

    CALL parse('R_star', R_star,input_file)
    WRITE(help,fmt='(ES15.6)') R_star
    WRITE(unit=3,fmt='(A)') 'R_star = {'//TRIM(help)//&
    '}           radius of star'

    CALL parse('M_star', M_star,input_file)
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


    CALL parse('no_photon',no_photon,new_input_file)
    WRITE(help,fmt='(ES9.2)') real(no_photon,kind=r2)
    WRITE(unit=3,fmt='(A)') 'no_photon = {'//TRIM(help)// &
    '}              number of photons per wavelength'

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


    CALL parse('zoom_map',zoom_map(1),new_input_file)
    WRITE(help,fmt='(F6.2)') zoom_map(1)
    WRITE(unit=3,fmt='(A)') 'zoom_map = {'//TRIM(help)//&
    '}                  factor to zoom inside the model >= 1 !'
    WRITE(unit=3,fmt='(A)') ''
    
    
    !acceptance angle (MC3D style...), max 5 deg
    al_map(1)   = MIN(5.0_r2, MAX(th_map(1), 1.0_r2))
    
    ! additional model parameter
    ! you can define additional model parameter via "plist" array
    ! example:
    !ALLOCATE(plist(1:1))
    !CALL parse('keyword', plist(1), new_input_file)
    
    CALL InitModel(model,1, ref_u_str, ref_u, r_in, r_ou, mass_dust, T_star,    &
                   R_star, M_star, n_map, distance, no_photon, th_map,ph_map,   &
                   zoom_map, al_map, n_bin_map, plist)
    
    DEALLOCATE( th_map,   &
                ph_map,   &
                zoom_map, &
                al_map)

    !-----------------------------  Init Gas  ---------------------------------!

    n_tr = 1    ! Number of transitions, please do not change this value.
                ! Only calculate one transition at once, for the moment
    ALLOCATE( tr_cat(1:n_tr) )

    CALL parse('molratio',molratio,new_input_file)
    WRITE(help,fmt='(ES15.6)') molratio
    WRITE(unit=3,fmt='(A)') 'molratio = {'//TRIM(help)//&
    '}         gas-dust ratio'

    CALL parse('abundance',abundance,new_input_file)
    WRITE(help,fmt='(ES15.6)') abundance
    WRITE(unit=3,fmt='(A)') 'abundance = {'//TRIM(help)//&
    '}        selected molecule number relative abundance to H2'

    !set the method to calculate the level populations
    ! 1 := LTE
    ! 2 := FEP
    ! 3 := LVG
    ! .
    ! 9 := read user file (TbD)
    CALL parse('line_mode',mode,new_input_file)
    WRITE(help,fmt='(I1)') mode
    WRITE(unit=3,fmt='(A)') 'line_mode = {'//TRIM(help)// &
    '}                      set the method to calculate the level populations,'
    WRITE(unit=3,fmt='(A)') '              &
    &                       1 := LTE  2 := FEP  3 := LVG'

    !tr_cat(1) = 4            !transition (tr) in molecular cat   !tbd
    CALL parse('line', tr_cat(1), new_input_file)
    WRITE(help,fmt='(I2)') tr_cat(1)
    WRITE(unit=3,fmt='(A)') 'line = {'//TRIM(help)// &
    '}                          line transition no (see gas_input file)  '
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

    CALL parse('v_turb', v_turb, new_input_file)
    WRITE(help,fmt='(F8.2)') v_turb
    WRITE(unit=3,fmt='(A)') 'v_turb = {'//TRIM(help)//&
    '}                    turbulent velocity (in m/s)'
    CALL InitGas(gas, mode, gas_cat_name, molratio, abundance, n_tr, tr_cat,   &
                 i_vel_chan, vel_max, v_turb)

    DEALLOCATE( tr_cat )

    !----------------------------  Init Dust  ---------------------------------!
    ! TbD: similar to the source type, generate a function to add dust
    !      -> no dust (n_dust=0) can be assumed (no temp calculation then)
    n_dust = 1                 ! Number of dust species"
    ALLOCATE( dens_dust(1:n_dust) )

    dens_dust(:) = 2.5      ! Density [g/cm^3] of dust grain species 1"
    sizexp = -3.0           ! Exponent of the dust grain size distribution
                            ! Please ALWAYS check the resulting mass of
                            ! each dust species to ensure the correct usage

    aniso = 1               ! Scattering 1) Anisotropic (Mie, prefered)
                            !            2) Isotropic
                            !            3) Henyey-Greenstein (not much tested)

    ALLOCATE(dust_cat(1:n_dust))

    CALL parse('dust_cat_name',dust_cat(1), input_file)
    WRITE(unit=3,fmt='(A)') 'dust_cat_name = {'//TRIM(dust_cat(1))// &
    '}           dust catalog name (in input/dust) (in mc3d style)'
    WRITE(unit=3,fmt='(A)') ''

    !dust_cat(2) = 'mrn---v4'

    CALL InitDust(dust, basics, gas, 1, 'mrn', n_dust, dens_dust,       &
                  sizexp, aniso, dust_cat, n_scatt_th, nrndpt)
    DEALLOCATE( dens_dust, &
                dust_cat)
    print '(A,I3)', " Number of dust components included  : ", dust%n_dust

    !----------------------------  Init Sources  ------------------------------!

    CALL InitSources(sources_in, 1, 'sources', dust)

    ! -------------------------------------
    ! Mode 1: give R_star and T_star, e.g.,
    !   CALL AddSources(sources_in,1 ,(/0.0_r2,0.0_r2,0.0_r2/),                &
    !                   T_star=model%t_star, R_star=model%r_star)
    ! Mode 2: give T_star and Luminosity, e.g.,
    !   CALL AddSources(sources_in,2,(/0.0_r2,0.0_r2,0.0_r2/),                 &
    !                   T_star=model%T_star, L_star=REAL(L_sun*1.9,kind=r1))
    ! -------------------------------------

    ! 1 source (used as the primary source, defined in the input file)
    CALL AddSources(sources_in, 1, (/0.0_r2,0.0_r2,0.0_r2/),                   &
                    T_star=model%t_star, R_star=model%r_star)
    ! 2 source (added by hand, could be an embedded planet or whatever)
    ! CALL AddSources(sources_in,2,(/sf*4.95_r2,-sf*0.72_r2,0.0_r2/),          &
    !                 T_star=1000.0, L_star=REAL(1e-4*L_sun,kind=r1))
    print '(A,I3)', " Number of radiation sources included: ",                 &
                      sources_in%n_sources
    print '(A,F5.2,A)'," Total Luminosity : ", sources_in%L_total /            &
                        L_sun, " L_sun"

    !----------------------------  Init Grid  ---------------------------------!
    ! we need to gerneralize the usage of this key TbD
    velo_type  = 1                    ! analytical velocity distribution
!~     velo_type  = 2                 ! lin. interpolated velocity distribution

    CALL parse('velo_type', velo_type, input_file)
    
    CALL parse('grid_type', grid_type, input_file)
    CALL parse('grid_name', grid_name, input_file)

    IF ( old_model ) THEN
        ! read from the old model boundary files
        file_a = TRIM(basics%path_results) //                                  &
                 TRIM(basics%pronam_old) // '_a_boundaries.dat'
        file_b = TRIM(basics%path_results) //                                  &
                 TRIM(basics%pronam_old) // '_b_boundaries.dat'
        file_c = TRIM(basics%path_results) //                                  &
                 TRIM(basics%pronam_old) // '_c_boundaries.dat'

        CALL read_no_cells(n_a, n_b, n_c, file_a, file_b, file_c)
        grid_type = 9
        sf = 1.0 ! not used
    ELSE IF (grid_type == 9) THEN
        ! grid is defined in input files see
        ! input/grid/a_coordinates.dat
        ! input/grid/b_coordinates.dat
        ! input/grid/c_coordinates.dat
        CALL read_no_cells(n_a, n_b, n_c)
        sf = 1 ! it is not needed
        velo_type  = 2
    ELSE
        ! this is the normal way, read the grid properties from the input file
        ! possible values: spherical (well tested)
        !                  cylindrical, cartesian (working but there might
        !                  be some bugs, test for yourself)
        !
        ! EXAMPLE:
        !
        ! spherical :  r = / 1 = rho  \
        !                  | 2 = theta|
        !                  \ 3 = phi  /
        ! n(1) = n_r :  number of point in r direction
        ! n(2) = n_th : number of point in th direction, must be odd
        ! n(3) = n_ph : number of point in oh direction, must be 1 or even
        ! sf   = step factor for logarithmic r scale (MC3D style)

        CALL parse('n_a',n_a, new_input_file)
        CALL parse('n_b',n_b, new_input_file)
        CALL parse('n_c',n_c, new_input_file)
        CALL parse('sf',sf, new_input_file)

    END IF

    ! write informations about the grid in the log/input file

    WRITE(unit=3,fmt='(A)') 'grid_name = {'//TRIM(grid_name)// &
        '}              possible values: spherical (OK), cylindrical (OK), cartesian (testing)'

    WRITE(help,fmt='(I1)') grid_type
    WRITE(unit=3,fmt='(A)') 'grid_type = {'//TRIM(help)// &
    '}                      default = 1, 2 and higher values for user defined grids'
    WRITE(help,fmt='(I1)') grid_type
    WRITE(unit=3,fmt='(A)') 'velo_type = {'//TRIM(help)// &
    '}                      1 : use analytical definition in model_mod.f90 (default)'
    WRITE(unit=3,fmt='(A)') '                          &
    &           2 : lin. interpolation is used, handy for full 3D data'
    WRITE(unit=3,fmt='(A)') ''
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
    '}                  stepfactor'
    WRITE(unit=3,fmt='(A)') ''

    ! give the type of the specified grid
    ! 1 == analytical version (prefered)
    !   Some examples here:
    ! SELECT CASE(grid_type)
    !
    ! CASE('cylindrical')
    !    grid_type = 1
    !    n_a  = 100
    !    sf   = 1.00001
    !    n_b  = 70
    !    n_c  = 101
    !
    ! CASE('spherical')
    !    grid_type = 1
    !    n_a  = 300
    !    sf   = 1.01
    !    n_b  = 257
    !    n_c  = 1
    !
    ! END SELECT
    CALL InitGrid(grid,grid_type,grid_name, n_a, sf, n_b, n_c,dust%n_dust, &
                  gas%egy_lvl)

    !---------------------------  Init Fluxes  --------------------------------!
    ! we always use Jy/pixel from now on, other units (e.g. T_b)
    ! can be calculated in the post process
    flux_unit = 'Jy_pix'

    CALL InitFluxes(fluxes,1, flux_unit, basics,                               &
                    model%n_bin_map,gas%n_tr, gas%i_vel_chan, dust%n_lam)

    print *, "Code initialization done !"
    Close(unit=3)

    !------------------------------  done!  -----------------------------------!
END SUBROUTINE inimol

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
    !--------------------------------------------------------------------------!

    CALL CloseBasic(basics)
    CALL CloseDust(dust)
    CALL CloseFluxes(fluxes)
    CALL CloseGrid(grid)
    CALL CloseModel(model)
    CALL CloseGas(gas)
    CALL CloseSources(sources_in)

END SUBROUTINE cleanup

END MODULE initiate
