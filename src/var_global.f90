module var_global
  use datatype

  implicit none
  PUBLIC
 
  ! ---
  ! parameters
  ! ---
  ! Msun      ... solar mass [kg]
  ! con_gamma ... grav. constant          [m^3 / (kg * s^2)]
  ! con_K     ... Boltzmann constant      [J/K]
  ! con_sigma ... stefan/boltzmann const. [W / (m^2 * K^4)] 
  ! con_mu    ... atomic mass unit        [kg]
  ! con_ c2   ... con_h*con_c/con_k
  ! real(kind=r3), public, parameter :: &

    real(kind=r2),parameter :: &
        PI        = 3.14159265358979324_r2, &
        con_h     = 6.626176e-34_r2, &
        con_k     = 1.380662e-23_r2, &
        con_Na    = 6.02214129e23_r2, &
        con_c     = 299792458.0_r2, &
        con_c2    = 1.43878631423230881E-2_r2, &
        con_AU    = 1.496e+11_r2, &
        con_pc    = 3.0856776e+16_r2,&
        con_gamma = 6.67384e-11_r2, &
        con_sigma = 5.6704e-8_r2, &
        con_mu    = 1.660538e-27_r2, &
        L_sun     = 3.85e+26_r2, &
        M_sun     = 1.9891e+30_r2, &
        R_sun     = 0.6960e+9, &
        SBK       = 5.67e-8_r2, &
        w2erg     = 1.0e+07_r2, &
        rel_err   = 1.0e-8, &   ! controls the rel error in the raytracing alg
        abs_err   = 1.0e-20     ! controls the abs error in the raytracing alg
    REAL(kind=r2),PARAMETER,DIMENSION(1:5)    :: &
        ! mass of collision partners (not final and needs an update)
        ! 1 = H2 = 2.01588 u
        ! 2 = He
        ! 3 = Electron
        !
        col_p_weight = (/2.01588_r2,0.0_r2,0.0_r2,0.0_r2,0.0_r2/)
        
    REAL(kind=r2),PARAMETER,DIMENSION(1:6)       ::  &
        RK_c  = (/0.0_r2, 1.0_r2/4.0_r2, 3.0_r2/8.0_r2, 12.0_r2/13.0_r2, 1.0_r2, 0.5_r2/), &
        RK_b1 = (/16.0_r2/135.0_r2, 0.0_r2, 6656.0_r2/12825.0_r2, 28561.0_r2/56430.0_r2, &
                -9.0_r2/50.0_r2, 2.0_r2/55.0_r2/), &
        RK_b2 = (/25.0_r2/216.0_r2, 0.0_r2, 1408.0_r2/2565.0_r2, 2197.0_r2/4104.0_r2, &
                -1.0_r2/5.0_r2, 0.0_r2/)
    
        ! Runge Kutta constants
    REAL(kind=r2),PARAMETER,DIMENSION(1:6,1:6)       ::  &   
        RK_a  = RESHAPE( (/ &
                0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, & 
                1.0_r2/4.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                3.0_r2/32.0_r2, 9.0_r2/32.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                1932.0_r2/2197.0_r2, -7200.0_r2/2197.0_r2, 7296.0_r2/2197.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                439.0_r2/216.0_r2, -8.0_r2, 3680.0_r2/513.0_r2, -845.0_r2/4104.0_r2, 0.0_r2, 0.0_r2, &
                -8.0_r2/27.0_r2, 2.0_r2, -3544.0_r2/2565.0_r2, 1859.0_r2/4104.0_r2, -11.0_r2/40.0_r2, &
                0.0_r2/), (/6,6/))
!~         RK_a  = RESHAPE( (/ &
!~                 0.0_r2, 1.0_r2/4.0_r2, 3.0_r2/32.0_r2, 1932.0_r2/2197.0_r2, 439.0_r2/216.0_r2, -8.0_r2/27.0_r2, &
!~                  
!~                 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                 , 9.0_r2/32.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                  -7200.0_r2/2197.0_r2, 7296.0_r2/2197.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                  -8.0_r2, 3680.0_r2/513.0_r2, -845.0_r2/4104.0_r2, 0.0_r2, 0.0_r2, &
!~                  2.0_r2, -3544.0_r2/2565.0_r2, 1859.0_r2/4104.0_r2, -11.0_r2/40.0_r2, &
!~                 0.0_r2/), (/6,6/)))

  character(len=75), public, parameter :: &
       hrule = "---------------------------------------------------------------------------"
  character(len= 5), public, parameter :: path_misc        = "misc/"
  character(len=11), public, parameter :: path_dust_cat    = "input/dust/"
  character(len=18), public, parameter :: path_dust_single = "input/dust/single/"
  character(len=21), public, parameter :: path_dust_nk     = "input/dust/tables_nk/"
  character(len=21), public, parameter :: path_dust_wv     = "input/dust/tables_wv/"
  character(len=15), public, parameter :: path_dust_tmp    = "input/dust/tmp/"
  character(len=12), public, parameter :: path_model       = "input/model/"
  character(len=13), public, parameter :: path_extern      = "input/extern/"
  character(len=10), public, parameter :: path_mol         = "input/mol/"
  
  ! ---
  ! strings
  ! ---
  character(len=8), public, parameter :: mol3d_version    = "Feb 2015"

  ! ---
  ! integer
  ! ---
  integer, public :: velo_type

  integer, public :: n_interact_max

  ! ---
  ! logical
  ! ---
  logical, public :: show_error


end module var_global
